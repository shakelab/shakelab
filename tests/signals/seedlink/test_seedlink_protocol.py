# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
Tests of :mod:`shakelab.signals.seedlink.protocol`.
"""

from __future__ import annotations

import socket
import struct

import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.seedlink import protocol
from shakelab.signals.seedlink.protocol import Connection, SeedLinkError

OGS_HELLO = (
    "SeedLink v3.1 (2017.052 RingServer) :: SLPROTO:3.1 CAP EXTREPLY "
    "NSWILDCARD BATCH WS:13"
)


# ----------------------------------------------------------------------------
# Negotiation
# ----------------------------------------------------------------------------


def test_parse_hello_v3():
    versions, capabilities = protocol.parse_hello(OGS_HELLO)
    assert versions == ["3.1"]
    assert {"NSWILDCARD", "EXTREPLY", "BATCH"} <= capabilities


def test_parse_hello_v4_and_banner_only():
    versions, capabilities = protocol.parse_hello(
        "SeedLink v4.0 (RingServer) :: SLPROTO:4.0 SLPROTO:3.1 TIME"
    )
    assert versions == ["4.0", "3.1"]
    assert capabilities == {"TIME"}
    assert protocol.parse_hello("SeedLink v3.0 (old server)")[0] == ["3.0"]


@pytest.mark.parametrize("offered, requested, expected", [
    (["4.0", "3.1"], "auto", "4.0"),
    (["3.1"], "auto", "3.1"),
    (["3.0"], "auto", "3.1"),
    (["4.0", "3.1"], "3.1", "3.1"),
    (["4.1"], "4.0", "4.0"),
])
def test_choose_version(offered, requested, expected):
    assert protocol.choose_version(offered, requested) == expected


@pytest.mark.parametrize("offered, requested", [
    (["3.1"], "4.0"),
    (["2.5"], "auto"),
    ([], "auto"),
])
def test_choose_version_unsupported(offered, requested):
    with pytest.raises(SeedLinkError) as error:
        protocol.choose_version(offered, requested)
    assert error.value.code == "UNSUPPORTED"


def test_normalize_protocol():
    assert protocol.normalize_protocol("AUTO") == "auto"
    assert protocol.normalize_protocol("4") == "4.0"
    assert protocol.normalize_protocol(3) == "3.1"
    with pytest.raises(ValueError):
        protocol.normalize_protocol("2.5")
    with pytest.raises(TypeError):
        protocol.normalize_protocol(None)


# ----------------------------------------------------------------------------
# Commands
# ----------------------------------------------------------------------------


def test_station_command():
    assert protocol.station_command("3.1", "IV", "MODE") == (
        "STATION MODE IV"
    )
    assert protocol.station_command("4.0", "IV", "MO*") == (
        "STATION IV_MO*"
    )


@pytest.mark.parametrize("selectors, v3, v4", [
    ([("*", "HHZ")], ["SELECT HHZ.D"], ["SELECT *_H_H_Z"]),
    ([("00", "HH?")], ["SELECT 00HH?.D"], ["SELECT 00_H_H_?"]),
    ([("", "HNZ")], ["SELECT HNZ.D"], ["SELECT _H_N_Z"]),
    ([("00", "*")], ["SELECT 00???.D"], ["SELECT 00_*"]),
    ([("*", "*")], [], []),
    ([("*", "H*")], [], []),
    ([("*", "HHZ"), ("*", "*")], [], []),
    ([("*", "HHZ"), ("*", "HHZ")], ["SELECT HHZ.D"], ["SELECT *_H_H_Z"]),
])
def test_select_commands(selectors, v3, v4):
    assert protocol.select_commands("3.1", selectors) == v3
    assert protocol.select_commands("4.0", selectors) == v4


def test_data_command_v3():
    start = Date("2026-11-05T08:00:00.25").to_unix_ns()
    end = Date("2026-11-05T08:10:00.5").to_unix_ns()
    assert protocol.data_command("3.1") == "DATA"
    assert protocol.data_command("3.1", fetch=True) == "FETCH"
    assert protocol.data_command("3.1", sequence=0x1A2B) == "DATA 001A2B"
    assert protocol.data_command("3.1", sequence=2 ** 24 + 1) == (
        "DATA 000001"
    )
    assert protocol.data_command("3.1", start_ns=start, end_ns=end) == (
        "TIME 2026,11,05,08,00,00 2026,11,05,08,10,01"
    )
    assert protocol.data_command(
        "3.1", sequence=5, start_ns=start, end_ns=end, fetch=True,
    ) == "FETCH 000005 2026,11,05,08,00,00"


def test_data_command_v4():
    start = Date("2026-11-05T08:00:00.25").to_unix_ns()
    end = Date("2026-11-05T08:10:00").to_unix_ns()
    assert protocol.data_command("4.0") == "DATA"
    assert protocol.data_command("4.0", sequence=123456789012) == (
        "DATA 123456789012"
    )
    assert protocol.data_command("4.0", start_ns=start, end_ns=end) == (
        "DATA ALL 2026-11-05T08:00:00.250000Z 2026-11-05T08:10:00.000000Z"
    )
    assert protocol.data_command("4.0", sequence=7, start_ns=start) == (
        "DATA 7 2026-11-05T08:00:00.250000Z"
    )


def test_start_command():
    assert protocol.start_command("3.1", fetch=True) == "END"
    assert protocol.start_command("4.0", fetch=True) == "ENDFETCH"
    assert protocol.start_command("4.0") == "END"


def test_sequence_order():
    assert protocol.next_sequence("3.1", 0xFFFFFF) == 0
    assert protocol.next_sequence("4.0", 0xFFFFFF) == 0x1000000
    assert protocol.is_newer("3.1", 0, 0xFFFFFF)
    assert not protocol.is_newer("3.1", 0xFFFFFF, 0)
    assert not protocol.is_newer("3.1", 5, 5)
    assert protocol.is_newer("4.0", 6, 5)
    assert not protocol.is_newer("4.0", 5, 5)


# ----------------------------------------------------------------------------
# Packets
# ----------------------------------------------------------------------------


def _connection(version):
    """Connection on one end of a socket pair, and the other end."""
    left, right = socket.socketpair()
    connection = Connection("localhost")
    connection._socket = left
    connection.version = version
    return connection, right


def test_read_packet_v3():
    connection, peer = _connection("3.1")
    record = bytes(range(256)) * 2
    peer.sendall(b"SL00002A" + record + b"SLINFO *" + record + b"END")

    packet = connection.read_packet()
    assert (packet.kind, packet.sequence, packet.format) == (
        "data", 42, "2"
    )
    assert packet.payload == record

    info = connection.read_packet()
    assert info.kind == "info" and not info.last
    assert connection.read_packet() is None
    peer.close()
    with pytest.raises(ConnectionError):
        connection.read_packet()


def test_read_packet_v4():
    connection, peer = _connection("4.0")
    payload = b"x" * 100
    header = struct.pack("<2sccIQB", b"SE", b"2", b"D", 100, 2 ** 40, 7)
    peer.sendall(header + b"IV_MODE" + payload)
    json_info = b'{"software": "x"}'
    peer.sendall(
        struct.pack("<2sccIQB", b"SE", b"J", b"I", len(json_info), 0, 0)
        + json_info
    )

    packet = connection.read_packet()
    assert packet.kind == "data"
    assert packet.sequence == 2 ** 40
    assert packet.station == "IV_MODE"
    assert packet.payload == payload
    info = connection.read_packet()
    assert info.kind == "info" and info.payload == json_info


def test_partial_packet_survives_timeout():
    connection, peer = _connection("3.1")
    connection.set_poll(0.05)
    record = b"r" * 512
    peer.sendall(b"SL000001" + record[:100])
    with pytest.raises(TimeoutError):
        connection.read_packet()
    peer.sendall(record[100:])
    assert connection.read_packet().payload == record


def test_invalid_packet():
    connection, peer = _connection("3.1")
    peer.sendall(b"XX000001" + b"\0" * 512)
    with pytest.raises(SeedLinkError):
        connection.read_packet()


def test_command_replies():
    connection, peer = _connection("4.0")
    peer.sendall(b"OK\r\nERROR ARGUMENTS bad pattern\r\nWHAT\r\n")
    connection.command("SELECT *")
    with pytest.raises(SeedLinkError) as error:
        connection.command("SELECT !")
    assert error.value.code == "ARGUMENTS"
    assert "bad pattern" in str(error.value)
    with pytest.raises(SeedLinkError):
        connection.command("SELECT x")
    with pytest.raises(ValueError):
        connection.send("SELECT " + "x" * 300)
