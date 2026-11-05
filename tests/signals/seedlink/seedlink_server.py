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
Local SeedLink server for the tests of :mod:`shakelab.signals.seedlink`
(fixture ``make_server`` in ``conftest.py``).

:class:`FakeServer` speaks version 3.1 or 4.0 of the protocol on a local
port, with a fixed set of MiniSEED packets per station. It records the
commands received, honours sequence numbers and time windows, ignores
``SELECT`` (the client filters the records), and can drop the first
connection after some packets to test reconnection.
"""

from __future__ import annotations

import fnmatch
import io
import json
import socket
import struct
import threading

import numpy as np
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.libio import mseed

START = Date("2026-11-05T08:00:00")
RATE = 100.0

HELLO_V3 = (
    "SeedLink v3.1 (2017.052 RingServer) :: SLPROTO:3.1 CAP EXTREPLY "
    "NSWILDCARD BATCH WS:13"
)
HELLO_V3_PLAIN = "SeedLink v3.1 (2020.075 SeisComP) :: SLPROTO:3.1"
HELLO_V4 = (
    "SeedLink v4.0 (2025.001 RingServer) :: SLPROTO:4.0 SLPROTO:3.1 TIME"
)


# ----------------------------------------------------------------------------
# Test data
# ----------------------------------------------------------------------------


def waveform(sid, nsamp, offset=0.0, seed=0):
    """Integer test signal of one channel."""
    rng = np.random.default_rng(seed)
    data = np.cumsum(rng.integers(-500, 500, nsamp)).astype(np.int32)
    return Record(
        data=data, delta=1.0 / RATE, time=START + offset, sid=sid
    )


def mseed_packets(records):
    """512-byte MiniSEED 2 records of a list of Records."""
    collection = StreamCollection()
    for record in records:
        collection.append(record)
    output = io.BytesIO()
    mseed.mswrite(output, collection, reclen=512)
    content = output.getvalue()
    return [content[i:i + 512] for i in range(0, len(content), 512)]


def log_records(text):
    """MiniSEED 2 log records (encoding 0) carrying ``text``."""
    blocks = []
    data = text.encode("ascii")
    size = 512 - 56
    chunks = [data[i:i + size] for i in range(0, len(data), size)] or [b""]
    for chunk in chunks:
        header = struct.pack(
            ">6sss5s2s3s2sHHBBBBHHhhBBBBiHH",
            b"000001", b"D", b" ", b"INFO ", b"  ", b"LOG", b"SL",
            2026, 309, 8, 0, 0, 0, 0, len(chunk), 0, 0, 0, 0, 0, 1, 0,
            56, 48,
        )
        blockette = struct.pack(">HHBBBB", 1000, 0, 0, 1, 9, 0)
        record = header + blockette + chunk
        blocks.append(record + b"\0" * (512 - len(record)))
    return blocks


class Station:
    """Packets of one station: (sequence, start_ns, end_ns, bytes)."""

    def __init__(self, network, station, records, first_sequence=100):
        self.network = network
        self.station = station
        self.packets = []
        sequence = first_sequence
        for block in mseed_packets(records):
            msr = mseed.MSRecord.from_bytes(block)
            start = msr.start_ns
            end = start + int(msr.nsamp * 1e9 / msr.rate)
            self.packets.append((sequence, start, end, block))
            sequence += 1

    @property
    def key(self):
        return f"{self.network}_{self.station}"


# ----------------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------------


class FakeServer:
    """
    SeedLink server on a local port, run in a background thread.

    Parameters
    ----------
    stations : list of Station
    version : {"3.1", "4.0"}
    hello : str, optional
        First line of the HELLO response.
    drop_after : int, optional
        Close the first connection after sending this many packets.
    reject : str, optional
        Command prefix answered with an error.
    info_error : bool, optional
        Answer the INFO requests with an error.
    """

    def __init__(self, stations, version="3.1", hello=None,
                 drop_after=None, reject=None, info_error=False):
        self.stations = {item.key: item for item in stations}
        self.version = version
        self.hello = hello or (HELLO_V4 if version == "4.0" else HELLO_V3)
        self.drop_after = drop_after
        self.reject = reject
        self.info_error = info_error
        self.commands = []
        self.connections = 0
        self._listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._listener.bind(("127.0.0.1", 0))
        self._listener.listen(5)
        self._listener.settimeout(0.2)
        self.port = self._listener.getsockname()[1]
        self._closing = threading.Event()
        self._thread = threading.Thread(target=self._serve, daemon=True)
        self._thread.start()

    def close(self):
        self._closing.set()
        self._thread.join(5)
        self._listener.close()

    # -- serving ---------------------------------------------------------

    def _serve(self):
        while not self._closing.is_set():
            try:
                conn, _ = self._listener.accept()
            except socket.timeout:
                continue
            self.connections += 1
            try:
                self._session(conn, self.connections)
            except OSError:
                pass
            finally:
                conn.close()

    def _session(self, conn, number):
        conn.settimeout(0.2)
        reader = _LineReader(conn, self._closing)
        requests = []
        current = None
        fetch = False
        v4 = self.version == "4.0"

        while True:
            line = reader.line()
            if line is None:
                return
            if not line:
                continue
            self.commands.append(line)
            word = line.split()[0].upper()
            args = line.split()[1:]

            if self.reject and line.startswith(self.reject):
                conn.sendall(b"ERROR ARGUMENTS rejected for the test\r\n")
            elif word == "HELLO":
                conn.sendall(
                    f"{self.hello}\r\nTest Seismic Network\r\n".encode()
                )
            elif word in ("SLPROTO", "USERAGENT", "SELECT"):
                conn.sendall(b"OK\r\n")
            elif word == "STATION":
                if v4:
                    network, _, station = args[0].partition("_")
                else:
                    station = args[0]
                    network = args[1] if len(args) > 1 else "*"
                current = {
                    "network": network, "station": station,
                    "sequence": None, "start": None, "end": None,
                }
                requests.append(current)
                conn.sendall(b"OK\r\n")
            elif word in ("DATA", "FETCH", "TIME"):
                fetch = fetch or word == "FETCH"
                _parse_data(word, args, current, v4)
                conn.sendall(b"OK\r\n")
            elif word == "INFO":
                self._send_info(conn, args[0] if args else "ID")
            elif word in ("END", "ENDFETCH"):
                fetch = fetch or word == "ENDFETCH"
                self._transfer(conn, reader, requests, fetch, number)
                return
            elif word == "BYE":
                return
            else:
                conn.sendall(b"ERROR UNSUPPORTED\r\n")

    def _transfer(self, conn, reader, requests, fetch, number):
        sent = 0
        done = set()
        for request in requests:
            for key, station in self.stations.items():
                if key in done:
                    continue
                if not (
                    fnmatch.fnmatchcase(station.network, request["network"])
                    and fnmatch.fnmatchcase(station.station,
                                            request["station"])
                ):
                    continue
                done.add(key)
                for sequence, start, end, block in station.packets:
                    if (
                        request["sequence"] is not None
                        and sequence < request["sequence"]
                    ):
                        continue
                    if request["start"] is not None and end <= (
                        request["start"]
                    ):
                        continue
                    if request["end"] is not None and start >= (
                        request["end"]
                    ):
                        continue
                    if (
                        number == 1
                        and self.drop_after is not None
                        and sent >= self.drop_after
                    ):
                        return
                    conn.sendall(self._packet(station, sequence, block))
                    sent += 1

        if fetch:
            conn.sendall(b"END")
        # Real time: keep the connection until the client leaves,
        # answering keep-alive requests.
        while True:
            line = reader.line()
            if line is None:
                return
            self.commands.append(line)
            if line.upper().startswith("BYE"):
                return
            if line.upper().startswith("INFO"):
                self._send_info(conn, "ID")

    def _packet(self, station, sequence, block):
        if self.version == "4.0":
            key = station.key.encode()
            header = struct.pack(
                "<2sccIQB", b"SE", b"2", b"D", len(block), sequence,
                len(key),
            )
            return header + key + block
        return f"SL{sequence % 2 ** 24:06X}".encode() + block

    def _streams(self, station):
        """(location, channel, start ns, end ns) of a station's data."""
        spans = {}
        for _, start, end, block in station.packets:
            msr = mseed.MSRecord.from_bytes(block)
            key = (
                msr.header["LOCATION_IDENTIFIER"].strip(),
                msr.header["CHANNEL_IDENTIFIER"].strip(),
            )
            first, last = spans.get(key, (start, end))
            spans[key] = (min(first, start), max(last, end))
        return [key + span for key, span in sorted(spans.items())]

    def _send_info(self, conn, item):
        if self.info_error:
            item = "ERROR"
        if self.version == "4.0":
            self._send_info_v4(conn, item)
        else:
            self._send_info_v3(conn, item)

    def _send_info_v4(self, conn, item):
        content = {"software": "fake", "organization": "test"}
        if item == "ERROR":
            content["error"] = {"code": "UNAUTHORIZED",
                                "message": "not allowed"}
        elif item in ("STATIONS", "STREAMS"):
            content["station"] = []
            for key, station in self.stations.items():
                entry = {
                    "id": key, "description": f"Station {station.station}",
                    "start_seq": station.packets[0][0],
                    "end_seq": station.packets[-1][0] + 1,
                }
                if item == "STREAMS":
                    entry["stream"] = [
                        {"id": f"{loc}_{'_'.join(cha)}", "format": "2",
                         "subformat": "D", "start_time": _iso(start),
                         "end_time": _iso(end)}
                        for loc, cha, start, end in self._streams(station)
                    ] + [{"id": "_L_O_G", "format": "2", "subformat": "L",
                          "start_time": _iso(0), "end_time": _iso(0)}]
                content["station"].append(entry)
        payload = json.dumps(content).encode()
        subformat = b"E" if item == "ERROR" else b"I"
        conn.sendall(
            struct.pack("<2sccIQB", b"SE", b"J", subformat, len(payload),
                        0, 0)
            + payload
        )

    def _send_info_v3(self, conn, item):
        if item == "ERROR":
            text = '<?xml version="1.0"?><seedlink><error message=' \
                '"not allowed"/></seedlink>'
        else:
            body = []
            for station in self.stations.values():
                streams = ""
                if item == "STREAMS":
                    streams = "".join(
                        f'<stream location="{loc or "  "}" '
                        f'seedname="{cha}" type="D" '
                        f'begin_time="{_v3_info_time(start)}" '
                        f'end_time="{_v3_info_time(end)}"/>'
                        for loc, cha, start, end in self._streams(station)
                    ) + '<stream location="  " seedname="LOG" type="L"/>'
                body.append(
                    f'<station name="{station.station}" '
                    f'network="{station.network}" '
                    f'description="Station {station.station}" '
                    f'begin_seq="{station.packets[0][0]:06X}" '
                    f'end_seq="{station.packets[-1][0]:06X}">'
                    f"{streams}</station>"
                )
            text = (
                '<?xml version="1.0"?><seedlink software="fake">'
                + "".join(body) + "</seedlink>"
            )
        blocks = log_records(text)
        for index, block in enumerate(blocks):
            flag = b"*" if index < len(blocks) - 1 else b" "
            conn.sendall(b"SLINFO " + flag + block)


def _iso(ns):
    return Date.from_unix_ns(ns).iso8601


def _v3_info_time(ns):
    return Date.from_unix_ns(ns).strformat("%Y/%m/%d %H:%M:%S")


class _LineReader:
    """Line reader of a server socket with a short timeout."""

    def __init__(self, conn, closing):
        self.conn = conn
        self.closing = closing
        self.buffer = b""

    def line(self):
        while b"\n" not in self.buffer and b"\r" not in self.buffer:
            if self.closing.is_set():
                return None
            try:
                chunk = self.conn.recv(4096)
            except socket.timeout:
                continue
            if not chunk:
                return None
            self.buffer += chunk
        for terminator in (b"\r\n", b"\n", b"\r"):
            if terminator in self.buffer:
                line, self.buffer = self.buffer.split(terminator, 1)
                return line.decode().strip()
        return None


def _parse_data(word, args, request, v4):
    """Sequence number and time window of a DATA/FETCH/TIME command."""
    if word == "TIME":
        request["start"] = _v3_ns(args[0])
        if len(args) > 1:
            request["end"] = _v3_ns(args[1])
        return
    if not args:
        return
    if v4:
        if args[0] != "ALL":
            request["sequence"] = int(args[0])
        if len(args) > 1:
            request["start"] = Date(args[1]).to_unix_ns()
        if len(args) > 2:
            request["end"] = Date(args[2]).to_unix_ns()
        return
    request["sequence"] = int(args[0], 16)
    if len(args) > 1:
        request["start"] = _v3_ns(args[1])


def _v3_ns(text):
    year, month, day, hour, minute, second = (
        int(value) for value in text.split(",")
    )
    return Date([year, month, day, hour, minute, second]).to_unix_ns()
