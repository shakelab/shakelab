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
Tests of the selection from files and metadata, the server information,
the channel status and the command line of
:mod:`shakelab.signals.seedlink`, against the local server of
``conftest.py``.
"""

from __future__ import annotations

import pytest

from seedlink_server import START, Station, waveform
from shakelab.libutils.timeN import Date
from shakelab.signals.libio.mseed import msread
from shakelab.signals.metadata import ChannelEpoch, MetadataCollection
from shakelab.signals.seedlink import (
    ChannelStatus,
    Client,
    SeedLinkError,
    StationInfo,
)
from shakelab.signals.seedlink.cli import main

VERSIONS = ["3.1", "4.0"]


@pytest.fixture
def stations():
    """IV.MODE (two channels, 60 s), OX.ABC and OX.BAD (one, 30 s)."""
    mode = [
        waveform("IV.MODE..HNZ", 6000, seed=1),
        waveform("IV.MODE..HNE", 6000, seed=2),
    ]
    abc = [waveform("OX.ABC.00.HHZ", 3000, seed=3)]
    bad = [waveform("OX.BAD..HHZ", 3000, seed=4)]
    return [
        Station("IV", "MODE", mode),
        Station("OX", "ABC", abc, first_sequence=7),
        Station("OX", "BAD", bad, first_sequence=0xFFFFF0),
    ]


def _client(server, **kwargs):
    kwargs.setdefault("timeout", 5.0)
    return Client("127.0.0.1", port=server.port, **kwargs)


# ----------------------------------------------------------------------------
# Selection
# ----------------------------------------------------------------------------


def test_select_braces(make_server, stations):
    server = make_server(stations)
    client = _client(server)
    client.select("{IV.MODE..HNZ,OX.{ABC,BAD}}")
    collection = client.run(fetch=True, duration=10)
    assert collection.sid == [
        "IV.MODE..HNZ", "OX.ABC.00.HHZ", "OX.BAD..HHZ",
    ]
    assert "STATION ABC OX" in server.commands


def test_select_file(make_server, stations, tmp_path):
    path = tmp_path / "channels.txt"
    path.write_text(
        "# test channels\nIV.MODE..HNE   OX.BAD # inline comment\n\n"
    )
    server = make_server(stations)
    client = _client(server)
    client.select(file=path)
    assert client.selection == ["IV.MODE..HNE", "OX.BAD.*.*"]
    collection = client.run(fetch=True, duration=10)
    assert collection.sid == ["IV.MODE..HNE", "OX.BAD..HHZ"]

    (tmp_path / "empty.txt").write_text("# nothing\n")
    with pytest.raises(ValueError):
        client.select(file=tmp_path / "empty.txt")


def test_select_metadata():
    metadata = MetadataCollection()
    metadata.append(ChannelEpoch(sid="IV.MODE..HNZ", starttime=START))
    metadata.append(ChannelEpoch(
        sid="IV.MODE..HHZ", starttime="2000-01-01T00:00:00",
        endtime="2010-01-01T00:00:00",
    ))
    client = Client("localhost")
    client.select(metadata)
    assert client.selection == ["IV.MODE..HNZ"]

    closed = MetadataCollection()
    closed.append(ChannelEpoch(
        sid="IV.MODE..HHZ", starttime="2000-01-01T00:00:00",
        endtime="2010-01-01T00:00:00",
    ))
    with pytest.raises(ValueError, match="valid now"):
        client.select(closed)
    with pytest.raises(ValueError):
        client.select()
    with pytest.raises(TypeError):
        client.select({"IV.MODE"})


# ----------------------------------------------------------------------------
# Server information
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
def test_stations(make_server, stations, version):
    server = make_server(stations, version=version)
    client = _client(server)

    result = client.stations()
    assert [station.code for station in result] == [
        "IV.MODE", "OX.ABC", "OX.BAD",
    ]
    assert result[1] == StationInfo(
        "OX", "ABC", "Station ABC", 7,
        stations[1].packets[-1][0] + (version == "4.0"),
    )
    assert [s.code for s in client.stations("OX.{BAD,XYZ}")] == ["OX.BAD"]
    assert client.stations("XX") == []


@pytest.mark.parametrize("version", VERSIONS)
def test_streams(make_server, stations, version):
    server = make_server(stations, version=version)
    client = _client(server)

    streams = client.streams()
    assert [stream.sid for stream in streams] == [
        "IV.MODE..HNE", "IV.MODE..HNZ", "OX.ABC.00.HHZ", "OX.BAD..HHZ",
    ]
    abc = client.streams("OX.ABC")[0]
    assert abc.start == START
    assert abs(abc.end - (START + 30.0)) < 1e-3
    assert [s.sid for s in client.streams("*.*.*.HNZ")] == [
        "IV.MODE..HNZ"
    ]


@pytest.mark.parametrize("version", VERSIONS)
def test_info_error(make_server, stations, version):
    server = make_server(stations, version=version, info_error=True)
    with pytest.raises(SeedLinkError, match="not allowed"):
        _client(server).streams()


# ----------------------------------------------------------------------------
# Channel status
# ----------------------------------------------------------------------------


def test_status(make_server, capsys):
    records = [
        waveform("IV.MODE..HNZ", 2000),
        waveform("IV.MODE..HNZ", 2000, offset=30.0, seed=1),
    ]
    server = make_server([Station("IV", "MODE", records)])
    client = _client(server)
    client.select("IV.MODE")
    assert client.status() == []

    client.run(fetch=True, duration=10)

    status, = client.status()
    assert isinstance(status, ChannelStatus)
    assert status.sid == "IV.MODE..HNZ"
    assert status.last_sample == START + 30.0 + 19.99
    assert status.gaps == 1 and status.overlaps == 0
    assert status.packets == len(
        server.stations["IV_MODE"].packets
    )
    assert status.latency == pytest.approx(
        Date.now() - status.last_sample, abs=5.0
    )
    assert abs(status.last_packet - Date.now()) < 60
    assert client.status("OX") == []

    client.info()
    output = capsys.readouterr().out
    assert "IV.MODE..HNZ" in output and "latency" in output


# ----------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------


def test_cli_info(make_server, stations, capsys):
    server = make_server(stations)
    target = f"127.0.0.1:{server.port}"
    assert main(["info", target]) == 0
    assert "OX.ABC     Station ABC" in capsys.readouterr().out
    assert main(["info", target, "OX.ABC", "--streams"]) == 0
    output = capsys.readouterr().out
    assert output.startswith("OX.ABC.00.HHZ")
    assert "2026-11-05T08:00:00" in output


def test_cli_monitor(make_server, stations, capsys):
    server = make_server(stations)
    target = f"127.0.0.1:{server.port}"
    code = main([
        "monitor", target, "OX.ABC", "--duration", "1", "--timeout", "5",
    ])
    assert code == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("OX.ABC.00.HHZ    2026-11-05T08:00:00")
    assert len(lines) == len(stations[1].packets)

    code = main([
        "monitor", target, "OX.ABC", "--duration", "1.2", "--status",
        "0.5",
    ])
    assert code == 0
    assert "latency" in capsys.readouterr().out


def test_cli_get_and_archive(make_server, stations, tmp_path):
    server = make_server(stations)
    target = f"127.0.0.1:{server.port}"
    output = tmp_path / "abc.mseed"
    code = main([
        "get", target, "OX.ABC", "--start", "2026-11-05T08:00:10",
        "--end", "2026-11-05T08:00:20", "-o", str(output),
    ])
    assert code == 0
    record = msread(str(output))["OX.ABC.00.HHZ"].record[0]
    assert record.time <= START + 10.0
    assert record.endtime >= START + 20.0

    assert main([
        "get", target, "OX.ABC", "--duration", "1", "-o", str(output),
    ]) == 1

    code = main([
        "archive", target, "OX.ABC", "--sds", str(tmp_path / "sds"),
        "--state", str(tmp_path / "state"), "--duration", "1",
    ])
    assert code == 0
    assert (tmp_path / "sds/2026/OX/ABC/HHZ.D/OX.ABC.00.HHZ.D.2026.309"
            ).exists()
    assert (tmp_path / "state").exists()

    code = main([
        "archive", target, "OX.ABC", "--path",
        str(tmp_path / "rot" / "{sta}.{time:%H%M%S}.mseed"),
        "--rotate", "10", "--duration", "1",
    ])
    assert code == 0
    assert len(list((tmp_path / "rot").iterdir())) == 3


def test_cli_errors(capsys):
    assert main(["info", "127.0.0.1:x"]) == 1
    assert "HOST[:PORT]" in capsys.readouterr().err
    assert main(["monitor", "127.0.0.1:1", "--duration", "1"]) == 1
    assert "patterns" in capsys.readouterr().err
    assert main(["get", "127.0.0.1:1", "IV", "-o", "x.mseed"]) == 1
    assert "--end or --duration" in capsys.readouterr().err
    assert main(["info", "127.0.0.1:1"]) == 1
    assert "127.0.0.1:1" in capsys.readouterr().err
