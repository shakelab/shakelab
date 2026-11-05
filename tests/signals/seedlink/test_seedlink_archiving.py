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
Tests of the archive, the state file and the :func:`archive` and
:func:`records` functions of :mod:`shakelab.signals.seedlink`, against
the local server of ``conftest.py``.
"""

from __future__ import annotations

import json
import threading

import numpy as np
import pytest

from seedlink_server import Station, waveform
from shakelab.signals.libio.mseed import msread
from shakelab.signals.seedlink import Client, archive, records

VERSIONS = ["3.1", "4.0"]


@pytest.fixture
def stations():
    """IV.MODE (two channels, 60 s) and OX.ABC (one channel, 30 s)."""
    mode = [
        waveform("IV.MODE..HNZ", 6000, seed=1),
        waveform("IV.MODE..HNE", 6000, seed=2),
    ]
    abc = [waveform("OX.ABC.00.HHZ", 3000, seed=3)]
    return (
        mode + abc,
        [Station("IV", "MODE", mode), Station("OX", "ABC", abc, 7)],
    )


def _client(server, **kwargs):
    kwargs.setdefault("timeout", 5.0)
    return Client("127.0.0.1", port=server.port, **kwargs)


def _packets(station, start=0):
    """Concatenated packets of a station."""
    return b"".join(packet[3] for packet in station.packets[start:])


# ----------------------------------------------------------------------------
# Archive
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
def test_single_file_archive(make_server, stations, tmp_path, version):
    records_, servers = stations
    server = make_server(servers, version=version)
    client = _client(server)
    client.select(["IV.MODE", "OX.ABC"])
    client.set_archive(tmp_path / "all.mseed")

    assert client.run(fetch=True, duration=10) is None

    content = (tmp_path / "all.mseed").read_bytes()
    assert content == _packets(servers[0]) + _packets(servers[1])
    collection = msread(str(tmp_path / "all.mseed"))
    for record in records_:
        stored = collection[record.sid].record
        assert len(stored) == 1
        np.testing.assert_array_equal(stored[0].data, record.data)


def test_archive_by_channel_and_selection(make_server, stations, tmp_path):
    server = make_server(stations[1])
    client = _client(server)
    client.select("IV.MODE..HNZ")
    client.set_archive(tmp_path / "{net}.{sta}.{loc}.{cha}.mseed")
    client.set_buffer()

    collection = client.run(fetch=True, duration=10)

    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "IV.MODE..HNZ.mseed"
    ]
    archived = msread(str(tmp_path / "IV.MODE..HNZ.mseed"))
    np.testing.assert_array_equal(
        archived["IV.MODE..HNZ"].record[0].data,
        collection["IV.MODE..HNZ"].record[0].data,
    )


def test_rotated_and_sds_archives(make_server, stations, tmp_path):
    server = make_server(stations[1])
    client = _client(server)
    client.select("OX.ABC")
    client.set_archive(
        tmp_path / "rot" / "{sta}.{time:%H%M%S}.mseed", rotate="10s",
    )
    client.run(fetch=True, duration=10)

    names = sorted(path.name for path in (tmp_path / "rot").iterdir())
    assert names[0] == "ABC.080000.mseed"
    assert len(names) == 3
    joined = b"".join(
        (tmp_path / "rot" / name).read_bytes() for name in names
    )
    assert joined == _packets(stations[1][1])

    client = _client(server)
    client.select("OX.ABC")
    client.set_archive(tmp_path / "sds", layout="sds")
    client.run(fetch=True, duration=10)
    path = tmp_path / "sds/2026/OX/ABC/HHZ.D/OX.ABC.00.HHZ.D.2026.309"
    assert path.read_bytes() == _packets(stations[1][1])


def test_archive_function(make_server, stations, tmp_path):
    server = make_server(stations[1])
    archive(
        "127.0.0.1", "OX.ABC", tmp_path / "{sta}.mseed", port=server.port,
        timeout=5.0, duration=1.0, state_file=tmp_path / "state.json",
    )
    assert (tmp_path / "ABC.mseed").read_bytes() == _packets(
        stations[1][1]
    )
    state = json.loads((tmp_path / "state.json").read_text())
    assert state["sequence"] == {"OX_ABC": stations[1][1].packets[-1][0]}


# ----------------------------------------------------------------------------
# State file
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
def test_state_file_resumes(make_server, stations, tmp_path, version):
    server = make_server(stations[1], version=version, drop_after=None)
    state = tmp_path / "run" / "state.json"
    mode = stations[1][0]

    first = _client(server, state_file=state)
    first.select("IV.MODE")
    first.set_buffer()
    first.run(fetch=True, duration=10)
    content = json.loads(state.read_text())
    assert content["server"] == f"127.0.0.1:{server.port}"
    assert content["protocol"] == version
    assert content["sequence"]["IV_MODE"] == mode.packets[-1][0]

    # A new client starts after the packets already received.
    content["sequence"]["IV_MODE"] = mode.packets[-4][0]
    state.write_text(json.dumps(content))
    second = _client(server, state_file=state)
    second.select("IV.MODE")
    second.set_archive(tmp_path / "rest.mseed")
    second.run(fetch=True, duration=10)
    assert (tmp_path / "rest.mseed").read_bytes() == _packets(mode, -3)
    resumed = mode.packets[-3][0]
    expected = (
        f"DATA {resumed}" if version == "4.0" else f"FETCH {resumed:06X}"
    )
    assert expected in server.commands


def test_state_file_of_other_server_ignored(make_server, stations,
                                            tmp_path, caplog):
    server = make_server(stations[1])
    state = tmp_path / "state.json"
    state.write_text(json.dumps({
        "format": "shakelab-seedlink-state", "version": 1,
        "server": "elsewhere:18000", "protocol": "3.1",
        "sequence": {"IV_MODE": 10 ** 6},
    }))
    client = _client(server, state_file=state)
    client.select("IV.MODE")
    collection = client.run(fetch=True, duration=10)
    assert "belongs to elsewhere:18000" in caplog.text
    assert collection["IV.MODE..HNZ"].nsamp == 6000


@pytest.mark.parametrize("text", ["{", "[]", '{"format": "other"}'])
def test_invalid_state_file(stations, tmp_path, text):
    state = tmp_path / "state.json"
    state.write_text(text)
    client = Client("127.0.0.1", port=1, state_file=state)
    client.select("IV.MODE")
    with pytest.raises(ValueError, match="invalid state file"):
        client.run(duration=1)


def test_state_file_type():
    with pytest.raises(TypeError):
        Client("localhost", state_file=5)


# ----------------------------------------------------------------------------
# Iterator
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
def test_records_iterator(make_server, stations, version):
    records_, servers = stations
    server = make_server(servers, version=version)

    received = list(records(
        "127.0.0.1", ["IV.MODE..HNZ", "OX.ABC"], port=server.port,
        timeout=5.0, fetch=True, max_queue=2,
    ))

    total = {}
    for record in received:
        total[record.sid] = total.get(record.sid, 0) + record.nsamp
    assert total == {"IV.MODE..HNZ": 6000, "OX.ABC.00.HHZ": 3000}


def test_records_iterator_left_early(make_server, stations):
    server = make_server(stations[1])
    iterator = records(
        "127.0.0.1", "IV.MODE", port=server.port, timeout=5.0,
        max_queue=1,
    )
    first = next(iterator)
    assert first.sid == "IV.MODE..HNZ"
    iterator.close()

    # The background client has ended and the iterator is exhausted.
    assert not any(
        thread.name == "seedlink-127.0.0.1"
        for thread in threading.enumerate()
    )
    with pytest.raises(StopIteration):
        next(iterator)


def test_records_invalid_queue():
    with pytest.raises(ValueError):
        next(records("localhost", "IV.MODE", max_queue=0))
