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
Tests of :mod:`shakelab.signals.seedlink.client`, against the local
server of ``conftest.py``.
"""

from __future__ import annotations

import socket
import warnings
from pathlib import Path

import numpy as np
import pytest

from seedlink_server import (
    HELLO_V3_PLAIN,
    RATE,
    START,
    Station,
    waveform,
)
from shakelab.signals.base import StreamCollection
from shakelab.signals.libio import mseed
from shakelab.signals.metadata import MetadataWarning
from shakelab.signals.seedlink import (
    Client,
    SeedLinkError,
    SeedLinkEvent,
    read,
)
from shakelab.signals.seedlink import client as client_module

DATA = Path(__file__).parents[1] / "data"
VERSIONS = ["3.1", "4.0"]
NSAMP = 6000


@pytest.fixture
def mode():
    """Three components of IV.MODE, 60 s at 100 Hz."""
    records = [
        waveform(f"IV.MODE..HN{component}", NSAMP, seed=index)
        for index, component in enumerate("ZNE")
    ]
    return records, Station("IV", "MODE", records)


@pytest.fixture
def other():
    """One channel of another station."""
    records = [waveform("OX.ABC.00.HHZ", 3000, seed=9)]
    return records, Station("OX", "ABC", records, first_sequence=5)


def _client(server, **kwargs):
    kwargs.setdefault("timeout", 5.0)
    return Client("127.0.0.1", port=server.port, **kwargs)


def _stop_when_complete(client, nsamp, channels):
    """Stop the client when the channels have nsamp samples each."""
    def check(client, record):
        streams = client.buffer
        if len(streams) < channels:
            return
        if all(stream.nsamp >= nsamp for stream in streams):
            client.stop()
    client.on_record(check)


def _assert_same(collection, records):
    """The collection holds the records, one contiguous record each."""
    assert sorted(collection.sid) == sorted(r.sid for r in records)
    for record in records:
        stream = collection[record.sid]
        assert len(stream.record) == 1
        assert stream.record[0].time == record.time
        np.testing.assert_array_equal(stream.record[0].data, record.data)


# ----------------------------------------------------------------------------
# Selection and configuration
# ----------------------------------------------------------------------------


def test_select_patterns():
    client = Client("localhost")
    client.select("IV.MODE")
    client.select(["IV.CAE*.*.HH[ZN]", "OX.ABC..HHZ", "IV.MODE"])
    assert client.selection == [
        "IV.MODE.*.*", "IV.CAE*.*.HH[ZN]", "OX.ABC..HHZ",
    ]
    client.select("OX.{CAE,BAD}.*.HHZ")
    assert client.selection[-2:] == ["OX.CAE.*.HHZ", "OX.BAD.*.HHZ"]


@pytest.mark.parametrize("pattern", ["", "A.B.C.D.E", "IV..*.HHZ",
                                     "IV.MO DE"])
def test_select_invalid(pattern):
    with pytest.raises(ValueError):
        Client("localhost").select(pattern)


def test_configuration_types():
    client = Client("localhost")
    with pytest.raises(TypeError):
        client.select(5)
    with pytest.raises(TypeError):
        client.set_buffer([])
    with pytest.raises(ValueError):
        client.set_buffer(StreamCollection(), max_duration=10)
    with pytest.raises(TypeError):
        client.on_record("not callable")
    with pytest.raises(ValueError):
        client.every(10, print, clock="cpu")
    with pytest.raises(ValueError):
        client.every(0, print)
    with pytest.raises(TypeError):
        Client("localhost", port="18000")
    with pytest.raises(ValueError):
        client.run(duration=1)  # nothing selected


def test_callback_methods_return_function():
    client = Client("localhost")

    @client.on_record
    def handler(client, record):
        pass

    assert handler is not None
    assert client.every(1.0, handler) is handler


# ----------------------------------------------------------------------------
# Data transfer
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
def test_run_unlimited_buffer(make_server, mode, other, version):
    records, station = mode
    server = make_server([station, other[1]], version=version)
    client = _client(server)
    client.select("IV.MODE..HN?")
    client.set_buffer()
    _stop_when_complete(client, NSAMP, 3)

    collection = client.run(duration=10)

    assert client.protocol == version
    assert client.organization == "Test Seismic Network"
    _assert_same(collection, records)
    assert client._sequence["IV_MODE"] == station.packets[-1][0]
    assert server.received("BYE"), server.commands[-1:]


@pytest.mark.parametrize("version", VERSIONS)
def test_fetch_ends_with_server(make_server, mode, version):
    records, station = mode
    server = make_server([station], version=version)
    client = _client(server)
    client.select("IV.MODE")
    events = []
    client.on_event(lambda client, event: events.append(event.kind))

    collection = client.run(fetch=True, duration=10)

    _assert_same(collection, records)
    assert events == ["connected", "end"]


@pytest.mark.parametrize("version", VERSIONS)
def test_past_time_window(make_server, mode, version):
    records, station = mode
    server = make_server([station], version=version)
    start, end = START + 20.0, START + 30.0

    collection = read(
        "127.0.0.1", "IV.MODE..HNZ", start=start, end=end,
        port=server.port, timeout=5.0,
    )

    record = collection["IV.MODE..HNZ"].record[0]
    assert record.time == start
    assert record.nsamp == 1001
    first = int(20 * RATE)
    np.testing.assert_array_equal(
        record.data, records[0].data[first:first + 1001]
    )
    if version == "4.0":
        assert any(command.startswith("DATA ALL 2026-11-05T08:00:20")
                   for command in server.commands)
        assert "ENDFETCH" in server.commands
    else:
        assert any(command.startswith("TIME 2026,11,05,08,00,20")
                   for command in server.commands)


def test_circular_buffer(make_server, mode):
    records, station = mode
    server = make_server([station])
    client = _client(server)
    client.select("IV.MODE..HNZ")
    client.set_buffer(max_duration=10.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        collection = client.run(fetch=True, duration=10)

    stream = collection["IV.MODE..HNZ"]
    assert stream.max_duration == 10.0
    end = stream.record[-1].endtime
    assert end == records[0].endtime
    assert (end - stream.record[0].time) <= 10.0 + 4.0


def test_get_data_joins_circular_buffer(make_server, mode):
    records, station = mode
    server = make_server([station])
    client = _client(server)
    client.select("IV.MODE")
    client.set_buffer(max_duration=120.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        client.run(fetch=True, duration=10)
        data = client.get_data()

    # The buffer keeps the packets; the copy joins them.
    assert all(len(stream.records()) > 1 for stream in client.buffer)
    assert data.max_duration is None
    for stream in data:
        assert type(stream).__name__ == "Stream"
        assert len(stream.records()) == 1
    _assert_same(data, records)
    data["IV.MODE..HNZ"].record[0].data[:] = 0
    assert client.buffer["IV.MODE..HNZ"].record[0].data.any()

    window = client.get_data("*.*.*.*Z", start=START + 1.0, end=START + 2.0)
    assert window["IV.MODE..HNZ"].nsamp == 101


def test_get_data_after_trimming(make_server, mode):
    records, station = mode
    server = make_server([station], interleave=True)
    client = _client(server)
    client.select("IV.MODE")
    client.set_buffer(max_duration=10.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        client.run(fetch=True, duration=10)
    data = client.get_data()

    for record in records:
        stream = data[record.sid]
        assert len(stream.records()) == 1
        joined = stream.record[0]
        assert joined.endtime == record.endtime
        assert 10.0 <= joined.endtime - joined.time <= 14.0
        first = round((joined.time - record.time) * RATE)
        np.testing.assert_array_equal(joined.data, record.data[first:])


@pytest.mark.parametrize("version", VERSIONS)
def test_past_window_of_several_stations(make_server, mode, other, version):
    """Dial-up: the run waits for the server, which sends one station
    (and one channel) after the other."""
    server = make_server([mode[1], other[1]], version=version)
    client = _client(server)
    client.select(["IV.MODE", "OX.ABC"])
    collection = client.run(start=START, end=START + 25.0)

    assert sorted(collection.sid) == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ", "OX.ABC.00.HHZ",
    ]
    for stream in collection:
        assert stream.endtime >= START + 25.0


@pytest.mark.parametrize("interleave", [False, True])
def test_interleaved_data_clock(make_server, mode, interleave):
    """With the packets in time order, the channels of the stations
    are aligned at every tick of the data clock."""
    records = [
        waveform(f"OX.ABC..HN{component}", NSAMP, seed=10 + index)
        for index, component in enumerate("ZNE")
    ]
    stations = [mode[1], Station("OX", "ABC", records)]
    server = make_server(stations, interleave=interleave)
    client = _client(server)
    client.select(["IV.MODE", "OX.ABC"])
    client.set_buffer()
    spreads = []

    def tick(client):
        ends = [stream.endtime for stream in client.buffer]
        spreads.append((len(ends), max(ends) - min(ends)))

    client.every(10.0, tick, clock="data")
    client.run(fetch=True, duration=10)

    assert len(spreads) == 5
    if interleave:
        assert all(count == 6 and spread < 5.0 for count, spread in spreads)
    else:
        # One station after the other: the second one is missing.
        assert all(count == 3 for count, _ in spreads)


@pytest.mark.parametrize("version", VERSIONS)
def test_reconnection_resumes(make_server, mode, version, monkeypatch):
    monkeypatch.setattr(client_module, "RECONNECT_DELAY", (0.05, 0.1))
    records, station = mode
    server = make_server([station], version=version, drop_after=10)
    client = _client(server)
    client.select("IV.MODE")
    events = []
    client.on_event(lambda client, event: events.append(event.kind))

    collection = client.run(fetch=True, duration=10)

    _assert_same(collection, records)
    assert server.connections == 2
    assert events == ["connected", "disconnected", "connected", "end"]
    resumed = station.packets[10][0]
    expected = (
        f"DATA {resumed}" if version == "4.0" else f"FETCH {resumed:06X}"
    )
    assert expected in server.commands


def test_repeated_packets_skipped(make_server, mode):
    records, station = mode
    server = make_server([station])
    client = _client(server)
    client.select("IV.MODE")
    client.run(fetch=True, duration=10)
    # A second run resumes after the last packet: nothing new arrives.
    client._sequence["IV_MODE"] -= 3
    client.set_buffer()
    collection = client.run(fetch=True, duration=10)
    last = [
        mseed.MSRecord.from_bytes(packet[3])
        for packet in station.packets[-3:]
    ]
    assert sorted(collection.sid) == sorted({msr.sid for msr in last})
    assert sum(stream.nsamp for stream in collection) == sum(
        msr.nsamp for msr in last
    )


def test_v3_without_station_wildcards(make_server, mode, other):
    server = make_server([mode[1], other[1]], hello=HELLO_V3_PLAIN)
    client = _client(server)
    client.select("I?.MO*.*.HNZ")
    collection = client.run(fetch=True, duration=10)
    assert collection.sid == ["IV.MODE..HNZ"]
    assert "INFO STATIONS" in server.commands
    assert "STATION MODE IV" in server.commands


def test_server_error(make_server, mode):
    server = make_server([mode[1]], reject="SELECT")
    client = _client(server)
    client.select("IV.MODE..HNZ")
    with pytest.raises(SeedLinkError) as error:
        client.run(duration=5)
    assert error.value.code == "ARGUMENTS"
    assert not client.running


def test_first_connection_failure():
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    client = Client("127.0.0.1", port=port)
    client.select("IV.MODE")
    with pytest.raises(OSError, match=f"127.0.0.1:{port}"):
        client.run(duration=5)


# ----------------------------------------------------------------------------
# Callbacks, background and metadata
# ----------------------------------------------------------------------------


def test_record_and_data_timers(make_server, mode):
    records, station = mode
    server = make_server([station])
    client = _client(server)
    client.select("IV.MODE..HNZ")
    received = []
    ticks = []
    client.on_record(lambda client, record: received.append(record.nsamp))
    client.every(
        20.0,
        lambda client: ticks.append(client.buffer["IV.MODE..HNZ"].nsamp),
        clock="data",
    )
    client.set_buffer()

    client.run(fetch=True, duration=10)

    assert sum(received) == NSAMP
    assert len(ticks) == 2
    assert ticks[0] >= 20 * RATE


def test_wall_timer_and_stop_from_callback(make_server, mode):
    server = make_server([mode[1]])
    client = _client(server)
    client.select("IV.MODE..HNZ")
    ticks = []

    def tick(client):
        ticks.append(1)
        if len(ticks) == 2:
            client.stop()

    client.every(0.2, tick)
    assert client.run(duration=10) is None
    assert len(ticks) == 2


def test_gap_event(make_server):
    records = [
        waveform("IV.MODE..HNZ", 2000),
        waveform("IV.MODE..HNZ", 2000, offset=30.0, seed=1),
    ]
    server = make_server([Station("IV", "MODE", records)])
    client = _client(server)
    client.select("IV.MODE")
    events = []
    client.on_event(lambda client, event: events.append(event))

    collection = client.run(fetch=True, duration=10)

    gaps = [event for event in events if event.kind == "gap"]
    assert len(gaps) == 1
    assert isinstance(gaps[0], SeedLinkEvent)
    assert gaps[0].sid == "IV.MODE..HNZ"
    assert gaps[0].start == START + 20.0
    assert gaps[0].end == START + 30.0
    assert len(collection["IV.MODE..HNZ"].record) == 2


def test_background_and_get_data(make_server, mode):
    records, station = mode
    server = make_server([station])
    client = _client(server)
    client.select("IV.MODE")
    client.set_buffer()

    client.start()
    assert client.running
    with pytest.raises(RuntimeError):
        client.select("IV.CAE")
    for _ in range(100):
        data = client.get_data("IV.MODE..HN?")
        if len(data) == 3 and all(s.nsamp == NSAMP for s in data):
            break
        client.wait(0.05)
    client.stop()

    assert not client.running
    _assert_same(data, records)
    window = client.get_data("*.*.*.*Z", start=START + 1.0, end=START + 2.0)
    assert window.sid == ["IV.MODE..HNZ"]
    assert window["IV.MODE..HNZ"].nsamp == 101


def test_background_error_reraised(make_server, mode):
    server = make_server([mode[1]])
    client = _client(server)
    client.select("IV.MODE")

    def fail(client, record):
        raise ZeroDivisionError("callback failure")

    client.on_record(fail)
    client.start(fetch=True)
    with pytest.raises(ZeroDivisionError):
        client.wait(10)
    assert not client.running


def test_metadata_assigned(make_server, mode, other):
    server = make_server([mode[1], other[1]])
    client = _client(server, metadata=DATA / "MODE.xml")
    client.select(["IV.MODE", "OX.ABC"])

    with pytest.warns(MetadataWarning, match="OX.ABC.00.HHZ"):
        collection = client.run(fetch=True, duration=10)

    record = collection["IV.MODE..HNZ"].record[0]
    assert len(collection["IV.MODE..HNZ"].record) == 1
    assert record.head.meta.sample_rate == 100.0
    assert record.location is not None
    assert collection["OX.ABC.00.HHZ"].record[0].head.meta.empty


def test_read_requires_stop_condition():
    with pytest.raises(ValueError):
        read("localhost", "IV.MODE")
