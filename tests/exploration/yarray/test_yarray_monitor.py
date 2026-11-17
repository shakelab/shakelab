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
Unit tests for :mod:`shakelab.exploration.yarray.monitor`.

The monitors are driven by a minimal client (an object with
``get_data``) for the synchronous and threaded processing, the queue,
the callbacks and the queries; one integration test runs the
beamforming monitor within the SeedLink client on the local test
server of ``tests/signals/seedlink/seedlink_server.py``, with a data
clock. The FDD monitor is tested on the same synthetic records, for
the stations it keeps, the gaps and the joined results.
"""

from __future__ import annotations

import importlib.util
import json
import math
import threading
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.beamforming import Beamformer
from shakelab.exploration.yarray.config import BeamformingConfig, FDDConfig
from shakelab.exploration.yarray.fdd import FDDResult
from shakelab.exploration.yarray.io import read_result
from shakelab.exploration.yarray.monitor import (
    ArrayMonitor,
    BeamformingMonitor,
    FDDMonitor,
)
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection

START = Date("2026-11-05T08:00:00")
DELTA = 0.01
VELOCITY = 400.0
AZIMUTH = 60.0
COORDINATES = np.random.default_rng(1).uniform(-100.0, 100.0, (6, 2))
SIDS = [f"XX.S{i:02d}." for i in range(6)]
CONFIG = BeamformingConfig(
    frequencies=[4.0], grid_min=200.0, grid_max=800.0, grid_size=31,
    grid_scale="linear", azimuth_size=72, window_type="fixed",
    window_length=2.0, windows_per_block=20, method="bartlett",
)


def _records(duration=60.0, scale=1.0, integer=False):
    """Vertical records of a 4 Hz plane wave at the stations."""
    angle = math.radians(AZIMUTH)
    slowness = np.array([math.sin(angle), math.cos(angle)]) / VELOCITY
    time_axis = np.arange(int(duration / DELTA)) * DELTA
    rng = np.random.default_rng(3)
    records = []

    for sid, xy in zip(SIDS, COORDINATES):
        data = scale * (
            np.cos(2 * np.pi * 4.0 * (time_axis - xy @ slowness))
            + 0.2 * rng.standard_normal(time_axis.size)
        )

        if integer:
            data = np.round(data).astype(np.int32)

        records.append(
            Record(data=data, delta=DELTA, time=START, sid=sid + ".HHZ")
        )

    return records


class _Client:
    """Minimal stand-in of the SeedLink client: a fixed buffer."""

    def __init__(self, records):
        self.buffer = StreamCollection()

        for record in records:
            self.buffer.append(record)

        self.calls = []

    def get_data(self, sid="*", *, starttime=None, endtime=None):
        self.calls.append(sid)
        return self.buffer.copy()


def _monitor(**kwargs):
    """Monitor with the coordinates of the synthetic stations."""
    options = {
        "window": 30.0,
        "coordinates": dict(zip(SIDS, COORDINATES)),
        "background": False,
    }
    options.update(kwargs)
    return BeamformingMonitor(CONFIG, **options)


def test_synchronous_analysis():
    """A call processes the last window of the buffer."""
    monitor = _monitor()
    client = _Client(_records())
    monitor(client)

    status = monitor.status()
    result = monitor.latest()
    best = result.maxima_at(rank=1)[0]

    assert status["analyses"] == 1 and status["calls"] == 1
    assert status["stations"] == 6
    assert status["last_window"][1] == START + 59.99
    assert status["last_window"][0] == START + 29.99
    assert best.velocity == pytest.approx(VELOCITY, abs=20.0)
    assert best.azimuth == pytest.approx(AZIMUTH)
    assert monitor.array().nsta == 6
    assert client.calls == ["*"]


def test_callbacks_and_queries(tmp_path):
    """on_result, results in an interval, joined result, output."""
    path = tmp_path / "monitor.jsonl"
    log = tmp_path / "log" / "monitor.jsonl"
    monitor = _monitor(output=path, keep=2, log=log)
    seen = []
    monitor.on_result(lambda monitor, result: seen.append(result))
    client = _Client(_records())

    for _ in range(3):
        monitor(client)

    assert len(seen) == 3
    assert len(monitor.results()) == 2
    assert len(read_result(path)[1]) == 3
    events = [json.loads(line) for line in log.read_text().splitlines()]
    assert [event["event"] for event in events] == ["result"] * 3
    assert events[0]["stations"] == 6
    assert Date(events[0]["window"][1]) == START + 59.99
    assert monitor.results(endtime=START + 10.0) == []

    joined = monitor.result()
    assert len(joined.block_times[4.0]) == 2 * len(
        seen[0].block_times[4.0]
    )

    latest = monitor.latest()
    latest.maxima.clear()
    assert monitor.latest().maxima


def test_skipped_windows(tmp_path):
    """Empty buffers and too few stations are skipped and reported."""
    reasons = []
    log = tmp_path / "monitor.jsonl"
    monitor = _monitor(min_stations=7, log=log)
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))

    monitor(_Client([]))
    monitor(_Client(_records()))

    assert monitor.status()["skipped"] == 2
    assert "no data" in reasons[0]
    assert "at least 7" in reasons[1]
    assert monitor.latest() is None

    events = [json.loads(line) for line in log.read_text().splitlines()]
    assert [event["event"] for event in events] == ["skip", "skip"]
    assert events[0]["window"] is None
    assert events[1]["reason"] == reasons[1]
    assert Date(events[1]["window"][0]) == START + 29.99


def test_too_short_buffer_is_skipped():
    """
    A buffer shorter than the window is skipped, quietly.

    At the start of the acquisition the buffer holds less than one
    window: no array is built (no warnings about stations) and the
    reason says how much data there is.
    """
    monitor = _monitor(window=100.0)
    reasons = []
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        monitor(_Client(_records()))

    assert monitor.status()["skipped"] == 1
    assert reasons[0] == "the buffer covers 59.99 s, the window is 100 s"


def test_window_without_blocks_is_skipped():
    """A window shorter than one covariance block is skipped."""
    monitor = _monitor(window=10.0)
    reasons = []
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        monitor(_Client(_records()))

    assert monitor.latest() is None
    assert "covariance block" in reasons[0]


def test_background_processing():
    """The worker processes the queued windows; close waits."""
    monitor = _monitor(background=True)
    client = _Client(_records())
    monitor(client)
    monitor.close()

    assert monitor.status()["analyses"] == 1
    assert not monitor.running


class _SlowBeamformer(Beamformer):
    """Beamformer that waits for an event before processing."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.release = threading.Event()

    def process(self, *args, **kwargs):
        self.release.wait(5.0)
        return super().process(*args, **kwargs)


def test_late_windows_are_dropped():
    """
    With the worker busy, a new call replaces the pending window.

    The first window is being processed, the second is queued; the
    third replaces it and the second is reported as dropped.
    """
    beamformer = _SlowBeamformer(config=CONFIG, cache=True)
    monitor = BeamformingMonitor(
        CONFIG, window=30.0, coordinates=dict(zip(SIDS, COORDINATES))
    )
    monitor._processor = beamformer
    reasons = []
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))
    client = _Client(_records())

    monitor(client)
    deadline = time.monotonic() + 5.0

    while monitor.status()["pending"] and time.monotonic() < deadline:
        time.sleep(0.01)

    monitor(client)
    monitor(client)
    beamformer.release.set()
    monitor.close()

    assert monitor.status()["analyses"] == 2
    assert monitor.status()["skipped"] == 1
    assert "slower" in reasons[0]


class _FailingBeamformer(Beamformer):
    """Beamformer with an unexpected error."""

    def process(self, *args, **kwargs):
        raise KeyError("unexpected")


def test_unexpected_error_stops_the_worker():
    """Unexpected errors stop the worker and are raised by close."""
    monitor = BeamformingMonitor(
        CONFIG, window=30.0, coordinates=dict(zip(SIDS, COORDINATES)),
    )
    monitor._processor = _FailingBeamformer(config=CONFIG)
    monitor(_Client(_records()))
    deadline = time.monotonic() + 5.0

    while monitor.running and time.monotonic() < deadline:
        time.sleep(0.01)

    assert "unexpected" in monitor.status()["error"]

    with pytest.raises(RuntimeError, match="stopped"):
        monitor.close()
    with pytest.raises(RuntimeError):
        monitor(_Client(_records()))


def test_invalid_arguments():
    """Wrong configurations and windows are rejected."""
    with pytest.raises(TypeError):
        BeamformingMonitor(42, window=10.0)
    with pytest.raises(ValueError, match="frequencies"):
        BeamformingMonitor(BeamformingConfig(), window=10.0)
    with pytest.raises(ValueError):
        BeamformingMonitor(CONFIG, window=0.0)
    with pytest.raises(ValueError):
        BeamformingMonitor(CONFIG, window=10.0, keep=0)
    with pytest.raises(ValueError, match="jsonl"):
        BeamformingMonitor(CONFIG, window=10.0, output="result.json")
    with pytest.raises(ValueError, match="jsonl"):
        FDDMonitor(FDD_CONFIG, window=10.0, log="log.json")

    monitor = BeamformingMonitor(CONFIG, window=10.0)
    assert monitor.config == CONFIG and monitor.config is not CONFIG


def _seedlink_test_server():
    """The module of the local SeedLink test server."""
    path = (
        Path(__file__).parents[2] / "signals" / "seedlink"
        / "seedlink_server.py"
    )
    spec = importlib.util.spec_from_file_location("_seedlink_server", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_within_the_seedlink_client():
    """
    The monitor runs as a function of the SeedLink client.

    Six stations stream 60 s of a plane wave from the local test
    server, which sends the stations one after the other: when the data
    clock reaches 30 s only the first station has arrived, and the
    monitor, called in the receiving thread, skips the window (too few
    stations). Called again on the complete buffer, it finds the wave.
    """
    from shakelab.signals.seedlink import Client

    server_module = _seedlink_test_server()
    records = _records(scale=1000.0, integer=True)
    stations = [
        server_module.Station("XX", sid.split(".")[1], [record])
        for sid, record in zip(SIDS, records)
    ]
    server = server_module.FakeServer(stations)
    threads = []

    try:
        client = Client("127.0.0.1", port=server.port, timeout=5.0)
        client.select("XX.*..HHZ")
        client.set_buffer()
        monitor = BeamformingMonitor(
            CONFIG, window=30.0,
            coordinates=dict(zip(SIDS, COORDINATES)),
        )
        monitor.on_skip(
            lambda monitor, reason: threads.append(threading.current_thread())
        )
        client.every(30.0, monitor, clock="data")
        client.run(fetch=True, duration=20)

        assert monitor.status()["calls"] == 1

        monitor(client)
        monitor.close()

    finally:
        server.close()

    status = monitor.status()
    assert status["skipped"] == 1 and status["analyses"] == 1
    assert threads[0].name == "BeamformingMonitor"
    best = monitor.latest().maxima_at(rank=1)[0]
    assert best.velocity == pytest.approx(VELOCITY, abs=20.0)
    assert best.azimuth == pytest.approx(AZIMUTH)


def test_circular_buffer_of_the_seedlink_client():
    """
    The monitor processes the circular buffer of the client.

    The test server sends the packets of the six stations in time order,
    as a real-time server; the buffer keeps the last 40 s, as separate
    packets. Every 15 s of data the monitor processes the last 30 s:
    the first call is skipped (15 s of data), the others are analysed.
    """
    from shakelab.signals.seedlink import Client

    server_module = _seedlink_test_server()
    records = _records(scale=1000.0, integer=True)
    stations = [
        server_module.Station("XX", sid.split(".")[1], [record])
        for sid, record in zip(SIDS, records)
    ]
    server = server_module.FakeServer(stations, interleave=True)
    reasons = []

    try:
        client = Client("127.0.0.1", port=server.port, timeout=5.0)
        client.select("XX.*..HHZ")
        client.set_buffer(max_duration=40.0)
        monitor = BeamformingMonitor(
            CONFIG, window=30.0, background=False,
            coordinates=dict(zip(SIDS, COORDINATES)),
        )
        monitor.on_skip(lambda monitor, reason: reasons.append(reason))
        client.every(15.0, monitor, clock="data")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            client.run(fetch=True, duration=20)

        assert all(
            len(stream.records()) > 1 for stream in client.buffer
        )

    finally:
        server.close()

    status = monitor.status()
    assert status["calls"] == 3
    assert status["analyses"] == 2
    assert len(reasons) == 1 and "the buffer covers" in reasons[0]
    for result in monitor.results():
        best = result.maxima_at(rank=1)[0]
        assert best.velocity == pytest.approx(VELOCITY, abs=20.0)
        assert best.azimuth == pytest.approx(AZIMUTH)


# ----------------------------------------------------------------------------
# Base class and FDD monitor
# ----------------------------------------------------------------------------

FDD_CONFIG = FDDConfig(
    fmin=2.0, fmax=8.0, window_type="fixed", window_length=2.0,
    windows_per_block=5, components="vertical",
)


def _noise(sid, duration=60.0, start=START, seed=0):
    """Vertical record of white noise."""
    data = np.random.default_rng(seed).standard_normal(
        int(duration / DELTA)
    )
    return Record(data=data, delta=DELTA, time=start, sid=sid + ".HHZ")


def test_base_monitor_is_not_used_directly():
    """ArrayMonitor is only the base of the monitors."""
    with pytest.raises(TypeError, match="BeamformingMonitor or FDDMonitor"):
        ArrayMonitor(CONFIG, window=30.0)

    with pytest.raises(TypeError):
        FDDMonitor(CONFIG, window=30.0)
    with pytest.raises(TypeError):
        BeamformingMonitor(FDD_CONFIG, window=30.0)


def test_fdd_monitor_keeps_its_stations(tmp_path):
    """
    The first result fixes the stations of the FDD monitor.

    A station added to the buffer later is ignored; a window without
    one of the first stations is skipped. The results kept and the
    output file join into one FDDResult.
    """
    output = tmp_path / "fdd.jsonl"
    monitor = FDDMonitor(FDD_CONFIG, window=30.0, background=False,
                         output=output, arrays=False)
    reasons = []
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))
    records = _records()
    client = _Client(records)

    assert monitor.stations is None

    monitor(client)

    assert monitor.stations == tuple(SIDS)
    assert monitor.latest().channels == tuple(
        f"{sid} vertical" for sid in SIDS
    )

    client.buffer.append(_noise("XX.S99."))
    monitor(client)

    assert monitor.latest().stations == monitor.results()[0].stations

    monitor(_Client(records[1:]))

    assert monitor.status()["analyses"] == 2
    assert "XX.S00. missing" in reasons[0]

    joined = monitor.result()
    first = monitor.results()[0]

    assert isinstance(joined, FDDResult)
    assert joined.nblocks == 2 * first.nblocks
    assert FDDResult(output).nblocks == joined.nblocks
    assert not (tmp_path / "fdd_arrays").exists()
    assert np.isnan(joined.coordinates).all()


def test_fdd_monitor_given_stations_and_gaps():
    """
    Stations given in advance, in their order; gaps skip the window.

    With the default ``gaps="raise"`` a selected station that does not
    cover the window skips it; stations not processed do not matter.
    """
    order = [SIDS[2], SIDS[0], SIDS[1]]
    records = _records()
    late = _noise(SIDS[3], duration=45.0, start=Date(START) + 15.0)
    client = _Client(records[:3] + [late])
    reasons = []
    monitor = FDDMonitor(FDD_CONFIG, window=50.0, stations=order,
                         background=False)
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))
    monitor(client)

    assert monitor.latest().stations == tuple(order)

    everything = FDDMonitor(FDD_CONFIG, window=50.0, background=False)
    everything.on_skip(lambda monitor, reason: reasons.append(reason))
    everything(client)

    assert everything.latest() is None
    assert "XX.S03." in reasons[-1]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dropped = FDDMonitor(FDD_CONFIG, window=50.0, gaps="drop",
                             background=False)
        dropped(client)

    assert dropped.stations == tuple(SIDS[:3])

    for stations in ("", ["XX.S00.", "XX.S00."]):
        with pytest.raises(ValueError):
            FDDMonitor(FDD_CONFIG, window=30.0, stations=stations)
    with pytest.raises(TypeError):
        FDDMonitor(FDD_CONFIG, window=30.0, stations=[1, 2])


def test_beamforming_monitor_needs_positions():
    """Without positions the beamforming windows are skipped."""
    reasons = []
    monitor = BeamformingMonitor(CONFIG, window=30.0, background=False)
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))
    monitor(_Client(_records()))

    assert monitor.status()["skipped"] == 1
    assert "needs the station positions" in reasons[0]
