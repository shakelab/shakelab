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

The monitor is driven by a minimal client (an object with
``get_data``) for the synchronous and threaded processing, the queue,
the callbacks and the queries; one integration test runs it within the
SeedLink client on the local test server of
``tests/signals/seedlink/seedlink_server.py``, with a data clock.
"""

from __future__ import annotations

import importlib.util
import math
import threading
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.beamforming import Beamformer
from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.io import read_results
from shakelab.exploration.yarray.monitor import ArrayMonitor
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection

START = Date("2026-11-05T08:00:00")
DELTA = 0.01
VELOCITY = 400.0
AZIMUTH = 60.0
COORDINATES = np.random.default_rng(1).uniform(-100.0, 100.0, (6, 2))
SIDS = [f"XX.S{i:02d}." for i in range(6)]
CONFIG = BeamformingConfig(
    [4.0],
    SearchGrid(200.0, 800.0, 31, scale="linear", azimuth_size=72),
    window=WindowConfig(type="fixed", length=2.0, windows_per_block=20),
    method="bartlett",
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

    def get_data(self, sid="*", start=None, end=None):
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
    return ArrayMonitor(Beamformer(CONFIG, cache=True), **options)


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
    path = tmp_path / "monitor.pkl"
    monitor = _monitor(output=path, keep=2)
    seen = []
    monitor.on_result(lambda monitor, result: seen.append(result))
    client = _Client(_records())

    for _ in range(3):
        monitor(client)

    assert len(seen) == 3
    assert len(monitor.results()) == 2
    assert len(read_results(path)) == 3
    assert monitor.results(endtime=START + 10.0) == []

    joined = monitor.result()
    assert len(joined.block_times[4.0]) == 2 * len(
        seen[0].block_times[4.0]
    )


def test_skipped_windows():
    """Empty buffers and too few stations are skipped and reported."""
    reasons = []
    monitor = _monitor(min_stations=7)
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))

    monitor(_Client([]))
    monitor(_Client(_records()))

    assert monitor.status()["skipped"] == 2
    assert "no data" in reasons[0]
    assert "at least 7" in reasons[1]
    assert monitor.latest() is None


def test_too_short_window_is_skipped():
    """A window without stations covering it is skipped, not raised."""
    monitor = _monitor(window=100.0)
    reasons = []
    monitor.on_skip(lambda monitor, reason: reasons.append(reason))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        monitor(_Client(_records()))

    assert monitor.status()["skipped"] == 1
    assert "no station covers" in reasons[0]


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

    def process(self, array):
        self.release.wait(5.0)
        return super().process(array)


def test_late_windows_are_dropped():
    """
    With the worker busy, a new call replaces the pending window.

    The first window is being processed, the second is queued; the
    third replaces it and the second is reported as dropped.
    """
    beamformer = _SlowBeamformer(CONFIG)
    monitor = ArrayMonitor(
        beamformer, window=30.0, coordinates=dict(zip(SIDS, COORDINATES))
    )
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

    def process(self, array):
        raise KeyError("unexpected")


def test_unexpected_error_stops_the_worker():
    """Unexpected errors stop the worker and are raised by close."""
    monitor = ArrayMonitor(
        _FailingBeamformer(CONFIG), window=30.0,
        coordinates=dict(zip(SIDS, COORDINATES)),
    )
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
    """Wrong beamformer and windows are rejected."""
    with pytest.raises(TypeError):
        ArrayMonitor("capon", window=10.0)
    with pytest.raises(ValueError):
        ArrayMonitor(CONFIG, window=0.0)
    with pytest.raises(ValueError):
        ArrayMonitor(CONFIG, window=10.0, keep=0)

    assert isinstance(ArrayMonitor(CONFIG, window=10.0).beamformer,
                      Beamformer)


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
        monitor = ArrayMonitor(
            Beamformer(CONFIG, cache=True), window=30.0,
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
    assert threads[0].name == "ArrayMonitor"
    best = monitor.latest().maxima_at(rank=1)[0]
    assert best.velocity == pytest.approx(VELOCITY, abs=20.0)
    assert best.azimuth == pytest.approx(AZIMUTH)
