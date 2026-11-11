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
Real-time array monitoring with the SeedLink client.

An :class:`ArrayMonitor` is a function of the SeedLink client: called
at regular intervals (``client.every``) in the receiving thread, it
copies the latest window of the buffer and hands it to its worker
thread, so that the reception is never delayed by the processing. The
monitor can be queried at any time (``status``, ``latest``,
``results``, ``result``, ``array``) while the client runs.

Without arguments the example replays the nine SESAME stations (60 s)
through a small stand-in of the client, which delivers the data in
packets in a receiving thread and calls the monitor every 10 s of data,
while the main thread queries the monitor. With a server it runs the
same monitor on live data::

    python examples/exploration/yarray_realtime.py \\
        --server rtserve.example.org --select "XX.A*..HH?" \\
        --metadata array.xml --duration 600

(positions and orientations from the metadata; a SeedLink server of
an array is needed). Run from the root of the repository::

    python examples/exploration/yarray_realtime.py
"""

from __future__ import annotations

import argparse
import threading
import time
import warnings

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayMonitor,
    Beamformer,
    BeamformingConfig,
    FrequencyConfig,
    SearchGrid,
    WindowConfig,
    read_array,
)
from shakelab.exploration.yarray.statistics import maximum_density
from shakelab.graphics.yarray import plot_density, plot_dispersion
from shakelab.signals.base import StreamCollection

ARRAY = "tests/exploration/yarray/data/sesame_array_9.json"

warnings.filterwarnings("ignore", message="STLO")

# Windows of 20 s analysed every 10 s: 19 windows of 2 s (overlap 0.5)
# per analysis, one covariance block of 18.
CONFIG = BeamformingConfig(
    frequencies=FrequencyConfig(minimum=3.0, maximum=10.0, size=8),
    grid=SearchGrid(100.0, 1500.0, 50, azimuth_size=72),
    window=WindowConfig(type="fixed", length=2.0, windows_per_block=18),
    method="capon",
    transform="rfft",
)
WINDOW = 20.0
INTERVAL = 10.0


class ReplayClient:
    """
    Stand-in of :class:`shakelab.signals.seedlink.Client` for files.

    Delivers the records of an array in packets of ``packet`` seconds
    of data, in a receiving thread, and calls the functions registered
    with :meth:`every` when the data time passes a multiple of their
    interval (as ``clock="data"`` of the client). ``get_data`` returns
    the data received so far, like the buffer of the client.
    """

    def __init__(self, records, packet=1.0, speed=50.0):
        self._records = list(records)
        self._packet = packet
        self._speed = speed
        self._timers = []
        self._now = None
        self._lock = threading.Lock()
        self._thread = None

    def every(self, interval, func, clock="data"):
        """Call ``func(client)`` every ``interval`` s of data."""
        self._timers.append([interval, func, None])
        return func

    def get_data(self, select="*"):
        """Copy of the data received so far."""
        collection = StreamCollection()

        with self._lock:
            now = self._now

        if now is None:
            return collection

        for record in self._records:
            if now > record.starttime:
                collection.append(record.extract(None, now))

        return collection

    def start(self):
        """Start receiving in a background thread."""
        self._thread = threading.Thread(target=self._receive,
                                        name="ReplayClient")
        self._thread.start()

    def join(self):
        """Wait for the end of the data."""
        self._thread.join()

    def _receive(self):
        """Advance the data time packet by packet."""
        start = min(record.starttime for record in self._records)
        end = min(record.endtime for record in self._records)
        now = start

        while now < end:
            now = min(now + self._packet, end)

            with self._lock:
                self._now = now

            for timer in self._timers:
                interval, func, last = timer
                step = int((now - start) // interval)

                if step > 0 and step != last:
                    timer[2] = step
                    func(self)

            time.sleep(self._packet / self._speed)


def replay(array):
    """Monitor the replayed SESAME records; return the monitor."""
    records = [record for station in array for record in station.records]
    client = ReplayClient(records)

    # Positions by NET.STA.LOC of the records (".S1001." in the files).
    monitor = ArrayMonitor(
        Beamformer(CONFIG, cache=True), window=WINDOW,
        coordinates={
            record.sid.rsplit(".", 1)[0]: station.coordinates
            for station in array for record in station.records
        },
    )

    @monitor.on_result
    def report(monitor, result):
        """Called in the worker thread after each analysis."""
        best = result.maxima_at(rank=1, frequency=5.0)
        print(f"  analysis of {monitor.status()['last_window'][1]}: "
              f"5 Hz {best[0].velocity:.0f} m/s" if best else
              "  analysis without maxima at 5 Hz")

    @monitor.on_skip
    def skipped(monitor, reason):
        """Called when a window cannot be analysed."""
        print(f"  skipped: {reason}")

    client.every(INTERVAL, monitor, clock="data")
    client.start()

    # The main thread is free: query the monitor while data arrive.
    while client._thread.is_alive():
        status = monitor.status()
        print(f"calls {status['calls']}, analyses {status['analyses']}, "
              f"pending {status['pending']}")
        time.sleep(0.4)

    client.join()
    monitor.close()          # waits for the last analysis
    return monitor


def live(arguments):
    """Monitor a SeedLink server (not run by the test suite)."""
    from shakelab.signals import seedlink
    from shakelab.signals.io import metadata_reader

    client = seedlink.Client(
        arguments.server,
        metadata=metadata_reader(arguments.metadata),
    )
    client.select(arguments.select)
    client.set_buffer(max_duration=3 * WINDOW)
    monitor = ArrayMonitor(Beamformer(CONFIG, cache=True), window=WINDOW,
                           output="monitor.pkl")
    monitor.on_result(lambda monitor, result: print(result.summary()))
    client.on_event(monitor.on_event)
    client.every(INTERVAL, monitor)
    client.start()

    try:
        for _ in range(int(arguments.duration // INTERVAL)):
            time.sleep(INTERVAL)
            print(monitor.status())
    finally:
        client.stop()
        monitor.close()

    return monitor


parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
parser.add_argument("--server", help="SeedLink server (host[:port])")
parser.add_argument("--select", default="*", help="channel pattern")
parser.add_argument("--metadata", help="metadata of the array")
parser.add_argument("--duration", type=float, default=600.0,
                    help="seconds of monitoring")
arguments, _ = parser.parse_known_args()

if arguments.server:
    monitor = live(arguments)
else:
    monitor = replay(read_array(ARRAY))

print(monitor.status())
joined = monitor.result()            # all the windows, joined
joined.info()

if joined.maxima:
    density = maximum_density(joined, coordinate="log_velocity")
    figure, axes = plt.subplots(1, 2, figsize=(14, 5))
    plot_dispersion(joined, ax=axes[0])
    plot_density(density, peaks=density.peaks(2), ax=axes[1])
    figure.tight_layout()
    plt.show()
