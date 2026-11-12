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

The client receives the data in a background thread and calls an
ArrayMonitor at regular intervals; the monitor processes the last
seconds of data in its own thread, so that the reception is never
delayed, and can be queried meanwhile.
"""

from shakelab.exploration.yarray import ArrayMonitor, BeamformingConfig
from shakelab.signals.io import metadata_reader
from shakelab.signals.seedlink import Client

config = BeamformingConfig(
    fmin=2.0, fmax=10.0, nfreq=15,
    grid_min=100.0, grid_max=2000.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=5.0, transform="rfft",
)

client = Client("rtserve.example.org", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=600)

monitor = ArrayMonitor(config, window=300, output="monitor.pkl")


@monitor.on_result
def show(monitor, result):
    print(result.maxima_at(rank=1, frequency=5.0))


@monitor.on_skip
def skipped(monitor, reason):           # windows that cannot be processed
    print("skipped:", reason)


client.every(60, monitor)               # clock="data" to replay an archive
client.on_event(monitor.on_event)       # connection events in the log
client.start()

# From any thread, while it runs:
monitor.status()
monitor.latest()
monitor.result("2026-11-01T10:00", "2026-11-01T11:00")

client.stop()
monitor.close()
