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
ArrayMonitor every INTERVAL seconds; the monitor processes the last
WINDOW seconds in its own thread, so that the reception is never
delayed, and can be queried meanwhile.

Set SERVER, SELECT and METADATA to an array streamed by a SeedLink
server, then run (from any folder; results go to ./output)::

    python examples/exploration/yarray_12_realtime.py
"""

import time

from shakelab.exploration.yarray import ArrayMonitor, BeamformingConfig
from shakelab.libutils.timeN import Date
from shakelab.signals.io import metadata_reader
from shakelab.signals.seedlink import Client

SERVER = "rtserve.example.org"      # SeedLink server of the array
SELECT = "XX.A*..HH?"               # channels of the array stations
METADATA = "array.xml"              # positions and orientations
WINDOW = 300.0                      # seconds analysed at each call
INTERVAL = 60.0                     # seconds between two analyses

config = BeamformingConfig(
    fmin=2.0, fmax=10.0, nfreq=15,
    grid_min=100.0, grid_max=2000.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=5.0, transform="rfft",
)

client = Client(SERVER, metadata=metadata_reader(METADATA))
client.select(SELECT)
client.set_buffer(max_duration=2 * WINDOW)

monitor = ArrayMonitor(config, window=WINDOW, output="output/monitor.pkl")


@monitor.on_result
def show(monitor, result):
    best = result.maxima_at(rank=1, frequency=5.0)[0]
    print(f"{best.time}  5 Hz: {best.velocity:.0f} m/s, "
          f"azimuth {best.azimuth:.0f}")


# clock="data": intervals on the time of the data (also for old data).
client.every(INTERVAL, monitor, clock="data")
client.start(start=Date.now() - 1800.0, duration=900.0)

while client.running:
    time.sleep(INTERVAL)
    print(monitor.status())

client.stop()
monitor.close()
monitor.result().dispersion().info()
