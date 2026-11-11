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

A SeedLink client receives the data of an array in a background
thread. Every ``INTERVAL`` seconds it calls an ArrayMonitor, which
copies the last ``WINDOW`` seconds of the buffer and processes them
in its own worker thread, so that the reception is never delayed.
Meanwhile the main program is free and can query the monitor.

Set SERVER, SELECT and METADATA to an array streamed by a SeedLink
server, then run (from any folder; results go to ./output)::

    python examples/exploration/yarray_12_realtime.py
"""

import time

from shakelab.exploration import yarray as ya
from shakelab.libutils.timeN import Date
from shakelab.signals import seedlink
from shakelab.signals.io import metadata_reader

SERVER = "rtserve.example.org"      # SeedLink server of the array
SELECT = "XX.A*..HH?"               # channels of the array stations
METADATA = "array.xml"              # positions and orientations
WINDOW = 300.0                      # seconds analysed at each call
INTERVAL = 60.0                     # seconds between two analyses

# 1. Processing applied to each window: vertical Capon beamforming.
config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=2.0, maximum=10.0, size=15),
    grid=ya.SearchGrid(100.0, 2000.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=5.0,
                           windows_per_block=40),
    method="capon",
    transform="rfft",
)

# 2. The SeedLink client keeps the last two windows in its buffer.
client = seedlink.Client(SERVER, metadata=metadata_reader(METADATA))
client.select(SELECT)
client.set_buffer(max_duration=2 * WINDOW)

# 3. The monitor: one analysis at each call; each result is also
#    appended to a file. cache=True keeps the steering vectors between
#    analyses (same stations, same grid).
monitor = ya.ArrayMonitor(
    ya.Beamformer(config, cache=True),
    window=WINDOW,
    output="output/monitor.pkl",
)


# 4. Optional callbacks, run in the worker thread of the monitor.
def on_result(monitor, result):
    best = result.maxima_at(rank=1, frequency=5.0)[0]
    print(f"{best.time}  5 Hz: {best.velocity:.0f} m/s, "
          f"azimuth {best.azimuth:.0f}")


def on_skip(monitor, reason):
    print(f"window skipped: {reason}")


monitor.on_result(on_result)
monitor.on_skip(on_skip)

# 5. The client calls the monitor at regular intervals. With
#    clock="data" the interval is counted on the time of the data, which
#    also works when old data arrive faster than real time.
client.every(INTERVAL, monitor, clock="data")

# 6. Receive in a background thread: the last 30 minutes still held by
#    the server, then the live data, for 15 minutes in all.
client.start(start=Date.now() - 1800.0, duration=900.0)

# 7. The main program is free: query the monitor while it runs.
while client.running:
    time.sleep(INTERVAL)
    status = monitor.status()
    print(f"analyses {status['analyses']}, skipped {status['skipped']}, "
          f"stations {status['stations']}")

# 8. Stop: the monitor waits for the last analysis.
client.stop()
monitor.close()

# 9. All the analyses joined into one result.
result = monitor.result()
result.info()
