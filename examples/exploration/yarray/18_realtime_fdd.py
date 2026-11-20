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
Real-time frequency domain decomposition with the SeedLink client.

An FDDMonitor decomposes, at regular intervals, the last window of
data of the client buffer. Its results join into one FDDResult, since
the monitor always processes the same stations: those of its first
result, or those given in advance.
"""

from shakelab.exploration.yarray import FDDConfig, FDDMonitor
from shakelab.signals.seedlink import Client

config = FDDConfig(
    fmin=0.5, fmax=20.0,
    window_type="fixed", window_length=20.0,    # seconds
    windows_per_block=10,
    components="three_component",               # or "vertical", ...
)

client = Client("rtserve.example.org")          # no metadata needed
client.select("XX.B*..HN?")
client.set_buffer(max_duration=2400)

monitor = FDDMonitor(
    config, window=1800,
    output="fdd.jsonl",                         # arrays=False: no .npz
    log="monitor.jsonl",                        # analyses and skips
    stations=["XX.B01.", "XX.B02.", "XX.B03."],  # default: the first result
    # gaps="drop" to exclude stations with gaps, default "raise"
)


@monitor.on_skip
def skipped(monitor, reason):           # e.g. a station is missing
    print("skipped:", reason)


client.every(300, monitor)
client.on_event(monitor.on_event)
client.start()

# From any thread, while it runs:
monitor.status()
result = monitor.result()               # all the results, joined
result.info()
result.mode(result.peaks()[0]).info()

client.stop()
monitor.close()
