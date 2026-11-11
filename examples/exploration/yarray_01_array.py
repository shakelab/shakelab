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
Building the array data.

An ArrayData holds the records of the stations on a common time axis,
with their coordinates (east, north, up, in metres). This example
builds the same array of the SESAME simulation (nine three-component
stations, 60 s) in three ways.

Run (from any folder)::

    python examples/exploration/yarray_01_array.py
"""

import warnings
from pathlib import Path

from shakelab.exploration import yarray as ya
from shakelab.signals.io import reader

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

# The SESAME files keep local coordinates in STLA/STLO: the SAC reader
# warns that they are not geographic.
warnings.filterwarnings("ignore", message="STLO")

# 1. From an array file: a JSON file listing, for each station, its
#    local coordinates and one waveform file per component.
array = ya.read_array(DATA / "sesame_array_9.json")
array.info()

# 2. From waveforms read by ShakeLab and a table of coordinates.
#    The keys are NET.STA.LOC of the records (here ".S1001." and so
#    on: the SESAME files have no network code). With geographic
#    positions in the metadata (e.g. StationXML) the table is not
#    needed: the coordinates are computed in a local ENU frame.
collection = reader(str(DATA / "sesame" / "*.sac"))

coordinates = {}
for stream in collection:
    header = stream.records()[0].head.extra["sac"]
    coordinates[stream.sid.rsplit(".", 1)[0]] = (header["USER7"],
                                                 header["USER8"])

array = ya.ArrayData.from_collection(collection, coordinates=coordinates)
print(array.sids)

#    The channels and the time window can be selected.
array = ya.ArrayData.from_collection(
    collection,
    coordinates=coordinates,
    select="*.S10[01]*..*",
    starttime="2003-01-01T00:01:50",
    endtime="2003-01-01T00:02:30",
)
print(array.sids, array.starttime, array.nsamp)

# 3. From records already in memory (e.g. processed by the user).
records = [record for stream in collection for record in stream.records()]
array = ya.ArrayData.from_records(records, coordinates=coordinates)

# The array is immutable: later changes of the records do not
# affect it. Its main properties:
print(array.delta, array.units)          # sampling, units (None here)
print(array.coordinates[:3])             # (n_stations, 3) in metres
station = array.stations[0]
print(station.sid, station.ncomp, station.vertical()[:5])
