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

An ArrayData holds the stations of an array: records, coordinates
(east, north, up, in metres) and orientation of the components. It is
filled from an array file, from waveforms or one station at a time.

Run (from any folder)::

    python examples/exploration/yarray_01_array.py
"""

import warnings
from pathlib import Path

from shakelab.exploration.yarray import ArrayData
from shakelab.signals.io import reader

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

# From an array file (local coordinates and one file per component).
array = ArrayData(DATA / "sesame_array.json")
array.info()

# From waveforms read by ShakeLab. With geographic positions in the
# metadata (e.g. StationXML) no coordinates are needed.
collection = reader(str(DATA / "sesame" / "*.sac"))
coordinates = {
    stream.sid.rsplit(".", 1)[0]: (
        stream.records()[0].head.extra["sac"]["USER7"],
        stream.records()[0].head.extra["sac"]["USER8"],
    )
    for stream in collection
}
array = ArrayData(collection, coordinates=coordinates,
                  select="*.S100*..*", starttime="2003-01-01T00:01:50")
print(array.sids, array.starttime, array.nsamp)

# One station at a time, from any source.
array = ArrayData()

for station in ["S1001", "S1002", "S1003"]:
    array.add(collection, select=f"*.{station}..*", coordinates=coordinates)

array.remove(".S1003.")
print(array)
print(array.coordinates, array.delta, array.units)
