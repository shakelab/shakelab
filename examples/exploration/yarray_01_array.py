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
"""

from shakelab.exploration.yarray import ArrayData
from shakelab.signals.io import metadata_reader, reader

# From an array file (local coordinates and one file per component).
array = ArrayData("array.json")
array.info()

# From waveforms with metadata (geographic positions and orientations).
collection = reader("data/*.mseed", metadata=metadata_reader("array.xml"))
array = ArrayData(collection, select="XX.*..HH?",
                  starttime="2026-11-01T10:00:00")

# Without positions in the metadata: local coordinates (east, north in m).
coordinates = {"XX.A01.": (0.0, 0.0), "XX.A02.": (15.0, 0.0),
               "XX.A03.": (-7.5, 13.0), "XX.A04.": (-7.5, -13.0)}
array = ArrayData(reader("data/*.sac"), coordinates=coordinates)

# One station at a time, from any source.
array = ArrayData()
array.add(collection, select="XX.A01..*")
array.add(collection, select="XX.A02..*")
array.remove("XX.A02.")
array.info()
print(array.sids, array.coordinates, array.delta, array.units)
