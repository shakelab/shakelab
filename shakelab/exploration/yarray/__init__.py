# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Seismic-array processing (YArray).

Frequency-wavenumber analysis of synchronized seismic arrays, usable
from Python, from the command line and within the callbacks of the
SeedLink client. The package is under development: so far it contains
the array data model (:mod:`.data`).
"""

from __future__ import annotations

from shakelab.exploration.yarray.data import (
    ArrayData,
    ArrayDataWarning,
    StationData,
)

__all__ = [
    "ArrayData",
    "ArrayDataWarning",
    "StationData",
]
