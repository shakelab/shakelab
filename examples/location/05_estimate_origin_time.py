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
Estimate origin time once the hypocenter is known.

The Cartesian positions are in metres. Two synthetic P picks have
standard errors of 0.1 and 0.2 s, so their weights are 10 and 5.
"""

from shakelab.libutils.geodeticN import MetricPoint
from shakelab.libutils.timeN import Date
from shakelab.location import (
    CartHypocenter,
    CartStation,
    HomogeneousEngine,
    PhasePick,
)
from shakelab.location.relative_locator import estimate_origin_time

source = CartHypocenter(MetricPoint(0, 0, -4000), t0=0)
stations = {
    "A": CartStation("A", MetricPoint(3000, 0, 0)),
    "B": CartStation("B", MetricPoint(0, 0, 0)),
}
origin = Date("2024-03-27T21:19:37Z")
picks = [
    PhasePick("A", "P", origin + 1.1, weight=10),
    PhasePick("B", "P", origin + 0.8, weight=5),
]
result = estimate_origin_time(
    picks, stations, HomogeneousEngine(5000, 3000), source,
)
print(f"Estimated origin: {result.t0}")
print(f"Shift from synthetic origin: {result.t0 - origin:+.2f} s")
print(f"Arrival residuals (s): {result.residuals}")
