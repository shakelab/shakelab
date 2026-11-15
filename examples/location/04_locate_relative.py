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
Locate a source from arrival-time differences between stations.

Synthetic P-wave travel times come from straight-line distances in a
homogeneous medium. Each observation is arrival at i minus arrival at j.
"""

from itertools import combinations
from math import dist

from shakelab.libutils.geodeticN import WgsPoint, geo_to_ecef
from shakelab.libutils.timeN import Date
from shakelab.location import HomogeneousEngine, Hypocenter, Station
from shakelab.location.relative_locator import (
    RelativePhaseObservation,
    locate_event_relative,
)

stations = {
    "A": Station("A", WgsPoint(13.00, 46.00, 200)),
    "B": Station("B", WgsPoint(13.10, 46.00, 220)),
    "C": Station("C", WgsPoint(13.00, 46.10, 210)),
    "D": Station("D", WgsPoint(13.10, 46.10, 230)),
    "E": Station("E", WgsPoint(13.05, 46.05, 250)),
}
truth = WgsPoint(13.05, 46.05, -8000)
source = geo_to_ecef(truth)
vp = 6000.0
travel_times = {
    code: dist(source, geo_to_ecef(station.location)) / vp
    for code, station in stations.items()
}
differences = [
    RelativePhaseObservation(i, j, "P", travel_times[i] - travel_times[j])
    for i, j in combinations(stations, 2)
]

initial = Hypocenter(WgsPoint(13.06, 46.04, -9000), Date())
event, result, _ = locate_event_relative(
    differences, stations, HomogeneousEngine(vp, 3500), initial,
    damping=0.0, tol_space_m=1e-4,
)
print(f"Converged: {result.converged}")
print(f"Longitude {event.lon:.4f}, latitude {event.lat:.4f}")
print(f"Depth {event.depth_km:.2f} km (target 8.00 km)")
print(f"RMS {result.rms:.3e} s for {result.num_observations} differences")
