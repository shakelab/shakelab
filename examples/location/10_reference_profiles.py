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
Choose a published 1D velocity profile and calculate a P arrival.

The catalog works offline. Depths and offsets are in metres.
"""

from shakelab.location.profiles import list_models, load_model, model_info
from shakelab.location.ttsolver import solve_traveltime

print("Available profiles:", ", ".join(list_models()))

name = "ogs-friuli"
info = model_info(name)
print(info["title"])
print(info["reference"], info["source_url"])
print("Depth reference:", info["depth_reference"])

model = load_model(name)
print("Layer tops (km):", model.z_top / 1000)
print("Vp (km/s):", model.vp / 1000)
print("Vs (km/s):", model.vs / 1000)
times, _ = solve_traveltime(
    model, [10000.0, 25000.0], source_depth=9000.0, phase="P",
    include_diving=False,
)
print("P times at 10 and 25 km (s):", times)
