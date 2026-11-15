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
Plot P and S velocities from the Friuli 1D model file.

The CSV stores layer thickness in metres and velocity in metres per second.
Layer boundaries follow from the cumulative thickness.
"""

from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.location.profiles import read_model

profile = Path("data/Model/friuli_velocity_profile_custom.csv")
model = read_model(profile)
edges = [0, *model.z_bottom / 1000]

_, ax = plt.subplots(layout="constrained")
ax.stairs(model.vp / 1000, edges, orientation="horizontal",
          baseline=None, label="P")
ax.stairs(model.vs / 1000, edges, orientation="horizontal",
          baseline=None, label="S")
ax.set(xlabel="Velocity (km/s)", ylabel="Depth (km)")
ax.invert_yaxis()
ax.legend()
plt.show()
