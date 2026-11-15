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
Trace a P–P reflection from a specified horizontal interface.

Interface index 1 denotes the base of the first layer. The source
and receiver here are both at the surface.
"""

import matplotlib.pyplot as plt
import numpy as np

from shakelab.location.ttsolver import LayeredModel1D, reflected_path

model = LayeredModel1D([4000, np.inf], [4500, 6000], [2500, 3500])
ray = reflected_path(
    model, offset=6000, source_depth=0, receiver_depth=0,
    interface_index=1, phase="P",
)
segments = ray["segments"]
distance = np.cumsum([0, *(segment.x for segment in segments)])
depth = [segments[0].z_start, *(segment.z_end for segment in segments)]
print(f"Reflection time {ray['t']:.3f} s")
print(f"Reflection depth {ray['reflection_depth'] / 1000:.1f} km")

_, ax = plt.subplots(layout="constrained")
ax.plot(distance / 1000, np.asarray(depth) / 1000, "o-", label="P–P")
ax.axhline(ray["reflection_depth"] / 1000, color="0.6", ls="--")
ax.set(xlabel="Offset (km)", ylabel="Depth (km)")
ax.invert_yaxis()
ax.legend()
plt.show()
