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
Compare direct arrivals with first arrivals over a faster half-space.

All offsets and depths are in metres. The faster lower layer makes
the head wave the first arrival at sufficiently large offsets.
"""

import matplotlib.pyplot as plt
import numpy as np

from shakelab.location.ttsolver import LayeredModel1D, solve_traveltime

model = LayeredModel1D([4000, np.inf], [4500, 6000], [2500, 3500])
offsets = np.linspace(0, 50000, 101)
direct, _ = solve_traveltime(model, offsets, 1000, arrivals="direct")
first, branches = solve_traveltime(
    model, offsets, 1000, include_diving=False, return_branches=True,
)

# One list of candidate branches is returned for each offset.
first_kinds = [min(rays, key=lambda ray: ray["t"])["kind"]
               for rays in branches]
head = np.asarray(first_kinds) == "head"
print(f"Head waves are first at {head.sum()} of {offsets.size} offsets")

_, ax = plt.subplots(layout="constrained")
ax.plot(offsets / 1000, direct, label="Direct P")
ax.plot(offsets[head] / 1000, first[head], label="Head wave (first)")
ax.set(xlabel="Offset (km)", ylabel="Travel time (s)")
ax.legend()
plt.show()
