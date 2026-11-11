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
Geometric response of an array.

The response depends only on the station positions: it tells which
wavevectors the array resolves (main lobe) and which it confuses
(aliases), hence the phase velocities it can measure at each
frequency. It is computed before a survey, from the planned layout.
"""

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration.yarray import ArrayData, ArrayResponse
from shakelab.graphics.yarray import plot_array_response, plot_velocity_limits

array = ArrayData("array.json")
response = ArrayResponse(array)          # or coordinates (n, 2|3); grid_max
response.info()                          # geometry and grid

analysis = response.analyze()            # resolution_level, ambiguity_level
analysis.info()
limits = analysis.velocity_limits([2.0, 5.0, 10.0, 20.0])
response.search_domain(0.2).info()       # is |k| <= 0.2 rad/m safe?

# A planned layout (east, north in m), positions known within 0.5 m.
coordinates = [[0, 0], [0, 10], [9, -5], [-9, -5], [16, 19], [9, -25],
               [-25, 5], [59, 10], [-39, 46], [-20, -56]]
covariances = np.repeat([np.eye(2) * 0.25], len(coordinates), axis=0)
uncertain = ArrayResponse(coordinates, coordinate_covariances=covariances)
uncertain.analyze(expected=True).info()

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_array_response(analysis, ax=axes[0])
plot_velocity_limits(limits, ax=axes[1])
figure.tight_layout()
plt.show()
