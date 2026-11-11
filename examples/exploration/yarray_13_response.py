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

Run (from any folder)::

    python examples/exploration/yarray_13_response.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import (
    plot_array_geometry,
    plot_array_response,
    plot_resolution,
    plot_search_domain,
    plot_velocity_limits,
)

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

# 1. Station coordinates (east, north, in metres): a centre and three
#    rings of three stations at 10, 25 and 60 m. The response accepts
#    plain coordinates, or an ArrayData.
coordinates = np.array([
    [0.0, 0.0],
    [0.0, 10.0], [8.7, -5.0], [-8.7, -5.0],
    [16.1, 19.2], [8.7, -24.6], [-24.6, 5.4],
    [59.1, 10.4], [-38.6, 46.0], [-20.5, -56.4],
])

# 2. Response on a wavenumber grid chosen from the geometry.
response = ya.array_response(coordinates)
response.info()

# 3. Main lobe, ambiguities and limits, at half power (default).
analysis = ya.analyze_response(response)
analysis.info()

limits = analysis.limits
print(f"resolved from {limits.resolution_radius:.3f} rad/m, "
      f"without aliases up to {limits.aliasing_radius:.3f} rad/m")

# 4. Phase velocities measurable at each frequency.
velocity = analysis.velocity_limits([2.0, 5.0, 10.0, 20.0])
for frequency, low, high in zip(velocity.frequencies, velocity.minimum,
                                velocity.maximum):
    print(f"{frequency:4.0f} Hz: {low:5.0f} to {high:5.0f} m/s")

# 5. Is a search up to 0.2 rad/m free of aliases?
domain = ya.analyze_search_domain(response, 0.2)
domain.info()

# 6. Uncertain positions (0.5 m of standard deviation): the expected
#    response loses coherence at high wavenumbers.
covariances = np.repeat([np.eye(2) * 0.25], len(coordinates), axis=0)
uncertain = ya.array_response(coordinates,
                              coordinate_covariances=covariances)
ya.analyze_response(uncertain, expected=True).info()

# 7. The array of the SESAME test data.
array = ya.ArrayData(DATA / "sesame_array.json")
ya.analyze_response(ya.array_response(array)).info()

figure, axes = plt.subplots(2, 3, figsize=(17, 10))
plot_array_geometry(response, ax=axes[0, 0])
plot_array_response(analysis, ax=axes[0, 1])
plot_resolution(analysis, ax=axes[0, 2])
plot_search_domain(analysis, domain, ax=axes[1, 0])
plot_velocity_limits(analysis.velocity_limits(np.geomspace(1, 30, 30)),
                     ax=axes[1, 1])
plot_array_geometry(array, labels=True, ax=axes[1, 2])
figure.tight_layout()
plt.show()
