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
Operational assessment of an array.

Where will the f-k analysis be reliable? The performance experiment is
repeated over a grid of frequencies, velocities and azimuths, and each
point is marked usable when it meets explicit criteria; otherwise it
gets the factor that limits it. Recommendations summarize the result.

Run (from any folder)::

    python examples/exploration/yarray_15_operational.py
"""

import warnings
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_operational

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.ArrayData(DATA / "sesame_array.json")
response = ya.array_response(array)

# 1. The processing to assess (as in yarray_14_performance.py).
config = ya.PerformanceConfig(
    snr_db=-3.0,
    grid=ya.SearchGrid(0.0, 0.24, 41, quantity="wavenumber",
                       azimuth_size=72),
    estimator=ya.EstimatorConfig("capon"),
    monte_carlo=ya.MonteCarloConfig(trials=50, realizations=20, seed=1),
)

# 2. The domain of interest and the criteria of a usable point.
domain = ya.OperationalConfig(
    frequencies=[3.0, 5.0, 8.0, 12.0],
    velocities=[150.0, 250.0, 400.0, 600.0, 1000.0],
    azimuths=[0.0, 60.0, 120.0],        # default: 12 azimuths
    criteria=ya.OperationalCriteria(
        minimum_local_probability=0.8,
        maximum_alias_probability=0.05,
        maximum_relative_velocity_rmse=0.15,
    ),
)

assessment = ya.assess_array(response, config, domain)
assessment.info()

# 3. Usable velocities at each frequency (worst azimuth).
for index, frequency in enumerate(assessment.frequencies):
    print(frequency, assessment.usable_velocity_range(index))

# 4. Why the other points are not usable.
print(Counter(factor.value for factor in assessment.limiting_factor.flat))

for recommendation in assessment.recommendations:
    print(recommendation.severity.value, recommendation.message)

plot_operational(assessment)
plt.show()
