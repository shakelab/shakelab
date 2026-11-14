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
Statistical performance of the estimators.

A plane wave in noise is simulated many times (Monte Carlo) and each
simulation is processed as real data would be. Each estimate is LOCAL
(near the true wavevector), ALIAS (in an ambiguity of the array),
OTHER (elsewhere, noise) or BOUNDARY (at the edge of the search grid).
The probabilities and the errors tell how reliable the analysis is for
a given array, estimator, signal-to-noise ratio and number of windows.
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    ArrayResponse,
    PerformanceConfig,
)
from shakelab.graphics.yarray import plot_performance, plot_performance_sweep

response = ArrayResponse(ArrayData("array.json"))

config = PerformanceConfig(snr_db=-8.0, method="capon",
                           grid_min=0.0, grid_max=0.24, grid_size=41,
                           azimuth_size=72,      # wavenumbers (rad/m)
                           trials=200, realizations=20, seed=1)
# or PerformanceConfig("planning.json", snr_db=-3.0)

result = response.performance(config, velocity=250.0, frequency=6.0,
                              azimuth=60.0)      # or wavenumber=0.15
result.info()
print(result.velocity_statistics())    # the frequency is in the source

for method in ("bartlett", "capon", "music"):
    test = response.performance(config.copy(method=method), wavenumber=0.15)
    print(method, test.outcome.local_probability)

sweep = response.performance_sweep(config, velocity=250.0, frequency=6.0,
                                   snr_db=[-12, -6, 0, 6],
                                   realizations=[10, 20, 40])
sweep.info()

figure, axes = plt.subplots(1, 2, figsize=(14, 6))
plot_performance(result, ax=axes[0])
plot_performance_sweep(sweep, ax=axes[1])
figure.tight_layout()
plt.show()
