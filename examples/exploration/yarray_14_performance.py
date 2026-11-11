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

Run (from any folder)::

    python examples/exploration/yarray_14_performance.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_performance, plot_performance_sweep

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array.json")
response = ya.array_response(array)

# 1. The experiment: the search grid (of wavenumbers), the estimator,
#    the SNR per station and the Monte Carlo sampling (trials, windows
#    per covariance matrix, seed for reproducible results).
config = ya.PerformanceConfig(
    snr_db=-8.0,
    grid=ya.SearchGrid(0.0, 0.24, 41, quantity="wavenumber",
                       azimuth_size=72),
    estimator=ya.EstimatorConfig("capon"),
    monte_carlo=ya.MonteCarloConfig(trials=200, realizations=20, seed=1),
)

# 2. The true wave: 250 m/s at 6 Hz, travelling towards N60E.
source = ya.PlaneWaveSource.from_velocity(250.0, 6.0, azimuth=60.0)

result = ya.evaluate_performance(response, source, config)
result.info()

outcome = result.outcome
print(outcome.local_probability, outcome.wilson_interval("local"))

velocity = result.velocity_statistics(6.0)
print(f"velocity {velocity.median:.0f} m/s "
      f"({velocity.lower_quantile:.0f} - {velocity.upper_quantile:.0f}), "
      f"relative RMSE {velocity.relative_rmse:.3f}")

# 3. The three estimators on the same simulations (same seed).
for method in ["bartlett", "capon", "music"]:
    test = ya.evaluate_performance(
        response, source, config.replace(estimator=ya.EstimatorConfig(method))
    )
    print(f"{method:9s} local {test.outcome.local_probability:.2f}")

# 4. How much SNR and how many windows are needed?
sweep = ya.performance_sweep(
    response, source,
    config.replace(monte_carlo=ya.MonteCarloConfig(trials=50, seed=2)),
    snr_db=[-12.0, -6.0, 0.0, 6.0],
    realizations=[10, 20, 40],
)
print(sweep.local_probability)

figure, axes = plt.subplots(1, 2, figsize=(14, 6))
plot_performance(result, ax=axes[0])
plot_performance_sweep(sweep, ax=axes[1])
figure.tight_layout()
plt.show()
