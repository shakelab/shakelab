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
Array design and survey planning.

Before deploying an array (or before trusting its results) three
questions are answered, from the geometry alone to the processing:

1. geometric response (:mod:`.response`): which wavenumbers the layout
   resolves without ambiguity, hence which phase velocities at each
   frequency; the effect of uncertain station positions; whether a
   search domain is free of aliases. Two layouts are compared: the
   nine SESAME stations and a planned layout of three rings;
2. statistical performance of the estimator (:mod:`.performance`):
   Monte Carlo simulations of a plane wave in noise, processed as the
   data would be, classify each estimate as LOCAL, ALIAS, OTHER or
   BOUNDARY; a sweep shows the SNR and the number of windows needed;
3. operational assessment (:mod:`.operational`): the envelope of
   frequencies, velocities and azimuths where the analysis meets
   explicit criteria, with the limiting factors and recommendations.

The number of Monte Carlo trials is small so that the example runs in
a few seconds; use several hundreds for real planning.

Run from the root of the repository::

    python examples/exploration/yarray_planning.py
"""

from __future__ import annotations

import warnings
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration.yarray import (
    EstimatorConfig,
    MonteCarloConfig,
    OperationalConfig,
    OperationalCriteria,
    PerformanceConfig,
    PlaneWaveSource,
    SearchGrid,
    analyze_response,
    analyze_search_domain,
    array_response,
    assess_array,
    evaluate_performance,
    performance_sweep,
    read_array,
)
from shakelab.graphics.yarray import (
    plot_array_geometry,
    plot_array_response,
    plot_operational,
    plot_performance,
    plot_performance_sweep,
    plot_resolution,
    plot_velocity_limits,
)

warnings.filterwarnings("ignore", message="STLO")

FREQUENCIES = np.geomspace(2.0, 20.0, 12)

# ----------------------------------------------------------------------------
# 1. Geometric response
# ----------------------------------------------------------------------------

# a) The nine SESAME stations (the array file gives the coordinates;
#    plain coordinates work the same).
sesame = read_array("tests/exploration/yarray/data/sesame_array_9.json")
response = array_response(sesame)
analysis = analyze_response(response)
analysis.info()
limits = analysis.velocity_limits(FREQUENCIES)

# b) A planned layout: a centre and three rings of three stations at
#    10, 25 and 60 m, rotated by 40 degrees from ring to ring.
angles = np.radians(np.arange(3) * 120.0)
rings = [np.zeros((1, 2))] + [
    radius * np.column_stack((np.sin(angles + np.radians(40.0 * step)),
                              np.cos(angles + np.radians(40.0 * step))))
    for step, radius in enumerate((10.0, 25.0, 60.0))
]
planned = np.vstack(rings)
planned_analysis = analyze_response(array_response(planned))
planned_analysis.info()

# c) Uncertain positions (0.5 m standard deviation, e.g. hand-held
#    GNSS): the expected response loses coherence at high wavenumbers.
uncertain = array_response(
    planned, coordinate_covariances=np.repeat(
        np.eye(2)[None] * 0.25, len(planned), axis=0
    ),
)
expected = analyze_response(uncertain, expected=True)
expected.info()

# d) Is a search up to 0.2 rad/m free of aliases? (the differences of
#    the searched wavevectors reach 0.4 rad/m).
for name, item in (("SESAME", response), ("planned", uncertain)):
    domain = analyze_search_domain(item, 0.2)
    print(f"{name}: ", end="")
    domain.info()

# ----------------------------------------------------------------------------
# 2. Statistical performance of the estimators
# ----------------------------------------------------------------------------

# A wavenumber grid matching the alias-free domain of the SESAME
# layout; 20 windows per block, -8 dB of SNR per station.
grid = SearchGrid(0.0, 0.24, 41, quantity="wavenumber", azimuth_size=72)
source = PlaneWaveSource.from_velocity(250.0, 6.0, azimuth=60.0)
template = PerformanceConfig(
    snr_db=-8.0,
    grid=grid,
    estimator=EstimatorConfig("capon"),
    monte_carlo=MonteCarloConfig(trials=60, realizations=20, seed=1),
)

for method in ("bartlett", "capon", "music"):
    experiment = evaluate_performance(
        response, source, template.replace(estimator=EstimatorConfig(method))
    )
    velocity = experiment.velocity_statistics(6.0)
    print(f"{method:9s} local {experiment.outcome.local_probability:.2f}"
          f", alias {experiment.outcome.alias_probability:.2f}, velocity "
          f"{velocity.median:.0f} m/s ({velocity.lower_quantile:.0f} - "
          f"{velocity.upper_quantile:.0f}) for {velocity.true_velocity:.0f}")

capon = evaluate_performance(response, source, template)
capon.info()

#    SNR and number of windows: where does the estimate become local?
sweep = performance_sweep(
    response, source,
    template.replace(monte_carlo=template.monte_carlo.replace(
        trials=30, seed=2)),
    snr_db=[-12.0, -6.0, 0.0, 6.0], realizations=[12, 20, 40],
)
sweep.info()

# ----------------------------------------------------------------------------
# 3. Operational envelope
# ----------------------------------------------------------------------------

# At -3 dB; out of the radial range of the grid a point is BOUNDARY.
assessment = assess_array(
    response,
    template.replace(snr_db=-3.0, monte_carlo=template.monte_carlo.replace(
        trials=20, seed=3)),
    OperationalConfig(
        frequencies=[3.0, 5.0, 8.0, 12.0],
        velocities=np.geomspace(150.0, 1200.0, 7),
        azimuths=[0.0, 60.0, 120.0],
        criteria=OperationalCriteria(minimum_local_probability=0.8,
                                     maximum_relative_velocity_rmse=0.15),
    ),
)
assessment.info()
print(Counter(item.value for item in assessment.limiting_factor.ravel()))

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

figure, axes = plt.subplots(2, 4, figsize=(22, 10))
plot_array_geometry(sesame, labels=True, ax=axes[0, 0])
plot_array_response(analysis, ax=axes[0, 1])
plot_velocity_limits(limits, ax=axes[0, 2])
plot_array_geometry(uncertain, ax=axes[0, 3])
axes[0, 3].set_title("Planned layout")
plot_performance(capon, ax=axes[1, 0])
plot_performance_sweep(sweep, ax=axes[1, 1])
plot_operational(assessment, ax=axes[1, 2])
plot_resolution(planned_analysis, ax=axes[1, 3])
plot_resolution(analysis, ax=axes[1, 3])
axes[1, 3].legend(["planned", "SESAME"])
figure.tight_layout()
plt.show()
