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
Unit tests for :mod:`shakelab.exploration.yarray.operational`.

Configuration checks, envelope shape and consistency, points outside
the search grid, recommendations, conservative products and a
regression against the August 2026 code.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.exploration.yarray.config import SearchGrid
from shakelab.exploration.yarray.operational import (
    LimitingFactor,
    OperationalConfig,
    OperationalCriteria,
    RecommendationSeverity,
    RecommendationType,
    assess_array,
)
from shakelab.exploration.yarray.performance import (
    EstimatorConfig,
    MonteCarloConfig,
    PerformanceConfig,
)
from shakelab.exploration.yarray.response import array_response

SEVEN = np.array([[0, 0], [30, 0], [0, 40], [25, 35], [-20, 10],
                  [10, -30], [-35, -20]], dtype=float)
GRID = SearchGrid(0.0, 0.2, 41, quantity="wavenumber", azimuth_size=72)
TEMPLATE = PerformanceConfig(
    snr_db=-3.0, grid=GRID, estimator=EstimatorConfig("capon"),
    monte_carlo=MonteCarloConfig(trials=30, realizations=20, seed=5),
)
AXES = {"frequencies": [3.0, 10.0], "velocities": [200.0, 400.0, 800.0],
        "azimuths": [0.0, 90.0]}


@pytest.fixture(scope="module")
def assessment():
    """Assessment of the seven-station array."""
    return assess_array(array_response(SEVEN), TEMPLATE,
                        OperationalConfig(**AXES))


def test_configuration_checks():
    config = OperationalConfig([1.0, 2.0], [100.0, 200.0])
    assert config.azimuths.size == 12

    with pytest.raises(ValueError):
        OperationalConfig([2.0, 1.0], [100.0])

    with pytest.raises(ValueError):
        OperationalConfig([1.0], [100.0], azimuths=[0.0, 360.0])

    with pytest.raises(ValueError):
        OperationalCriteria(minimum_local_probability=1.5)

    with pytest.raises(ValueError):
        OperationalCriteria(maximum_relative_velocity_rmse=-0.1)

    with pytest.raises(TypeError):
        assess_array(array_response(SEVEN), TEMPLATE, AXES)


def test_regression_against_august_code(assessment):
    local = [[[1.0, 0.9333333333333333], [1.0, 0.9],
              [0.9666666666666667, 0.9333333333333333]],
             [[0.0, 0.0], [1.0, 0.9666666666666667],
              [0.9333333333333333, 0.9666666666666667]]]
    factors = ["none", "aliasing", "resolution", "search_boundary",
               "search_boundary", "aliasing", "search_boundary",
               "search_boundary", "none", "none", "aliasing", "none"]
    assert np.allclose(assessment.local_probability, local)
    assert [item.value for item in assessment.limiting_factor.ravel()] == \
        factors


def test_envelope_consistency(assessment):
    shape = (2, 3, 2)

    for name in ("local_probability", "alias_probability",
                 "other_probability", "boundary_probability", "usable",
                 "limiting_factor", "relative_velocity_rmse"):
        assert getattr(assessment, name).shape == shape

    total = (assessment.local_probability + assessment.alias_probability
             + assessment.other_probability
             + assessment.boundary_probability)
    assert np.allclose(total, 1.0)
    usable = np.vectorize(lambda item: item is LimitingFactor.NONE)(
        assessment.limiting_factor
    )
    assert np.array_equal(usable, assessment.usable)
    assert assessment.resolution_radius > 0.0


def test_points_outside_the_grid(assessment):
    # 10 Hz at 200 m/s is k = 0.314, beyond the grid maximum of 0.2.
    assert np.all(assessment.boundary_probability[1, 0] == 1.0)
    assert np.all(np.isnan(assessment.grid_quantization_error[1, 0]))
    assert all(item is LimitingFactor.SEARCH_BOUNDARY
               for item in assessment.limiting_factor[1, 0])


def test_conservative_products(assessment):
    usable = assessment.conservative_usable
    assert np.array_equal(usable, np.all(assessment.usable, axis=2))
    assert np.array_equal(assessment.conservative_local_probability,
                          np.min(assessment.local_probability, axis=2))
    rmse = assessment.conservative_relative_velocity_rmse
    assert np.isnan(rmse[1, 0]) and np.isfinite(rmse[0, 0])

    for index in range(2):
        interval = assessment.usable_velocity_range(index)

        if np.any(usable[index]):
            assert interval[0] <= interval[1]
        else:
            assert interval is None


def test_recommendations(assessment):
    kinds = {item.category for item in assessment.recommendations}
    assert RecommendationType.SEARCH_RANGE in kinds
    severe = assess_array(
        array_response(SEVEN), TEMPLATE,
        OperationalConfig([30.0], [100.0], azimuths=[0.0]),
    )
    assert severe.recommendations[0].severity is \
        RecommendationSeverity.CRITICAL


def test_info(assessment, capsys):
    assessment.info()
    text = capsys.readouterr().out
    assert "OperationalAssessment" in text and "3 Hz" in text
