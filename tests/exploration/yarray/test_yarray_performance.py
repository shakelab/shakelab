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
Unit tests for :mod:`shakelab.exploration.yarray.performance`.

Configuration checks, reproducibility, outcome statistics (Wilson
intervals), behaviour with SNR, classification at the grid boundary,
velocity statistics, sweeps and a regression against the outcomes of
the August 2026 code (same seeds, same random draws).
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.exploration.yarray.config import SearchGrid
from shakelab.exploration.yarray.performance import (
    ClassificationConfig,
    EstimatorConfig,
    MonteCarloConfig,
    OutcomeClass,
    OutcomeSummary,
    PerformanceConfig,
    PlaneWaveSource,
    evaluate_performance,
    performance_sweep,
    theoretical_cross_spectrum,
)
from shakelab.exploration.yarray.response import array_response

SEVEN = np.array([[0, 0], [30, 0], [0, 40], [25, 35], [-20, 10],
                  [10, -30], [-35, -20]], dtype=float)
GRID = SearchGrid(0.0, 0.2, 41, quantity="wavenumber", azimuth_size=72)
SOURCE = PlaneWaveSource(0.15, 60.0)


@pytest.fixture(scope="module")
def response():
    """Response of the seven-station array."""
    return array_response(SEVEN)


def _config(method="bartlett", snr=-6.0, trials=100, realizations=20,
            seed=11, **options):
    """Performance configuration."""
    return PerformanceConfig(
        snr_db=snr, grid=GRID, estimator=EstimatorConfig(method),
        monte_carlo=MonteCarloConfig(trials, realizations, seed), **options,
    )


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------


def test_source():
    source = PlaneWaveSource(0.1, 450.0)
    assert source.azimuth == 90.0
    assert source.kx == pytest.approx(0.1) and abs(source.ky) < 1e-15
    wave = PlaneWaveSource.from_velocity(400.0, 8.0, 30.0)
    assert wave.wavenumber == pytest.approx(2 * np.pi * 8.0 / 400.0)

    with pytest.raises(ValueError):
        PlaneWaveSource(-0.1, 0.0)

    with pytest.raises(ValueError):
        PlaneWaveSource(0.1, 0.0, signal_power=0.0)


def test_configuration_checks():
    assert EstimatorConfig("capon").method.value == "capon"

    with pytest.raises(ValueError, match="unknown method"):
        EstimatorConfig("rtbf")

    with pytest.raises(ValueError):
        MonteCarloConfig(trials=1)

    with pytest.raises(ValueError):
        ClassificationConfig(resolution_level=1.0)

    with pytest.raises(ValueError, match="wavenumbers"):
        PerformanceConfig(0.0, SearchGrid(100.0, 1000.0, 10))

    with pytest.raises(ValueError, match="interior"):
        PerformanceConfig(0.0, GRID, classification=ClassificationConfig(
            boundary_margin=21))

    config = _config()
    assert config.snr == pytest.approx(10 ** -0.6)
    changed = config.replace(snr_db=3.0)
    assert changed.snr_db == 3.0 and config.snr_db == -6.0


def test_input_checks(response):
    with pytest.raises(ValueError, match="outside"):
        evaluate_performance(response, PlaneWaveSource(0.5, 0.0), _config())

    with pytest.raises(ValueError, match="realizations"):
        evaluate_performance(response, SOURCE,
                             _config("capon", realizations=5))

    line = array_response([[0, 0], [10, 0], [30, 0]])

    with pytest.raises(ValueError, match="two-dimensional"):
        evaluate_performance(line, SOURCE, _config())


def test_theoretical_cross_spectrum(response):
    matrix = theoretical_cross_spectrum(response, SOURCE, _config(snr=0.0))
    assert np.allclose(matrix, matrix.conj().T)
    assert np.allclose(np.diag(matrix).real, 2.0)
    assert np.allclose(np.abs(matrix[0, 1]), 1.0)


# ----------------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------------


def test_regression_against_august_code(response):
    reference = {
        "bartlett": ((90, 5, 3, 2), 0.011284288316733603),
        "capon": ((82, 11, 5, 2), 0.012120859244919237),
        "music": ((87, 6, 5, 2), 0.012214157714528547),
    }

    for method, (counts, rmse) in reference.items():
        result = evaluate_performance(response, SOURCE, _config(method))
        outcome = result.outcome
        assert (outcome.local_count, outcome.alias_count,
                outcome.other_count, outcome.boundary_count) == counts
        assert result.local_statistics.rmse == pytest.approx(rmse,
                                                             rel=1e-9)
        assert result.alias_counts.sum() == outcome.alias_count


def test_reproducibility(response):
    first = evaluate_performance(response, SOURCE, _config(trials=20))
    second = evaluate_performance(response, SOURCE, _config(trials=20))
    other = evaluate_performance(response, SOURCE,
                                 _config(trials=20, seed=12))
    assert np.array_equal(first.samples.power, second.samples.power)
    assert not np.array_equal(first.samples.power, other.samples.power)


def test_high_snr_is_local(response):
    result = evaluate_performance(response, SOURCE,
                                  _config(snr=20.0, trials=30))
    assert result.outcome.local_probability == 1.0
    assert result.local_statistics.rmse < 0.5 * \
        result.analysis.limits.resolution_radius
    assert result.global_rmse == pytest.approx(result.local_statistics.rmse)
    assert not result.reached_alias_families


def test_low_snr_degrades(response):
    good = evaluate_performance(response, SOURCE, _config(snr=0.0))
    bad = evaluate_performance(response, SOURCE, _config(snr=-15.0))
    assert bad.outcome.local_probability < good.outcome.local_probability
    assert bad.outcome.alias_probability + bad.outcome.other_probability \
        > good.outcome.alias_probability + good.outcome.other_probability


def test_boundary_outcome(response):
    edge = PlaneWaveSource(0.2, 0.0)
    result = evaluate_performance(
        response, edge, _config(snr=20.0, trials=10,
                                classification=ClassificationConfig(
                                    boundary_margin=1)),
    )
    assert result.outcome.boundary_probability == 1.0
    assert result.local_statistics is None
    assert np.all(result.samples.mask("boundary"))


def test_samples_and_grid_quantization(response):
    result = evaluate_performance(response, PlaneWaveSource(0.1012, 37.0),
                                  _config(snr=10.0, trials=10))
    samples = result.samples
    assert samples.size == 10
    assert np.allclose(samples.wavenumber, np.hypot(samples.kx, samples.ky))
    assert 0.0 < result.grid_quantization_error < GRID.values[1]
    assert np.allclose(np.linalg.norm(result.nearest_search_wavevector),
                       0.1)
    assert result.effective_realizations == 20.0
    assert result.noise_power == pytest.approx(0.1)


def test_outcome_summary():
    summary = OutcomeSummary(100, 90, 5, 5, 0)
    assert summary.count("alias") == 5
    assert summary.probability(OutcomeClass.LOCAL) == 0.9
    assert summary.standard_error("local") == pytest.approx(0.03)
    low, high = summary.wilson_interval("local")
    assert low < 0.9 < high
    low, high = summary.wilson_interval("boundary")
    assert low == 0.0 and 0.0 < high < 0.05
    low, high = OutcomeSummary(20, 20, 0, 0, 0).wilson_interval("local")
    assert 0.8 < low < 1.0 and high == 1.0

    with pytest.raises(ValueError):
        OutcomeSummary(10, 5, 5, 5, 0)

    with pytest.raises(ValueError):
        summary.count("unknown")


def test_velocity_statistics(response):
    result = evaluate_performance(response, SOURCE, _config(snr=0.0))
    slow = result.velocity_statistics(2.0)
    fast = result.velocity_statistics(8.0)
    assert slow.true_velocity == pytest.approx(2 * np.pi * 2.0 / 0.15)
    assert fast.true_velocity == pytest.approx(4 * slow.true_velocity)
    assert fast.relative_rmse == pytest.approx(slow.relative_rmse)
    assert slow.lower_quantile <= slow.median <= slow.upper_quantile
    assert slow.count == result.outcome.local_count

    with pytest.raises(ValueError):
        result.velocity_statistics(0.0)


def test_sweep(response):
    config = _config(trials=20, seed=4)
    sweep = performance_sweep(response, SOURCE, config,
                              snr_db=[-10.0, 0.0, 10.0],
                              realizations=[10, 20], seed_stride=7)
    assert sweep.shape == (3, 2)
    assert sweep.local_probability.shape == (3, 2)
    point = sweep.result(2, 1)
    assert point.config.snr_db == 10.0
    assert point.config.monte_carlo.realizations == 20
    assert point.config.monte_carlo.seed == 4 + 5 * 7
    assert sweep.local_probability[2, 1] >= sweep.local_probability[0, 0]
    single = evaluate_performance(response, SOURCE, point.config)
    assert np.array_equal(single.samples.power, point.samples.power)

    with pytest.raises(ValueError):
        performance_sweep(response, SOURCE, config, snr_db=[0.0],
                          realizations=[1.5])


def test_info(response, capsys):
    result = evaluate_performance(response, SOURCE, _config(trials=10))
    result.info()
    performance_sweep(response, SOURCE, _config(trials=10), [0.0],
                      [20]).info()
    text = capsys.readouterr().out
    assert "PerformanceResult" in text and "PerformanceSweepResult" in text
