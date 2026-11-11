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

Sources and input checks, reproducibility, outcome statistics (Wilson
intervals), behaviour with SNR, classification at the grid boundary,
velocity statistics, sweeps and a regression against the outcomes of
the August 2026 code (same seeds, same random draws). The
configuration is tested in ``test_yarray_config.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.exploration.yarray.config import PerformanceConfig
from shakelab.exploration.yarray.performance import (
    OutcomeClass,
    OutcomeSummary,
    PlaneWaveSource,
    _cross_spectrum,
)
from shakelab.exploration.yarray.response import ArrayResponse

SEVEN = np.array([[0, 0], [30, 0], [0, 40], [25, 35], [-20, 10],
                  [10, -30], [-35, -20]], dtype=float)
GRID = {"grid_min": 0.0, "grid_max": 0.2, "grid_size": 41,
        "azimuth_size": 72}
SOURCE = {"wavenumber": 0.15, "azimuth": 60.0}


@pytest.fixture(scope="module")
def response():
    """Response of the seven-station array."""
    return ArrayResponse(SEVEN)


def _config(method="bartlett", snr=-6.0, trials=100, realizations=20,
            seed=11, **options):
    """Performance configuration."""
    return PerformanceConfig(snr_db=snr, method=method, trials=trials,
                             realizations=realizations, seed=seed,
                             **GRID, **options)


# ----------------------------------------------------------------------------
# Sources and checks
# ----------------------------------------------------------------------------


def test_source():
    source = PlaneWaveSource(0.1, 450.0)
    assert source.azimuth == 90.0
    assert source.kx == pytest.approx(0.1) and abs(source.ky) < 1e-15
    assert source.frequency is None and source.velocity is None
    assert PlaneWaveSource(0.1, 0.0, frequency=2.0).velocity == \
        pytest.approx(4 * np.pi / 0.1)

    with pytest.raises(ValueError):
        PlaneWaveSource(-0.1, 0.0)

    with pytest.raises(ValueError):
        PlaneWaveSource(0.1, 0.0, signal_power=0.0)

    with pytest.raises(ValueError):
        PlaneWaveSource(0.1, 0.0, frequency=-1.0)


def test_source_parameters(response):
    config = _config(trials=5)
    wave = response.performance(config, velocity=400.0, frequency=8.0,
                                azimuth=30.0).source
    assert wave.wavenumber == pytest.approx(2 * np.pi * 8.0 / 400.0)
    assert wave.frequency == 8.0 and wave.azimuth == 30.0
    assert response.performance(config, wavenumber=0.1).source.azimuth == 0

    with pytest.raises(TypeError, match="either"):
        response.performance(config)

    with pytest.raises(TypeError, match="either"):
        response.performance(config, wavenumber=0.1, velocity=400.0,
                             frequency=8.0)

    with pytest.raises(TypeError, match="frequency"):
        response.performance(config, velocity=400.0)

    with pytest.raises(ValueError):
        response.performance(config, velocity=-400.0, frequency=8.0)


def test_input_checks(response):
    with pytest.raises(ValueError, match="outside"):
        response.performance(_config(), wavenumber=0.5)

    with pytest.raises(ValueError, match="realizations"):
        response.performance(_config("capon", realizations=5), **SOURCE)

    with pytest.raises(ValueError, match="needs snr_db"):
        response.performance(PerformanceConfig(**GRID), **SOURCE)

    with pytest.raises(ValueError, match="needs grid"):
        response.performance(PerformanceConfig(snr_db=0.0), **SOURCE)

    with pytest.raises(TypeError):
        response.performance({"snr_db": 0.0}, **SOURCE)

    line = ArrayResponse([[0, 0], [10, 0], [30, 0]])

    with pytest.raises(ValueError, match="two-dimensional"):
        line.performance(_config(), **SOURCE)


def test_config_is_copied(response):
    config = _config(trials=5)
    result = response.performance(config, **SOURCE)
    config.update(trials=7)
    assert result.config.trials == 5 and result.config is not config


def test_cross_spectrum(response):
    matrix = _cross_spectrum(response, PlaneWaveSource(0.15, 60.0),
                             _config(snr=0.0))
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
        result = response.performance(_config(method), **SOURCE)
        outcome = result.outcome
        assert (outcome.local_count, outcome.alias_count,
                outcome.other_count, outcome.boundary_count) == counts
        assert result.local_statistics.rmse == pytest.approx(rmse,
                                                             rel=1e-9)
        assert result.alias_counts.sum() == outcome.alias_count


def test_reproducibility(response):
    first = response.performance(_config(trials=20), **SOURCE)
    second = response.performance(_config(trials=20), **SOURCE)
    other = response.performance(_config(trials=20, seed=12), **SOURCE)
    assert np.array_equal(first.samples.power, second.samples.power)
    assert not np.array_equal(first.samples.power, other.samples.power)


def test_high_snr_is_local(response):
    result = response.performance(_config(snr=20.0, trials=30), **SOURCE)
    assert result.outcome.local_probability == 1.0
    assert result.local_statistics.rmse < 0.5 * \
        result.analysis.limits.resolution_radius
    assert result.global_rmse == pytest.approx(result.local_statistics.rmse)
    assert not result.reached_alias_families


def test_low_snr_degrades(response):
    good = response.performance(_config(snr=0.0), **SOURCE)
    bad = response.performance(_config(snr=-15.0), **SOURCE)
    assert bad.outcome.local_probability < good.outcome.local_probability
    assert bad.outcome.alias_probability + bad.outcome.other_probability \
        > good.outcome.alias_probability + good.outcome.other_probability


def test_boundary_outcome(response):
    result = response.performance(
        _config(snr=20.0, trials=10, boundary_margin=1), wavenumber=0.2
    )
    assert result.outcome.boundary_probability == 1.0
    assert result.local_statistics is None
    assert np.all(result.samples.mask("boundary"))


def test_samples_and_grid_quantization(response):
    result = response.performance(_config(snr=10.0, trials=10),
                                  wavenumber=0.1012, azimuth=37.0)
    samples = result.samples
    assert samples.size == 10
    assert np.allclose(samples.wavenumber, np.hypot(samples.kx, samples.ky))
    assert 0.0 < result.grid_quantization_error < 0.2 / 40
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
    result = response.performance(_config(snr=0.0), **SOURCE)
    slow = result.velocity_statistics(2.0)
    fast = result.velocity_statistics(8.0)
    assert slow.true_velocity == pytest.approx(2 * np.pi * 2.0 / 0.15)
    assert fast.true_velocity == pytest.approx(4 * slow.true_velocity)
    assert fast.relative_rmse == pytest.approx(slow.relative_rmse)
    assert slow.lower_quantile <= slow.median <= slow.upper_quantile
    assert slow.count == result.outcome.local_count

    with pytest.raises(ValueError):
        result.velocity_statistics(0.0)

    with pytest.raises(ValueError, match="frequency"):
        result.velocity_statistics()

    known = response.performance(_config(snr=0.0), velocity=250.0,
                                 frequency=6.0, azimuth=60.0)
    assert known.velocity_statistics().true_velocity == pytest.approx(250.0)


def test_sweep(response):
    config = _config(trials=20, seed=4)
    sweep = response.performance_sweep(config, snr_db=[-10.0, 0.0, 10.0],
                                       realizations=[10, 20], **SOURCE)
    assert sweep.shape == (3, 2)
    assert sweep.local_probability.shape == (3, 2)
    point = sweep.result(2, 1)
    assert point.config.snr_db == 10.0
    assert point.config.realizations == 20
    assert point.config.seed == 4 + 5 * 1000003
    assert config.realizations == 20 and config.seed == 4
    assert sweep.local_probability[2, 1] >= sweep.local_probability[0, 0]
    single = response.performance(point.config, **SOURCE)
    assert np.array_equal(single.samples.power, point.samples.power)

    with pytest.raises(ValueError):
        response.performance_sweep(config, [0.0], [1.5], **SOURCE)

    template = PerformanceConfig(trials=5, **GRID)
    assert response.performance_sweep(template, [0.0], [20],
                                      **SOURCE).shape == (1, 1)


def test_info(response, capsys):
    result = response.performance(_config(trials=10), **SOURCE)
    result.info()
    response.performance_sweep(_config(trials=10), [0.0], [20],
                               **SOURCE).info()
    text = capsys.readouterr().out
    assert "PerformanceResult" in text and "PerformanceSweepResult" in text
