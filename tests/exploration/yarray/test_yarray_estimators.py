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
Unit tests for :mod:`shakelab.exploration.yarray.estimators`.

The estimators are checked on the exact matrix of one plane wave in
white noise, ``C = S a0 a0^H + s I``, for which Bartlett gives ``S`` at
the true wavevector and Capon gives ``S + s / n`` (Capon, 1969); the
empirical standard errors are checked against Monte Carlo simulations.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.exploration.yarray.estimators import (
    BeamformingMethod,
    bartlett_power,
    bartlett_power_statistics,
    beam_power,
    capon_power,
    capon_power_statistics,
    check_realizations,
    music_power,
    regularized_matrix,
    steering_vectors,
)

RNG = np.random.default_rng(7)
COORDINATES = RNG.uniform(-100.0, 100.0, (8, 2))
K0 = (0.02, -0.01)


def _plane_wave_matrix(signal=2.0, noise=0.5):
    """Exact cross-spectral matrix of one plane wave in white noise."""
    a0 = steering_vectors(COORDINATES, *K0)
    return signal * np.outer(a0, a0.conj()) + noise * np.eye(8), a0


def _grid():
    """Wavevector grid containing the true wavevector."""
    kx, ky = np.meshgrid(np.linspace(-0.04, 0.04, 41),
                         np.linspace(-0.04, 0.04, 41), indexing="ij")
    return kx, ky


def _realizations(count, signal=2.0, noise=0.5, rng=RNG):
    """Random Fourier coefficients of one plane wave in white noise."""
    a0 = steering_vectors(COORDINATES, *K0)
    amplitude = np.sqrt(signal / 2.0) * (
        rng.standard_normal(count) + 1j * rng.standard_normal(count)
    )
    white = np.sqrt(noise / 2.0) * (
        rng.standard_normal((count, 8)) + 1j * rng.standard_normal((count, 8))
    )
    return amplitude[:, None] * a0 + white


def test_steering_vectors():
    """Phases -k.r, shape kx.shape + (n,)."""
    kx, ky = _grid()
    steering = steering_vectors(COORDINATES, kx, ky)

    assert steering.shape == (41, 41, 8)
    np.testing.assert_allclose(
        steering[30, 15],
        np.exp(-1j * (kx[30, 15] * COORDINATES[:, 0]
                      + ky[30, 15] * COORDINATES[:, 1])),
    )

    with pytest.raises(ValueError):
        steering_vectors(COORDINATES, kx, ky[:3])


def test_bartlett_and_capon_at_the_true_wavevector():
    """Both give S + s / n at the true wavevector (no loading)."""
    matrix, a0 = _plane_wave_matrix()

    assert bartlett_power(matrix, a0) == pytest.approx(2.0 + 0.5 / 8)
    assert capon_power(matrix, a0, 0.0) == pytest.approx(2.0 + 0.5 / 8)


def test_maxima_at_the_true_wavevector():
    """All estimators peak at the true wavevector."""
    matrix, _ = _plane_wave_matrix()
    steering = steering_vectors(COORDINATES, *_grid())

    for power in (
        bartlett_power(matrix, steering),
        capon_power(matrix, steering, 0.001),
        music_power(matrix, steering, 1),
    ):
        assert np.unravel_index(np.argmax(power), power.shape) == (30, 15)


def test_capon_is_sharper_than_bartlett():
    """Capon is never above Bartlett (Lacoss et al., 1969)."""
    matrix, _ = _plane_wave_matrix()
    steering = steering_vectors(COORDINATES, *_grid())

    assert np.all(
        capon_power(matrix, steering, 0.0)
        <= bartlett_power(matrix, steering) * (1 + 1e-12)
    )


def test_music_definition():
    """MUSIC is 1 / (a^H E E^H a) with the noise eigenvectors E."""
    matrix, _ = _plane_wave_matrix()
    steering = steering_vectors(COORDINATES, *_grid())[3, 4]
    noise = np.linalg.eigh(matrix)[1][:, :7]

    assert music_power(matrix, steering, 1) == pytest.approx(
        1.0 / np.real(steering.conj() @ noise @ noise.conj().T @ steering)
    )

    with pytest.raises(ValueError):
        music_power(matrix, steering, 8)


def test_regularized_matrix():
    """Loading relative to the mean diagonal; Hermitian part."""
    matrix = np.array([[2.0, 1.0j], [0.0, 4.0]])
    loaded = regularized_matrix(matrix, 0.1)

    np.testing.assert_allclose(loaded, [[2.3, 0.5j], [-0.5j, 4.3]])

    with pytest.raises(ValueError):
        regularized_matrix(np.zeros((2, 2)), 0.1)
    with pytest.raises(ValueError):
        regularized_matrix(matrix, -1.0)


def test_bartlett_attenuation():
    """
    Unit coherence changes nothing; zero coherence leaves the diagonal.

    With g = 0 the expected power is trace(C) / n**2 at every
    wavevector.
    """
    matrix, _ = _plane_wave_matrix()
    steering = steering_vectors(COORDINATES, *_grid())

    np.testing.assert_allclose(
        bartlett_power(matrix, steering, np.ones(steering.shape)),
        bartlett_power(matrix, steering),
    )
    np.testing.assert_allclose(
        bartlett_power(matrix, steering, np.zeros(steering.shape)),
        np.real(np.trace(matrix)) / 64,
    )


def test_statistics_means_equal_estimators():
    """The statistics functions return the ordinary powers."""
    values = _realizations(50)
    matrix = values.T @ values.conj() / 50
    steering = steering_vectors(COORDINATES, *_grid())

    power, error = bartlett_power_statistics(values, steering, 50.0)
    np.testing.assert_allclose(power, bartlett_power(matrix, steering))
    assert np.all(error > 0.0)

    power, error = capon_power_statistics(values, steering, 50.0, 0.01)
    np.testing.assert_allclose(power, capon_power(matrix, steering, 0.01))
    assert np.all(error > 0.0)


@pytest.mark.parametrize("method", ["bartlett", "capon"])
def test_standard_error_against_monte_carlo(method):
    """
    The empirical standard error matches the spread of independent
    estimates at the true wavevector.

    400 Monte Carlo sets of 60 realizations; the mean predicted error
    must match the observed standard deviation within 20% (the
    estimate of a standard deviation from 400 samples has a relative
    error of about 4%; the delta method is a first-order
    approximation).
    """
    rng = np.random.default_rng(11)
    a0 = steering_vectors(COORDINATES, *K0)
    powers, errors = [], []

    for _ in range(400):
        values = _realizations(60, rng=rng)

        if method == "bartlett":
            power, error = bartlett_power_statistics(values, a0, 60.0)
        else:
            power, error = capon_power_statistics(values, a0, 60.0, 0.01)

        powers.append(power)
        errors.append(error)

    assert np.mean(errors) == pytest.approx(np.std(powers), rel=0.2)


def test_check_realizations():
    """Capon and MUSIC need more realizations than the dimension."""
    check_realizations(5, 8, "bartlett")
    check_realizations(5, 8, "capon", strict=False)

    with pytest.raises(ValueError, match="more realizations"):
        check_realizations(8, 8, BeamformingMethod.MUSIC)


def test_beam_power_dispatch():
    """beam_power selects the estimator by name."""
    matrix, a0 = _plane_wave_matrix()

    assert beam_power(matrix, a0, 20, "capon", diagonal_loading=0.0) == (
        pytest.approx(capon_power(matrix, a0, 0.0))
    )
    assert beam_power(matrix, a0, 2, "bartlett") == pytest.approx(
        bartlett_power(matrix, a0)
    )

    with pytest.raises(ValueError, match="Bartlett"):
        beam_power(matrix, a0, 20, "capon", attenuation=np.ones(8))
    with pytest.raises(ValueError, match="method"):
        beam_power(matrix, a0, 20, "beam")
