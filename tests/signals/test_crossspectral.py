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
Unit tests for :mod:`shakelab.signals.crossspectral`.

The coefficients are checked against the power spectral density of
white noise (``2 * sigma**2 * delta``, one-sided), against
``wft_transform(scaling="psd")`` of :mod:`shakelab.signals.spectrogram`,
between the direct and RFFT transforms, and the cross-spectral matrix
against the phase delay of a shifted sinusoid.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.signals.crossspectral import (
    coefficients,
    cross_spectral_matrix,
)
from shakelab.signals.spectrogram import wft_transform
from shakelab.signals.windowing import window_starts

DELTA = 0.01


def _noise(shape, sigma=2.0, seed=0):
    """White Gaussian noise."""
    return sigma * np.random.default_rng(seed).standard_normal(shape)


def test_shapes_and_window_starts():
    """Coefficients are (n_freq, n_windows, channels)."""
    data = _noise((3, 2000))
    freq, values = coefficients(data, DELTA, [2.0, 5.0], 400, overlap=0.5)

    np.testing.assert_array_equal(freq, [2.0, 5.0])
    assert values.shape == (2, window_starts(2000, 400, 0.5).size, 3)

    _, single = coefficients(data[0], DELTA, 2.0, 400)
    assert single.shape == (1, 5)


def test_white_noise_psd():
    """
    The diagonal is the one-sided PSD of white noise, 2 sigma**2 delta.

    Each estimate averages about 500 windows (relative standard error
    about 5%); the mean over 78 frequencies and two channels has a
    relative standard error below 0.5%: tolerance 2%.
    """
    data = _noise((2, 100000))
    freq = np.arange(1.0, 40.0, 0.5)
    _, values = coefficients(data, DELTA, freq, 400, overlap=0.5)
    psd = np.array([
        np.real(np.diag(cross_spectral_matrix(item))) for item in values
    ])

    assert np.mean(psd) == pytest.approx(2 * 4.0 * DELTA, rel=0.02)
    assert np.std(psd) / np.mean(psd) < 0.1


def test_agrees_with_wft_psd():
    """The PSD equals that of wft_transform with the same windows."""
    data = _noise(6000)
    freq, wtime, spectrum = wft_transform(
        data, DELTA, 400, step=200, window=("tukey", 0.1), centred=False,
        scaling="psd",
    )
    k = np.argmin(np.abs(freq - 3.0))
    _, values = coefficients(data, DELTA, freq[k], 400, step=200)

    np.testing.assert_allclose(
        np.abs(values[0]) ** 2, np.abs(spectrum[k]) ** 2, rtol=1e-10
    )


def test_rfft_nearest_equals_direct_at_bins():
    """At FFT bins the two transforms give the same coefficients."""
    data = _noise((4, 3000))
    direct = coefficients(data, DELTA, [2.5, 7.5], 400, overlap=0.5)[1]
    freq, rfft = coefficients(
        data, DELTA, [2.51, 7.49], 400, overlap=0.5, transform="rfft"
    )

    np.testing.assert_allclose(freq, [2.5, 7.5])
    np.testing.assert_allclose(rfft, direct, atol=1e-12)


def test_rfft_nearest_removes_repeated_bins():
    """Requested frequencies on the same bin give one coefficient."""
    freq, values = coefficients(
        _noise(3000), DELTA, [2.49, 2.51, 5.0], 400, transform="rfft"
    )

    np.testing.assert_allclose(freq, [2.5, 5.0])
    assert values.shape[0] == 2


def test_rfft_interpolation():
    """
    Interpolated coefficients approximate the direct ones.

    With the default eight-fold padding the cubic interpolation error
    is below 1% of the largest coefficient.
    """
    data = _noise((3, 3000))
    freq = [2.13, 6.71, 30.02]
    direct = coefficients(data, DELTA, freq, 400, overlap=0.5)[1]
    out, interp = coefficients(
        data, DELTA, freq, 400, overlap=0.5, transform="rfft",
        selection="interp",
    )

    np.testing.assert_array_equal(out, freq)
    assert np.max(np.abs(interp - direct)) < 0.01 * np.max(np.abs(direct))


def test_unscaled_coefficients():
    """scaling=None gives delta * DFT of the tapered window."""
    data = _noise(400)
    _, values = coefficients(data, DELTA, 2.5, 400, window="boxcar",
                             scaling=None)

    assert values[0, 0] == pytest.approx(DELTA * np.fft.rfft(data)[10])


def test_cross_spectrum_phase_of_a_delay():
    """
    A delay of tau between two channels gives the phase -2 pi f tau.

    Channel 1 is channel 0 delayed by 5 samples:
    C[0, 1] = S exp(+2 pi i f tau) with the kernel exp(-2 pi i f t).
    """
    base = _noise(20005)
    data = np.stack([base[5:], base[:-5]])
    _, values = coefficients(data, DELTA, 4.0, 500, overlap=0.5)
    matrix = cross_spectral_matrix(values[0])

    assert np.angle(matrix[0, 1]) == pytest.approx(
        2 * np.pi * 4.0 * 5 * DELTA, abs=0.05
    )
    np.testing.assert_allclose(matrix, matrix.conj().T)


def test_cross_spectral_matrix_variance():
    """The element variance matches an explicit computation."""
    rng = np.random.default_rng(1)
    x = rng.standard_normal((50, 3)) + 1j * rng.standard_normal((50, 3))
    y = rng.standard_normal((50, 2)) + 1j * rng.standard_normal((50, 2))

    matrix, spread = cross_spectral_matrix(x, y, variance=True)
    products = x[:, :, None] * y[:, None, :].conj()

    np.testing.assert_allclose(matrix, products.mean(axis=0))
    np.testing.assert_allclose(
        spread, np.var(products, axis=0, ddof=1), rtol=1e-12
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"freq": 0.0},
        {"freq": 50.0},
        {"freq": 2.0, "transform": "fft"},
        {"freq": 2.0, "selection": "linear"},
        {"freq": 2.0, "scaling": "amplitude"},
        {"freq": 2.0, "transform": "rfft", "padding": 0},
    ],
)
def test_invalid_arguments(kwargs):
    """Unknown options and frequencies outside (0, Nyquist) fail."""
    with pytest.raises(ValueError):
        coefficients(_noise(1000), DELTA, size=200, **kwargs)


def test_invalid_data_and_matrices():
    """Non-finite data, short series and wrong shapes fail."""
    with pytest.raises(ValueError):
        coefficients([1.0, np.nan, 2.0], DELTA, 2.0, 2)
    with pytest.raises(ValueError):
        coefficients(_noise(100), DELTA, 2.0, 200)
    with pytest.raises(ValueError):
        cross_spectral_matrix(np.ones(5))
    with pytest.raises(ValueError):
        cross_spectral_matrix(np.ones((1, 2)), variance=True)
    with pytest.raises(ValueError):
        cross_spectral_matrix(np.ones((3, 2)), np.ones((4, 2)))
