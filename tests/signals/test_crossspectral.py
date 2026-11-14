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
against the phase delay of a shifted sinusoid. :class:`CrossSpectra` is
checked against the matrix of mixed white noise (``A diag(2 sigma**2
delta) A^T``), against :func:`cross_spectral_matrix` block by block,
and against a rank-one wavefield for the decomposition.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.signals.crossspectral import (
    CrossSpectra,
    coefficients,
    cross_spectral_matrix,
)
from shakelab.signals.spectrogram import wft_transform
from shakelab.signals.windowing import window_size, window_starts

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


# ----------------------------------------------------------------------------
# CrossSpectra
# ----------------------------------------------------------------------------


def test_cross_spectra_blocks_match_cross_spectral_matrix():
    """Each block is the matrix of its windows; times are block centres."""
    data = _noise((3, 6000), seed=1)
    spectra = CrossSpectra.from_data(
        data, DELTA, [2.0, 5.0], duration=4.0, overlap=0.5,
        windows_per_block=10, block_overlap=0.5,
    )
    starts = window_starts(6000, 400, overlap=0.5)
    blocks = window_starts(starts.size, 10, overlap=0.5)
    _, values = coefficients(data, DELTA, [2.0, 5.0], 400, overlap=0.5)

    assert spectra.nchannels == 3
    np.testing.assert_array_equal(spectra.nblocks, [blocks.size] * 2)
    np.testing.assert_array_equal(spectra.realizations(1), 10)

    for k in range(2):
        matrices = spectra.matrices(k)

        for row, first in enumerate(blocks):
            expected = cross_spectral_matrix(values[k, first:first + 10])
            np.testing.assert_allclose(matrices[row], expected,
                                       rtol=1e-12)

    centre = (starts[blocks[1]] + starts[blocks[1] + 9] + 399) / 2 * DELTA
    assert spectra.times(0)[1] == pytest.approx(centre)


def test_cross_spectra_mixed_white_noise():
    """The mean matrix of mixed white noise is A diag(2 s**2 dt) A^T."""
    sigma = np.array([1.0, 2.0, 0.5])
    mixing = np.array([[1.0, 0.5, 0.0], [0.2, 1.0, 0.3], [0.0, -0.4, 1.0]])
    sources = sigma[:, None] * np.random.default_rng(2).standard_normal(
        (3, 200000)
    )
    spectra = CrossSpectra.from_data(
        mixing @ sources, DELTA, np.linspace(2.0, 40.0, 20), duration=5.0,
        overlap=0.5, window="hann",
    )
    expected = mixing @ np.diag(2.0 * sigma ** 2 * DELTA) @ mixing.T
    scale = np.abs(expected).max()

    # About 4% scatter at each frequency: compare the frequency mean.
    np.testing.assert_allclose(spectra.mean().mean(axis=0), expected,
                               atol=0.03 * scale)


def test_cross_spectra_cycles_and_subfrequencies():
    """Windows in cycles change with frequency; sub-frequencies average."""
    data = _noise((2, 8000), seed=3)
    spectra = CrossSpectra.from_data(
        data, DELTA, [1.0, 4.0], cycles=10, overlap=0.5,
        windows_per_block=5, subfrequencies=3, bandwidth=0.2,
    )
    sizes = [window_size(DELTA, cycles=10, frequency=f) for f in (1.0, 4.0)]
    counts = [
        window_starts(window_starts(8000, n, 0.5).size, 5).size
        for n in sizes
    ]

    np.testing.assert_array_equal(spectra.nblocks, counts)
    np.testing.assert_array_equal(spectra.realizations(0), 15)

    targets = 4.0 * np.exp(np.array([-0.1, 0.0, 0.1]))
    _, values = coefficients(data, DELTA, targets, sizes[1], overlap=0.5)
    series = values[:, :5, :].reshape(-1, 2)
    np.testing.assert_allclose(spectra.matrices(1)[0],
                               cross_spectral_matrix(series), rtol=1e-12)


def test_cross_spectra_window_longer_than_data():
    """A frequency whose window does not fit has no blocks."""
    spectra = CrossSpectra.from_data(
        _noise((2, 1000), seed=4), DELTA, [0.5, 5.0], cycles=10,
    )

    np.testing.assert_array_equal(spectra.nblocks, [0, 1])
    assert spectra.matrices(0).shape == (0, 2, 2)
    assert np.all(np.isnan(spectra.mean()[0]))


def test_cross_spectra_coherence():
    """Unit coherence of scaled copies, times the mean power."""
    base = _noise(4000, seed=5)
    data = np.vstack([base, 3.0 * base, _noise(4000, seed=6)])
    spectra = CrossSpectra.from_data(data, DELTA, [2.0, 6.0], duration=2.0)
    coherence = spectra.coherence()
    original = spectra.matrices(0)[0]
    result = coherence.matrices(0)[0]
    mean_power = np.real(np.trace(original)) / 3

    np.testing.assert_allclose(np.real(np.diag(result)), mean_power)
    assert abs(result[0, 1]) == pytest.approx(mean_power)
    assert abs(result[0, 2]) < mean_power
    np.testing.assert_array_equal(coherence.times(1), spectra.times(1))

    with pytest.raises(ValueError):
        CrossSpectra.from_data(
            np.vstack([base, np.zeros(4000)]), DELTA, 2.0, duration=2.0
        ).coherence()


def test_cross_spectra_decompose_rank_one():
    """A rank-one wavefield gives one dominant eigenvalue and its shape."""
    shape = np.array([0.2, 0.6, 1.0, 0.7])
    signal = _noise(20000, seed=7)
    noise = 0.05 * _noise((4, 20000), seed=8)
    spectra = CrossSpectra.from_data(
        shape[:, None] * signal + noise, DELTA, [2.0, 8.0], duration=5.0,
        overlap=0.5, windows_per_block=20,
    )
    values, vectors = spectra.decompose(1)
    matrices = spectra.matrices(1)

    assert values.shape == (spectra.nblocks[1], 4)
    assert np.all(np.diff(values, axis=1) <= 0.0)
    assert np.all(values[:, 1] < 0.01 * values[:, 0])

    first = vectors[:, :, 0]
    mac = np.abs(first @ shape) ** 2 / (shape @ shape)
    assert np.all(mac > 0.999)

    rebuilt = np.einsum("bij,bj,bkj->bik", vectors, values, vectors.conj())
    np.testing.assert_allclose(rebuilt, matrices,
                               atol=1e-10 * np.abs(matrices).max())


def test_cross_spectra_summary_and_copy(capsys):
    """Summary, info, labels and independent copies."""
    spectra = CrossSpectra.from_data(
        _noise((2, 2000), seed=9), DELTA, [2.0, 3.0], duration=2.0,
        channels=["XX.A01..HHZ", "XX.A02..HHZ"],
    )
    copy = spectra.copy()

    assert "2 channels" in spectra.summary()
    spectra.info()
    assert "XX.A01..HHZ" in capsys.readouterr().out
    assert copy.channels == spectra.channels
    np.testing.assert_array_equal(copy.matrices(0), spectra.matrices(0))
    assert spectra.matrices(-1).shape == spectra.matrices(1).shape

    with pytest.raises(IndexError):
        spectra.matrices(2)
    with pytest.raises(TypeError):
        spectra.times(1.0)


@pytest.mark.parametrize(
    "kwargs, error",
    [
        ({}, ValueError),
        ({"duration": 2.0, "cycles": 5}, ValueError),
        ({"duration": 2.0, "subfrequencies": 3}, ValueError),
        ({"duration": 2.0, "bandwidth": 0.1}, ValueError),
        ({"duration": 2.0, "subfrequencies": 2.0, "bandwidth": 0.1},
         TypeError),
        ({"duration": 2.0, "windows_per_block": 0}, ValueError),
        ({"duration": 2.0, "windows_per_block": 2.5}, TypeError),
        ({"duration": 2.0, "block_overlap": 1.0}, ValueError),
        ({"duration": 2.0, "channels": ["a"]}, ValueError),
        ({"duration": 2.0, "channels": "ab"}, TypeError),
    ],
)
def test_cross_spectra_invalid_arguments(kwargs, error):
    """Inconsistent windows, blocks, sub-frequencies and labels fail."""
    with pytest.raises(error):
        CrossSpectra.from_data(_noise((2, 2000)), DELTA, 2.0, **kwargs)


def test_cross_spectra_invalid_construction():
    """Wrong data and inconsistent shapes fail."""
    with pytest.raises(ValueError):
        CrossSpectra.from_data(_noise(2000), DELTA, 2.0, duration=2.0)
    with pytest.raises(ValueError):
        CrossSpectra([1.0], [np.zeros((1, 2, 2))], [[0.0, 1.0]], [[1]])
    with pytest.raises(ValueError):
        CrossSpectra([1.0, 2.0], [np.zeros((1, 2, 2))], [[0.0]], [[1]])
    with pytest.raises(ValueError):
        CrossSpectra([1.0], [np.zeros((1, 2, 3))], [[0.0]], [[1]])
    with pytest.raises(ValueError):
        CrossSpectra([1.0], [np.zeros((1, 2, 2))], [[0.0]], [[0]])
