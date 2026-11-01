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
Unit tests for :mod:`shakelab.signals.fourier`.

The tests check the frequency axes, the scaling of the transforms
against analytical amplitudes, the Butterworth responses against their
defining formula and the causality of the minimum-phase version, the
regularized inversion and the smoothing kernels.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.signals import fourier


# ----------------------------------------------------------------------------
# Axes and transforms
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("nsamp, expected", [(10, 6), (11, 6), (1, 1)])
def test_rfft_length(nsamp, expected):
    """The one-sided length is nsamp // 2 + 1."""
    assert fourier.rfft_length(nsamp) == expected
    assert fourier.frequency_axis(0.1, nsamp).size == expected


def test_axes():
    """Frequency axis, inverse length and frequency ranges."""
    np.testing.assert_allclose(
        fourier.frequency_axis(0.1, 10), [0, 1, 2, 3, 4, 5]
    )
    assert fourier.irfft_length(0.01, 1.0) == 100
    np.testing.assert_allclose(
        fourier.frequency_range(1, 100, 3), [1, 10, 100]
    )
    np.testing.assert_allclose(
        fourier.frequency_range(1, 3, 3, log=False), [1, 2, 3]
    )


def test_transform_scaling_and_round_trip():
    """
    X(f) = delta * DFT approximates the continuous transform.

    A sine of unit amplitude lasting T seconds has |X| = T / 2 at its
    frequency.
    """
    delta = 0.01
    t = np.arange(1000) * delta
    x = np.sin(2.0 * np.pi * 5.0 * t)

    spectrum = fourier.fft(x, delta)
    k = np.argmin(np.abs(fourier.frequency_axis(delta, x.size) - 5.0))

    assert np.abs(spectrum[k]) == pytest.approx(5.0)
    np.testing.assert_allclose(
        fourier.ifft(spectrum, delta, x.size), x, atol=1e-12
    )


def test_shift_time():
    """An integer shift moves a pulse; the shift is circular."""
    x = np.zeros(16)
    x[3] = 1.0

    np.testing.assert_allclose(
        fourier.shift_time(x, 1.0, 2.0), np.roll(x, 2), atol=1e-12
    )
    np.testing.assert_allclose(
        fourier.shift_time(x, 1.0, -5.0), np.roll(x, -5), atol=1e-12
    )


# ----------------------------------------------------------------------------
# Butterworth
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("order", [1, 2, 3, 4, 8])
def test_butterworth_zero_phase_amplitude(order):
    """
    The zero-phase response is 1 / sqrt(1 + (f / fc) ** (2 n)).

    The amplitude is 1 / sqrt(2) at the corner for every order, odd
    orders included, and the response is real.
    """
    freq = np.array([0.0, 0.5, 1.0, 2.0, 10.0])
    lowpass = fourier.butterworth(freq, 1.0, order)
    highpass = fourier.butterworth(freq, 1.0, order, btype="highpass")

    np.testing.assert_allclose(
        lowpass, 1.0 / np.sqrt(1.0 + freq ** (2 * order))
    )
    np.testing.assert_allclose(highpass[1:], np.sqrt(
        1.0 - np.abs(lowpass[1:]) ** 2
    ))
    assert highpass[0] == 0.0
    assert np.abs(lowpass[2]) == pytest.approx(1.0 / np.sqrt(2.0))


@pytest.mark.parametrize("btype", ["lowpass", "highpass"])
def test_butterworth_minphase_same_amplitude(btype):
    """The causal response has the same amplitude as the zero-phase."""
    freq = np.linspace(0.0, 20.0, 201)
    causal = fourier.butterworth(freq, 2.0, 3, btype, minphase=True)
    zero_phase = fourier.butterworth(freq, 2.0, 3, btype)

    np.testing.assert_allclose(np.abs(causal), np.abs(zero_phase),
                               atol=1e-12)
    assert np.any(np.abs(np.angle(causal[1:])) > 0.1)


def test_butterworth_minphase_is_causal():
    """
    The causal filter leaves no energy before an impulse.

    The zero-phase filter instead gives a response symmetric around the
    impulse.
    """
    delta = 0.01
    pulse = np.zeros(1024)
    pulse[100] = 1.0 / delta
    freq = fourier.frequency_axis(delta, pulse.size)
    spectrum = fourier.fft(pulse, delta)

    causal = fourier.ifft(
        spectrum * fourier.butterworth(freq, 5.0, 4, minphase=True),
        delta, pulse.size,
    )
    symmetric = fourier.ifft(
        spectrum * fourier.butterworth(freq, 5.0, 4), delta, pulse.size
    )

    assert np.max(np.abs(causal[:100])) < 1e-3 * np.max(np.abs(causal))
    np.testing.assert_allclose(
        symmetric[80:100], symmetric[101:121][::-1], atol=1e-9
    )


def test_butterworth_errors():
    """Invalid corners, orders and types are rejected."""
    with pytest.raises(ValueError):
        fourier.butterworth([1.0], 0.0)
    with pytest.raises(ValueError):
        fourier.butterworth([1.0], 1.0, order=0)
    with pytest.raises(ValueError):
        fourier.butterworth([1.0], 1.0, order=2.5)
    with pytest.raises(ValueError):
        fourier.butterworth([1.0], 1.0, btype="bandpass")


# ----------------------------------------------------------------------------
# Inversion
# ----------------------------------------------------------------------------


def test_inverse_spectrum():
    """Values above the water level are inverted, the rest zeroed."""
    spectrum = np.array([2.0, 1e-9, 1j])

    sharp = fourier.inverse_spectrum(spectrum, 120, "sharp")
    np.testing.assert_allclose(sharp, [0.5, 0.0, -1j])

    smooth = fourier.inverse_spectrum(spectrum, 120, "smooth")
    np.testing.assert_allclose(smooth, [0.5, 0.0, -1j], rtol=1e-5)

    with pytest.raises(ValueError):
        fourier.inverse_spectrum(spectrum, 120, "other")


# ----------------------------------------------------------------------------
# Smoothing
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["gaussian", "konno-ohmachi"])
def test_smooth_preserves_constants(method):
    """Normalized kernels leave a constant unchanged."""
    freq = np.linspace(0.0, 10.0, 200)
    values = np.full(freq.size, 3.0)

    np.testing.assert_allclose(fourier.smooth(freq, values, method), 3.0)


@pytest.mark.parametrize("method", ["gaussian", "konno-ohmachi"])
def test_smooth_reduces_noise_and_keeps_zero_frequency(method):
    """Smoothing reduces scatter; the zero-frequency value is kept."""
    rng = np.random.default_rng(0)
    freq = np.linspace(0.0, 50.0, 2001)
    values = 1.0 + 0.5 * rng.standard_normal(freq.size)

    smoothed = fourier.smooth(freq, values, method)

    assert smoothed[0] == values[0]
    assert np.std(smoothed[200:]) < 0.25 * np.std(values[200:])


def test_konno_ohmachi_window():
    """
    The Konno-Ohmachi window follows its definition.

    Two values at f = fc and at b * log10(f / fc) = pi / 2 have weights
    1 and (2 / pi) ** 4.
    """
    b = 40.0
    fc = 1.0
    f2 = fc * 10.0 ** (np.pi / 2.0 / b)
    freq = np.array([fc, f2])
    values = np.array([1.0, 0.0])

    smoothed = fourier.smooth(
        freq, values, "konno-ohmachi", center=[fc], bandwidth=b
    )
    weight = (2.0 / np.pi) ** 4

    assert smoothed[0] == pytest.approx(1.0 / (1.0 + weight))


def test_smooth_at_center_frequencies():
    """Smoothed values can be computed at chosen frequencies."""
    freq = np.linspace(0.0, 10.0, 501)
    smoothed = fourier.smooth(
        freq, freq ** 0, "konno-ohmachi", center=[0.0, 1.0, 5.0]
    )

    assert np.isnan(smoothed[0])
    np.testing.assert_allclose(smoothed[1:], 1.0)


def test_smooth_errors():
    """Unknown methods and invalid parameters are rejected."""
    freq = np.linspace(0.0, 1.0, 5)

    with pytest.raises(ValueError):
        fourier.smooth(freq, freq, "boxcar")
    with pytest.raises(ValueError):
        fourier.smooth(freq, freq, sigma=0.0)
    with pytest.raises(ValueError):
        fourier.smooth(freq, freq[:3])


def test_smooth_amplitude_keeps_phase():
    """Amplitude smoothing leaves the phase unchanged."""
    rng = np.random.default_rng(1)
    freq = np.linspace(0.0, 10.0, 300)
    spectrum = (1.0 + 0.3 * rng.random(300)) * np.exp(
        1j * rng.uniform(-np.pi, np.pi, 300)
    )

    for log in (True, False):
        smoothed = fourier.smooth_amplitude(freq, spectrum, log=log)
        np.testing.assert_allclose(
            np.angle(smoothed), np.angle(spectrum), atol=1e-12
        )


def test_smooth_amplitude_phase_unwraps_phase():
    """
    A linear phase (a pure delay) survives smoothing.

    The wrapped phase jumps by 2 pi; averaging it without unwrapping
    would distort the delay.
    """
    freq = np.linspace(0.0, 10.0, 1001)
    spectrum = 2.0 * np.exp(-2j * np.pi * freq * 0.7)

    smoothed = fourier.smooth_amplitude_phase(freq, spectrum, sigma=0.01)

    np.testing.assert_allclose(np.abs(smoothed), 2.0)
    inner = slice(100, 900)
    np.testing.assert_allclose(
        np.unwrap(np.angle(smoothed))[inner],
        np.unwrap(np.angle(spectrum))[inner],
        atol=1e-2,
    )
