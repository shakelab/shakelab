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
Fourier-domain primitives.

This module collects the functions on which the spectral objects of
:mod:`shakelab.signals` are built: frequency axes, the scaled discrete
Fourier transform, time shifts, Butterworth transfer functions,
regularized inversion and spectral smoothing. All functions work on
NumPy arrays and do not depend on other ShakeLab classes.

Conventions
-----------
- Spectra are one-sided (``numpy.fft.rfft`` layout) for real signals
  of ``nsamp`` samples at interval ``delta``: frequencies
  ``k / (nsamp * delta)``, ``k = 0, ..., rfft_length(nsamp) - 1``.
- The forward transform is scaled by the sampling interval,
  ``X(f) = delta * DFT(x)``, approximating the continuous Fourier
  transform: a signal in m/s gives a spectrum in m/s * s. The inverse
  divides by ``delta``, so the two are exact inverses.
- The kernel is ``exp(-2j * pi * f * t)``: a causal filter with
  Laplace transfer function ``H(s)`` acts through ``H(2j * pi * f)``.
- Frequencies are in Hz.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import signal as _signal

__all__ = [
    "rfft_length",
    "irfft_length",
    "frequency_axis",
    "frequency_range",
    "fft",
    "ifft",
    "shift_time",
    "butterworth",
    "inverse_spectrum",
    "smooth",
    "smooth_amplitude",
    "smooth_amplitude_phase",
]


# ----------------------------------------------------------------------------
# Axes and lengths
# ----------------------------------------------------------------------------


def rfft_length(nsamp: int) -> int:
    """
    Return the number of frequencies of the one-sided spectrum.

    Parameters
    ----------
    nsamp
        Number of samples of the time series.

    Returns
    -------
    int
        ``nsamp // 2 + 1``.

    Examples
    --------
    >>> rfft_length(10), rfft_length(11)
    (6, 6)
    """
    return int(nsamp) // 2 + 1


def irfft_length(delta: float, dfreq: float) -> int:
    """
    Return the number of samples for a given time and frequency step.

    Parameters
    ----------
    delta
        Sampling interval in seconds.
    dfreq
        Frequency step in Hz.

    Returns
    -------
    int
        ``round(1 / (dfreq * delta))``.

    Examples
    --------
    >>> irfft_length(0.01, 1.0)
    100
    """
    return int(round(1.0 / (dfreq * delta)))


def frequency_axis(delta: float, nsamp: int) -> np.ndarray:
    """
    Return the frequencies of the one-sided spectrum.

    Parameters
    ----------
    delta
        Sampling interval in seconds.
    nsamp
        Number of samples of the time series.

    Returns
    -------
    numpy.ndarray
        Frequencies in Hz, from 0 to the Nyquist frequency (included
        only for even ``nsamp``).
    """
    return np.fft.rfftfreq(int(nsamp), delta)


def frequency_range(
    fmin: float,
    fmax: float,
    fnum: int,
    log: bool = True,
) -> np.ndarray:
    """
    Return a linearly or logarithmically spaced frequency range.

    Parameters
    ----------
    fmin, fmax
        First and last frequency in Hz (both strictly positive for a
        logarithmic range).
    fnum
        Number of frequencies.
    log
        Logarithmic spacing (default) or linear spacing.

    Returns
    -------
    numpy.ndarray
        Frequencies in Hz.
    """
    if log:
        return np.logspace(np.log10(fmin), np.log10(fmax), fnum)

    return np.linspace(fmin, fmax, fnum)


# ----------------------------------------------------------------------------
# Transforms
# ----------------------------------------------------------------------------


def fft(
    data: np.ndarray,
    delta: float,
    nsamp: Optional[int] = None,
) -> np.ndarray:
    """
    Return the scaled one-sided Fourier transform of a real signal.

    Parameters
    ----------
    data
        Real samples.
    delta
        Sampling interval in seconds.
    nsamp
        Transform length; the signal is zero-padded or truncated to it.
        Default is the number of samples.

    Returns
    -------
    numpy.ndarray
        ``delta * numpy.fft.rfft(data, nsamp)``.
    """
    return np.fft.rfft(data, nsamp) * delta


def ifft(spectrum: np.ndarray, delta: float, nsamp: int) -> np.ndarray:
    """
    Return the real signal of a scaled one-sided spectrum.

    Parameters
    ----------
    spectrum
        One-sided spectrum, as returned by :func:`fft`.
    delta
        Sampling interval in seconds.
    nsamp
        Number of samples of the signal (needed to distinguish even and
        odd lengths).

    Returns
    -------
    numpy.ndarray
        ``numpy.fft.irfft(spectrum, nsamp) / delta``.
    """
    return np.fft.irfft(spectrum, int(nsamp)) / delta


def shift_time(data: np.ndarray, delta: float, shift: float) -> np.ndarray:
    """
    Shift a real signal in time by a phase rotation of its spectrum.

    Parameters
    ----------
    data
        Real samples.
    delta
        Sampling interval in seconds.
    shift
        Time shift in seconds (positive delays the signal); fractional
        samples are allowed.

    Returns
    -------
    numpy.ndarray
        Shifted samples. The shift is circular: samples moved past one
        end reappear at the other (pad with zeros to avoid it).
    """
    nsamp = len(data)
    freq = frequency_axis(delta, nsamp)
    rotation = np.exp(-2j * np.pi * shift * freq)

    return np.fft.irfft(np.fft.rfft(data) * rotation, nsamp)


# ----------------------------------------------------------------------------
# Filters and inversion
# ----------------------------------------------------------------------------


def butterworth(
    freq: np.ndarray,
    corner: float,
    order: int = 4,
    btype: str = "lowpass",
    minphase: bool = False,
) -> np.ndarray:
    """
    Return the transfer function of a Butterworth filter.

    Parameters
    ----------
    freq
        Frequencies in Hz (non-negative).
    corner
        Corner frequency in Hz, where the amplitude is ``1 / sqrt(2)``
        (-3 dB).
    order
        Filter order (strictly positive).
    btype
        ``"lowpass"`` or ``"highpass"``.
    minphase
        If False (default), return the zero-phase response, i.e. the
        real amplitude ``1 / sqrt(1 + (f / corner) ** (2 * order))``
        (low-pass) or ``1 / sqrt(1 + (corner / f) ** (2 * order))``
        (high-pass): the filtered signal is not delayed. If True,
        return the complex response of the causal analog Butterworth
        filter, which has the same amplitude and minimum phase, as a
        physical filter.

    Returns
    -------
    numpy.ndarray
        Complex transfer function at ``freq``.

    Raises
    ------
    ValueError
        If the corner, the order or the filter type are invalid.

    References
    ----------
    Butterworth, S. (1930). On the theory of filter amplifiers.
    Wireless Engineer, 7, 536-541.
    """
    if not np.isfinite(corner) or corner <= 0.0:
        raise ValueError("corner must be finite and strictly positive")

    if int(order) != order or order <= 0:
        raise ValueError("order must be a strictly positive integer")

    if btype not in ("lowpass", "highpass"):
        raise ValueError("btype must be 'lowpass' or 'highpass'")

    freq = np.asarray(freq, dtype=float)
    order = int(order)

    if minphase:
        b, a = _signal.butter(
            order, 2.0 * np.pi * corner, btype=btype, analog=True
        )
        _, response = _signal.freqs(b, a, worN=2.0 * np.pi * freq)
        return response.astype(complex)

    if btype == "lowpass":
        ratio = freq / corner
        amplitude = 1.0 / np.sqrt(1.0 + ratio ** (2 * order))
    else:
        amplitude = np.zeros(freq.shape)
        positive = freq > 0.0
        ratio = corner / freq[positive]
        amplitude[positive] = 1.0 / np.sqrt(1.0 + ratio ** (2 * order))

    return amplitude.astype(complex)


def inverse_spectrum(
    spectrum: np.ndarray,
    waterlevel: float = 120.0,
    method: str = "smooth",
) -> np.ndarray:
    """
    Return a regularized inverse of a spectrum.

    Parameters
    ----------
    spectrum
        Complex spectrum ``S``.
    waterlevel
        Level in dB defining the threshold ``w = 10 ** (-waterlevel /
        20)``.
    method
        ``"sharp"``: ``1 / S`` where ``|S| >= w``, zero elsewhere;
        ``"smooth"``: ``conj(S) / (|S| ** 2 + w)``, set to zero where
        ``|S| <= w``.

    Returns
    -------
    numpy.ndarray
        Inverse spectrum.

    Raises
    ------
    ValueError
        If the method is unknown.

    Notes
    -----
    The threshold is absolute, so its effect depends on the units of
    the spectrum; its definition is to be reviewed together with the
    instrument response correction.
    """
    spectrum = np.asarray(spectrum, dtype=complex)
    amplitude = np.abs(spectrum)
    level = 10.0 ** (-waterlevel / 20.0)

    if method == "sharp":
        inverse = np.zeros(spectrum.shape, dtype=complex)
        valid = amplitude >= level
        inverse[valid] = 1.0 / spectrum[valid]
        return inverse

    if method == "smooth":
        inverse = spectrum.conj() / (amplitude ** 2 + level)
        inverse[amplitude <= level] = 0.0
        return inverse

    raise ValueError("method must be 'sharp' or 'smooth'")


# ----------------------------------------------------------------------------
# Smoothing
# ----------------------------------------------------------------------------


def smooth(
    freq: np.ndarray,
    values: np.ndarray,
    method: str = "gaussian",
    center: Optional[np.ndarray] = None,
    sigma: float = 0.2,
    truncation: float = 3.0,
    bandwidth: float = 40.0,
) -> np.ndarray:
    """
    Smooth real values defined on a frequency axis.

    The smoothed value at a center frequency ``fc`` is the weighted
    mean of the values at the positive frequencies ``f``, with weights
    given by a kernel of ``f / fc``.

    Parameters
    ----------
    freq
        Frequencies in Hz.
    values
        Real values at ``freq``.
    method
        Smoothing kernel:

        - ``"gaussian"``: Gaussian in natural-log frequency,
          ``exp(-(ln(f / fc) / sigma) ** 2 / 2)``, truncated at
          ``truncation * sigma``;
        - ``"konno-ohmachi"``: the window of Konno and Ohmachi (1998),
          ``(sin(b * log10(f / fc)) / (b * log10(f / fc))) ** 4`` with
          ``b = bandwidth``.
    center
        Frequencies at which the smoothed values are computed. Default
        is ``freq``.
    sigma
        Standard deviation of the Gaussian kernel, in natural-log
        frequency units.
    truncation
        Truncation of the Gaussian kernel, in units of ``sigma``.
    bandwidth
        Bandwidth coefficient ``b`` of the Konno-Ohmachi window (larger
        values give a narrower window).

    Returns
    -------
    numpy.ndarray
        Smoothed values at ``center``. Values at non-positive center
        frequencies are not smoothed: they are copied from ``values``
        when ``center`` is ``freq``, and are NaN otherwise.

    Raises
    ------
    ValueError
        If the method or the kernel parameters are invalid.

    References
    ----------
    Konno, K., and Ohmachi, T. (1998). Ground-motion characteristics
    estimated from spectral ratio between horizontal and vertical
    components of microtremor. Bulletin of the Seismological Society of
    America, 88(1), 228-241.
    """
    freq = np.asarray(freq, dtype=float)
    values = np.asarray(values, dtype=float)

    if freq.shape != values.shape or freq.ndim != 1:
        raise ValueError("freq and values must be 1-D arrays of equal size")

    if method == "gaussian":
        if sigma <= 0.0 or truncation <= 0.0:
            raise ValueError("sigma and truncation must be positive")
    elif method == "konno-ohmachi":
        if bandwidth <= 0.0:
            raise ValueError("bandwidth must be positive")
    else:
        raise ValueError("method must be 'gaussian' or 'konno-ohmachi'")

    same_axis = center is None
    center = freq if same_axis else np.asarray(center, dtype=float)

    positive = freq > 0.0
    log_freq = np.log(freq[positive])
    data = values[positive]

    output = np.full(center.shape, np.nan)
    targets = np.flatnonzero(center > 0.0)

    # Process the centers in blocks to bound the memory of the
    # (centers x frequencies) weight matrix.
    block = max(1, 2 ** 22 // max(log_freq.size, 1))

    for start in range(0, targets.size, block):
        index = targets[start:start + block]
        distance = log_freq[None, :] - np.log(center[index])[:, None]

        if method == "gaussian":
            weights = np.exp(-0.5 * (distance / sigma) ** 2)
            weights[np.abs(distance) > truncation * sigma] = 0.0
        else:
            x = bandwidth * distance / np.log(10.0)
            with np.errstate(invalid="ignore", divide="ignore"):
                weights = (np.sin(x) / x) ** 4
            weights[x == 0.0] = 1.0

        total = weights.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            output[index] = weights @ data / total

    if same_axis:
        output[~(center > 0.0)] = values[~(center > 0.0)]

    return output


def smooth_amplitude(
    freq: np.ndarray,
    spectrum: np.ndarray,
    method: str = "gaussian",
    log: bool = True,
    **options,
) -> np.ndarray:
    """
    Smooth the amplitude of a complex spectrum, keeping its phase.

    Parameters
    ----------
    freq
        Frequencies in Hz.
    spectrum
        Complex spectrum at ``freq``.
    method
        Smoothing kernel (see :func:`smooth`).
    log
        If True (default), average the logarithm of the amplitude
        (geometric mean); otherwise average the amplitude itself.
    **options
        Kernel parameters passed to :func:`smooth` (``sigma``,
        ``truncation``, ``bandwidth``).

    Returns
    -------
    numpy.ndarray
        Complex spectrum with smoothed amplitude and the original phase.
    """
    spectrum = np.asarray(spectrum, dtype=complex)
    amplitude = np.abs(spectrum)
    phase = np.angle(spectrum)

    if log:
        with np.errstate(divide="ignore"):
            smoothed = np.exp(
                smooth(freq, np.log(amplitude), method, **options)
            )
    else:
        smoothed = smooth(freq, amplitude, method, **options)

    return smoothed * np.exp(1j * phase)


def smooth_amplitude_phase(
    freq: np.ndarray,
    spectrum: np.ndarray,
    method: str = "gaussian",
    **options,
) -> np.ndarray:
    """
    Smooth the log-amplitude and the unwrapped phase of a spectrum.

    Parameters
    ----------
    freq
        Frequencies in Hz.
    spectrum
        Complex spectrum at ``freq``.
    method
        Smoothing kernel (see :func:`smooth`).
    **options
        Kernel parameters passed to :func:`smooth`.

    Returns
    -------
    numpy.ndarray
        Complex spectrum ``exp(smoothed log-amplitude + 1j * smoothed
        phase)``. The phase is unwrapped along frequency before
        smoothing, so that jumps of 2 pi do not bias the average.
    """
    spectrum = np.asarray(spectrum, dtype=complex)

    with np.errstate(divide="ignore"):
        log_amplitude = np.log(np.abs(spectrum))

    phase = np.unwrap(np.angle(spectrum))

    smoothed = (
        smooth(freq, log_amplitude, method, **options)
        + 1j * smooth(freq, phase, method, **options)
    )

    return np.exp(smoothed)
