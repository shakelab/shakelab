# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Fourier coefficients of windows and cross-spectral matrices.

Methods based on the coherence between channels (frequency-wavenumber
analysis of arrays, frequency domain decomposition, spatial
autocorrelation) estimate, at each frequency, the cross-spectral matrix
of the channels by averaging the products of the Fourier coefficients
of many windows (Welch's method extended to several channels). This
module provides the two steps:

- :func:`coefficients`: the Fourier coefficients of the windows of a
  set of channels at selected frequencies, computed directly at any
  frequency (``transform="direct"``) or from one real FFT per window
  (``transform="rfft"``, at the nearest bins or interpolated);
- :func:`cross_spectral_matrix`: the average of the outer products of
  the coefficients, optionally with the variance of each element.

Scaling
-------
With ``scaling="psd"`` (default) the coefficients are scaled so that
their squared modulus is the one-sided power spectral density of the
window (units of the data squared per Hz, as ``scaling="psd"`` of
:func:`shakelab.signals.spectrogram.wft_transform`), and the
cross-spectral matrix is the one-sided cross power spectral density
matrix. The data may be in any units, also uncorrected (counts): the
matrices are then in those units squared per Hz. With ``scaling=None``
the coefficients are the transform of the tapered window with the
convention of :mod:`shakelab.signals.fourier`
(``delta * DFT``).

The phase of the coefficients refers to the first sample of each
window; it cancels in the cross-spectral matrix, which depends only on
the phase differences between channels.

References
----------
Welch, P. D. (1967). The use of fast Fourier transform for the
estimation of power spectra: a method based on time averaging over
short, modified periodograms. IEEE Transactions on Audio and
Electroacoustics, 15(2), 70-73.

Bendat, J. S., and Piersol, A. G. (2010). Random Data: Analysis and
Measurement Procedures, 4th edition. Wiley, Chapter 11.
"""

from __future__ import annotations

from typing import Any, Optional, Union

import numpy as np

from shakelab.signals import fourier
from shakelab.signals.windowing import segments, taper, window_starts

__all__ = [
    "SCALINGS",
    "SELECTIONS",
    "TRANSFORMS",
    "coefficients",
    "cross_spectral_matrix",
]

# Accepted values of the arguments of :func:`coefficients`.
TRANSFORMS = ("direct", "rfft")
SELECTIONS = ("nearest", "interp")
SCALINGS = ("psd", None)

# Default zero padding of the RFFT, by frequency selection.
_DEFAULT_PADDING = {"nearest": 1, "interp": 8}

# Number of windows transformed at once by the RFFT (limits memory).
_RFFT_CHUNK = 64


# ----------------------------------------------------------------------------
# Fourier coefficients
# ----------------------------------------------------------------------------


def coefficients(
    data: Any,
    delta: float,
    freq: Any,
    size: int,
    overlap: Optional[float] = None,
    step: Optional[int] = None,
    window: Any = ("tukey", 0.1),
    transform: str = "direct",
    selection: str = "nearest",
    padding: Optional[int] = None,
    scaling: Optional[str] = "psd",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the Fourier coefficients of the windows of a set of channels.

    Parameters
    ----------
    data : array-like
        Real samples, time along the last axis: ``(nsamp,)`` for one
        channel, ``(n_channels, nsamp)`` or any shape ``(...,
        nsamp)``.
    delta : float
        Sampling interval in seconds.
    freq : float or array-like
        Frequencies in Hz, strictly between zero and the Nyquist
        frequency.
    size : int
        Number of samples of each window (see
        :func:`shakelab.signals.windowing.window_size`).
    overlap : float, optional
        Fractional overlap between consecutive windows, in ``[0, 1)``.
        Default is no overlap.
    step : int, optional
        Step between windows in samples; alternative to ``overlap``.
    window : str, tuple or array-like, optional
        Taper (see :func:`shakelab.signals.windowing.taper`). Default is
        a Tukey window with 10% tapered edges.
    transform : {"direct", "rfft"}, optional
        ``"direct"`` (default) evaluates the transform at the requested
        frequencies (:func:`shakelab.signals.fourier.dft`). ``"rfft"``
        computes one real FFT per window and takes the requested
        frequencies from it, as selected by ``selection``: faster for
        many frequencies.
    selection : {"nearest", "interp"}, optional
        With ``transform="rfft"``: ``"nearest"`` (default) returns the
        FFT bins nearest to the requested frequencies (repeated bins
        once), ``"interp"`` interpolates the complex spectrum at the
        requested frequencies with local cubic polynomials.
    padding : int, optional
        With ``transform="rfft"``, zero-padding factor of the windows
        (FFT length ``padding * size``): it refines the frequency grid,
        not the resolution. Default is 1 for ``"nearest"`` and 8 for
        ``"interp"``.
    scaling : {"psd", None}, optional
        ``"psd"`` (default): the squared modulus is the one-sided power
        spectral density. ``None``: ``delta * DFT`` of the tapered
        window.

    Returns
    -------
    freq : numpy.ndarray
        Frequencies of the coefficients: the requested ones, or the
        selected FFT bins with ``selection="nearest"``.
    values : numpy.ndarray
        Complex coefficients with shape ``(n_freq, n_windows) + shape``,
        where ``shape`` is the shape of ``data`` without the time axis.
        The windows start at the samples given by
        :func:`shakelab.signals.windowing.window_starts` with the same
        arguments.

    Raises
    ------
    ValueError
        If the data are not finite, an option is unknown, a frequency is
        outside ``(0, Nyquist)``, no window fits in the data, or the
        window arguments are invalid.

    Examples
    --------
    >>> freq, values = coefficients(data, 0.01, [2.0, 4.0], 500,
    ...                             overlap=0.5)
    >>> csm = cross_spectral_matrix(values[0])   # at 2 Hz
    """
    data = np.asarray(data, dtype=float)

    if data.ndim == 0 or not np.all(np.isfinite(data)):
        raise ValueError("data must be a finite array of samples")

    if transform not in TRANSFORMS:
        raise ValueError(f"transform must be one of {TRANSFORMS}")

    if selection not in SELECTIONS:
        raise ValueError(f"selection must be one of {SELECTIONS}")

    if scaling not in SCALINGS:
        raise ValueError(f"scaling must be one of {SCALINGS}")

    delta = float(delta)

    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError("delta must be positive")

    freq = np.atleast_1d(np.asarray(freq, dtype=float))
    nyquist = 0.5 / delta

    if (
        freq.ndim != 1
        or not np.all(np.isfinite(freq))
        or np.any(freq <= 0.0)
        or np.any(freq >= nyquist)
    ):
        raise ValueError(
            f"frequencies must be in (0, {nyquist:g}) Hz (Nyquist)"
        )

    starts = window_starts(data.shape[-1], size, overlap=overlap, step=step)

    if starts.size == 0:
        raise ValueError(
            f"no window of {size} samples fits in {data.shape[-1]} samples"
        )

    weights = taper(size, window)
    windows = segments(data, size, overlap=overlap, step=step)

    if transform == "direct":
        values = fourier.dft(windows * weights, delta, freq)
    else:
        freq, values = _rfft_coefficients(
            windows, weights, delta, freq, selection, padding
        )

    if scaling == "psd":
        values = values * np.sqrt(2.0 / (delta * np.dot(weights, weights)))

    # (n_windows, ..., n_freq) -> (n_freq, n_windows, ...)
    return freq, np.moveaxis(values, -1, 0)


def _rfft_coefficients(
    windows: np.ndarray,
    weights: np.ndarray,
    delta: float,
    freq: np.ndarray,
    selection: str,
    padding: Optional[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Coefficients from one real FFT per window, ``delta * DFT``."""
    if padding is None:
        padding = _DEFAULT_PADDING[selection]

    if (
        isinstance(padding, bool)
        or not isinstance(padding, (int, np.integer))
        or padding < 1
    ):
        raise ValueError("padding must be a positive integer")

    size = weights.size
    nfft = size * int(padding)
    bin_width = 1.0 / (nfft * delta)

    if selection == "nearest":
        index = np.rint(freq / bin_width).astype(np.int64)
        _, first = np.unique(index, return_index=True)
        index = index[np.sort(first)]
        freq = index * bin_width
    else:
        position = freq / bin_width

    nwin = windows.shape[0]
    shape = windows.shape[1:-1]
    values = np.empty((nwin,) + shape + (freq.size,), dtype=complex)

    for first in range(0, nwin, _RFFT_CHUNK):
        chunk = slice(first, first + _RFFT_CHUNK)
        spectrum = fourier.fft(windows[chunk] * weights, delta, nfft)

        if selection == "nearest":
            values[chunk] = spectrum[..., index]
        else:
            values[chunk] = _cubic_interpolation(spectrum, position)

    return freq, values


def _cubic_interpolation(spectrum: np.ndarray, position: np.ndarray) -> Any:
    """
    Interpolate a complex spectrum at fractional bin positions.

    Four-point Lagrange cubic polynomials in the interior, linear
    interpolation within one bin of the ends, exact values at the bins.
    The last axis of ``spectrum`` holds the bins.
    """
    last = spectrum.shape[-1] - 1
    base = np.clip(np.floor(position).astype(np.int64), 0, last - 1)
    x = position - base
    interior = (base >= 1) & (base + 2 <= last)
    result = np.empty(spectrum.shape[:-1] + position.shape, dtype=complex)

    if np.any(interior):
        b = base[interior]
        t = x[interior]
        weights = (
            -t * (t - 1.0) * (t - 2.0) / 6.0,
            (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0,
            -(t + 1.0) * t * (t - 2.0) / 2.0,
            (t + 1.0) * t * (t - 1.0) / 6.0,
        )
        result[..., interior] = sum(
            spectrum[..., b + offset] * weight
            for offset, weight in zip((-1, 0, 1, 2), weights)
        )

    edge = ~interior

    if np.any(edge):
        b = base[edge]
        t = x[edge]
        result[..., edge] = (
            spectrum[..., b] * (1.0 - t) + spectrum[..., b + 1] * t
        )

    return result


# ----------------------------------------------------------------------------
# Cross-spectral matrices
# ----------------------------------------------------------------------------


def cross_spectral_matrix(
    values: Any,
    other: Any = None,
    variance: bool = False,
) -> Union[np.ndarray, tuple[np.ndarray, np.ndarray]]:
    """
    Return the cross-spectral matrix averaged over realizations.

    ``C[i, j] = mean_m X[m, i] * conj(Y[m, j])``, with ``Y = X`` by
    default: a Hermitian matrix whose diagonal holds the power spectral
    densities of the channels (with coefficients scaled as
    ``scaling="psd"`` of :func:`coefficients`).

    Parameters
    ----------
    values : array-like
        Complex coefficients ``(n_realizations, n_channels)`` at one
        frequency (e.g. ``coefficients(...)[1][k]``).
    other : array-like, optional
        Coefficients ``(n_realizations, n_other)`` of a second set of
        channels, for the cross terms between two sets.
    variance : bool, optional
        If True, also return the sample variance of the products of
        each element, ``L / (L - 1) * (mean |x y*|**2 - |C|**2)``; the
        variance of the average is this value divided by the effective
        number of independent realizations (see
        :func:`shakelab.signals.windowing.effective_realizations`).

    Returns
    -------
    numpy.ndarray or tuple of numpy.ndarray
        Matrix ``(n_channels, n_other)``, and its element variance if
        ``variance`` is True.

    Raises
    ------
    ValueError
        If the coefficients are not two-dimensional, the two sets have
        different numbers of realizations, there are no realizations,
        or fewer than two with ``variance``.
    """
    x = np.asarray(values, dtype=complex)
    y = x if other is None else np.asarray(other, dtype=complex)

    if x.ndim != 2 or y.ndim != 2:
        raise ValueError(
            "coefficients must be (n_realizations, n_channels) arrays"
        )

    if x.shape[0] != y.shape[0]:
        raise ValueError("the two sets have different realizations")

    count = x.shape[0]

    if count == 0:
        raise ValueError("at least one realization is required")

    matrix = x.T @ y.conj() / count

    if not variance:
        return matrix

    if count < 2:
        raise ValueError("the variance needs at least two realizations")

    squares = (np.abs(x) ** 2).T @ (np.abs(y) ** 2) / count
    spread = count / (count - 1) * (squares - np.abs(matrix) ** 2)

    return matrix, np.maximum(spread, 0.0)
