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
Division of time series into windows.

Spectral estimates averaged over many windows (Welch spectra, spectral
ratios, cross-spectral matrices for array processing, FDD, SPAC) all
divide a series into windows of equal length, possibly overlapping, and
often group consecutive windows into blocks. This module collects these
operations in one place:

- :func:`window_size`: number of samples of a window given as a
  duration or as a number of cycles of a frequency;
- :func:`window_starts`: first index of every complete window, or of
  every complete block of windows;
- :func:`segments`: the windows of an array, as a view without copies;
- :func:`taper`: a taper, optionally normalized;
- :func:`effective_realizations`: the number of independent windows
  equivalent to a set of overlapping tapered windows;
- :func:`split_record`: the windows of a
  :class:`~shakelab.signals.base.Record`, as new records.

Conventions
-----------
Only complete windows are returned: a window is never padded and the
samples after the last complete window are left out. The step between
consecutive windows is given either in samples (``step``) or as a
fractional overlap (``overlap``), in which case it is
``floor(size * (1 - overlap))``, at least one sample.

References
----------
Welch, P. D. (1967). The use of fast Fourier transform for the
estimation of power spectra: a method based on time averaging over
short, modified periodograms. IEEE Transactions on Audio and
Electroacoustics, 15(2), 70-73.

Percival, D. B., and Walden, A. T. (1993). Spectral Analysis for
Physical Applications. Cambridge University Press, Section 6.17.
"""

from __future__ import annotations

from numbers import Integral, Real
from typing import Any, Optional

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import signal

from shakelab.signals.base import Record

__all__ = [
    "effective_realizations",
    "segments",
    "split_record",
    "taper",
    "window_size",
    "window_starts",
]

# Relative tolerance below which a duration is taken as an exact
# multiple of the sampling interval (5.0 / 0.01 is 499.99999999999994).
_SIZE_TOLERANCE = 1e-9

# Normalizations accepted by :func:`taper`.
_TAPER_NORMS = (None, "rms", "sum")


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _check_positive_real(value: Any, name: str) -> float:
    """Return a finite, strictly positive real number."""
    if (
        not isinstance(value, Real)
        or isinstance(value, (bool, np.bool_))
    ):
        raise TypeError(f"{name} must be a real number")

    value = float(value)

    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive")

    return value


def _check_integer(value: Any, name: str, minimum: int) -> int:
    """Return an integer not smaller than ``minimum``."""
    if (
        not isinstance(value, Integral)
        or isinstance(value, (bool, np.bool_))
    ):
        raise TypeError(f"{name} must be an integer")

    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")

    return int(value)


def _resolve_step(size: int, overlap: Any, step: Any) -> int:
    """Return the step between windows from ``overlap`` or ``step``."""
    if step is not None:
        if overlap is not None:
            raise ValueError("give either overlap or step, not both")

        return _check_integer(step, "step", 1)

    if overlap is None:
        overlap = 0.0

    if (
        not isinstance(overlap, Real)
        or isinstance(overlap, (bool, np.bool_))
    ):
        raise TypeError("overlap must be a real number")

    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must be in the range [0, 1)")

    return max(1, int(np.floor(size * (1.0 - float(overlap)))))


# ----------------------------------------------------------------------------
# Window layout
# ----------------------------------------------------------------------------


def window_size(
    delta: float,
    duration: Optional[float] = None,
    cycles: Optional[float] = None,
    frequency: Optional[float] = None,
) -> int:
    """
    Return the number of samples of a window.

    The window is given either by its duration or by a number of cycles
    of a frequency (duration ``cycles / frequency``), as for the
    frequency-dependent windows of array processing.

    Parameters
    ----------
    delta : float
        Sampling interval in seconds.
    duration : float, optional
        Window duration in seconds.
    cycles : float, optional
        Number of cycles of ``frequency`` in the window.
    frequency : float, optional
        Frequency in Hz, required with ``cycles``.

    Returns
    -------
    int
        Number of samples: the duration divided by ``delta``, rounded
        down (durations that are multiples of ``delta`` up to rounding
        errors give the exact multiple).

    Raises
    ------
    TypeError
        If an argument is not a real number.
    ValueError
        If neither or both of ``duration`` and ``cycles`` are given, if
        ``frequency`` is missing with ``cycles`` or given without it, if
        a value is not positive, or if the window has fewer than two
        samples.

    Examples
    --------
    >>> window_size(0.01, duration=5.0)
    500
    >>> window_size(0.01, cycles=20, frequency=4.0)
    500
    """
    delta = _check_positive_real(delta, "delta")

    if (duration is None) == (cycles is None):
        raise ValueError("give either duration or cycles")

    if cycles is not None:
        if frequency is None:
            raise ValueError("frequency is required with cycles")

        cycles = _check_positive_real(cycles, "cycles")
        frequency = _check_positive_real(frequency, "frequency")
        duration = cycles / frequency

    else:
        if frequency is not None:
            raise ValueError("frequency is used only with cycles")

        duration = _check_positive_real(duration, "duration")

    ratio = duration / delta
    nearest = np.round(ratio)

    if abs(ratio - nearest) <= _SIZE_TOLERANCE * max(1.0, ratio):
        size = int(nearest)
    else:
        size = int(np.floor(ratio))

    if size < 2:
        raise ValueError(
            f"a window of {duration:g} s contains fewer than two samples"
        )

    return size


def window_starts(
    nsamp: int,
    size: int,
    overlap: Optional[float] = None,
    step: Optional[int] = None,
) -> np.ndarray:
    """
    Return the first index of every complete window.

    The same function gives the windows of a series of samples and the
    blocks of a series of windows (``nsamp`` windows grouped ``size``
    at a time).

    Parameters
    ----------
    nsamp : int
        Number of elements of the series.
    size : int
        Number of elements of each window.
    overlap : float, optional
        Fractional overlap between consecutive windows, in ``[0, 1)``.
        Default is no overlap.
    step : int, optional
        Step between consecutive windows, in elements; alternative to
        ``overlap``.

    Returns
    -------
    numpy.ndarray
        Start indices (int64), empty if ``size`` exceeds ``nsamp``.

    Raises
    ------
    TypeError
        If an argument has the wrong type.
    ValueError
        If ``nsamp`` is negative, ``size`` or ``step`` is smaller than
        one, ``overlap`` is outside ``[0, 1)``, or both ``overlap`` and
        ``step`` are given.

    Examples
    --------
    >>> window_starts(10, 4, overlap=0.5)
    array([0, 2, 4, 6])
    """
    nsamp = _check_integer(nsamp, "nsamp", 0)
    size = _check_integer(size, "size", 1)
    step = _resolve_step(size, overlap, step)

    if size > nsamp:
        return np.empty(0, dtype=np.int64)

    return np.arange(0, nsamp - size + 1, step, dtype=np.int64)


def segments(
    data: Any,
    size: int,
    overlap: Optional[float] = None,
    step: Optional[int] = None,
    axis: int = -1,
) -> np.ndarray:
    """
    Return the complete windows of an array.

    Parameters
    ----------
    data : array-like
        Input array; the windows are taken along ``axis``.
    size : int
        Number of samples of each window.
    overlap : float, optional
        Fractional overlap between consecutive windows, in ``[0, 1)``.
        Default is no overlap.
    step : int, optional
        Step between consecutive windows, in samples; alternative to
        ``overlap``.
    axis : int, optional
        Axis along which the windows are taken. Default is the last.

    Returns
    -------
    numpy.ndarray
        Read-only view with shape ``(n_windows,) + shape + (size,)``,
        where ``shape`` is the shape of ``data`` without ``axis``. The
        window ``i`` starts at sample ``window_starts(...)[i]``. When no
        window fits, ``n_windows`` is zero.

    Raises
    ------
    TypeError
        If an argument has the wrong type.
    ValueError
        If ``data`` has no dimensions, ``axis`` is out of bounds, or the
        window arguments are invalid (see :func:`window_starts`).

    Notes
    -----
    No data are copied: the result shares memory with ``data`` (with
    its array conversion, for other inputs) and must be copied before
    being modified.
    """
    data = np.asarray(data)

    if data.ndim == 0:
        raise ValueError("data must have at least one dimension")

    if (
        not isinstance(axis, Integral)
        or isinstance(axis, (bool, np.bool_))
    ):
        raise TypeError("axis must be an integer")

    if not -data.ndim <= axis < data.ndim:
        raise ValueError(
            f"axis {axis} is out of bounds for an array of dimension "
            f"{data.ndim}"
        )

    axis = int(axis) % data.ndim
    size = _check_integer(size, "size", 1)
    step = _resolve_step(size, overlap, step)

    if size > data.shape[axis]:
        shape = (0,) + data.shape[:axis] + data.shape[axis + 1:]
        return np.empty(shape + (size,), dtype=data.dtype)

    # sliding_window_view puts the window index at ``axis`` and the
    # samples of each window last.
    view = sliding_window_view(data, size, axis=axis)
    index = [slice(None)] * view.ndim
    index[axis] = slice(None, None, step)

    return np.moveaxis(view[tuple(index)], axis, 0)


# ----------------------------------------------------------------------------
# Tapers
# ----------------------------------------------------------------------------


def taper(
    size: int,
    window: Any = "hann",
    norm: Optional[str] = None,
) -> np.ndarray:
    """
    Return a taper for windows of a given length.

    Parameters
    ----------
    size : int
        Number of samples.
    window : str, tuple or array-like, optional
        Window specification of ``scipy.signal.get_window`` (name, or
        tuple of name and parameters, e.g. ``("tukey", 0.1)``) or an
        array of ``size`` values. Default is ``"hann"``. Named windows
        are periodic (``fftbins=True``), as appropriate for spectral
        analysis.
    norm : {None, "rms", "sum"}, optional
        ``None`` (default): no normalization. ``"rms"``: unit mean
        square value, so that tapering preserves the power of a
        stationary signal. ``"sum"``: unit sum of the values, so that a
        tapered sinusoid of amplitude ``A`` has a DFT of modulus
        ``A / 2`` at its frequency.

    Returns
    -------
    numpy.ndarray
        Taper values.

    Raises
    ------
    TypeError
        If ``size`` is not an integer.
    ValueError
        If ``size`` is smaller than two, an array window has the wrong
        length or is not finite, ``norm`` is unknown, or the taper has
        zero energy.
    """
    size = _check_integer(size, "size", 2)

    if norm not in _TAPER_NORMS:
        names = ", ".join(repr(value) for value in _TAPER_NORMS)
        raise ValueError(f"norm must be one of {names}")

    if isinstance(window, (str, tuple)):
        values = signal.get_window(window, size, fftbins=True)
    else:
        values = np.asarray(window, dtype=float)

        if values.shape != (size,):
            raise ValueError(
                f"window must contain {size} values, got shape "
                f"{values.shape}"
            )

    values = np.array(values, dtype=float)

    if not np.all(np.isfinite(values)):
        raise ValueError("window values must be finite")

    if norm == "rms":
        scale = np.sqrt(np.mean(values ** 2))
    elif norm == "sum":
        scale = np.sum(values)
    else:
        return values

    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("the taper cannot be normalized")

    return values / scale


def effective_realizations(
    taper_values: Any,
    step: int,
    count: int,
) -> float:
    """
    Return the effective number of independent windows.

    Overlapping tapered windows of a stationary series are correlated,
    so that the average of ``count`` periodograms has a larger variance
    than the average of ``count`` independent ones. The effective
    number is the number of independent windows with the same
    variance (Percival and Walden, 1993, Section 6.17):

    ``count / (1 + 2 * sum_q (1 - q / count) * rho_q**2)``,

    where ``rho_q = sum_n w[n] w[n + q * step] / sum_n w[n]**2`` is the
    overlap of the taper ``w`` with itself shifted by ``q`` windows.

    Parameters
    ----------
    taper_values : array-like
        Taper of the windows (any normalization).
    step : int
        Step between consecutive windows, in samples.
    count : int
        Number of windows averaged.

    Returns
    -------
    float
        Effective number of independent windows, between one and
        ``count``.

    Raises
    ------
    TypeError
        If ``step`` or ``count`` is not an integer.
    ValueError
        If the taper is not one-dimensional, has fewer than two values,
        is not finite or has zero energy, or if ``step`` or ``count``
        is smaller than one.

    Notes
    -----
    The formula assumes a series whose spectrum varies little over the
    bandwidth of the taper (white noise in the derivation); it accounts
    only for the correlation introduced by the overlap.
    """
    values = np.asarray(taper_values, dtype=float)

    if values.ndim != 1 or values.size < 2:
        raise ValueError(
            "taper_values must be one-dimensional with at least two "
            "values"
        )

    if not np.all(np.isfinite(values)):
        raise ValueError("taper_values must be finite")

    step = _check_integer(step, "step", 1)
    count = _check_integer(count, "count", 1)

    energy = float(np.dot(values, values))

    if energy <= 0.0:
        raise ValueError("the taper has zero energy")

    denominator = 1.0

    for lag in range(1, count):
        shift = lag * step

        if shift >= values.size:
            break

        rho = float(np.dot(values[:-shift], values[shift:])) / energy
        denominator += 2.0 * (1.0 - lag / count) * rho ** 2

    return count / denominator


# ----------------------------------------------------------------------------
# Records
# ----------------------------------------------------------------------------


def split_record(
    record: Record,
    duration: float,
    overlap: Optional[float] = None,
) -> list[Record]:
    """
    Split a record into windows of equal duration.

    Parameters
    ----------
    record : Record
        Record to split; it is not modified.
    duration : float
        Duration of each window, in seconds (see :func:`window_size`).
    overlap : float, optional
        Fractional overlap between consecutive windows, in ``[0, 1)``.
        Default is no overlap.

    Returns
    -------
    list of Record
        Independent records, one per complete window, with the header
        of ``record`` and the start time of their first sample. Empty if
        the record is shorter than one window.

    Raises
    ------
    TypeError
        If ``record`` is not a Record.
    ValueError
        If the record has no sampling interval, or the window arguments
        are invalid.
    """
    if not isinstance(record, Record):
        raise TypeError("record must be a Record")

    if record.delta is None:
        raise ValueError("the record has no sampling interval")

    size = window_size(record.delta, duration=duration)
    starts = window_starts(record.nsamp, size, overlap=overlap)
    pieces = []

    for start in starts:
        piece = Record()
        piece.head = record.head.copy(parent=piece)
        piece.data = record.data[start:start + size]

        if record.time is not None:
            piece.time = record.time + int(start) * record.delta

        pieces.append(piece)

    return pieces
