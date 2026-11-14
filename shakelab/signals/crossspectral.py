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

The class :class:`CrossSpectra` combines the two steps for a set of
frequencies: it averages the matrices by blocks of consecutive windows
(their evolution in time), optionally also over sub-frequencies around
each frequency (a frequency smoothing that adds realizations), with
windows of fixed duration or of a fixed number of cycles of each
frequency, and provides the coherence matrices and the decomposition
in eigenvalues and eigenvectors used by the frequency domain
decomposition.

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

Brincker, R., Zhang, L., and Andersen, P. (2001). Modal identification
of output-only systems using frequency domain decomposition. Smart
Materials and Structures, 10(3), 441-445.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Optional, Union

import numpy as np

from shakelab.signals import fourier
from shakelab.signals.windowing import (
    segments,
    taper,
    window_size,
    window_starts,
)

__all__ = [
    "SCALINGS",
    "SELECTIONS",
    "TRANSFORMS",
    "CrossSpectra",
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


# ----------------------------------------------------------------------------
# Cross-spectral matrices by frequency and block
# ----------------------------------------------------------------------------


class CrossSpectra:
    """
    Cross-spectral matrices of a set of channels, by frequency and block.

    At each frequency the matrices average the outer products of the
    Fourier coefficients of blocks of consecutive windows (and of the
    sub-frequencies around the frequency, if any). With windows of a
    fixed number of cycles the window length, and so the number and the
    times of the blocks, change with the frequency: the matrices, times
    and realizations are therefore given frequency by frequency
    (:meth:`matrices`, :meth:`times`, :meth:`realizations`).

    Usually created with :meth:`from_data`.

    Parameters
    ----------
    frequency : array-like
        Frequencies in Hz, ``(n_freq,)``.
    matrices : sequence of array-like
        For each frequency, the Hermitian matrices of its blocks,
        ``(n_blocks, n_channels, n_channels)`` (``n_blocks`` may be
        zero).
    times : sequence of array-like
        For each frequency, the centres of its blocks in seconds from
        the first sample, ``(n_blocks,)``.
    realizations : sequence of array-like
        For each frequency, the number of products averaged in each
        block (windows times sub-frequencies), ``(n_blocks,)``.
    channels : sequence of str, optional
        Labels of the channels (e.g. station identifiers and
        components).

    Raises
    ------
    TypeError
        If ``channels`` is not a sequence of strings.
    ValueError
        If the shapes are inconsistent.

    Examples
    --------
    >>> spectra = CrossSpectra.from_data(data, 0.01, [1.0, 2.0],
    ...                                  cycles=20, overlap=0.5,
    ...                                  windows_per_block=30)
    >>> values, vectors = spectra.decompose(0)
    """

    def __init__(
        self,
        frequency: Any,
        matrices: Sequence,
        times: Sequence,
        realizations: Sequence,
        channels: Optional[Sequence[str]] = None,
    ) -> None:
        frequency = np.atleast_1d(np.asarray(frequency, dtype=float))

        if frequency.ndim != 1:
            raise ValueError("frequency must be one-dimensional")

        nfreq = frequency.size

        if not len(matrices) == len(times) == len(realizations) == nfreq:
            raise ValueError(
                "matrices, times and realizations need one item for each "
                "frequency"
            )

        size = None
        stored_matrices = []
        stored_times = []
        stored_counts = []

        for index in range(nfreq):
            matrix = np.asarray(matrices[index], dtype=complex)
            time = np.asarray(times[index], dtype=float).reshape(-1)
            count = np.asarray(realizations[index]).reshape(-1)

            if matrix.ndim != 3 or matrix.shape[1] != matrix.shape[2]:
                raise ValueError(
                    "matrices must be (n_blocks, n_channels, n_channels)"
                )

            if size is None:
                size = matrix.shape[1]

            elif matrix.shape[1] != size:
                raise ValueError("the number of channels changes")

            if not time.size == count.size == matrix.shape[0]:
                raise ValueError(
                    "times and realizations need one value for each block"
                )

            if count.size and (
                not np.issubdtype(count.dtype, np.integer)
                or np.any(count < 1)
            ):
                raise ValueError("realizations must be positive integers")

            stored_matrices.append(matrix)
            stored_times.append(time)
            stored_counts.append(count.astype(np.int64))

        if channels is not None:
            if isinstance(channels, str) or not all(
                isinstance(name, str) for name in channels
            ):
                raise TypeError("channels must be a sequence of strings")

            channels = tuple(channels)

            if size is not None and len(channels) != size:
                raise ValueError(
                    f"{len(channels)} channel labels for {size} channels"
                )

        self._frequency = frequency
        self._matrices = stored_matrices
        self._times = stored_times
        self._realizations = stored_counts
        self._size = 0 if size is None else size
        self._channels = channels

    def __repr__(self) -> str:
        return (
            f"CrossSpectra({self._size} channels, {self.nfreq} "
            f"frequencies, {int(self.nblocks.sum())} blocks)"
        )

    # -- construction --------------------------------------------------

    @classmethod
    def from_data(
        cls,
        data: Any,
        delta: float,
        freq: Any,
        *,
        duration: Optional[float] = None,
        cycles: Optional[float] = None,
        overlap: Optional[float] = None,
        window: Any = ("tukey", 0.1),
        windows_per_block: Optional[int] = None,
        block_overlap: float = 0.0,
        subfrequencies: int = 1,
        bandwidth: float = 0.0,
        scaling: Optional[str] = "psd",
        channels: Optional[Sequence[str]] = None,
    ) -> CrossSpectra:
        """
        Return the cross-spectral matrices of a set of channels.

        Parameters
        ----------
        data : array-like
            Real samples ``(n_channels, nsamp)``, synchronized.
        delta : float
            Sampling interval in seconds.
        freq : float or array-like
            Frequencies in Hz.
        duration : float, optional
            Duration of the windows in seconds (fixed windows).
        cycles : float, optional
            Number of cycles of each frequency in its windows;
            alternative to ``duration`` (one of the two is required).
        overlap : float, optional
            Fractional overlap between consecutive windows, in
            ``[0, 1)``. Default is no overlap.
        window : str, tuple or array-like, optional
            Taper (see :func:`shakelab.signals.windowing.taper`).
            Default is a Tukey window with 10% tapered edges.
        windows_per_block : int, optional
            Consecutive windows averaged in each block. Default: all
            the windows in one block.
        block_overlap : float, optional
            Fractional overlap between consecutive blocks, in windows,
            in ``[0, 1)``. Default is 0.
        subfrequencies : int, optional
            Number of frequencies averaged around each frequency,
            equally spaced in logarithm between ``f * exp(-bandwidth /
            2)`` and ``f * exp(bandwidth / 2)`` (about ``f * (1 +/-
            bandwidth / 2)``). Default is 1 (no averaging).
        bandwidth : float, optional
            Relative width of the band of the sub-frequencies; positive
            if and only if ``subfrequencies`` is larger than one.
        scaling : {"psd", None}, optional
            Scaling of the coefficients (see :func:`coefficients`):
            with ``"psd"`` (default) the matrices are one-sided cross
            power spectral densities.
        channels : sequence of str, optional
            Labels of the channels.

        Returns
        -------
        CrossSpectra
            Matrices of every frequency; a frequency whose window does
            not fit in the data, or whose windows do not fill a block,
            has no blocks.

        Raises
        ------
        TypeError
            If an argument has the wrong type.
        ValueError
            If the data are not two-dimensional and finite, the window
            or the sub-frequencies are inconsistent, or a frequency
            (or sub-frequency) is outside ``(0, Nyquist)``.

        Notes
        -----
        The realizations of a block are the number of products
        averaged (windows times sub-frequencies); overlapping windows
        and close sub-frequencies are correlated, so the effective
        number of independent realizations is smaller (see
        :func:`shakelab.signals.windowing.effective_realizations`).
        """
        data = np.asarray(data, dtype=float)

        if data.ndim != 2 or not np.all(np.isfinite(data)):
            raise ValueError(
                "data must be a finite (n_channels, nsamp) array"
            )

        freq = np.atleast_1d(np.asarray(freq, dtype=float))

        if freq.ndim != 1 or freq.size == 0:
            raise ValueError("freq must be one or more frequencies")

        offsets = _subfrequency_offsets(subfrequencies, bandwidth)

        if windows_per_block is not None:
            if (
                isinstance(windows_per_block, (bool, np.bool_))
                or not isinstance(windows_per_block, (int, np.integer))
            ):
                raise TypeError("windows_per_block must be an integer")

            if windows_per_block < 1:
                raise ValueError("windows_per_block must be at least 1")

        # Check the block overlap once (any number of windows).
        window_starts(1, 1, overlap=block_overlap)

        nsamp = data.shape[-1]
        sizes = [
            window_size(delta, duration=duration, cycles=cycles,
                        frequency=None if cycles is None else value)
            for value in freq
        ]

        matrices: list = [None] * freq.size
        times: list = [None] * freq.size
        counts: list = [None] * freq.size
        empty = np.zeros((0, data.shape[0], data.shape[0]), dtype=complex)

        # One transform for the frequencies sharing a window length.
        for size in sorted(set(sizes)):
            members = [k for k, value in enumerate(sizes) if value == size]
            starts = window_starts(nsamp, size, overlap=overlap)

            if starts.size == 0:
                for k in members:
                    matrices[k], times[k], counts[k] = empty, [], []

                continue

            targets = (freq[members][:, None] * offsets[None, :]).ravel()
            _, values = coefficients(
                data, delta, targets, size, overlap=overlap,
                window=window, scaling=scaling,
            )
            # (n_members, n_sub, n_windows, n_channels)
            values = values.reshape(
                len(members), offsets.size, starts.size, data.shape[0]
            )
            per_block = (
                starts.size if windows_per_block is None
                else windows_per_block
            )
            blocks = window_starts(
                starts.size, per_block, overlap=block_overlap
            )
            centres = (
                (starts[blocks] + starts[blocks + per_block - 1] + size)
                * 0.5 * float(delta)
            )

            for row, k in enumerate(members):
                matrices[k] = _block_matrices(values[row], blocks, per_block)
                times[k] = centres
                counts[k] = np.full(
                    blocks.size, per_block * offsets.size, dtype=np.int64
                )

        return cls(freq, matrices, times, counts, channels=channels)

    # -- properties ----------------------------------------------------

    @property
    def frequency(self) -> np.ndarray:
        """Frequencies in Hz (a copy)."""
        return self._frequency.copy()

    @property
    def nfreq(self) -> int:
        """Number of frequencies."""
        return self._frequency.size

    @property
    def nchannels(self) -> int:
        """Number of channels."""
        return self._size

    @property
    def channels(self) -> Optional[tuple[str, ...]]:
        """Labels of the channels, or None."""
        return self._channels

    @property
    def nblocks(self) -> np.ndarray:
        """Number of blocks at each frequency."""
        return np.array([m.shape[0] for m in self._matrices], dtype=int)

    # -- access --------------------------------------------------------

    def matrices(self, index: int) -> np.ndarray:
        """
        Return the matrices of the blocks of a frequency.

        Parameters
        ----------
        index : int
            Index of the frequency.

        Returns
        -------
        numpy.ndarray
            Copy of the matrices, ``(n_blocks, n_channels,
            n_channels)``.
        """
        return self._matrices[self._index(index)].copy()

    def times(self, index: int) -> np.ndarray:
        """
        Return the centres of the blocks of a frequency.

        Parameters
        ----------
        index : int
            Index of the frequency.

        Returns
        -------
        numpy.ndarray
            Seconds from the first sample, ``(n_blocks,)``.
        """
        return self._times[self._index(index)].copy()

    def realizations(self, index: int) -> np.ndarray:
        """
        Return the number of products averaged in each block.

        Parameters
        ----------
        index : int
            Index of the frequency.

        Returns
        -------
        numpy.ndarray
            Windows times sub-frequencies, ``(n_blocks,)``.
        """
        return self._realizations[self._index(index)].copy()

    def mean(self) -> np.ndarray:
        """
        Return the matrices averaged over the blocks.

        Returns
        -------
        numpy.ndarray
            ``(n_freq, n_channels, n_channels)``, blocks weighted by
            their realizations; NaN at frequencies without blocks.
        """
        out = np.full((self.nfreq, self._size, self._size), np.nan,
                      dtype=complex)

        for index, matrix in enumerate(self._matrices):
            weights = self._realizations[index]

            if weights.size:
                out[index] = np.tensordot(weights, matrix, axes=1) / (
                    weights.sum()
                )

        return out

    def coherence(self) -> CrossSpectra:
        """
        Return the coherence matrices scaled to the mean power.

        Each matrix becomes ``C[i, j] / sqrt(C[i, i] * C[j, j])`` (unit
        diagonal, off-diagonal moduli up to one) multiplied by the mean
        of its diagonal, ``trace(C) / n``: the differences of gain
        between the channels are removed, the overall power is kept.

        Returns
        -------
        CrossSpectra
            New object with the same frequencies, times and
            realizations.

        Raises
        ------
        ValueError
            If a channel has no power in a block.
        """
        matrices = []

        for matrix in self._matrices:
            power = np.real(np.einsum("bii->bi", matrix))

            if np.any(power <= 0.0):
                raise ValueError(
                    "the coherence needs channels with positive power"
                )

            scale = np.sqrt(power)
            normalized = matrix / (scale[:, :, None] * scale[:, None, :])
            matrices.append(normalized * power.mean(axis=1)[:, None, None])

        return CrossSpectra(
            self._frequency, matrices, self._times, self._realizations,
            channels=self._channels,
        )

    def decompose(self, index: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the eigenvalues and eigenvectors of the blocks.

        For the Hermitian, positive semi-definite matrices the
        eigendecomposition coincides with the singular value
        decomposition of the frequency domain decomposition.

        Parameters
        ----------
        index : int
            Index of the frequency.

        Returns
        -------
        values : numpy.ndarray
            Eigenvalues in decreasing order, ``(n_blocks, n_channels)``
            (rounding errors below zero are set to zero).
        vectors : numpy.ndarray
            Unit eigenvectors in the columns, ``(n_blocks, n_channels,
            n_channels)``, in the order of ``values``.
        """
        matrix = self._matrices[self._index(index)]
        values, vectors = np.linalg.eigh(matrix)

        return (
            np.maximum(values[:, ::-1], 0.0),
            vectors[:, :, ::-1],
        )

    def copy(self) -> CrossSpectra:
        """
        Return an independent copy.

        Returns
        -------
        CrossSpectra
            Copy of the matrices, times and realizations.
        """
        return CrossSpectra(
            self._frequency, [m.copy() for m in self._matrices],
            [t.copy() for t in self._times],
            [c.copy() for c in self._realizations],
            channels=self._channels,
        )

    # -- summary -------------------------------------------------------

    def summary(self) -> str:
        """
        Return a summary as text.

        Returns
        -------
        str
            Channels, frequencies and blocks.
        """
        lines = [repr(self)]

        if self._channels is not None:
            lines.append("  channels: " + ", ".join(self._channels))

        if self.nfreq:
            lines.append(
                f"  frequencies: {self._frequency.min():g} - "
                f"{self._frequency.max():g} Hz"
            )
            nblocks = self.nblocks
            lines.append(
                f"  blocks per frequency: {nblocks.min()} - "
                f"{nblocks.max()}"
            )

        return "\n".join(lines)

    def info(self) -> None:
        """Print the summary (see :meth:`summary`)."""
        print(self.summary())

    # -- helpers -------------------------------------------------------

    def _index(self, index: Any) -> int:
        """Return a valid frequency index."""
        if (
            isinstance(index, (bool, np.bool_))
            or not isinstance(index, (int, np.integer))
        ):
            raise TypeError("the frequency index must be an integer")

        if not -self.nfreq <= index < self.nfreq:
            raise IndexError(f"frequency index {index} out of range")

        return int(index) % self.nfreq


def _subfrequency_offsets(count: Any, bandwidth: Any) -> np.ndarray:
    """Return the relative sub-frequencies ``exp(linspace(-b/2, b/2))``."""
    if (
        isinstance(count, (bool, np.bool_))
        or not isinstance(count, (int, np.integer))
    ):
        raise TypeError("subfrequencies must be an integer")

    if count < 1:
        raise ValueError("subfrequencies must be at least 1")

    try:
        bandwidth = float(bandwidth)
    except (TypeError, ValueError):
        raise TypeError("bandwidth must be a real number") from None

    if not np.isfinite(bandwidth) or bandwidth < 0.0:
        raise ValueError("bandwidth must be zero or positive")

    if (count > 1) != (bandwidth > 0.0):
        raise ValueError(
            "a positive bandwidth needs more than one sub-frequency, "
            "and more sub-frequencies a positive bandwidth"
        )

    return np.exp(np.linspace(-0.5, 0.5, count) * bandwidth)


def _block_matrices(
    values: np.ndarray,
    blocks: np.ndarray,
    per_block: int,
) -> np.ndarray:
    """Return the matrices of the blocks from ``(n_sub, n_win, n_ch)``."""
    nsub, _, nch = values.shape
    out = np.empty((blocks.size, nch, nch), dtype=complex)

    for row, first in enumerate(blocks):
        series = values[:, first:first + per_block, :].reshape(-1, nch)
        out[row] = cross_spectral_matrix(series)

    return out
