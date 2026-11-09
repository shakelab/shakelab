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
Time-frequency analysis.

The module has two levels:

- functions on arrays, along any axis, used also by higher modules
  such as :mod:`shakelab.signals.wavenumber`: the windowed Fourier
  transform (:func:`wft_transform`), the continuous wavelet transform
  with a complex Morlet wavelet (:func:`cwt_transform`) and the
  Stockwell transform (:func:`st_transform`);
- the :class:`Spectrogram` of a :class:`~shakelab.signals.base.Record`,
  with the metadata of the record, time and frequency axes and the edge
  effects of each method (``Record.spectrogram`` builds it).

All transforms return the frequency as the first axis; the time axis
keeps the position of the transformed axis of the input array.

Scaling
-------
By default (``scaling="amplitude"``) the three transforms are
calibrated in the units of the data: a sinusoid of amplitude ``A``
gives coefficients of modulus ``A`` at its frequency, so that the
methods can be compared directly. The windowed Fourier transform also
offers ``scaling="psd"``, whose squared modulus is a one-sided power
spectral density (units squared per Hz, as in Welch's method).

References
----------
Torrence, C., and Compo, G. P. (1998). A practical guide to wavelet
analysis. Bulletin of the American Meteorological Society, 79(1),
61-78.

Stockwell, R. G., Mansinha, L., and Lowe, R. P. (1996). Localization of
the complex spectrum: the S transform. IEEE Transactions on Signal
Processing, 44(4), 998-1001.
"""

from __future__ import annotations

from numbers import Integral, Real
from typing import Any, Optional

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import signal

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.header import Header

__all__ = [
    "METHODS",
    "Spectrogram",
    "cwt_coi",
    "cwt_transform",
    "morlet_wavelet",
    "power_to_db",
    "st_transform",
    "wft_transform",
]

# Methods of :class:`Spectrogram`, by name.
METHODS = ("wft", "cwt", "st")

# Half support of the Morlet kernel, in Gaussian standard deviations.
_MORLET_SUPPORT = 5.0


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _check_delta(delta: Any) -> float:
    """Return a positive sampling interval."""
    if (
        not isinstance(delta, Real)
        or isinstance(delta, (bool, np.bool_))
    ):
        raise TypeError("delta must be a real number")

    delta = float(delta)

    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError("delta must be positive")

    return delta


def _check_axis(axis: Any, ndim: int) -> int:
    """Return the non-negative index of an axis of an array."""
    if (
        not isinstance(axis, Integral)
        or isinstance(axis, (bool, np.bool_))
    ):
        raise TypeError("axis must be an integer")

    if not -ndim <= axis < ndim:
        raise ValueError(
            f"axis {axis} is out of bounds for an array of dimension "
            f"{ndim}"
        )

    return int(axis) % ndim


def _check_data(data: Any) -> np.ndarray:
    """Return the input data as an array with at least one dimension."""
    data = np.asarray(data)

    if data.ndim == 0:
        raise ValueError("data must have at least one dimension")

    return data


def _check_scaling(scaling: Any, allowed: tuple) -> Any:
    """Return a valid scaling name."""
    if scaling not in allowed:
        names = ", ".join(repr(value) for value in allowed)
        raise ValueError(f"scaling must be one of {names}")

    return scaling


def _one_sided_factor(index: np.ndarray, nsamp: int) -> np.ndarray:
    """
    Return the one-sided amplitude factor of DFT frequency indices.

    The factor is 2 for the frequencies with a negative counterpart
    and 1 for zero and for the Nyquist frequency of an even length.
    """
    index = np.asarray(index)
    factor = np.full(index.shape, 2.0)
    factor[index == 0] = 1.0

    if nsamp % 2 == 0:
        factor[index == nsamp // 2] = 1.0

    return factor


def _frequency_first(values: np.ndarray, axis: int) -> np.ndarray:
    """
    Move the frequency axis from the end to the front.

    ``values`` has the transformed (time) axis in the next to last
    position and frequency last; the time axis is moved back to
    ``axis``, the frequency axis to position 0.
    """
    values = np.moveaxis(values, -1, 0)
    return np.moveaxis(values, -1, axis + 1)


# ----------------------------------------------------------------------------
# Windowed Fourier transform
# ----------------------------------------------------------------------------


def wft_transform(
    data: Any,
    delta: float,
    window_length: int,
    step: Optional[int] = None,
    window: Any = "hann",
    nfft: Optional[int] = None,
    centred: bool = True,
    axis: int = -1,
    scaling: Optional[str] = "amplitude",
    power: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute a windowed (short-time) Fourier transform.

    Parameters
    ----------
    data : array-like
        Input data array.
    delta : float
        Sampling interval along the selected axis, in seconds.
    window_length : int
        Length of the moving window, in samples.
    step : int, optional
        Step between consecutive windows, in samples. Default is half
        the window length (50% overlap).
    window : str, tuple or array-like, optional
        Taper: a window specification of ``scipy.signal.get_window``
        (name, or tuple of name and parameters, e.g.
        ``("tukey", 0.25)``) or an array of ``window_length`` values.
        Default is ``"hann"``.
    nfft : int, optional
        Length of each transform (zero padding of the windows); at
        least ``window_length``. Default is ``window_length``.
    centred : bool, optional
        If ``True`` (default), the data are padded with zeros at both
        ends so that the window centres fall on the samples
        ``0, step, 2 * step, ...`` and cover the whole series. If
        ``False``, only windows fully inside the data are used.
    axis : int, optional
        Axis along which the transform is computed. Default is the
        last axis.
    scaling : {"amplitude", "psd", None}, optional
        ``"amplitude"`` (default): a sinusoid of amplitude ``A`` gives
        coefficients of modulus ``A``. ``"psd"``: the squared modulus
        is the one-sided power spectral density, in units squared per
        Hz. ``None``: the plain DFT of the tapered windows.
    power : bool, optional
        If ``True``, return the squared modulus of the coefficients.
        Default is ``False``.

    Returns
    -------
    freq : numpy.ndarray
        Frequencies, in Hz.
    time : numpy.ndarray
        Times of the window centres, in seconds from the first sample.
    spectrum : numpy.ndarray
        Coefficients with shape ``(n_freq,) + data.shape``, where the
        transformed axis has ``n_windows`` elements.

    Raises
    ------
    TypeError
        If an argument has the wrong type.
    ValueError
        If an argument is out of range or the window is longer than the
        data (with ``centred=False``).
    """
    data = _check_data(data)
    delta = _check_delta(delta)
    axis = _check_axis(axis, data.ndim)
    scaling = _check_scaling(scaling, ("amplitude", "psd", None))

    if (
        not isinstance(window_length, Integral)
        or isinstance(window_length, (bool, np.bool_))
    ):
        raise TypeError("window_length must be an integer")

    window_length = int(window_length)

    if window_length < 1:
        raise ValueError("window_length must be positive")

    if step is None:
        step = max(window_length // 2, 1)

    if not isinstance(step, Integral) or isinstance(step, bool):
        raise TypeError("step must be an integer")

    step = int(step)

    if step < 1:
        raise ValueError("step must be positive")

    if nfft is None:
        nfft = window_length

    if not isinstance(nfft, Integral) or isinstance(nfft, bool):
        raise TypeError("nfft must be an integer")

    nfft = int(nfft)

    if nfft < window_length:
        raise ValueError("nfft must not be smaller than window_length")

    if not isinstance(centred, bool):
        raise TypeError("centred must be a boolean")

    if isinstance(window, (str, tuple)):
        taper = signal.get_window(window, window_length)
    else:
        taper = np.asarray(window, dtype=float)

        if taper.shape != (window_length,):
            raise ValueError(
                "window array must have window_length values"
            )

    moved = np.moveaxis(data, axis, -1)
    nsamp = moved.shape[-1]
    half = window_length // 2

    if centred:
        nwin = (nsamp - 1) // step + 1
        right = max((nwin - 1) * step + window_length - nsamp - half, 0)
        pad = [(0, 0)] * (moved.ndim - 1) + [(half, right)]
        moved = np.pad(moved, pad)
        centres = np.arange(nwin) * step
    else:
        if window_length > nsamp:
            raise ValueError("window_length exceeds the data length")

        nwin = (nsamp - window_length) // step + 1
        centres = np.arange(nwin) * step + half

    segments = sliding_window_view(moved, window_length, axis=-1)
    segments = segments[..., :nwin * step:step, :] * taper

    coeff = np.fft.rfft(segments, n=nfft, axis=-1)
    index = np.arange(coeff.shape[-1])

    # Phase referred to the window centre (the time of the
    # coefficients), as for the wavelet and Stockwell transforms.
    coeff *= np.exp(2j * np.pi * index * half / nfft)

    if scaling == "amplitude":
        coeff *= _one_sided_factor(index, nfft) / np.sum(taper)
    elif scaling == "psd":
        density = np.sum(taper ** 2) / delta
        coeff *= np.sqrt(_one_sided_factor(index, nfft) / density)

    freq = np.fft.rfftfreq(nfft, d=delta)
    time = centres * delta
    spectrum = _frequency_first(coeff, axis)

    if power:
        spectrum = np.abs(spectrum) ** 2

    return freq, time, spectrum


# ----------------------------------------------------------------------------
# Continuous wavelet transform
# ----------------------------------------------------------------------------


def morlet_wavelet(
    length: int,
    width: float,
    omega0: float = 6.0,
    norm: str = "energy",
) -> np.ndarray:
    """
    Generate a complex Morlet wavelet.

    Parameters
    ----------
    length : int
        Number of samples of the wavelet; the wavelet is centred on
        sample ``(length - 1) / 2`` (an odd length centres it on a
        sample).
    width : float
        Wavelet scale (standard deviation of the Gaussian envelope), in
        samples.
    omega0 : float, optional
        Non-dimensional central angular frequency. Default is 6.
    norm : {"energy", "amplitude"}, optional
        ``"energy"`` (default): unit energy. ``"amplitude"``: the
        Gaussian envelope sums to one, so that the correlation with a
        complex exponential of unit modulus at the central frequency
        has unit modulus.

    Returns
    -------
    numpy.ndarray
        Complex wavelet samples.

    Raises
    ------
    ValueError
        If ``length`` or ``width`` is not positive or ``norm`` is
        unknown.
    """
    if length <= 0:
        raise ValueError("length must be positive")

    if width <= 0:
        raise ValueError("width must be positive")

    x = np.arange(length) - (length - 1) / 2.0
    envelope = np.exp(-(x ** 2) / (2.0 * width ** 2))
    wavelet = envelope * np.exp(1j * omega0 * x / width)

    if norm == "energy":
        return wavelet / np.sqrt(np.sum(envelope ** 2))

    if norm == "amplitude":
        return wavelet / np.sum(envelope)

    raise ValueError("norm must be 'energy' or 'amplitude'")


def _cwt_width_transform(
    data: Any,
    widths: Any,
    omega0: float = 6.0,
    axis: int = -1,
    power: bool = False,
    scaling: str = "energy",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a Morlet continuous wavelet transform from wavelet widths.

    Parameters
    ----------
    data : array-like
        Input data array.
    widths : array-like
        Wavelet scales, in samples.
    omega0 : float, optional
        Non-dimensional central angular frequency.
    axis : int, optional
        Axis along which the transform is computed.
    power : bool, optional
        If ``True``, return the squared modulus.
    scaling : {"energy", "amplitude"}, optional
        ``"energy"``: wavelets of unit energy. ``"amplitude"``: a
        sinusoid of amplitude ``A`` gives coefficients of modulus
        ``A`` at its frequency.

    Returns
    -------
    widths : numpy.ndarray
        Wavelet scales, in samples.
    transform : numpy.ndarray
        Coefficients with shape ``(n_widths,) + data.shape``.
    """
    data = _check_data(data)
    axis = _check_axis(axis, data.ndim)
    scaling = _check_scaling(scaling, ("energy", "amplitude"))
    widths = np.atleast_1d(np.asarray(widths, dtype=float))

    if widths.ndim != 1:
        raise ValueError("widths must be one-dimensional")

    if np.any(~np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("all widths must be positive")

    moved = np.moveaxis(data, axis, -1)
    nsamp = moved.shape[-1]
    transform = np.zeros(
        moved.shape + (widths.size,), dtype=np.complex128
    )

    for iwidth, width in enumerate(widths):
        # Odd kernel centred on a sample, truncated where the Gaussian
        # envelope is negligible and never longer than needed.
        half = min(int(np.ceil(_MORLET_SUPPORT * width)), nsamp)
        wavelet = morlet_wavelet(
            2 * half + 1, width, omega0=omega0, norm=scaling
        )

        # Correlation with the wavelet = convolution with the
        # time-reversed conjugate.
        kernel = np.conj(wavelet[::-1])
        kernel = kernel.reshape((1,) * (moved.ndim - 1) + (-1,))

        transform[..., iwidth] = signal.fftconvolve(
            moved, kernel, mode="same", axes=-1
        )

    if scaling == "amplitude":
        transform *= 2.0

    transform = _frequency_first(transform, axis)

    if power:
        transform = np.abs(transform) ** 2

    return widths, transform


def cwt_transform(
    data: Any,
    delta: float,
    freq: Any,
    omega0: float = 6.0,
    axis: int = -1,
    scaling: str = "amplitude",
    power: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a continuous wavelet transform with a Morlet wavelet.

    Parameters
    ----------
    data : array-like
        Input data array.
    delta : float
        Sampling interval along the selected axis, in seconds.
    freq : array-like
        Frequencies to analyse, in Hz, up to the Nyquist frequency.
    omega0 : float, optional
        Non-dimensional central angular frequency of the wavelet.
        Default is 6.
    axis : int, optional
        Axis along which the transform is computed. Default is the
        last axis.
    scaling : {"amplitude", "energy"}, optional
        ``"amplitude"`` (default): a sinusoid of amplitude ``A`` gives
        coefficients of modulus ``A`` at its frequency. ``"energy"``:
        wavelets of unit energy (the modulus of a sinusoid grows with
        the scale, as ``1 / sqrt(f)``).
    power : bool, optional
        If ``True``, return the squared modulus. Default is ``False``.

    Returns
    -------
    freq : numpy.ndarray
        Analysed frequencies, in Hz.
    transform : numpy.ndarray
        Coefficients with shape ``(n_freq,) + data.shape``.

    Raises
    ------
    ValueError
        If a frequency is not positive or exceeds the Nyquist frequency.

    Notes
    -----
    The wavelet width (Gaussian standard deviation, in samples) of the
    frequency ``f`` is ``omega0 / (2 * pi * f * delta)``, so that the
    central frequency of the wavelet is ``f``. The kernel is truncated
    at five standard deviations; the edges of the series are affected
    within :func:`cwt_coi` of each end.
    """
    delta = _check_delta(delta)
    freq = np.atleast_1d(np.asarray(freq, dtype=float))

    if np.any(~np.isfinite(freq)) or np.any(freq <= 0):
        raise ValueError("all frequencies must be positive")

    if np.any(freq > 0.5 / delta):
        raise ValueError("frequencies must not exceed the Nyquist frequency")

    widths = omega0 / (2.0 * np.pi * freq * delta)

    _, transform = _cwt_width_transform(
        data, widths, omega0=omega0, axis=axis, power=power,
        scaling=scaling,
    )

    return freq, transform


def cwt_coi(freq: Any, omega0: float = 6.0) -> np.ndarray:
    """
    Return the cone of influence of the Morlet transform.

    Parameters
    ----------
    freq : array-like
        Frequencies, in Hz.
    omega0 : float, optional
        Non-dimensional central angular frequency. Default is 6.

    Returns
    -------
    numpy.ndarray
        For each frequency, the time from either end of the series
        within which the coefficients are affected by the edges, in
        seconds.

    Notes
    -----
    The e-folding time of the wavelet power of an edge discontinuity,
    ``sqrt(2)`` times the Gaussian standard deviation of the wavelet
    (Torrence and Compo, 1998).
    """
    freq = np.asarray(freq, dtype=float)
    return np.sqrt(2.0) * omega0 / (2.0 * np.pi * freq)


# ----------------------------------------------------------------------------
# Stockwell transform
# ----------------------------------------------------------------------------


def st_transform(
    data: Any,
    delta: float,
    freq: Any = None,
    axis: int = -1,
    scaling: Optional[str] = "amplitude",
    power: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the Stockwell transform.

    Parameters
    ----------
    data : array-like
        Input data array.
    delta : float
        Sampling interval along the selected axis, in seconds.
    freq : array-like, optional
        Frequencies to compute, in Hz, from 0 to the Nyquist frequency;
        each is replaced by the nearest DFT frequency ``k / (n *
        delta)`` and duplicates are removed. Default is every DFT
        frequency (memory grows as ``n ** 2``).
    axis : int, optional
        Axis along which the transform is computed. Default is the
        last axis.
    scaling : {"amplitude", None}, optional
        ``"amplitude"`` (default): a sinusoid of amplitude ``A`` gives
        coefficients of modulus ``A`` at its frequency. ``None``: the
        definition of Stockwell et al. (1996), half of it except at
        zero and Nyquist frequency.
    power : bool, optional
        If ``True``, return the squared modulus. Default is ``False``.

    Returns
    -------
    freq : numpy.ndarray
        Frequencies of the transform, in Hz.
    transform : numpy.ndarray
        Coefficients with shape ``(n_freq,) + data.shape``.

    Raises
    ------
    ValueError
        If a frequency is negative or exceeds the Nyquist frequency.

    Notes
    -----
    Row ``k > 0`` is the inverse DFT of the spectrum shifted by ``k``
    and multiplied by the Gaussian ``exp(-2 pi^2 m^2 / k^2)``; row 0 is
    the mean of the series. The transform is circular: the edges of the
    series interact within about ``sqrt(2) / f`` seconds.
    """
    data = _check_data(data)
    delta = _check_delta(delta)
    axis = _check_axis(axis, data.ndim)
    scaling = _check_scaling(scaling, ("amplitude", None))

    moved = np.moveaxis(data, axis, -1)
    nsamp = moved.shape[-1]
    nyquist = nsamp // 2

    if freq is None:
        index = np.arange(nyquist + 1)
    else:
        freq = np.atleast_1d(np.asarray(freq, dtype=float))

        if np.any(~np.isfinite(freq)) or np.any(freq < 0):
            raise ValueError("frequencies must not be negative")

        if np.any(freq > 0.5 / delta * (1.0 + 1e-9)):
            raise ValueError(
                "frequencies must not exceed the Nyquist frequency"
            )

        index = np.unique(np.rint(freq * nsamp * delta).astype(int))
        index = index[index <= nyquist]

    spectrum = np.fft.fft(moved, axis=-1)
    offset = np.fft.fftfreq(nsamp) * nsamp
    transform = np.zeros(
        moved.shape + (index.size,), dtype=np.complex128
    )

    for irow, k in enumerate(index):
        if k == 0:
            transform[..., irow] = np.mean(moved, axis=-1)[..., None]
            continue

        gauss = np.exp(-2.0 * np.pi ** 2 * offset ** 2 / k ** 2)
        shifted = np.roll(spectrum, -k, axis=-1)
        transform[..., irow] = np.fft.ifft(shifted * gauss, axis=-1)

    if scaling == "amplitude":
        transform *= _one_sided_factor(index, nsamp)

    transform = _frequency_first(transform, axis)

    if power:
        transform = np.abs(transform) ** 2

    return index / (nsamp * delta), transform


# ----------------------------------------------------------------------------
# Decibels
# ----------------------------------------------------------------------------


def power_to_db(
    power: Any,
    reference: float = 1.0,
    floor: Optional[float] = None,
) -> np.ndarray:
    """
    Convert power values to decibels.

    Parameters
    ----------
    power : array-like
        Power values (squared modulus), non-negative.
    reference : float, optional
        Power of 0 dB. Default is 1 (decibels relative to one unit
        squared).
    floor : float, optional
        Dynamic range in dB: values more than ``floor`` dB below the
        maximum are raised to that level (zeros included). Default is
        no floor (zeros give ``-inf``).

    Returns
    -------
    numpy.ndarray
        ``10 * log10(power / reference)``.

    Raises
    ------
    ValueError
        If ``reference`` or ``floor`` is not positive.
    """
    power = np.asarray(power, dtype=float)

    if reference <= 0:
        raise ValueError("reference must be positive")

    if floor is not None:
        if floor <= 0:
            raise ValueError("floor must be positive")

        peak = np.max(power) if power.size else 0.0

        if peak > 0:
            power = np.maximum(power, peak * 10.0 ** (-floor / 10.0))

    with np.errstate(divide="ignore"):
        return 10.0 * np.log10(power / reference)


# ----------------------------------------------------------------------------
# Spectrogram
# ----------------------------------------------------------------------------


class Spectrogram:
    """
    Time-frequency representation of a record.

    Parameters
    ----------
    data : array-like, optional
        Complex coefficients with shape ``(n_freq, n_time)``.
    frequency : array-like, optional
        Frequencies, in Hz.
    time : array-like, optional
        Times of the coefficients, in seconds from ``head.time``.
    head : Header, optional
        Metadata of the record; a copy is stored.
    method : str, optional
        Transform: ``"wft"``, ``"cwt"`` or ``"st"``.
    scaling : str, optional
        Scaling of the coefficients (see the transforms).
    coi : array-like, optional
        For each frequency, the time from either end of the record
        within which the coefficients are affected by the edges, in
        seconds.
    duration : float, optional
        Duration of the record, in seconds (end of the edge effects).
        Default is the last time.
    parameters : dict, optional
        Parameters of the transform.

    Attributes
    ----------
    head : Header
        Metadata of the record (identifier, start time, units).
    data : numpy.ndarray
        Complex coefficients, ``(n_freq, n_time)``.
    frequency : numpy.ndarray
        Frequencies, in Hz.
    time : numpy.ndarray
        Times, in seconds from the start of the record.
    method : str or None
        Transform.
    scaling : str or None
        Scaling of the coefficients.
    coi : numpy.ndarray or None
        Edge effects for each frequency, in seconds.
    duration : float or None
        Duration of the record, in seconds.
    parameters : dict
        Parameters of the transform.

    Raises
    ------
    ValueError
        If the shapes of the data and of the axes do not match.
    """

    def __init__(
        self,
        data: Any = None,
        frequency: Any = None,
        time: Any = None,
        head: Optional[Header] = None,
        method: Optional[str] = None,
        scaling: Optional[str] = None,
        coi: Any = None,
        duration: Optional[float] = None,
        parameters: Optional[dict] = None,
    ) -> None:
        self.head = Header() if head is None else head.copy()
        self.method = method
        self.scaling = scaling
        self.parameters = {} if parameters is None else dict(parameters)

        data = np.zeros((0, 0), dtype=complex) if data is None else data
        self.data = np.array(data, dtype=complex, ndmin=2)

        nfreq, ntime = self.data.shape
        self.frequency = self._axis(frequency, nfreq, "frequency")
        self.time = self._axis(time, ntime, "time")
        self.coi = None if coi is None else self._axis(coi, nfreq, "coi")

        if duration is None and ntime:
            duration = float(self.time[-1])

        self.duration = duration

    @staticmethod
    def _axis(values: Any, size: int, name: str) -> np.ndarray:
        """Return an axis of the given size."""
        if values is None:
            if size:
                raise ValueError(f"{name} values are required")
            return np.zeros(0)

        values = np.array(values, dtype=float, ndmin=1)

        if values.shape != (size,):
            raise ValueError(
                f"{values.size} {name} values do not match {size} "
                "coefficients"
            )

        return values

    # -- construction ---------------------------------------------------

    @classmethod
    def from_record(
        cls,
        record: Record,
        method: str = "wft",
        *,
        window: Optional[float] = None,
        overlap: float = 0.5,
        taper: Any = "hann",
        nfft: Optional[int] = None,
        fmin: Optional[float] = None,
        fmax: Optional[float] = None,
        nfreq: Optional[int] = None,
        spacing: Optional[str] = None,
        omega0: Optional[float] = None,
        scaling: Optional[str] = "amplitude",
    ) -> Spectrogram:
        """
        Return the spectrogram of a Record, with a copy of its metadata.

        Parameters
        ----------
        record : Record
            Record with a sampling interval and at least two samples.
        method : {"wft", "cwt", "st"}, optional
            Windowed Fourier transform (default), Morlet wavelet
            transform or Stockwell transform.
        window : float, optional
            ``"wft"``: window length in seconds. Default is one
            twentieth of the record (at least 16 samples).
        overlap : float, optional
            ``"wft"``: overlap of consecutive windows, as a fraction of
            the window, in [0, 1). Default is 0.5.
        taper : str, tuple or array-like, optional
            ``"wft"``: taper of the windows (see :func:`wft_transform`).
            Default is ``"hann"``.
        nfft : int, optional
            ``"wft"``: transform length (zero padding of the windows).
        fmin, fmax : float, optional
            Frequency range, in Hz. ``"cwt"``: default ``fmin`` is the
            frequency whose edge effects leave half of the record,
            ``fmax`` 80% of the Nyquist frequency. ``"st"``: default
            from the lowest non-zero DFT frequency to the Nyquist
            frequency. ``"wft"``: the frequencies of the transform are
            cut to the range.
        nfreq : int, optional
            ``"cwt"``, ``"st"``: number of frequencies. Default is 100.
        spacing : {"log", "linear"}, optional
            ``"cwt"``, ``"st"``: spacing of the frequencies. Default is
            ``"log"``.
        omega0 : float, optional
            ``"cwt"``: central angular frequency of the wavelet.
            Default is 6.
        scaling : str, optional
            Scaling of the coefficients: ``"amplitude"`` (default) for
            every method, ``"psd"`` or None for ``"wft"``, ``"energy"``
            for ``"cwt"``, None for ``"st"``.

        Returns
        -------
        Spectrogram
            Spectrogram of the record samples.

        Raises
        ------
        TypeError
            If ``record`` is not a Record, or an option is not used by
            the method.
        ValueError
            If the record has no sampling interval or fewer than two
            samples, or an option is out of range.
        """
        if not isinstance(record, Record):
            raise TypeError("record must be a Record instance")

        if method not in METHODS:
            names = ", ".join(repr(name) for name in METHODS)
            raise ValueError(f"method must be one of {names}")

        if record.delta is None:
            raise ValueError(
                "sampling interval is required to compute a spectrogram"
            )

        if record.nsamp < 2:
            raise ValueError(
                "at least two samples are required to compute a "
                "spectrogram"
            )

        given = {
            "window": window is not None,
            "overlap": overlap != 0.5,
            "taper": not (isinstance(taper, str) and taper == "hann"),
            "nfft": nfft is not None,
            "nfreq": nfreq is not None,
            "spacing": spacing is not None,
            "omega0": omega0 is not None,
        }
        used = {
            "wft": ("window", "overlap", "taper", "nfft"),
            "cwt": ("nfreq", "spacing", "omega0"),
            "st": ("nfreq", "spacing"),
        }[method]
        unused = [name for name, value in given.items()
                  if value and name not in used]

        if unused:
            raise TypeError(
                f"option(s) {', '.join(unused)} not used by method "
                f"{method!r}"
            )

        delta = record.delta
        data = np.asarray(record.data, dtype=float)
        duration = record.nsamp * delta
        nyquist = 0.5 / delta

        for name, value in (("fmin", fmin), ("fmax", fmax)):
            if value is not None and not 0 <= value <= nyquist:
                raise ValueError(
                    f"{name} must be within 0 and the Nyquist frequency"
                )

        if fmin is not None and fmax is not None and fmax < fmin:
            raise ValueError("fmax must not be smaller than fmin")

        if method == "wft":
            if window is None:
                length = max(record.nsamp // 20, 16)
            else:
                length = int(round(window / delta))

            length = min(length, record.nsamp)

            if length < 2:
                raise ValueError("window must contain at least 2 samples")

            if not 0 <= overlap < 1:
                raise ValueError("overlap must be in [0, 1)")

            step = max(int(round(length * (1.0 - overlap))), 1)
            freq, time, coeff = wft_transform(
                data, delta, length, step=step, window=taper, nfft=nfft,
                scaling=scaling,
            )
            keep = np.ones(freq.size, dtype=bool)

            if fmin is not None:
                keep &= freq >= fmin
            if fmax is not None:
                keep &= freq <= fmax

            freq, coeff = freq[keep], coeff[keep]
            coi = np.full(freq.size, length // 2 * delta)
            parameters = {
                "window_length": length,
                "step": step,
                "taper": taper,
                "nfft": length if nfft is None else nfft,
            }

        else:
            nfreq = 100 if nfreq is None else nfreq
            spacing = "log" if spacing is None else spacing

            if not isinstance(nfreq, Integral) or nfreq < 1:
                raise ValueError("nfreq must be a positive integer")

            if spacing not in ("log", "linear"):
                raise ValueError("spacing must be 'log' or 'linear'")

            time = np.arange(record.nsamp) * delta

            if method == "cwt":
                omega0 = 6.0 if omega0 is None else float(omega0)

                if fmin is None:
                    # cwt_coi(fmin) is one quarter of the record.
                    fmin = 4.0 * np.sqrt(2.0) * omega0 / (
                        2.0 * np.pi * duration
                    )
                    fmin = min(fmin, nyquist)

                if fmax is None:
                    fmax = 0.8 * nyquist

                fmin = max(fmin, 1.0 / duration)
                freq = _frequency_grid(fmin, fmax, nfreq, spacing)
                freq, coeff = cwt_transform(
                    data, delta, freq, omega0=omega0, scaling=scaling
                )
                coi = cwt_coi(freq, omega0)
                parameters = {"omega0": omega0}

            else:
                fmin = 1.0 / duration if fmin is None else fmin
                fmax = nyquist if fmax is None else fmax
                freq = _frequency_grid(
                    max(fmin, 1.0 / duration), fmax, nfreq, spacing
                )

                if fmin == 0:
                    freq = np.concatenate(([0.0], freq))

                freq, coeff = st_transform(
                    data, delta, freq, scaling=scaling
                )

                with np.errstate(divide="ignore"):
                    coi = np.sqrt(2.0) / freq

                parameters = {}

            parameters.update({"nfreq": nfreq, "spacing": spacing})

        return cls(
            coeff, freq, time, head=record.head, method=method,
            scaling=scaling, coi=coi,
            duration=(record.nsamp - 1) * delta, parameters=parameters,
        )

    def copy(self) -> Spectrogram:
        """Return an independent copy of the spectrogram."""
        return Spectrogram(
            self.data.copy(), self.frequency, self.time, head=self.head,
            method=self.method, scaling=self.scaling, coi=self.coi,
            duration=self.duration, parameters=self.parameters,
        )

    # -- properties -----------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"Spectrogram(sid={self.head.sid!r}, method={self.method!r}, "
            f"nfreq={self.nfreq}, ntime={self.ntime})"
        )

    def __str__(self) -> str:
        if self.nfreq:
            band = (
                f"{self.frequency.min():g}-{self.frequency.max():g} Hz"
            )
        else:
            band = "no frequency"

        return (
            f"source identifier {self.head.sid}, method {self.method}, "
            f"scaling {self.scaling}, {self.nfreq} frequencies "
            f"({band}), {self.ntime} times from {self.head.time}"
        )

    def info(self) -> None:
        """Print a human-readable description of the spectrogram."""
        print(self)

    @property
    def sid(self) -> Optional[str]:
        """Source identifier of the record."""
        return self.head.sid

    @property
    def units(self) -> Optional[str]:
        """Units of the record samples."""
        return self.head.units

    @property
    def starttime(self) -> Optional[Date]:
        """Start time of the record (time origin of :attr:`time`)."""
        return self.head.time

    @property
    def nfreq(self) -> int:
        """Number of frequencies."""
        return self.data.shape[0]

    @property
    def ntime(self) -> int:
        """Number of times."""
        return self.data.shape[1]

    @property
    def amplitude(self) -> np.ndarray:
        """Modulus of the coefficients."""
        return np.abs(self.data)

    @property
    def power(self) -> np.ndarray:
        """Squared modulus of the coefficients."""
        return np.abs(self.data) ** 2

    @property
    def phase(self) -> np.ndarray:
        """Phase of the coefficients, in radians, in (-pi, pi]."""
        return np.angle(self.data)

    @property
    def coi_mask(self) -> np.ndarray:
        """
        Coefficients affected by the edges of the record.

        Returns
        -------
        numpy.ndarray
            Boolean array ``(n_freq, n_time)``, True within :attr:`coi`
            of either end of the record; all False without ``coi``.
        """
        mask = np.zeros(self.data.shape, dtype=bool)

        if self.coi is None or not self.ntime:
            return mask

        edge = self.coi[:, None]
        return (self.time < edge) | (self.time > self.duration - edge)

    # -- operations -----------------------------------------------------

    def db(
        self,
        reference: float = 1.0,
        floor: Optional[float] = None,
    ) -> np.ndarray:
        """
        Return the power in decibels.

        Parameters
        ----------
        reference : float, optional
            Power of 0 dB, in units squared. Default is 1.
        floor : float, optional
            Dynamic range in dB below the maximum (see
            :func:`power_to_db`). Default is no floor.

        Returns
        -------
        numpy.ndarray
            ``10 * log10(power / reference)``.
        """
        return power_to_db(self.power, reference=reference, floor=floor)

    def cut(
        self,
        starttime: Any = None,
        endtime: Any = None,
        fmin: Optional[float] = None,
        fmax: Optional[float] = None,
        inplace: bool = True,
    ) -> Optional[Spectrogram]:
        """
        Cut the spectrogram to a time window and a frequency band.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the window. Numerical values are seconds from the
            start of the record. Default is the first time.
        endtime : Date, str or float, optional
            End of the window, included. Default is the last time.
        fmin : float, optional
            Lowest frequency, in Hz, included. Default is no limit.
        fmax : float, optional
            Highest frequency, in Hz, included. Default is no limit.
        inplace : bool, optional
            If ``True`` (default), modify the spectrogram; otherwise
            return a cut copy.

        Returns
        -------
        Spectrogram or None
            The cut spectrogram when ``inplace=False``.

        Raises
        ------
        TypeError
            If ``inplace`` is not a boolean or a time has the wrong
            type.
        ValueError
            If a bound is reversed or nothing is left.
        """
        if not isinstance(inplace, bool):
            raise TypeError("inplace must be a boolean")

        start = self._relative_time(starttime, "starttime")
        end = self._relative_time(endtime, "endtime")

        if start is not None and end is not None and end < start:
            raise ValueError("endtime precedes starttime")

        if fmin is not None and fmax is not None and fmax < fmin:
            raise ValueError("fmax must not be smaller than fmin")

        # Rounding tolerance of the time comparisons.
        tol = 1e-9 * max(np.max(np.abs(self.time), initial=0.0), 1.0)
        keep_time = np.ones(self.ntime, dtype=bool)
        keep_freq = np.ones(self.nfreq, dtype=bool)

        if start is not None:
            keep_time &= self.time >= start - tol
        if end is not None:
            keep_time &= self.time <= end + tol
        if fmin is not None:
            keep_freq &= self.frequency >= fmin
        if fmax is not None:
            keep_freq &= self.frequency <= fmax

        if not keep_time.any() or not keep_freq.any():
            raise ValueError(
                "no coefficient falls within the requested window"
            )

        target = self if inplace else self.copy()
        target.data = self.data[np.ix_(keep_freq, keep_time)]
        target.frequency = self.frequency[keep_freq]
        target.time = self.time[keep_time]

        if self.coi is not None:
            target.coi = self.coi[keep_freq]

        return None if inplace else target

    def _relative_time(self, value: Any, name: str) -> Optional[float]:
        """Return a time in seconds from the start of the record."""
        if value is None:
            return None

        if isinstance(value, (Date, str)):
            if self.head.time is None:
                raise ValueError(
                    f"{name} as a date requires the record start time"
                )

            return Date(value) - self.head.time

        if isinstance(value, Real) and not isinstance(value, bool):
            return float(value)

        raise TypeError(f"{name} must be a Date, a string or a number")


def _frequency_grid(
    fmin: float,
    fmax: float,
    nfreq: int,
    spacing: str,
) -> np.ndarray:
    """Return ``nfreq`` frequencies from ``fmin`` to ``fmax``."""
    if fmax < fmin:
        raise ValueError("fmax must not be smaller than fmin")

    if nfreq == 1:
        return np.array([fmin])

    if spacing == "log":
        return np.geomspace(fmin, fmax, nfreq)

    return np.linspace(fmin, fmax, nfreq)
