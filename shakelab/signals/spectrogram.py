# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Time-frequency analysis tools.

This module provides basic algorithms to compute time-frequency
representations of one-dimensional signals or arrays along a selected
axis.

Implemented transforms
----------------------
wft_transform
    Windowed Fourier transform.

cwt_transform
    Continuous wavelet transform using a complex Morlet wavelet.

st_transform
    Stockwell transform.

The functions are designed as generic building blocks that can also be
used by higher-level modules, such as frequency-wavenumber analysis tools.
"""

import numpy as np
from scipy import signal


def morlet_wavelet(length, width, omega0=6.0):
    """
    Generate a complex Morlet wavelet.

    Parameters
    ----------
    length : int
        Number of samples of the wavelet.
    width : float
        Wavelet scale, expressed in samples.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.

    Returns
    -------
    wavelet : ndarray
        Complex Morlet wavelet with unit energy.

    Notes
    -----
    The wavelet is normalised to unit energy.
    """
    if length <= 0:
        raise ValueError("length must be positive")

    if width <= 0:
        raise ValueError("width must be positive")

    x = np.arange(length) - (length - 1) / 2.0

    wavelet = np.exp(1j * omega0 * x / width)
    wavelet *= np.exp(-(x ** 2) / (2.0 * width ** 2))

    norm = np.sqrt(np.sum(np.abs(wavelet) ** 2))

    return wavelet / norm


def wft_transform(
        data, delta, window_length, step=None, window="hann",
        axis=-1, power=False):
    """
    Compute a windowed Fourier transform.

    Parameters
    ----------
    data : array_like
        Input data array.
    delta : float
        Sampling interval along the selected axis.
    window_length : int
        Length of the moving window, in samples.
    step : int, optional
        Step between consecutive windows, in samples. If None, half
        window overlap is used.
    window : str or tuple or array_like, optional
        Window specification passed to ``scipy.signal.get_window``. If an
        array is provided, it must have length ``window_length``.
    axis : int, optional
        Axis along which the transform is computed.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.

    Returns
    -------
    freq : ndarray
        Frequency axis, in Hz.
    time : ndarray
        Window centre times, in seconds.
    spectrum : ndarray
        Windowed Fourier spectrum. The first two axes are frequency and
        window index; the remaining axes correspond to the non-transformed
        dimensions of the input array.

    Notes
    -----
    The selected axis is moved internally to the last position. The
    returned spectrum has shape ``(n_freq, n_windows, ...)``.
    """
    data = np.asarray(data)

    if delta <= 0:
        raise ValueError("delta must be positive")

    if window_length <= 0:
        raise ValueError("window_length must be positive")

    if step is None:
        step = window_length // 2

    if step <= 0:
        raise ValueError("step must be positive")

    data = np.moveaxis(data, axis, -1)
    n_samples = data.shape[-1]

    if window_length > n_samples:
        raise ValueError("window_length exceeds selected axis length")

    starts = np.arange(0, n_samples - window_length + 1, step)
    centres = starts + (window_length - 1) / 2.0

    freq = np.fft.rfftfreq(window_length, d=delta)
    time = centres * delta

    if np.ndim(window) == 0:
        taper = signal.get_window(window, window_length)
    else:
        taper = np.asarray(window)

        if taper.shape[0] != window_length:
            raise ValueError("window length does not match window_length")

    out_shape = (
        len(freq),
        len(starts),
    ) + data.shape[:-1]

    spectrum = np.zeros(out_shape, dtype=np.complex128)

    for index, start in enumerate(starts):
        stop = start + window_length
        segment = data[..., start:stop] * taper
        coeff = np.fft.rfft(segment, axis=-1)
        spectrum[:, index, ...] = np.moveaxis(coeff, -1, 0)

    if power:
        spectrum = np.abs(spectrum) ** 2

    return freq, time, spectrum


def _cwt_width_transform(
        data, widths, omega0=6.0, axis=-1, power=False):
    """
    Compute a continuous wavelet transform from wavelet widths.

    Parameters
    ----------
    data : array_like
        Input data array.
    widths : array_like
        Wavelet scales, expressed in samples.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.
    axis : int, optional
        Axis along which the transform is computed.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.

    Returns
    -------
    widths : ndarray
        Wavelet scales, expressed in samples.
    transform : ndarray
        CWT coefficients. The first axis corresponds to width; the
        remaining axes preserve the original input array shape.
    """
    data = np.asarray(data)
    widths = np.asarray(widths, dtype=float)

    if np.any(widths <= 0):
        raise ValueError("all widths must be positive")

    axis = np.core.numeric.normalize_axis_index(axis, data.ndim)

    moved = np.moveaxis(data, axis, -1)
    n_samples = moved.shape[-1]

    out_shape = (len(widths),) + moved.shape
    transform = np.zeros(out_shape, dtype=np.complex128)

    for iwidth, width in enumerate(widths):
        wavelet = morlet_wavelet(
            n_samples,
            width,
            omega0=omega0,
        )

        kernel = np.conj(wavelet[::-1])

        kernel_shape = [1] * moved.ndim
        kernel_shape[-1] = len(kernel)
        kernel = kernel.reshape(kernel_shape)

        transform[iwidth] = signal.fftconvolve(
            moved,
            kernel,
            mode="same",
            axes=-1,
        )

    transform = np.moveaxis(transform, -1, axis + 1)

    if power:
        transform = np.abs(transform) ** 2

    return widths, transform


def cwt_transform(
        data, delta, freq, omega0=6.0, axis=-1,
        power=False):
    """
    Compute a continuous wavelet transform with a Morlet wavelet.

    Parameters
    ----------
    data : array_like
        Input data array.
    delta : float
        Sampling interval along the selected axis.
    freq : array_like
        Frequencies to analyse, in Hz.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.
    axis : int, optional
        Axis along which the transform is computed.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.

    Returns
    -------
    freq : ndarray
        Analysed frequencies, in Hz.
    transform : ndarray
        Complex CWT coefficients. The first axis corresponds to
        frequency; the remaining axes follow the input array, with the
        transformed axis preserved as the last dimension.

    Notes
    -----
    The transform is internally computed from wavelet widths derived
    from the analysed frequencies:

    :contentReference[oaicite:0]{index=0}

    where:
    - ``a`` is the wavelet width in samples,
    - ``omega0`` is the Morlet central frequency,
    - ``f`` is the analysed frequency,
    - ``delta`` is the sampling interval.

    The selected axis is moved internally to the last position. The
    output has shape ``(n_freq, ..., n_samples)``.
    """
    if delta <= 0:
        raise ValueError("delta must be positive")

    freq = np.asarray(freq, dtype=float)

    if np.any(freq <= 0):
        raise ValueError("all frequencies must be positive")

    widths = omega0 / (
        2.0 * np.pi * freq * delta
    )

    _, transform = _cwt_width_transform(
        data,
        widths,
        omega0=omega0,
        axis=axis,
        power=power,
    )

    return freq, transform


def st_transform(data, delta, axis=-1, power=False):
    """
    Compute the Stockwell transform.

    Parameters
    ----------
    data : array_like
        Input data array.
    delta : float
        Sampling interval along the selected axis.
    axis : int, optional
        Axis along which the transform is computed.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.

    Returns
    -------
    freq : ndarray
        Positive frequency axis, in Hz.
    transform : ndarray
        Stockwell coefficients. The first axis corresponds to frequency;
        the remaining axes follow the input array, with the transformed
        axis preserved as the last dimension.

    Notes
    -----
    This implementation computes the one-sided Stockwell transform,
    including the zero-frequency row.
    """
    data = np.asarray(data)

    if delta <= 0:
        raise ValueError("delta must be positive")

    data = np.moveaxis(data, axis, -1)
    original_shape = data.shape
    n_samples = original_shape[-1]

    flat_data = data.reshape((-1, n_samples))
    n_series = flat_data.shape[0]

    n_half = int(np.fix(n_samples / 2))

    is_odd = 1
    if n_half * 2 == n_samples:
        is_odd = 0

    f_left = np.arange(0, n_half + 1) / (delta * n_samples)
    f_right = np.arange(-n_half + 1 - is_odd, 0) / (
        delta * n_samples
    )

    freq_full = np.concatenate([f_left, f_right])
    freq = f_left

    fft_data = np.fft.fft(flat_data, axis=1)

    transform = np.zeros(
        (len(freq), n_series, n_samples),
        dtype=np.complex128,
    )

    transform[0, :, :] = np.mean(flat_data, axis=1)[:, np.newaxis]

    if n_half > 0:
        pos_freq = freq_full[1:n_half + 1]

        freq_matrix = freq_full[np.newaxis, :]
        inv_freq = (1.0 / pos_freq)[:, np.newaxis]

        omega = 2.0 * np.pi * freq_matrix * inv_freq
        gauss = np.exp(-(omega ** 2) / 2.0)

        for ifreq in range(1, n_half + 1):
            shifted = np.roll(fft_data, -ifreq, axis=1)
            windowed = shifted * gauss[ifreq - 1][np.newaxis, :]
            transform[ifreq, :, :] = np.fft.ifft(windowed, axis=1)

    output_shape = (len(freq),) + original_shape
    transform = transform.reshape(output_shape)

    if power:
        transform = np.abs(transform) ** 2

    return freq, transform


def power_to_db(power, eps=1e-12):
    """
    Convert a power spectrum to decibel scale.

    Parameters
    ----------
    power : array_like
        Input power spectrum.
    eps : float, optional
        Small stabilisation value added to avoid logarithm singularities.

    Returns
    -------
    output : ndarray
        Spectrum expressed in decibels.
    """
    power = np.asarray(power)

    return 10.0 * np.log10(power + eps)

