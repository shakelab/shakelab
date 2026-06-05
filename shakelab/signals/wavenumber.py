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
Frequency-wavenumber analysis tools for DAS data.

This module provides functions for frequency-wavenumber analysis of DAS
records arranged as a 2-D matrix with shape
``(n_channels, n_samples)``.

The implemented methods combine Fourier, wavelet and spatial wavelet
transforms to analyse the temporal and spatial spectral content of DAS
signals at different levels of localisation.

Implemented transforms
----------------------
1) Fourier(f) --> Fourier(k)
   Classical frequency-wavenumber spectrum:
   ``f, k``

2) Fourier(f) --> Wavelet(k)
   Localised spatial wavenumber spectrum:
   ``f, x, k``

3) Wavelet(f) --> Fourier(k)
   Time-frequency-wavenumber spectrum:
   ``t, f, k``

4) Wavelet(f) --> Wavelet(k)
   Fully localised time-frequency-space-wavenumber spectrum:
   ``t, f, x, k``

The module currently uses complex Morlet wavelets for both temporal and
spatial continuous wavelet transforms.
"""

import numpy as np
from scipy import signal

from shakelab.signals.spectrogram import morlet_wavelet
from shakelab.signals.spectrogram import cwt_transform, _cwt_width_transform


def fk_transform(data, dx, dt, power=False, shift=True, angular=False):
    """
    Compute the 2-D Fourier spectrum in frequency-wavenumber domain.

    Parameters
    ----------
    data : array_like
        Input DAS data matrix with shape ``(n_channels, n_samples)``.
    dx : float
        Spatial sampling interval between channels, in metres.
    dt : float
        Temporal sampling interval, in seconds.
    power : bool, optional
        If True, return the squared amplitude spectrum. If False, return
        the complex Fourier spectrum. Default is False.
    shift : bool, optional
        If True, shift the wavenumber axis so that zero wavenumber is at
        the centre. Default is True.
    angular : bool, optional
        If True, return angular wavenumber in rad/m. If False, return
        cyclic wavenumber in cycles/m. Default is False.

    Returns
    -------
    freq : ndarray
        Frequency axis, in Hz.
    wavenumber : ndarray
        Wavenumber axis, in cycles/m or rad/m.
    spectrum : ndarray
        Frequency-wavenumber spectrum with shape
        ``(n_k, n_freq)``.

    Notes
    -----
    The transform is computed using a real FFT along the time axis and a
    complex FFT along the spatial/channel axis.
    """
    data = np.asarray(data)

    if data.ndim != 2:
        raise ValueError("data must be a 2-D array")

    if dx <= 0:
        raise ValueError("dx must be positive")

    if dt <= 0:
        raise ValueError("dt must be positive")

    n_channels, n_samples = data.shape

    freq = np.fft.rfftfreq(n_samples, d=dt)
    wavenumber = np.fft.fftfreq(n_channels, d=dx)

    spectrum = np.fft.rfft(data, axis=1)
    spectrum = np.fft.fft(spectrum, axis=0)

    if angular:
        wavenumber = 2.0 * np.pi * wavenumber

    if shift:
        wavenumber = np.fft.fftshift(wavenumber)
        spectrum = np.fft.fftshift(spectrum, axes=0)

    if power:
        spectrum = np.abs(spectrum) ** 2

    return freq, wavenumber, spectrum


def fk_cwt_transform(
        data, dx, dt, wavenumber=None, n_k=64, omega0=6.0,
        include_negative=True, freq_min=None, freq_max=None,
        power=False):
    """
    Compute a frequency-space-wavenumber spectrum.

    The function first applies a real Fourier transform along the time
    axis. Then, for each selected temporal frequency, it applies a
    complex Morlet continuous wavelet transform along the spatial axis.

    Parameters
    ----------
    data : array_like
        Input DAS data matrix with shape ``(n_channels, n_samples)``.
    dx : float
        Spatial sampling interval between channels, in metres.
    dt : float
        Temporal sampling interval, in seconds.
    wavenumber : array_like, optional
        Target cyclic wavenumbers, in cycles/m. Zero is not allowed.
        If None, logarithmically spaced positive wavenumbers are used.
    n_k : int, optional
        Number of positive wavenumbers used when ``wavenumber`` is None.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.
    include_negative : bool, optional
        If True, include both positive and negative wavenumbers.
    freq_min : float, optional
        Minimum temporal frequency to include, in Hz.
    freq_max : float, optional
        Maximum temporal frequency to include, in Hz.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.

    Returns
    -------
    freq : ndarray
        Selected temporal frequencies, in Hz.
    x : ndarray
        Spatial coordinate of the channels, in metres.
    k : ndarray
        Cyclic wavenumbers, in cycles/m.
    spectrum : ndarray
        Spectrum with shape ``(n_freq, n_k, n_channels)``.

    Notes
    -----
    The output can be interpreted as ``spectrum[f, k, x]``.
    """
    data = np.asarray(data)

    if data.ndim != 2:
        raise ValueError("data must be a 2-D array")

    if dx <= 0:
        raise ValueError("dx must be positive")

    if dt <= 0:
        raise ValueError("dt must be positive")

    if n_k <= 0:
        raise ValueError("n_k must be positive")

    n_channels, n_samples = data.shape

    freq = np.fft.rfftfreq(n_samples, d=dt)
    spectrum_f = np.fft.rfft(data, axis=1)

    freq_mask = np.ones(freq.shape, dtype=bool)

    if freq_min is not None:
        freq_mask &= freq >= freq_min

    if freq_max is not None:
        freq_mask &= freq <= freq_max

    freq = freq[freq_mask]
    spectrum_f = spectrum_f[:, freq_mask]

    if wavenumber is None:
        k_min = 1.0 / (n_channels * dx)
        k_max = 1.0 / (2.0 * dx)

        k_pos = np.logspace(
            np.log10(k_min),
            np.log10(k_max),
            n_k,
        )

        if include_negative:
            k = np.concatenate((-k_pos[::-1], k_pos))
        else:
            k = k_pos
    else:
        k = np.asarray(wavenumber, dtype=float)

    if np.any(k == 0.0):
        raise ValueError("wavenumber values must be different from zero")

    x = np.arange(n_channels) * dx

    output = np.zeros(
        (len(freq), len(k), n_channels),
        dtype=np.complex128,
    )

    for ik, kval in enumerate(k):
        width = omega0 / (2.0 * np.pi * abs(kval) * dx)
        wavelet_omega0 = np.sign(kval) * omega0

        for ifreq in range(len(freq)):
            spatial_signal = spectrum_f[:, ifreq]

            coeff = spatial_cwt(
                spatial_signal,
                [width],
                omega0=wavelet_omega0,
            )

            output[ifreq, ik, :] = coeff[0]

    if power:
        output = np.abs(output) ** 2

    return freq, x, k, output


def cwt_fk_transform(
        data, dx, dt, freq=None, n_freq=64, freq_min=None,
        freq_max=None, omega0=6.0, power=False, shift=True,
        angular=False):
    """
    Compute a time-frequency-wavenumber spectrum.

    The function first applies a Morlet continuous wavelet transform
    along the time axis for each channel. Then, for each frequency and
    time sample, it applies a Fourier transform along the spatial axis.

    Parameters
    ----------
    data : array_like
        Input DAS data matrix with shape ``(n_channels, n_samples)``.
    dx : float
        Spatial sampling interval between channels, in metres.
    dt : float
        Temporal sampling interval, in seconds.
    freq : array_like, optional
        Frequencies to analyse, in Hz. If None, frequencies are generated
        automatically between ``freq_min`` and ``freq_max``.
    n_freq : int, optional
        Number of frequencies used when ``freq`` is None.
    freq_min : float, optional
        Minimum frequency, in Hz.
    freq_max : float, optional
        Maximum frequency, in Hz.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.
    power : bool, optional
        If True, return squared amplitudes. If False, return complex
        coefficients.
    shift : bool, optional
        If True, shift the wavenumber axis so that zero wavenumber is at
        the centre.
    angular : bool, optional
        If True, return angular wavenumber in rad/m. If False, return
        cyclic wavenumber in cycles/m.

    Returns
    -------
    freq : ndarray
        Analysed frequencies, in Hz.
    k : ndarray
        Wavenumber axis, in cycles/m or rad/m.
    time : ndarray
        Time axis, in seconds.
    spectrum : ndarray
        Spectrum with shape ``(n_freq, n_k, n_samples)``.

    Notes
    -----
    The output can be interpreted as ``spectrum[f, k, t]``.
    """
    data = np.asarray(data)

    if data.ndim != 2:
        raise ValueError("data must be a 2-D array")

    if dx <= 0:
        raise ValueError("dx must be positive")

    if dt <= 0:
        raise ValueError("dt must be positive")

    if n_freq <= 0:
        raise ValueError("n_freq must be positive")

    n_channels, n_samples = data.shape

    if freq is None:
        if freq_min is None:
            freq_min = 1.0 / (n_samples * dt)

        if freq_max is None:
            freq_max = 1.0 / (2.0 * dt)

        freq = np.logspace(
            np.log10(freq_min),
            np.log10(freq_max),
            n_freq,
        )
    else:
        freq = np.asarray(freq, dtype=float)

    if np.any(freq <= 0):
        raise ValueError("all frequencies must be positive")

    time = np.arange(n_samples) * dt
    k = np.fft.fftfreq(n_channels, d=dx)

    tf_data = time_cwt(
        data,
        dt=dt,
        freq=freq,
        omega0=omega0,
    )

    spectrum = np.fft.fft(tf_data, axis=1)

    if angular:
        k = 2.0 * np.pi * k

    if shift:
        k = np.fft.fftshift(k)
        spectrum = np.fft.fftshift(spectrum, axes=1)

    if power:
        spectrum = np.abs(spectrum) ** 2

    return freq, k, time, spectrum


def cwt_cwt_transform(
        data, dx, dt, freq=None, wavenumber=None, n_freq=32,
        n_k=32, freq_min=None, freq_max=None, omega0_time=6.0,
        omega0_space=6.0, include_negative=True, power=False):
    """
    Compute a time-frequency-space-wavenumber spectrum.

    The function first applies a Morlet CWT along the time axis for each
    channel. Then, for each temporal frequency and time sample, it applies
    a Morlet CWT along the spatial axis.

    Parameters
    ----------
    data : array_like
        Input DAS data matrix with shape ``(n_channels, n_samples)``.
    dx : float
        Spatial sampling interval between channels, in metres.
    dt : float
        Temporal sampling interval, in seconds.
    freq : array_like, optional
        Frequencies to analyse, in Hz. If None, they are generated
        automatically.
    wavenumber : array_like, optional
        Cyclic wavenumbers to analyse, in cycles/m. If None, they are
        generated automatically.
    n_freq : int, optional
        Number of frequencies used when ``freq`` is None.
    n_k : int, optional
        Number of positive wavenumbers used when ``wavenumber`` is None.
    freq_min : float, optional
        Minimum temporal frequency, in Hz.
    freq_max : float, optional
        Maximum temporal frequency, in Hz.
    omega0_time : float, optional
        Morlet central frequency for the temporal CWT.
    omega0_space : float, optional
        Morlet central frequency for the spatial CWT.
    include_negative : bool, optional
        If True, include positive and negative wavenumbers.
    power : bool, optional
        If True, return squared amplitudes.

    Returns
    -------
    freq : ndarray
        Analysed temporal frequencies, in Hz.
    x : ndarray
        Spatial coordinate of the channels, in metres.
    k : ndarray
        Analysed cyclic wavenumbers, in cycles/m.
    time : ndarray
        Time axis, in seconds.
    spectrum : ndarray
        Spectrum with shape ``(n_freq, n_k, n_channels, n_samples)``.

    Notes
    -----
    The output can be interpreted as ``spectrum[f, k, x, t]`` where:
    - ``f`` is temporal frequency,
    - ``k`` is spatial wavenumber,
    - ``x`` is channel position,
    - ``t`` is time.
    """
    data = np.asarray(data)

    if data.ndim != 2:
        raise ValueError("data must be a 2-D array")

    if dx <= 0:
        raise ValueError("dx must be positive")

    if dt <= 0:
        raise ValueError("dt must be positive")

    if n_freq <= 0:
        raise ValueError("n_freq must be positive")

    if n_k <= 0:
        raise ValueError("n_k must be positive")

    n_channels, n_samples = data.shape

    if freq is None:
        if freq_min is None:
            freq_min = 1.0 / (n_samples * dt)

        if freq_max is None:
            freq_max = 1.0 / (2.0 * dt)

        freq = np.logspace(
            np.log10(freq_min),
            np.log10(freq_max),
            n_freq,
        )
    else:
        freq = np.asarray(freq, dtype=float)

    if np.any(freq <= 0):
        raise ValueError("all frequencies must be positive")

    if wavenumber is None:
        k_min = 1.0 / (n_channels * dx)
        k_max = 1.0 / (2.0 * dx)

        k_pos = np.logspace(
            np.log10(k_min),
            np.log10(k_max),
            n_k,
        )

        if include_negative:
            k = np.concatenate((-k_pos[::-1], k_pos))
        else:
            k = k_pos
    else:
        k = np.asarray(wavenumber, dtype=float)

    if np.any(k == 0.0):
        raise ValueError("wavenumber values must be different from zero")

    x = np.arange(n_channels) * dx
    time = np.arange(n_samples) * dt

    tf_data = time_cwt(
        data,
        dt=dt,
        freq=freq,
        omega0=omega0_time,
    )

    output = np.zeros(
        (len(freq), len(k), n_channels, n_samples),
        dtype=np.complex128,
    )

    for ik, kval in enumerate(k):
        width = omega0_space / (
            2.0 * np.pi * abs(kval) * dx
        )

        wavelet_omega0 = np.sign(kval) * omega0_space

        for ifreq in range(len(freq)):
            coeff = spatial_cwt_matrix(
                tf_data[ifreq],
                [width],
                omega0=wavelet_omega0,
            )

            output[ifreq, ik, :, :] = coeff[0]

    if power:
        output = np.abs(output) ** 2

    return freq, x, k, time, output


def time_cwt(data, dt, freq, omega0=6.0):
    """
    Compute a Morlet CWT along the time axis.

    Parameters
    ----------
    data : array_like
        Input data matrix with shape ``(n_channels, n_samples)``.
    dt : float
        Temporal sampling interval, in seconds.
    freq : array_like
        Frequencies to analyse, in Hz.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.

    Returns
    -------
    output : ndarray
        Complex CWT coefficients with shape
        ``(n_freq, n_channels, n_samples)``.
    """
    _, output = cwt_transform(
        data,
        delta=dt,
        freq=freq,
        omega0=omega0,
        axis=1,
        power=False,
    )

    return output


def spatial_cwt(data, widths, omega0=6.0):
    """
    Compute a continuous wavelet transform along one spatial axis.

    Parameters
    ----------
    data : array_like
        One-dimensional spatial signal.
    widths : array_like
        Wavelet scales, expressed in samples.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.

    Returns
    -------
    output : ndarray
        Complex CWT coefficients with shape ``(n_widths, n_samples)``.
    """
    _, output = _cwt_width_transform(
        data,
        widths,
        omega0=omega0,
        axis=0,
        power=False,
    )

    return output


def spatial_cwt_matrix(data, widths, omega0=6.0):
    """
    Compute a spatial continuous wavelet transform on a 2-D matrix.

    Parameters
    ----------
    data : array_like
        Input matrix with shape ``(n_channels, n_samples)``.
    widths : array_like
        Wavelet scales, expressed in samples.
    omega0 : float, optional
        Non-dimensional central frequency of the Morlet wavelet.

    Returns
    -------
    output : ndarray
        Complex CWT coefficients with shape
        ``(n_widths, n_channels, n_samples)``.

    Notes
    -----
    The transform is applied along the spatial/channel axis while
    preserving the temporal dimension.
    """
    _, output = _cwt_width_transform(
        data,
        widths,
        omega0=omega0,
        axis=0,
        power=False,
    )

    return output

