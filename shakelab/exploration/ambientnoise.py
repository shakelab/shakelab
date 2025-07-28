# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
"""
import numpy as np

from shakelab.signals.base import Record


def average_spectral_ratio(numerator, denominator,
                            window_length,
                            overlapping=0.0,
                            waterlevel=1e-6,
                            highpass=None, lowpass=None,
                            filter_order=4,
                            taper_time=None,
                            logsmooth_sigma=None,
                            return_individual=False):
    """
    Compute the average spectral ratio across multiple time windows using
    lognormal statistics in the complex domain.

    Parameters:
        numerator (Record): Signal at the numerator (e.g., top).
        denominator (Record): Signal at the denominator (e.g., base).
        window_length (float): Length of each window in seconds.
        overlapping (float): Fractional overlap between windows.
        waterlevel (float): Stabilization level for spectral division.
        highpass, lowpass (float): Optional filter corner frequencies (Hz).
        filter_order (int): Butterworth filter order.
        taper_time (float): Optional taper time (s).
        logsmooth_sigma (float): Optional log-smoothing of final spectrum.
        return_individual (bool): If True, return also list of
                                  individual spectral ratios.

    Returns:
        Spectrum: Complex log-mean spectrum (average).
        np.ndarray: Complex log-standard deviation spectrum.
        list[Spectrum] (optional): List of spectral ratios per window.
    """
    # Copy and pre-filter full signals before splitting
    num_filtered = numerator.copy()
    den_filtered = denominator.copy()

    if highpass is not None or lowpass is not None:
        num_filtered.filter(highpass=highpass, lowpass=lowpass,
                            order=filter_order)
        den_filtered.filter(highpass=highpass, lowpass=lowpass,
                            order=filter_order)

    # Split into time windows
    num_segments = split_record(num_filtered, window_length, overlapping)
    den_segments = split_record(den_filtered, window_length, overlapping)

    if len(num_segments) != len(den_segments):
        raise ValueError("Unequal number of segments")

    spectra = []
    for num_win, den_win in zip(num_segments, den_segments):
        spec = spectral_ratio(
            num_win, den_win,
            waterlevel=waterlevel,
            taper_time=taper_time
        )
        spectra.append(spec)

    spec_array = np.array([s.data for s in spectra])
    log_spec = np.log(spec_array)

    logmean = np.mean(log_spec, axis=0)
    logstd = np.std(log_spec, axis=0)

    mean_spec = spectra[0].copy()
    mean_spec.data = np.exp(logmean)

    if logsmooth_sigma is not None:
        mean_spec.logsmooth(sigma=logsmooth_sigma)

    if return_individual:
        return mean_spec, np.exp(logstd), spectra
    else:
        return mean_spec, np.exp(logstd)


def spectral_ratio(numerator, denominator,
                   waterlevel=1e-6,
                   highpass=None, lowpass=None,
                   filter_order=4,
                   taper_time=None,
                   logsmooth_sigma=None):
    """
    Compute the spectral ratio of two Record objects, with optional
    filtering, tapering and log-smoothing.

    Parameters:
        numerator (Record): Signal at the numerator (e.g., top).
        denominator (Record): Signal at the denominator (e.g., base).
        waterlevel (float): Minimum relative power in denominator spectrum.
        highpass (float): Optional highpass frequency (Hz).
        lowpass (float): Optional lowpass frequency (Hz).
        filter_order (int): Butterworth filter order (default=4).
        taper_time (float): Optional taper time in seconds.
        logsmooth_sigma (float): Std dev in log-freq for smoothing.

    Returns:
        Spectrum: Spectral ratio (complex), smoothed if requested.
    """
    if numerator.nsamp != denominator.nsamp:
        raise ValueError("Records must have the same number of samples")

    if numerator.delta != denominator.delta:
        raise ValueError("Sampling intervals must match")

    if numerator.time != denominator.time:
        raise ValueError("Records must be synchronized")

    num = numerator.copy()
    den = denominator.copy()

    if taper_time is not None:
        num.taper(taper_time)
        den.taper(taper_time)

    if highpass is not None or lowpass is not None:
        num.filter(highpass=highpass, lowpass=lowpass,
                   order=filter_order)
        den.filter(highpass=highpass, lowpass=lowpass,
                   order=filter_order)

    spec1 = num.to_spectrum()
    spec2 = den.to_spectrum()

    if waterlevel in [0., None]:
        ratio  = spec1 / spec2
    else:
        S = spec2.data
        S_star = np.conj(S)
        power = np.abs(S)**2
        epsilon = waterlevel * np.max(power)
        inv = S_star / (power + epsilon)
    
        ratio = spec1.copy()
        ratio.data = spec1.data * inv

    if logsmooth_sigma is not None:
        ratio.logsmooth(sigma=logsmooth_sigma)

    return ratio


def split_record(record, window_length, overlapping=0.0):
    """
    Split a Record object into multiple fixed-length Record windows.

    Parameters:
        record (Record): Input Record object.
        window_length (float): Window length in seconds.
        overlapping (float): Overlap ratio (0.0 to <1.0).

    Returns:
        list[Record]: List of Record segments.
    """
    dt = record.delta
    nsamp = record.nsamp

    start_indices = split_window(nsamp, dt, window_length, overlapping)
    win_size = int(window_length / dt)

    segments = []
    for i0 in start_indices:
        i1 = i0 + win_size
        data_segment = record.data[i0:i1]

        # Create a new record and copy header
        sub = Record(data=data_segment)
        sub.head = record.head.copy()
        sub.head._parent = sub  # update parent reference
        sub.head.time = record.time + i0 * dt

        segments.append(sub)

    return segments


def split_window(number_of_samples, dt, window_length, overlapping=0.0):
    """
    Split a time series into windows and return the start indices
    of each window.

    Parameters:
        number_of_samples (int): Total number of samples in the time series.
        dt (float): Sampling interval in seconds.
        window_length (float): Window length in seconds.
        overlapping (float): Overlap between windows (0.0 to <1.0).
        Default is 0.

    Returns:
        list[int]: Indices of the first sample of each window.
    """
    if not (0.0 <= overlapping < 1.0):
        raise ValueError("Overlapping must be in the range [0.0, 1.0)")

    win_size = int(window_length / dt)
    if win_size <= 0:
        raise ValueError("Window size must be at least one sample")

    step = int(win_size * (1.0 - overlapping))
    if step <= 0:
        raise ValueError("Step between windows must be positive")

    indices = list(range(0, number_of_samples - win_size + 1, step))
    return indices


def combine_spectral_directions(spectrum_x, spectrum_y, method='energy'):
    """
    Combine two orthogonal spectra along their dominant horizontal direction.

    Parameters
    ----------
    spectrum_x : Spectrum
        Spectrum of the component along X (e.g., HN).
    spectrum_y : Spectrum
        Spectrum of the component along Y (e.g., HE).
    method : str, optional
        Method to estimate dominant direction:
            - 'energy' : maximizes total energy (default)
            - 'cross'  : uses complex cross-correlation
            - 'phase'  : uses average phase angle

    Returns
    -------
    Spectrum
        Complex spectrum projected along the estimated dominant direction.
    float
        Estimated azimuth in radians (theta, counterclockwise from X).
    """
    Sx = spectrum_x.data
    Sy = spectrum_y.data

    if method == 'energy':
        # Eigenvector of 2D spectral covariance matrix
        cross = np.sum(Sx * Sy.conj())
        dx = np.sum(np.abs(Sx)**2)
        dy = np.sum(np.abs(Sy)**2)
        theta = 0.5 * np.arctan2(2 * np.real(cross), dx - dy)

    elif method == 'cross':
        # Direction from cross-spectral phase
        cross = Sx * Sy.conj()
        theta = np.angle(np.sum(cross))

    elif method == 'phase':
        # Direction from phase angle of imaginary parts
        theta = np.arctan2(
            np.sum(np.imag(Sy)),
            np.sum(np.imag(Sx))
        )

    else:
        raise ValueError(f"Unknown method '{method}'")

    # Project the signal onto the estimated direction
    combined = Sx * np.cos(theta) + Sy * np.sin(theta)

    result = spectrum_x.copy()
    result.data = combined

    return result, theta
