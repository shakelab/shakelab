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
# You should have received a copy of the GNU General Public License with
# this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Engineering signal parameters and response functions.
"""

import numpy as np
from scipy import signal
from scipy.integrate import trapezoid, cumulative_trapezoid

from shakelab.engineering.response import (
    sdof_response_spectrum,
    sdof_interdrift,
    newmark_integration
)
from shakelab.libutils.constants import PI, GRAVITY


def hilbert_transform(record):
    """
    Compute the analytic signal using the Hilbert transform.

    The analytic signal is a complex representation of the input
    signal where the real part is the original signal and the
    imaginary part is the Hilbert transform.

    Parameters
    ----------
    record : Record
        Input signal record.

    Returns
    -------
    ndarray
        Complex analytic signal, same length as input.
    """
    return signal.hilbert(record.data)


def amplitude_envelope(record, method='hilbert', window=0.1):
    """
    Compute the amplitude envelope of the signal.

    Parameters
    ----------
    record : Record
        Input signal record.
    method : str, default 'hilbert'
        Envelope computation method: 'hilbert' or 'rms'.
    window : float, optional
        Window size in seconds for RMS method.

    Returns
    -------
    ndarray
        Amplitude envelope of the signal.
    """
    if method == 'hilbert':
        return np.abs(hilbert_transform(record))

    elif method == 'rms':
        win_len = max(1, int(round(window / record.head.delta)))
        squared = record.data**2
        kernel = np.ones(win_len) / win_len
        smoothed = np.convolve(squared, kernel, mode='same')
        return np.sqrt(smoothed)

    else:
        raise ValueError("Unsupported method. Use 'hilbert' or 'rms'.")


def instantaneous_phase(record):
    """
    Compute the instantaneous phase of the signal from the Hilbert transform.

    The instantaneous phase is obtained by unwrapping the angle of the
    analytic (Hilbert) signal.

    Parameters
    ----------
    record : Record
        Input signal record.

    Returns
    -------
    ndarray
        Unwrapped instantaneous phase in radians.
    """
    analytic = hilbert_transform(record)
    return np.unwrap(np.angle(analytic))


def instantaneous_frequency(record, method='gradient', pad=False):
    """
    Compute the instantaneous frequency from the unwrapped phase.

    Parameters
    ----------
    record : Record
        Input signal record.
    method : str, default 'diff'
        Method to compute time derivative of the phase:
        - 'diff': first-order forward difference (N-1 samples)
        - 'gradient': centered finite difference (N samples)
    pad : bool, default False
        If True, pads the result to match the input length (for 'diff' method)

    Returns
    -------
    ndarray
        Instantaneous frequency in Hz.
    """
    phase = instantaneous_phase(record)

    if method == 'diff':
        freq = np.diff(phase) / (2 * PI) * record.head.rate
        if pad:
            freq = np.append(freq, freq[-1])

    elif method == 'gradient':
        freq = np.gradient(phase, record.head.delta) / (2 * PI)

    else:
        raise ValueError("Invalid method. Use 'diff' or 'gradient'.")

    return freq


def peak_amplitude(record):
    """
    Compute the peak absolute amplitude of the signal.

    Parameters
    ----------
    record : Record
        Input signal record.

    Returns
    -------
    float or tuple
        (Peak, time).
    """
    abs_data = np.abs(record.data)
    peak = np.max(abs_data)

    index = np.argmax(abs_data)
    time = record.starttime + index * record.delta
    return peak, time


def arias_intensity(record, normalize=False):
    """
    Compute the Arias Intensity of the input signal.

    Parameters
    ----------
    record : Record
        Input signal record (should be in acceleration).
    normalize : bool, default False
        If True, return the cumulative normalized Arias intensity.

    Returns
    -------
    float or ndarray
        Arias intensity (float), or normalized cumulative array if normalize.
    """
    acc_squared = record.data**2
    cumulative = cumulative_trapezoid(
        acc_squared,
        dx=record.head.delta,
        initial=0.0
        )
    total = cumulative[-1]
    scale = PI / (2 * GRAVITY)

    if normalize:
        return cumulative / total
    else:
        return scale * total


def cumulative_absolute_velocity(record, normalize=False):
    """
    Compute the Cumulative Absolute Velocity (CAV) of the signal.

    Parameters
    ----------
    record : Record
        Input signal record (usually acceleration).
    normalize : bool, default False
        If True, return cumulative normalized CAV array.

    Returns
    -------
    float or ndarray
        CAV value (float) or normalized cumulative array if normalize=True.
    """
    abs_data = np.abs(record.data)
    cumulative = cumulative_trapezoid(
        abs_data,
        dx=record.head.delta,
        initial=0.0
        )

    if normalize:
        return cumulative / cumulative[-1]
    else:
        return cumulative[-1]


def bracketed_duration(record, threshold=0.05, relative=False,
                       return_times=False):
    """
    Compute the bracketed duration of the signal.

    Parameters
    ----------
    record : Record
        Input signal record.
    threshold : float, default 0.05
        Absolute threshold, or fraction of peak if relative=True.
    relative : bool, default False
        If True, interpret threshold as a fraction of peak amplitude.
    return_times : bool, default False
        If True, return also the start and end time of the interval.

    Returns
    -------
    float or tuple
        Duration in seconds, or (duration, t0, t1) if return_times=True.
        Returns None if no values exceed the threshold.
    """
    abs_data = np.abs(record.data)

    if relative:
        threshold *= np.max(abs_data)

    idx = np.argwhere(abs_data >= threshold)

    if idx.size == 0:
        return None

    i0, i1 = idx[0][0], idx[-1][0]
    t0 = record.starttime + i0 * record.delta
    t1 = record.starttime + i1 * record.delta
    duration = t1 - t0

    if return_times:
        return duration, t0, t1
    else:
        return duration


def significant_duration(record, threshold=(0.05, 0.95),
                         return_times=False):
    """
    Compute the significant duration based on cumulative Arias intensity.

    Parameters
    ----------
    record : Record
        Input signal record.
    threshold : tuple of float, default (0.05, 0.95)
        Lower and upper normalized energy bounds (e.g., 5%–95%).
    return_times : bool, default False
        If True, also return the start and end time.

    Returns
    -------
    float or tuple
        Duration in seconds, or (duration, t0, t1) if return_times=True.
        Returns None if energy is zero or indices not found.
    """
    energy = record.data**2
    cumulative = cumulative_trapezoid(
        energy,
        dx=record.head.delta,
        initial=0.0
        )
    total = cumulative[-1]

    if total == 0:
        return None

    norm = cumulative / total

    try:
        i0 = np.argmax(norm >= threshold[0])
        i1 = np.argmax(norm >= threshold[1])
    except IndexError:
        return None

    t0 = record.starttime + i0 * record.delta
    t1 = record.starttime + i1 * record.delta
    duration = t1 - t0

    if return_times:
        return duration, t0, t1
    else:
        return duration


def root_mean_square(record, window=None):
    """
    Compute the root mean square (RMS) of the signal.

    Parameters
    ----------
    record : Record
        Input signal record.
    window : float or None, default None
        RMS window length in seconds. If None, compute global RMS.
        If > 0, return moving RMS with centered window.

    Returns
    -------
    float or ndarray
        RMS value (scalar if window=None, else moving RMS array).
    """
    data_sq = record.data**2

    if window is None:
        integral = trapezoid(data_sq, dx=record.head.delta)
        return np.sqrt(integral / record.duration)

    else:
        nwin = max(1, int(round(window / record.head.delta)))
        kernel = np.ones(nwin) / nwin
        smoothed = np.convolve(data_sq, kernel, mode='same')
        return np.sqrt(smoothed)


def sdof_peak_response_spectrum(record, periods, zeta=0.05):
    """
    Compute the SDOF peak response spectrum over a range of oscillator periods.

    For each period, computes the maximum values of displacement, velocity,
    acceleration, and their pseudo forms, assuming a linear elastic SDOF
    system.

    Parameters
    ----------
    record : Record
        Input signal record (should be in acceleration).
    periods : array_like
        Oscillator periods (in seconds).
    zeta : float, default 0.05
        Damping ratio (5% critical by default).

    Returns
    -------
    dict of ndarray
        Dictionary containing arrays of peak response values:
        - 'sd'  : maximum displacement [m]
        - 'sv'  : maximum velocity [m/s]
        - 'sa'  : maximum acceleration [m/s²]
        - 'psv' : pseudo velocity [m/s]
        - 'psa' : pseudo acceleration [m/s²]
    """
    periods = np.atleast_1d(periods)
    rssp = sdof_response_spectrum(
        record.data,
        record.head.delta,
        periods,
        zeta=zeta
    )
    return {
        'sd': rssp[0],
        'sv': rssp[1],
        'sa': rssp[2],
        'psv': rssp[3],
        'psa': rssp[4]
    }


def sdof_time_response(record, period, zeta=0.05):
    """
    Compute the time history response of a SDOF system to the input motion.

    Solves the dynamic equilibrium equation using Newmark-beta integration,
    returning displacement, velocity, and acceleration time series.

    Parameters
    ----------
    record : Record
        Input signal record (usually acceleration in m/s²).
    period : float
        Natural period of the oscillator (in seconds).
    zeta : float, default 0.05
        Damping ratio (dimensionless, default is 5%).

    Returns
    -------
    dict of ndarray
        Time series response of the SDOF system:
        - 'd' : displacement [m]
        - 'v' : velocity [m/s]
        - 'a' : acceleration [m/s²]
    """
    resp = newmark_integration(
        record.data, record.head.delta, period, zeta=zeta
    )
    return {
        'd': resp[0],
        'v': resp[1],
        'a': resp[2]
    }


def sdof_interstory_drift(record, period, zeta=0.05):
    """
    Compute the interstory drift time history for an equivalent 2DOF system.

    The drift is computed as the relative displacement between two
    mass levels of a simplified structure responding to the input motion.

    Parameters
    ----------
    record : Record
        Input signal record (acceleration, in m/s²).
    period : float
        Natural period of the equivalent SDOF system (in seconds).
    zeta : float, default 0.05
        Damping ratio (dimensionless, default is 5%).

    Returns
    -------
    ndarray
        Interstory drift time series [m].
    """
    return sdof_interdrift(
        record.data,
        record.head.delta,
        period,
        zeta=zeta
    )


class ResponseSpectrum():
    """
    Compute and store the peak response spectrum of a SDOF system.

    Parameters
    ----------
    record : Record
        The input signal, typically in acceleration [m/s²].
    damping : float, default 0.05
        Damping ratio (5% by default).

    Attributes
    ----------
    record : Record
        Original input record (used for PGA/PGV/PGD).
    damping : float
        Damping ratio used for the computation.
    periods : ndarray
        Oscillator periods [s].
    sd, sv, sa, psv, psa : ndarray
        Peak displacement, velocity, acceleration, pseudo-velocity
        and pseudo-acceleration.
    """

    def __init__(self, record, damping=0.05):
        self.record = record
        self.damping = damping
        self.periods = None
        self.sd = None
        self.sv = None
        self.sa = None
        self.psv = None
        self.psa = None

    def compute(self, periods):
        """
        Compute the peak response spectrum for given oscillator periods.

        Parameters
        ----------
        periods : array_like
            List or array of oscillator periods [s].
        """
        self.periods = np.atleast_1d(periods)

        sd, sv, sa, psv, psa = sdof_response_spectrum(
            self.record.data,
            self.record.head.delta,
            self.periods,
            zeta=self.damping
        )

        self.sd = sd
        self.sv = sv
        self.sa = sa
        self.psv = psv
        self.psa = psa

    def __getitem__(self, key):
        if key not in ('sd', 'sv', 'sa', 'psv', 'psa'):
            raise KeyError(f"{key} is not a valid spectral key.")
        return getattr(self, key)

    def to_dict(self):
        """
        Export the response spectrum as a dictionary.

        Returns
        -------
        dict
            Dictionary with keys: 'periods', 'sd', 'sv', 'sa', 'psv', 'psa'.
        """
        return {
            "periods": self.periods,
            "sd": self.sd,
            "sv": self.sv,
            "sa": self.sa,
            "psv": self.psv,
            "psa": self.psa
        }

    @property
    def pga(self):
        """
        Peak ground acceleration [m/s²].

        Returns
        -------
        float
            Maximum absolute value of input acceleration.
        """
        return np.max(np.abs(self.record.data))

    @property
    def pgv(self):
        """
        Peak ground velocity [m/s].

        Returns
        -------
        float
            Maximum absolute velocity obtained by integration.
        """
        vel = np.gradient(self.record.data, self.record.head.delta)
        return np.max(np.abs(vel))

    @property
    def pgd(self):
        """
        Peak ground displacement [m].

        Returns
        -------
        float
            Maximum absolute displacement by double integration.
        """
        disp = cumulative_trapezoid(
            self.record.data, dx=self.record.head.delta, initial=0.0
        )
        return np.max(np.abs(disp))

    def __repr__(self):
        n = len(self.periods) if self.periods is not None else 0
        return (f"<ResponseSpectrum: {n} periods, "
                f"damping={self.damping:.2f}>")


def convolve_1d_site(record, model1d, component='sh', angle=0.0):
    """
    Placeholder for 1D site response convolution.

    Parameters
    ----------
    record : Record
        Input signal record.
    model1d : object
        1D soil profile model.
    component : str, default 'sh'
        Wave type ('sh', 'sv', or 'p').
    angle : float, default 0.0
        Incidence angle in degrees.

    Returns
    -------
    NotImplementedError
        Currently not implemented.
    """
    raise NotImplementedError("soil1d_convolve is not yet implemented.")
