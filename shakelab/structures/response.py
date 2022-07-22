# ****************************************************************************
#
# Copyright (C) 2019-2022, ShakeLab Developers.
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

from shakelab.libutils.constants import PI
from shakelab.signals import fourier


def sdof_response_spectrum(accg, delta, periods, zeta=0.05, method='newmark'):
    """
    """
    tlen = len(periods)

    sd = np.zeros(tlen, dtype=float)
    sv = np.zeros(tlen, dtype=float)
    sa = np.zeros(tlen, dtype=float)
    psv = np.zeros(tlen, dtype=float)
    psa = np.zeros(tlen, dtype=float)

    integration = {'newmark' : newmark_integration,
                   'fourier' : fourier_integration}

    for j in range(tlen):

        d, v, a = integration[method](accg, delta, periods[j], zeta)

        sd[j] = _peak_value(d)
        sv[j] = _peak_value(v)

        omega0 = 2.*PI/periods[j]

        if omega0 == 0.:
            sa[j] = _peak_value(accg)
        else:
            sa[j] = _peak_value(a)

        psv[j] = sd[j]*omega0
        psa[j] = sd[j]*(omega0**2)

    return sd, sv, sa, psv, psa

def sdof_interdrift(accg, delta, period, zeta=0.05, norm=1.):
    """
    """
    top = newmark_integration(accg, delta, period, zeta=zeta)
    bot = newmark_integration(accg, delta, 1e3, zeta=zeta)

    return (top[0] - bot[0]) / norm

def newmark_integration(accg, delta, period, zeta=0.05, beta=0.25, gamma=0.5):
    """
    accg: input ground acceleration
    """
    alen = len(accg)

    d = np.zeros(alen, dtype=float)
    v = np.zeros(alen, dtype=float)
    a = np.zeros(alen, dtype=float)

    a[0] = accg[0]

    omega0 = 2.*PI/period

    K = omega0**2
    C = 2.*zeta*(K**0.5)

    B = 1./(beta*delta*delta) + (gamma*C)/(beta*delta)
    A = B + K

    E = 1./(beta*delta) + (gamma/beta-1.)*C
    G = 1./(2.*beta) + C*(delta/2.)*(gamma/beta-2.) - 1

    for t in range(alen-1):

        d[t+1] = (accg[t+1] + B*d[t] + E*v[t] + G*a[t])/A

        a[t+1] = (d[t+1] - d[t] - delta*v[t])/(beta*delta*delta)
        a[t+1] += a[t]*(1. - 1./(2.*beta))

        v[t+1] = v[t] + delta*a[t] + delta*gamma*(a[t+1] - a[t])

    return -d, -v, -a+accg

def fourier_integration(accg, delta, period, zeta=0.05, nc=1):
    """
    Fourier domain solution (circular convolution)
    nc: number trace length cycles padded with zero
    """

    alen = len(accg)

    # 0-padding to avoid circular reverberation
    trace = np.pad(accg, (0, nc*alen), 'constant')

    # Natural frequency
    omega0 = 2.*PI/period

    # Convolution axis
    omega = 2.*PI*fourier.frequency(len(accg), delta)

    # Input Fourier spectrum
    trace_fft = fourier.fft(accg);

    # Harmonic oscillator spectrum
    sdof_fft = sdof_transfer_function(omega, omega0, zeta)

    # Convolution
    d = fourier.ifft(trace_fft * sdof_fft)
    v = fourier.ifft(trace_fft * sdof_fft * 1j*omega)
    a = fourier.ifft(trace_fft * sdof_fft * -(omega**2))

    return -d[:alen], -v[:alen], -a[:alen]+accg[:alen]

def omegaD(omega0, zeta):
    """
    Damped natural frequency
    """
    return omega0*np.sqrt(1.-zeta**2)

def omega0(m, k):
    """
    Compute the natural angular frequency of a SDOF
    from given mass and elastic constant.
    """
    return np.sqrt(k / m)

def sdof_transfer_function(omega, omega0, zeta=0.05):
    """
    """
    return 1./(omega0**2 + 2.*1j*zeta*omega0*omega - omega**2)

def _peak_value(trace):
    """
    """
    return np.max(np.abs(trace))
