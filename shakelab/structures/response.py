# ****************************************************************************
#
# Copyright (C) 2019-2020, ShakeLab Developers.
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

from shakelab.libutils.constants import GRAVITY, PI

import numpy as np


def peak_value(trace):
    """
    """
    return np.max(np.abs(trace))


class ResponseSDOF():
    """
    """
    def __init__(self, period):
        self.period = period

    def compute_spectrum(self, trace, dt):
        """
        """

        tlen = len(self.period)
        self.sd = np.zeros(tlen)
        self.sv = np.zeros(tlen)
        self.sa = np.zeros(tlen)
        self.psv = np.zeros(tlen)
        self.psa = np.zeros(tlen)

        for j in range(tlen):

            omega = 2.*PI/self.period[j]

            d, v, a = newmark_integration(trace, dt, omega)
            # d, v, a = fft_integration(trace, dt, omega)

            self.sd[j] = peak_value(d)
            self.sv[j] = peak_value(v)
            if omega == 0.:
                self.sa[j] = peak_value(trace)
            else:
                self.sa[j] = peak_value(a)

            self.psv[j] = self.sd[j]*omega
            self.psa[j] = self.sd[j]*(omega**2)


def newmark_integration(trace, dt, omega, zeta=0.05, beta=0.25, gamma=0.5):
    """
    if input in cm/s2, output is cm (sd) and cm/s2 (sa, psa)
    """

    alen = len(trace)
    d = np.zeros(alen)
    v = np.zeros(alen)
    a = np.zeros(alen)
    z = np.zeros(alen)

    a[0] = trace[0]

    K = omega**2
    C = 2.*zeta*(K**0.5)

    B = 1./(beta*dt*dt) + (gamma*C)/(beta*dt)
    A = B + K

    E = 1./(beta*dt) + (gamma/beta-1.)*C
    G = 1./(2.*beta) + C*(dt/2.)*(gamma/beta-2.) - 1

    for t in range(alen-1):

        d[t+1] = (trace[t+1] + B*d[t] + E*v[t] + G*a[t])/A

        a[t+1] = (d[t+1] - d[t] - dt*v[t])/(beta*dt*dt)
        a[t+1] += a[t]*(1. - 1./(2.*beta))

        v[t+1] = v[t] + dt*a[t] + dt*gamma*(a[t+1] - a[t])

    z = a - trace

    return d, v, z


def fft_integration(trace, dt, omega0, zeta=0.05):

    alen = len(trace)

    # Padding (20 cycles of the period)
    npad = int(20.*2.*PI/(omega0*dt))
    if npad > alen:
        trace = np.pad(trace, (0, npad-alen))
        alen += (npad-alen)

    amax = np.max(np.abs(trace))

    # Input spectrum
    afft = np.fft.fft(trace);

    # Signal duration
    #time = (alen-1)*dt

    # Convolution axis
    # omegaT = 2.*PI*(np.arange(alen))/time
    omegaT = 2.*PI*np.fft.fftfreq(alen, dt)

    # Damped natural frequency
    # omegaD = omega0*np.sqrt(1.-zeta**2)

    # Harmonic oscillator
    sdof = 1./(omega0**2 + 2.*1j*zeta*omega0*omegaT - omegaT**2)

    omegaT = 2.*PI*fft_positive_axis(alen, dt)

    # Convolution
    d = np.fft.ifft(afft*sdof)
    v = np.fft.ifft(afft*sdof*omegaT)
    a = np.fft.ifft(afft*sdof*omegaT**2) - amax

    return np.real(d), np.real(v), np.real(a)


def sdof_transfer_function(fmax, fnum, omega0, zeta=0.05):
    """
    """
    omega = np.linspace(0., fmax, fnum)

    tf = 1./(omega0**2 + 2.*1j*zeta*omega0*omega - omega**2)

    omega = np.linspace(0., 2*fmax, 2*fnum)
    tf = np.concatenate((tf, np.flip(tf)))

    return tf
    
