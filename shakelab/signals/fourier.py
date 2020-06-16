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

import numpy as np


def fft_axis(snum, dt):
    """
    Equivalent to numpy.fft.fftfreq()
    Only used as reference.
    """

    # Check for odd/even number of samples
    if np.mod(snum, 2) == 0:
        pax = np.arange(0, snum/2)
        nax = np.arange(-snum/2, 0)
    else:
        pax = np.arange(0, (snum+1)/2)
        nax = np.arange(-(snum-1)/2, 0)
        
    return np.concatenate((pax, nax))/(dt*snum)

def fft_positive_axis(snum, dt):
    """
    Compute the positive frequency axis of an fft.
    """

    # return np.arange(snum)*(1./(snum*dt))
    return np.linspace(0., (snum-1.)/(dt*snum), snum)

def shift_time(signal, dt, time):
    """
    Shift a signal in time by using fft-based circular convolution.
    No zero-padding is assumed.
    """

    frax = np.fft.fftfreq(len(signal), dt)
    expt = np.exp(-2*1j*np.pi*time*frax)
    shift = np.fft.ifft(np.fft.fft(signal)*expt)

    return np.real(shift)
