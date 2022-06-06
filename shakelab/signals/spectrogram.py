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
Collection of algorithms to compute time-frequency spectrograms.
"""

import numpy as np
from numpy.matlib import repmat
from scipy.linalg import toeplitz


def wft(signal, delta):
    """
    """
    pass

def cwt(signal, delta, wavelet):
    """
    """
    pass

def stransform(signal, delta):

    # Compute S-Transform without for loops
    # Coded by Kalyan S. Dash
    # IIT Bhubaneswar, India

    signal = np.reshape(signal, (1, -1))

    (_, n) = signal.shape

    nhaf = int(np.fix(n/2))

    odvn = 1
    if nhaf*2 == n:
        odvn = 0

    f_l = np.arange(0, nhaf + 1)/(delta*n)
    f_r = np.arange(-nhaf + 1 - odvn, 0)/(delta*n)
    f = np.concatenate([f_l, f_r])

    hft = np.fft.fft(signal)

    # Compute all frequency domain Gaussians as one matrix
    invfk = np.array([1.0 / f[1 : nhaf + 1]]).T
    w = 2 * np.pi * repmat(f, nhaf, 1) * repmat(invfk, 1, n)
    gw = np.exp((-w**2)/2)

    # Compute Toeplitz matrix with the shifted fft(h)
    hw = toeplitz(np.conj(hft.T[0 : nhaf + 1]), hft)

    # Exclude the first row, corresponding to zero frequency
    hw = hw[1 : nhaf + 1, :]

    # Compute Stockwell Transform
    st = np.fft.ifft(hw*gw)

    # Add the zero freq row
    st0 = np.mean(signal) * np.ones((1, n))
    st = np.concatenate([st0, st])

    return f_l, st

def istransform(h, delta):
    """
    Not yet implemented
    """
    pass
