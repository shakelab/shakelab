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

from shakelab.signals import base
from shakelab.libutils.time import Date


def fft(data):
    """
    Note: better using numpy.fft or scipy.fftpack?
    """
    return np.fft.fft(data) * (2/len(data))

def ifft(data):
    """
    """
    return np.fft.ifft(data) / (2/len(data))

def fft_axis_double_sided(snum, dt):
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

def fft_axis_single_sided(snum, dt):
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
    shift = ifft(fft(signal)*expt)

    return np.real(shift)


class Spectrum():
    """
    Fourier spectrum base class
    """

    def __init__(self):
        self.data = []
        self.time = Date()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, sliced):
        return self.data[sliced]

    @property
    def amplitude(self):
        """
        """
        return np.abs(self.data)

    @property
    def phase(self, unwrap=False):
        """
        """
        phase = np.angle(self.data)
        if unwrap:
            return np.unwrap(phase)
        else:
            return phase


class DiscreteSpectrum(Spectrum):
    """
    Discrete Fourier spectrum
    """

    def __init__(self):
        pass


class FFTSpectrum(Spectrum):
    """
    Fourier spectrum in FFT format
    """

    def __init__(self, record=None):
        self.df = None

        if record is not None:
            self.fft(record)

    def fft(self, record, norm=False):
        """
        """
        self.df = 1./(record.dt * len(record))
        self.data = fft(record.data)
        if norm:
            self.data = self.data / record.dt
        self.time = record.time

    def ifft(self, norm=False):
        """
        """
        record = base.Record()
        record.dt = 1./(self.df * len(self))
        record.data = np.real(ifft(self.data))
        if norm:
            record.data = record.data / self.df
        record.time = self.time
        return record



