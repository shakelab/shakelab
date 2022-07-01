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
from scipy import interpolate
from copy import deepcopy

import shakelab.signals.base as base

def _halflen(snum):
    """
    Computing half-length of FFT, accounting for
    odd/even number of samples
    """
    return np.rint(snum/2)

def fft(data):
    """
    Amplitude scaling should be optional...
    """
    return np.fft.rfft(data) * 2 / len(data)

def ifft(data):
    """
    """
    return np.fft.irfft(data) * len(data)

def frequency(snum, delta):
    """
    """
    return np.fft.rfftfreq(snum, delta)

def shift_time(signal, delta, shift):
    """
    Shift a signal in time by using fft-based circular convolution.
    No zero-padding is assumed.
    """
    freq = frequency(len(signal), delta)
    expt = np.exp(-2*1j*np.pi*shift*freq)

    return ifft(fft(signal)*expt)


class Spectrum():
    """
    Discrete Fourier spectrum class in FFT format
    """

    def __init__(self, record=None):
        self.head = base.Header()
        self.data = np.array([], dtype="complex")

        if record is not None:
            self.fft(record)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, sliced):
        return self.data[sliced]

    def fft(self, record, norm=False):
        """
        """
        self.head = deepcopy(record.head)
        self.data = fft(record.data)

        if norm:
            self.data = self.data / record.dt

    def ifft(self, norm=False):
        """
        """
        record = base.Record()
        record.head = deepcopy(self.head)
        record.data = ifft(self.data)

        if norm:
            record.data = record.data / self.df

        return record

    @property
    def df(self):
        """
        """
        return 1/(self.head.delta * (len(self)-1))

    @property
    def amplitude(self):
        """
        """
        return np.abs(self.data)

    @property
    def phase(self):
        """
        """
        phase = np.angle(self.data)

    @property
    def unwrap(self):
        """
        """
        return np.unwrap(self.phase)

    @property
    def frequency(self):
        """
        """
        return frequency(2*len(self)-1, self.head.delta)

    def filter(self, highpass=None, lowpass=None):
        """
        """
        if (highpass is not None):
            self.data[self.frequency < highpass] = 0.

        if (lowpass is not None):
            self.data[self.frequency > lowpass] = 0.

    def resample(self, frequency):
        """
        Extract spectrum samples at specific frequencies.
        """
        f = interpolate.interp1d(self.frequency, self.data)
        return f(frequency)

    def logsmooth(self, sigma, memsafe=False):
        """
        Logarithm smoothing of (complex) spectra.
        Note: 0-frequency is preserved
        TO-DO: algorithm is slow and memory consuming, must be optimized.
        """
        slen = len(self) - 1
        freq = np.log(self.frequency[1:]) / (np.sqrt(2) * sigma)
        data = np.log(self.data[1:])
        s0 = self.data[0]

        if not memsafe:
            # Fast vectorial version, although memory consuming

            # Gaussian weighting window
            g = np.exp(-(np.tile(freq, slen) - np.repeat(freq, slen))**2)
            g = g.reshape(slen, slen)
            
            sumg = np.matmul(g, np.ones(slen))
            smat = g / sumg

            return np.insert(np.exp(np.matmul(data, smat)), 0, s0)

        else:
            # Slower version, but memory saving for large time series.
            sdata = np.zeros(slen, dtype='complex')

            for i,f0 in enumerate(freq):
                # Gaussian weighting window
                g = np.exp(-(freq - f0)**2)

                sdata[i]= np.matmul(data, g/sum(g))

            return np.insert(np.exp(sdata), 0, s0)
