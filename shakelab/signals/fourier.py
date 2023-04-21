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
from scipy import interpolate
from copy import deepcopy

import shakelab.signals.base as base

def _fnum(tnum):
    """
    Computing the length of RFFT (positive) frequency axis,
    accounting for odd/even number of samples
    """
    if (tnum % 2) == 0:
        return int((tnum/2)+1)
    else:
        return int((tnum+1)/2)

def _tnum(delta, dfreq):
    """
    """
    return int((1 / dfreq) / delta)

def delta_to_dfreq(delta, tnum):
    """
    """
    return 1 / (delta * tnum)

def dfreq_to_delta(dfreq, tnum):
    """
    """
    return 1 / (dfreq * tnum)

def fft(data, nsamp):
    """
    wrapper
    """
    return np.fft.rfft(data, nsamp)

def ifft(data, nsamp):
    """
    wrapper
    """
    return np.fft.irfft(data, nsamp)

def frequency_axis(delta, tnum):
    """
    """
    return np.fft.rfftfreq(tnum, delta)

def frequency_range(fmin, fmax, fnum, log=True):
    """
    Compute a linear or logarithmic frequency range

    :param float fmin:
        Minimum frequency

    :param float fmax:
        Maximum frequency

    :param int fnum:
        Number of frequencies

    :param boolean log:
        Select between linear or logarithmic spacing
        (default is logarithmic)

    :return numpy.array freq:
        The frequency axis
    """

    if log:
        freq = np.logspace(np.log10(fmin), np.log10(fmax), fnum)
    else:
        freq = np.linspace(fmin, fmax, fnum)

    return freq

def _to_complex(amplitude, phase):
    """
    """
    return amplitude * np.exp(1j * phase)

def shift_time(signal, delta, shift):
    """
    Shift a signal in time by using fft-based circular convolution.
    No zero-padding is assumed.
    """
    nsamp = len(signal)
    freq = frequency_axis(delta, len(signal))
    expt = np.exp(-2*1j*np.pi*shift*freq)

    return ifft(fft(signal, nsamp)*expt, nsamp)


class Spectrum():
    """
    Discrete Fourier spectrum class in RFFT format (one sided positive)
    """

    def __init__(self, record=None, delta=None, nsamp=None, data=None):

        self.head = base.Header(self)
        self.data = np.array([], dtype="complex")

        if record is not None:
            self.fft(record)
        else:
            if delta is not None:
                self.head.delta = delta

            if nsamp is not None:
                self.head.nsamp = nsamp

            if data is not None:
                self.data = data        

    def __len__(self):
        return self.nfreq

    def __getitem__(self, sliced):
        return self.data[sliced]

    @property
    def nfreq(self):
        return _fnum(self.head.nsamp)

    @property
    def df(self):
        """
        """
        return delta_to_dfreq(self.head.delta, self.head.nsamp)

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
    def phase_unwrap(self):
        """
        """
        return np.unwrap(self.phase)

    @property
    def frequency(self):
        """
        """
        return frequency_axis(self.head.delta, self.head.nsamp)

    def fft(self, record, norm=False):
        """
        """
        self.head = deepcopy(record.head)
        self.data = fft(record.data, self.head.nsamp)

        if norm:
            # to check
            self.data =  self.data / record.dt / len(self.head.nsamp)

    def ifft(self, norm=False):
        """
        """
        record = base.Record()
        record.head = deepcopy(self.head)
        record.data = ifft(self.data, self.head.nsamp)

        if norm:
            # to check
            record.data = record.data / self.df * len(self.head.nsamp)

        return record

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

    def logsmooth(self, sigma=0.2, memsafe=False):
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

