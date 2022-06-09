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
Module for basic waveform analysis
"""

import scipy as sp
import numpy as np

from scipy import signal, fftpack, integrate
from copy import deepcopy

from shakelab.libutils.time import Date
from shakelab.libutils.geodetic import WgsPoint
import shakelab.signals.fourier as fourier


class Header(object):
    """
    """
    def __init__(self):
        self._rate = None
        self._delta = None
        self.time = Date()
        self.location = WgsPoint(None, None)
        self.sid = StreamId()
        self.response = None
        self.meta = {}

    @property
    def delta(self):
        return self._delta

    @property
    def rate(self):
        return self._rate

    @delta.setter
    def delta(self, value):
        self._delta = value
        self._rate = 1./value

    @rate.setter
    def rate(self, value):
        self._rate = value
        self._delta = 1./value


class StreamId(object):
    """
    """
    def __init__(self, code=None):
        self.network = ''
        self.station = ''
        self.location = ''
        self.channel = ''

        if code is not None:
            self.set(code)

    def __repr__(self):
        return '{0}.{1}.{2}.{3}'.format(self.network,
                                        self.station,
                                        self.location,
                                        self.channel)

    def set(self, code):
        """
        """
        if isinstance(code, (list, tuple)):
            self.network = code[0]
            self.station = code[1]
            self.location = code[2]
            self.channel = code[3]

        if isinstance(code, str):
            code = code.split('.')
            self.network = code[0]
            self.station = code[1]
            self.location = code[2]
            self.channel = code[3]

        if isinstance(code, dict):
            self.network = code['network']
            self.station = code['station']
            self.location = code['location']
            self.channel = code['channel']


class Record(object):
    """
    Individual (continuos) recording block.
    """
    def __init__(self):
        self.head = Header()
        self.data = np.array([])

    def __len__(self):
        return len(self.data)

    def __getitem__(self, sliced):
        return self.data[sliced]

    def duration(self):
        """
        """
        return (len(self) - 1) * self.head.delta

    def time_axis(self, reference='relative', shift=0.):
        """
        to do: add reference
        """
        tax = np.arange(0., len(self)) * self.head.delta
        if reference in ['a', 'absolute']:
            tax += self.head.time.to_seconds()
        return tax + shift

    @property
    def time(self):
        return self.time_axis()

    def remove_mean(self):
        """
        """
        self.data -= np.mean(self.data)

    def filter(self, highpass=None, lowpass=None, order=2, minphase=False):
        """
        zero-phase and min-phase are allowed
        """
        # Corner frequencies
        corners = []

        if (highpass is not None):
            corners.append(2. * highpass * self.head.delta)
            filter_type = 'high'

        if (lowpass is not None):
            corners.append(2. * lowpass * self.head.delta)
            filter_type = 'low'

        if (highpass is not None) and (lowpass is not None):
            filter_type = 'band'

        if len(corners) > 0:
            # Butterworth filter
            sos = signal.butter(order, corners, analog=False,
                                btype=filter_type, output='sos')

        if minphase:
	        self.data = signal.sosfilt(sos, self.data)
        else:
            self.data = signal.sosfiltfilt(sos, self.data)

    def cut(self, starttime=None, endtime=None):
        """
        Cut the signal to the nearest time sample.
        """
        i0 = 0
        t0 = 0.
        i1 = len(self)
        t1 = self.duration()

        if (starttime is not None):
            if isinstance(starttime, Date):
                t0 = starttime - self.time
            elif isinstance(starttime, (int, float)):
                t0 = starttime

        if (endtime is not None):
            if isinstance(endtime, Date):
                t1 = endtime - self.time
            elif isinstance(endtime, (int, float)):
                t1 = endtime

        if (0. < t0 < self.duration()):
            i0 = int(np.argwhere(self.taxis() > t0)[0])

        if (0. < t1 < self.duration()):
            i1 = int(np.argwhere(self.taxis() > t1)[0])

        if (i1 > i0):
            self.data = self.data[i0:i1]
            self.time += t0
        else:
            print('Error: endtime before starttime')

    def taper(self, time=0.1):
        """
        time is in seconds.
        negative time means the whole window (cosine taper)
        """
        tnum = len(self)
        if time < 0:
            alpha = 1
        else:
            alpha = min(2 * float(time)/(self.head.delta * tnum), 1)
        self.data = (self.data * sp.signal.tukey(tnum, alpha))

    def zero_padding(self, time):
        """
        """
        zeros = np.zeros(round(time/self.head.delta))
        self.data = np.concatenate((self.data, zeros))

    def shift(self, time, padding=True):
        """
        Shift a signal in time by using fft-based circular convolution.
        """
        if padding:
            zeros = np.zeros(len(self))
            data = np.concatenate((self.data, zeros))
        else:
            data = self.data

        data = fourier.shift_time(data, self.head.delta, time)
        self.data = data[0:len(self.data)]

    def fft(self):
        """
        """
        return fourier.Spectrum(self)

    def ifft(self, spectrum):
        """
        """
        record = spectrum.ifft()
        self.head = record.head
        self.data = record.data

    def integrate(self, method='fft'):
        """
        """
        if method == 'cum':
            self.data = integrate.cumtrapz(self.data, dx=self.head.delta,
                                           initial=0)
        elif method == 'fft':
            self.data = fftpack.diff(self.data, order=-1,
                                     period=self.duration())
        else:
            raise NotImplementedError('method not implemented')

    def differentiate(self, method='fft'):
        """
        """
        if method == 'grad':
            self.data = np.gradient(self.data, self.head.delta)

        elif method == 'fft':
            self.data = fftpack.diff(self.data, order=1,
                                     period=self.duration())
        else:
            raise NotImplementedError('method not implemented')

    def convolve(self, record, mode='full', method='fft'):
        """
        """
        self.data = signal.convolve(self.data, record.data,
                                    mode=mode, method=method)

    def correlate(self, record, mode='full', method='fft'):
        """
        """
        self.data = signal.correlate(self.data, record.data,
                                     mode=mode, method=method)

    def copy(self):
        """
        """
        return deepcopy(self)

    def peak_amplitude(self):
        """
        """
        return np.max(np.abs(self.data))

    def significant_duration(self):
        """
        """
        pass

    def bracketed_duration(self):
        """
        """
        pass

    def integral_amplitude(self):
        """
        """
        pass


class Channel(object):
    """
    """

    def __init__(self):
        self.record = []

    def __len__(self):
        return len(self.record)

    def __getitem__(self, sliced):
        return self.record[sliced]


class Station(object):
    """
    Single measuring location in space with one or more recordings.
    It can be multichannel, but channels must be synchronous and
    related to the same measurement location.
    """

    keymap = {0 : 0, 'e' : 0, 'E' : 0, 'ew' : 0, 'EW' : 0, 'x' : 0, 'X' : 0,
              1 : 1, 'n' : 1, 'N' : 1, 'ns' : 1, 'NS' : 1, 'y' : 1, 'Y' : 1,
              2 : 2, 'u' : 2, 'U' : 2, 'ud' : 2, 'UD' : 2, 'z' : 2, 'Z' : 2}

    def __init__(self):
        self.id = None
        self.latitude = None
        self.longitude = None
        self.channel = []

    def __len__(self):
        return len(self.channel)

    def __getitem__(self, item):
        return self.channel[keymap[item]]


class Array(object):
    """
    Array of measuring locations.
    """

    def __init__(self):
        self.id = None
        self.station = []

    def __len__(self):
        return len(self.station)

    def __getitem__(self, sliced):
        return self.station[sliced]
