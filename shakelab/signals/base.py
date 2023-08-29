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

import numpy as np

from scipy import signal, fftpack, integrate
from copy import deepcopy

from shakelab.signals import fourier
from shakelab.signals.libio import mseed
from shakelab.libutils.time import Date
from shakelab.libutils.constants import PI, GRAVITY
from shakelab.libutils.geodetic import WgsPoint
from shakelab.structures.response import (sdof_response_spectrum,
                                          sdof_interdrift,
                                          newmark_integration)


class Header(object):
    """
    """
    def __init__(self, parent_record=None):

        self._rate = None
        self._delta = None
        self._nsamp = None

        self.sid = None
        self.eid = None
        self.time = Date()
        self.location = WgsPoint(None, None)

        self.units = None
        self.response = None
        self.meta = {}

        # Reference to the host record is passed to access
        # record's methods and data within the header
        # (e.g. the recording length)
        self._parent = parent_record

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

    @property
    def nsamp(self):
        if self._parent is not None:
            return self._parent.nsamp
        else:
            return None

    def copy(self):
        return deepcopy(self)

class Record(object):
    """
    Individual (continuos) recording block.
    """
    def __init__(self, time=None, delta=None, data=None):
        self.head = Header(self)
        self.data = np.array([])

        if data is not None:
            self.data = np.array(data)

        if delta is not None:
            self.head.delta = delta

        if time is not None:
            self.head.time = time

    def __len__(self):
        return self.nsamp

    def __getitem__(self, sliced):
        return self.data[sliced]

    def __sum__(self, value):
        if isinstance(value, (int, float)):
            self.data += value
        elif isinstance(value, Record):
            self.append(value)

    def __mul__(self, value):
        if isinstance(value, (int, float)):
            self.data *= value

    @property
    def nsamp(self):
        """
        """
        return len(self.data)

    @property
    def delta(self):
        """
        """
        return self.head.delta

    @property
    def duration(self, precision=9):
        """
        """
        d = (len(self) - 1) * self.head.delta
        return round(d, precision)

    @property
    def time(self):
        return self.head.time

    @property
    def starttime(self):
        return self.head.time

    @property
    def endtime(self):
        return self.head.time + self.duration

    def time_axis(self, reference='relative', shift=0.):
        """
        to do: add reference
        """
        tax = np.arange(0., len(self)) * self.head.delta
        if reference in ['a', 'absolute']:
            tax += self.head.time.to_seconds()
        return tax + shift

    def append(self, record, enforce=False, fillvalue=0., precision=9):
        """
        """
        if self.delta != record.delta:
            print('Samping rate must be uniform. Not merging.')
            return False

        if self.time > record.time:
            print('New record starts before previous record. Not merging.')
            return False

        d0 = round(self.duration + self.delta, precision)
        d1 = round(record.time - self.time, precision)

        if (d1 == d0):
            self.data = np.concatenate((self.data, record.data))
            return True

        else:
            (q, r) = divmod(round(d1 - d0, precision), self.delta)

            if (q > 0) and enforce:
                if r == 0:
                    infill = np.ones(int(q)) * fillvalue
                    segments = (self.data, infill, record.data)
                    self.data = np.concatenate(segments)
                    return True
                else:
                    print('Sampling mismatch')
                    return False

            elif (q < 0) and enforce:
                if r == 0:
                    segments = (self.data[0:int(q)], record.data)
                    self.data = np.concatenate(segments)
                    return True
                else:
                    print('Sampling mismatch')
                    return False

            else:
                print('Not contiguous data')
                return False

    def remove_mean(self):
        """
        """
        self.data = self.data - np.mean(self.data)

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

    def cut(self, starttime=None, endtime=None, inplace=True):
        """
        Cut the signal in place to the nearest time sample.
        """
        i0 = 0
        t0 = 0.
        i1 = len(self) - 1
        t1 = self.duration

        if (starttime is not None):
            if isinstance(starttime, Date):
                t0 = starttime - self.head.time
            elif isinstance(starttime, (int, float)):
                t0 = starttime

        if (endtime is not None):
            if isinstance(endtime, Date):
                t1 = endtime - self.head.time
            elif isinstance(endtime, (int, float)):
                t1 = endtime

        if (0. < t0 < self.duration):
            i0 = int(np.argwhere(self.time_axis() >= t0)[0])

        if (0. < t1 < self.duration):
            i1 = int(np.argwhere(self.time_axis() <= t1)[-1])

        if (i1 > i0):
            if inplace:
                self.data = self.data[i0:i1+1]
                self.head.time += t0
            else:
                rec = self.copy()
                rec.data = self.data[i0:i1+1]
                rec.head.time += t0
                return rec

        else:
            print('Warning: endtime before starttime')
            return None

    def extract(self, starttime=None, endtime=None):
        """
        Return a new record with the requested time window
        """
        return self.cut(starttime, endtime, inplace=False)

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
        self.data = (self.data * signal.tukey(tnum, alpha))

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
                                     period=self.duration)
        else:
            raise NotImplementedError('method not implemented')

    def differentiate(self, method='fft'):
        """
        """
        if method == 'grad':
            self.data = np.gradient(self.data, self.head.delta)

        elif method == 'fft':
            self.data = fftpack.diff(self.data, order=1,
                                     period=self.duration)
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

    @property
    def analytic_signal(self):
        """
        Compute the analytic signal using the Hilbert transform.
        The Hilbert transformed signal will be the imaginary part of
        the analytical signal.
        """
        return signal.hilbert(self.data)

    @property
    def amplitude_envelope(self):
        """
        """
        return np.abs(self.analytic_signal)

    @property
    def instantaneous_phase(self):
        """
        """
        return np.unwrap(np.angle(self.analytic_signal))

    @property
    def instantaneous_frequency(self):
        """
        """
        return np.diff(self.instantaneous_phase) / (2*PI) * self.head.rate

    def peak_amplitude(self):
        """
        """
        return np.max(np.abs(self.data))

    def arias_intesity(self):
        """
        """
        integral = integrate.trapz(self.data**2, dx=self.head.delta)
        return PI * integral / (2*GRAVITY)

    def cumulative_absolute_velocity(self):
        """
        """
        return integrate.trapz(np.abs(self.data), dx=self.head.delta)

    def bracketed_duration(self, threshold=0.05):
        """
        """
        data = np.argwhere(np.abs(self.data) >= threshold)

        if data.size != 0:
            i0 = data[0][0]
            i1 = data[-1][0]

            time = self.time_axis()
            return time[i1] - time[i0]

        else:
            return None

    def significant_duration(self, threshold=(0.05,0.95)):
        """
        """
        num = integrate.cumtrapz(self.data**2, dx=self.head.delta)
        den = integrate.trapz(self.data**2, dx=self.head.delta)
        cum_arias = num/den

        i0 = np.argwhere(cum_arias >= threshold[0])[0][0]
        i1 = np.argwhere(cum_arias <= threshold[1])[-1][0]

        time = self.time_axis()
        return time[i1] - time[i0]

    def root_mean_square(self):
        """
        """
        integral = integrate.trapz(self.data**2, dx=self.head.delta)

        return np.sqrt(integral/self.duration)

    def sdof_response_spectrum(self, periods, zeta=0.05):
        """
        """
        periods = np.array(periods, ndmin=1)

        rssp = sdof_response_spectrum(self.data, self.head.delta,
                                      periods, zeta=zeta)

        return {'sd' : rssp[0],
                'sv' : rssp[1],
                'sa' : rssp[2],
                'psv' : rssp[3],
                'psa' : rssp[4]}

    def sdof_convolve(self, period, zeta=0.05):
        """
        """
        resp = newmark_integration(self.data, self.head.delta,
                                   period, zeta=zeta)

        return {'d' : resp[0],
                'v' : resp[1],
                'a' : resp[2]}

    def sdof_interdrift(self, period, zeta=0.05):
        """
        """
        return sdof_interdrift(self.data, self.head.delta,
                               period, zeta=zeta)

    def soil1d_convolve(self, model1d, component='sh', angle=0.):
        """
        """
        pass

    def copy(self):
        """
        """
        return deepcopy(self)


class Stream(object):
    """
    Representation of a single stream (or channel), which could be 
    continuos or with gaps.
    """
    def __init__(self, id):
        self.sid = id
        self.record = []

    def __len__(self):
        return len(self.record)

    def __getitem__(self, id):
        """
        """
        if isinstance(id, str):
            return self.record[self._idx(id)]
        else:
            return self.record[id]

    def _idx(self, id):
        """
        """
        eid_list = self.eid
        if id in eid_list:
            return eid_list.index(id)
        else:
            print('Id not found.')
            return None

    @property
    def eid(self):
        """
        """
        return [record.head.eid for record in self.record]

    def append(self, record, enforce=False):
        """
        """
        if record.head.sid != self.sid:
            assert ValueError('Record ID mismatching')

        if not self.record:
            self.record = [record]
        else:
            if not self.record[-1].append(record, enforce=enforce):
                self.record.append(record)

    def remove(self):
        pass

    def get(self, eid=None, starttime=None, endtime=None):
        """
        Return selected record, force merge if record
        not contiguous.
        """
        if eid is not None:
            return self[eid].extract(starttime, endtime)

        else:
            out = None
            for rec in self.record:
                sel = rec.extract(starttime, endtime)
                if sel is not None:
                    if out is None:
                        out = sel
                    else:
                        out.append(sel, enforce=True)
            return out
        

    def sort(self):
        """
        """
        def get_time(rec):
            return rec.head.time.seconds

        self.record.sort(key=get_time)

    def fix(self):
        """
        Utility to fix gaps and overlaps between records
        """
        pass


class StreamCollection():
    """
    """
    def __init__(self):
        self.stream = []

    def __len__(self):
        return len(self.stream)

    def __getitem__(self, id):
        """
        """
        if isinstance(id, str):
            return self.stream[self._idx(id)]
        else:
            return self.stream[id]

    def _idx(self, id):
        """
        """
        sid_list = self.sid
        if id in sid_list:
            return sid_list.index(id)
        else:
            print('Id not found.')
            return None

    @property
    def sid(self):
        """
        """
        return [stream.sid for stream in self.stream]

    def add(self, record):
        """
        """
        sid = record.head.sid
        if sid not in self.sid:
            self.stream.append(Stream(sid))

        self[sid].append(record)

    def remove(self):
        pass

    def get(self, id, starttime=None, endtime=None):
        """
        """
        return self[id].get(starttime, endtime)

    def merge(self, stream_collection):
        """
        """
        for stream in self.stream:
            for record in stream.record:
                self.add(record)

    def read(self, byte_stream, ftype='mseed', byte_order='be'):
        """
        """
        if ftype == 'mseed':
            mseed.msread(byte_stream, byte_order, self)

    def write(self, file):
        """
        """
        pass


# ---------------------
"""

class Location(object):
    def __init__(self, id):
        self.id = id
        self.channel = {}

    def __len__(self):
        return len(self.channel)

    def __getitem__(self, item):
        return self.channel[item]

    def add(self):
        pass

    def remove(self):
        pass

class Station(object):

    def __init__(self, id):
        self.id = id
        self.location = {}

    def __len__(self):
        return len(self.location)

    def __getitem__(self, item):
        return self.location[item]

    def add(self):
        pass

    def remove(self):
        pass

class StationCollection(object):

    def __init__(self):
        self.id = None
        self.station = {}  # Must be list! implement iterators

    def __len__(self):
        return len(self.station)

    def __getitem__(self, sliced):
        return self.station[sliced]

    def add(self):
        pass

    def remove(self):
        pass

    def filter(self, streamid=None, station=None, channel=None):
        pass
"""