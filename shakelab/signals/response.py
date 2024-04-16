# ****************************************************************************
#
# Copyright (C) 2019-2023, ShakeLab Developers.
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

import os as _os
import json as _json
import numpy as np
from copy import deepcopy

from shakelab.libutils.utils import cast_value
from shakelab.signals import fourier
from shakelab.signals import base
from shakelab.signals import stationxml
from shakelab.libutils.time import Date


class ResponseCollection():
    """
    TO DO: It should be derived from a base class, which includes
           also StreamCollection.
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

    def append(self, response):
        """
        TO DO: must include also the possibility to add single stages.
        """
        sid = response.sid
        if sid not in self.sid:
            self.stream.append(response)

    def remove(self):
        pass

    def get(self, id, time=None, stage_number=None):
        """
        """
        out = self[id]
        if time is not None:
            out = out.get(time)
        if stage_number is not None:
            out = out.get(stage_number)
        return out

    def read(self, byte_stream, ftype='sxml'):
        """
        """
        if ftype == 'sjson':
            pass

        if ftype == 'sxml':
            # TO IMPROVE: should not replace but append
            buf = stationxml.parse_sxml(byte_stream)
            self.stream = buf.stream

    def write(self, file):
        """
        """
        pass


class StreamResponse():
    """
    """
    def __init__(self, sid):
        self.sid = sid
        self.record = []

    def __len__(self):
        return len(self.record)

    def __getitem__(self, time):
        """
        """
        return self.get(time)

    def append(self, stage_record):
        """
        """
        if not isinstance(stage_record, list):
            stage_record = [stage_record]

        for record in stage_record:
            if isinstance(record, StageRecord):
                self.record.append(record)
            else:
                raise ValueError('Not a valid stage record')

    def remove(self):
        pass

    def get(self, time=None):
        """
        """
        for record in self.record:
            if record.match(time):
                return record
        return None


class StageRecord():
    """
    """
    def __init__(self, starttime=None, endtime=None):

        if starttime is not None:
            if not isinstance(starttime, Date):
                starttime = Date(starttime)

        if endtime is not None:
            if not isinstance(endtime, Date):
                endtime = Date(endtime)

        self.starttime = starttime
        self.endtime = endtime
        self.stage = []

    def __getitem__(self, idx):
        """
        """
        if isinstance(idx, int):
            return self.select_stage(idx)
        else:
            raise NotImplementedError('Not a valid stage number')

    def append(self, stage_response):
        """
        """
        if not isinstance(stage_response, list):
            stage_response = [stage_response]

        for stage in stage_response:
            if isinstance(stage, StageResponse):
                self.stage.append(stage)
            else:
                raise ValueError('Not a valid stage response')

    def match(self, time):
        """
        """
        if not isinstance(time, Date):
            if time == 'now' or time is None:
                time = Date('now')
            else:
                time = Date(time)

        if time >= self.starttime:
            if self.endtime is None:
                return True
            else:
                if time < self.endtime:
                    return True
                else:
                    return False
        else:
            return False

    def select_stage(self, stage_number):
        """
        """
        srec = StageRecord()
        srec.starttime = self.starttime
        srec.endtime = self.endttime

        for stage in self.stage:
            if stage.stage_number == stage_number:
                srec.stage.append(stage)

    #def get(self, stage_number=None):
    #    """
    #    Select specific stage sequence numbers
    #    """
    #    if stage_number is not None:
    #        return [s for s in self.stage if s.stage_number==stage_number+1]

    def convolve_record(self, record):
        """
        """
        if isinstance(record, base.Record):
            rec_mod = record.copy()
            for s in s in self.stage:
                rec_mod = s.convolve_record(rec_mod)
            return rec_mode

        elif isinstance(record, base.Stream):
            pass

        elif isinstance(record, base.StreamCollection):
            pass

    def deconvolve_record(self, record):
        """
        """
        if isinstance(record, base.Record):
            rec_mod = record.copy()
            for s in s in self.stage:
                rec_mod = s.deconvolve_record(rec_mod)
            return rec_mode

        elif isinstance(record, base.Stream):
            pass

        elif isinstance(record, base.StreamCollection):
            pass


    def copy(self):
        """
        """
        return deepcopy(self)

class StageResponse():
    """
    Base class for different stage responses
    """
    _KEYMAP = {}
    stage_type = None

    def __init__(self, data={}, **kwargs):

        # Initialise attributes to default value
        for key in self._KEYMAP:
            value = self._KEYMAP[key][1]
            exec('self.{0}={1}[0]'.format(key, [value]))

        if data is not None:
            self.set(data)

        for key in kwargs:
            self[key] = kwargs[key]

    def __setitem__(self, key, value):

        if key in self._KEYMAP:
            value = cast_value(value, self._KEYMAP[key][0])
            exec('from numpy import array')
            exec('self.{0}={1}[0]'.format(key, [value]))
        else:
            raise KeyError('{0}'.format(key))

    def __getitem__(self, key):

        if key in self._KEYMAP:
            return eval('self.{0}'.format(key))
        else:
            raise KeyError('{0}'.format(key))

    def set(self, data):
        """
        """
        if isinstance(data, dict):
            for (key, value) in data.items():
                self[key] = value
        else:
            raise ValueError('not a valid data format')

    def get(self):
        """
        """
        data = {}
        for key in self.keys:
            data[key] = self[key]
        return data


class StageGain(StageResponse):
    """
    """
    _KEYMAP = {
        'description': (str, None),
        'sensitivity': (float, None),
        'frequency' : (float, None),
        'stage_number' : (int, None)
        }
    stage_type = 'gain'

    def convolve_record(self, record):
        """
        """
        rec_mod = record.copy()
        rec_mod.data = rec_mod.data * self.sensitivity
        return rec_mod

    def deconvolve_record(self, record):
        """
        """
        rec_mod = record.copy()
        rec_mod.data = rec_mod.data / self.sensitivity
        return rec_mod


class StagePoleZero(StageResponse):
    """
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'pz_type' : (str, None),
        'normalization_factor': (float, None),
        'normalization_frequency': (float, None),
        'poles': (np.array, None),
        'zeros': (np.array, None),
        'stage_number' : (int, None)
        }
    stage_type = "paz"

    def load_from_file(self, sensor_id, json_file=None):
        """
        """
        self.paz = load_paz_from_file(sensor_id, json_file)

    def response_function(self, frequency):
        """
        """
        return paz_transfer_function(self.normalization_factor,
                                     self.poles, self.zeros,
                                     frequency)

    def to_spectrum(self, delta, nsamp):
        """
        """
        frequency = fourier.frequency_axis(delta, nsamp)
        rf = self.response_function(frequency)

        sp = fourier.Spectrum()
        sp.set_data(rf, delta, nsamp)
        return sp

    def convolve_record(self, record):
        """
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        spec.data *= resp

        return spec.to_record()

    def deconvolve_record(self, record, waterlevel=100):
        """
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        iresp = inverse_spectrum(resp.data, waterlevel)
        spec.data *= iresp

        return spec.to_record()


class StagePolynomial(StageResponse):
    """
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'numerator': (np.array, None),
        'denominator': (np.array, None),
        'stage_number' : (int, None)
        }

    def convolve_record(self, record):
        """
        """
        pass

    def deconvolve_record(self, record):
        """
        """
        pass


class StageFIR(StageResponse):
    """
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'simmetry' : (str, None),
        'normalization' : (str, None),
        'numerator': (np.array, 1),
        'denominator': (np.array, 1),
        'stage_number' : (int, None)
        }

    def to_spectrum(self, delta, nsamp):
        """
        """
        nsamp_half = fourier.rfft_length(nsamp)

        h, w = fir_transfer_function(self.numerator,
                                     self.denominator,
                                     nsamp_half, delta)

        sp = fourier.Spectrum()
        sp.set_data(h, delta, nsamp)
        return sp

    def convolve_record(self, record):
        """
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        spec.data *= resp

        return spec.to_record()

    def deconvolve_record(self, record, waterlevel=100):
        """
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        iresp = inverse_spectrum(resp.data, waterlevel)
        spec.data *= iresp

        return spec.to_record()


#class StageDecimation(StageResponse):


def load_paz_from_file(sensor_id, paz_file=None):
    """
    """
    if paz_file is None:
        full_path = _os.path.dirname(__file__)
        json_file = _os.path.join(full_path, 'data', 'sensor_paz.json')

    with open(paz_file) as jf:
        paz = _json.load(jf)

    # Converting to complex number format
    for k in paz.keys():
        paz[k]['poles'] = [p[0]+p[1]*1j for p in paz[k]['poles']]
        paz[k]['zeros'] = [p[0]+p[1]*1j for p in paz[k]['zeros']]

    return paz[sensor_id]

def paz_transfer_function(a0, poles, zeros, freq):
    """
    Note: poles and zeros must be in radians/seconds

    a0 = normalisation factor
    """
    omega = 2*np.pi*freq

    zo = np.ones_like(omega, dtype=np.complex128)
    po = np.ones_like(omega, dtype=np.complex128)

    for zn in zeros:
        zo *= (1j*omega - zn)

    for pn in poles:
        po *= (1j*omega - pn)

    return a0 * (zo/po)

def polynomial_transfer_function(ncoeff, dcoeff, freq):
    """
    TO CHECK!
    """
    omega = 2*np.pi*freq

    if isinstance(ncoeff, (int, float)):
        ncoeff = [ncoeff]

    if isinstance(dcoeff, (int, float)):
        dcoeff = [dcoeff]

    nlen = len(ncoeff)
    dlen = len(dcoeff)

    num = np.zeros_like(omega, dtype=np.complex128)
    den = np.zeros_like(omega, dtype=np.complex128)

    for n, nc in enumerate(ncoeff):
        num += (nc*(1j*omega)**(nlen-n))

    for n, dc in enumerate(dcoeff):
        den += (dc*(1j*omega)**(dlen-n))

    return num/den

def fir_transfer_function(ncoeff, dcoeff, nsamp, delta=1):
    """
    """
    if isinstance(ncoeff, (int, float)):
        ncoeff = [ncoeff]

    if isinstance(dcoeff, (int, float)):
        dcoeff = [dcoeff]

    omega = np.linspace(0, np.pi, nsamp)
    
    num = np.zeros_like(omega, dtype=np.complex128)
    den = np.zeros_like(omega, dtype=np.complex128)
    
    for n, nc in enumerate(ncoeff):
        num += nc * np.exp(-1j * omega * n)

    for n, dc in enumerate(dcoeff):
        den += dc * np.exp(-1j * omega * n)

    return num / den, omega / (2*np.pi*delta)

def inverse_spectrum(spectrum, waterlevel=100, method='smooth'):
    """
    THIS COULD BE MOVED TO FOURIER MODULE
    """
    if not isinstance(spectrum, np.ndarray):
        spectrum = np.array(spectrum)

    abs_spec = np.abs(spectrum)

    # Conversion from Db to actual spectrum level (to check!)
    #wlev_db = abs_spec.max() * 10.0 ** (-waterlevel / 20.0)
    wlev_db = 10.0 ** (-waterlevel / 20.0)

    match method:
        case 'sharp':
            # Preallocation of the inverse spectrum
            inv_spec = np.array([0+1j*0] * len(spectrum))

            # Identification of the spectrum above waterlevel
            i0 = (abs_spec >= wlev_db)

            inv_spec[i0] = 1/spectrum[i0]

        case 'smooth':
            inv_spec = spectrum.conj()/(spectrum*spectrum.conj() + wlev_db)

            # Removing values below waterlevel (to check!)
            i0 = (abs_spec <= wlev_db)
            inv_spec[i0] = 0.

        case _:
            raise ValueError('Not a valid method')

    return inv_spec