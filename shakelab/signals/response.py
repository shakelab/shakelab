# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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

from shakelab.libutils.utils import (cast_value, 
                                     serialize_ndarray,
                                     deserialize_complex)
from shakelab.signals import fourier
from shakelab.signals import base
from shakelab.signals import stationxml
from shakelab.libutils.time import Date

ACTIVE_STAGES = {
    'gain' : True,
    'paz' : True,
    'polynomial' : False,
    'fir' : False
    }


class ResponseCollection:
    """
    A container of StreamResponse objects indexed by FDSN codes.
    """
    def __init__(self, byte_stream=None, ftype='sxml'):
        self.response = []
        if byte_stream is not None:
            self.read(byte_stream, ftype=ftype)

    def __len__(self):
        return len(self.response)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.response[self._idx(key)]
        return self.response[key]

    def _idx(self, sid):
        sid_list = self.sid
        if sid in sid_list:
            return sid_list.index(sid)
        print('Stream ID not found.')
        return None

    @property
    def sid(self):
        return [resp.sid for resp in self.response]

    def append(self, obj):
        if isinstance(obj, StreamResponse):
            if obj.sid not in self.sid:
                self.response.append(obj)
        elif isinstance(obj, StageSet):
            sid = obj.sid
            for resp in self.response:
                if resp.sid == sid:
                    resp.append(obj)
                    return
            sr = StreamResponse(sid)
            sr.append(obj)
            self.response.append(sr)
        else:
            raise ValueError('Unsupported type to append')

    def extract(self, sid, time=None, stage_number=None):
        """
        Return the stage or StageSet matching the stream ID and time.
        """
        resp = self[sid]
        result = resp.extract(time)
        if stage_number is not None:
            return result[stage_number]
        return result

    def read(self, byte_stream, ftype=None, append=False):
        """
        Read from an XML or JSON file. If ftype is None, infer it from
        file extension (.xml or .json). If append is False, replace content.
        """
        if ftype is None:
            _, ext = _os.path.splitext(str(byte_stream))
            ftype = ext.lstrip('.').lower()
    
        if ftype == 'xml':
            new = stationxml.parse_sxml(byte_stream)
            if append:
                for r in new.response:
                    self.append(r)
            else:
                self.response = new.response
    
        elif ftype == 'json':
            with open(byte_stream, 'r') as f:
                obj = _json.load(f)
    
            for entry in obj:
                sid = entry['sid']
                sset = StageSet(entry['starttime'], entry['endtime'])
    
                for s in entry['stages']:
                    stage_type = s['stage_type']
                    data = deserialize_complex(s['data'])
                    data['stage_number'] = s['stage_number']
    
                    if stage_type == 'gain':
                        stage = StageGain(data)
                    elif stage_type == 'paz':
                        stage = StagePoleZero(data)
                    elif stage_type == 'fir':
                        stage = StageFIR(data)
                    elif stage_type == 'coefficients':
                        stage = StageCoefficients(data)
                    elif stage_type == 'polynomial':
                        stage = StagePolynomial(data)
                    else:
                        errstr = f"Unsupported stage_type: {stage_type}"
                        raise ValueError(errstr)
    
                    sset.append(stage)
    
                if append and sid in self.sid:
                    self[sid].append(sset)
                else:
                    sr = StreamResponse(sid)
                    sr.append(sset)
                    self.append(sr)
    
        else:
            raise ValueError(f"Unknown file type: {ftype}")

    def write(self, file, ftype='json'):
        """
        Write the response data to a file in json format.
        """
        if ftype != 'json':
            raise NotImplementedError('Only sjson format is supported...')
    
        out = []
    
        for resp in self.response:
            for sset in resp.stageset:
                stages = []
                for stage in sset.stage:
                    stages.append({
                        'stage_number': stage.stage_number,
                        'stage_type': stage.stage_type,
                        'data': serialize_ndarray(stage.to_dict())
                    })
                out.append({
                    'sid': resp.sid,
                    'starttime': str(sset.starttime),
                    'endtime': str(sset.endtime) if sset.endtime else None,
                    'stages': stages
                })
    
        with open(file, 'w') as f:
            _json.dump(out, f, indent=2)

    def info(self):
        """
        Print info for all StreamResponses.
        """
        for resp in self.response:
            resp.info()


class StreamResponse():
    """
    Represents all stage sets for a given stream (by FDSN code).
    """
    def __init__(self, sid):
        self.sid = sid
        self.stageset = []

    def __len__(self):
        return len(self.stageset)

    def __getitem__(self, time):
        return self.extract(time)

    def append(self, stage_set):
        if not isinstance(stage_set, list):
            stage_set = [stage_set]
        for sset in stage_set:
            if not isinstance(sset, StageSet):
                raise ValueError('Not a valid stage set')
            self.stageset.append(sset)

    def extract(self, time=None):
        """
        Return the StageSet active at the specified time.
        """
        for sset in self.stageset:
            if sset.match(time):
                return sset
        return None

    def info(self):
        """
        Print info for each StageSet in this stream.
        """
        print(f"Stream ID: {self.sid}")
        for sset in self.stageset:
            sset.info()


class StageSet():
    """
    A time-bound set of stage responses valid between starttime and endtime.
    """
    def __init__(self, starttime=None, endtime=None):
        if starttime is not None and not isinstance(starttime, Date):
            starttime = Date(starttime)
        if endtime is not None and not isinstance(endtime, Date):
            endtime = Date(endtime)

        self.starttime = starttime
        self.endtime = endtime
        self.stage = []

    def __getitem__(self, idx):
        """
        Return the stage with the given stage_number (int).
        """
        if isinstance(idx, int):
            return self.get_stage(idx)
        raise NotImplementedError('Not a valid stage index')

    def info(self):
        """
        Print the validity interval and stage list of the StageSet.
        """
        start_str = str(self.starttime) if self.starttime else 'unknown'
        end_str = str(self.endtime) if self.endtime else 'now'
        print(f"Response valid from: {start_str} to: {end_str}")
    
        if not self.stage:
            print("Stages: none")
            return
    
        stage_str = "; ".join(
            f"({s.stage_number}) {s.stage_type}"
            for s in self.stage
        )
        print(f"Stages: {stage_str}")

    def append(self, stage_response):
        """
        Append one or more stage responses.
        """
        if not isinstance(stage_response, list):
            stage_response = [stage_response]

        for stage in stage_response:
            if not isinstance(stage, StageResponse):
                raise ValueError('Not a valid stage response')
            self.stage.append(stage)

    def match(self, time):
        """
        Check if a given time falls within the stage set's validity period.
        """
        if not isinstance(time, Date):
            time = Date('now') if time in ['now', None] else Date(time)

        return (
            (self.starttime is None or time >= self.starttime) and
            (self.endtime is None or time < self.endtime)
        )

    def get_stage(self, stage_number):
        """
        Return the stage with the specified stage_number.
        """
        for stage in self.stage:
            if stage.stage_number == stage_number:
                return stage
        raise IndexError(f'Stage {stage_number} not found')

    def convolve_record(self, record):
        """
        Apply all stages to the input record (single Record only).
        """
        if isinstance(record, base.Record):
            rec_mod = record.copy()
            for stage in self.stage:
                if ACTIVE_STAGES(stage.stage_type):
                    rec_mod = stage.convolve_record(rec_mod)
            return rec_mod
        raise NotImplementedError('Unsupported record type')

    def deconvolve_record(self, record):
        """
        Remove all stages from the input record (single Record only).
        """
        if isinstance(record, base.Record):
            rec_mod = record.copy()
            for stage in self.stage:
                if ACTIVE_STAGES(stage.stage_type):
                    rec_mod = stage.deconvolve_record(rec_mod)
            return rec_mod
        raise NotImplementedError('Unsupported record type')

    def copy(self):
        """
        Return a deep copy of the StageSet.
        """
        return deepcopy(self)


class StageResponse:
    """
    Base class for different stage responses.

    Attributes are defined in _KEYMAP with a structure:
    _KEYMAP = {
        "attribute_name": (expected_type, default_value)
    }

    Parameters can be set via a dictionary at initialization,
    additional keyword arguments, or loaded from a JSON file.
    """
    _KEYMAP = {}
    stage_type = None

    def __init__(self, data=None, **kwargs):
        # Initialise attributes to default values
        for key in self._KEYMAP:
            default_value = self._KEYMAP[key][1]
            setattr(self, key, default_value)

        # Apply data dictionary if provided
        if data is not None:
            self.from_dict(data)

        # Apply keyword arguments
        for key in kwargs:
            self[key] = kwargs[key]

    def __setitem__(self, key, value):
        """
        Set a parameter by key using type conversion defined in _KEYMAP.
        """
        if key in self._KEYMAP:
            dtype = self._KEYMAP[key][0]
            value = cast_value(value, dtype)
            setattr(self, key, value)
        else:
            raise KeyError(f'Invalid key: {key}')

    def __getitem__(self, key):
        """
        Get a parameter value by key.
        """
        if key in self._KEYMAP:
            return getattr(self, key)
        else:
            raise KeyError(f'Invalid key: {key}')

    def from_dict(self, data):
        """
        Set multiple parameters using a dictionary.
        """
        if isinstance(data, dict):
            for key, value in data.items():
                self[key] = value
        else:
            raise ValueError("Input data must be a dictionary.")

    def to_dict(self):
        """
        Get a dictionary of all parameters currently set.
        """
        return {key: getattr(self, key) for key in self._KEYMAP}

    def from_json(self, file_name):
        """
        Load parameters from a JSON file and update the object attributes.

        Parameters:
        - file_name (str): path to the JSON file

        Raises:
        - ValueError: if the file cannot be read or is not valid JSON
        """
        try:
            with open(file_name, 'r') as f:
                data = _json.load(f)
            self.from_dict(data)
        except (IOError, _json.JSONDecodeError) as e:
            raise ValueError(f"Failed to load JSON file: {e}")

    def to_json(self, file_name):
        """
        Save current parameters to a JSON file.

        Parameters:
        - file_name (str): path to the output JSON file
        """
        with open(file_name, 'w') as f:
            _json.dump(self.to_dict(), f, indent=2)


class StageGain(StageResponse):
    """
    Stage representing a simple gain stage with sensitivity.
    """
    _KEYMAP = {
        'description': (str, None),
        'sensitivity': (float, 1.0),
        'frequency': (float, None),
        'stage_number': (int, None)
    }
    stage_type = 'gain'

    def convolve_record(self, record):
        """
        Apply the gain (sensitivity) to the record.
        """
        if self.sensitivity is None:
            raise ValueError("Sensitivity is not set.")
        rec_mod = record.copy()
        rec_mod.data *= self.sensitivity
        return rec_mod

    def deconvolve_record(self, record):
        """
        Remove the gain (sensitivity) from the record.
        """
        if self.sensitivity is None or self.sensitivity == 0:
            raise ValueError("Sensitivity must be set and non-zero.")
        rec_mod = record.copy()
        rec_mod.data = rec_mod.data.astype(float) / self.sensitivity
        return rec_mod


class StagePoleZero(StageResponse):
    """
    Stage representing a response in pole-zero format.
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'pz_type': (str, None),
        'normalization_factor': (float, None),
        'normalization_frequency': (float, None),
        'poles': (complex, None),
        'zeros': (complex, None),
        'stage_number': (int, None)
    }
    stage_type = "paz"

    def response_function(self, frequency):
        """
        Evaluate the complex frequency response.
        """
        if self.poles is None:
            raise ValueError("Poles must be defined.")
        
        if self.zeros is None:
            self.zeros = np.array([], dtype=complex)

        return paz_transfer_function(
            self.normalization_factor,
            self.poles,
            self.zeros,
            frequency
        )

    def to_spectrum(self, delta, nsamp):
        """
        Return the complex frequency response as a Spectrum object.
        """
        frequency = fourier.frequency_axis(delta, nsamp)
        rf = self.response_function(frequency)

        sp = fourier.Spectrum()
        sp.set_data(rf, delta, nsamp)
        return sp

    def convolve_record(self, record):
        """
        Apply the response function to a record
        (frequency domain multiplication).
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        spec.data = spec.data.astype(complex) * resp.data
        return spec.to_record()

    def deconvolve_record(self, record, waterlevel=100):
        """
        Remove the response from a record using spectral division.
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        iresp = inverse_spectrum(resp.data, waterlevel)
        spec.data = spec.data.astype(complex) * iresp
        return spec.to_record()


class StagePolynomial(StageResponse):
    """
    Stage representing a polynomial transfer function of the form:
    H(jω) = (sum n_k * (jω)^k) / (sum d_k * (jω)^k)
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'numerator': (float, [1.0]),
        'denominator': (float, [1.0]),
        'stage_number': (int, None)
    }
    stage_type = "polynomial"

    def to_spectrum(self, delta, nsamp):
        """
        Compute the complex transfer function spectrum H(jω)
        for the given sampling interval and number of samples.
        """
        freq = fourier.frequency_axis(delta, nsamp)
        h = polynomial_transfer_function(self.numerator,
                                         self.denominator,
                                         freq)

        sp = fourier.Spectrum()
        sp.set_data(h, delta, nsamp)
        return sp

    def convolve_record(self, record):
        """
        Apply the polynomial transfer function to a record.
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        spec.data = spec.data.astype(complex) * resp.data
        return spec.to_record()

    def deconvolve_record(self, record, waterlevel=100):
        """
        Remove the polynomial response from a record.
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        iresp = inverse_spectrum(resp.data, waterlevel)
        spec.data = spec.data.astype(complex) * iresp
        return spec.to_record()


class StageFIR(StageResponse):
    """
    Stage representing a Finite Impulse Response (FIR) filter.
    """
    _KEYMAP = {
        'description': (str, None),
        'input_units': (str, None),
        'output_units': (str, None),
        'simmetry': (str, None),
        'normalization': (str, None),
        'numerator': (float, [1.0]),
        'denominator': (float, [1.0]),
        'stage_number': (int, None)
    }
    stage_type = "fir"

    def to_spectrum(self, delta, nsamp):
        """
        Return the FIR filter frequency response as a Spectrum object.
        """
        nsamp_half = fourier.rfft_length(nsamp)

        h, w = fir_transfer_function(
            self.numerator,
            self.denominator,
            nsamp_half,
            delta
        )

        sp = fourier.Spectrum()
        sp.set_data(h, delta, nsamp)
        return sp

    def convolve_record(self, record):
        """
        Apply the FIR filter to the record in the frequency domain.
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        spec.data = spec.data.astype(complex) * resp.data
        return spec.to_record()

    def deconvolve_record(self, record, waterlevel=100):
        """
        Remove the FIR response from the record using spectral division.
        """
        spec = record.to_spectrum()
        resp = self.to_spectrum(record.head.delta, record.nsamp)
        iresp = inverse_spectrum(resp.data, waterlevel)
        spec.data = spec.data.astype(complex) * iresp
        return spec.to_record()

#class StageDecimation(StageResponse):


def paz_transfer_function(a0, poles, zeros, freq):
    """
    Compute the complex transfer function H(jω) of a PAZ system.

    Parameters
    ----------
    a0 : float
        Normalization factor.
    poles : array-like of complex
        Poles of the transfer function (in radians/sec).
    zeros : array-like of complex
        Zeros of the transfer function (in radians/sec).
    freq : array-like
        Frequencies at which to evaluate the response (in Hz).

    Returns
    -------
    H : np.ndarray of complex
        Transfer function H(jω) evaluated at each frequency.
    """
    omega = 2 * np.pi * np.array(freq)
    s = 1j * omega

    zo = np.ones_like(s, dtype=np.complex128)
    po = np.ones_like(s, dtype=np.complex128)

    if zeros is not None and len(zeros) > 0:
        for z in zeros:
            zo *= (s - z)

    if poles is not None and len(poles) > 0:
        for p in poles:
            po *= (s - p)

    po = np.where(np.abs(po) < 1e-12, 1e-12, po)

    return a0 * (zo / po)

def polynomial_transfer_function(ncoeff, dcoeff, freq):
    """
    Compute the frequency response of a polynomial transfer function,
    where numerator and denominator are polynomials in s = jω.

    H(jω) = (sum n_k * (jω)^k) / (sum d_k * (jω)^k)

    Parameters
    ----------
    ncoeff : list or array-like
        Coefficients of the numerator polynomial (ascending order).
    dcoeff : list or array-like
        Coefficients of the denominator polynomial (ascending order).
    freq : array-like
        Frequencies in Hz.

    Returns
    -------
    H : np.ndarray
        Complex transfer function evaluated at each frequency.
    """
    omega = 2 * np.pi * np.array(freq)

    # Convert scalars to list
    if isinstance(ncoeff, (int, float)):
        ncoeff = [ncoeff]
    if isinstance(dcoeff, (int, float)):
        dcoeff = [dcoeff]

    s = 1j * omega

    num = sum(nc * s**k for k, nc in enumerate(ncoeff))
    den = sum(dc * s**k for k, dc in enumerate(dcoeff))

    # Prevent division by zero
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)

    return num / den

def fir_transfer_function(ncoeff, dcoeff, nsamp, delta=1.0):
    """
    Compute the complex frequency response of an FIR filter.

    Parameters
    ----------
    ncoeff : list or array-like
        FIR numerator coefficients (impulse response).
    dcoeff : list or array-like
        FIR denominator coefficients (typically [1.0]).
    nsamp : int
        Number of frequency points to compute.
    delta : float, optional
        Sampling interval in seconds.

    Returns
    -------
    H : np.ndarray
        Complex frequency response.
    freqs : np.ndarray
        Frequencies in Hz corresponding to the response.
    """
    # Handle None and scalar input
    if ncoeff is None or isinstance(ncoeff, (int, float)):
        ncoeff = [ncoeff or 1.0]

    if dcoeff is None or isinstance(dcoeff, (int, float)):
        dcoeff = [dcoeff or 1.0]

    omega = np.linspace(0, np.pi, nsamp)  # angular frequency (radians)
    jn = np.arange(len(ncoeff)).reshape(-1, 1)  # for broadcasting
    jd = np.arange(len(dcoeff)).reshape(-1, 1)

    # Numerator and denominator using vectorized matrix multiplication
    num = np.sum(np.array(ncoeff)[:, None] * np.exp(-1j * omega * jn), axis=0)
    den = np.sum(np.array(dcoeff)[:, None] * np.exp(-1j * omega * jd), axis=0)

    # Prevent division by zero
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)

    freqs = omega / (2 * np.pi * delta)  # Hz
    return num / den, freqs

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

    if method == 'sharp':
        # Preallocation of the inverse spectrum
        inv_spec = np.array([0+1j*0] * len(spectrum))

        # Identification of the spectrum above waterlevel
        i0 = (abs_spec >= wlev_db)

        inv_spec[i0] = 1/spectrum[i0]

    elif method == 'smooth':
        inv_spec = spectrum.conj()/(spectrum*spectrum.conj() + wlev_db)

        # Removing values below waterlevel (to check!)
        i0 = (abs_spec <= wlev_db)
        inv_spec[i0] = 0.

    else:
        raise ValueError('Not a valid method')

    return inv_spec

