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

import os as _os
import json as _json
import numpy as np

from shakelab.signals.fourier import fft, ifft, fft_axis


paz_map = {
    'description': None,
    'input_units': None,
    'output_units': None,
    'corner_frequency': None,
    'damping': None,
    'sensitivity': None,
    'normalization_factor': None,
    'normalization_frequency': None,
    'poles': None,
    'zeros': None
    }

class SensorResponse():
    """
    """

    def __init__(self, paz=None, **kwargs):
        self.paz = paz_map
        self.freq = None
        self.resp = None

        if paz is not None:
            self.paz = paz

        for key in kwargs:
            self.paz[key] = kwargs[key]

    def load_paz(self, sensor_id, json_file=None):
        """
        """
        self.paz = load_paz_from_file(sensor_id, json_file)

    def compute_response(self, frequency):
        """
        """
        self.freq = frequency
        self.resp = sensor_transfer_function(frequency, self.paz)


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

def sensor_transfer_function(frequency, paz):
    """
    Note: poles and zeros must be in radians/seconds
    """
    omega = 2*np.pi*np.array(frequency)

    sg = paz['sensitivity']
    a0 = paz['normalization_factor']

    zo = 1.
    for zn in paz['zeros']:
        zo *= (1j*omega - zn)

    po = 1.
    for pn in paz['poles']:
        po *= (1j*omega - pn)

    return sg * a0 *(zo/po)

def remove_sensor_response(record, polezeros, wlev=0.1):
    """
    """
    spec = fft(record.data)
    freq = fft_axis(len(record), record.dt)
    resp = sensor_response(freq, polezeros)

    # Compute inverse, including waterlevel
    i0 = (resp != 0.)
    resp[i0] = np.conj(resp[i0])/(resp[i0]*np.conj(resp[i0]) + wlev)

    return np.real(ifft(spec * resp))




