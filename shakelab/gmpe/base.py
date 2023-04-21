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

import abc as _abc
import os as _os
import json as _json

class GMPE(metaclass=_abc.ABCMeta):
    """
    Base class for ground motion prediction equation models (GMPEs).
    """

    _COEFF_FILE = None
    _COEFF_SET = 'default'

    def __init__(self):
        self.import_coeff_from_json(self._COEFF_FILE,
                                    self._COEFF_SET)

    @property
    @_abc.abstractmethod
    def REFERENCE_VELOCITY(self):
        pass

    @property
    @_abc.abstractmethod
    def DISTANCE_METRIC(self):
        pass

    @property
    @_abc.abstractmethod
    def MAGNITUDE_TYPE(self):
        pass

    @_abc.abstractmethod
    def ground_motion(self, imt, magnitude, distance):
        """
        """
        pass

    def import_coeff_from_json(self, json_file, coeff_set='default'):
        """
        Loads the coefficient from a separate file in json format.
        File has to be stored in the same directory of the GMPE class.
        """
        full_path = _os.path.dirname(__file__)
        path_file = _os.path.join(full_path, 'data', json_file)

        self.coeff = {}
        with open(path_file) as jf:
            data = _json.load(jf)['coefficients'][coeff_set]
            self.keys = data['keys']
            for kp in data['type']:
                self.coeff[kp] = [float(c) for c in data['type'][kp]]

    def get_coefficients(self, imt):
        """
        Extract coefficients corresponding to a given intensity
        measure type.
        """
        if imt in self.coeff:
            return dict(zip(self.keys, self.coeff[imt]))
        else:
            raise KeyError('Not a valid intensity measure type')

    def list_imts(self):
        """
        List the available intensity measure types for the gmpe.
        """
        return [k for k in self.coeff.keys()]

