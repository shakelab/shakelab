# ============================================================================
#
# Copyright (C) 2019 Valerio Poggi.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ============================================================================
"""
"""

import abc as _abc
import os as _os
import json as _json

class GMPE(metaclass=_abc.ABCMeta):
    """
    Base class for ground motion prediction equation models (GMPEs).
    """

    def __init__(self, json_file, version='default'):
        self.import_coeff_from_json(json_file, version)

    def import_coeff_from_json(self, json_file, version='default'):
        """
        Loads the coefficient from a separate file in json format.
        File has to be stored in the same directory of the GMPE class.
        """
        full_path = _os.path.dirname(__file__)
        path_file = _os.path.join(full_path, 'data', json_file)

        self.coeff = {}
        with open(path_file) as jf:
            data = _json.load(jf)['coefficients'][version]
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
            print('Error: Not a valid intensity measure type')

    @_abc.abstractmethod
    def ground_motion(self, imt, magnitude, distance):
        """
        """
        pass



