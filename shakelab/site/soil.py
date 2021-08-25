# ****************************************************************************
#
# Copyright (C) 2019-2021, ShakeLab Developers.
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

from shakelab.site import engpar


_KEYMAP = ['hl', 'vp', 'vs', 'dn', 'qp', 'qs']


class Layer():

    def __init__(self, data=None):

        # Initialise attributes to default value
        for key in _KEYMAP:
            exec('self.{0}={1}[0]'.format(key, [None]))

        if data is not None:
            self.set(data)

    def __setitem__(self, key, value):

        if key in _KEYMAP:
            exec('self.{0}=float({1})'.format(key, value))
        else:
            raise KeyError('{0}'.format(key))

    def __getitem__(self, key):

        if key in _KEYMAP:
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

    @property
    def keys(self):
        return _KEYMAP

class Model1D():

    def __init__(self, ascii_file=None):
        self.layer = []

        if ascii_file is not None:
            self.read(ascii_file)

    def __getitem__(self, key):

        if key in _KEYMAP:
            return eval('[l.{0} for l in self.layer]'.format(key))
        else:
            raise KeyError('{0}'.format(key))

    def __len__(self):
        return len(self.layer)

    def add_layer(self, data=None, idx=-1):
        """
        Method to add a single (horizontal) layer to the 1d
        soil profile at arbitrary location.

        :param list or dictionary data:
            data can be a list of values sorted according to key list,
            or a dictionary with the corresponding keys

        :param int index:
            index of the position along the profile where the layer
            should be added. Use -1 for the last layer (default).
        """
        if idx < 0:
            idx = len(self.layer) + idx + 1

        self.layer.insert(idx, Layer(data))

    def del_layer(self, idx=-1):
        """
        Method to remove a single layer from te 1d soil profile
        ar arbitrary location

        :param int index:
            index of the position along the profile where the layer
            should be removed. Use -1 for the last layer (default).
        """
        del self.layer[idx]

    def set(self, key, data):
        print('set')

    def get(self, key):
        return self[key]

    def __set_hl(self, data):
        self.set('hl', data)

    def __get_hl(self):
        return self.get('hl')

    def __set_vp(self, data):
        self.set('vp', data)

    def __get_vp(self):
        return self.get('vp')

    def __set_vs(self, data):
        self.set('vs', data)

    def __get_vs(self):
        return self.get('vs')

    def __set_dn(self, data):
        self.set('dn', data)

    def __get_dn(self):
        return self.get('dn')

    def __set_qp(self, data):
        self.set('qp', data)

    def __get_qp(self):
        return self.get('qp')

    def __set_qs(self, data):
        self.set('qs', data)

    def __get_qs(self):
        return self.get('qs')

    hl = property(__get_hl, __set_hl)
    vp = property(__get_vp, __set_vp)
    vs = property(__get_vs, __set_vs)
    dn = property(__get_dn, __set_dn)
    qp = property(__get_qp, __set_qp)
    qs = property(__get_qs, __set_qs)

    def read(self, ascii_file, header=None, skipline=0,
                   comment='#', delimiter=','):
        """
        Method to parse soil properties from tabular ascii file;
        arbitrary formatting is allowed

        :param string ascii_file:
            input model file. Default format is:
                hl,vp,vs,dn
                10,300,200,1900
                10,500,300,1900

        :param list header:
            list of header keys, to be used when not
            available within the input file

        :param int skipline:
            number of intitial lines to skip;
            default value is 0

        :param char or string comment:
            string to mark comments (which are not parsed);
            default value is the hash character

        :param char delimiter:
            character separator between data fields;
            default value is comma
        """

        # Open input ascii file
        with open(ascii_file, 'r') as f:

            # Ignore initial line(s) if necessary
            for _ in range(0, skipline):
                next(f)

            for line in f:
                line = line.strip()

                # Skip comments
                if line[0] != comment:
                    data = line.split(delimiter)

                    # Import header and data
                    if not header:
                        header = data
                    else:
                        layer = {k: float(d) for k, d in zip(header, data)}
                        self.add_layer(layer)

            f.close()
            return

        print('Error: Wrong file or file path')

    def traveltime_velocity(self, key='vs', depth=30.):
        """
        Compute and store travel-time average velocity at a given depth.
        Multiple depths are also allowed (as list).
        Default depth is Vs30.

        :param string key:
            velocity to compute ('vp', 'vs')

        :param float depth:
            calculation depth
        """

        if key not in ['vp', 'vs']:
            raise KeyError('{0}'.format(key))

        return engpar.traveltime_velocity(np.array(self['hl']),
                                          np.array(self[key]),
                                          depth)

    @property
    def vs30(self):
        """
        """
        return self.traveltime_velocity('vs', 30)
