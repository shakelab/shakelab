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

from shakelab.site import engpar
from shakelab.site.cps.surf96 import surf96
from shakelab.site.response import soil_response
from shakelab.site.engpar import compute_site_kappa
from shakelab.signals.fourier import frequency_axis, Spectrum
from shakelab.libutils import utils

# Precision for decimal rounding
DECIMALS = 6

# Map of soil parameters
_KEYMAP = ['hl', 'vp', 'vs', 'dn', 'qp', 'qs']


class Layer():
    """
    """

    def __init__(self, data=None):

        # Initialise attributes to default value
        for key in _KEYMAP:
            exec('self.{0}={1}[0]'.format(key, [None]))

        if data is not None:
            self.set(data)

    def __repr__(self):
        return repr(self.get())

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
        Import data in dictionary or list format.
        """
        if isinstance(data, dict):
            for (key, value) in data.items():
                self[key] = value
        elif isinstance(data, list):
            for (key, value) in zip(_KEYMAP, data):
                self[key] = value
        else:
            raise ValueError('not a valid data format')

    def get(self):
        """
        Export data in dictionary format.
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

    def __setitem__(self, key, data):

        if key in _KEYMAP:
            for i, value in enumerate(data):
                exec('self.layer[{0}].{1}=float({2})'.format(i, key, value))
        else:
            raise KeyError('{0}'.format(key))

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

    def set_key(self, key, data):
        self[key] = data

    def get_key(self, key):
        return self[key]

    def _set_hl(self, data):
        self.set_key('hl', data)

    def _get_hl(self):
        return self.get_key('hl')

    def _set_vp(self, data):
        self.set_key('vp', data)

    def _get_vp(self):
        return self.get_key('vp')

    def _set_vs(self, data):
        self.set_key('vs', data)

    def _get_vs(self):
        return self.get_key('vs')

    def _set_dn(self, data):
        self.set_key('dn', data)

    def _get_dn(self):
        return self.get_key('dn')

    def _set_qp(self, data):
        self.set_key('qp', data)

    def _get_qp(self):
        return self.get_key('qp')

    def _set_qs(self, data):
        self.set_key('qs', data)

    def _get_qs(self):
        return self.get_key('qs')

    hl = property(_get_hl, _set_hl)
    vp = property(_get_vp, _set_vp)
    vs = property(_get_vs, _set_vs)
    dn = property(_get_dn, _set_dn)
    qp = property(_get_qp, _set_qp)
    qs = property(_get_qs, _set_qs)

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

    def write(self, ascii_file):
        """
        """
        pass

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

    def _dispersion(self, frequency, mode, itype, ifunc):
        """
        """
        hl = np.array(self.hl) / 1e3
        vp = np.array(self.vp) / 1e3
        vs = np.array(self.vs) / 1e3
        dn = np.array(self.dn) / 1e3

        return surf96(1./frequency, hl, vp, vs, dn,
                      mode=mode, itype=itype, ifunc=ifunc,
                      dc=0.0001, dt=0.025)

    def love_phase_dispersion(self, frequency, mode=0):
        """
        """
        return self._dispersion(frequency, mode, 0, 1)

    def rayleigh_phase_dispersion(self, frequency, mode=0):
        """
        """
        return self._dispersion(frequency, mode, 0, 2)

    def love_group_dispersion(self, frequency, mode=0):
        """
        """
        return self._dispersion(frequency, mode, 1, 1)

    def rayleigh_group_dispersion(self, frequency, mode=0):
        """
        """
        return self._dispersion(frequency, mode, 1, 2)

    def soil_response(self, delta, nsamp, iwave='sh', iangle=0.,
                      ilayer=0, elastic=False):
        """
        Compute the soil (complex) response at the
        surface for outcropping rock reference conditions.
        Calculation can be done for an arbitrary angle of
        incidence (0-90), elastic or anelastic.

        ilayer = measuring layer
        """
        freq = frequency_axis(delta, nsamp)
        out = soil_response(freq, self, iwave=iwave,
                            iangle=iangle, elastic=elastic)

        if iwave == 'sh':
            sp_h = Spectrum(delta=delta, nsamp=nsamp, data=out[ilayer])
            return sp_h

        else:
            sp_h = Spectrum(delta=delta, nsamp=nsamp, data=out[0][ilayer])
            sp_v = Spectrum(delta=delta, nsamp=nsamp, data=out[1][ilayer])
            return sp_h, sp_v


    def site_kappa(self, depth=None):
        """
        Compute the Kappa parameter directly from the site model
        using harmonic averaging. Caculation depth can be specified,
        otherwise depth of the last layer interface is assumed.

        :param float depth:
            averaging depth in meters (optional)
        """

        # Compute kappa attenuation
        kappa = compute_site_kappa(np.array(self.hl),
                                   np.array(self.vs),
                                   np.array(self.qs),
                                   depth)

        return utils.a_round(kappa, DECIMALS)

    def soil_class(self, code='EC8'):
        """
        Compute geotechnical classification according to specified
        building code. Default is EC8 (missing special classes).

        :param string code:
            the reference building code for the classification
            (default EC8)
        """

        return engpar.soil_class(self.vs30, code)
