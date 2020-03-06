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

import json as _json
import numpy as _np
import scipy.stats as _stat
import scipy.interpolate as _ipol

import matplotlib.pyplot as plt
from abc import ABCMeta, abstractmethod
import xml.etree.cElementTree as xet

GMI_DEFAULT = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


class FragilityModel(metaclass=ABCMeta):
    """
    Base class for a generic fragility model.
    """
    def __init__(self, model_name, gmi_type=None, bounds=(0., 10.)):
        self.id = model_name
        self.gmt = gmi_type
        self.bounds = bounds
        self.damage_state = {}

    @abstractmethod
    def get_poes(self, dsl, gmi):
        """
        Compute the probability of exceedance of a given damage state
        for a specific ground motion intensity level.

        :param dsl:
            Damage state label
        :param gmi:
            Ground motion intensity
        """

    @abstractmethod
    def add_damage_state(self, dsl):
        """
        Add a fragility model (parametric or discrete) for a given
        damage state.
        """


class FragilityModelParametric(FragilityModel):
    """
    Subclass of :class:`FagilityModel` with paramteric format (mean, stdv)
    of the damage states.
    """

    def add_damage_state(self, dsl, mean, stdv):
        self.damage_state[dsl] = (mean, stdv)

    def get_poes(self, dsl, gmi):
        """
        """
        beta, teta = normal_to_lognormal(self.damage_state[dsl][0],
                                         self.damage_state[dsl][1])
        poes = _stat.lognorm.cdf(gmi, beta, scale=teta)
        return poes

    def to_discrete(self, gmi=GMI_DEFAULT):
        """
        """
        fmd = FragilityModelDiscrete(self.id)
        fmd.gmt = self.gmt
        fmd.gmi = _np.array(gmi)
        for dsl in self.damage_state.keys():
            fmd.add_damage_state(dsl, self.get_poes(dsl, gmi))
        return fmd


class FragilityModelDiscrete(FragilityModel):
    """
    Subclass of :class:`FagilityModel` with discrete format of the
    damage states (array of PoE values)
    """

    def add_damage_state(self, dsl, poes):
        self.damage_state[dsl] = _np.array(poes)

    def get_poes(self, dsl, gmi, sampling='lin'):
        if sampling == 'lin':
            scale = lambda x: x
        if sampling == 'log':
            scale = lambda x: _np.log(x)
        fi = _ipol.interp1d(scale(self.gmi), self.damage_state[dsl])
        return fi(scale(gmi))


def normal_to_lognormal(mean, stdv):
    """
    NOTE: This conversion is actually not clear to me and must be
    investigated further.
    """
    beta = _np.sqrt(_np.log(1 + (stdv / mean) ** 2))
    teta = mean / _np.sqrt(1 + (stdv / mean) ** 2)
    return beta, teta


def plot_fragility_model(fm, gmi, file=None):
    """
    """
    plt.figure(figsize=(6, 4), dpi=150)
    for dsl in fm.damage_state.keys():
        x = gmi
        y = fm.get_poes(dsl, gmi)
        plt.plot(x, y, label=dsl)
    plt.title(fm.id)
    plt.legend(loc='lower right')
    plt.xlabel('Ground Motion Intensity')
    plt.ylabel('Probability of Exeedance')
    plt.gca().set_xlim([_np.min(x), _np.max(x)])
    plt.gca().set_ylim([0., 1.])
    plt.grid('on')
    plt.show(block=False)

    if file is not None:
        plt.savefig(file, dpi=150)


class FragilityCollection(object):
    """
    """
    def __init__(self, json_file=None):
        self.model = []
        if json_file is not None:
            self.import_from_json(json_file)

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            model = self.model[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return model

    def add_model(self, fragility_model):
        self.model.append(fragility_model)

    def add_from_dict(self, fmd):
        """
        :param fmd:
            Fragility model as json dictionary (string variables).
        """
        if fmd['format'] == 'parametric':
            fm = FragilityModelParametric(fmd['id'])
            fm.gmt = fmd['gmt']
            fm.bounds = (float(fmd['bounds']['min']),
                         float(fmd['bounds']['max']))
            for dsl in fmd['damage_states']:
                mean = float(dsl['mean'])
                stdv = float(dsl['stdv'])
                fm.add_damage_state(dsl['id'], mean, stdv)

        if fmd['format'] == 'discrete':
            fm = FragilityModelDiscrete(fmd['id'])
            fm.gmt = fmd['gmt']
            fm.bounds = (float(fmd['bounds']['min']),
                         float(fmd['bounds']['max']))
            fm.gmi = _np.array([float(f) for f in fmd['intensity']])
            for dsl in fmd['damage_states']:
                poes = _np.array([float(f) for f in dsl['poes']])
                fm.add_damage_state(dsl['id'], poes)

        self.model.append(fm)

    def import_from_json(self, json_file):
        with open(json_file) as jf:
            data = _json.load(jf)
            for fmd in data['fragility_collection']:
                self.add_from_dict(fmd)

    def to_discrete(self, gmi=GMI_DEFAULT):

        fc = FragilityCollection()
        for fm in self.model:
            if isinstance(fm, FragilityModelParametric):
                fc.add_model(fm.to_discrete(gmi))
        return fc


class TreeItem(object):
    """
    """
    def __init__(self, id=None):
        self.id = id
        self.branch = {}


class TaxonomyTree(object):
    """
    """
    def __init__(self, json_file=None):
        self.tree = []
        if json_file is not None:
            self.import_from_json(json_file)

    def import_from_json(self, json_file):
        with open(json_file) as jf:
            data = _json.load(jf)
            for tte in data['taxonomy_list']:
                self.add_from_dict(tte)

    def add_from_dict(self, tte):
        """
        """
        ti = TreeItem(tte['id'])
        for fm in tte['fragility']:
            ti.branch[fm['id']] = float(fm['weight'])
        self.tree.append(ti)

    def get_element(self, id):
        """
        """
        return next((x for x in self.tree if x.id  == id), None)


class TaxonomyItem(object):
    """
    """
    def __init__(self, id=None, number_of_buildings=None):
        self.id = id
        self.number_of_buildings = number_of_buildings
        self.occupants = {'day': None,
                          'transit': None,
                          'night': None}
        self.cost = {'structural': None,
                     'nonstructural': None,
                     'content': None,
                     'bi': None}


class LocationItem(object):
    """
    """
    def __init__(self, id=None, code=None,
                 latitude=None, longitude=None, area=None):
        self.id = id
        self.code = code
        self.latitude = latitude
        self.longitude = longitude
        self.area = area
        self.taxonomy = []


class ExposureDatabase(object):
    """
    """
    def __init__(self, json_file=None):
        self.location = []
        if json_file is not None:
            self.import_from_json(json_file)

    def add_from_dict(self, exp):
        """
        """
        li = LocationItem(exp['id'],
                          exp['code'],
                          float(exp['latitude']),
                          float(exp['longitude']),
                          float(exp['area']))

        for tax in exp['taxonomy']:

            ti = TaxonomyItem()
            ti.id = tax['id']
            ti.number_of_buildings = float(tax['number_of_buildings'])

            for key, value in tax['occupants'].items():
                ti.occupants[key] = float(value)
            for key, value in tax['cost'].items():
                ti.cost[key] = float(value)

            li.taxonomy.append(ti)

        self.location.append(li)

    def import_from_json(self, json_file):

        with open(json_file) as jf:
            data = _json.load(jf)
            for exp in data['exposure']:
                self.add_from_dict(exp)

