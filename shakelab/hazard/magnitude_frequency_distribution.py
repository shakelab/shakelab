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

import numpy as np


def moment_to_magnitude(moment):

    # Using Hanks & Kanamori equation
    return ((2./3.) * np.log10(moment)) - 10.7

def magnitude_to_moment(magnitude):

    # Using Hanks & Kanamori equation
    return 10**((3./2.) * magnitude + 10.7)

class MagnitudeOccurrencePair():

    def __init__(self, magnitude=None, occurrence=None):
        self.magnitude = magnitude
        self.occurrence = occurrence

    def get_total_moment(self):
        # Using Hanks & Kanamori equation
        return magnitude_to_moment(self.magnitude) * self.occurrence

class MagnitudeDistribution():

    def __init__(self):
        self.mfd = []

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            magnitude_rate_pair = self.mfd[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return magnitude_rate_pair

    def add(self, magnitude_rate_pair):
        self.mfd.append(magnitude_rate_pair)

    def get_magnitude_list(self):
        return [m[0] for m in self.mfd]

    def get_rate_list(self):
        return [m[1] for m in self.mfd]

    def get_total_moment(self):
        pass

class BoundedGutembergRichter():

    def __init__(self, a_value, b_value, min_mag, max_mag, ref_mag=0.):
        self.a_value = a_value
        self.b_value = b_value
        self.min_mag = min_mag
        self.max_mag = max_mag

    def cumulative_rates(self, magnitude=None):

        if magnitude is None:
            magnitude = self.min_mag

        p1 = 10.**self.a_value
        p2 = 10.**(-self.b_value * magnitude)
        p3 = 10.**(-self.b_value * self.max_mag)

        cumrts = p1 * (p2 - p3)

        cumrts[magnitude > self.max_mag] = 0.
        return cumrts

    def inverse_sampling(self, snum):
        """
        Samples are derived using the inverse transform sampling (also known
        as inversion sampling, inverse probability integral transform)
        """

        um = np.random.rand(snum)
        cf = 1.-(10.**(-self.b_value*(self.max_mag-self.min_mag)))
        return self.min_mag-np.log10(1.-(um*cf))/self.b_value





