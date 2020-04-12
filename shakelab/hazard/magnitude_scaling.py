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

class wells_coppersmith_1994():
    """
    Wells and Coppersmith (1994, Bull. Seism. Soc. Am.)
    magnitude-scaling relationship.
    """

    def __init__(self):
        """
        """

    def area_to_magnitude(area, type='generic')
        """
        note: uncertainty on magnitude
        """
        if type in ['g', 'generic']:
            mag = 4.07 + 0.98*np.log(area)
            std = 0.24

        if type in ['n', 'normal']:
            mag = 3.93 + 1.02*np.log(area)
            std = 0.25

        if type in ['r', 'reverse']:
            mag = 4.33 + 0.90*np.log(area)
            std = 0.25

        if type in ['s', 'strike-slip']:
            mag = 3.98 + 1.02*np.log(area)
            std = 0.23

        else:
            print('type not recognized')

        return mag, std

    def magnitude_to_area(area, type='generic'):
        """
        note: uncertainty on log area
        """
        if type in ['g', 'generic']:
            area = 10**(-3.49 + 0.91*mag)
            std = 0.24

        if type in ['n', 'normal']:
            area = 10**(-3.42 + 0.90*mag)
            std = 0.22

        if type in ['r', 'reverse']:
            area = 10**(-3.99 + 0.98*mag)
            std = 0.26

        if type in ['s', 'strike-slip']:
            area = 10**(-2.87 + 0.82*mag)
            std = 0.22

        return mag, std