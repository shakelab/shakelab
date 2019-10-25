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
Module for basic waveform analysis
"""

import numpy as _np
from shakelab.utils.time import Date

keymap = {0 : 0, 'e' : 0, 'E' : 0, 'ew' : 0, 'EW' : 0, 'x' : 0, 'X' : 0,
          1 : 1, 'n' : 1, 'N' : 1, 'ns' : 1, 'NS' : 1, 'y' : 1, 'Y' : 1,
          2 : 2, 'u' : 2, 'U' : 2, 'ud' : 2, 'UD' : 2, 'z' : 2, 'Z' : 2}


class Record(object):
    """
    Single recording. It can be multichannel, but channels must
    be synchronous and related to the same measurement location.
    """

    def __init__(self):
        self.id = None
        self.delta = None
        self.scale = None
        self.time = Date()
        self.data = None

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[keymap[item]]

class Station(object):
    """
    Single measuring location in space with one or more recordings.
    """

    def __init__(self):
        self.id = None
        self.latitude = None
        self.longitude = None
        self.record = []

    def __len__(self):
        return len(self.recording)

    def __getitem__(self, sliced):
        return self.recording[sliced]

class Array(object):
    """
    Array of measuring locations.
    """

    def __init__(self):
        self.id = None
        self.station = []

    def __len__(self):
        return len(self.station)

    def __getitem__(self, sliced):
        return self.station[sliced]

