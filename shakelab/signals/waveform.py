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


class Recording(object):
    """
    Single recording, split or not into multiple windows.
    It can be multichannel, but channels must be synchronous,
    sampled uniformly and related to the same measurement location.
    """

    def __init__(self):
        self.id = None
        self.dt = None
        self.units = None
        self.time = Date()
        self.channel = []

    def __len__(self):
        return len(self.channel)

    def __getitem__(self, sliced):
        return self.channel[sliced]

class Station(object):
    """
    Single measuring location in space with one or more recordings.
    """

    def __init__(self):
        self.id = None
        self.name = None
        self.latitude = None
        self.longitude = None
        self.recording = []

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

