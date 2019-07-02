# =============================================================================
#
# Copyright (C) 2019 Valerio Poggi.
# This file is part of QuakeLab.
#
# QuakeLab is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# QuakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# =============================================================================
"""
Module for basic waveform analysis
"""

import numpy as _np


class Recording(object):
    """
    """

    def __init__(self):
        self.head = {'id' : '',
                     'dt' : -1,
                     'time' : [0,0,0,0,0,0.],
                     'channels' : 1,
                     'units' : ''}
        self.data = {0 : None}

