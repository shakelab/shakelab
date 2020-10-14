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
Seedlink implementation
"""

class RingBuffer:
    """
    Modified from the example of Sébastien Keim
    """
    
    def __init__(self, size):
        self.size = size
        self.data = [None] * size
        self._i0 = 0

    def append_value(self, data):
        """
        Append an element overwriting the oldest one.
        """
        self.data[self._i0] = data
        self._i0 = (self._i0+1) % self.size

    def append_array(self, data):
        """
        """
        for d in data:
            self.append_value(d)

    def get_buffer(self):
        """
        Return list of elements in correct order
        """
        return self.data[self._i0:]+self.data[:self._i0]

