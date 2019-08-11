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
An simple Python library for MiniSeeed file manipulation
"""


class MiniSeed(object):

    def __init__(self, file=[]):
        """
        """

        if file:
            # Import miniSEED file
            self.read(file)


    #--------------------------------------------------------------------------

    def read(self, file):
        """
        """

        with open(file, 'rb') as fid:

            rec = read_record(fid, offset)


# =============================================================================

def read_record(fid, offset):
    """
    """
    return []

def read_header(fid):
    """
    """
    a=1

def read_blockette(fid):
    """
    """
    a=1

def read_data():
    """
    """
    a=1



