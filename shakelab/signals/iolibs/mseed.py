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

from struct import pack, unpack


class MiniSeed(object):

    def __init__(self, file=[], byte_order='le'):
        """
        """

        # Variable initialisation
        self.head = {}
        self.data = []

        # Set byte order
        self.byte = byte_order

        if file:
            # Import miniSEED file
            self.read(file)


    #--------------------------------------------------------------------------

    def read(self, file):
        """
        """

        with open(file, 'rb') as fid:

            rec = self.read_record(fid, 0)

    def read_record(self, fid, offset):
        """
        """
        # Record sequence number
        print(_fread(fid, 6, 's', self.byte))

        # Data header/quality indicator
        print(_fread(fid, 1, 's', self.byte))

        # Reserved
        print(_fread(fid, 1, 's', self.byte))

        return []

    def read_header(self, fid):
        """
        """
        a=1

    def read_blockette(self, fid):
        """
        """
        a=1

    def read_data(self, fid):
        """
        """
        a=1


# =============================================================================
# INTERNAL: bytewise read

def _fread(fid, bnum, bkey, bord):

    hex = fid.read(bnum)

    if bkey == 's': bkey = str(bnum) + bkey
    if bord == 'be': bkey = '>' + bkey
    if bord == 'le': bkey = '<' + bkey

    data = unpack(bkey, hex)[0]

    return data


# =============================================================================
# INTERNAL: bytewise write

def _fwrite(fid, data, bnum, bkey, bord):

    if bkey == 's': bkey = str(bnum) + bkey
    if bord == 'be': bkey = '>' + bkey
    if bord == 'le': bkey = '<' + bkey

    hex = pack(bkey, data)

    fid.write(hex)

