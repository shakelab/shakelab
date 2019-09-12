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

    def __init__(self, file=None, byte_order='le'):
        """
        """

        # Variable initialisation
        self.head = {}
        self.block = {}
        self.data = []


        # Set byte order
        self.byte = byte_order

        self.offset = 0.

        if file is not None:
            # Import miniSEED file
            self.read(file)

    #--------------------------------------------------------------------------

    def read(self, file):
        """
        """

        with open(file, 'rb') as fid:

            rec = self._read_record(fid)

    def _read_record(self, fid):
        """
        """

        # Importing header structure
        for H in _HdrStruc:
            self.head[H[0]] = _fread(fid, H[1], H[2], self.byte)

        # Importing blockettes
        for nb in range(self.head['NUMBER_OF_BLOCKETTES_TO_FOLLOW']):
            block_type = _fread(fid, 2, 'H', self.byte)
            offset = _fread(fid, 2, 'H', self.byte)

            self.block[block_type] = {}
            for B in BlkDict[block_type]:
                self.block[block_type][B[0]] = _fread(fid, B[1], B[2], self.byte)

        # Preallocate data array
        self.data = [None] * self.head['NUMBER_OF_SAMPLES']

        # Importing data
        ef = self.block[1000]['ENCODING_FORMAT']

        if ef in [0, 1, 3, 4, 5]:
            for s in range(self.head['NUMBER_OF_SAMPLES']):
                self.data[s] = _fread(fid, T[ef][0], T[ef][1], self.byte)

        if ef is 10:

            # for s in range(self.head['NUMBER_OF_SAMPLES']):
            #     self.data[s] = _fread(fid, 4, 'I', self.byte)

            frame = [None] * 16
            for w in range(16):
                frame[w] = _fread(fid, 4, 'I', self.byte)

            cn = [_nibble(frame[0], n) for n in range(16)]
            x0 = frame[1]
            xn = frame[2]

            for i in range(16):
                if cn[i] is 0:
                    print('no difference')
                if cn[i] is 1:
                    print('1byte difference')
                if cn[i] is 2:
                    print('2byte difference')
                if cn[i] is 3:
                    print('4byte difference')

            """
            frame_number = self.head['NUMBER_OF_SAMPLES'] // 16
            reminder = self.head['NUMBER_OF_SAMPLES'] % 16

            for f in range(frame_number):
                frame = [None] * 16
                for w in range(16):
                    frame[w] = _fread(fid, 4, 'I', self.byte)
                if f is 0:
                    x0 = frame[1]
                    xn = frame[2]
            """

        #c0 = [_nibble(self.data[0], n) for n in range(16)]
        #print([bin(c) for c in c0])

        print(self.head)
        print(self.block)
        print(self.data)
        print(fid.tell())

# =============================================================================
# INTERNAL: Extract 2bits nibble from integer

def _nibble(value, position):

    return (value >> 2*position) & 3

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

# =============================================================================
# INTERNAL: Header Structure

_HdrStruc = [('SEQUENCE_NUMBER', 6, 's'),
             ('DATA_HEADER_QUALITY_INDICATOR', 1, 's'),
             ('RESERVED_BYTE', 1, 's'),
             ('STATION_CODE', 5, 's'),
             ('LOCATION_IDENTIFIER', 2, 's'),
             ('CHANNEL_IDENTIFIER', 3, 's'),
             ('NETWORK_CODE', 2, 's'),
             ('YEAR', 2, 'H'),
             ('DAY', 2, 'H'),
             ('HOURS', 1, 'B'),
             ('MINUTES', 1, 'B'),
             ('SECONDS', 1, 'B'),
             ('UNUSED', 1, 'B'),
             ('MSECONDS', 2, 'H'),
             ('NUMBER_OF_SAMPLES', 2, 'H'),
             ('SAMPLE_RATE_FACTOR', 2, 'h'),
             ('SAMPLE_RATE_MULTIPLIER', 2, 'h'),
             ('ACTIVITY_FLAGS', 1, 'B'),
             ('IO_FLAGS', 1, 'B'),
             ('DATA_QUALITY_FLAGS', 1, 'B'),
             ('NUMBER_OF_BLOCKETTES_TO_FOLLOW', 1, 'B'),
             ('TIME_CORRECTION', 4, 'l'),
             ('OFFSET_TO_BEGINNING_OF_DATA', 2, 'H'),
             ('OFFSET_TO_BEGINNING_OF_BLOCKETTE', 2, 'H')]

BlkDict =  {1000 : [('ENCODING_FORMAT', 1, 'B'),
                    ('WORD_ORDER', 1, 'B'),
                    ('DATA_RECORD_LENGTH', 1, 'B'),
                    ('RESERVED', 1, 'B')],
            1001 : [('TIMING_QUALITY', 1, 'B'),
                    ('MICRO_SEC', 1, 'B'),
                    ('RESERVED', 1, 'B'),
                    ('FRAME_COUNT', 1, 'B')]}

T = {0 : (1, 's'),
     1 : (1, 'h'),
     3 : (4, 'i'),
     4 : (4, 'f')}
