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
import sys
import matplotlib.pyplot as plt


class MiniSeed(object):

    def __init__(self, file=None, byte_order='le'):
        """
        """

        # Variable initialisation
        self.header = {}
        self.blockette = {}
        self.data = []
        self._isfirst = True

        # Set byte order
        self.byte_order = byte_order

        # Import miniSEED file
        if file is not None:
            self.read(file)

    def read(self, file):
        """
        """

        with open(file, 'rb') as fid:

            while True:
                # Reading the fixed header (48 bytes)
                byte_stream = ByteStream(self.byte_order)

                if not byte_stream.read(fid, 48):
                    break
                self._read_header(byte_stream)

                # Reading the blockettes
                data_offset = self.header['OFFSET_TO_BEGINNING_OF_DATA']
                byte_stream.read(fid, data_offset - 48)
                self._read_blockette(byte_stream)

                # Reading data
                data_len = self.blockette[1000]['DATA_RECORD_LENGTH']
                byte_stream.read(fid, 2**data_len - data_offset)

                self._read_data(byte_stream)

                #print(self.header)
                #print(self.blockette)

            #print(self.data)
            plt.plot(self.data, '-')
            plt.show(block=False)

            #sys.exit()

    def _read_header(self, byte_stream):
        """
        Importing header structure
        """

        head_struc = [('SEQUENCE_NUMBER', 's', 6),
                      ('DATA_HEADER_QUALITY_INDICATOR', 's', 1),
                      ('RESERVED_BYTE', 's', 1),
                      ('STATION_CODE', 's', 5),
                      ('LOCATION_IDENTIFIER', 's', 2),
                      ('CHANNEL_IDENTIFIER', 's', 3),
                      ('NETWORK_CODE', 's', 2),
                      ('YEAR', 'H', 2),
                      ('DAY', 'H', 2),
                      ('HOURS', 'B', 1),
                      ('MINUTES', 'B', 1),
                      ('SECONDS', 'B', 1),
                      ('UNUSED', 'B', 1),
                      ('MSECONDS', 'H', 2),
                      ('NUMBER_OF_SAMPLES', 'H', 2),
                      ('SAMPLE_RATE_FACTOR', 'h', 2),
                      ('SAMPLE_RATE_MULTIPLIER', 'h', 2),
                      ('ACTIVITY_FLAGS', 'B', 1),
                      ('IO_FLAGS', 'B', 1),
                      ('DATA_QUALITY_FLAGS', 'B', 1),
                      ('NUMBER_OF_BLOCKETTES_TO_FOLLOW', 'B', 1),
                      ('TIME_CORRECTION', 'l', 4),
                      ('OFFSET_TO_BEGINNING_OF_DATA', 'H', 2),
                      ('OFFSET_TO_BEGINNING_OF_BLOCKETTE', 'H', 2)]

        self.header = {}
        for hs in head_struc:
            self.header[hs[0]] = byte_stream.get(hs[1], hs[2])

    def _read_blockette(self, byte_stream):
        """
        Importing blockettes
        """

        block_struc =  {1000 : [('ENCODING_FORMAT', 'B', 1),
                                ('WORD_ORDER', 'B', 1),
                                ('DATA_RECORD_LENGTH', 'B', 1),
                                ('RESERVED', 'B', 1)],
                        1001 : [('TIMING_QUALITY', 'B', 1),
                                ('MICRO_SEC', 'B', 1),
                                ('RESERVED', 'B', 1),
                                ('FRAME_COUNT', 'B', 1)]}

        for nb in range(self.header['NUMBER_OF_BLOCKETTES_TO_FOLLOW']):

            # Blockette code
            block_type = byte_stream.get('H', 2)

            # Offset to the beginning of the next blockette
            offset_next = byte_stream.get('H', 2)

            blockette = {}
            for bs in block_struc[block_type]:
                blockette[bs[0]] = byte_stream.get(bs[1], bs[2])
                self.blockette[block_type] = blockette

    def _read_data(self, byte_stream):
        """
        Importing data
        """

        data_struc = {0 : ('s', 1),
                      1 : ('h', 2),
                      3 : ('i', 4),
                      4 : ('f', 4)}

        nos = self.header['NUMBER_OF_SAMPLES']
        enc = self.blockette[1000]['ENCODING_FORMAT']

        data = [None] * nos

        if enc in [0, 1, 3, 4]:

            for ds in range(nos):
                data[ds] = byte_stream.get(data_struc[enc][0],
                                           data_struc[enc][1])

        if enc is 10:

            cnt = 0
            for fn in range(byte_stream.len()//64):

                word = [None] * 16
                for wn in range(16):
                    word[wn] = byte_stream.get('i', 4)

                if fn is 0:
                    first = word[1]
                    last = word[2]

                    if self._isfirst:
                        current = first
                        self._isfirst = False
                    else:
                        current = self.data[-1]

                cn = [_nibble(word[0], 15-n) for n in range(16)]

                for i in range(16):
                    if cn[i] in [1, 2, 3]:
                        for sample in _splitI32(word[i], cn[i]):
                            current += sample
                            data[cnt] = current
                            cnt += 1

            if current != last:
                raise ValueError('Sample mismatch in record')

        self.data += data[:nos]

# =============================================================================
# INTERNAL: binary operations

class ByteStream(object):
    """
    """

    def __init__(self, byte_order='le'):

        self.data = []
        self.offset = 0
        self.byte_order = byte_order

    def read(self, fid, byte_num):
        """
        """
        self.data = fid.read(byte_num)
        self.offset = 0

        if self.data:
            return True
        else:
            return False

    def get(self, byte_key, byte_num=1):
        """
        """
        byte_pack = self.data[self.offset:self.offset + byte_num]
        self.offset += byte_num

        if byte_key == 's':
            byte_key = str(byte_num) + byte_key

        byte_map = {'be' : '>', 'le' : '<'}
        byte_key = byte_map[self.byte_order] + byte_key

        return unpack(byte_key, byte_pack)[0]

    def len(self):
        """
        """
        return len(self.data)


def _nibble(value, position):
    """
    Extract 2bits nibble from integer
    NOTE: 3 corresponds to the mask int('00000011', 2)
    """

    return (value >> 2*position) & 3


def _splitI32(value, order):
    """
    Split a long-word (32 bits) into 2 or 4 integers of
    respectively 16 and 8 bits.
    """

    if order is 1:
        out = [None] * 4
        for i in [0, 1, 2, 3]:
            s = (value >> 8*i) & 255
            out[3-i] = s if s < 2**7 else s - 2**8

    if order is 2:
        out = [None] * 2
        for i in [0, 1]:
            s = (value >> 16*i) & 65535
            out[1-i] = s if s < 2**15 else s - 2**16

    if order is 3:
        out = [value]

    return out



