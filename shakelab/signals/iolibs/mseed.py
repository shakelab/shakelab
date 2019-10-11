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
    """
    Non-continuous recordings are not (yet) supported
    """

    def __init__(self, file=None, byte_order='le'):

        # Record list initialisation
        self.record = []

        # Set byte order
        self.byte_order = byte_order

        # Import miniSEED file
        if file is not None:
            self.read(file)

    def read(self, file):
        """
        """

        with open(file, 'rb') as fid:

            # Setting last sample from previous record
            # (used in STEIN compression)
            _last = None

            while True:
                # Reading the fixed header (48 bytes)
                byte_stream = ByteStream(self.byte_order)

                if not byte_stream.read(fid, 48):
                    break

                # Initialise new record
                record = Record(_last)

                # Reading header information
                _read_header(record, byte_stream)

                # Reading the blockettes
                data_offset = record.header['OFFSET_TO_BEGINNING_OF_DATA']
                byte_stream.read(fid, data_offset - 48)
                _read_blockette(record, byte_stream)

                # Reading data
                data_len = record.blockette[1000]['DATA_RECORD_LENGTH']
                byte_stream.read(fid, 2**data_len - data_offset)
                _read_data(record, byte_stream)

                _last = record.data[-1]
                self.record.append(record)

                #print(record.header)
                #print(record.blockette)

            data = []
            for r in self.record:
                data += r.data
            plt.plot(data, '-')
            plt.show(block=False)
            #sys.exit()


class Record(object):
    """
    """
    def __init__(self, first=None):
    
        self.header = {}
        self.blockette = {}
        self.data = []
        self.first = first


def _read_header(record, byte_stream):
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

    record.header = {}
    for hs in head_struc:
        record.header[hs[0]] = byte_stream.get(hs[1], hs[2])

def _read_blockette(record, byte_stream):
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

    for nb in range(record.header['NUMBER_OF_BLOCKETTES_TO_FOLLOW']):

        # Blockette code
        block_type = byte_stream.get('H', 2)

        if block_type in block_struc:
            # Offset to the beginning of the next blockette
            blockette = {'OFFSET_NEXT' : byte_stream.get('H', 2)}

            # Loop through blockette specific keys
            for bs in block_struc[block_type]:
                blockette[bs[0]] = byte_stream.get(bs[1], bs[2])
                record.blockette[block_type] = blockette
        else:
            print('Blockette type {0} not supported'.format(block_type))

def _read_data(record, byte_stream):
    """
    Importing data
    """

    data_struc = {0 : ('s', 1),
                  1 : ('h', 2),
                  3 : ('i', 4),
                  4 : ('f', 4)}

    nos = record.header['NUMBER_OF_SAMPLES']
    enc = record.blockette[1000]['ENCODING_FORMAT']

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
                if record.first is None:
                    current = first
                else:
                    current = record.first

            cn = [_nibble(word[0], 15-n) for n in range(16)]

            for i in range(16):
                if cn[i] in [1, 2, 3]:
                    for sample in _splitI32(word[i], cn[i]):
                        current += sample
                        data[cnt] = current
                        cnt += 1

        if current != last:
            raise ValueError('Sample mismatch in record')

    record.data = data[:nos]


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



