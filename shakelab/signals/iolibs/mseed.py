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

from shakelab.utils.time import Date

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

            # Last sample from previous record
            # (used in STEIN compression)
            last_sample = None

            # Loop over records
            while True:

                # Initialise stream object
                byte_stream = ByteStream(self.byte_order)

                # Reading the fixed header (48 bytes)
                if not byte_stream.read(fid, 48):
                    break

                # Initialise new record
                record = Record(last_sample)

                # Reading header information
                record.read_header(byte_stream)

                # Reading the blockettes
                data_offset = record.get_data_offset()
                byte_stream.read(fid, data_offset - 48)
                record.read_blockette(byte_stream)

                # Reading data
                data_length = record.get_data_length()
                byte_stream.read(fid, 2**data_length - data_offset)
                record.read_data(byte_stream)

                last_sample = record.data[-1]
                self.record.append(record)

                #print(record.header)
                #print(record.blockette)

            data = []
            for r in self.record:
                data += r.data
            plt.plot(data, '-')
            plt.show(block=False)
            #sys.exit()

# =============================================================================
# Record related methods

class Record(object):
    """
    """
    def __init__(self, first=None):
    
        self.header = {}
        self.blockette = {}
        self.data = []
        self.first = first


    def read_header(self, byte_stream):
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

    def read_blockette(self, byte_stream):
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

            if block_type in block_struc:
                # Offset to the beginning of the next blockette
                blockette = {'OFFSET_NEXT' : byte_stream.get('H', 2)}

                # Loop through blockette specific keys
                for bs in block_struc[block_type]:
                    blockette[bs[0]] = byte_stream.get(bs[1], bs[2])
                    self.blockette[block_type] = blockette
            else:
                print('Blockette type {0} not supported'.format(block_type))

    def read_data(self, byte_stream):
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

        if enc in [10, 11]:

            cnt = 0
            for fn in range(byte_stream.len()//64):

                word = [None] * 16
                for wn in range(16):
                    word[wn] = byte_stream.get('i', 4)

                if fn is 0:
                    first = word[1]
                    last = word[2]
                    if self.first is None:
                        current = first
                    else:
                        current = self.first

                cn = [_binmask(word[0], 2, 15-n) for n in range(16)]

                for i in range(16):
                    if cn[i] in [1, 2, 3]:
                        for diff in _w32split(word[i], cn[i], enc):
                            current += diff
                            data[cnt] = current
                            cnt += 1

            if current != last:
                raise ValueError('Sample mismatch in record')

        self.data = data[:nos]

    def get_data_offset(self):
        """
        """
        return self.header['OFFSET_TO_BEGINNING_OF_DATA']

    def get_data_length(self):
        """
        """
        return self.blockette[1000]['DATA_RECORD_LENGTH']

    def get_time_start(self):
        """
        """
        year = self.header['YEAR']
        day = self.header['DAY']
        hour = self.header['HOURS']
        minute = self.header['MINUTES']
        second = self.header['SECONDS']

    def get_time_end(self):
        """
        """
        a = 1

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


def _binmask(word, bits, position):
    """
    Extract N-bits nibble from long word and convert it to integer
    """
    return (word >> bits * position) & (2**(bits) - 1)

def _getdiff(word, bits, dnum):
    """
    """
    out = [None] * dnum
    for i in range(dnum):
        s = _binmask(word, bits, i)
        out[(dnum-1)-i] = s if s < 2**(bits-1) else s - 2**(bits)
    return out

def _w32split(word, order, scheme):
    """
    Split a long-word (32 bits) into integers of different
    lenght in bits (depending on STEIM1 or STEIM2 scheme)
    """

    if order is 1:
        out = _getdiff(word, 8, 4)

    # STEIM 1
    if scheme is 10:
        if order is 2:
            out = _getdiff(word, 16, 2)

        if order is 3:
            out = [word]

    # STEIM 2
    if scheme is 11:
        dnib = _binmask(word, 2, 15)

        if order is 2:
            if dnib is 1:
                out = _getdiff(word, 30, 1)

            if dnib is 2:
                out = _getdiff(word, 15, 2)

            if dnib is 3:
                out = _getdiff(word, 10, 3)

        if order is 3:
            if dnib is 0:
                out = _getdiff(word, 6, 5)

            if dnib is 1:
                out = _getdiff(word, 5, 6)

            if dnib is 2:
                # TO CHECK!
                out = _getdiff(word, 4, 7)

    return out



