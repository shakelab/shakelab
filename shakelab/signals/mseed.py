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
An simple Python library for MiniSeeed file manipulation
"""

from struct import pack, unpack
from shakelab.libutils.time import Date
from shakelab.libutils.time import day_to_month

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

block_struc = {1000: [('ENCODING_FORMAT', 'B', 1),
                      ('WORD_ORDER', 'B', 1),
                      ('DATA_RECORD_LENGTH', 'B', 1),
                      ('RESERVED', 'B', 1)],
               1001: [('TIMING_QUALITY', 'B', 1),
                      ('MICRO_SEC', 'B', 1),
                      ('RESERVED', 'B', 1),
                      ('FRAME_COUNT', 'B', 1)]}

data_struc = {0: ('s', 1),
              1: ('h', 2),
              3: ('i', 4),
              4: ('f', 4)}


class MiniSeed(object):
    """
    """

    def __init__(self, file=None, byte_order='le'):

        # Record list initialisation
        self.record = []

        # Set byte order
        self.byte_order = byte_order

        # Import miniSEED file
        if file is not None:
            self.read(file)

    def read(self, file, byte_order=None):
        """
        """
        self.record = []

        # Set byte order
        if byte_order is not None:
            self.byte_order = byte_order

        # Initialise stream object
        # (This approach will pre-load the whole file into memory)
        byte_stream = ByteStream(byte_order=self.byte_order)
        byte_stream.read(file)

        # Loop over records
        while True:

            # Initialise new record
            record = Record()
            record.read_stream(byte_stream)

            # Split record in case of multiplexing
            # and non-contiguous data
            if not self.record:
                hs_0 = record.header_set()
                tm_0 = record.time_seconds()
                dt_0 = record.duration()

                self.record.append(record)

            else:
                hs_1 = record.header_set()
                tm_1 = record.time_seconds()
                dt_1 = record.duration()

                if (hs_0 == hs_1) and (tm_0 + dt_0 == tm_1):
                    self.record[-1].data += record.data
                else:
                    self.record.append(record)

                hs_0 = hs_1
                tm_0 = tm_1
                dt_0 = dt_1

            # Check if last record, otherwise exit
            if byte_stream.offset >= byte_stream.length:
                break

    def write(self, file, byte_order=None):
        """
        """
        pass

class Record(object):
    """
    """
    def __init__(self, byte_stream=None):

        self.header = {}
        self.blockette = {}
        self.data = []

        if byte_stream is not None:
            self.read_stream(byte_stream)

    def read_stream(self, byte_stream):
        """
        """
        # Reading header information
        self.read_header(byte_stream)

        # Reading the blockettes
        self.read_blockette(byte_stream)

        # Reading data
        self.read_data(byte_stream)

    def read_header(self, byte_stream):
        """
        Importing header structure
        """
        self._record_offset = byte_stream.offset

        self.header = {}
        for hs in head_struc:
            self.header[hs[0]] = byte_stream.get(hs[1], hs[2])

    def read_blockette(self, byte_stream):
        """
        Importing blockettes
        """
        byte_stream.offset = (self._record_offset +
                              self.header['OFFSET_TO_BEGINNING_OF_BLOCKETTE'])

        for nb in range(self.header['NUMBER_OF_BLOCKETTES_TO_FOLLOW']):

            # Blockette code
            block_type = byte_stream.get('H', 2)

            # Offset to the beginning of the next blockette
            offset_next = byte_stream.get('H', 2)

            if block_type in block_struc:
                # Blockette initialisation
                blockette = {'OFFSET_NEXT': offset_next}

                # Loop through blockette specific keys
                for bs in block_struc[block_type]:
                    blockette[bs[0]] = byte_stream.get(bs[1], bs[2])
                    self.blockette[block_type] = blockette
            else:
                print('Blockette type {0} not supported'.format(block_type))
                byte_stream.offset += offset_next

    def read_data(self, byte_stream):
        """
        Importing data
        """
        byte_stream.offset = (self._record_offset +
                              self.header['OFFSET_TO_BEGINNING_OF_DATA'])

        nos = self.header['NUMBER_OF_SAMPLES']
        enc = self.blockette[1000]['ENCODING_FORMAT']

        if enc in [0, 1, 3, 4]:

            bnum = self.data_length()//data_struc[enc][1]
            data = [None] * bnum

            # Reading all data bytes, including zeros
            for ds in range(bnum):
                data[ds] = byte_stream.get(data_struc[enc][0],
                                           data_struc[enc][1])

        elif enc in [10, 11]:

            cnt = 0
            data = [None] * nos

            for fn in range(self.data_length()//64):

                word = [None] * 16
                for wn in range(16):
                    word[wn] = byte_stream.get('i', 4)

                # First and last sample
                if fn == 0:
                    first = word[1]
                    last = word[2]

                # Extract nibbles
                cn = [_binmask(word[0], 2, 15-n) for n in range(16)]

                # Collecting differences
                for i in range(16):
                    if cn[i] in [1, 2, 3]:
                        for diff in _w32split(word[i], cn[i], enc):
                            data[cnt] = diff
                            cnt += 1

            # Last sample from previous record (for check)
            self.previous = first - data[0]

            # Computing full samples from differences
            data[0] = first
            for i in range(1, nos):
                data[i] += data[i-1]

            if data[-1] != last:
                raise ValueError('Sample mismatch in record')

        else:
            raise ValueError('Not recognized data format: ', enc)

        # Store data
        self.data = data[:nos]

    def data_length(self):
        """
        """
        return (2**self.blockette[1000]['DATA_RECORD_LENGTH'] -
                self.header['OFFSET_TO_BEGINNING_OF_DATA'])

    def time_date(self):
        """
        """
        year = self.header['YEAR']
        day = self.header['DAY']
        hour = self.header['HOURS']
        minute = self.header['MINUTES']
        second = self.header['SECONDS']
        msecond = self.header['MSECONDS'] * 1e-4

        # Convert total days to month/day
        (month, day) = day_to_month(year, day)

        return Date([year, month, day, hour, minute, second + msecond])

    def time_seconds(self):
        """
        """
        date = self.time_date()
        return round(date.to_seconds(), 4)

    def sampling_rate(self):
        """
        """
        srate = self.header['SAMPLE_RATE_FACTOR']
        rmult = self.header['SAMPLE_RATE_MULTIPLIER']

        if srate < 0:
            srate = -1./srate
        if rmult < 0:
            rmult = 1./rmult
        srate *= rmult

        return srate

    def duration(self):
        """
        """
        nsamp = self.header['NUMBER_OF_SAMPLES']
        srate = self.sampling_rate()

        return (nsamp/srate)

    def header_set(self):
        """
        List the header itmes which must be equal between
        consecutive records from the same stream
        """
        items = [3, 4, 5, 6]
        return set([self.header[head_struc[i][0]] for i in items])

    def decode(self):
        """
        """
        return "".join([d.decode() for d in self.data])


class ByteStream(object):
    """
    """

    def __init__(self, data=[], byte_order='le'):

        self.offset = 0
        self.data = data
        self.byte_order = byte_order
        self.length = len(self.data)

    def read(self, file):
        """
        """
        self.offset = 0
        with open(file, 'rb') as fid:
            self.data = fid.read()
        self.length = len(self.data)

    def get(self, byte_key, byte_num=1):
        """
        """
        byte_pack = self.data[self.offset:self.offset + byte_num]
        self.offset += byte_num

        if byte_key == 's':
            byte_key = str(byte_num) + byte_key

        byte_map = {'be': '>', 'le': '<'}
        byte_key = byte_map[self.byte_order] + byte_key

        return unpack(byte_key, byte_pack)[0]

    def reminder(self):
        """
        """
        return (self.length - self.offset)


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
    if order == 1:
        out = _getdiff(word, 8, 4)

    # STEIM 1
    if scheme == 10:
        if order == 2:
            out = _getdiff(word, 16, 2)

        if order == 3:
            out = _getdiff(word, 32, 1)

    # STEIM 2
    if scheme == 11:
        dnib = _binmask(word, 2, 15)

        if order == 2:
            if dnib == 0:
                raise ValueError('Nibble not recognized')

            elif dnib == 1:
                out = _getdiff(word, 30, 1)

            elif dnib == 2:
                out = _getdiff(word, 15, 2)

            elif dnib == 3:
                out = _getdiff(word, 10, 3)

        if order == 3:
            if dnib == 0:
                out = _getdiff(word, 6, 5)

            elif dnib == 1:
                out = _getdiff(word, 5, 6)

            elif dnib == 2:
                # TO CHECK!
                out = _getdiff(word, 4, 7)

            elif dnib == 3:
                raise ValueError('Nibble not recognized')

    return out
