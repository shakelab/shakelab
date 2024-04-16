# ****************************************************************************
#
# Copyright (C) 2019-2023, ShakeLab Developers.
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
import numpy as np
from copy import deepcopy

from shakelab.libutils.time import Date
from shakelab.signals.binutils import ByteStream
from shakelab.signals import base

DEFAULT_BYTE_ORDER = 'be'
ADMITTED_RECORD_LENGTH = [256, 1024, 2048, 4096]
ADMITTED_ENCODING = [0, 1, 3, 4, 10, 11]

def msread(input_data_source, stream_collection=None,
           byte_order=DEFAULT_BYTE_ORDER):
    """
    """
    if stream_collection is None:
        stream_collection = base.StreamCollection()

    record_list = msrawread(input_data_source, byte_order=byte_order)
    for record in record_list:
        stream_collection.append(record.to_shakelab())

    return stream_collection

def msrawread(input_data_source, byte_order=DEFAULT_BYTE_ORDER):
    """
    """
    if isinstance(input_data_source, ByteStream):
        byte_stream = input_data_source
    else:
        byte_stream = ByteStream(byte_order=byte_order)
        byte_stream.ropen(input_data_source)

    record_list = []

    # Loop over records
    while True:
        try:
            record = MSRecord()
            record.read(byte_stream)
            record_list.append(record)
        except:
            print('Invalid record found. Stop reading.')
            break

        # Check if end of stream, otherwise exit
        if byte_stream.offset >= byte_stream.length:
            break

    byte_stream.close()

    return record_list

def msrawwrite(record_list, output_data_source,
               byte_order=DEFAULT_BYTE_ORDER):
    """
    """
    if isinstance(output_data_source, ByteStream):
        byte_stream = output_data_source
    else:
        byte_stream = ByteStream(byte_order=byte_order)
        byte_stream.wopen(output_data_source)

    for record in record_list:
        record.write(byte_stream)

    byte_stream.close()


class MSRecord(object):
    """
    MiniSeed record class
    """
    def __init__(self, byte_stream=None, **kwargs):

        self._header_init()
        self._blockette_init()
        self.data = []

        for k,v in kwargs.items():
            if k in self.header:
                self.header.update({k : v})

        if byte_stream is not None:
            self.read(byte_stream)

    def __len__(self):

        return self.nsamp

    def _header_init(self):
        """
        """
        self.header = {}
        for hs in head_struc:
            self.header[hs[0]] = None

    def _blockette_init(self):
        """
        """
        self.blockette = {1000 : {}}
        for bs in block_struc[1000]:
            self.blockette[1000][bs[0]] = None

    @property
    def nsamp(self):
        """
        """
        return self.header['NUMBER_OF_SAMPLES']

    @property
    def delta(self):
        """
        """
        srate = self.header['SAMPLE_RATE_FACTOR']
        rmult = self.header['SAMPLE_RATE_MULTIPLIER']

        if srate < 0:
            srate = -1./srate
        if rmult < 0:
            rmult = 1./rmult
        srate *= rmult

        return 1./srate

    @property
    def seqn(self):
        """
        Sequence number
        """
        return int(self.header['SEQUENCE_NUMBER'])

    @property
    def time(self):
        """
        """
        date = '{0:04d}-'.format(self.header['YEAR'])
        date += '{0:03d}T'.format(self.header['DAY'])
        date += '{0:02d}:'.format(self.header['HOURS'])
        date += '{0:02d}:'.format(self.header['MINUTES'])
        date += '{0:02d}.'.format(self.header['SECONDS'])
        date += '{0:04d}'.format(self.header['MSECONDS'])

        return Date(date)

    @property
    def code(self):
        """
        Stream identifier (FDSN code)
        """
        net = self.header['NETWORK_CODE'].strip()
        sta = self.header['STATION_CODE'].strip()
        loc = self.header['LOCATION_IDENTIFIER'].strip()
        chn = self.header['CHANNEL_IDENTIFIER'].strip()

        return '{0}.{1}.{2}.{3}'.format(net, sta, loc, chn)

    @property
    def duration(self):
        """
        """
        return (self.nsamp - 1) * self.delta

    @property
    def _bytelen(self):
        """
        """
        return (2**self.blockette[1000]['DATA_RECORD_LENGTH'] -
                self.header['OFFSET_TO_BEGINNING_OF_DATA'])

    def write(self, byte_stream, sequence_number, 
                    record_length=None, encoding=None):
        """
        Write the MSRecord information to the given byte stream.
        Note: Data length might not fill exactly in one binary record,
              therefore two options are possible:
              1) last record is filled with zeros for the remaining bytes
              2) exceeding data is returned to a record for subsequent use
        """
        rec = MSRecord()
        rec.header = deepcopy(self.header)
        rec.blockette = deepcopy(self.blockette)

        if record_length is None:
            record_length = 2**rec.blockette[1000]['DATA_RECORD_LENGTH']
        else:
            base2len = int(np.log2(record_length))
            rec.blockette[1000]['DATA_RECORD_LENGTH'] = base2len        

        if encoding is None:
            encoding = rec.blockette[1000]['ENCODING_FORMAT']
        else:
            if encoding in ADMITTED_ENCODING:
                rec.blockette[1000]['ENCODING_FORMAT'] = int(encoding)
            else:
                raise ValueError('Encoding format not recognized')

        header_size = 48

        total_blockette_size = 0
        for block_type in rec.blockette:
            total_blockette_size += blockette_size(block_type)

        data_offset = header_size + total_blockette_size
        rec.header['OFFSET_TO_BEGINNING_OF_DATA'] = data_offset

        data_length = record_length - data_offset

        record_offset = byte_stream.offset

        data_struc = {0: ('s', 1),
                      1: ('h', 2),
                      3: ('i', 4),
                      4: ('f', 4)}

        if encoding in [0, 1, 3, 4]:
            bnum = data_length//data_struc[encoding][1]
            print(bnum)
            #if bnum > len(self.data):
            #    bnum = len(self.data)
            bnum_left = max([len(self.data)-bnum, 0])
            rec.header['NUMBER_OF_SAMPLES'] = bnum

        # Write header information
        for hs in head_struc:
            byte_stream.put(rec.header[hs[0]], hs[1], hs[2])

        # Absolute offset should be added
        block_offset = rec.header['OFFSET_TO_BEGINNING_OF_BLOCKETTE']

        # Write blockette information
        for block_type, blockette in rec.blockette.items():

            byte_stream.goto(record_offset + block_offset)

            block_offset += blockette_size(block_type)

            # Blockette code
            byte_stream.put(block_type, 'H', 2)

            # Offset to the beginning of the next blockette
            byte_stream.put(block_offset, 'H', 2)

            for bs in block_struc[block_type]:
                byte_stream.put(blockette[bs[0]], bs[1], bs[2])

        # Write data
        if encoding in [0, 1, 3, 4]:

            for data in self.data[0:bnum]:
                byte_stream.put(data,
                                data_struc[encoding][0],
                                data_struc[encoding][1])

            for data in range(0, bnum_left):
                byte_stream.put(0.,
                                data_struc[encoding][0],
                                data_struc[encoding][1])

        elif encoding in [10, 11]:
            # Handle STEIM1/2 encoding
            # Your implementation for STEIM1/2 encoding goes here
            pass

        else:
            raise ValueError('Not recognized data format: ', encoding)

    def read(self, byte_stream):
        """
        """
        # Set record initial offset
        self._record_offset = byte_stream.offset

        # Reading header information
        self._get_header(byte_stream)

        # Reading the blockettes
        self._get_blockette(byte_stream)

        # Reading data
        self._get_data(byte_stream)

    def _get_header(self, byte_stream):
        """
        Importing header structure
        """
        self.header = {}
        for hs in head_struc:
            self.header[hs[0]] = byte_stream.get(hs[1], hs[2])

    def _get_blockette(self, byte_stream):
        """
        Importing blockettes
        """
        block_offset = self.header['OFFSET_TO_BEGINNING_OF_BLOCKETTE']

        for nb in range(self.header['NUMBER_OF_BLOCKETTES_TO_FOLLOW']):

            byte_stream.goto(self._record_offset + block_offset)

            # Blockette code
            block_type = byte_stream.get('H', 2)
            print(block_type)

            # Offset to the beginning of the next blockette
            block_offset = byte_stream.get('H', 2)

            if block_type in block_struc:
                # Blockette initialisation
                blockette = {'OFFSET_NEXT': block_offset}

                # Loop through blockette specific keys
                for bs in block_struc[block_type]:
                    blockette[bs[0]] = byte_stream.get(bs[1], bs[2])
                    self.blockette[block_type] = blockette

            else:
                print('Blockette type {0} not supported'.format(block_type))

    def _get_data(self, byte_stream):
        """
        Importing data
        """
        data_struc = {0: ('s', 1),
                      1: ('h', 2),
                      3: ('i', 4),
                      4: ('f', 4)}

        offset = (self._record_offset +
                  self.header['OFFSET_TO_BEGINNING_OF_DATA'])

        byte_stream.goto(offset)

        nos = self.header['NUMBER_OF_SAMPLES']
        enc = self.blockette[1000]['ENCODING_FORMAT']

        if enc in [0, 1, 3, 4]:

            bnum = self._bytelen//data_struc[enc][1]
            data = [None] * bnum

            # Reading all data bytes, including zeros
            for ds in range(bnum):
                data[ds] = byte_stream.get(data_struc[enc][0],
                                           data_struc[enc][1])
                print(data[ds])

            # Decode ASCII data (e.g. logs)
            if enc == 0:
                data = "".join([d.decode() for d in data])

        elif enc in [10, 11]:

            if byte_stream.byte_order == 'le':
                raise ValueError('STEIM1/2 only defined for Big-Endian')

            cnt = 0
            data = [None] * nos

            for fn in range(self._bytelen//64):

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

    def append(self, record):
        """
        Append the data from a recording to the current.
        Note that sequence number is updated.

        TODO: add header consistency check
        """
        self.header['SEQUENCE_NUMBER'] = record.header['SEQUENCE_NUMBER']
        self.header['NUMBER_OF_SAMPLES'] += record.header['NUMBER_OF_SAMPLES']
        self.data += record.data

    def to_shakelab(self):
        """
        Convert MiniSeed record to Shakelab record object
        """
        record = base.Record()

        record.head.sid = self.code
        record.head.delta = self.delta
        record.head.time = self.time
        record.data = np.array(self.data)

        return record


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

def blockette_size(block_type):
    """
    Return the size a specific blockette type by summing
    the size of individual variables
    """
    if block_type in block_struc:
        size = sum([i[2] for i in block_struc[block_type]])
    else:
        raise ValueError('Not a valid blockette type')

    # Add standard four bytes common to all blockettes.
    size += 4

    return size


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
