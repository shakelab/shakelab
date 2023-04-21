# ****************************************************************************
#
# Copyright (C) 2019-2021, ShakeLab Developers.
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
Module to import / export data formats
"""

import os
import numpy as np
from io import BytesIO
from struct import pack, unpack

#from shakelab.signals.base import Record
#from shakelab.signals.libio import mseed, sac, smdb

from shakelab.libutils.time import Date
from shakelab.libutils.time import days_to_month


class ByteStream(object):
    """
    This class allows reading binary data from a file or from a
    byte buffer using the same format.
    """
    def __init__(self, byte_stream=None, byte_order='be'):

        self.buffer = b''
        self.length = 0

        if byte_stream is not None:
            self.open(byte_stream)

        self.byte_order = byte_order

    def open(self, byte_stream):
        """
        """
        if isinstance(byte_stream, bytes):
            self.buffer = BytesIO(byte_stream)
        else:
            self.buffer = open(byte_stream, 'rb')

        self.goto(0, 2)
        self.length = self.offset
        self.goto(0)

    def goto(self, offset, whence=0):
        """
        offset − position of the read/write pointer within the file.
        whence − 0 for absolute file positioning, 1 for relative to
        the current position and 2 seek relative to the file's end.
        """
        self.buffer.seek(offset, whence)

    def shift(self, places):
        """
        """
        self.goto(places, 1)

    @property
    def offset(self):
        """
        """
        return self.buffer.tell()

    def get(self, byte_format, byte_num=1, offset=None):
        """
        """
        if offset is not None:
            self.goto(offset)

        byte_buffer = self.buffer.read(byte_num)

        if byte_format == 's':
            byte_format = str(byte_num) + byte_format

        byte_map = {'be': '>', 'le': '<'}
        byte_format = byte_map[self.byte_order] + byte_format

        value = unpack(byte_format, byte_buffer)[0]

        if isinstance(value, bytes):
            value = value.decode()

        return value

    def close(self):
        """
        """
        self.buffer.close()


def reader(file, ftype=None, byte_order='be'):
    """
    """

    if ftype is None:
        # Try to identify file from extension
        fext = os.path.splitext(file)[1]

        if fext in ['.ms', '.mseed', '.miniseed', '.seed']:
            ftype = 'mseed'

        elif fext in ['.sac', '.SAC']:
            ftype = 'sac'

        else:
            raise NotImplementedError('format not recognized')

    # Import recordings

    if ftype == 'mseed':
        sc = mseed.read(file, byte_order=byte_order)

    elif ftype == 'sac':

        sc = sac.Sac(file, byte_order=byte_order)
        record = Record()
        record.head.delta= sc.delta
        record.head.time = Date(sc.time, format='julian')
        record.data = np.array(sc.data[0])
        rec_list.append(record)

    elif ftype == 'itaca':

        it = smdb.Itaca(file)
        record = Record()
        record.head.rate = it.sampling_rate()
        record.head.time = Date(it.time)
        record.data = np.array(it.data)
        rec_list.append(record)

    elif ftype == 'ascii':
        raise NotImplementedError('format not yet implemented')

    elif ftype == 'seisan':
        raise NotImplementedError('format not yet implemented')

    elif ftype == 'seg2':
        raise NotImplementedError('format not yet implemented')

    elif ftype == 'dat':
        raise NotImplementedError('format not yet implemented')

    elif ftype == 'gse':
        raise NotImplementedError('format not yet implemented')

    elif ftype == 'reftek':
        raise NotImplementedError('format not yet implemented')

    else:
        pass

    return rec_list

