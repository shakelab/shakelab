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
Utilities for binary file handling
"""

import numpy as np
from io import BytesIO
from struct import pack, unpack

DEFAULT_BYTE_ORDER = 'be'


class ByteStream(object):
    """
    This class allows reading binary data from a file or from a
    byte buffer using the same format.
    """
    def __init__(self, byte_order=DEFAULT_BYTE_ORDER):

        self.buffer = b''
        self.byte_order = byte_order
        self.mode = None

    def ropen(self, byte_stream):
        """
        """
        self.mode = 'r'
        if isinstance(byte_stream, bytes):
            self.buffer = BytesIO(byte_stream)
        else:
            self.buffer = open(byte_stream, 'rb')

    def wopen(self, byte_stream=None):
        """
        """
        self.mode = 'w'
        if byte_stream is None:
            self.buffer = BytesIO()
        else:
            self.buffer = open(byte_stream, 'wb')

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

    @property
    def length(self):
        """
        """
        offset = self.offset
        self.goto(0, 2)
        length = self.offset
        self.goto(offset)

        return length

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

    def put(self, value, byte_format, byte_num=1, offset=None):
        """
        Write data to the byte stream.
        """
        # If offset is provided, move the stream pointer to that position
        if offset is not None:
            self.goto(offset)

        # If the byte_format is a string, prepend it with
        # the number of bytes to read
        if byte_format == 's':
            byte_format = str(byte_num) + byte_format
            value = str(value).encode('utf-8')
        if byte_format == 'i':
            value = int(value)
        if byte_format == 'h':
            value = int(value)
        if byte_format == 'f':
            value = float(value)

        # Determine the byte order for packing based on
        # the byte_order attribute
        byte_map = {'be': '>', 'le': '<'}
        byte_format = byte_map[self.byte_order] + byte_format

        # Pack the value into bytes
        packed_bytes = pack(byte_format, value)

        # Write the packed bytes into the buffer
        self.buffer.write(packed_bytes)

    def flush(self):
        """
        """
        self.buffer.flush()

    def close(self):
        """
        """
        self.buffer.close()
        self.mode = None
