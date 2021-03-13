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
An simple Python library for SAC file manipulation
"""

from struct import pack, unpack
from os.path import isfile

from shakelab.libutils.time import Date
from shakelab.libutils.time import day_to_month

class Sac(object):

    def __init__(self, file=[], byte_order='le'):
        """
        Base class for sac file I/O

        Usage:
            s = SacLib.Sac
            s = SacLib.Sac('MyFile.sac')

        Optional parameters:
            byte_order - 'le' (Little Endian, default)
                       - 'be' (Big Endian)

        Attributes:
            S.head - SAC header information
            S.data[0] - SAC trace
            S.data[1] - SAC trace (optional block)
            S.byte - byte order

        Methods:
            S.read - Read SAC file from disk
            S.write - write SAC file to disk
            S.copy - Create copy of SAC object
            S.print - Print header information
        """

        # Variable initialisation
        self.head = {}
        self.data = [[], []]

        # Set byte order
        self.byte = byte_order

        if file:
            # Import SAC file
            self.read(file)

        else:
            # Set header defaults
            for H in _HdrStruc:
                self.head[H[0]] = H[3]

    def read(self, file, byte_order='le'):
        """
        Read SAC file from disk

        Usage:
            s.read('MyFile.sac')
            s.read('MyFile.sac', byte_order='be')
        """

        # Set byte order
        self.byte = byte_order

        # Open input SAC file
        with open(file, 'rb') as fid:

            # Import header
            for H in _HdrStruc:
                self.head[H[0]] = _fread(fid, H[1], H[2], self.byte)

            # Import first block of data
            for N in range(0, self.head['NPTS']):
                self.data[0].append(_fread(fid, 4, 'f', self.byte))

            # Import second block of data
            if self.head['LEVEN'] != 1 or self.head['IFTYPE'] in (2, 3):
                for N in range(0, self.head['NPTS']):
                    self.data[1].append(_fread(fid, 4, 'f', self.byte))

            fid.close()
            return

        # Warn user if model file does not exist
        print('Error: File not found')

    def write(self, file, byte_order='le', owrite=False):
        """
        Write SAC file to disk

        Usage:
            s.write('MyFile.sac')
            s.write('MyFile.sac', byte_order='be')
            s.write('MyFile.sac', owrite=True)
        """

        if isfile(file) and not owrite:
            print('Error: file exists, not overwriting....')
            return

        if byte_order:
            # Set byte order (le/be)
            self.byte = byte_order

        # Open output SAC file
        with open(file, 'wb') as fid:

            # Export header
            for H in _HdrStruc:
                _fwrite(fid, self.head[H[0]], H[1], H[2], self.byte)

            # Export first block of data
            for D in self.data[0]:
                _fwrite(fid, D, 4, 'f', self.byte)

            # Export second block of data
            if self.data[1]:
                for D in self.data[1]:
                    _fwrite(fid, D, 4, 'f', self.byte)

            fid.close()
            return

        # Warn user if model file does not exist
        print('Error: File not found')

    def info(self):
        """
        Print header details

        Usage:
            s.print()
        """

        print('------------')
        for H in _HdrStruc:
            data = self.head[H[0]]

            if data != H[3]:
                print('{0:>12} = {1}'.format(H[0], data))

    @property
    def time(self):
        """
        """
        year = self.head['NZYEAR']
        day = self.head['NZJDAY']
        hour = self.head['NZHOUR']
        minute = self.head['NZMIN']
        second = self.head['NZSEC']
        msecond = self.head['NZMSEC'] * 1e-4

        # Convert total days to month/day
        (month, day) = day_to_month(year, day)

        return Date([year, month, day, hour, minute, second + msecond])

    @property
    def seconds(self):
        """
        Rounding to millisecond precision
        """
        return round(self.time.to_second(), 4)

    @property
    def delta(self):
        """
        """
        return self.head['DELTA']

    @property
    def duration(self):
        """
        """
        nsamp = self.header['NPTS']
        srate = self.sampling_rate

        return (nsamp * delta)


def _fread(fid, bnum, bkey, bord):
    """
    Bytewise read.
    """

    hex = fid.read(bnum)

    if bkey == 's':
        bkey = str(bnum) + bkey
    if bord == 'be':
        bkey = '>' + bkey
    if bord == 'le':
        bkey = '<' + bkey

    data = unpack(bkey, hex)[0]

    return data


def _fwrite(fid, data, bnum, bkey, bord):
    """
    Bytewise write.
    """

    if bkey == 's':
        bkey = str(bnum) + bkey
    if bord == 'be':
        bkey = '>' + bkey
    if bord == 'le':
        bkey = '<' + bkey

    hex = pack(bkey, data)

    fid.write(hex)


# INTERNAL: Header Structure (sorted)
#   [0] - Field Name
#   [1] - Length in Bytes
#   [2] - Variable Type
#   [3] - Default Value

_HdrStruc = [('DELTA', 4, 'f', -12345),
             ('DEPMIN', 4, 'f', -12345),
             ('DEPMAX', 4, 'f', -12345),
             ('SCALE', 4, 'f', -12345),
             ('ODELTA', 4, 'f', -12345),
             ('B', 4, 'f', -12345),
             ('E', 4, 'f', -12345),
             ('O', 4, 'f', -12345),
             ('A', 4, 'f', -12345),
             ('INTERNAL1', 4, 'f', -12345),
             ('T0', 4, 'f', -12345),
             ('T1', 4, 'f', -12345),
             ('T2', 4, 'f', -12345),
             ('T3', 4, 'f', -12345),
             ('T4', 4, 'f', -12345),
             ('T5', 4, 'f', -12345),
             ('T6', 4, 'f', -12345),
             ('T7', 4, 'f', -12345),
             ('T8', 4, 'f', -12345),
             ('T9', 4, 'f', -12345),
             ('F', 4, 'f', -12345),
             ('RESP0', 4, 'f', -12345),
             ('RESP1', 4, 'f', -12345),
             ('RESP2', 4, 'f', -12345),
             ('RESP3', 4, 'f', -12345),
             ('RESP4', 4, 'f', -12345),
             ('RESP5', 4, 'f', -12345),
             ('RESP6', 4, 'f', -12345),
             ('RESP7', 4, 'f', -12345),
             ('RESP8', 4, 'f', -12345),
             ('RESP9', 4, 'f', -12345),
             ('STLA', 4, 'f', -12345),
             ('STLO', 4, 'f', -12345),
             ('STEL', 4, 'f', -12345),
             ('STDP', 4, 'f', -12345),
             ('EVLA', 4, 'f', -12345),
             ('EVLO', 4, 'f', -12345),
             ('EVEL', 4, 'f', -12345),
             ('EVDP', 4, 'f', -12345),
             ('MAG', 4, 'f', -12345),
             ('USER0', 4, 'f', -12345),
             ('USER1', 4, 'f', -12345),
             ('USER2', 4, 'f', -12345),
             ('USER3', 4, 'f', -12345),
             ('USER4', 4, 'f', -12345),
             ('USER5', 4, 'f', -12345),
             ('USER6', 4, 'f', -12345),
             ('USER7', 4, 'f', -12345),
             ('USER8', 4, 'f', -12345),
             ('USER9', 4, 'f', -12345),
             ('DIST', 4, 'f', -12345),
             ('AZ', 4, 'f', -12345),
             ('BAZ', 4, 'f', -12345),
             ('GCARC', 4, 'f', -12345),
             ('INTERNAL2', 4, 'f', -12345),
             ('INTERNAL3', 4, 'f', -12345),
             ('DEPMEN', 4, 'f', -12345),
             ('CMPAZ', 4, 'f', -12345),
             ('CMPINC', 4, 'f', -12345),
             ('XMINIMUM', 4, 'f', -12345),
             ('XMAXIMUM', 4, 'f', -12345),
             ('YMINIMUM', 4, 'f', -12345),
             ('YMAXIMUM', 4, 'f', -12345),
             ('UNUSED1', 4, 'f', -12345),
             ('UNUSED2', 4, 'f', -12345),
             ('UNUSED3', 4, 'f', -12345),
             ('UNUSED4', 4, 'f', -12345),
             ('UNUSED5', 4, 'f', -12345),
             ('UNUSED6', 4, 'f', -12345),
             ('UNUSED7', 4, 'f', -12345),
             ('NZYEAR', 4, 'i', -12345),
             ('NZJDAY', 4, 'i', -12345),
             ('NZHOUR', 4, 'i', -12345),
             ('NZMIN', 4, 'i', -12345),
             ('NZSEC', 4, 'i', -12345),
             ('NZMSEC', 4, 'i', -12345),
             ('NVHDR', 4, 'i', 6),
             ('NORID', 4, 'i', -12345),
             ('NEVID', 4, 'i', -12345),
             ('NPTS', 4, 'i', -12345),
             ('INTERNAL4', 4, 'i', -12345),
             ('NWFID', 4, 'i', -12345),
             ('NXSIZE', 4, 'i', -12345),
             ('NYSIZE', 4, 'i', -12345),
             ('UNUSED8', 4, 'i', -12345),
             ('IFTYPE', 4, 'i', 1),
             ('IDEP', 4, 'i', -12345),
             ('IZTYPE', 4, 'i', -12345),
             ('UNUSED9', 4, 'i', -12345),
             ('IINST', 4, 'i', -12345),
             ('ISTREG', 4, 'i', -12345),
             ('IEVREG', 4, 'i', -12345),
             ('IEVTYP', 4, 'i', -12345),
             ('IQUAL', 4, 'i', -12345),
             ('ISYNTH', 4, 'i', -12345),
             ('IMAGTYP', 4, 'i', -12345),
             ('IMAGSRC', 4, 'i', -12345),
             ('UNUSED10', 4, 'i', -12345),
             ('UNUSED11', 4, 'i', -12345),
             ('UNUSED12', 4, 'i', -12345),
             ('UNUSED13', 4, 'i', -12345),
             ('UNUSED14', 4, 'i', -12345),
             ('UNUSED15', 4, 'i', -12345),
             ('UNUSED16', 4, 'i', -12345),
             ('UNUSED17', 4, 'i', -12345),
             ('LEVEN', 4, 'i', 1),
             ('LPSPOL', 4, 'i', 0),
             ('LOVROK', 4, 'i', 1),
             ('LCALDA', 4, 'i', 1),
             ('UNUSED18', 4, 'i', 0),
             ('KSTNM', 8, 's', '-12345  '),
             ('KEVNM', 16, 's', '-12345  -12345  '),
             ('KHOLE', 8, 's', '-12345  '),
             ('KO', 8, 's', '-12345  '),
             ('KA', 8, 's', '-12345  '),
             ('KT0', 8, 's', '-12345  '),
             ('KT1', 8, 's', '-12345  '),
             ('KT2', 8, 's', '-12345  '),
             ('KT3', 8, 's', '-12345  '),
             ('KT4', 8, 's', '-12345  '),
             ('KT5', 8, 's', '-12345  '),
             ('KT6', 8, 's', '-12345  '),
             ('KT7', 8, 's', '-12345  '),
             ('KT8', 8, 's', '-12345  '),
             ('KT9', 8, 's', '-12345  '),
             ('KF', 8, 's', '-12345  '),
             ('KUSER0', 8, 's', '-12345  '),
             ('KUSER1', 8, 's', '-12345  '),
             ('KUSER2', 8, 's', '-12345  '),
             ('KCMPNM', 8, 's', '-12345  '),
             ('KNETWK', 8, 's', '-12345  '),
             ('KDATRD', 8, 's', '-12345  '),
             ('KINST', 8, 's', '-12345  ')]
