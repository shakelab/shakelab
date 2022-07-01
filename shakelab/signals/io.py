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

from shakelab.signals.base import Record
from shakelab.signals.libio import mseed, sac, smdb

from shakelab.libutils.time import Date
from shakelab.libutils.time import days_to_month


def reader(file, ftype=None, byte_order='be'):
    """
    """
    # Initialise an empty trace
    rec_list = []

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

        ms = mseed.MiniSeed(file, byte_order=byte_order)

        for (code, stream) in ms.stream.items():
            for msrec in stream:

                record = Record()
                record.head.delta = msrec.delta
                record.head.time = Date(msrec.time, format='julian')
                record.head.sid.set(msrec.sid)
                record.data = np.array(msrec.data)
                rec_list.append(record)

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


def RecDatabase(object):
    """
    """

    def __init__(self, **kwargs):
        self.record = []


