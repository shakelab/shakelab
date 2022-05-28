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
Module to import / export data formats
"""

import numpy as np

from shakelab.signals.base import Record
from shakelab.signals import mseed, sac, smdb, fourier

from shakelab.libutils.time import Date
from shakelab.libutils.time import days_to_month


def reader(file, format='mseed', path=None, byte_order='be',
                     **kwargs):
    """
    """
    # Setting path
    if path is not None:
        file = path + file

    # Initialise an empty trace
    rec_list = []

    # Import recordings from file

    if format in ('ms', 'mseed', 'miniseed'):

        ms = mseed.MiniSeed(file, byte_order=byte_order)

        for (code, stream) in ms.stream.items():
            for msrec in stream:

                record = Record()
                record.head.delta = msrec.delta
                record.head.time = Date(msrec.time, format='julian')
                record.head.sid.set(msrec.sid)
                record.data = np.array(msrec.data)

                rec_list.append(record)

    elif format == 'sac':

        sc = sac.Sac(file, byte_order=byte_order)
        record = Record()
        record.head.rate = sc.sampling_rate()
        record.head.time = Date(sc.time)
        record.data = np.array(sc.data[0])
        rec_list.append(record)

    elif format == 'itaca':

        it = smdb.Itaca(file)
        rec = Record()
        rec.dt = it.sampling_rate()
        rec.data = it.time_date()
        rec.data = it.data
        rec_list.append(rec)

    elif format == 'ascii':
        raise NotImplementedError('format not yet implemented')

    elif format == 'seisan':
        raise NotImplementedError('format not yet implemented')

    elif format == 'seg2':
        raise NotImplementedError('format not yet implemented')

    elif format == 'dat':
        raise NotImplementedError('format not yet implemented')

    elif format == 'gse':
        raise NotImplementedError('format not yet implemented')

    else:
        raise NotImplementedError('format not recognized')

    return rec_list

