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


def reader(file, format='sac', path=None, byte_order='le',
                     **kwargs):
    """
    """

    # Initialise an empty trace
    rec_list = []

    # Import recordings from file
    if format == 'mseed':
        ms = mseed.MiniSeed(file, byte_order=byte_order)
        for mr in ms.record:
            rec = Record()
            rec.dt = 1./mr.sampling_rate
            rec.time = mr.time
            rec.data = np.array(mr.data)
            rec_list.append(rec)

    elif format == 'sac':
        sc = sac.Sac(file, byte_order=byte_order)
        rec = Record()
        rec.dt = sc.sampling_rate()
        rec.time = sc.time_date()
        rec.data = np.array(sc.data[0])
        rec_list.append(rec)

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

