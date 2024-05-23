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

from shakelab.signals import base
from shakelab.signals.libio import mseed, sac, smdb
from shakelab.libutils.time import Date


def reader(file, ftype=None, stream_collection=None, byte_order='be'):
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
            raise NotImplementedError('file tyep not recognized')

    if stream_collection is None:
        stream_collection = base.StreamCollection()

    # Import recordings
    if ftype == 'mseed':
        stream_collection = mseed.msread(file,
                                         stream_collection=stream_collection,
                                         byte_order=byte_order)

    elif ftype == 'sac':

        record = sac.sacread(file, byte_order=byte_order)
        stream_collection.append(record)

    elif ftype == 'itaca':

        #it = smdb.Itaca(file)
        #record = Record()
        #record.head.rate = it.sampling_rate()
        #record.head.time = Date(it.time)
        #record.data = np.array(it.data)
        #rec_list.append(record)
        pass

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

    return stream_collection

