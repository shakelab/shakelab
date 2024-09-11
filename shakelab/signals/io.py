# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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
from glob import glob

from shakelab.signals import base
from shakelab.signals.libio import sac, smdb

from shakelab.libutils.time import Date

USE_LIBMSEED = True

def get_mseed_module(use_libmseed=USE_LIBMSEED):
    """
    Load the appropriate module conditionally.
    """
    if use_libmseed:
        from shakelab.signals.libio import cymseed as mseed
    else:
        from shakelab.signals.libio import mseed
    return mseed


def reader(file_path, stream_collection=None, format=None, byte_order='be'):
    """
    """
    if os.path.isdir(file_path):
       file_path = os.path.join(file_path, '*')

    file_list = glob(file_path, recursive=False)

    if stream_collection is None:
        stream_collection = base.StreamCollection()

    for file in file_list:

        if format is None:
            # Try to identify file from extension
            fext = os.path.splitext(file)[1]
    
            if fext in ['.ms', '.mseed', '.miniseed', '.seed']:
                format = 'mseed'
    
            elif fext in ['.sac', '.SAC']:
                format = 'sac'
    
            else:
                #raise NotImplementedError('file type not recognized')
                print('File type not recognized. Defualt to mseed')
                format = 'mseed'
    
        # Import recordings
        if format in ['miniseed', 'mseed', 'ms']:
            mseed = get_mseed_module(USE_LIBMSEED)
            stream_collection = mseed.msread(file,
                                     stream_collection=stream_collection)
    
        elif format == 'sac':
            record = sac.sacread(file, byte_order=byte_order)
            stream_collection.append(record)
    
        elif format == 'itaca':
    
            #it = smdb.Itaca(file)
            #record = Record()
            #record.head.rate = it.sampling_rate()
            #record.head.time = Date(it.time)
            #record.data = np.array(it.data)
            #rec_list.append(record)
            pass
    
        elif format == 'ascii':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        elif format == 'seisan':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        elif format == 'seg2':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        elif format == 'dat':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        elif format == 'gse':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        elif format == 'reftek':
            raise NotImplementedError(f'{format}: format not yet implemented')
    
        else:
            raise ValueError(f'{format}: format not recognized')

    return stream_collection


def writer(file_path, stream_collection, format=None, byte_order='be'):
    """
    """
    if format is None:
        # Try to identify file from extension
        fext = os.path.splitext(file_path)[1]

        if fext in ['.ms', '.mseed', '.miniseed', '.seed']:
            format = 'mseed'

        elif fext in ['.sac', '.SAC']:
            format = 'sac'

        else:
            #raise NotImplementedError('file type not recognized')
            print('File type not recognized. Defualt to mseed')
            format = 'mseed'

    if format in ['miniseed', 'mseed', 'ms']:
        mseed = get_mseed_module(USE_LIBMSEED)
        mseed.mswrite(file_path, stream_collection)

    elif format == 'sac':
        for stream in stream_collection:
            for record in stream:
                # TO CHANGE THE FILE NAME
                sac.sacwrite(file, record, byte_order=byte_order)

    elif format == 'itaca':
        pass

    elif format == 'ascii':
        raise NotImplementedError(f'{format}: format not yet implemented')

    elif format == 'seisan':
        raise NotImplementedError(f'{format}: format not yet implemented')

    elif format == 'seg2':
        raise NotImplementedError(f'{format}: format not yet implemented')

    elif format == 'dat':
        raise NotImplementedError(f'{format}: format not yet implemented')

    elif format == 'gse':
        raise NotImplementedError(f'{format}: format not yet implemented')

    elif format == 'reftek':
        raise NotImplementedError(f'{format}: format not yet implemented')

    else:
        raise ValueError(f'{format}: format not recognized')
