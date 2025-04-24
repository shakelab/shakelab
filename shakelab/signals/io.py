# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
from shakelab.signals.libio import sac, dyna

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


def reader(file_path, stream_collection=None, format=None,
           byte_order='be', is_db=False):
    """
    Read seismic waveform data from file(s) into a StreamCollection object.

    Parameters
    ----------
    file_path : str or list
        Path to a waveform file, wildcard pattern, list of files, or a text
        file containing a list of file paths (if `is_db=True`).
    stream_collection : StreamCollection, optional
        Existing StreamCollection to which the data will be appended.
        If None, a new StreamCollection is created.
    format : str, optional
        File format. If None, it is inferred from file extensions.
        Supported: 'mseed', 'sac', 'dyna'.
    byte_order : str, default='be'
        Byte order for reading SAC files ('be' or 'le').
    is_db : bool, default=False
        If True, interpret `file_path` as a text file listing paths to data
        files instead of as a waveform file.

    Returns
    -------
    StreamCollection
        A StreamCollection object containing the loaded waveform records.

    Raises
    ------
    ValueError
        If the specified or inferred format is not recognized.
    NotImplementedError
        If the format is known but not yet implemented.
    """
    if stream_collection is None:
        stream_collection = base.StreamCollection()

    file_list = []

    # Handle list of file paths
    if isinstance(file_path, list):
        file_list = file_path

    # Handle a single string path
    elif isinstance(file_path, str):

        # Explicit text file list mode
        if is_db:
            if not os.path.isfile(file_path):
                raise FileNotFoundError(
                    f"List file not found: {file_path}"
                )
            with open(file_path, 'r') as f:
                file_list = [
                    line.strip() for line in f
                    if line.strip() and not line.startswith('#')
                ]

        # Wildcards
        elif '*' in file_path or '?' in file_path or '[' in file_path:
            file_list = glob(file_path)

        # Directory
        elif os.path.isdir(file_path):
            file_list = glob(os.path.join(file_path, '*'))

        # Single file
        else:
            file_list = [file_path]

    else:
        raise TypeError("file_path must be str or list")

    for file in file_list:
        current_format = format
        if current_format is None:
            fext = os.path.splitext(file)[1].lower()
            if fext in ['.ms', '.mseed', '.miniseed', '.seed']:
                current_format = 'mseed'
            elif fext == '.sac':
                current_format = 'sac'
            elif fext == '.txt':
                current_format = 'dyna'
            else:
                print(f"Unrecognized extension for '{file}'. "
                      f"Defaulting to 'mseed'.")
                current_format = 'mseed'

        # Load the waveform
        if current_format in ['miniseed', 'mseed', 'ms']:
            mseed = get_mseed_module(USE_LIBMSEED)
            stream_collection = mseed.msread(
                file, stream_collection=stream_collection
            )

        elif current_format == 'sac':
            record = sac.sacread(file, byte_order=byte_order)
            stream_collection.append(record)

        elif current_format == 'dyna':
            record = dyna.dynaread(file)
            stream_collection.append(record)

        elif current_format in ['ascii', 'seisan', 'seg2',
                                'dat', 'gse', 'reftek']:
            raise NotImplementedError(f'{current_format}: format '
                                      'not yet implemented')

        else:
            raise ValueError(f'{current_format}: format not recognized')

    return stream_collection


def writer(input_data, file_path, format=None,
            byte_order='be', owrite=False):
    """
    Write one or more seismic records to disk.

    Parameters
    ----------
    input_data : StreamCollection or Stream or Record
        Data to be written. Accepts any shakelab waveform container.
    file_path : str
        Path to the output file (MiniSEED) or directory (SAC).
    format : str, optional
        Output format. If None, inferred from file extension.
        Supported: 'mseed', 'sac'.
    byte_order : str, default='be'
        Byte order for SAC files ('be' or 'le').
    owrite : bool, default=False
        If True, allows overwriting existing files.

    Raises
    ------
    ValueError
        If the format is unknown.
    NotImplementedError
        If the format is not implemented.
    """
    if isinstance(input_data, (base.Record, base.Stream)):
        stream_collection = base.StreamCollection()
        stream_collection.append(input_data)

    elif isinstance(input_data, base.StreamCollection):
        stream_collection = input_data

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
        if not owrite and os.path.exists(file_path):
            print(f"File already exists. Skipping write.")
            return
    
        mseed = get_mseed_module(USE_LIBMSEED)
        mseed.mswrite(file_path, stream_collection)

    elif format == 'sac':
        if not os.path.isdir(file_path):
            os.makedirs(file_path, exist_ok=True)

        for stream in stream_collection:
            for record in stream:
                sid = record.head.sid
                t0 = record.head.time.iso8601
                filename = f"{sid}_{t0}.sac"
                outpath = os.path.join(file_path, filename)
                sac.sacwrite(outpath, record,
                             byte_order=byte_order, owrite=owrite)

    elif format in ['ascii', 'seisan', 'seg2',
                    'dat', 'gse', 'reftek', 'itaca']:
        raise NotImplementedError(f'{format}: format not yet implemented')

    else:
        raise ValueError(f'{format}: format not recognized')
