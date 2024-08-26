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
An Cython I/O wrapper to the EarthScope MiniSeed library Version 3.
"""
from shakelab.signals.libio.cython import libmseed
from shakelab.signals import base
from shakelab.signals.fdsnws import FDSNCode

import numpy as np
from io import BytesIO

MSEED_VERSION = 2


def msread(input_data_source, stream_collection=None):
    """
    Read MiniSEED data from a file or byte stream and populate a
    StreamCollection.

    Parameters:
        input_data_source (str or bytes):
            The path to the MiniSEED file or a byte stream of MiniSEED data.
        stream_collection (StreamCollection, optional):
            A collection to store the read records. 
            If None, a new StreamCollection is created.

    Returns:
        StreamCollection:
            A collection of records containing the read MiniSEED data.
    """
    if stream_collection is None:
        stream_collection = base.StreamCollection()

    # Determine if the input is a file path or a byte stream
    if isinstance(input_data_source, bytes):
        fid = BytesIO(input_data_source)
    else:
        fid = open(input_data_source, 'rb')

    # Read the MiniSEED data from the file or byte stream
    buffer = fid.read()

    if buffer:
        # Initialize the MiniSeed object and read the buffer data
        ms = libmseed.MiniSeed()
        ms.read(buffer)

        # Export the read records into a list of dictionaries
        record_list = ms.export_records()

        # Process each record dictionary and convert it into a Record object
        if record_list:
            for record_dict in record_list:

                record = base.Record()

                # Format the FDSN code
                fdsn_code = '{0}.{1}.{2}.{3}'.format(
                    record_dict['network'],
                    record_dict['station'],
                    record_dict['location'],
                    record_dict['channel']
                    )

                # Populate the Record object's metadata and data
                record.head.sid = fdsn_code
                record.head.rate = record_dict['rate']
                record.head.time = record_dict['starttime']
                record.data = np.array(record_dict['data'])

                # Append the Record to the StreamCollection
                stream_collection.append(record)

    return stream_collection

def mswrite(output_data_source, stream_collection,
            encoding=11, reclen=512, msformat=None):
    """
    Write a StreamCollection to a MiniSEED file.

    Parameters:
        output_data_source (str):
            The path to the output MiniSEED file.
        stream_collection (StreamCollection):
            A collection of records to be written to the MiniSEED file.
        encoding (int, optional):
            The encoding format to use for the MiniSEED data.
            Default is 11 (Steim-2).
        reclen (int, optional):
            The record length to use for the MiniSEED data.
            Default is 512 bytes.
        msformat (int, optional):
            The miniseed format version of the output file (2 or 3).
            Default is 2.

    Returns:
        None
    """
    if msformat is None:
        msformat = MSEED_VERSION

    # Initialize a new MiniSeed object
    ms = libmseed.MiniSeed()

    # Iterate through each stream/record in the StreamCollection
    for stream in stream_collection:
        code = FDSNCode(stream.sid)

        for record in stream:
            # Create a record dictionary to store the metadata and data
            record_dict = {
                'network' : code.network,
                'station' : code.station,
                'location' : code.location,
                'channel' : code.channel,
                'starttime' : record.time.get_date('iso8601'),
                'rate' : record.rate,
                'nsamp' : record.nsamp,
                'data' : record.data
            }
            # Import the record dictionary into the MiniSeed object
            ms.import_record(record_dict, encoding=encoding, reclen=reclen)

    # Write the MiniSEED data to a byte buffer
    buffer = ms.write(msformat, reclen, encoding)

    # Write the buffer to the specified output file
    with open(output_data_source, 'wb') as fid:
        fid.write(buffer)
