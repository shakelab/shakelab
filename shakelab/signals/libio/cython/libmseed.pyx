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
MiniSEED Processing Module

This module provides the `MiniSeed` class and associated utility functions 
for processing MiniSEED data using the `libmseed` library. It allows for 
the reading, writing, and manipulation of MiniSEED data, as well as the 
conversion of time formats and encoding/decoding of records.

Note on Encoding:
    0 -> Text encoding (UTF-8)
    1 -> 16-bit integer
    3 -> 32-bit integer
    4 -> 32-bit float (IEEE)
    5 -> 64-bit float (IEEE)
    10 -> Steim-1 compressed integers
    11 -> Steim-2 compressed integers
    12 -> [Legacy] GEOSCOPE 24-bit integer
    13 -> [Legacy] GEOSCOPE 16-bit gain ranged, 3-bit exponent
    14 -> [Legacy] GEOSCOPE 16-bit gain ranged, 4-bit exponent
    16 -> [Legacy] CDSN 16-bit gain ranged
    30 -> [Legacy] SRO 16-bit gain ranged
    32 -> [Legacy] DWWSSN 16-bit gain ranged
"""

cimport libmseed
from cython.operator import dereference

from libc.string cimport memset, memcpy
from libc.stdlib cimport malloc, realloc, free

# -------------------------------------------------------------------------- #
# TIME CONVERSION FUNCTIONS
# -------------------------------------------------------------------------- #

# Define constants for time format and subsecond precision
TIMEFORMAT_ISO = 1 
SUBSECOND_MICRO = 3

def convert_components_to_nstime(int year, int yday, int hour,
                                 int min, int sec, uint32_t nsec):
    """
    Convert time components to nstime_t (nanoseconds since epoch).
    
    Parameters:
        year (int): The year component.
        yday (int): The day of the year component.
        hour (int): The hour component.
        min (int): The minute component.
        sec (int): The second component.
        nsec (uint32_t): The nanosecond component.
    
    Returns:
        nstime_t: The time in nanoseconds since epoch.
    """
    return ms_time2nstime(year, yday, hour, min, sec, nsec)

def convert_timestr_to_nstime(str timestr):
    """
    Convert a formatted time string to nstime_t.
    
    Parameters:
        timestr (str): The ISO formatted time string.
    
    Returns:
        nstime_t: The time in nanoseconds since epoch.
    """
    return ms_timestr2nstime(timestr.encode('utf-8'))

def convert_nstime_to_components(int64_t nstime):
    """
    Convert nstime_t to year, day of year, hour,
    minute, second, and nanosecond components.
    
    Parameters:
        nstime_t (int64_t): The time in nanoseconds since epoch.
    
    Returns:
        tuple: A tuple containing year, day of year, hour,
               minute, second, nanosecond.
    """
    cdef uint16_t year, yday
    cdef uint8_t hour, min, sec
    cdef uint32_t nsec
    
    ms_nstime2time(nstime, &year, &yday, &hour, &min, &sec, &nsec)
    
    return year, yday, hour, min, sec, nsec

def convert_nstime_to_timestr(int64_t nstime):
    """
    Convert nstime_t to a formatted ISO time string.
    
    Parameters:
        nstime_t (int64_t): The time in nanoseconds since epoch.
    
    Returns:
        str: The ISO formatted time string.
    """
    cdef char timestr[30]
    
    cdef int timeformat_iso = TIMEFORMAT_ISO
    cdef int subsecond_micro = SUBSECOND_MICRO
    
    ms_nstime2timestr(nstime, timestr, timeformat_iso, subsecond_micro)
    
    return timestr.decode('utf-8')

# -------------------------------------------------------------------------- #

cdef void record_handler(char *record,
                         int reclen,
                         void *handlerdata) noexcept:
    """
    Handle each MiniSEED record during the packing process.

    This function is called by `mstl3_pack` for each MiniSEED record.
    It appends the record data to a binary buffer (bytearray)
    passed as `handlerdata`.

    Parameters:
        record (char*): Pointer to the MiniSEED record data.
        reclen (int): Length of the MiniSEED record.
        handlerdata (void*): Pointer to the handler data,
                             which is expected to be a bytearray.

    Notes:
        - The `noexcept` specifier indicates that this function
          does not raise Python exceptions.
        - The handler assumes that `handlerdata` is a bytearray
          that can be safely extended.
    """
    # For debug
    # cdef MS3Record* msr = NULL
    # if not msr3_parse (record, reclen, &msr, 0, 0):
    #     msr3_print (msr, 0)
    # msr3_free(&msr)

    # Cast the handlerdata to a bytearray to append the record data
    buffer = <bytearray> handlerdata

    # Ensure the record is valid and has a positive length before appending
    if record != NULL and reclen > 0:
        buffer.extend(<bytes>record[:reclen])

# -------------------------------------------------------------------------- #
# MINISEED MAIN CLASS
# -------------------------------------------------------------------------- #

cdef class MiniSeed:
    """
    MiniSeed class for handling MiniSEED data.
    
    This class provides functionality to initialize, read, write, export, 
    and import MiniSEED records using the libmseed library.
    """
    cdef MS3TraceList* mstl

    def __cinit__(self):
        """
        Initialize the MiniSeed class by creating an MS3TraceList.
        """
        self.mstl = mstl3_init(NULL)

    def __dealloc__(self):
        """
        Free the MS3TraceList when the MiniSeed object is deallocated.
        """
        mstl3_free(&self.mstl, 0)

    @property
    def numtraceids(self):
        return self.mstl.numtraceids

    def print(self):
        pass

    def read(self, bytes byte_stream,
                   uint32_t flags=0,
                   int8_t verbose=0):
        """
        Read and parse MiniSEED data from a byte stream.

        Parameters:
            byte_stream (bytes): The MiniSEED data to be parsed.
            flags (uint32_t): Flags to control parsing behavior.
            verbose (int8_t): Verbosity level for logging.
        """
        cdef MS3Record* msr
        cdef MS3RecordPtr* ptr

        cdef uint64_t bufferlength = len(byte_stream)
        cdef uint64_t offset = 0

        cdef int8_t splitversion = 0
        cdef int parsevalue

        # Set flags for validation, unpacking, and record listing
        flags |= MSF_VALIDATECRC
        flags |= MSF_UNPACKDATA
        flags |= MSF_SKIPNOTDATA
        flags |= MSF_RECORDLIST

        # Parse through the byte stream, creating and adding MS3Records
        while int(bufferlength - offset) > int(MINRECLEN):

            msr = msr3_init(NULL)
            if not msr:
                raise MemoryError("Failed to allocate memory for MS3Record")

            parsevalue = msr3_parse(byte_stream[offset:],
                                    bufferlength,
                                    &msr,
                                    flags,
                                    verbose)

            # Verbose logging of parsed record details
            if verbose > 0:
                iso_starttime = convert_nstime_to_timestr(msr.starttime)
                strout = f"SID: {msr.sid.decode('utf-8')}, "
                strout += f"Start Time: {iso_starttime}, "
                strout += f"Sample Count: {msr.samplecnt}"
                print(strout)

            if parsevalue < 0:
                if msr:
                    msr3_free(&msr)
                raise RuntimeError(f"Error parsing MiniSEED data")

            if parsevalue > 0:
                break

            # Add the parsed MS3Record to the trace list
            mstl3_addmsr_recordptr(self.mstl,
                                   msr,
                                   &ptr,
                                   splitversion,
                                   1,
                                   flags,
                                   NULL)

            offset += msr.reclen

    def write(self, int msformat=2,
                    int reclen=512,
                    int8_t encoding=11,
                    uint32_t flags=0,
                    int8_t verbose=0):
        """
        Write the trace list to a MiniSEED byte stream.

        Parameters:
            msformat (int): MiniSEED format version (2 or 3).
            reclen (int): Record length for the MiniSEED data.
            encoding (int8_t): Encoding format for the data.
            flags (uint32_t): Flags to control writing behavior.
            verbose (int8_t): Verbosity level for logging.

        Returns:
            bytes: The packed MiniSEED data.
        """
        cdef int64_t psamples
        cdef int precords

        byte_stream = bytearray()
    
        # Set flags for validation, unpacking, maintaining trace list,
        # and flushing data
        flags |= MSF_VALIDATECRC
        flags |= MSF_UNPACKDATA
        flags |= MSF_MAINTAINMSTL
        flags |= MSF_FLUSHDATA

        # Set miniSEED version 2 if requested
        if msformat == 2:
            flags |= MSF_PACKVER2

        # Pack the trace list into MiniSEED format
        precords = mstl3_pack(self.mstl,
                              record_handler,
                              <void *> byte_stream,
                              reclen,
                              encoding,
                              &psamples,
                              flags,
                              verbose,
                              NULL)

        return bytes(byte_stream)

    def export_records(self):
        """
        Convert the parsed MiniSEED data into a list of dictionaries.

        Each dictionary contains the NSLC codes, start time, sample rate,
        number of samples, and the data array for a segment.

        Returns:
            list: A list of dictionaries representing the MiniSEED records.
        """
        cdef list record_list = []
        traceid = self.mstl.traces

        # Prepare buffers for NSLC components
        cdef char network[10]
        cdef char station[10]
        cdef char location[10]
        cdef char channel[10]

        cdef int *int_ptr
        cdef float *float_ptr
        cdef double *double_ptr
        cdef char *char_ptr

        while traceid.next[0] != NULL:
            traceid = dereference(traceid.next[0])

            seg = traceid.first
            while seg:

                rng_cnt = range(seg.samplecnt)

                # Parse the SID into FSDN codes
                ms_sid2nslc(traceid.sid, network, station, location, channel)

                # Process the data samples based on the sample type
                if seg.sampletype == 'i':
                    int_ptr = <int *>seg.datasamples
                    data = [int_ptr[i] for i in rng_cnt]

                elif seg.sampletype == 'f':
                    float_ptr = <float *>seg.datasamples
                    data = [float_ptr[i] for i in rng_cnt]

                elif seg.sampletype == 'd':
                    double_ptr = <double *>seg.datasamples
                    data = [double_ptr[i] for i in rng_cnt]

                elif seg.sampletype == 't':
                    char_ptr = <char *>seg.datasamples
                    data = [char_ptr[i].decode('utf-8') for i in rng_cnt]

                else:
                    raise ValueError(f"Unsupported sampletype")

                # Create a record dictionary
                record_dict = {
                    'network': network.decode('utf-8'),
                    'station': station.decode('utf-8'),
                    'location': location.decode('utf-8'),
                    'channel': channel.decode('utf-8'),
                    'starttime': convert_nstime_to_timestr(seg.starttime),
                    'rate': seg.samprate,
                    'nsamp': seg.numsamples,
                    'data' : data
                }
                record_list.append(record_dict)
                seg = seg.next

        return record_list

    def import_record(self, record_dict,
                      int16_t encoding=11,
                      int32_t reclen=512,
                      int8_t verbose=0):
        """
        Import a record from a dictionary and add it to the trace list.

        Parameters:
            record_dict (dict): Dictionary containing record details
                (NSLC codes, start time, sample rate, number of samples).
            encoding (int16_t): Encoding format for the data.
            reclen (int32_t): Record length for the MiniSEED data.
            verbose (int8_t): Verbosity level for logging.
        """
        cdef MS3Record *msr
        cdef MS3RecordPtr* ptr
        cdef char sid[64]
        
        cdef int8_t splitversion = 0
        cdef uint32_t flags=0

        # Set flags for validation, unpacking, and record listing
        flags |= MSF_VALIDATECRC
        flags |= MSF_UNPACKDATA
        flags |= MSF_SKIPNOTDATA
        flags |= MSF_RECORDLIST

        # Initialize a new MiniSEED record
        msr = msr3_init(NULL)
        if not msr:
            raise MemoryError("Failed to allocate memory for MS3Record")

        # Define FDSN codes from the record dictionary
        network = record_dict['network'].encode('utf-8')
        station = record_dict['station'].encode('utf-8')
        location = record_dict['location'].encode('utf-8')
        channel = record_dict['channel'].encode('utf-8')

        cdef char* net = <char*> network  # Network code
        cdef char* sta = <char*> station  # Station code
        cdef char* loc = <char*> location # Location code
        cdef char* chn = <char*> channel  # Channel code

        # Generate SID using the FDSN codes
        ms_nslc2sid(sid, sizeof(sid), 0, net, sta, loc, chn)
        
        # Assign the generated SID to msr.sid
        msr.sid = sid

        # Set the record properties from the dictionary
        msr.starttime = convert_timestr_to_nstime(record_dict['starttime'])
        msr.samprate = <double> record_dict['rate']
        msr.samplecnt = <int64_t> record_dict['nsamp']
        msr.numsamples = msr.samplecnt

        msr.encoding = encoding
        msr.reclen = reclen

        # Allocate memory for data samples
        cdef int32_t* data = <int32_t*> malloc(msr.samplecnt * sizeof(int32_t))
        if not data:
            msr3_free(&msr)
            raise MemoryError("Failed to allocate memory for data samples")

        # Generate example data
        for i in range(msr.samplecnt):
            data[i] = record_dict['data'][i]

        # Set the sample data and its size in the MS3Record structure
        msr.datasamples = <void*> data

        # Determine the sample type and data size based on encoding
        if encoding in [0]:
            msr.sampletype = b't'
            msr.datasize = msr.samplecnt * sizeof(char)

        elif encoding in [1]:
            msr.sampletype = b'i'
            msr.datasize = msr.samplecnt * sizeof(int16_t)

        elif encoding in [4]:
            msr.sampletype = b'f'
            msr.datasize = msr.samplecnt * sizeof(float)

        elif encoding in [5]:
            msr.sampletype = b'd'
            msr.datasize = msr.samplecnt * sizeof(double)

        elif encoding in [3, 10, 11]:
            msr.sampletype = b'i'
            msr.datasize = msr.samplecnt * sizeof(int32_t)

        else:
            raise ValueError(f"Unsupported encoding")

        # NOTE: is it needed to pack data?

        # Add the record to the trace list
        mstl3_addmsr_recordptr(self.mstl,
                               msr,
                               &ptr,
                               splitversion,
                               1,
                               flags,
                               NULL)

        # Free the MS3Record
        msr3_free(&msr)

# -------------------------------------------------------------------------- #

