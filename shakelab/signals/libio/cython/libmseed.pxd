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

cdef extern from "stdint.h":

    ctypedef signed char int8_t
    ctypedef unsigned char uint8_t
    ctypedef signed short int int16_t
    ctypedef unsigned short int uint16_t
    ctypedef signed int int32_t
    ctypedef unsigned int uint32_t
    ctypedef signed long long int64_t
    ctypedef unsigned long long uint64_t
    ctypedef int64_t nstime_t

cdef extern from "libmseed.h":

    cdef int MINRECLEN

    cdef int MS_ENDOFFILE        # End of file reached return value
    cdef int MS_NOERROR          # No error
    cdef int MS_GENERROR         # Generic unspecified error
    cdef int MS_NOTSEED          # Data not SEED
    cdef int MS_WRONGLENGTH      # Length of data read was not correct
    cdef int MS_OUTOFRANGE       # SEED record length out of range
    cdef int MS_UNKNOWNFORMAT    # Unknown data encoding format
    cdef int MS_STBADCOMPFLAG    # Steim, invalid compression flag(s)
    cdef int MS_INVALIDCRC       # Invalid CRC

    cdef int MSF_UNPACKDATA      # Unpack data samples
    cdef int MSF_SKIPNOTDATA     # Skip input that cannot be identified as miniSEED
    cdef int MSF_VALIDATECRC     # Validate CRC (if version 3)
    cdef int MSF_PNAMERANGE      # Parse and utilize byte range from path name suffix
    cdef int MSF_ATENDOFFILE     # Reading routine is at the end of the file
    cdef int MSF_SEQUENCE        # UNSUPPORTED: Maintain a record-level sequence number
    cdef int MSF_FLUSHDATA       # Pack all available data even if final record would not be filled
    cdef int MSF_PACKVER2        # Pack as miniSEED version 2 instead of 3
    cdef int MSF_RECORDLIST      # Build a ::MS3RecordList for each ::MS3TraceSeg
    cdef int MSF_MAINTAINMSTL    # Do not modify a trace list when packing

    cdef int DE_TEXT             # Text encoding (UTF-8)
    cdef int DE_INT16            # 16-bit integer
    cdef int DE_INT32            # 32-bit integer
    cdef int DE_FLOAT32          # 32-bit float (IEEE)
    cdef int DE_FLOAT64          # 64-bit float (IEEE)
    cdef int DE_STEIM1           # Steim-1 compressed integers
    cdef int DE_STEIM2           # Steim-2 compressed integers
    cdef int DE_GEOSCOPE24       # [Legacy] GEOSCOPE 24-bit integer
    cdef int DE_GEOSCOPE163      # [Legacy] GEOSCOPE 16-bit gain ranged, 3-bit exponent
    cdef int DE_GEOSCOPE164      # [Legacy] GEOSCOPE 16-bit gain ranged, 4-bit exponent
    cdef int DE_CDSN             # [Legacy] CDSN 16-bit gain ranged
    cdef int DE_SRO              # [Legacy] SRO 16-bit gain ranged
    cdef int DE_DWWSSN           # [Legacy] DWWSSN 16-bit gain ranged

    ctypedef struct MS3Record:
        const char     *record
        int32_t         reclen
        uint8_t         swapflag
        char            sid[64]
        uint8_t         formatversion
        uint8_t         flags
        nstime_t        starttime
        double          samprate
        int16_t         encoding
        uint8_t         pubversion
        int64_t         samplecnt
        uint32_t        crc
        uint16_t        extralength
        uint32_t        datalength
        char           *extra
        void           *datasamples
        uint64_t        datasize
        int64_t         numsamples
        char            sampletype

    ctypedef struct MS3RecordPtr:
        char           *bufferptr
        #FILE          *fileptr
        char           *filename
        int64_t         fileoffset
        MS3Record      *msr
        nstime_t        endtime
        uint32_t        dataoffset
        void           *prvtptr
        MS3RecordPtr   *next

    ctypedef struct MS3RecordList:
        uint64_t        recordcnt
        MS3RecordPtr   *first
        MS3RecordPtr   *last

    ctypedef struct MS3TraceSeg:
        nstime_t        starttime
        nstime_t        endtime
        double          samprate
        int64_t         samplecnt
        void           *datasamples
        uint64_t        datasize
        int64_t         numsamples
        char            sampletype
        void           *prvtptr
        MS3RecordList  *recordlist
        MS3TraceSeg    *prev
        MS3TraceSeg    *next

    ctypedef struct MS3TraceID:
        char            sid[64]
        uint8_t         pubversion
        nstime_t        earliest
        nstime_t        latest
        void           *prvtptr
        uint32_t        numsegments
        MS3TraceSeg    *first
        MS3TraceSeg    *last
        MS3TraceID     *next[8]
        uint8_t         height

    ctypedef struct MS3TraceList:
        uint32_t        numtraceids
        MS3TraceID      traces
        uint64_t        prngstate

    ctypedef struct MS3Tolerance:
        pass

    MS3Record* msr3_init(
        MS3Record *msr
        )

    void msr3_free(
        MS3Record **ppmsr
        )

    MS3TraceList* mstl3_init (
        MS3TraceList *mstl
        )

    void mstl3_free (
        MS3TraceList **ppmstl,
        int8_t freeprvtptr
        )

    int msr3_parse(
        char* buffer,
        uint64_t recbuflen,
        MS3Record **ppmsr,
        uint32_t dataflag,
        int8_t verbose
        )

    #int msr3_unpack_data(
    #    MS3Record *msr,
    #    int8_t verbose
    #    )

    #int msr3_pack(
    #    const MS3Record *msr,
    #    void (*record_handler) (char *, int, void *),
    #    void *handlerdata,
    #    int64_t *packedsamples,
    #    uint32_t flags,
    #    int8_t verbose
    #    )

    void msr3_print(
        const MS3Record *msr,
        int8_t details
        )

    int64_t mstl3_pack(
        MS3TraceList *mstl,
        void (*record_handler) (char *, int, void *),
        void *handlerdata,
        int reclen,
        int8_t encoding,
        int64_t *packedsamples,
        uint32_t flags,
        int8_t verbose,
        char *extra
        )

    MS3TraceSeg* mstl3_addmsr_recordptr(
        MS3TraceList *mstl,
        const MS3Record *msr,
        MS3RecordPtr **recordptr,
        int8_t splitversion,
        int8_t mergegaps,
        uint32_t flags,
        const MS3Tolerance *tolerance
        )

    int ms_nstime2time(
        int64_t nstime,
        uint16_t *year,
        uint16_t *yday,
        uint8_t *hour,
        uint8_t *min,
        uint8_t *sec,
        uint32_t *nsec
        )

    char* ms_nstime2timestr(
        int64_t nstime,
        char *timestr,
        int timeformat,
        int subsecond
        )

    int64_t ms_time2nstime(
        int year,
        int yday,
        int hour,
        int min,
        int sec,
        uint32_t nsec
        )

    int64_t ms_timestr2nstime(
        const char *timestr
        )

    int ms_sid2nslc(
        const char *sid,
        char *net,
        char *sta,
        char *loc,
        char *chan
        )

    int ms_nslc2sid(
        char *sid,
        int sidlen,
        uint16_t flags,
        const char *net,
        const char *sta,
        const char *loc,
        const char *chan
        )
