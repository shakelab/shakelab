# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
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
Lightweight MiniSEED reader and writer.

This module provides a minimal pure-Python MiniSEED implementation used by
ShakeLab.  It is intended mainly as a transparent reader/writer and as a
fallback implementation.  High-performance decoding of compressed MiniSEED
records, especially STEIM1/2, should preferably be delegated to a libmseed
backend by the calling layer.

The implementation supports:

- fixed section data header;
- blockettes 1000 and 1001;
- uncompressed encodings 0, 1, 3, 4;
- STEIM1 and STEIM2 encoding and decoding in pure Python.

The public API includes:

- msread()
- mswrite()
- msrawread()
- msrawiter()
- msrawwrite()
- MSRecord
"""

import logging
from collections.abc import Callable
from copy import deepcopy
from fractions import Fraction

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.signals.binutils import ByteStream


LOGGER = logging.getLogger(__name__)

# MiniSEED fixed-header and payload byte order used by default.
DEFAULT_BYTE_ORDER = "be"

# MiniSEED data encoding identifiers.
ENCODING_ASCII = 0
ENCODING_INT16 = 1
ENCODING_INT32 = 3
ENCODING_FLOAT32 = 4
ENCODING_STEIM1 = 10
ENCODING_STEIM2 = 11

ADMITTED_RECORD_LENGTH = (256, 512, 1024, 2048, 4096, 8192)
ADMITTED_ENCODING = (
    ENCODING_ASCII,
    ENCODING_INT16,
    ENCODING_INT32,
    ENCODING_FLOAT32,
    ENCODING_STEIM1,
    ENCODING_STEIM2,
)

# Fixed section data header size in bytes.
HEADER_SIZE = 48

# One STEIM frame contains sixteen 32-bit words.
FRAME_SIZE = 64

# Signed 32-bit sample limits required by STEIM1 and STEIM2.
INT32_MIN = -(2 ** 31)
INT32_MAX = 2 ** 31 - 1

SIMPLE_DATA_FORMAT = {
    ENCODING_ASCII: ("S1", 1),
    ENCODING_INT16: ("i2", 2),
    ENCODING_INT32: ("i4", 4),
    ENCODING_FLOAT32: ("f4", 4),
}

STRUCT_DATA_FORMAT = {
    ENCODING_ASCII: ("s", 1),
    ENCODING_INT16: ("h", 2),
    ENCODING_INT32: ("i", 4),
    ENCODING_FLOAT32: ("f", 4),
}


def msread(input_data_source, stream_collection=None,
           byte_order=DEFAULT_BYTE_ORDER):
    """
    Read a MiniSEED file and return a ShakeLab StreamCollection.

    Parameters
    ----------
    input_data_source : str or ByteStream
        Input MiniSEED file path or an already opened ByteStream.
    stream_collection : StreamCollection, optional
        Existing collection to which decoded records are appended.
    byte_order : str, optional
        Byte order used by ByteStream. Default is ``"be"``.

    Returns
    -------
    StreamCollection
        Collection containing the decoded ShakeLab records.
    """
    from shakelab.signals.base import StreamCollection

    if stream_collection is None:
        stream_collection = StreamCollection()

    for record in msrawiter(input_data_source, byte_order=byte_order):
        stream_collection.append(record.to_shakelab())

    return stream_collection


def mswrite(output_data_source, stream_collection,
            encoding=ENCODING_STEIM2, reclen=512,
            msformat=None):
    """
    Write a ShakeLab StreamCollection to a MiniSEED 2 file.

    Parameters
    ----------
    output_data_source : str or ByteStream
        Output MiniSEED file path or an already opened ByteStream.
    stream_collection : StreamCollection
        ShakeLab waveform collection to write.
    encoding : int, optional
        MiniSEED encoding. Default is STEIM2.
    reclen : int, optional
        MiniSEED record length in bytes. Default is 512.
    msformat : int, optional
        MiniSEED format version. Only version 2 is supported. ``None``
        is interpreted as version 2.

    Raises
    ------
    TypeError
        If the input collection contains unsupported objects.
    ValueError
        If the encoding, record length, format version or waveform
        metadata are invalid.
    """
    if msformat is None:
        msformat = 2

    if msformat != 2:
        raise ValueError(
            "The pure-Python backend supports only MiniSEED 2"
        )

    _validate_encoding(encoding)
    _validate_record_length(reclen)

    close_stream = not isinstance(output_data_source, ByteStream)

    if close_stream:
        byte_stream = ByteStream(byte_order=DEFAULT_BYTE_ORDER)
        byte_stream.wopen(output_data_source)
    else:
        byte_stream = output_data_source

    sequence_number = 1

    try:
        for stream in stream_collection:
            for record in stream:
                data = np.asarray(record.data)

                if data.ndim != 1:
                    raise ValueError(
                        "MiniSEED waveform data must be one-dimensional"
                    )

                sample_index = 0

                while sample_index < len(data):
                    if sequence_number > 999999:
                        raise ValueError(
                            "MiniSEED sequence number exceeds six digits"
                        )
                
                    raw_record = _from_shakelab_record(
                        record,
                        sample_index,
                        encoding,
                        reclen,
                    )
                
                    raw_record.data = data[sample_index:]
                
                    previous_sample = None
                
                    if (
                        sample_index > 0
                        and encoding in (
                            ENCODING_STEIM1,
                            ENCODING_STEIM2,
                        )
                    ):
                        previous_sample = int(
                            data[sample_index - 1]
                        )
                
                    sample_count = raw_record.write(
                        byte_stream,
                        sequence_number=sequence_number,
                        record_length=reclen,
                        encoding=encoding,
                        allow_partial=True,
                        previous_sample=previous_sample,
                    )
                
                    if sample_count <= 0:
                        raise ValueError(
                            "No waveform samples could be written"
                        )
                
                    sample_index += sample_count
                    sequence_number += 1

    finally:
        if close_stream:
            byte_stream.close()


def msrawread(input_data_source, byte_order=DEFAULT_BYTE_ORDER):
    """
    Read a MiniSEED file and return a list of raw MSRecord objects.

    Parameters
    ----------
    input_data_source : str or ByteStream
        Input MiniSEED file path or an already opened ByteStream.
    byte_order : str, optional
        Byte order used by ByteStream. Default is ``"be"``.

    Returns
    -------
    list of MSRecord
        Raw MiniSEED records.
    """
    return list(msrawiter(input_data_source, byte_order=byte_order))


def msrawiter(input_data_source, byte_order=DEFAULT_BYTE_ORDER):
    """
    Iterate over raw MiniSEED records.

    This generator avoids materialising the full raw record list before
    conversion, reducing memory usage for long files.

    Parameters
    ----------
    input_data_source : str or ByteStream
        Input MiniSEED file path or an already opened ByteStream.
    byte_order : str, optional
        Byte order used by ByteStream. Default is ``"be"``.

    Yields
    ------
    MSRecord
        One decoded raw MiniSEED record at a time.
    """
    close_stream = not isinstance(input_data_source, ByteStream)

    if close_stream:
        byte_stream = ByteStream(byte_order=byte_order)
        byte_stream.ropen(input_data_source)
    else:
        byte_stream = input_data_source

    try:
        while byte_stream.offset < byte_stream.length:
            offset = byte_stream.offset

            try:
                record = MSRecord()
                record.read(byte_stream)
            except EOFError as exc:
                msg = "Truncated MiniSEED record at byte offset {0}"
                raise EOFError(msg.format(offset)) from exc
            except Exception as exc:
                msg = "Invalid MiniSEED record at byte offset {0}"
                raise ValueError(msg.format(offset)) from exc

            yield record

    finally:
        if close_stream:
            byte_stream.close()


def msrawwrite(record_list, output_data_source,
               byte_order=DEFAULT_BYTE_ORDER,
               record_length=None,
               encoding=None):
    """
    Write a list of raw MSRecord objects to MiniSEED.

    Parameters
    ----------
    record_list : sequence of MSRecord
        Records to write.
    output_data_source : str or ByteStream
        Output file path or an already opened ByteStream.
    byte_order : str, optional
        Byte order used by ByteStream. Default is ``"be"``.
    record_length : int, optional
        Output record length in bytes. If omitted, each record retains
        the value declared in blockette 1000.
    encoding : int, optional
        Output MiniSEED encoding. If omitted, each record retains the
        encoding declared in blockette 1000.
    """
    close_stream = not isinstance(output_data_source, ByteStream)

    if close_stream:
        byte_stream = ByteStream(byte_order=byte_order)
        byte_stream.wopen(output_data_source)
    else:
        byte_stream = output_data_source

    try:
        for index, record in enumerate(record_list, start=1):
            record.write(
                byte_stream,
                sequence_number=index,
                record_length=record_length,
                encoding=encoding,
            )
    finally:
        if close_stream:
            byte_stream.close()


def _from_shakelab_record(record, sample_index, encoding,
                          record_length):
    """
    Create an MSRecord from a ShakeLab waveform record.

    Parameters
    ----------
    record : Record
        ShakeLab waveform record.
    sample_index : int
        Index of the first sample represented by the MiniSEED record.
    encoding : int
        MiniSEED encoding.
    record_length : int
        MiniSEED record length in bytes.

    Returns
    -------
    MSRecord
        Initialized raw MiniSEED record.
    """
    if not isinstance(sample_index, (int, np.integer)):
        raise TypeError("sample_index must be an integer")

    if sample_index < 0:
        raise ValueError("sample_index cannot be negative")

    if record.head.delta is None:
        raise ValueError("Waveform sampling interval is not defined")

    delta = float(record.head.delta)

    if not np.isfinite(delta) or delta <= 0:
        raise ValueError(
            "Waveform sampling interval must be positive and finite"
        )

    start_time = record.head.time + sample_index * delta
    rate_factor, rate_multiplier = _sampling_rate_factors(delta)

    network, station, location, channel = _split_fdsn_code(
        record.head.sid
    )

    year, day, hour, minute, second, fraction = _mseed_time_fields(
        start_time
    )

    raw_record = MSRecord()

    raw_record.header.update({
        "SEQUENCE_NUMBER": "000001",
        "DATA_HEADER_QUALITY_INDICATOR": "D",
        "RESERVED_BYTE": " ",
        "STATION_CODE": _format_ascii_field(
            station,
            5,
            "station code",
        ),
        "LOCATION_IDENTIFIER": _format_ascii_field(
            location,
            2,
            "location identifier",
        ),
        "CHANNEL_IDENTIFIER": _format_ascii_field(
            channel,
            3,
            "channel identifier",
        ),
        "NETWORK_CODE": _format_ascii_field(
            network,
            2,
            "network code",
        ),
        "YEAR": year,
        "DAY": day,
        "HOURS": hour,
        "MINUTES": minute,
        "SECONDS": second,
        "UNUSED": 0,
        "MSECONDS": fraction,
        "NUMBER_OF_SAMPLES": 0,
        "SAMPLE_RATE_FACTOR": rate_factor,
        "SAMPLE_RATE_MULTIPLIER": rate_multiplier,
        "ACTIVITY_FLAGS": 0,
        "IO_FLAGS": 0,
        "DATA_QUALITY_FLAGS": 0,
        "NUMBER_OF_BLOCKETTES_TO_FOLLOW": 1,
        "TIME_CORRECTION": 0,
        "OFFSET_TO_BEGINNING_OF_DATA": 0,
        "OFFSET_TO_BEGINNING_OF_BLOCKETTE": HEADER_SIZE,
    })

    raw_record.blockette[1000].update({
        "ENCODING_FORMAT": int(encoding),
        "WORD_ORDER": 1,
        "DATA_RECORD_LENGTH": int(np.log2(record_length)),
        "RESERVED": 0,
    })

    return raw_record


class MSRecord:
    """
    MiniSEED record.

    The class stores fixed-header fields in ``header``, supported blockettes
    in ``blockette`` and decoded data in ``data``.  For numeric encodings,
    ``data`` is a NumPy array.  For encoding 0, ``data`` is a string.
    """

    def __init__(self, byte_stream=None, **kwargs):
        self._record_offset = 0
        self.header = _empty_header()
        self.blockette = _empty_blockettes()
        self.data = np.array([], dtype=np.float64)

        for key, value in kwargs.items():
            if key in self.header:
                self.header[key] = value

        if byte_stream is not None:
            self.read(byte_stream)

    def __len__(self):
        return self.nsamp

    @property
    def nsamp(self):
        """Number of samples declared in the fixed header."""
        return self.header["NUMBER_OF_SAMPLES"]

    @property
    def delta(self):
        """Sampling interval in seconds."""
        srate = self.header["SAMPLE_RATE_FACTOR"]
        rmult = self.header["SAMPLE_RATE_MULTIPLIER"]

        if srate < 0:
            srate = -1.0 / srate

        if rmult < 0:
            rmult = 1.0 / rmult

        srate *= rmult

        return 1.0 / srate

    @property
    def seqn(self):
        """Sequence number."""
        return int(_decode_ascii(self.header["SEQUENCE_NUMBER"]))

    @property
    def time(self):
        """Record start time as a ShakeLab Date object."""
        date = "{0:04d}-".format(self.header["YEAR"])
        date += "{0:03d}T".format(self.header["DAY"])
        date += "{0:02d}:".format(self.header["HOURS"])
        date += "{0:02d}:".format(self.header["MINUTES"])
        date += "{0:02d}.".format(self.header["SECONDS"])
        date += "{0:04d}".format(self.header["MSECONDS"])

        return Date(date)

    @property
    def code(self):
        """FDSN stream identifier."""
        net = _decode_ascii(self.header["NETWORK_CODE"]).strip()
        sta = _decode_ascii(self.header["STATION_CODE"]).strip()
        loc = _decode_ascii(self.header["LOCATION_IDENTIFIER"]).strip()
        chn = _decode_ascii(self.header["CHANNEL_IDENTIFIER"]).strip()

        return "{0}.{1}.{2}.{3}".format(net, sta, loc, chn)

    @property
    def duration(self):
        """Record duration in seconds."""
        return (self.nsamp - 1) * self.delta

    @property
    def record_length(self):
        """MiniSEED record length in bytes."""
        try:
            return 2 ** self.blockette[1000]["DATA_RECORD_LENGTH"]
        except KeyError as exc:
            raise ValueError("Missing required blockette 1000") from exc

    @property
    def data_byte_length(self):
        """Available data payload length in bytes."""
        return self.record_length - self.header["OFFSET_TO_BEGINNING_OF_DATA"]

    def read(self, byte_stream):
        """
        Read the record from a ByteStream.

        Parameters
        ----------
        byte_stream : ByteStream
            Open binary stream positioned at the beginning of a MiniSEED
            record.
        """
        self._record_offset = byte_stream.offset

        self._get_header(byte_stream)
        self._get_blockettes(byte_stream)
        self._validate()
        self._get_data(byte_stream)

        byte_stream.goto(self._record_offset + self.record_length)

    def write(self, byte_stream, sequence_number=None,
              record_length=None, encoding=None,
              allow_partial=False, previous_sample=None):
        """
        Write the record to a ByteStream.

        Parameters
        ----------
        byte_stream : ByteStream
            Open output stream.
        sequence_number : int, optional
            Sequence number to write. If omitted, the existing header value is
            preserved.
        record_length : int, optional
            Output MiniSEED record length in bytes. If omitted, the value
            declared in blockette 1000 is retained.
        encoding : int, optional
            Output MiniSEED encoding. If omitted, the value declared in
            blockette 1000 is retained.
        allow_partial : bool, optional
            If ``True``, write as many samples as fit and return their
            number. If ``False``, raise an error unless all samples fit.
        previous_sample : int, optional
            Sample immediately preceding this record. Used only by STEIM
            encodings to calculate the first difference.

        Returns
        -------
        int
            Number of waveform samples written.

        """
        rec = MSRecord()
        rec.header = deepcopy(self.header)
        rec.blockette = deepcopy(self.blockette)

        rec.header["NUMBER_OF_BLOCKETTES_TO_FOLLOW"] = len(
            rec.blockette
        )

        if sequence_number is not None:
            rec.header["SEQUENCE_NUMBER"] = "{0:06d}".format(
                sequence_number
            )

        if record_length is None:
            record_length = rec.record_length
        else:
            _validate_record_length(record_length)
            rec.blockette[1000]["DATA_RECORD_LENGTH"] = int(
                np.log2(record_length)
            )

        if encoding is None:
            encoding = rec.blockette[1000]["ENCODING_FORMAT"]
        else:
            _validate_encoding(encoding)
            rec.blockette[1000]["ENCODING_FORMAT"] = int(encoding)

        record_offset = byte_stream.offset
        blockette_end = HEADER_SIZE + _total_blockette_size(
            rec.blockette
        )

        rec.header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"] = HEADER_SIZE

        if encoding in (
            ENCODING_STEIM1,
            ENCODING_STEIM2,
        ):
            data_offset = _next_frame_boundary(blockette_end)
            rec.blockette[1000]["WORD_ORDER"] = 1
        else:
            data_offset = blockette_end
            rec.blockette[1000]["WORD_ORDER"] = _mseed_word_order(
                byte_stream
            )

        if data_offset >= record_length:
            raise ValueError(
                "MiniSEED record has no space available for data"
            )

        rec.header["OFFSET_TO_BEGINNING_OF_DATA"] = data_offset

        if encoding in STRUCT_DATA_FORMAT:
            data = _prepare_output_data(self.data, encoding)

            _, sample_size = STRUCT_DATA_FORMAT[encoding]
            max_samples = (
                record_length - data_offset
            ) // sample_size

            sample_count = min(len(data), max_samples)

            if sample_count != len(data) and not allow_partial:
                raise ValueError(
                    "Data do not fit in the selected MiniSEED record "
                    "length: {0} of {1} samples can be written using "
                    "encoding {2}".format(
                        sample_count,
                        len(data),
                        encoding,
                    )
                )

            payload = None
            frame_count = None

        elif encoding in (
            ENCODING_STEIM1,
            ENCODING_STEIM2,
        ):
            (
                data,
                payload,
                sample_count,
                frame_count,
            ) = _encode_steim_payload(
                self.data,
                encoding,
                record_length,
                data_offset,
                require_all=not allow_partial,
                previous_sample=previous_sample,
            )

        else:
            raise ValueError(
                "Unsupported MiniSEED encoding: {0}".format(
                    encoding
                )
            )

        rec.header["NUMBER_OF_SAMPLES"] = sample_count

        if 1001 in rec.blockette:
            if frame_count is None:
                rec.blockette[1001]["FRAME_COUNT"] = 0
            else:
                rec.blockette[1001]["FRAME_COUNT"] = frame_count

        self._write_header(byte_stream, rec)
        self._write_blockettes(
            byte_stream,
            rec,
            record_offset,
        )

        byte_stream.goto(record_offset + data_offset)

        if encoding in STRUCT_DATA_FORMAT:
            self._write_data(
                byte_stream,
                data,
                sample_count,
                encoding,
            )
        else:
            byte_stream.buffer.write(payload)

        padding = record_offset + record_length - byte_stream.offset

        if padding < 0:
            raise ValueError(
                "Encoded data exceed the selected MiniSEED record length"
            )

        if padding:
            byte_stream.buffer.write(b"\x00" * padding)

        return sample_count

    def append(self, record):
        """
        Append another MSRecord to this record.

        Parameters
        ----------
        record : MSRecord
            Record to append.

        Notes
        -----
        This method assumes compatible stream id, sampling rate and timing.
        Consistency checks should be performed by the calling code.
        """
        self.header["SEQUENCE_NUMBER"] = record.header["SEQUENCE_NUMBER"]
        self.header["NUMBER_OF_SAMPLES"] += record.header[
            "NUMBER_OF_SAMPLES"
        ]

        if isinstance(self.data, np.ndarray):
            self.data = np.concatenate((self.data, record.data))
        else:
            self.data += record.data

    def to_shakelab(self):
        """
        Convert the MiniSEED record to a ShakeLab Record.

        Returns
        -------
        Record
            ShakeLab signal record.
        """
        from shakelab.signals.base import Record

        record = Record()

        record.head.sid = self.code
        record.head.delta = self.delta
        record.head.time = self.time

        if isinstance(self.data, np.ndarray):
            record.data = self.data.copy()
        else:
            record.data = np.array(self.data)

        return record

    def _get_header(self, byte_stream):
        """Read the fixed section data header."""
        self.header = {}

        for name, dtype, size in HEAD_STRUCT:
            self.header[name] = byte_stream.get(dtype, size)

    def _get_blockettes(self, byte_stream):
        """Read supported blockettes."""
        self.blockette = {}
        block_offset = self.header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"]

        for _ in range(self.header["NUMBER_OF_BLOCKETTES_TO_FOLLOW"]):
            byte_stream.goto(self._record_offset + block_offset)

            block_type = byte_stream.get("H", 2)
            next_offset = byte_stream.get("H", 2)

            if block_type not in BLOCK_STRUCT:
                LOGGER.debug("Unsupported MiniSEED blockette %s", block_type)
                block_offset = next_offset
                continue

            blockette = {"OFFSET_NEXT": next_offset}

            for name, dtype, size in BLOCK_STRUCT[block_type]:
                blockette[name] = byte_stream.get(dtype, size)

            self.blockette[block_type] = blockette
            block_offset = next_offset

            if block_offset == 0:
                break

    def _get_data(self, byte_stream):
        """Read and decode the record payload."""
        offset = (
            self._record_offset
            + self.header["OFFSET_TO_BEGINNING_OF_DATA"]
        )

        byte_stream.goto(offset)

        nsamp = self.header["NUMBER_OF_SAMPLES"]
        encoding = self.blockette[1000]["ENCODING_FORMAT"]

        if encoding in SIMPLE_DATA_FORMAT:
            self.data = self._read_simple_data(byte_stream, encoding, nsamp)

        elif encoding in (
            ENCODING_STEIM1,
            ENCODING_STEIM2,
        ):
            self.data = self._read_steim_data(byte_stream, encoding, nsamp)

        else:
            raise ValueError("Unsupported MiniSEED encoding: {0}".format(
                encoding
            ))

    def _read_simple_data(self, byte_stream, encoding, nsamp):
        """Read uncompressed data using NumPy from a single byte block."""
        dtype_code, sample_size = SIMPLE_DATA_FORMAT[encoding]
        byte_count = self.data_byte_length

        raw = _read_bytes(byte_stream, byte_count)

        if encoding == ENCODING_ASCII:
            return raw[:nsamp].decode("ascii", errors="replace")

        dtype = np.dtype(dtype_code).newbyteorder(_numpy_byte_order(
            byte_stream
        ))

        available = byte_count // sample_size
        count = min(nsamp, available)

        data = np.frombuffer(raw, dtype=dtype, count=count)

        return data.copy()

    def _read_steim_data(self, byte_stream, encoding, nsamp):
        """Read and decode STEIM1 or STEIM2 compressed data."""
        word_order = self.blockette[1000]["WORD_ORDER"]

        if word_order != 1:
            raise ValueError(
                "STEIM1/2 payload must use big-endian word order"
            )

        raw = _read_bytes(
            byte_stream,
            self.data_byte_length,
        )

        return _decode_steim(raw, encoding, nsamp)

    def _validate(self):
        """Validate minimal record consistency."""
        if 1000 not in self.blockette:
            raise ValueError("Missing required MiniSEED blockette 1000")

        encoding = self.blockette[1000]["ENCODING_FORMAT"]
        record_length = self.record_length

        _validate_encoding(encoding)
        _validate_record_length(record_length)

        data_offset = self.header["OFFSET_TO_BEGINNING_OF_DATA"]

        if data_offset < HEADER_SIZE:
            raise ValueError("Invalid MiniSEED data offset")

        if data_offset > record_length:
            raise ValueError("MiniSEED data offset exceeds record length")

    @staticmethod
    def _write_header(byte_stream, record):
        """Write fixed section data header."""
        for name, dtype, size in HEAD_STRUCT:
            byte_stream.put(record.header[name], dtype, size)

    @staticmethod
    def _write_blockettes(byte_stream, record, record_offset):
        """Write supported blockettes."""
        block_offset = record.header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"]
        block_types = list(record.blockette)

        for index, block_type in enumerate(block_types):
            byte_stream.goto(record_offset + block_offset)

            size = blockette_size(block_type)
            next_offset = block_offset + size

            if index == len(block_types) - 1:
                next_offset = 0

            byte_stream.put(block_type, "H", 2)
            byte_stream.put(next_offset, "H", 2)

            blockette = record.blockette[block_type]

            for name, dtype, field_size in BLOCK_STRUCT[block_type]:
                byte_stream.put(blockette[name], dtype, field_size)

            block_offset += size

    @staticmethod
    def _write_data(byte_stream, data, sample_count, encoding):
        """Write uncompressed data samples."""
        dtype, size = STRUCT_DATA_FORMAT[encoding]

        if encoding == ENCODING_ASCII:
            byte_stream.buffer.write(bytes(data[:sample_count]))
            return

        for value in data[:sample_count]:
            byte_stream.put(value.item(), dtype, size)


def _read_bytes(byte_stream, size):
    """
    Read a raw byte block from ByteStream.

    ByteStream.get("s", n) decodes byte strings to Python strings.
    For MiniSEED binary payloads this is unsafe, so we read directly from
    the underlying buffer.
    """
    raw = byte_stream.buffer.read(size)

    if len(raw) != size:
        raise EOFError("Unexpected end of MiniSEED byte stream")

    return raw


def _prepare_output_data(data, encoding):
    """Convert record data to a writeable sequence."""
    if encoding == ENCODING_ASCII:
        if isinstance(data, str):
            return data.encode("ascii")
        return bytes(data)

    if isinstance(data, np.ndarray):
        return data

    return np.asarray(data)


def _decode_steim(data, encoding, nsamp):
    """
    Decode a STEIM1 or STEIM2 payload.

    Parameters
    ----------
    data : bytes
        Complete STEIM payload composed of 64-byte frames.
    encoding : int
        MiniSEED encoding code for STEIM1 or STEIM2.
    nsamp : int
        Number of samples declared for the record.

    Returns
    -------
    numpy.ndarray
        Decoded signed 32-bit samples.

    Raises
    ------
    ValueError
        If the payload or encoding is invalid, or if the decoded samples
        are inconsistent with the integration constants.
    """
    if encoding not in (
        ENCODING_STEIM1,
        ENCODING_STEIM2,
    ):
        raise ValueError(
            "Unsupported STEIM encoding: {0}".format(encoding)
        )

    if nsamp < 0:
        raise ValueError("STEIM sample count cannot be negative")

    if nsamp == 0:
        return np.array([], dtype=np.int32)

    if len(data) % FRAME_SIZE != 0:
        raise ValueError(
            "STEIM payload length must be a multiple of 64 bytes"
        )

    frame_count = len(data) // FRAME_SIZE
    samples = np.empty(nsamp, dtype=np.int32)

    difference_count = 0
    first = None
    last = None

    for frame_index in range(frame_count):
        frame_offset = frame_index * FRAME_SIZE
        frame = data[frame_offset:frame_offset + FRAME_SIZE]

        words = np.frombuffer(
            frame,
            dtype=">i4",
            count=16,
        )

        if frame_index == 0:
            first = int(words[1])
            last = int(words[2])

        control = int(words[0])
        start_word = 3 if frame_index == 0 else 1

        for word_index in range(start_word, 16):
            code = _binmask(
                control,
                2,
                15 - word_index,
            )

            if code == 0:
                continue

            differences = _w32split(
                int(words[word_index]),
                code,
                encoding,
            )

            for difference in differences:
                if difference_count >= nsamp:
                    break

                samples[difference_count] = difference
                difference_count += 1

            if difference_count >= nsamp:
                break

        if difference_count >= nsamp:
            break

    if first is None or last is None:
        raise ValueError("Empty STEIM data section")

    if difference_count != nsamp:
        raise ValueError(
            "Insufficient differences in STEIM payload: "
            "expected {0}, decoded {1}".format(
                nsamp,
                difference_count,
            )
        )

    # The first encoded difference is d0 = X0 - X(-1). X0 is stored
    # explicitly in the first frame, so d0 is replaced by X0 before
    # integrating the remaining differences.
    samples[0] = first

    for index in range(1, nsamp):
        value = int(samples[index - 1]) + int(samples[index])

        if value < INT32_MIN or value > INT32_MAX:
            raise ValueError(
                "STEIM integration exceeds signed 32-bit range"
            )

        samples[index] = value

    if int(samples[-1]) != last:
        raise ValueError(
            "Sample mismatch in STEIM payload: "
            "expected XN={0}, decoded XN={1}".format(
                last,
                int(samples[-1]),
            )
        )

    return samples


def _prepare_steim_data(data):
    """
    Validate and convert input samples for STEIM encoding.

    Parameters
    ----------
    data : array-like
        Input integer samples.

    Returns
    -------
    numpy.ndarray
        One-dimensional array with dtype int64.

    Raises
    ------
    ValueError
        If the input is not one-dimensional, contains non-integer values,
        or exceeds the signed 32-bit range.
    """
    samples = np.asarray(data)

    if samples.ndim != 1:
        raise ValueError("STEIM data must be one-dimensional")

    if not np.issubdtype(samples.dtype, np.integer):
        raise ValueError("STEIM encoding requires integer samples")

    samples = samples.astype(np.int64, copy=False)

    if samples.size:
        if np.any(samples < INT32_MIN) or np.any(samples > INT32_MAX):
            raise ValueError(
                "STEIM samples must fit in the signed 32-bit range"
            )

    return samples


def _signed_fits(value: int, bits: int) -> bool:
    """
    Return True if a value fits in a signed integer field.

    Parameters
    ----------
    value : int
        Integer value.
    bits : int
        Number of available bits.

    Returns
    -------
    bool
        True if the value can be represented.
    """
    minimum = -(1 << (bits - 1))
    maximum = (1 << (bits - 1)) - 1

    return minimum <= int(value) <= maximum


def _pack_signed_fields(
    values: list[int],
    bits: int,
) -> int:
    """
    Pack signed integer fields into one 32-bit word.

    Values are packed from the most significant field to the least
    significant field.

    Parameters
    ----------
    values : list of int
        Values to pack.
    bits : int
        Number of bits assigned to each value.

    Returns
    -------
    int
        Unsigned 32-bit packed word.
    """
    mask = (1 << bits) - 1
    word = 0

    for value in values:
        if not _signed_fits(value, bits):
            raise ValueError(
                "Value {0} does not fit in {1} signed bits".format(
                    value,
                    bits,
                )
            )

        word = (word << bits) | (int(value) & mask)

    return word & 0xFFFFFFFF


def _set_steim_control(
    control: int,
    word_index: int,
    code: int,
) -> int:
    """
    Set the two-bit control code for one STEIM frame word.

    Parameters
    ----------
    control : int
        Current control word.
    word_index : int
        Frame word index, from 1 to 15.
    code : int
        Two-bit STEIM control code.

    Returns
    -------
    int
        Updated control word.
    """
    if word_index < 1 or word_index > 15:
        raise ValueError("STEIM word index must be between 1 and 15")

    if code < 0 or code > 3:
        raise ValueError("STEIM control code must be between 0 and 3")

    shift = (15 - word_index) * 2

    return control | (code << shift)


def _pack_steim1_word(
    differences: np.ndarray,
    index: int,
) -> tuple[int, int, int]:
    """
    Pack the next STEIM1 difference word.

    Packing priority is:

    - four signed 8-bit differences;
    - two signed 16-bit differences;
    - one signed 32-bit difference.

    Incomplete groups at the end are padded with zeros.

    Parameters
    ----------
    differences : numpy.ndarray
        Array of integer differences.
    index : int
        Index of the first unencoded difference.

    Returns
    -------
    tuple
        Packed word, control code and number of real differences consumed.
    """
    remaining = len(differences) - index

    if remaining <= 0:
        raise ValueError("No STEIM1 differences available for packing")

    count = min(4, remaining)
    values = [
        int(value)
        for value in differences[index:index + count]
    ]

    if all(_signed_fits(value, 8) for value in values):
        packed_values = values + [0] * (4 - count)

        return (
            _pack_signed_fields(packed_values, 8),
            1,
            count,
        )

    count = min(2, remaining)
    values = [
        int(value)
        for value in differences[index:index + count]
    ]

    if all(_signed_fits(value, 16) for value in values):
        packed_values = values + [0] * (2 - count)

        return (
            _pack_signed_fields(packed_values, 16),
            2,
            count,
        )

    value = int(differences[index])

    if not _signed_fits(value, 32):
        raise ValueError(
            "Difference {0} cannot be represented by STEIM1".format(
                value
            )
        )

    return value & 0xFFFFFFFF, 3, 1


def _pack_steim2_word(
    differences: np.ndarray,
    index: int,
) -> tuple[int, int, int]:
    """
    Pack the next STEIM2 difference word.

    Packing priority is:

    - seven signed 4-bit differences;
    - six signed 5-bit differences;
    - five signed 6-bit differences;
    - four signed 8-bit differences;
    - three signed 10-bit differences;
    - two signed 15-bit differences;
    - one signed 30-bit difference.

    Incomplete groups at the end are padded with zeros.

    Parameters
    ----------
    differences : numpy.ndarray
        Array of integer differences.
    index : int
        Index of the first unencoded difference.

    Returns
    -------
    tuple
        Packed word, control code and number of real differences consumed.

    Raises
    ------
    ValueError
        If no differences are available or the next difference cannot be
        represented by STEIM2.
    """
    remaining = len(differences) - index

    if remaining <= 0:
        raise ValueError("No STEIM2 differences available for packing")

    packings = (
        # capacity, bits, control code, DNIB
        (7, 4, 3, 2),
        (6, 5, 3, 1),
        (5, 6, 3, 0),
        (4, 8, 1, None),
        (3, 10, 2, 3),
        (2, 15, 2, 2),
        (1, 30, 2, 1),
    )

    for capacity, bits, control, dnib in packings:
        count = min(capacity, remaining)

        values = [
            int(value)
            for value in differences[index:index + count]
        ]

        if not all(_signed_fits(value, bits) for value in values):
            continue

        packed_values = values + [0] * (capacity - count)
        word = _pack_signed_fields(packed_values, bits)

        if dnib is not None:
            word |= dnib << 30

        return word & 0xFFFFFFFF, control, count

    value = int(differences[index])

    raise ValueError(
        "Difference {0} cannot be represented by STEIM2".format(
            value
        )
    )


SteimWordPacker = Callable[
    [np.ndarray, int],
    tuple[int, int, int],
]


def _encode_steim(
    data,
    frame_count: int,
    pack_word: SteimWordPacker,
    scheme_name: str,
    previous_sample: int | None = None,
    prepared: bool = False,
) -> tuple[bytes, int, int]:
    """
    Encode integer samples as a STEIM payload.

    Parameters
    ----------
    data : array-like
        Input signed integer samples.
    frame_count : int
        Maximum number of 64-byte frames available.
    pack_word : callable
        Function used to pack differences into one STEIM word.
    scheme_name : str
        Human-readable encoding name used in error messages.
    previous_sample : int, optional
        Sample immediately preceding the record. If omitted, zero is used
        as the integration reference for the first difference.
    prepared : bool, optional
        If ``True``, input data are assumed to have already been validated
        and converted by ``_prepare_steim_data()``.

    Returns
    -------
    tuple
        Encoded payload bytes, number of samples encoded and number of
        frames used.

    Raises
    ------
    ValueError
        If the data, frame count or differences are invalid.

    Notes
    -----
    The returned payload contains only the frames actually used. Padding
    to the full MiniSEED record length remains the responsibility of the
    record writer.
    """
    if frame_count < 1:
        raise ValueError("At least one STEIM frame is required")

    if prepared:
        samples = data
    else:
        samples = _prepare_steim_data(data)

    if samples.size == 0:
        return b"", 0, 0

    if previous_sample is None:
        previous_sample = 0

    previous_sample = int(previous_sample)

    if not _signed_fits(previous_sample, 32):
        raise ValueError(
            "Previous STEIM sample must fit in signed 32 bits"
        )

    differences = np.empty(samples.size, dtype=np.int64)

    differences[0] = int(samples[0]) - previous_sample

    if samples.size > 1:
        differences[1:] = samples[1:] - samples[:-1]

    frames = []
    difference_index = 0

    for frame_index in range(frame_count):
        words = [0] * 16
        control = 0

        if frame_index == 0:
            words[1] = int(samples[0]) & 0xFFFFFFFF
            word_index = 3
        else:
            word_index = 1

        while (
            word_index < 16
            and difference_index < len(differences)
        ):
            word, code, consumed = pack_word(
                differences,
                difference_index,
            )

            words[word_index] = word

            control = _set_steim_control(
                control,
                word_index,
                code,
            )

            difference_index += consumed
            word_index += 1

        words[0] = control
        frames.append(words)

        if difference_index >= len(differences):
            break

    sample_count = difference_index

    if sample_count == 0:
        raise ValueError(
            "No samples could be encoded as {0}".format(
                scheme_name
            )
        )

    # XN is the last sample actually stored in this payload.
    frames[0][2] = (
        int(samples[sample_count - 1])
        & 0xFFFFFFFF
    )

    payload = np.asarray(
        frames,
        dtype=">u4",
    ).tobytes()

    return payload, sample_count, len(frames)


def _encode_steim1(
    data,
    frame_count: int,
    previous_sample: int | None = None,
    prepared: bool = False,
) -> tuple[bytes, int, int]:
    """
    Encode integer samples as a STEIM1 payload.

    Parameters
    ----------
    data : array-like
        Input signed integer samples.
    frame_count : int
        Maximum number of 64-byte frames available.
    previous_sample : int, optional
        Sample immediately preceding the record.
    prepared : bool, optional
        If ``True``, input data are assumed to have already been validated
        and converted by ``_prepare_steim_data()``.

    Returns
    -------
    tuple
        Encoded payload bytes, number of samples encoded and number of
        frames used.
    """
    return _encode_steim(
        data=data,
        frame_count=frame_count,
        pack_word=_pack_steim1_word,
        scheme_name="STEIM1",
        previous_sample=previous_sample,
        prepared=prepared,
    )


def _encode_steim2(
    data,
    frame_count: int,
    previous_sample: int | None = None,
    prepared: bool = False,
) -> tuple[bytes, int, int]:
    """
    Encode integer samples as a STEIM2 payload.

    Parameters
    ----------
    data : array-like
        Input signed integer samples.
    frame_count : int
        Maximum number of 64-byte frames available.
    previous_sample : int, optional
        Sample immediately preceding the record.
    prepared : bool, optional
        If ``True``, input data are assumed to have already been validated
        and converted by ``_prepare_steim_data()``.

    Returns
    -------
    tuple
        Encoded payload bytes, number of samples encoded and number of
        frames used.
    """
    return _encode_steim(
        data=data,
        frame_count=frame_count,
        pack_word=_pack_steim2_word,
        scheme_name="STEIM2",
        previous_sample=previous_sample,
        prepared=prepared,
    )


def _encode_steim_payload(
    data,
    encoding: int,
    record_length: int,
    data_offset: int,
    require_all: bool = True,
    previous_sample: int | None = None,
) -> tuple[np.ndarray, bytes, int, int]:
    """
    Encode a complete STEIM payload for one MiniSEED record.

    Parameters
    ----------
    data : array-like
        Input integer samples.
    encoding : int
        MiniSEED encoding code for STEIM1 or STEIM2.
    record_length : int
        Total MiniSEED record length in bytes.
    data_offset : int
        Offset to the beginning of the STEIM payload.
    require_all : bool, optional
        If ``True``, raise an error unless all samples fit in the
        available frames.
    previous_sample : int, optional
        Sample immediately preceding the record.

    Returns
    -------
    tuple
        Prepared integer data, encoded payload bytes, number of samples
        encoded and number of frames used.

    Raises
    ------
    ValueError
        If the encoding is unsupported, no complete frame is available,
        or the samples do not fit in the selected record length.
    """
    available_bytes = record_length - data_offset
    max_frame_count = available_bytes // FRAME_SIZE

    if max_frame_count < 1:
        raise ValueError(
            "MiniSEED record has no complete STEIM frame"
        )

    if encoding == ENCODING_STEIM1:
        encoder = _encode_steim1
        scheme_name = "STEIM1"
        max_differences_per_word = 4

    elif encoding == ENCODING_STEIM2:
        encoder = _encode_steim2
        scheme_name = "STEIM2"
        max_differences_per_word = 7

    else:
        raise ValueError(
            "Unsupported MiniSEED STEIM encoding: {0}".format(
                encoding
            )
        )

    input_data = np.asarray(data)

    if input_data.ndim != 1:
        raise ValueError("STEIM data must be one-dimensional")

    input_count = len(input_data)

    if not require_all:
        data_word_count = (
            13
            + 15 * (max_frame_count - 1)
        )

        max_sample_count = (
            data_word_count
            * max_differences_per_word
        )

        input_data = input_data[:max_sample_count]

    samples = _prepare_steim_data(input_data)

    payload, sample_count, frame_count = encoder(
        samples,
        frame_count=max_frame_count,
        previous_sample=previous_sample,
        prepared=True,
    )

    if sample_count != input_count and require_all:
        raise ValueError(
            "{0} data do not fit in the selected MiniSEED "
            "record length: {1} of {2} samples were "
            "encoded".format(
                scheme_name,
                sample_count,
                input_count,
            )
        )

    return samples, payload, sample_count, frame_count


def _numpy_byte_order(byte_stream):
    """Return NumPy byte-order marker from ByteStream byte order."""
    if byte_stream.byte_order == "le":
        return "<"

    return ">"


def _decode_ascii(value):
    """Decode ASCII fields stored as bytes or return strings unchanged."""
    if isinstance(value, bytes):
        return value.decode("ascii", errors="replace")

    return str(value)


def _empty_header():
    """Return an empty fixed-header dictionary."""
    return {name: None for name, _, _ in HEAD_STRUCT}


def _empty_blockettes():
    """Return default blockette dictionary."""
    blockettes = {1000: {}}

    for name, _, _ in BLOCK_STRUCT[1000]:
        blockettes[1000][name] = None

    return blockettes


def _validate_encoding(encoding):
    """Validate MiniSEED encoding code."""
    if encoding not in ADMITTED_ENCODING:
        raise ValueError("Unsupported MiniSEED encoding: {0}".format(
            encoding
        ))


def _validate_record_length(record_length):
    """Validate MiniSEED record length."""
    if record_length not in ADMITTED_RECORD_LENGTH:
        raise ValueError("Unsupported MiniSEED record length: {0}".format(
            record_length
        ))


def _total_blockette_size(blockettes):
    """Return total encoded size of the supplied blockettes."""
    return sum(blockette_size(block_type) for block_type in blockettes)


def _next_frame_boundary(offset):
    """
    Return the first 64-byte frame boundary at or after an offset.

    Parameters
    ----------
    offset : int
        Byte offset relative to the beginning of the record.

    Returns
    -------
    int
        Offset aligned to a STEIM frame boundary.
    """
    remainder = offset % FRAME_SIZE

    if remainder == 0:
        return offset

    return offset + FRAME_SIZE - remainder


def _binmask(
    word: int,
    bits: int,
    position: int,
) -> int:
    """
    Extract an unsigned bit field from a 32-bit word.

    Parameters
    ----------
    word : int
        Input word.
    bits : int
        Number of bits to extract.
    position : int
        Field position, counted from the least significant side.

    Returns
    -------
    int
        Extracted unsigned integer.
    """
    return (word >> bits * position) & (2 ** bits - 1)


def _getdiff(
    word: int,
    bits: int,
    diff_count: int,
) -> list[int]:
    """
    Split a 32-bit word into signed differences.

    Parameters
    ----------
    word : int
        Input 32-bit word.
    bits : int
        Number of bits per difference.
    diff_count : int
        Number of differences packed in the word.

    Returns
    -------
    list of int
        Signed integer differences.
    """
    output = [0] * diff_count

    for index in range(diff_count):
        value = _binmask(word, bits, index)

        if value >= 2 ** (bits - 1):
            value -= 2 ** bits

        output[(diff_count - 1) - index] = value

    return output


def _w32split(
    word: int,
    order: int,
    scheme: int,
) -> list[int]:
    """
    Split a 32-bit STEIM word into integer differences.

    Parameters
    ----------
    word : int
        Input 32-bit word.
    order : int
        Two-bit STEIM control code.
    scheme : int
        STEIM scheme: 10 for STEIM1, 11 for STEIM2.

    Returns
    -------
    list of int
        Decoded differences.
    """
    if order == 1:
        return _getdiff(word, 8, 4)

    if scheme == ENCODING_STEIM1:
        return _split_steim1(word, order)

    if scheme == ENCODING_STEIM2:
        return _split_steim2(word, order)

    raise ValueError(
        "Unsupported MiniSEED STEIM encoding: {0}".format(
            scheme
        )
    )


def _split_steim1(
    word: int,
    order: int,
) -> list[int]:
    """Decode differences from one STEIM1 data word."""
    if order == 2:
        return _getdiff(word, 16, 2)

    if order == 3:
        return _getdiff(word, 32, 1)

    raise ValueError(
        "Unsupported STEIM1 control code: {0}".format(
            order
        )
    )


def _split_steim2(
    word: int,
    order: int,
) -> list[int]:
    """Decode differences from one STEIM2 data word."""
    dnib = _binmask(word, 2, 15)

    if order == 2:
        if dnib == 1:
            return _getdiff(word, 30, 1)

        if dnib == 2:
            return _getdiff(word, 15, 2)

        if dnib == 3:
            return _getdiff(word, 10, 3)

    if order == 3:
        if dnib == 0:
            return _getdiff(word, 6, 5)

        if dnib == 1:
            return _getdiff(word, 5, 6)

        if dnib == 2:
            return _getdiff(word, 4, 7)

    raise ValueError(
        "Unsupported STEIM2 control code and DNIB combination"
    )


def _mseed_word_order(byte_stream):
    """Return MiniSEED blockette 1000 word-order code."""
    if byte_stream.byte_order == "be":
        return 1

    if byte_stream.byte_order == "le":
        return 0

    raise ValueError("Unsupported byte order: {0}".format(
        byte_stream.byte_order
    ))


def _split_fdsn_code(code):
    """
    Split an FDSN source identifier into MiniSEED header fields.

    Parameters
    ----------
    code : str
        Identifier in ``NET.STA.LOC.CHA`` form.

    Returns
    -------
    tuple of str
        Network, station, location and channel codes.
    """
    if not isinstance(code, str):
        raise TypeError("Waveform source identifier must be a string")

    parts = code.split(".")

    if len(parts) != 4:
        raise ValueError(
            "Waveform source identifier must use NET.STA.LOC.CHA form"
        )

    network, station, location, channel = parts

    if not station:
        raise ValueError("Station code cannot be empty")

    if not channel:
        raise ValueError("Channel code cannot be empty")

    return network, station, location, channel


def _format_ascii_field(value, width, field_name):
    """
    Validate and pad a fixed-width MiniSEED ASCII field.
    """
    value = str(value)

    try:
        value.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(
            "{0} must contain ASCII characters only".format(
                field_name
            )
        ) from exc

    if len(value) > width:
        raise ValueError(
            "{0} exceeds the MiniSEED field width of {1}".format(
                field_name,
                width,
            )
        )

    return value.ljust(width)


def _mseed_time_fields(date):
    """
    Convert a Date object to MiniSEED fixed-header time fields.

    The fractional-second field is expressed in units of 0.0001 seconds.
    """
    second_value = float(date.second)
    second = int(second_value)
    fraction = int(round(
        (second_value - second) * 10000
    ))

    if fraction == 10000:
        date = date + (1.0 - (second_value - second))
        second_value = float(date.second)
        second = int(second_value)
        fraction = 0

    return (
        int(date.year),
        int(date.ordinal_day),
        int(date.hour),
        int(date.minute),
        second,
        fraction,
    )


def _sampling_rate_factors(delta):
    """
    Convert a sampling interval to MiniSEED rate factor and multiplier.

    Parameters
    ----------
    delta : float
        Sampling interval in seconds.

    Returns
    -------
    tuple of int
        MiniSEED sample-rate factor and multiplier.

    Raises
    ------
    ValueError
        If the sampling rate cannot be represented by two signed
        16-bit MiniSEED fields.
    """
    delta = float(delta)

    if not np.isfinite(delta) or delta <= 0:
        raise ValueError(
            "Sampling interval must be positive and finite"
        )

    rate = 1.0 / delta
    fraction = Fraction(str(rate)).limit_denominator(32767)

    numerator = fraction.numerator
    denominator = fraction.denominator

    if numerator > 32767 or denominator > 32767:
        raise ValueError(
            "Sampling rate cannot be represented by MiniSEED "
            "factor and multiplier fields"
        )

    if denominator == 1:
        factor = numerator
        multiplier = 1
    else:
        factor = numerator
        multiplier = -denominator

    reconstructed_rate = factor

    if multiplier < 0:
        reconstructed_rate /= abs(multiplier)
    else:
        reconstructed_rate *= multiplier

    tolerance = max(abs(rate) * 1e-12, 1e-15)

    if abs(reconstructed_rate - rate) > tolerance:
        raise ValueError(
            "Sampling rate cannot be represented accurately by "
            "MiniSEED factor and multiplier fields"
        )

    return int(factor), int(multiplier)


def blockette_size(block_type):
    """
    Return the encoded size of a MiniSEED blockette.

    Parameters
    ----------
    block_type : int
        Blockette type.

    Returns
    -------
    int
        Blockette size in bytes, including the common 4-byte header.
    """
    if block_type not in BLOCK_STRUCT:
        raise ValueError("Unsupported MiniSEED blockette: {0}".format(
            block_type
        ))

    size = sum(item[2] for item in BLOCK_STRUCT[block_type])

    return size + 4


HEAD_STRUCT = [
    ("SEQUENCE_NUMBER", "s", 6),
    ("DATA_HEADER_QUALITY_INDICATOR", "s", 1),
    ("RESERVED_BYTE", "s", 1),
    ("STATION_CODE", "s", 5),
    ("LOCATION_IDENTIFIER", "s", 2),
    ("CHANNEL_IDENTIFIER", "s", 3),
    ("NETWORK_CODE", "s", 2),
    ("YEAR", "H", 2),
    ("DAY", "H", 2),
    ("HOURS", "B", 1),
    ("MINUTES", "B", 1),
    ("SECONDS", "B", 1),
    ("UNUSED", "B", 1),
    ("MSECONDS", "H", 2),
    ("NUMBER_OF_SAMPLES", "H", 2),
    ("SAMPLE_RATE_FACTOR", "h", 2),
    ("SAMPLE_RATE_MULTIPLIER", "h", 2),
    ("ACTIVITY_FLAGS", "B", 1),
    ("IO_FLAGS", "B", 1),
    ("DATA_QUALITY_FLAGS", "B", 1),
    ("NUMBER_OF_BLOCKETTES_TO_FOLLOW", "B", 1),
    ("TIME_CORRECTION", "l", 4),
    ("OFFSET_TO_BEGINNING_OF_DATA", "H", 2),
    ("OFFSET_TO_BEGINNING_OF_BLOCKETTE", "H", 2),
]

BLOCK_STRUCT = {
    1000: [
        ("ENCODING_FORMAT", "B", 1),
        ("WORD_ORDER", "B", 1),
        ("DATA_RECORD_LENGTH", "B", 1),
        ("RESERVED", "B", 1),
    ],
    1001: [
        ("TIMING_QUALITY", "B", 1),
        ("MICRO_SEC", "B", 1),
        ("RESERVED", "B", 1),
        ("FRAME_COUNT", "B", 1),
    ],
}