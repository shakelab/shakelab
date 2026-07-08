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
- STEIM1 and STEIM2 decoding in pure Python.

The public API is intentionally kept compatible with the previous module:

- msread()
- msrawread()
- msrawwrite()
- MSRecord
"""

import logging
from copy import deepcopy

import numpy as np

from shakelab.libutils.time import Date
from shakelab.signals.binutils import ByteStream


LOGGER = logging.getLogger(__name__)

DEFAULT_BYTE_ORDER = "be"

ADMITTED_RECORD_LENGTH = (256, 512, 1024, 2048, 4096, 8192)
ADMITTED_ENCODING = (0, 1, 3, 4, 10, 11)

HEADER_SIZE = 48
FRAME_SIZE = 64

SIMPLE_DATA_FORMAT = {
    0: ("S1", 1),
    1: ("i2", 2),
    3: ("i4", 4),
    4: ("f4", 4),
}

STRUCT_DATA_FORMAT = {
    0: ("s", 1),
    1: ("h", 2),
    3: ("i", 4),
    4: ("f", 4),
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
              record_length=None, encoding=None):
        """
        Write the record to a ByteStream.

        Parameters
        ----------
        byte_stream : ByteStream
            Open output stream.
        sequence_number : int, optional
            Sequence number to write.  If omitted, the existing header value is
            preserved.
        record_length : int, optional
            Output MiniSEED record length in bytes.
        encoding : int, optional
            Output encoding.  Currently only uncompressed encodings are
            supported for writing.
        """
        rec = MSRecord()
        rec.header = deepcopy(self.header)
        rec.blockette = deepcopy(self.blockette)

        if sequence_number is not None:
            rec.header["SEQUENCE_NUMBER"] = "{0:06d}".format(sequence_number)

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
        
        rec.blockette[1000]["WORD_ORDER"] = _mseed_word_order(byte_stream)

        if encoding not in STRUCT_DATA_FORMAT:
            raise NotImplementedError(
                "Writing compressed MiniSEED is not implemented"
            )

        record_offset = byte_stream.offset
        data_offset = HEADER_SIZE + _total_blockette_size(rec.blockette)

        rec.header["OFFSET_TO_BEGINNING_OF_DATA"] = data_offset
        rec.header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"] = HEADER_SIZE

        sample_type, sample_size = STRUCT_DATA_FORMAT[encoding]
        max_samples = (record_length - data_offset) // sample_size

        data = _prepare_output_data(self.data, encoding)
        sample_count = min(len(data), max_samples)

        rec.header["NUMBER_OF_SAMPLES"] = sample_count

        self._write_header(byte_stream, rec)
        self._write_blockettes(byte_stream, rec, record_offset)
        self._write_data(byte_stream, data, sample_count, encoding)

        padding = record_offset + record_length - byte_stream.offset

        for _ in range(padding):
            byte_stream.put(0, "B", 1)

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

        elif encoding in (10, 11):
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

        if encoding == 0:
            return raw[:nsamp].decode("ascii", errors="replace")

        dtype = np.dtype(dtype_code).newbyteorder(_numpy_byte_order(
            byte_stream
        ))

        available = byte_count // sample_size
        count = min(nsamp, available)

        data = np.frombuffer(raw, dtype=dtype, count=count)

        return data.copy()

    def _read_steim_data(self, byte_stream, encoding, nsamp):
        """Decode STEIM1 or STEIM2 compressed data."""
        if nsamp == 0:
            return np.array([], dtype=np.int32)

        if byte_stream.byte_order == "le":
            raise ValueError("STEIM1/2 are only defined as big-endian")

        frame_count = self.data_byte_length // FRAME_SIZE
        data = np.empty(nsamp, dtype=np.int32)

        count = 0
        first = None
        last = None

        for frame_index in range(frame_count):
            raw = _read_bytes(byte_stream, FRAME_SIZE)
            words = np.frombuffer(raw, dtype=">i4", count=16)

            if frame_index == 0:
                first = int(words[1])
                last = int(words[2])

            control = int(words[0])
            codes = [_binmask(control, 2, 15 - idx) for idx in range(16)]

            start_word = 3 if frame_index == 0 else 1

            for word_index in range(start_word, 16):
                code = codes[word_index]

                if code == 0:
                    continue

                diffs = _w32split(int(words[word_index]), code, encoding)

                for diff in diffs:
                    if count >= nsamp:
                        break

                    data[count] = diff
                    count += 1

                if count >= nsamp:
                    break

            if count >= nsamp:
                break

        if first is None or last is None:
            raise ValueError("Empty STEIM data section")

        if nsamp == 0:
            return np.array([], dtype=np.int32)

        data[0] = first

        for idx in range(1, nsamp):
            data[idx] += data[idx - 1]

        if data[-1] != last:
            raise ValueError("Sample mismatch in MiniSEED STEIM record")

        return data

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

        if encoding == 0:
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
    if encoding == 0:
        if isinstance(data, str):
            return data.encode("ascii")
        return bytes(data)

    if isinstance(data, np.ndarray):
        return data

    return np.asarray(data)


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


def _binmask(word, bits, position):
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


def _getdiff(word, bits, diff_count):
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


def _w32split(word, order, scheme):
    """
    Split a 32-bit STEIM word into integer differences.

    Parameters
    ----------
    word : int
        Input 32-bit word.
    order : int
        Control nibble value.
    scheme : int
        STEIM scheme: 10 for STEIM1, 11 for STEIM2.

    Returns
    -------
    list of int
        Decoded differences.
    """
    if order == 1:
        return _getdiff(word, 8, 4)

    if scheme == 10:
        return _split_steim1(word, order)

    if scheme == 11:
        return _split_steim2(word, order)

    raise ValueError("Unsupported STEIM scheme: {0}".format(scheme))


def _split_steim1(word, order):
    """Split one STEIM1 word."""
    if order == 2:
        return _getdiff(word, 16, 2)

    if order == 3:
        return _getdiff(word, 32, 1)

    raise ValueError("Unsupported STEIM1 order: {0}".format(order))


def _split_steim2(word, order):
    """Split one STEIM2 word."""
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

    raise ValueError("Unsupported STEIM2 nibble/order combination")


def _mseed_word_order(byte_stream):
    """Return MiniSEED blockette 1000 word-order code."""
    if byte_stream.byte_order == "be":
        return 1

    if byte_stream.byte_order == "le":
        return 0

    raise ValueError("Unsupported byte order: {0}".format(
        byte_stream.byte_order
    ))


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