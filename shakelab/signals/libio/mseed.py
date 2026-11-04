# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
Pure-Python MiniSEED 2 reader and writer.

This module reads and writes MiniSEED 2 data records (SEED manual v2.4,
chapter 8) without compiled extensions. It is the transparent fallback
of the libmseed backend (:mod:`shakelab.signals.libio.cymseed`); both
return the same ShakeLab objects.

Reading
-------
- Byte order of each record detected from its header (as libmseed
  does), or imposed with ``byte_order``.
- Encodings 0 (text), 1 (int16), 3 (int32), 4 (float32), 5 (float64),
  10 (STEIM1) and 11 (STEIM2), in either word order. Records with
  other encodings are skipped by :func:`msread` with a
  :class:`MSeedWarning`.
- Start time with the microseconds of blockette 1001 and the time
  correction of the fixed header, when not already applied.
- Sampling rate from blockette 100, if present, otherwise from the rate
  factor and multiplier.
- Records without samples or with a null rate (log, detection, timing
  records) are not waveforms: :func:`msread` skips them, the raw
  functions return them.
- Contiguous records of a channel are joined into one
  :class:`~shakelab.signals.base.Record`.

Writing
-------
- Big-endian records with blockette 1000, and blockette 1001 when the
  start time has microseconds beyond the 0.1 ms of the fixed header.
- The encoding follows the type of the data unless given: STEIM2 for
  integers (STEIM1 or int32 if the differences between samples exceed
  the STEIM2 words), 4 for float32, 5 for float64.

MiniSEED 3 is read and written by the libmseed backend only.
"""

from __future__ import annotations

import os
import struct
import warnings
from datetime import date as _date
from fractions import Fraction
from typing import TYPE_CHECKING, Any, Callable, Iterator, Optional

import numpy as np

from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import Record, StreamCollection

__all__ = [
    "ENCODING_ASCII",
    "ENCODING_INT16",
    "ENCODING_INT32",
    "ENCODING_FLOAT32",
    "ENCODING_FLOAT64",
    "ENCODING_STEIM1",
    "ENCODING_STEIM2",
    "MSRecord",
    "MSeedWarning",
    "msread",
    "mswrite",
    "msrawiter",
    "msrawread",
    "msrawwrite",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

# MiniSEED data encoding identifiers (SEED manual, blockette 1000).
ENCODING_ASCII = 0
ENCODING_INT16 = 1
ENCODING_INT32 = 3
ENCODING_FLOAT32 = 4
ENCODING_FLOAT64 = 5
ENCODING_STEIM1 = 10
ENCODING_STEIM2 = 11

STEIM_ENCODINGS = (ENCODING_STEIM1, ENCODING_STEIM2)

# Uncompressed encodings: NumPy type code and sample size in bytes.
SIMPLE_ENCODINGS = {
    ENCODING_ASCII: ("S1", 1),
    ENCODING_INT16: ("i2", 2),
    ENCODING_INT32: ("i4", 4),
    ENCODING_FLOAT32: ("f4", 4),
    ENCODING_FLOAT64: ("f8", 8),
}

SUPPORTED_ENCODINGS = tuple(SIMPLE_ENCODINGS) + STEIM_ENCODINGS

# Record lengths: powers of two from 128 to 65536 bytes.
RECORD_LENGTHS = tuple(2 ** k for k in range(7, 17))

# Fixed section of the data header (48 bytes).
HEADER_SIZE = 48

# One STEIM frame holds sixteen 32-bit words.
FRAME_SIZE = 64

INT32_MIN = -(2 ** 31)
INT32_MAX = 2 ** 31 - 1

# Bit 1 of the activity flags: time correction already applied.
TIME_CORRECTION_APPLIED = 0x02

# Names of the fixed-header fields, in file order, and their layout.
HEADER_FIELDS = (
    "SEQUENCE_NUMBER",
    "DATA_HEADER_QUALITY_INDICATOR",
    "RESERVED_BYTE",
    "STATION_CODE",
    "LOCATION_IDENTIFIER",
    "CHANNEL_IDENTIFIER",
    "NETWORK_CODE",
    "YEAR",
    "DAY",
    "HOURS",
    "MINUTES",
    "SECONDS",
    "UNUSED",
    "MSECONDS",
    "NUMBER_OF_SAMPLES",
    "SAMPLE_RATE_FACTOR",
    "SAMPLE_RATE_MULTIPLIER",
    "ACTIVITY_FLAGS",
    "IO_FLAGS",
    "DATA_QUALITY_FLAGS",
    "NUMBER_OF_BLOCKETTES_TO_FOLLOW",
    "TIME_CORRECTION",
    "OFFSET_TO_BEGINNING_OF_DATA",
    "OFFSET_TO_BEGINNING_OF_BLOCKETTE",
)
HEADER_FORMAT = "6sss5s2s3s2sHHBBBBHHhhBBBBiHH"
TEXT_FIELDS = HEADER_FIELDS[:7]
TEXT_WIDTH = dict(zip(TEXT_FIELDS, (6, 1, 1, 5, 2, 3, 2)))

# Supported blockettes: field names and layout after the 4-byte
# type/next header. Other blockettes are skipped.
BLOCKETTE_FIELDS = {
    100: (("SAMPLE_RATE", "FLAGS", "RESERVED"), "fB3s"),
    1000: (
        ("ENCODING_FORMAT", "WORD_ORDER", "DATA_RECORD_LENGTH",
         "RESERVED"),
        "BBBB",
    ),
    1001: (
        ("TIMING_QUALITY", "MICRO_SEC", "RESERVED", "FRAME_COUNT"),
        "BbBB",
    ),
}

_UNIX_DAY0 = _date(1970, 1, 1).toordinal()
_NS_PER_DAY = 86_400 * 10 ** 9
_NS_PER_TENTH_MS = 100_000


class MSeedWarning(UserWarning):
    """Warning about MiniSEED records that are skipped or suspicious."""


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def msread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
    byte_order: Optional[str] = None,
) -> StreamCollection:
    """
    Read MiniSEED waveforms into a StreamCollection.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        MiniSEED file, its content, or a file opened in binary mode.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.
    byte_order : {"be", "le"}, optional
        Byte order of all records. If None (default), it is detected
        from each record header.

    Returns
    -------
    StreamCollection
        Collection with one record per contiguous segment of each
        channel.

    Warns
    -----
    MSeedWarning
        Once per file, if records with unsupported encodings are
        skipped; for each record whose last decoded sample differs from
        the one declared in the STEIM frames (the data are kept).

    Notes
    -----
    The records of each channel are sorted by time and joined when they
    have the same sampling rate and the gap from the end of the
    previous one is within
    :data:`~shakelab.signals.base.MERGE_TOLERANCE` samples, whatever
    their order in the file.
    Records without samples or sampling rate (log, detection and timing
    records) are skipped.
    """
    from shakelab.signals.base import (
        MERGE_TOLERANCE,
        Record,
        StreamCollection,
    )

    if stream_collection is None:
        stream_collection = StreamCollection()

    pieces: dict = {}
    skipped: dict = {}

    for raw in msrawiter(source, byte_order=byte_order):
        if raw.is_waveform:
            pieces.setdefault(raw.sid, []).append(
                (raw.start_ns, raw.rate, raw.data)
            )
        elif raw.nsamp and raw.rate and raw.data is None:
            skipped[raw.encoding] = skipped.get(raw.encoding, 0) + 1

    if skipped:
        counts = ", ".join(
            f"{count} with encoding {enc}"
            for enc, count in sorted(skipped.items())
        )
        warnings.warn(
            f"MiniSEED records with unsupported encodings skipped: "
            f"{counts}",
            MSeedWarning,
            stacklevel=2,
        )

    segments = []
    for sid, records in pieces.items():
        records.sort(key=lambda item: item[0])
        segment = None
        for start, rate, data in records:
            if segment is not None:
                step = 1e9 / rate
                end = segment["start"] + segment["nsamp"] * step
                if (
                    abs(segment["rate"] - rate) <= 1e-9 * rate
                    and abs(start - end) <= MERGE_TOLERANCE * step
                ):
                    segment["chunks"].append(data)
                    segment["nsamp"] += data.size
                    continue
            segment = {
                "sid": sid,
                "start": start,
                "rate": rate,
                "nsamp": data.size,
                "chunks": [data],
            }
            segments.append(segment)

    for segment in segments:
        chunks = segment["chunks"]
        data = chunks[0] if len(chunks) == 1 else np.concatenate(chunks)
        stream_collection.append(
            Record(
                data=data,
                delta=1.0 / segment["rate"],
                time=Date.from_unix_ns(segment["start"]),
                sid=segment["sid"],
            )
        )

    return stream_collection


def mswrite(
    target: Any,
    stream_collection: Any,
    encoding: Optional[int] = None,
    reclen: int = 512,
    msformat: Optional[int] = None,
) -> None:
    """
    Write waveforms to a MiniSEED 2 file.

    Parameters
    ----------
    target : str, path-like or binary file object
        Output file, or a file opened in binary mode.
    stream_collection : StreamCollection
        Waveforms to write; each record is split into as many MiniSEED
        records as needed.
    encoding : int, optional
        MiniSEED encoding. If None (default), it follows the type of
        each record's data: STEIM2 for integers, 4 for float32 and 5
        for other floating-point types.
    reclen : int, optional
        Record length in bytes, a power of two from 128 to 65536.
        Default is 512.
    msformat : int, optional
        MiniSEED version; only 2 (the default) is supported here.

    Raises
    ------
    ValueError
        If the version, encoding, record length or a source identifier
        is invalid, or if the data cannot be represented with the
        encoding (e.g. floating-point data with STEIM).
    """
    if msformat not in (None, 2):
        raise ValueError(
            "The pure-Python MiniSEED backend writes MiniSEED 2 only"
        )
    if encoding is not None:
        _validate_encoding(encoding)
    _validate_record_length(reclen)

    with _BinaryOutput(target) as fid:
        sequence_number = 1
        for stream in stream_collection:
            for record in stream:
                sequence_number = _write_shakelab_record(
                    fid, record, encoding, reclen, sequence_number
                )


def msrawiter(
    source: Any,
    byte_order: Optional[str] = None,
) -> Iterator[MSRecord]:
    """
    Iterate over the raw records of a MiniSEED file.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        MiniSEED file, its content, or a file opened in binary mode.
    byte_order : {"be", "le"}, optional
        Byte order of all records; detected from each header if None.

    Yields
    ------
    MSRecord
        Decoded records, in file order. Records with unsupported
        encodings have ``data`` set to None.

    Raises
    ------
    ValueError
        If a record is invalid (the byte offset is in the message).

    Warns
    -----
    MSeedWarning
        If fewer bytes than a header remain after the last record and
        they are not all zero (zero padding is ignored silently).
    """
    buffer = memoryview(_read_source(source))
    offset = 0
    size = len(buffer)

    while offset + HEADER_SIZE <= size:
        if not any(buffer[offset:]):
            return
        try:
            record = MSRecord.from_bytes(
                buffer, offset=offset, byte_order=byte_order
            )
        except (ValueError, struct.error) as error:
            raise ValueError(
                f"Invalid MiniSEED record at byte offset {offset}: {error}"
            ) from error
        yield record
        offset += record.record_length

    if offset < size and any(buffer[offset:]):
        warnings.warn(
            f"{size - offset} trailing bytes after the last MiniSEED "
            "record ignored",
            MSeedWarning,
            stacklevel=2,
        )


def msrawread(source: Any, byte_order: Optional[str] = None) -> list:
    """
    Read all raw records of a MiniSEED file (see :func:`msrawiter`).

    Returns
    -------
    list of MSRecord
        Decoded records, in file order.
    """
    return list(msrawiter(source, byte_order=byte_order))


def msrawwrite(
    record_list: Any,
    target: Any,
    byte_order: str = "be",
    record_length: Optional[int] = None,
    encoding: Optional[int] = None,
) -> None:
    """
    Write raw records to a MiniSEED file.

    Parameters
    ----------
    record_list : sequence of MSRecord
        Records to write; sequence numbers are renumbered from 1.
    target : str, path-like or binary file object
        Output file, or a file opened in binary mode.
    byte_order : {"be", "le"}, optional
        Byte order of headers and uncompressed data. STEIM frames are
        always big-endian. Default is ``"be"``.
    record_length : int, optional
        Record length in bytes; the declared one of each record if
        None.
    encoding : int, optional
        Encoding; the declared one of each record if None.

    Raises
    ------
    ValueError
        If a record does not fit in the record length.
    """
    with _BinaryOutput(target) as fid:
        for index, record in enumerate(record_list, start=1):
            content, _ = record.to_bytes(
                sequence_number=index,
                record_length=record_length,
                encoding=encoding,
                byte_order=byte_order,
            )
            fid.write(content)


# ----------------------------------------------------------------------------
# Raw record
# ----------------------------------------------------------------------------


class MSRecord:
    """
    MiniSEED 2 data record.

    Parameters
    ----------
    **fields
        Initial values of fixed-header fields (names as in
        :data:`HEADER_FIELDS`).

    Attributes
    ----------
    header : dict
        Fixed-header fields; text fields are strings.
    blockette : dict
        Supported blockettes (100, 1000, 1001) by type, each a dict of
        fields.
    data : numpy.ndarray, str or None
        Decoded samples (native byte order), text for encoding 0, None
        if the encoding is not supported.
    byte_order : str
        Byte order of the header as read (``"be"`` or ``"le"``).

    Notes
    -----
    Records are built from bytes with :meth:`from_bytes` and converted
    back with :meth:`to_bytes`; unsupported blockettes are not kept.
    """

    def __init__(self, **fields: Any) -> None:
        self.header: dict = {name: None for name in HEADER_FIELDS}
        self.blockette: dict = {
            1000: {
                name: None for name in BLOCKETTE_FIELDS[1000][0]
            }
        }
        self.data: Any = np.array([], dtype=np.int32)
        self.byte_order = "be"

        for key, value in fields.items():
            if key not in self.header:
                raise KeyError(f"unknown MiniSEED header field: {key}")
            self.header[key] = value

    def __len__(self) -> int:
        return self.nsamp

    def __repr__(self) -> str:
        return (
            f"MSRecord({self.sid}, {self.time}, nsamp={self.nsamp}, "
            f"rate={self.rate}, encoding={self.encoding})"
        )

    # -- properties -----------------------------------------------------

    @property
    def nsamp(self) -> int:
        """Number of samples declared in the fixed header."""
        return self.header["NUMBER_OF_SAMPLES"] or 0

    @property
    def rate(self) -> float:
        """
        Sampling rate in hertz (0 for records without sampling).

        Blockette 100 takes precedence over the rate factor and
        multiplier of the fixed header.
        """
        b100 = self.blockette.get(100)
        if b100 and b100["SAMPLE_RATE"] > 0:
            return float(b100["SAMPLE_RATE"])
        return _nominal_rate(
            self.header["SAMPLE_RATE_FACTOR"],
            self.header["SAMPLE_RATE_MULTIPLIER"],
        )

    @property
    def delta(self) -> Optional[float]:
        """Sampling interval in seconds, None without sampling."""
        rate = self.rate
        return 1.0 / rate if rate > 0 else None

    @property
    def seqn(self) -> int:
        """Sequence number."""
        text = self.header["SEQUENCE_NUMBER"].strip()
        return int(text) if text.isdigit() else 0

    @property
    def sid(self) -> str:
        """Source identifier ``NET.STA.LOC.CHA``."""
        h = self.header
        return ".".join(
            h[name].strip() for name in (
                "NETWORK_CODE", "STATION_CODE", "LOCATION_IDENTIFIER",
                "CHANNEL_IDENTIFIER",
            )
        )

    @property
    def start_ns(self) -> int:
        """
        Start time in nanoseconds from the Unix epoch.

        It includes the microseconds of blockette 1001 and the time
        correction of the fixed header, unless the activity flags state
        that it is already applied.
        """
        h = self.header
        days = _date(h["YEAR"], 1, 1).toordinal() - _UNIX_DAY0
        days += h["DAY"] - 1
        ns = days * _NS_PER_DAY
        ns += (h["HOURS"] * 3600 + h["MINUTES"] * 60 + h["SECONDS"]) * 10**9
        ns += h["MSECONDS"] * _NS_PER_TENTH_MS

        b1001 = self.blockette.get(1001)
        if b1001:
            ns += b1001["MICRO_SEC"] * 1000

        correction = h["TIME_CORRECTION"] or 0
        if correction and not (h["ACTIVITY_FLAGS"] or 0) & (
            TIME_CORRECTION_APPLIED
        ):
            ns += correction * _NS_PER_TENTH_MS

        return ns

    @property
    def time(self) -> Date:
        """Start time (see :attr:`start_ns`)."""
        return Date.from_unix_ns(self.start_ns)

    @property
    def duration(self) -> float:
        """Time between the first and the last sample, in seconds."""
        delta = self.delta
        return (self.nsamp - 1) * delta if delta and self.nsamp else 0.0

    @property
    def encoding(self) -> Optional[int]:
        """Encoding declared in blockette 1000."""
        return self.blockette.get(1000, {}).get("ENCODING_FORMAT")

    @property
    def record_length(self) -> int:
        """Record length in bytes, from blockette 1000."""
        exponent = self.blockette.get(1000, {}).get("DATA_RECORD_LENGTH")
        if exponent is None:
            raise ValueError("Missing required blockette 1000")
        return 2 ** exponent

    @property
    def is_waveform(self) -> bool:
        """True for records with decoded numeric samples and a rate."""
        return (
            self.nsamp > 0
            and self.rate > 0
            and isinstance(self.data, np.ndarray)
        )

    # -- conversion -----------------------------------------------------

    @classmethod
    def from_bytes(
        cls,
        buffer: Any,
        offset: int = 0,
        byte_order: Optional[str] = None,
    ) -> MSRecord:
        """
        Decode one record.

        Parameters
        ----------
        buffer : bytes-like
            Content holding the record.
        offset : int, optional
            Byte offset of the record in ``buffer``.
        byte_order : {"be", "le"}, optional
            Byte order of the header; detected if None.

        Returns
        -------
        MSRecord
            Decoded record.

        Raises
        ------
        ValueError
            If the header, the blockettes or the data are invalid.
        """
        buffer = memoryview(buffer)
        record = cls()

        if byte_order is None:
            byte_order = _detect_byte_order(buffer, offset)
        prefix = _struct_prefix(byte_order)
        record.byte_order = byte_order

        values = struct.unpack_from(prefix + HEADER_FORMAT, buffer, offset)
        header = dict(zip(HEADER_FIELDS, values))
        for name in TEXT_FIELDS:
            header[name] = header[name].decode("ascii", errors="replace")
        record.header = header
        _check_header(header)

        record.blockette = _read_blockettes(
            buffer, offset, header, prefix
        )
        if 1000 not in record.blockette:
            raise ValueError("Missing required blockette 1000")

        length = record.record_length
        if length not in RECORD_LENGTHS:
            raise ValueError(f"Unsupported record length: {length}")
        if offset + length > len(buffer):
            raise ValueError("Truncated record")

        record.data = record._decode(buffer[offset:offset + length])
        return record

    def to_bytes(
        self,
        sequence_number: Optional[int] = None,
        record_length: Optional[int] = None,
        encoding: Optional[int] = None,
        allow_partial: bool = False,
        previous_sample: Optional[int] = None,
        byte_order: str = "be",
    ) -> tuple:
        """
        Encode the record.

        Parameters
        ----------
        sequence_number : int, optional
            Sequence number to write; the current one if None.
        record_length : int, optional
            Record length in bytes; the declared one if None.
        encoding : int, optional
            Encoding; the declared one if None.
        allow_partial : bool, optional
            If True, encode as many samples as fit; otherwise raise
            ``ValueError`` unless all fit.
        previous_sample : int, optional
            Sample preceding the record, for the first STEIM difference.
        byte_order : {"be", "le"}, optional
            Byte order of the header and of uncompressed data (STEIM
            frames are big-endian). Default is ``"be"``.

        Returns
        -------
        tuple of (bytes, int)
            Encoded record and number of samples it holds.
        """
        prefix = _struct_prefix(byte_order)
        header = dict(self.header)
        blockettes = {k: dict(v) for k, v in self.blockette.items()}

        if sequence_number is not None:
            if not 0 < sequence_number <= 999999:
                raise ValueError("sequence number exceeds six digits")
            header["SEQUENCE_NUMBER"] = f"{sequence_number:06d}"

        b1000 = blockettes[1000]
        if record_length is None:
            record_length = self.record_length
        _validate_record_length(record_length)
        b1000["DATA_RECORD_LENGTH"] = record_length.bit_length() - 1

        if encoding is None:
            encoding = b1000["ENCODING_FORMAT"]
        _validate_encoding(encoding)
        b1000["ENCODING_FORMAT"] = int(encoding)
        b1000["RESERVED"] = b1000.get("RESERVED") or 0

        steim = encoding in STEIM_ENCODINGS
        b1000["WORD_ORDER"] = 1 if steim or byte_order == "be" else 0

        # Blockettes: 1000 first, then the others in type order.
        order = [1000] + sorted(k for k in blockettes if k != 1000)
        blockette_end = HEADER_SIZE + sum(
            _blockette_size(k) for k in order
        )
        data_offset = (
            _next_frame_boundary(blockette_end) if steim else blockette_end
        )
        if data_offset >= record_length:
            raise ValueError("record has no space available for data")

        if steim:
            _, payload, count, frames = _encode_steim_payload(
                self.data,
                encoding,
                record_length,
                data_offset,
                require_all=not allow_partial,
                previous_sample=previous_sample,
            )
            if 1001 in blockettes:
                blockettes[1001]["FRAME_COUNT"] = frames
        else:
            payload, count = _encode_simple(
                self.data,
                encoding,
                record_length - data_offset,
                prefix,
                allow_partial,
            )
            if 1001 in blockettes:
                blockettes[1001]["FRAME_COUNT"] = 0

        header["NUMBER_OF_SAMPLES"] = count
        header["NUMBER_OF_BLOCKETTES_TO_FOLLOW"] = len(order)
        header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"] = HEADER_SIZE
        header["OFFSET_TO_BEGINNING_OF_DATA"] = data_offset

        content = bytearray(record_length)
        struct.pack_into(
            prefix + HEADER_FORMAT,
            content,
            0,
            *(_header_value(name, header[name]) for name in HEADER_FIELDS),
        )

        position = HEADER_SIZE
        for index, block_type in enumerate(order):
            size = _blockette_size(block_type)
            following = (
                0 if index == len(order) - 1 else position + size
            )
            names, layout = BLOCKETTE_FIELDS[block_type]
            values = [blockettes[block_type].get(name) or 0
                      for name in names]
            if block_type == 100:
                values[2] = b"\x00\x00\x00"
            struct.pack_into(
                prefix + "HH" + layout,
                content,
                position,
                block_type,
                following,
                *values,
            )
            position += size

        content[data_offset:data_offset + len(payload)] = payload
        return bytes(content), count

    def to_shakelab(self) -> Record:
        """
        Convert the record to a ShakeLab Record.

        Raises
        ------
        ValueError
            If the record holds no waveform (see :attr:`is_waveform`).
        """
        from shakelab.signals.base import Record

        if not self.is_waveform:
            raise ValueError("the MiniSEED record holds no waveform")
        return Record(
            data=self.data.copy(),
            delta=self.delta,
            time=self.time,
            sid=self.sid,
        )

    # -- decoding -------------------------------------------------------

    def _decode(self, content: memoryview) -> Any:
        """Decode the data section of the record content."""
        nsamp = self.nsamp
        encoding = self.encoding
        if nsamp == 0:
            return np.array([], dtype=np.int32)

        data_offset = self.header["OFFSET_TO_BEGINNING_OF_DATA"]
        if not HEADER_SIZE <= data_offset < len(content):
            raise ValueError(f"Invalid data offset: {data_offset}")
        payload = content[data_offset:]
        word_order = "be" if self.blockette[1000]["WORD_ORDER"] else "le"

        if encoding in STEIM_ENCODINGS:
            usable = len(payload) - len(payload) % FRAME_SIZE
            samples, xn = _steim_samples(
                payload[:usable], encoding, nsamp, word_order
            )
            if int(samples[-1]) != xn:
                warnings.warn(
                    f"{self.sid} {self.time}: last decoded sample "
                    f"{int(samples[-1])} differs from the declared "
                    f"{xn}; data kept",
                    MSeedWarning,
                    stacklevel=4,
                )
            return samples

        if encoding in SIMPLE_ENCODINGS:
            code, size = SIMPLE_ENCODINGS[encoding]
            if nsamp * size > len(payload):
                raise ValueError("data exceed the record length")
            if encoding == ENCODING_ASCII:
                return bytes(payload[:nsamp]).decode(
                    "ascii", errors="replace"
                )
            dtype = np.dtype(_struct_prefix(word_order) + code)
            data = np.frombuffer(payload, dtype=dtype, count=nsamp)
            return data.astype(dtype.newbyteorder("="))

        return None


# ----------------------------------------------------------------------------
# Reading helpers
# ----------------------------------------------------------------------------


def _read_source(source: Any) -> bytes:
    """Return the content of a path, bytes-like or binary file."""
    if isinstance(source, (bytes, bytearray, memoryview)):
        return source
    if isinstance(source, (str, os.PathLike)):
        with open(source, "rb") as fid:
            return fid.read()
    if hasattr(source, "read"):
        content = source.read()
        if not isinstance(content, (bytes, bytearray)):
            raise TypeError("the file must be opened in binary mode")
        return content
    raise TypeError(
        "source must be a path, bytes or a binary file object"
    )


def _struct_prefix(byte_order: str) -> str:
    """Return the struct/NumPy prefix of a byte order."""
    if byte_order == "be":
        return ">"
    if byte_order == "le":
        return "<"
    raise ValueError(f"byte order must be 'be' or 'le', not {byte_order!r}")


def _detect_byte_order(buffer: memoryview, offset: int) -> str:
    """
    Detect the byte order of a header from its year and day of year
    (the test used by libmseed).
    """
    for byte_order in ("be", "le"):
        year, day = struct.unpack_from(
            _struct_prefix(byte_order) + "HH", buffer, offset + 20
        )
        if 1900 <= year <= 2100 and 1 <= day <= 366:
            return byte_order
    raise ValueError("not a MiniSEED header (invalid year or day)")


def _check_header(header: dict) -> None:
    """Check the plausibility of a fixed header."""
    sequence = header["SEQUENCE_NUMBER"]
    if not all(c.isdigit() or c in " \x00" for c in sequence):
        raise ValueError(f"invalid sequence number {sequence!r}")
    if header["DATA_HEADER_QUALITY_INDICATOR"] not in "DRQM":
        raise ValueError(
            "invalid quality indicator "
            f"{header['DATA_HEADER_QUALITY_INDICATOR']!r}"
        )
    if (
        header["HOURS"] > 23
        or header["MINUTES"] > 59
        or header["SECONDS"] > 60
        or header["MSECONDS"] > 9999
    ):
        raise ValueError("invalid start time")


def _read_blockettes(
    buffer: memoryview,
    offset: int,
    header: dict,
    prefix: str,
) -> dict:
    """Read the supported blockettes of a record."""
    blockettes: dict = {}
    position = header["OFFSET_TO_BEGINNING_OF_BLOCKETTE"]
    size = len(buffer) - offset

    for _ in range(header["NUMBER_OF_BLOCKETTES_TO_FOLLOW"]):
        if position == 0:
            break
        if position < HEADER_SIZE or position + 4 > size:
            raise ValueError(f"invalid blockette offset {position}")
        block_type, following = struct.unpack_from(
            prefix + "HH", buffer, offset + position
        )
        if block_type in BLOCKETTE_FIELDS:
            names, layout = BLOCKETTE_FIELDS[block_type]
            values = struct.unpack_from(
                prefix + layout, buffer, offset + position + 4
            )
            fields = dict(zip(names, values))
            fields["OFFSET_NEXT"] = following
            blockettes[block_type] = fields
        if following and following <= position:
            raise ValueError(f"invalid blockette offset {following}")
        position = following

    return blockettes


def _nominal_rate(factor: Optional[int], multiplier: Optional[int]) -> float:
    """
    Sampling rate from the rate factor and multiplier (SEED manual,
    fixed header field 10).
    """
    if not factor or not multiplier:
        return 0.0
    if factor > 0 and multiplier > 0:
        return float(factor * multiplier)
    if factor > 0:
        return -factor / multiplier
    if multiplier > 0:
        return -multiplier / factor
    return 1.0 / (factor * multiplier)


# STEIM word types: (differences, bits) by control code (and, for
# STEIM2, by the two leading bits of the word, "dnib").
_STEIM1_TYPES = {1: (4, 8), 2: (2, 16), 3: (1, 32)}
_STEIM2_TYPES = {
    (1, 0): (4, 8), (1, 1): (4, 8), (1, 2): (4, 8), (1, 3): (4, 8),
    (2, 1): (1, 30), (2, 2): (2, 15), (2, 3): (3, 10),
    (3, 0): (5, 6), (3, 1): (6, 5), (3, 2): (7, 4),
}


def _steim_tables(types: dict, size: int) -> tuple:
    """Lookup tables of differences and bits per word key."""
    counts = np.zeros(size, dtype=np.int64)
    bits = np.zeros(size, dtype=np.int64)
    for key, (count, nbits) in types.items():
        if isinstance(key, tuple):
            key = key[0] * 4 + key[1]
        counts[key] = count
        bits[key] = nbits
    return counts, bits


_STEIM1_TABLES = _steim_tables(_STEIM1_TYPES, 4)
_STEIM2_TABLES = _steim_tables(_STEIM2_TYPES, 16)


def _steim_samples(
    payload: Any,
    encoding: int,
    nsamp: int,
    word_order: str = "be",
) -> tuple:
    """
    Decode STEIM1 or STEIM2 frames (vectorized).

    Parameters
    ----------
    payload : bytes-like
        Frames, a multiple of 64 bytes.
    encoding : int
        10 (STEIM1) or 11 (STEIM2).
    nsamp : int
        Number of samples declared for the record.
    word_order : {"be", "le"}, optional
        Byte order of the 32-bit words.

    Returns
    -------
    tuple of (numpy.ndarray, int)
        Samples (int32) and the last sample declared in the first frame
        (XN), for the consistency check.

    Raises
    ------
    ValueError
        If the frames hold fewer differences than ``nsamp``, contain an
        invalid word code, or integrate beyond the 32-bit range.
    """
    if encoding not in STEIM_ENCODINGS:
        raise ValueError(f"Unsupported STEIM encoding: {encoding}")
    if len(payload) % FRAME_SIZE:
        raise ValueError("STEIM payload length must be a multiple of 64")
    if len(payload) == 0:
        raise ValueError("Empty STEIM data section")

    words = np.frombuffer(
        payload, dtype=_struct_prefix(word_order) + "u4"
    ).astype(np.int64).reshape(-1, 16)

    x0 = _signed(int(words[0, 1]), 32)
    xn = _signed(int(words[0, 2]), 32)

    shifts = 2 * (15 - np.arange(16))
    codes = (words[:, :1] >> shifts) & 3
    codes[:, 0] = 0
    codes[0, 1:3] = 0

    flat = words.ravel()
    codes = codes.ravel()

    if word_order == "le":
        # As in libmseed: 8-bit differences are bytes in file order,
        # 16-bit ones (STEIM1) two little-endian halfwords in file
        # order; the other words are little-endian 32-bit words.
        flat = flat.copy()
        bytewise = codes == 1
        flat[bytewise] = (
            np.ascontiguousarray(flat[bytewise], dtype="<u4")
            .byteswap()
            .astype(np.int64)
        )
        if encoding == ENCODING_STEIM1:
            halves = codes == 2
            flat[halves] = ((flat[halves] & 0xFFFF) << 16) | (
                flat[halves] >> 16
            )

    if encoding == ENCODING_STEIM1:
        key = codes
        count_table, bits_table = _STEIM1_TABLES
    else:
        key = codes * 4 + (flat >> 30)
        count_table, bits_table = _STEIM2_TABLES

    counts = count_table[key]
    positions = np.cumsum(counts) - counts
    total = int(counts.sum())

    invalid = (codes != 0) & (counts == 0) & (positions < nsamp)
    if np.any(invalid):
        raise ValueError("invalid STEIM word code")
    if total < nsamp:
        raise ValueError(
            f"Insufficient differences in STEIM payload: expected "
            f"{nsamp}, decoded {total}"
        )

    differences = np.empty(total, dtype=np.int64)
    for value in np.unique(key[counts > 0]):
        count = int(count_table[value])
        bits = int(bits_table[value])
        selected = key == value
        word = flat[selected]
        start = positions[selected]
        mask = (1 << bits) - 1
        half = 1 << (bits - 1)
        for j in range(count):
            field = (word >> (bits * (count - 1 - j))) & mask
            differences[start + j] = np.where(
                field >= half, field - (1 << bits), field
            )

    differences = differences[:nsamp]
    samples = x0 + np.cumsum(differences) - differences[0]
    if samples.min() < INT32_MIN or samples.max() > INT32_MAX:
        raise ValueError("STEIM integration exceeds the signed 32-bit range")

    return samples.astype(np.int32), xn


def _decode_steim(
    data: Any,
    encoding: int,
    nsamp: int,
    word_order: str = "be",
) -> np.ndarray:
    """
    Decode a STEIM payload, checking the last sample.

    Raises
    ------
    ValueError
        As :func:`_steim_samples`, and if the last decoded sample
        differs from the declared one.
    """
    if nsamp < 0:
        raise ValueError("STEIM sample count cannot be negative")
    if nsamp == 0:
        return np.array([], dtype=np.int32)
    samples, xn = _steim_samples(data, encoding, nsamp, word_order)
    if int(samples[-1]) != xn:
        raise ValueError(
            f"Sample mismatch in STEIM payload: expected XN={xn}, "
            f"decoded XN={int(samples[-1])}"
        )
    return samples


def _signed(value: int, bits: int) -> int:
    """Interpret an unsigned integer as two's complement."""
    return value - (1 << bits) if value >= 1 << (bits - 1) else value


# ----------------------------------------------------------------------------
# Writing helpers
# ----------------------------------------------------------------------------


class _BinaryOutput:
    """Context manager for a path or an open binary file."""

    def __init__(self, target: Any) -> None:
        self._target = target
        self._file = None

    def __enter__(self) -> Any:
        if isinstance(self._target, (str, os.PathLike)):
            self._file = open(self._target, "wb")
            return self._file
        if hasattr(self._target, "write"):
            return self._target
        raise TypeError("target must be a path or a binary file object")

    def __exit__(self, *args: Any) -> None:
        if self._file is not None:
            self._file.close()


def _write_shakelab_record(
    fid: Any,
    record: Record,
    encoding: Optional[int],
    reclen: int,
    sequence_number: int,
) -> int:
    """
    Write one ShakeLab record as MiniSEED records; return the next
    sequence number.
    """
    data = np.asarray(record.data)
    if data.ndim != 1:
        raise ValueError("MiniSEED waveform data must be one-dimensional")
    if record.head.delta is None:
        raise ValueError("Waveform sampling interval is not defined")
    delta = float(record.head.delta)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError(
            "Waveform sampling interval must be positive and finite"
        )

    if encoding is None:
        encoding = _default_encoding(data)

    network, station, location, channel = _split_fdsn_code(record.head.sid)
    factor, multiplier = _sampling_rate_factors(delta)
    start_ns = record.head.time.to_unix_ns()

    raw = MSRecord(
        DATA_HEADER_QUALITY_INDICATOR="D",
        RESERVED_BYTE=" ",
        STATION_CODE=_format_ascii_field(station, 5, "station code"),
        LOCATION_IDENTIFIER=_format_ascii_field(
            location, 2, "location identifier"
        ),
        CHANNEL_IDENTIFIER=_format_ascii_field(
            channel, 3, "channel identifier"
        ),
        NETWORK_CODE=_format_ascii_field(network, 2, "network code"),
        UNUSED=0,
        SAMPLE_RATE_FACTOR=factor,
        SAMPLE_RATE_MULTIPLIER=multiplier,
        ACTIVITY_FLAGS=0,
        IO_FLAGS=0,
        DATA_QUALITY_FLAGS=0,
        TIME_CORRECTION=0,
    )
    raw.blockette[1000].update(
        ENCODING_FORMAT=int(encoding),
        WORD_ORDER=1,
        DATA_RECORD_LENGTH=reclen.bit_length() - 1,
        RESERVED=0,
    )

    index = 0
    while index < data.size:
        time_fields, micro = _mseed_time_fields(
            start_ns + round(index * delta * 1e9)
        )
        raw.header.update(time_fields)
        if micro:
            raw.blockette[1001] = {
                "TIMING_QUALITY": 0,
                "MICRO_SEC": micro,
                "RESERVED": 0,
                "FRAME_COUNT": 0,
            }
        else:
            raw.blockette.pop(1001, None)

        raw.data = data[index:]
        previous = (
            int(data[index - 1])
            if index > 0 and encoding in STEIM_ENCODINGS
            else None
        )
        content, count = raw.to_bytes(
            sequence_number=sequence_number,
            allow_partial=True,
            previous_sample=previous,
        )
        if count <= 0:
            raise ValueError("No waveform samples could be written")
        fid.write(content)
        index += count
        sequence_number += 1

    return sequence_number


def _default_encoding(data: np.ndarray) -> int:
    """
    Encoding following the type of the data.

    Integers use STEIM2 when the differences between samples fit its
    30-bit words, otherwise STEIM1 (32-bit differences) or, beyond that,
    uncompressed int32.
    """
    if np.issubdtype(data.dtype, np.integer):
        values = data.astype(np.int64)
        if values.size and (
            values.min() < INT32_MIN or values.max() > INT32_MAX
        ):
            raise ValueError("integer samples exceed the 32-bit range")
        jump = int(np.abs(np.diff(values)).max()) if values.size > 1 else 0
        if jump < 2**29:
            return ENCODING_STEIM2
        if jump < 2**31:
            return ENCODING_STEIM1
        return ENCODING_INT32
    if data.dtype == np.float32:
        return ENCODING_FLOAT32
    if np.issubdtype(data.dtype, np.floating):
        return ENCODING_FLOAT64
    raise ValueError(
        f"data of type {data.dtype} cannot be written as MiniSEED"
    )


def _encode_simple(
    data: Any,
    encoding: int,
    available: int,
    prefix: str,
    allow_partial: bool,
) -> tuple:
    """
    Encode uncompressed samples; return the bytes and their number.
    """
    code, size = SIMPLE_ENCODINGS[encoding]
    capacity = available // size

    if encoding == ENCODING_ASCII:
        text = data if isinstance(data, (bytes, bytearray)) else (
            str(data).encode("ascii")
        )
        count = min(len(text), capacity)
        if count < len(text) and not allow_partial:
            raise ValueError(
                f"Data do not fit in the record length: {count} of "
                f"{len(text)} characters"
            )
        return bytes(text[:count]), count

    values = np.asarray(data)
    count = min(values.size, capacity)
    if count < values.size and not allow_partial:
        raise ValueError(
            f"Data do not fit in the selected MiniSEED record length: "
            f"{count} of {values.size} samples can be written using "
            f"encoding {encoding}"
        )
    values = values[:count]
    dtype = np.dtype(prefix + code)

    if np.issubdtype(dtype, np.integer):
        if not np.issubdtype(values.dtype, np.integer):
            raise ValueError(
                f"encoding {encoding} requires integer samples"
            )
        info = np.iinfo(dtype)
        if values.size and (
            values.min() < info.min or values.max() > info.max
        ):
            raise ValueError(
                f"samples exceed the range of encoding {encoding}"
            )

    return values.astype(dtype).tobytes(), count


def _header_value(name: str, value: Any) -> Any:
    """Return a fixed-header value ready for packing."""
    if name in TEXT_WIDTH:
        text = "" if value is None else str(value)
        return text.ljust(TEXT_WIDTH[name]).encode("ascii")
    return 0 if value is None else int(value)


def _blockette_size(block_type: int) -> int:
    """Size of a blockette in bytes, its 4-byte header included."""
    if block_type not in BLOCKETTE_FIELDS:
        raise ValueError(f"Unsupported MiniSEED blockette: {block_type}")
    return 4 + struct.calcsize(">" + BLOCKETTE_FIELDS[block_type][1])


def _next_frame_boundary(offset: int) -> int:
    """First 64-byte frame boundary at or after an offset."""
    return -(-offset // FRAME_SIZE) * FRAME_SIZE


def _validate_encoding(encoding: Any) -> None:
    """Check that an encoding is supported."""
    if encoding not in SUPPORTED_ENCODINGS:
        raise ValueError(f"Unsupported MiniSEED encoding: {encoding}")


def _validate_record_length(record_length: Any) -> None:
    """Check that a record length is supported."""
    if record_length not in RECORD_LENGTHS:
        raise ValueError(
            f"Unsupported MiniSEED record length: {record_length}"
        )


def _split_fdsn_code(code: Any) -> tuple:
    """
    Split a source identifier ``NET.STA.LOC.CHA`` into its codes.

    Raises
    ------
    TypeError
        If the identifier is not a string.
    ValueError
        If it does not have four parts, or the station or channel code
        is empty.
    """
    if not isinstance(code, str):
        raise TypeError("Waveform source identifier must be a string")
    parts = code.split(".")
    if len(parts) != 4:
        raise ValueError(
            "Waveform source identifier must use NET.STA.LOC.CHA form"
        )
    if not parts[1]:
        raise ValueError("Station code cannot be empty")
    if not parts[3]:
        raise ValueError("Channel code cannot be empty")
    return tuple(parts)


def _format_ascii_field(value: Any, width: int, field_name: str) -> str:
    """Validate and pad a fixed-width ASCII field."""
    value = str(value)
    if not value.isascii():
        raise ValueError(f"{field_name} must contain ASCII characters only")
    if len(value) > width:
        raise ValueError(
            f"{field_name} exceeds the MiniSEED field width of {width}"
        )
    return value.ljust(width)


def _mseed_time_fields(unix_ns: int) -> tuple:
    """
    Fixed-header time fields of an instant, and the microseconds for
    blockette 1001.

    The instant is rounded to the microsecond; the fixed header holds
    it to 0.1 ms and the remainder, from -50 to 49 microseconds, goes
    to blockette 1001 (as libmseed does).

    Returns
    -------
    tuple of (dict, int)
        Time fields of the fixed header and the microsecond offset.
    """
    micro_total = (unix_ns + 500) // 1000
    tenths = (micro_total + 50) // 100
    micro = micro_total - tenths * 100

    days, rest = divmod(tenths, 864_000_000)
    day = _date.fromordinal(days + _UNIX_DAY0)
    seconds, fraction = divmod(rest, 10_000)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)

    fields = {
        "YEAR": day.year,
        "DAY": day.timetuple().tm_yday,
        "HOURS": hours,
        "MINUTES": minutes,
        "SECONDS": seconds,
        "MSECONDS": fraction,
    }
    return fields, micro


def _sampling_rate_factors(delta: float) -> tuple:
    """
    Rate factor and multiplier of a sampling interval.

    Raises
    ------
    ValueError
        If the rate cannot be represented by two signed 16-bit fields.
    """
    delta = float(delta)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError("Sampling interval must be positive and finite")

    rate = 1.0 / delta
    fraction = Fraction(str(rate)).limit_denominator(32767)
    numerator, denominator = fraction.numerator, fraction.denominator

    if numerator > 32767 or denominator > 32767:
        raise ValueError(
            "Sampling rate cannot be represented by MiniSEED factor and "
            "multiplier fields"
        )

    multiplier = 1 if denominator == 1 else -denominator
    if abs(_nominal_rate(numerator, multiplier) - rate) > max(
        rate * 1e-12, 1e-15
    ):
        raise ValueError(
            "Sampling rate cannot be represented accurately by MiniSEED "
            "factor and multiplier fields"
        )
    return int(numerator), int(multiplier)


# -- STEIM encoding -----------------------------------------------------------


def _signed_fits(value: int, bits: int) -> bool:
    """True if a value fits in a signed field of ``bits`` bits."""
    return -(1 << (bits - 1)) <= int(value) <= (1 << (bits - 1)) - 1


def _pack_signed_fields(values: list, bits: int) -> int:
    """Pack signed fields into one word, most significant first."""
    mask = (1 << bits) - 1
    word = 0
    for value in values:
        if not _signed_fits(value, bits):
            raise ValueError(f"Value {value} does not fit in {bits} bits")
        word = (word << bits) | (int(value) & mask)
    return word & 0xFFFFFFFF


def _pack_steim1_word(differences: np.ndarray, index: int) -> tuple:
    """
    Pack the next STEIM1 word: four 8-bit, two 16-bit or one 32-bit
    differences. Return the word, its control code and the number of
    differences used.
    """
    remaining = len(differences) - index
    if remaining <= 0:
        raise ValueError("No STEIM1 differences available for packing")

    for capacity, bits, code in ((4, 8, 1), (2, 16, 2)):
        count = min(capacity, remaining)
        values = [int(v) for v in differences[index:index + count]]
        if all(_signed_fits(v, bits) for v in values):
            word = _pack_signed_fields(
                values + [0] * (capacity - count), bits
            )
            return word, code, count

    value = int(differences[index])
    if not _signed_fits(value, 32):
        raise ValueError(
            f"Difference {value} cannot be represented by STEIM1"
        )
    return value & 0xFFFFFFFF, 3, 1


def _pack_steim2_word(differences: np.ndarray, index: int) -> tuple:
    """
    Pack the next STEIM2 word (from seven 4-bit to one 30-bit
    differences). Return the word, its control code and the number of
    differences used.
    """
    remaining = len(differences) - index
    if remaining <= 0:
        raise ValueError("No STEIM2 differences available for packing")

    packings = (
        # capacity, bits, control code, dnib
        (7, 4, 3, 2),
        (6, 5, 3, 1),
        (5, 6, 3, 0),
        (4, 8, 1, None),
        (3, 10, 2, 3),
        (2, 15, 2, 2),
        (1, 30, 2, 1),
    )
    for capacity, bits, code, dnib in packings:
        count = min(capacity, remaining)
        values = [int(v) for v in differences[index:index + count]]
        if not all(_signed_fits(v, bits) for v in values):
            continue
        word = _pack_signed_fields(values + [0] * (capacity - count), bits)
        if dnib is not None:
            word |= dnib << 30
        return word & 0xFFFFFFFF, code, count

    raise ValueError(
        f"Difference {int(differences[index])} cannot be represented by "
        "STEIM2"
    )


def _prepare_steim_data(data: Any) -> np.ndarray:
    """
    Check and convert samples for STEIM encoding (int64 copy).

    Raises
    ------
    ValueError
        If the samples are not one-dimensional integers within the
        signed 32-bit range.
    """
    samples = np.asarray(data)
    if samples.ndim != 1:
        raise ValueError("STEIM data must be one-dimensional")
    if not np.issubdtype(samples.dtype, np.integer):
        raise ValueError(
            "STEIM encoding requires integer samples (use encoding 4 or "
            "5, or None, for floating-point data)"
        )
    samples = samples.astype(np.int64, copy=False)
    if samples.size and (
        samples.min() < INT32_MIN or samples.max() > INT32_MAX
    ):
        raise ValueError("STEIM samples must fit in the signed 32-bit range")
    return samples


def _encode_steim(
    data: Any,
    frame_count: int,
    pack_word: Callable,
    scheme_name: str,
    previous_sample: Optional[int] = None,
    prepared: bool = False,
) -> tuple:
    """
    Encode samples as STEIM frames.

    Parameters
    ----------
    data : array-like
        Signed integer samples.
    frame_count : int
        Maximum number of 64-byte frames.
    pack_word : callable
        Word packer (:func:`_pack_steim1_word` or
        :func:`_pack_steim2_word`).
    scheme_name : str
        Name used in error messages.
    previous_sample : int, optional
        Sample preceding the record (0 if None), for the first
        difference.
    prepared : bool, optional
        True if ``data`` already went through
        :func:`_prepare_steim_data`.

    Returns
    -------
    tuple of (bytes, int, int)
        Frames actually used, number of samples encoded and number of
        frames.
    """
    if frame_count < 1:
        raise ValueError("At least one STEIM frame is required")

    samples = data if prepared else _prepare_steim_data(data)
    if samples.size == 0:
        return b"", 0, 0

    previous_sample = int(previous_sample or 0)
    if not _signed_fits(previous_sample, 32):
        raise ValueError("Previous STEIM sample must fit in signed 32 bits")

    differences = np.empty(samples.size, dtype=np.int64)
    differences[0] = int(samples[0]) - previous_sample
    differences[1:] = np.diff(samples)

    frames = []
    used = 0
    for frame_index in range(frame_count):
        words = [0] * 16
        control = 0
        word_index = 1
        if frame_index == 0:
            words[1] = int(samples[0]) & 0xFFFFFFFF
            word_index = 3

        while word_index < 16 and used < differences.size:
            word, code, consumed = pack_word(differences, used)
            words[word_index] = word
            control |= code << (2 * (15 - word_index))
            used += consumed
            word_index += 1

        words[0] = control
        frames.append(words)
        if used >= differences.size:
            break

    if used == 0:
        raise ValueError(f"No samples could be encoded as {scheme_name}")

    # XN is the last sample actually stored.
    frames[0][2] = int(samples[used - 1]) & 0xFFFFFFFF
    payload = np.asarray(frames, dtype=">u4").tobytes()
    return payload, used, len(frames)


def _encode_steim1(
    data: Any,
    frame_count: int,
    previous_sample: Optional[int] = None,
    prepared: bool = False,
) -> tuple:
    """Encode samples as STEIM1 frames (see :func:`_encode_steim`)."""
    return _encode_steim(
        data, frame_count, _pack_steim1_word, "STEIM1",
        previous_sample, prepared,
    )


def _encode_steim2(
    data: Any,
    frame_count: int,
    previous_sample: Optional[int] = None,
    prepared: bool = False,
) -> tuple:
    """Encode samples as STEIM2 frames (see :func:`_encode_steim`)."""
    return _encode_steim(
        data, frame_count, _pack_steim2_word, "STEIM2",
        previous_sample, prepared,
    )


def _encode_steim_payload(
    data: Any,
    encoding: int,
    record_length: int,
    data_offset: int,
    require_all: bool = True,
    previous_sample: Optional[int] = None,
) -> tuple:
    """
    Encode the STEIM payload of one record.

    Returns
    -------
    tuple of (numpy.ndarray, bytes, int, int)
        Prepared samples, frames, number of samples encoded and number
        of frames.

    Raises
    ------
    ValueError
        If the encoding is not STEIM, no frame fits, or (with
        ``require_all``) not all samples fit.
    """
    max_frames = (record_length - data_offset) // FRAME_SIZE
    if max_frames < 1:
        raise ValueError("MiniSEED record has no complete STEIM frame")

    if encoding == ENCODING_STEIM1:
        encoder, name, per_word = _encode_steim1, "STEIM1", 4
    elif encoding == ENCODING_STEIM2:
        encoder, name, per_word = _encode_steim2, "STEIM2", 7
    else:
        raise ValueError(f"Unsupported MiniSEED STEIM encoding: {encoding}")

    values = np.asarray(data)
    if values.ndim != 1:
        raise ValueError("STEIM data must be one-dimensional")
    input_count = values.size

    if not require_all:
        # Upper bound of the samples that fit: only these are prepared.
        values = values[:(13 + 15 * (max_frames - 1)) * per_word]

    samples = _prepare_steim_data(values)
    payload, count, frames = encoder(
        samples,
        frame_count=max_frames,
        previous_sample=previous_sample,
        prepared=True,
    )
    if require_all and count != input_count:
        raise ValueError(
            f"{name} data do not fit in the selected MiniSEED record "
            f"length: {count} of {input_count} samples were encoded"
        )
    return samples, payload, count, frames
