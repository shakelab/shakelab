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
MiniSEED reader and writer based on libmseed.

This backend wraps the EarthScope libmseed library (version 3) through
the ``cymseed3`` package. It reads MiniSEED 2 and 3 and writes either
version in every encoding. Its interface and results match those of
the pure-Python backend (:mod:`shakelab.signals.libio.mseed`), which
documents the behavior in detail.

Notes
-----
Samples keep the type decoded by libmseed: int32 for integer encodings
(also for int16 records, which the pure-Python backend returns as
int16), float32 and float64. Records without samples and text records
are skipped, as by the pure-Python backend.

When writing with ``encoding=None``, libmseed chooses the encoding of
each segment with the same rule as the pure-Python backend: STEIM2 for
integers whose differences fit in 30 bits, STEIM1 if they fit in 32
bits, otherwise int32; 4 for float32, 5 for other floating-point data.
Contiguous records of the same channel are written as one segment.
MiniSEED 2 stores times to the microsecond.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
from cymseed3 import MiniSeed, MiniSeedWarning

from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import warn_external
from shakelab.signals.libio.mseed import (
    INT32_MAX,
    INT32_MIN,
    MSeedWarning,
    _BinaryOutput,
    _split_fdsn_code,
    _validate_encoding,
    _validate_record_length,
)

if TYPE_CHECKING:
    from shakelab.signals.base import StreamCollection

__all__ = ["MSEED_VERSION", "msread", "mswrite"]

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

# Default MiniSEED version for writing.
MSEED_VERSION = 2

# Sample types of the encodings, as written by cymseed3.
_ENCODING_TYPES = {
    1: np.int16,
    3: np.int32,
    4: np.float32,
    5: np.float64,
    10: np.int32,
    11: np.int32,
}


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def msread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
) -> StreamCollection:
    """
    Read MiniSEED 2 or 3 waveforms into a StreamCollection.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        MiniSEED file (read through ``mmap``), its content (e.g. a
        SeedLink packet), or a file opened in binary mode.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.

    Returns
    -------
    StreamCollection
        Collection with one record per contiguous segment of each
        channel, as joined by libmseed.

    Warns
    -----
    MSeedWarning
        For invalid or truncated records skipped by libmseed.
    """
    from shakelab.signals.base import Record, StreamCollection

    if stream_collection is None:
        stream_collection = StreamCollection()

    ms = MiniSeed()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", MiniSeedWarning)
        ms.read(source)
    for warning in caught:
        if issubclass(warning.category, MiniSeedWarning):
            warn_external(str(warning.message), MSeedWarning)
        else:
            warnings.warn_explicit(warning.message, warning.category,
                                   warning.filename, warning.lineno)

    for segment in ms.segments():
        if segment.sampletype == "t" or segment.rate <= 0.0:
            continue
        stream_collection.append(
            Record(
                data=np.frombuffer(segment.data,
                                   dtype=segment.data.typecode),
                delta=1.0 / segment.rate,
                time=Date.from_unix_ns(segment.start_ns),
                sid=".".join(segment.codes),
            )
        )

    return stream_collection


def mswrite(
    target: Any,
    stream_collection: Any,
    encoding: Optional[int] = None,
    reclen: int = 512,
    msformat: Optional[int] = None,
    owrite: bool = False,
) -> None:
    """
    Write waveforms to a MiniSEED 2 or 3 file.

    Parameters
    ----------
    target : str, path-like or binary file object
        Output file, or a file opened in binary mode.
    stream_collection : StreamCollection, Stream or Record
        Waveforms to write.
    encoding : int, optional
        MiniSEED encoding of all the records: 1 (int16), 3 (int32), 4
        (float32), 5 (float64), 10 (STEIM1) or 11 (STEIM2). If None
        (default), it follows the type of the data of each segment
        (see the module notes).
    reclen : int, optional
        Record length in bytes. Default is 512.
    msformat : {2, 3}, optional
        MiniSEED version; :data:`MSEED_VERSION` if None.
    owrite : bool, optional
        If True, overwrite an existing file. Default is False.

    Raises
    ------
    TypeError
        If ``stream_collection`` is not a waveform container.
    FileExistsError
        If the file exists and ``owrite`` is False.
    ValueError
        If the version, encoding, record length or a source identifier
        is invalid, or if the data cannot be represented with the
        encoding (e.g. floating-point data with STEIM).
    """
    from shakelab.signals.base import Record, Stream, StreamCollection

    if not isinstance(stream_collection, (Record, Stream, StreamCollection)):
        raise TypeError(
            "expected a StreamCollection, a Stream or a Record, got "
            f"{type(stream_collection).__name__}"
        )
    if msformat is None:
        msformat = MSEED_VERSION
    if msformat not in (2, 3):
        raise ValueError(f"Unsupported MiniSEED version: {msformat}")
    if encoding is not None:
        _validate_encoding(encoding)
        if encoding not in _ENCODING_TYPES:
            raise ValueError(
                f"encoding {encoding} is not supported for waveforms"
            )
    _validate_record_length(reclen)

    ms = MiniSeed()
    for record in stream_collection.records():
        ms.add(
            _split_fdsn_code(record.sid),
            record.head.time.to_unix_ns(),
            float(record.rate),
            _samples(np.asarray(record.data), encoding),
            encoding=encoding,
        )

    with _BinaryOutput(target, owrite) as fid:
        fid.write(ms.write(version=msformat, reclen=reclen))


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _samples(data: np.ndarray, encoding: Optional[int]) -> np.ndarray:
    """
    Convert samples to the type that cymseed3 writes with an encoding
    (integers to int32 if ``encoding`` is None).

    Raises
    ------
    ValueError
        If integer encodings get non-integer data, or the samples
        exceed the range of the encoding.
    """
    if data.ndim != 1:
        raise ValueError("MiniSEED samples must be one-dimensional")
    if encoding is None:
        if np.issubdtype(data.dtype, np.integer):
            dtype = np.int32
        elif data.dtype == np.float32:
            dtype = np.float32
        elif np.issubdtype(data.dtype, np.floating):
            dtype = np.float64
        else:
            raise ValueError(
                f"data of type {data.dtype} cannot be written as MiniSEED"
            )
    else:
        dtype = _ENCODING_TYPES[encoding]

    if np.issubdtype(dtype, np.integer):
        if not np.issubdtype(data.dtype, np.integer):
            raise ValueError(
                f"encoding {encoding} requires integer samples (use "
                "encoding 4 or 5, or None, for floating-point data)"
            )
        info = np.iinfo(dtype)
        low, high = (info.min, info.max) if dtype == np.int16 else (
            INT32_MIN, INT32_MAX)
        if data.size and (data.min() < low or data.max() > high):
            raise ValueError(
                f"samples exceed the range of encoding {encoding}"
                if encoding is not None else
                "integer samples exceed the 32-bit range"
            )
    return np.ascontiguousarray(data, dtype=dtype)
