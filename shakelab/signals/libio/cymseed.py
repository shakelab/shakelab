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
version. Its interface and results match those of the pure-Python
backend (:mod:`shakelab.signals.libio.mseed`), which documents the
behavior in detail.

Notes
-----
``cymseed3`` returns samples as Python lists without their original
type: integers are converted to int32 and floating-point values to
float64 (float32 data read by the pure-Python backend stay float32).

``cymseed3`` (up to 0.1.x) packs samples only as int32: floating-point
and int16 data are corrupted, and one encoding is used for the whole
file. Until it is fixed, MiniSEED 2 files that need other encodings
(1, 4, 5) are written by the pure-Python backend, and MiniSEED 3 is
limited to integer encodings (3, 10, 11).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

import numpy as np
from cymseed3 import MiniSeed

from shakelab.libutils.timeN import Date
from shakelab.signals.libio import mseed as _pymseed
from shakelab.signals.libio.mseed import (
    _default_encoding,
    _read_source,
    _split_fdsn_code,
    _validate_encoding,
    _validate_record_length,
)

if TYPE_CHECKING:
    from shakelab.signals.base import StreamCollection

__all__ = ["MSEED_VERSION", "msread", "mswrite"]

# Default MiniSEED version for writing.
MSEED_VERSION = 2

# Encodings that cymseed3 packs correctly (samples held as int32).
_CYMSEED3_ENCODINGS = (3, 10, 11)


def msread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
) -> StreamCollection:
    """
    Read MiniSEED 2 or 3 waveforms into a StreamCollection.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        MiniSEED file, its content, or a file opened in binary mode.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.

    Returns
    -------
    StreamCollection
        Collection with one record per contiguous segment of each
        channel, as joined by libmseed.
    """
    from shakelab.signals.base import Record, StreamCollection

    if stream_collection is None:
        stream_collection = StreamCollection()

    content = bytes(_read_source(source))
    if not content:
        return stream_collection

    ms = MiniSeed()
    ms.read(content)

    for segment in ms.export_records() or []:
        data = segment["data"]
        if not data or isinstance(data[0], str):
            continue
        dtype = np.int32 if isinstance(data[0], int) else np.float64
        stream_collection.append(
            Record(
                data=np.asarray(data, dtype=dtype),
                delta=1.0 / segment["rate"],
                time=Date(segment["starttime"]),
                sid=".".join(
                    segment[code].strip()
                    for code in ("network", "station", "location", "channel")
                ),
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
    Write waveforms to a MiniSEED file.

    Parameters
    ----------
    target : str or path-like
        Output file.
    stream_collection : StreamCollection
        Waveforms to write.
    encoding : int, optional
        MiniSEED encoding. If None (default), it follows the type of
        each record's data: STEIM2 for integers, 4 for float32 and 5
        for other floating-point types.
    reclen : int, optional
        Record length in bytes. Default is 512.
    msformat : {2, 3}, optional
        MiniSEED version; :data:`MSEED_VERSION` if None.

    Raises
    ------
    ValueError
        If the version, encoding, record length or a source identifier
        is invalid.
    NotImplementedError
        For MiniSEED 3 with encodings other than 3, 10 and 11, or with
        records needing different encodings (see the module notes).
    """
    if msformat is None:
        msformat = MSEED_VERSION
    if msformat not in (2, 3):
        raise ValueError(f"Unsupported MiniSEED version: {msformat}")
    if encoding is not None:
        _validate_encoding(encoding)
    _validate_record_length(reclen)

    records = [record for stream in stream_collection for record in stream]
    if not records:
        return

    encodings = {
        _default_encoding(np.asarray(record.data))
        if encoding is None else encoding
        for record in records
    }

    if len(encodings) > 1 or not encodings <= set(_CYMSEED3_ENCODINGS):
        if msformat == 2:
            _pymseed.mswrite(target, stream_collection, encoding, reclen)
            return
        raise NotImplementedError(
            "MiniSEED 3 is written by cymseed3 with one integer encoding "
            "(3, 10 or 11) only"
        )

    (record_encoding,) = encodings
    ms = MiniSeed()

    for record in records:
        network, station, location, channel = _split_fdsn_code(record.sid)
        ms.import_record(
            {
                "network": network,
                "station": station,
                "location": location,
                "channel": channel,
                "starttime": record.time.get_date("iso8601"),
                "rate": record.rate,
                "nsamp": record.nsamp,
                "data": np.asarray(record.data),
            },
            encoding=record_encoding,
            reclen=reclen,
        )

    content = ms.write(msformat, reclen, record_encoding)
    with open(target, "wb") as fid:
        fid.write(content)
