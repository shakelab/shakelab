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
Tests of the libmseed MiniSEED backend (cymseed3).

Round trips in every encoding and in MiniSEED 2 and 3, files exchanged
with the pure-Python backend, types of the samples, file objects and
errors. The conformance with libmseed on its test files is checked for
both backends in ``test_mseed_reading.py``.
"""

from __future__ import annotations

import io as _io
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.libio import mseed
from shakelab.signals.libio.mseed import MSeedWarning

cymseed = pytest.importorskip("shakelab.signals.libio.cymseed")

LIBMSEED = Path(__file__).parent / "data" / "mseed" / "libmseed"
TIME = "2026-01-01T00:00:00.123456"
SID = "AA.NOA.00.HHZ"


def _record(data, sid=SID, time=TIME, delta=0.005):
    return Record(data=data, delta=delta, time=time, sid=sid)


def _collection(*records):
    collection = StreamCollection()
    for record in records:
        collection.append(record)
    return collection


def _harmonic(dtype, count=3000):
    t = np.arange(count) * 0.005
    values = 20000.0 * np.sin(2.0 * np.pi * t) * np.exp(-0.2 * t)
    if np.issubdtype(dtype, np.integer):
        values = np.rint(values)
    return values.astype(dtype)


# ----------------------------------------------------------------------------
# Round trips
# ----------------------------------------------------------------------------

ENCODINGS = {1: np.int16, 3: np.int32, 4: np.float32, 5: np.float64,
             10: np.int32, 11: np.int32}


@pytest.mark.parametrize("msformat", [2, 3])
@pytest.mark.parametrize("reclen", [256, 512, 4096])
@pytest.mark.parametrize("encoding", sorted(ENCODINGS))
def test_round_trip(encoding, reclen, msformat):
    """Every encoding, record length and version returns the samples."""
    data = _harmonic(ENCODINGS[encoding])
    buffer = _io.BytesIO()
    cymseed.mswrite(buffer, _record(data), encoding=encoding,
                    reclen=reclen, msformat=msformat)

    (stream,) = cymseed.msread(buffer.getvalue())
    (record,) = stream.records()

    np.testing.assert_array_equal(record.data, data)
    assert record.sid == SID
    assert record.rate == 200.0
    assert abs(record.time - Date(TIME)) < 1e-9


@pytest.mark.parametrize("dtype,expected", [
    (np.int16, np.int32), (np.int32, np.int32), (np.int64, np.int32),
    (np.float32, np.float32), (np.float64, np.float64),
])
def test_sample_types(dtype, expected):
    """Default encodings keep float32 and float64; integers are int32."""
    buffer = _io.BytesIO()
    cymseed.mswrite(buffer, _record(_harmonic(dtype)))

    (record,) = cymseed.msread(buffer.getvalue()).records()

    assert record.data.dtype == expected


def test_nanoseconds_in_miniseed3():
    """MiniSEED 3 keeps the nanoseconds of the start time."""
    time = Date.from_unix_ns(1_767_225_600_123_456_789)
    buffer = _io.BytesIO()
    cymseed.mswrite(buffer, _record(_harmonic(np.int32), time=time),
                    msformat=3)

    (record,) = cymseed.msread(buffer.getvalue()).records()

    assert record.head.time.to_unix_ns() == time.to_unix_ns()


# ----------------------------------------------------------------------------
# Exchange with the pure-Python backend
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("encoding", sorted(ENCODINGS))
def test_python_reads_libmseed_files(encoding):
    """MiniSEED 2 written by libmseed is read by the Python backend."""
    data = _harmonic(ENCODINGS[encoding])
    buffer = _io.BytesIO()
    cymseed.mswrite(buffer, _record(data), encoding=encoding, reclen=512)

    (record,) = mseed.msread(buffer.getvalue()).records()

    np.testing.assert_array_equal(record.data, data)
    assert abs(record.time - Date(TIME)) < 1e-9


@pytest.mark.parametrize("encoding", sorted(ENCODINGS))
def test_backends_write_the_same_samples(encoding):
    """Both backends write files decoded to the same samples."""
    data = _harmonic(ENCODINGS[encoding])
    python, libmseed = _io.BytesIO(), _io.BytesIO()
    mseed.mswrite(python, _record(data), encoding=encoding)
    cymseed.mswrite(libmseed, _record(data), encoding=encoding)

    for content in (python.getvalue(), libmseed.getvalue()):
        for read in (mseed.msread, cymseed.msread):
            (record,) = read(content).records()
            np.testing.assert_array_equal(record.data, data)


def test_several_channels_and_segments():
    """Channels and gaps are kept; contiguous records are joined."""
    first = _record(_harmonic(np.int32, 100))
    second = _record(_harmonic(np.int32, 50),
                     time="2026-01-01T00:00:00.623456")
    gap = _record(_harmonic(np.int32, 50), time="2026-01-01T00:01:00")
    other = _record(_harmonic(np.float32, 80), sid="AA.NOA.00.HHE")
    buffer = _io.BytesIO()
    cymseed.mswrite(buffer, _collection(first, second, gap, other))

    collection = cymseed.msread(buffer.getvalue())
    counts = {stream.sid: [rec.nsamp for rec in stream.records()]
              for stream in collection}

    assert counts == {SID: [150, 50], "AA.NOA.00.HHE": [80]}


# ----------------------------------------------------------------------------
# Sources, targets and errors
# ----------------------------------------------------------------------------


def test_sources_and_targets(tmp_path):
    """Paths, bytes and binary files; owrite."""
    path = tmp_path / "out.mseed"
    record = _record(_harmonic(np.int32))
    cymseed.mswrite(path, record)
    with pytest.raises(FileExistsError):
        cymseed.mswrite(path, record)
    cymseed.mswrite(str(path), record, msformat=3, owrite=True)

    with open(path, "rb") as fid:
        from_file = cymseed.msread(fid)
    for source in (path, str(path), path.read_bytes()):
        assert cymseed.msread(source).records()[0].nsamp == 3000
    assert from_file.records()[0].nsamp == 3000


def test_append_to_collection():
    """Records are appended to a given collection."""
    collection = StreamCollection()
    path = LIBMSEED / "reference-testdata-steim2.mseed2"
    result = cymseed.msread(path, collection)

    assert result is collection
    assert collection.records()[0].nsamp == 499


def test_skipped_records():
    """Text and detection records are not waveforms."""
    for name in ("reference-testdata-text.mseed2",
                 "testdata-detection.record.mseed2"):
        assert len(cymseed.msread(LIBMSEED / name)) == 0


def test_invalid_data_warns():
    """Invalid data are skipped with an MSeedWarning."""
    content = (LIBMSEED / "reference-testdata-steim2.mseed2").read_bytes()
    with pytest.warns(MSeedWarning, match="invalid data"):
        collection = cymseed.msread(b"x" * 300 + content)
    assert collection.records()[0].nsamp == 499


@pytest.mark.parametrize("data,encoding,match", [
    (_harmonic(np.float64), 11, "integer samples"),
    (_harmonic(np.float64), 1, "integer samples"),
    (np.array([0, 40000], dtype=np.int32), 1, "range of encoding 1"),
    (np.array([0, 2**40], dtype=np.int64), None, "32-bit range"),
    (np.array([0, 2**40], dtype=np.int64), 3, "range of encoding 3"),
    (np.array(["a", "b"]), None, "cannot be written"),
])
def test_data_errors(data, encoding, match):
    """Data that the encoding cannot represent."""
    with pytest.raises(ValueError, match=match):
        cymseed.mswrite(_io.BytesIO(), _record(data), encoding=encoding)


def test_argument_errors():
    """Invalid containers, versions, encodings and record lengths."""
    record = _record(_harmonic(np.int32))
    with pytest.raises(TypeError):
        cymseed.mswrite(_io.BytesIO(), [record])
    with pytest.raises(ValueError, match="version"):
        cymseed.mswrite(_io.BytesIO(), record, msformat=4)
    with pytest.raises(ValueError):
        cymseed.mswrite(_io.BytesIO(), record, encoding=0)
    with pytest.raises(ValueError):
        cymseed.mswrite(_io.BytesIO(), record, reclen=500)
    with pytest.raises(TypeError):
        cymseed.mswrite(12, record)


def test_explicit_float_encodings_convert():
    """Encodings 4 and 5 convert integers and other floats."""
    data = _harmonic(np.int32)
    for encoding, dtype in ((4, np.float32), (5, np.float64)):
        buffer = _io.BytesIO()
        cymseed.mswrite(buffer, _record(data), encoding=encoding)
        (record,) = cymseed.msread(buffer.getvalue()).records()
        assert record.data.dtype == dtype
        np.testing.assert_array_equal(record.data, data)


# ----------------------------------------------------------------------------
# SeedLink packets
# ----------------------------------------------------------------------------


def test_seedlink_miniseed3_packet():
    """A MiniSEED 3 SeedLink packet (one record) is decoded."""
    from cymseed3 import pack_records

    from shakelab.signals.seedlink.client import _decode
    from shakelab.signals.seedlink.protocol import Packet

    data = _harmonic(np.int32, 400)
    (raw,) = pack_records(("AA", "NOA", "00", "HHZ"),
                          Date(TIME).to_unix_ns(), 200.0, data,
                          version=3, reclen=4096)
    packet = Packet(kind="data", sequence=1, station="AA_NOA", format="3",
                    subformat="D", payload=raw.raw)

    records, station = _decode(packet)

    assert station == "AA_NOA"
    (record,) = records
    assert record.sid == SID
    np.testing.assert_array_equal(record.data, data)
