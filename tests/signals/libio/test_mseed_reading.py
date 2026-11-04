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
Tests of the pure-Python MiniSEED backend on reference files.

The files in ``data/mseed/libmseed`` are the test data of libmseed
(EarthScope); those in ``data/mseed/synth`` hold a Fibonacci series in
every supported encoding, byte order and record length, with the
expected values in ASCII. The tests also fix the behavior on
non-waveform records, unsupported encodings, time corrections,
blockettes 100/1001 and the choice of the writing encoding.
"""

from __future__ import annotations

import hashlib
import io as _io
import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.libio import mseed
from shakelab.signals.libio.mseed import (
    ENCODING_FLOAT32,
    ENCODING_FLOAT64,
    ENCODING_STEIM1,
    ENCODING_STEIM2,
    MSeedWarning,
    MSRecord,
    msrawread,
    msread,
    mswrite,
)

DATA = Path(__file__).parent / "data" / "mseed"
LIBMSEED = DATA / "libmseed"
SYNTH = sorted((DATA / "synth").glob("*.mseed"))
SIGNALS = Path(__file__).parents[1]
REFERENCE = json.loads((DATA / "libmseed_reference.json").read_text())


def _reference(name):
    """Data of a libmseed reference file (single record)."""
    collection = msread(LIBMSEED / f"reference-testdata-{name}.mseed2")
    assert len(collection) == 1 and len(collection[0]) == 1
    return collection[0][0]


# ----------------------------------------------------------------------------
# Conformance with libmseed
# ----------------------------------------------------------------------------


def _backend_reader(backend):
    """msread of a backend (the libmseed one needs cymseed3)."""
    if backend == "python":
        return msread
    pytest.importorskip("cymseed3")
    from shakelab.signals.libio import cymseed

    return cymseed.msread


@pytest.mark.parametrize("backend", ["python", "libmseed"])
@pytest.mark.parametrize(
    "path", sorted(REFERENCE["files"]), ids=lambda p: Path(p).name
)
def test_libmseed_reference(path, backend):
    """
    Both backends read what libmseed reads: identifiers, start times,
    rates and samples (compared through a hash, see the reference
    file).
    """
    read = _backend_reader(backend)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MSeedWarning)
        records = [rec for stream in read(SIGNALS / path) for rec in stream]
    expected = REFERENCE["files"][path]

    assert len(records) == len(expected)
    for record, segment in zip(
        sorted(records, key=lambda r: r.sid),
        sorted(expected, key=lambda s: s["sid"]),
    ):
        samples = np.asarray(record.data, dtype=">f8").tobytes()
        assert record.sid == segment["sid"]
        assert abs(record.time - Date(segment["starttime"])) < 1e-6
        assert record.rate == pytest.approx(segment["rate"], rel=1e-9)
        assert record.nsamp == segment["nsamp"]
        assert hashlib.sha256(samples).hexdigest() == segment["sha256"]


# ----------------------------------------------------------------------------
# Decoding
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("path", SYNTH, ids=lambda p: p.stem)
def test_synthetic_files(path):
    """Every encoding, byte order and record length gives the series."""
    expected = np.loadtxt(path.with_suffix(".ascii"))
    record = msread(path)[0][0]

    np.testing.assert_allclose(record.data, expected)
    assert record.data.dtype.isnative


@pytest.mark.parametrize(
    "name", ["int16", "float32", "float64", "steim1", "steim2", "defaults"]
)
def test_libmseed_encodings_agree(name):
    """
    The same series is decoded from every encoding (the integer files
    hold the floating-point series truncated).
    """
    expected = _reference("int32").data
    record = _reference(name)
    count = record.nsamp

    assert record.sid == "XX.TEST..BHZ"
    assert record.rate == 40.0
    if name.startswith("float"):
        np.testing.assert_array_equal(np.trunc(record.data), expected)
        np.testing.assert_allclose(
            record.data, _reference("float64").data, rtol=1e-6
        )
    else:
        np.testing.assert_array_equal(record.data, expected[:count])


@pytest.mark.parametrize("name", ["steim1", "steim2"])
def test_little_endian_steim(name):
    """Little-endian STEIM records decode as the big-endian ones."""
    with pytest.warns(MSeedWarning, match="trailing bytes"):
        little = msread(LIBMSEED / f"reference-testdata-{name}-LE.mseed2")
    little = little[0][0]
    big = _reference(name)

    np.testing.assert_array_equal(little.data, big.data)


def test_forced_byte_order():
    """A wrong imposed byte order is reported as an invalid record."""
    path = LIBMSEED / "reference-testdata-steim1.mseed2"

    assert msread(path, byte_order="be")[0][0].nsamp == 500
    with pytest.raises(ValueError, match="byte offset 0"):
        msread(path, byte_order="le")


def test_native_integer_type():
    """int16 samples are returned in native byte order."""
    record = _reference("int16")
    assert record.data.dtype == np.int16
    assert record.data.dtype.isnative


# ----------------------------------------------------------------------------
# Time and sampling rate
# ----------------------------------------------------------------------------


def test_unapplied_time_correction_and_negative_multiplier():
    """
    The time correction is added when not applied, and a negative rate
    multiplier divides (32760 / 819 = 40 Hz).
    """
    path = LIBMSEED / "testdata-unapplied-timecorrection.mseed2"
    raw = msrawread(path)[0]

    assert raw.header["SAMPLE_RATE_MULTIPLIER"] == -819
    assert raw.rate == pytest.approx(40.0)
    assert raw.header["TIME_CORRECTION"] == 10000
    assert raw.time == Date("2003-05-29T02:13:23.0434")

    # Once marked as applied, the correction is not added again.
    raw.header["ACTIVITY_FLAGS"] |= mseed.TIME_CORRECTION_APPLIED
    assert raw.time == Date("2003-05-29T02:13:22.0434")


def test_blockette_1001_microseconds():
    """The microseconds of blockette 1001 are part of the start time."""
    raw = msrawread(LIBMSEED / "testdata-3channel-signal.mseed2")[0]
    micro = raw.blockette[1001]["MICRO_SEC"]
    header_only = Date("2010-02-27T06:50:00.0695")

    assert micro != 0
    assert raw.time - header_only == pytest.approx(micro * 1e-6)


def test_blockette_100_rate(tmp_path):
    """Blockette 100 takes precedence over the nominal rate."""
    raw = msrawread(LIBMSEED / "reference-testdata-int32.mseed2")[0]
    raw.blockette[100] = {"SAMPLE_RATE": 39.999, "FLAGS": 0}

    content, _ = raw.to_bytes(record_length=4096)
    restored = MSRecord.from_bytes(content)

    assert restored.rate == pytest.approx(39.999, rel=1e-6)


@pytest.mark.parametrize(
    "factor, multiplier, rate",
    [(40, 1, 40.0), (20, 2, 40.0), (32760, -819, 40.0), (-10, 1, 0.1),
     (-10, -2, 0.05), (0, 1, 0.0)],
)
def test_nominal_rate(factor, multiplier, rate):
    """Rate factor and multiplier as in the SEED manual."""
    assert mseed._nominal_rate(factor, multiplier) == pytest.approx(rate)


# ----------------------------------------------------------------------------
# Records that are not waveforms, or not supported
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["reference-testdata-text", "testdata-detection.record"]
)
def test_non_waveform_records_skipped(name):
    """Log and detection records are skipped silently by msread."""
    path = LIBMSEED / f"{name}.mseed2"

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        collection = msread(path)
    raw = msrawread(path)

    assert len(collection) == 0
    assert raw and not any(record.is_waveform for record in raw)


def test_text_record_content():
    """Text records are decoded as strings by the raw reader."""
    raw = msrawread(LIBMSEED / "reference-testdata-text.mseed2")[0]
    assert isinstance(raw.data, str) and len(raw.data) == raw.nsamp


@pytest.mark.parametrize(
    "name, encoding",
    [("CDSN", 16), ("DWWSSN", 32), ("SRO", 30),
     ("GEOSCOPE-16bit-3exp-encoded", 14)],
)
def test_unsupported_encodings(name, encoding):
    """Records with unsupported encodings are skipped with one warning."""
    path = LIBMSEED / f"testdata-encoding-{name}.mseed2"

    with pytest.warns(MSeedWarning, match=f"encoding {encoding}") as caught:
        collection = msread(path)

    assert len(collection) == 0
    assert len(caught) == 1
    assert msrawread(path)[0].data is None


@pytest.mark.parametrize(
    "name", ["testdata-invalid-blockette-offsets",
             "testdata-no-blockette1000-steim1"],
)
def test_invalid_records(name):
    """Invalid records raise ValueError with their offset."""
    with pytest.raises(ValueError, match="byte offset 0"):
        msread(LIBMSEED / f"{name}.mseed2")


def test_steim_last_sample_mismatch():
    """A wrong last sample in the STEIM frames warns, data are kept."""
    content = bytearray(
        (LIBMSEED / "reference-testdata-steim1.mseed2").read_bytes()
    )
    raw = MSRecord.from_bytes(content)
    xn = raw.header["OFFSET_TO_BEGINNING_OF_DATA"] + 8
    content[xn:xn + 4] = (12345).to_bytes(4, "big", signed=True)

    with pytest.warns(MSeedWarning, match="differs from the declared"):
        record = msread(bytes(content))[0][0]

    np.testing.assert_array_equal(
        record.data[:raw.nsamp], _reference("steim1").data[:raw.nsamp]
    )


def test_trailing_bytes():
    """Zero padding is ignored; other trailing bytes with a warning."""
    content = (LIBMSEED / "reference-testdata-steim2.mseed2").read_bytes()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert msread(content + bytes(100))[0][0].nsamp == 499

    with pytest.warns(MSeedWarning, match="3 trailing bytes"):
        msread(content + b"abc")

    with pytest.raises(ValueError, match="byte offset"):
        msread(content + b"x" * 60)


# ----------------------------------------------------------------------------
# Joining records
# ----------------------------------------------------------------------------


def _collection(*records):
    """StreamCollection with the given records."""
    collection = StreamCollection()
    for record in records:
        collection.append(record)
    return collection


def test_records_joined_and_gaps_kept():
    """Contiguous records form one segment; a gap starts another."""
    data = np.arange(3000, dtype=np.int32)
    first = Record(data=data, delta=0.01, time="2026-01-01T00:00:00",
                   sid="XX.STA..HHZ")
    second = Record(data=data, delta=0.01, time="2026-01-01T00:01:00",
                    sid="XX.STA..HHZ")
    buffer = _io.BytesIO()
    mswrite(buffer, _collection(first, second), reclen=256)

    assert len(msrawread(buffer.getvalue())) > 10
    stream = msread(buffer.getvalue())[0]

    assert len(stream) == 2
    np.testing.assert_array_equal(stream[0].data, data)
    assert stream[1].time == Date("2026-01-01T00:01:00")


def test_interleaved_channels():
    """Records of different channels in one file are separated."""
    path = LIBMSEED / "testdata-3channel-signal.mseed2"
    collection = msread(path)

    assert collection.sid == ["IU.COLA.00.LH1", "IU.COLA.00.LH2",
                              "IU.COLA.00.LHZ"]
    assert all(len(stream) == 1 for stream in collection)


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "dtype, encoding",
    [(np.int32, ENCODING_STEIM2), (np.int16, ENCODING_STEIM2),
     (np.float32, ENCODING_FLOAT32), (np.float64, ENCODING_FLOAT64)],
)
def test_default_encoding(dtype, encoding):
    """The writing encoding follows the type of the data."""
    data = (np.arange(500) - 250).astype(dtype)
    record = Record(data=data, delta=0.005, time="2026-01-01T00:00:00",
                    sid="XX.STA.00.HNZ")
    buffer = _io.BytesIO()
    mswrite(buffer, _collection(record))

    raw = msrawread(buffer.getvalue())
    restored = msread(buffer.getvalue())[0][0]

    assert {item.encoding for item in raw} == {encoding}
    np.testing.assert_array_equal(restored.data, data)


def test_steim_rejects_floating_point():
    """STEIM chosen explicitly for float data gives a clear error."""
    record = Record(data=np.linspace(0.0, 1.0, 10), delta=0.01,
                    time="2026-01-01T00:00:00", sid="XX.STA..HHZ")

    with pytest.raises(ValueError, match="encoding 4 or 5"):
        mswrite(_io.BytesIO(), _collection(record),
                encoding=ENCODING_STEIM1)


def test_microsecond_start_time():
    """Microseconds go to blockette 1001 and are read back."""
    time = "2026-03-01T10:20:30.123457"
    record = Record(data=np.arange(2000, dtype=np.int32), delta=0.01,
                    time=time, sid="XX.STA..HHZ")
    buffer = _io.BytesIO()
    mswrite(buffer, _collection(record))

    raw = msrawread(buffer.getvalue())
    restored = msread(buffer.getvalue())[0][0]

    assert raw[0].header["MSECONDS"] == 1235
    assert raw[0].blockette[1001]["MICRO_SEC"] == -43
    assert restored.time == Date(time)
    assert len(msread(buffer.getvalue())[0]) == 1


def test_whole_tenths_without_blockette_1001():
    """Blockette 1001 is written only when needed."""
    record = Record(data=np.arange(10, dtype=np.int32), delta=0.01,
                    time="2026-03-01T10:20:30.1234", sid="XX.STA..HHZ")
    buffer = _io.BytesIO()
    mswrite(buffer, _collection(record), encoding=ENCODING_STEIM2)

    assert 1001 not in msrawread(buffer.getvalue())[0].blockette


@pytest.mark.parametrize(
    "unix_ns, fields, micro",
    [
        (0, (1970, 1, 0, 0, 0, 0), 0),
        (1_000_049_000, (1970, 1, 0, 0, 1, 0), 49),
        (1_000_050_000, (1970, 1, 0, 0, 1, 1), -50),
        (86_399_999_960_000, (1970, 2, 0, 0, 0, 0), -40),
    ],
)
def test_time_fields(unix_ns, fields, micro):
    """Header time to 0.1 ms, microseconds from -50 to 49."""
    header, offset = mseed._mseed_time_fields(unix_ns)
    names = ("YEAR", "DAY", "HOURS", "MINUTES", "SECONDS", "MSECONDS")

    assert tuple(header[name] for name in names) == fields
    assert offset == micro


def test_file_objects(tmp_path):
    """Paths, bytes and binary file objects are accepted."""
    path = LIBMSEED / "reference-testdata-steim2.mseed2"
    expected = msread(path)[0][0].data

    with open(path, "rb") as fid:
        from_file = msread(fid)[0][0].data
    from_bytes = msread(path.read_bytes())[0][0].data

    np.testing.assert_array_equal(from_file, expected)
    np.testing.assert_array_equal(from_bytes, expected)

    with pytest.raises(TypeError):
        msread(12)


def test_records_out_of_order():
    """Records in mixed order and lengths are sorted and joined."""
    path = LIBMSEED / "testdata-oneseries-mixedlengths-mixedorder.mseed2"
    raw = msrawread(path)
    stream = msread(path)[0]

    assert len({item.record_length for item in raw}) == 7
    assert len(stream) == 1
    assert stream[0].nsamp == sum(item.nsamp for item in raw)
    assert stream[0].time == min(item.time for item in raw)
