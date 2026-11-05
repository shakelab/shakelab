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
Tests of :mod:`shakelab.signals.libio.simpleascii`.

The files in ``data/simpleascii`` were written by ``mseed2ascii`` 2.6
(EarthScope) from the libmseed test files in ``data/mseed/libmseed``;
the samples read are compared with the libmseed reference, and the
files written from the MiniSEED records with those of ``mseed2ascii``.
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
from shakelab.signals.base import Record, Stream, StreamCollection
from shakelab.signals.io import reader, writer
from shakelab.signals.libio.mseed import msread
from shakelab.signals.libio.simpleascii import (
    SimpleAscii,
    SimpleAsciiWarning,
    simpleasciiread,
    simpleasciiwrite,
)

DATA = Path(__file__).parent / "data"
ASCII = DATA / "simpleascii"
MSEED = DATA / "mseed" / "libmseed"
REFERENCE = DATA / "mseed" / "libmseed_reference.json"

# Simple ASCII file: (MiniSEED source, style, columns).
FILES = {
    "int32.slist.txt": ("reference-testdata-int32", "slist", 1),
    "int32.slist6.txt": ("reference-testdata-int32", "slist", 6),
    "int32.tspair.txt": ("reference-testdata-int32", "tspair", 1),
    "float32.slist.txt": ("reference-testdata-float32", "slist", 1),
    "float64.slist.txt": ("reference-testdata-float64", "slist", 1),
    "float64.tspair.txt": ("reference-testdata-float64", "tspair", 1),
    "oddrate.slist.txt": ("reference-testdata-oddrate", "slist", 1),
    "3channel.slist.txt": ("testdata-3channel-signal", "slist", 1),
}

TEXT = (
    "TIMESERIES XX_TEST__LOG_D, 27 samples, 0 sps, "
    "2026-01-01T00:00:00.000000, SLIST, ASCII, Counts\n"
    "Station rebooted\nGPS locked\n"
)


def _records(collection):
    """Records of a collection, in identifier order."""
    return sorted(
        (rec for stream in collection for rec in stream),
        key=lambda rec: rec.sid,
    )


def _mseed_records(name):
    """Records of a libmseed test file, as written by mseed2ascii."""
    collection = msread(MSEED / f"{name}.mseed2")
    for rec in _records(collection):
        rec.head.units = "Counts"
    return collection


def _segment(n=4, **header):
    """Text of a small INTEGER segment."""
    fields = {
        "source": "XX_STA_00_HHZ", "samples": n, "sps": 100,
        "time": "2026-01-01T00:00:00.000000", "style": "SLIST",
    }
    fields.update(header)
    return (
        f"TIMESERIES {fields['source']}, {fields['samples']} samples, "
        f"{fields['sps']} sps, {fields['time']}, {fields['style']}, "
        "INTEGER, Counts\n" + "\n".join(str(k) for k in range(n)) + "\n"
    )


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(FILES))
def test_read_as_libmseed(name):
    """Identifiers, times and samples as read by libmseed."""
    with open(REFERENCE, encoding="utf-8") as fid:
        reference = json.load(fid)["files"]
    source = FILES[name][0]
    expected = sorted(
        reference[f"libio/data/mseed/libmseed/{source}.mseed2"],
        key=lambda seg: seg["sid"],
    )
    records = _records(simpleasciiread(ASCII / name))

    assert len(records) == len(expected)
    for rec, seg in zip(records, expected):
        assert rec.sid == seg["sid"]
        assert rec.nsamp == seg["nsamp"]
        assert rec.rate == pytest.approx(seg["rate"], rel=1e-12)
        assert abs(
            rec.time.to_unix_ns() - Date(seg["starttime"]).to_unix_ns()
        ) < 1000
        if rec.data.dtype.kind == "i":
            digest = hashlib.sha256(
                rec.data.astype(">f8").tobytes()
            ).hexdigest()
            assert digest == seg["sha256"]


@pytest.mark.parametrize(
    "name, rtol",
    [("float32.slist.txt", 1e-7), ("float64.slist.txt", 1e-9),
     ("float64.tspair.txt", 1e-9)],
)
def test_float_samples_within_rounding(name, rtol):
    """Floats agree with MiniSEED to the digits written by the format."""
    read = _records(simpleasciiread(ASCII / name))[0].data
    expected = _records(_mseed_records(FILES[name][0]))[0].data

    assert read.dtype == np.float64
    np.testing.assert_allclose(read, expected, rtol=rtol, atol=0)


def test_integer_types():
    """Integers are int32 when they fit, int64 otherwise."""
    assert simpleasciiread(_segment())[0][0].data.dtype == np.int32

    text = _segment(1).replace("\n0\n", "\n3000000000\n")
    rec = simpleasciiread(text)[0][0]
    assert rec.data.dtype == np.int64
    assert rec.data[0] == 3_000_000_000


def test_record_header():
    """Units as written, quality code in the extra fields."""
    rec = _records(simpleasciiread(ASCII / "3channel.slist.txt"))[0]

    assert rec.sid == "IU.COLA.00.LH1"
    assert rec.time == Date("2010-02-27T06:50:00.069539Z")
    assert rec.head.units == "Counts"
    assert rec.head.extra["simpleascii"] == {"QUALITY": "M"}


def test_source_name_without_quality():
    """The quality code is optional."""
    rec = simpleasciiread(_segment())[0][0]

    assert rec.sid == "XX.STA.00.HHZ"
    assert "simpleascii" not in rec.head.extra


def test_invalid_source_name():
    """A source name that is not Net_Sta_Loc_Chan gives no identifier."""
    segment = SimpleAscii.from_text(_segment(source="XX_STA"))[0]

    assert segment.sid is None
    assert segment.quality is None


def test_several_segments_of_one_channel():
    """Segments of the same channel are separate records."""
    text = _segment() + _segment(time="2026-01-01T00:00:10.000000")
    stream = simpleasciiread(text)["XX.STA.00.HHZ"]

    assert len(stream) == 2
    assert stream[1].time == Date("2026-01-01T00:00:10Z")


def test_append_to_collection():
    """Records are appended to a given collection."""
    collection = simpleasciiread(_segment())
    simpleasciiread(
        _segment(source="XX_STA_00_HHN"), stream_collection=collection
    )

    assert sorted(st.sid for st in collection) == [
        "XX.STA.00.HHN", "XX.STA.00.HHZ"
    ]


def test_read_sources():
    """Paths, content, bytes and file objects are accepted."""
    path = ASCII / "int32.slist.txt"
    text = path.read_text()
    expected = simpleasciiread(path)[0][0].data

    for source in (str(path), text, text.encode(), _io.StringIO(text),
                   _io.BytesIO(text.encode())):
        assert np.array_equal(simpleasciiread(source)[0][0].data, expected)
    with pytest.raises(TypeError):
        simpleasciiread(42)


def test_sample_count_mismatch_warns():
    """A number of samples different from the header gives a warning."""
    with pytest.warns(SimpleAsciiWarning, match="5 samples declared"):
        rec = simpleasciiread(_segment(4, samples=5))[0][0]
    assert rec.nsamp == 4


def test_tspair_start_mismatch_warns():
    """A first pair far from the start time gives a warning."""
    text = ASCII.joinpath("int32.tspair.txt").read_text().replace(
        "2012-05-12T00:00:00.000000  ", "2012-05-12T00:00:01.000000  ", 1
    )
    with pytest.warns(SimpleAsciiWarning, match="first time-sample pair"):
        simpleasciiread(text)


@pytest.mark.parametrize(
    "text, match",
    [
        ("1\n2\n", "no TIMESERIES"),
        ("junk\n" + _segment(), "no TIMESERIES"),
        ("TIMESERIES XX_A__HHZ, 1 samples, 1 sps\n1\n", "header fields"),
        (_segment(time="yesterday"), "invalid header"),
        (_segment(style="COLUMNS"), "unknown format"),
        (_segment(1, style="TSPAIR"), "incomplete"),
        (_segment(1) + "x\n", "invalid sample"),
    ],
)
def test_invalid_files(text, match):
    """Invalid files are rejected with the line of the problem."""
    with pytest.raises(ValueError, match=match):
        simpleasciiread(text)


def test_text_segment():
    """Text segments are read by SimpleAscii only."""
    segment = SimpleAscii.from_text(TEXT)[0]

    assert not segment.is_timeseries
    assert segment.data == "Station rebooted\nGPS locked"
    assert segment.quality == "D"
    assert segment.to_text() == TEXT
    with pytest.raises(ValueError, match="holds text"):
        simpleasciiread(TEXT)


def test_zero_rate_rejected():
    """A waveform needs a positive sampling rate."""
    with pytest.raises(ValueError, match="sampling rate"):
        simpleasciiread(_segment(sps=0))


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(FILES))
def test_write_as_mseed2ascii(name):
    """Records read from MiniSEED are written as mseed2ascii writes them."""
    source, style, columns = FILES[name]
    expected = (ASCII / name).read_text()
    quality = expected.split(",")[0].rsplit("_", 1)[1]
    collection = _mseed_records(source)
    for rec in _records(collection):
        rec.head.extra["simpleascii"] = {"QUALITY": quality}

    buffer = _io.StringIO()
    simpleasciiwrite(buffer, collection, style=style, columns=columns)
    assert buffer.getvalue() == expected


@pytest.mark.parametrize("style", ["slist", "tspair"])
def test_round_trip(style):
    """Read, write and read again."""
    collection = simpleasciiread(ASCII / "3channel.slist.txt")
    buffer = _io.StringIO()
    simpleasciiwrite(buffer, collection, style=style, columns=3)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        again = simpleasciiread(buffer.getvalue())
    for rec, other in zip(_records(collection), _records(again)):
        assert other.sid == rec.sid
        assert other.time == rec.time
        assert other.delta == rec.delta
        assert other.head.units == rec.head.units
        assert other.head.extra == rec.head.extra
        np.testing.assert_array_equal(other.data, rec.data)


def test_write_record_and_stream():
    """A Record or a Stream can be written directly."""
    rec = Record(data=np.arange(3), delta=0.5, time="2026-01-01T00:00:00Z",
                 sid="XX.A..HHZ")
    stream = Stream("XX.A..HHZ")
    stream.append(rec)
    for data in (rec, stream):
        buffer = _io.StringIO()
        simpleasciiwrite(buffer, data)
        assert buffer.getvalue().startswith(
            "TIMESERIES XX_A__HHZ, 3 samples, 2 sps, "
            "2026-01-01T00:00:00.000000, SLIST, INTEGER, Unknown\n"
        )
    with pytest.raises(TypeError):
        simpleasciiwrite(_io.StringIO(), [rec])


def test_rate_and_time_precision():
    """The rate is written exactly, the start time rounded to 1 us."""
    rec = Record(data=[0.5, 1.5], delta=60.0, sid="XX.A..LHZ",
                 time="2026-01-01T00:00:00.1234567Z", units="M/S")
    buffer = _io.StringIO()
    simpleasciiwrite(buffer, rec)
    text = buffer.getvalue()

    assert "0.016666666666666666 sps, 2026-01-01T00:00:00.123457, " \
        "SLIST, FLOAT" in text
    again = simpleasciiread(text)[0][0]
    assert again.delta == 60.0
    assert again.head.units == "M/S"


@pytest.mark.parametrize(
    "change, error, match",
    [
        ({"sid": None}, ValueError, "source identifier"),
        ({"time": None}, ValueError, "source identifier"),
        ({"units": "m per s"}, ValueError, "spaces"),
    ],
)
def test_write_invalid_record(change, error, match):
    """Records that cannot be written are rejected."""
    rec = Record(data=[1, 2], delta=0.01, sid="XX.A..HHZ",
                 time="2026-01-01T00:00:00Z")
    if "sid" in change:
        rec.sid = None
    if "time" in change:
        rec.head.time = None
    if "units" in change:
        rec.head.units = change["units"]
    with pytest.raises(error, match=match):
        simpleasciiwrite(_io.StringIO(), rec)


@pytest.mark.parametrize(
    "options, error",
    [({"style": "columns"}, ValueError), ({"columns": 0}, ValueError),
     ({"columns": 1.5}, TypeError)],
)
def test_write_invalid_options(options, error):
    """Unknown styles and invalid column numbers are rejected."""
    rec = simpleasciiread(_segment())[0][0]
    with pytest.raises(error):
        simpleasciiwrite(_io.StringIO(), rec, **options)


def test_write_file(tmp_path):
    """Files are not overwritten unless requested."""
    rec = simpleasciiread(_segment())[0][0]
    path = tmp_path / "out.txt"
    simpleasciiwrite(path, rec)

    with pytest.raises(FileExistsError):
        simpleasciiwrite(path, rec)
    simpleasciiwrite(path, rec, owrite=True)
    segment = SimpleAscii.read(path)[0]
    segment.write(tmp_path / "copy.txt", style="tspair")
    assert SimpleAscii.read(tmp_path / "copy.txt")[0].header["FORMAT"] \
        == "TSPAIR"


def test_header_must_be_a_dict():
    """The header of a segment is a dictionary."""
    with pytest.raises(TypeError):
        SimpleAscii(header=[("SOURCE", "XX_A__HHZ")])


# ----------------------------------------------------------------------------
# io
# ----------------------------------------------------------------------------


def test_io_reader_and_writer(tmp_path):
    """``io`` reads and writes the format by name."""
    collection = reader(str(ASCII / "3channel.slist.txt"),
                        format="simpleascii")
    path = tmp_path / "out.txt"
    writer(collection, str(path), format="simpleascii", style="tspair")

    again = reader(str(path), format="simpleascii")
    assert isinstance(again, StreamCollection)
    assert [st.sid for st in again] == [st.sid for st in collection]
    assert "TSPAIR" in path.read_text().splitlines()[0]
