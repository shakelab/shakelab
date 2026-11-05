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
Tests of :mod:`shakelab.signals.libio.ascii`.

The tables are written in the tests, in the layouts of files produced
by simple programs and spreadsheets; the MODE records come from
``tests/signals/data``.
"""

from __future__ import annotations

import io as _io
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.io import reader, writer
from shakelab.signals.libio.ascii import (
    AsciiSignalWarning,
    asciiread,
    asciiwrite,
)
from shakelab.signals.libio.mseed import msread

MODE = Path(__file__).resolve().parents[1] / "data" / "MODE.mseed"


def _records(collection):
    """Records of a collection, in identifier order."""
    return sorted(
        (rec for stream in collection for rec in stream),
        key=lambda rec: rec.sid,
    )


# ----------------------------------------------------------------------------
# Parameters from the file, from the arguments, from the data
# ----------------------------------------------------------------------------


def test_parameters_in_the_file():
    """Commented key-value lines give the parameters."""
    text = (
        "# network: IV\n# station: MODE\n# Sampling rate [Hz]: 100\n"
        "# start: 2012-05-20T02:03:23\n# units: M/S**2\n"
        "HNE,HNN,HNZ\n0.12,0.03,-0.4\n0.13,0.02,-0.3\n"
    )
    records = _records(asciiread(text))

    assert [rec.sid for rec in records] == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ"
    ]
    rec = records[2]
    assert rec.delta == 0.01
    assert rec.time == Date("2012-05-20T02:03:23Z")
    assert rec.head.units == "M/S**2"
    np.testing.assert_array_equal(rec.data, [-0.4, -0.3])
    assert rec.head.extra["ascii"]["column"] == "HNZ"
    assert rec.head.extra["ascii"]["metadata"]["station"] == "MODE"


def test_arguments_take_precedence():
    """Arguments and file use the same names; arguments win."""
    text = "# rate: 100\n# units: M/S\nZ\n1\n2\n"
    rec = asciiread(text, rate="200 Hz", units="CM/S")[0][0]

    assert rec.delta == 0.005
    assert rec.head.units == "CM/S"


@pytest.mark.parametrize(
    "line, delta",
    [("dt = 0.02 s", 0.02), ("Sample interval: 0.5", 0.5),
     ("fs: 50", 0.02), ("SPS = 4", 0.25)],
)
def test_parameter_synonyms(line, delta):
    """Synonyms, case, blanks and units are accepted."""
    assert asciiread(f"{line}\n1\n2\n")[0][0].delta == delta


def test_numbers_only_with_arguments():
    """A file of numbers only, with the parameters as arguments."""
    collection = asciiread(
        "  1  2  3\n  4  5  6\n", rate=200, start="2026-11-05 10:00:00",
        sid=["XX.A..HNE", "XX.A..HNN", "XX.A..HNZ"],
    )
    rec = collection["XX.A..HNN"][0]

    assert rec.delta == 0.005
    assert rec.time == Date("2026-11-05T10:00:00Z")
    np.testing.assert_array_equal(rec.data, [2.0, 5.0])


def test_only_the_sampling_is_required():
    """Without start time, codes and units the records are still read."""
    rec = asciiread("1\n2\n3\n", delta=0.1)[0][0]

    assert rec.sid == "...1"
    assert rec.time is None
    assert rec.head.units is None
    with pytest.raises(ValueError, match="sampling interval is unknown"):
        asciiread("1\n2\n3\n")


def test_contradictory_sampling():
    """A rate and an interval that disagree are rejected."""
    with pytest.raises(ValueError, match="disagree"):
        asciiread("1\n2\n", rate=100, delta=0.02)
    with pytest.raises(ValueError, match="positive"):
        asciiread("1\n2\n", rate=0)
    with pytest.raises(ValueError, match="invalid sampling rate"):
        asciiread("1\n2\n", rate="fast")
    with pytest.raises(ValueError, match="invalid start"):
        asciiread("1\n2\n", rate=1, start="yesterday")


# ----------------------------------------------------------------------------
# Columns
# ----------------------------------------------------------------------------


def test_column_of_seconds():
    """A column of seconds gives the sampling interval."""
    text = (
        "Sampling rate = 200 Hz\ntempo;X;Y\n"
        "0,000;1,5;-2,25\n0,005;1,6;-2,30\n"
    )
    records = _records(asciiread(
        text, sid="XX.A..HN{column}", start="2026-11-05 10:00:00"
    ))

    assert [rec.sid for rec in records] == ["XX.A..HNX", "XX.A..HNY"]
    assert records[1].delta == 0.005
    np.testing.assert_array_equal(records[1].data, [-2.25, -2.30])


def test_column_of_dates():
    """A column of dates gives the start time and the interval."""
    text = "date\tvalue\n2026-11-05 10:00:00.000\t1.5\n" \
        "2026-11-05 10:00:00.010\t1.6\n"
    rec = asciiread(text, network="XX", station="A")[0][0]

    assert rec.time == Date("2026-11-05T10:00:00Z")
    assert rec.delta == pytest.approx(0.01)
    assert rec.sid == "XX.A..value"


def test_unnamed_column_of_dates():
    """A first column of dates is the column of times."""
    text = "2026-11-05T10:00:00 1\n2026-11-05T10:00:01 2\n"
    rec = asciiread(text)[0][0]

    assert rec.delta == 1.0
    assert rec.time == Date("2026-11-05T10:00:00Z")


def test_seconds_from_start():
    """A column of seconds counts from the start time."""
    text = "t,a\n1.0,5\n1.5,6\n"
    rec = asciiread(text, start="2026-01-01T00:00:00")[0][0]
    assert rec.time == Date("2026-01-01T00:00:01Z")

    with pytest.warns(AsciiSignalWarning, match="start time is None"):
        rec = asciiread(text)[0][0]
    assert rec.time is None
    assert rec.delta == 0.5


def test_columns_without_line_of_names():
    """Names are given in file order; None skips a column."""
    rec = asciiread(
        "0.0 9 1\n0.5 9 2\n", columns=["t", None, "Z"], network="XX",
        station="B",
    )["XX.B..Z"][0]

    assert rec.delta == 0.5
    np.testing.assert_array_equal(rec.data, [1.0, 2.0])


def test_time_column_by_name_or_index():
    """The column of times can be named or indexed."""
    text = "a,b\n0,0\n2,4\n"
    assert asciiread(text, time_column="a")[0][0].delta == 2.0
    assert asciiread(text, time_column=1)[0][0].sid == "...a"
    with pytest.raises(ValueError, match="not found"):
        asciiread(text, time_column="c")
    with pytest.raises(ValueError, match="out of range"):
        asciiread(text, time_column=5)


def test_irregular_and_conflicting_times():
    """Irregular times and a different given interval give warnings."""
    with pytest.warns(AsciiSignalWarning, match="irregular"):
        asciiread("t,a\n0,1\n1,2\n2.5,3\n")
    with pytest.warns(AsciiSignalWarning, match="given one is used"):
        rec = asciiread("t,a\n0,1\n1,2\n", delta=0.5)[0][0]
    assert rec.delta == 0.5
    with pytest.raises(ValueError, match="do not increase"):
        asciiread("t,a\n0,1\n-1,2\n")


def test_identifiers():
    """Identifiers from patterns, lists, names or codes."""
    text = "IV.MODE..HNE,b\n1,2\n3,4\n"
    records = _records(asciiread(text, rate=1, station="IV.MODE"))
    assert [rec.sid for rec in records] == ["IV.MODE..HNE", "IV.MODE..b"]

    with pytest.raises(ValueError, match="pattern with"):
        asciiread(text, rate=1, sid="XX.A..HHZ")
    with pytest.raises(ValueError, match="1 sid values for 2"):
        asciiread(text, rate=1, sid=["XX.A..HHZ"])
    assert asciiread("1\n2\n", rate=1, sid="XX.A..HHZ")[0][0].sid \
        == "XX.A..HHZ"


def test_units_per_column():
    """Units can be given per data column."""
    text = "# units: M/S**2, M/S\na,b\n1,2\n"
    records = _records(asciiread(text, rate=1))
    assert [rec.head.units for rec in records] == ["M/S**2", "M/S"]

    records = _records(asciiread("a,b\n1,2\n", rate=1, units=["M", "CM"]))
    assert [rec.head.units for rec in records] == ["M", "CM"]


def test_no_data_column():
    """A table with only times has no data."""
    with pytest.raises(ValueError, match="no numeric data column"):
        asciiread("time\n0\n1\n")


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("time", ["relative", "absolute", None])
def test_round_trip(time):
    """Written tables are read back with the same records."""
    collection = msread(MODE)
    buffer = _io.StringIO()
    asciiwrite(buffer, collection, time=time)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        again = asciiread(buffer.getvalue())
    for rec, other in zip(_records(collection), _records(again)):
        assert other.sid == rec.sid
        assert other.time == rec.time
        assert other.delta == rec.delta
        np.testing.assert_array_equal(other.data, rec.data)


def test_written_text():
    """Parameters as comment lines, identifiers as names."""
    rec = Record(data=[0.5, -1.25], delta=0.01, sid="XX.A..HNZ",
                 time="2026-01-01T00:00:00Z", units="M/S**2")
    buffer = _io.StringIO()
    asciiwrite(buffer, rec, delimiter=";")

    assert buffer.getvalue().splitlines() == [
        "# delta: 0.01",
        "# start: 2026-01-01T00:00:00.000000Z",
        "# units: M/S**2",
        "time;XX.A..HNZ",
        "0;0.5",
        "0.01;-1.25",
    ]


def test_write_incompatible_records():
    """Records of different length or start cannot share a table."""
    a = Record(data=[1.0, 2.0], delta=1.0, sid="XX.A..HHZ",
               time="2026-01-01T00:00:00Z")
    b = Record(data=[1.0], delta=1.0, sid="XX.A..HHN",
               time="2026-01-01T00:00:00Z")
    from shakelab.signals.base import StreamCollection

    collection = StreamCollection()
    collection.append(a)
    collection.append(b)
    with pytest.raises(ValueError, match="separate files"):
        asciiwrite(_io.StringIO(), collection)
    with pytest.raises(ValueError, match="time must be"):
        asciiwrite(_io.StringIO(), a, time="iso")
    a.head.time = None
    with pytest.raises(ValueError, match="need the start time"):
        asciiwrite(_io.StringIO(), a, time="absolute")
    with pytest.raises(TypeError):
        asciiwrite(_io.StringIO(), [a])


def test_io_reader_and_writer(tmp_path):
    """``io`` reads and writes the format by name or ``.csv``."""
    path = tmp_path / "mode.csv"
    writer(msread(MODE), str(path), time="absolute")

    collection = reader(str(path))
    assert sorted(st.sid for st in collection) == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ"
    ]
    text = tmp_path / "x.txt"
    text.write_text("1\n2\n")
    rec = reader(str(text), format="ascii", rate=10,
                 start="2026-01-01T00:00:00")[0][0]
    assert rec.delta == 0.1
