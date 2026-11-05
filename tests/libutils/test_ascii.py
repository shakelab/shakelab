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
Unit tests for :mod:`shakelab.libutils.ascii`.

The tests cover row and column editing, broadcasting and validation
of column data, reading files with comments, blank lines, skipped
lines, unnamed columns and per-column types, the write/read round
trip, quoted fields, the selection helpers (``add_row``,
``to_array``, ``unique``, ``filter``, ``select``, ``sort``), the text
sources and the numeric tables of any layout (``read_columns``,
``write_columns``).
"""

from __future__ import annotations

import io

import numpy as np
import pytest

from shakelab.libutils.ascii import (
    AsciiTable,
    AsciiWarning,
    TextColumns,
    read_columns,
    read_text,
    write_columns,
    write_text,
)


# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------


@pytest.fixture
def table():
    """Two-row table with columns a and b."""
    tab = AsciiTable(["a", "b"])
    tab.add_element([1, "x"])
    tab.add_element([2, "y"])
    return tab


# ----------------------------------------------------------------------------
# Construction and container protocol
# ----------------------------------------------------------------------------


def test_header_is_copied():
    """The header list passed to the constructor is not shared."""
    header = ["a"]
    tab = AsciiTable(header)
    tab.add_key("b")
    assert header == ["a"]
    assert AsciiTable().header == []


def test_size_len_and_iteration(table):
    """size, len and iteration describe rows and columns."""
    assert table.size == (2, 2)
    assert len(table) == 2
    assert [row["a"] for row in table] == [1, 2]

    # Iteration is stateless: a second pass sees all rows again.
    assert [row["b"] for row in table] == ["x", "y"]


def test_getitem_and_extract(table):
    """Columns are extracted by name, optionally cast."""
    assert table["b"] == ["x", "y"]
    assert table.extract("a", dtype=float) == [1.0, 2.0]
    with pytest.raises(KeyError):
        table.extract("missing")


def test_setitem(table):
    """Scalars, strings included, are broadcast; sequences are checked."""
    table["b"] = "z"
    assert table["b"] == ["z", "z"]

    table["c"] = [10, 20]
    assert table.header == ["a", "b", "c"]
    assert table["c"] == [10, 20]

    with pytest.raises(ValueError):
        table["c"] = [1, 2, 3]


def test_add_element_checks_length(table):
    """Rows must have one value per column."""
    with pytest.raises(ValueError):
        table.add_element([3])


# ----------------------------------------------------------------------------
# Column editing
# ----------------------------------------------------------------------------


def test_add_key(table):
    """Columns are inserted at the given position."""
    table.add_key("id", 0, [7, 8])
    table.add_key("flag", data=True)
    table.add_key("note")

    assert table.header == ["id", "a", "b", "flag", "note"]
    assert table["id"] == [7, 8]
    assert table["flag"] == [True, True]
    assert table["note"] == [None, None]


def test_add_key_errors(table):
    """Duplicate keys and wrong data lengths raise ValueError."""
    with pytest.raises(ValueError):
        table.add_key("a")
    with pytest.raises(ValueError):
        table.add_key("c", data=[1])
    assert table.header == ["a", "b"]


def test_remove_and_rename_key(table):
    """Columns can be removed and renamed in place."""
    table.rename_key("a", "n")
    assert table.header == ["n", "b"]
    assert table["n"] == [1, 2]

    table.remove_key("b")
    assert table.header == ["n"]
    assert table.db == [{"n": 1}, {"n": 2}]

    with pytest.raises(KeyError):
        table.remove_key("b")

    table.add_key("m")
    with pytest.raises(ValueError):
        table.rename_key("n", "m")


def test_recast(table):
    """Column values are converted in place."""
    table.recast("a", str)
    assert table["a"] == ["1", "2"]


def test_concatenate(table):
    """Rows of a table with the same header are copied."""
    other = AsciiTable(["a", "b"])
    other.add_element([3, "z"])
    table.concatenate(other)

    assert table["a"] == [1, 2, 3]
    table.db[-1]["a"] = 99
    assert other["a"] == [3]

    with pytest.raises(ValueError):
        table.concatenate(AsciiTable(["a"]))


# ----------------------------------------------------------------------------
# Reading and writing
# ----------------------------------------------------------------------------


def test_read_with_comments_and_blank_lines(tmp_path):
    """Comments, blank lines and surrounding blanks are ignored."""
    path = tmp_path / "data.csv"
    path.write_text(
        "# catalogue\n"
        "\n"
        "Id, Magnitude\n"
        "  # comment\n"
        "1, 4.5\n"
        "\n"
        "2, \n"
    )

    tab = AsciiTable()
    tab.read(path, dtype=[int, float], empty=-1.0)

    assert tab.header == ["Id", "Magnitude"]
    assert tab.db == [
        {"Id": 1, "Magnitude": 4.5},
        {"Id": 2, "Magnitude": -1.0},
    ]


def test_read_skipline_header_and_unnamed_columns(tmp_path):
    """Initial lines can be skipped and unnamed columns are dropped."""
    path = tmp_path / "data.txt"
    path.write_text("title line\n1;2;3\n4;5;6\n")

    tab = AsciiTable()
    tab.read(
        path,
        header=["a", "", "c"],
        dtype=int,
        delimiter=";",
        skipline=1,
    )

    assert tab.header == ["a", "c"]
    assert tab.db == [{"a": 1, "c": 3}, {"a": 4, "c": 6}]


def test_read_errors(tmp_path):
    """Short lines, bad dtype lists and missing files raise."""
    path = tmp_path / "data.csv"
    path.write_text("a,b\n1\n")

    with pytest.raises(ValueError, match="line 2"):
        AsciiTable().read(path)

    path.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError):
        AsciiTable().read(path, dtype=[int])

    with pytest.raises(FileNotFoundError):
        AsciiTable().read(tmp_path / "missing.csv")


def test_write_read_round_trip(table, tmp_path):
    """A written table is read back identically."""
    table.add_key("c", data=[None, 3.5])
    path = tmp_path / "out.csv"
    table.write(path)

    assert path.read_text() == "a,b,c\n1,x,\n2,y,3.5\n"

    tab = AsciiTable()
    tab.read(path, dtype=[int, str, float])
    assert tab.header == table.header
    assert tab.db == table.db


def test_write_without_header(table, tmp_path):
    """The header line can be omitted."""
    path = tmp_path / "out.txt"
    table.write(path, header=False, delimiter=" ")
    assert path.read_text() == "1 x\n2 y\n"


# ----------------------------------------------------------------------------
# Quoted fields
# ----------------------------------------------------------------------------


def test_read_quoted_fields(tmp_path):
    """Quoted fields may contain delimiters, blanks and quotes."""
    path = tmp_path / "data.csv"
    path.write_text(
        'Id,Name,Note\n'
        '1, "Udine, Italy" ," keep "\n'
        '2,Gemona,"say ""hi"""\n'
        '3,5\'3",""\n'
    )

    tab = AsciiTable()
    tab.read(path, dtype=[int, str, str])

    assert tab["Name"] == ["Udine, Italy", "Gemona", "5'3\""]
    assert tab["Note"] == [" keep ", 'say "hi"', None]


def test_read_quote_errors(tmp_path):
    """Unterminated quotes and text after a closing quote raise."""
    path = tmp_path / "data.csv"
    path.write_text('a,b\n1,"open\n')
    with pytest.raises(ValueError, match="line 2: unterminated"):
        AsciiTable().read(path)

    path.write_text('a,b\n1,"x"y\n')
    with pytest.raises(ValueError, match="line 2: unexpected"):
        AsciiTable().read(path)


def test_read_quoting_disabled(tmp_path):
    """With quotechar=None quotes are ordinary characters."""
    path = tmp_path / "data.csv"
    path.write_text('a,b\n"x",y\n')

    tab = AsciiTable()
    tab.read(path, quotechar=None)
    assert tab["a"] == ['"x"']


def test_write_quotes_only_when_needed(tmp_path):
    """Fields are quoted only if they would not read back unchanged."""
    tab = AsciiTable(["name", "note"])
    tab.add_element(["Udine, Italy", 'say "hi"'])
    tab.add_element(["Gemona", " padded "])

    path = tmp_path / "out.csv"
    tab.write(path)
    assert path.read_text() == (
        'name,note\n'
        '"Udine, Italy","say ""hi"""\n'
        'Gemona," padded "\n'
    )

    back = AsciiTable()
    back.read(path)
    assert back.db == tab.db


# ----------------------------------------------------------------------------
# Selection helpers
# ----------------------------------------------------------------------------


@pytest.fixture
def catalogue():
    """Small catalogue with an empty magnitude."""
    tab = AsciiTable(["Id", "Agency", "Mag"])
    tab.add_element([1, "INGV", 3.2])
    tab.add_element([2, "OGS", None])
    tab.add_element([3, "INGV", 4.7])
    tab.add_element([4, "OGS", 2.1])
    return tab


def test_add_row(catalogue):
    """Rows given as mappings fill missing columns with None."""
    catalogue.add_row({"Id": 5, "Mag": 1.0})
    assert catalogue.db[-1] == {"Id": 5, "Agency": None, "Mag": 1.0}

    with pytest.raises(KeyError):
        catalogue.add_row({"Depth": 10.0})


def test_to_array(catalogue):
    """Columns become arrays with NaN for empty values."""
    mag = catalogue.to_array("Mag")
    assert mag.dtype == float
    np.testing.assert_array_equal(mag, [3.2, np.nan, 4.7, 2.1])

    ids = catalogue.to_array("Id", dtype=int)
    np.testing.assert_array_equal(ids, [1, 2, 3, 4])


def test_unique(catalogue):
    """Distinct values keep the order of first appearance."""
    assert catalogue.unique("Agency") == ["INGV", "OGS"]


def test_filter(catalogue):
    """Filtering returns a new table with copied rows."""
    strong = catalogue.filter(
        lambda r: r["Mag"] is not None and r["Mag"] >= 3.0
    )
    assert strong.header == catalogue.header
    assert strong["Id"] == [1, 3]

    strong.db[0]["Mag"] = 0.0
    assert catalogue.db[0]["Mag"] == 3.2


def test_select(catalogue):
    """Selecting columns returns a new table in the given order."""
    sub = catalogue.select(["Mag", "Id"])
    assert sub.header == ["Mag", "Id"]
    assert sub.db[0] == {"Mag": 3.2, "Id": 1}

    with pytest.raises(KeyError):
        catalogue.select(["Depth"])


def test_sort(catalogue):
    """Rows are sorted in place with empty values last."""
    catalogue.sort("Mag")
    assert catalogue["Id"] == [4, 1, 3, 2]

    catalogue.sort("Mag", reverse=True)
    assert catalogue["Id"] == [3, 1, 4, 2]


# ----------------------------------------------------------------------------
# Text sources
# ----------------------------------------------------------------------------


def test_read_text_sources(tmp_path):
    """Paths, content, bytes and file objects give the same text."""
    path = tmp_path / "x.txt"
    path.write_text("a\nb\n", encoding="utf-8")

    for source in (path, str(path), "a\nb\n", b"a\nb\n",
                   io.StringIO("a\nb\n"), io.BytesIO(b"a\nb\n")):
        assert read_text(source) == "a\nb\n"
    with pytest.raises(TypeError):
        read_text(42)


def test_write_text(tmp_path):
    """Files are not overwritten unless requested."""
    path = tmp_path / "x.txt"
    write_text(path, "one\n")
    with pytest.raises(FileExistsError):
        write_text(path, "two\n")
    write_text(path, "two\n", owrite=True)
    assert path.read_text() == "two\n"

    buffer = io.StringIO()
    write_text(buffer, "three\n")
    assert buffer.getvalue() == "three\n"
    with pytest.raises(TypeError):
        write_text(42, "x")


# ----------------------------------------------------------------------------
# Numeric tables of any layout
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text, delimiter, decimal",
    [
        ("a,b\n1.5,2\n3,4.25\n", ",", "."),
        ("a;b\n1,5;2\n3;4,25\n", ";", ","),
        ("a\tb\n1.5\t2\n3\t4.25\n", "\t", "."),
        ("a|b\n1.5|2\n3|4.25\n", "|", "."),
        ("a b\n 1.5   2\n 3   4.25\n", None, "."),
        ("a b\n1,5 2\n3 4,25\n", None, ","),
    ],
)
def test_read_columns_layouts(text, delimiter, decimal):
    """Delimiter and decimal separator are recognized."""
    table = read_columns(text)

    assert table.delimiter == delimiter
    assert table.decimal == decimal
    assert table.names == ["a", "b"]
    np.testing.assert_array_equal(table["a"], [1.5, 3.0])
    np.testing.assert_array_equal(table[1], [2.0, 4.25])
    assert table.nrows == 2 and len(table) == 2


def test_read_columns_metadata():
    """Key-value lines before the table, commented or not."""
    text = (
        "\ufeffExported by logger v2\n"
        "# station: IV.MODE\n"
        "Sampling rate = 200 Hz\n"
        "# start: 2012-05-20T02:03:23\n"
        "\n"
        "time,X\n"
        "0,1\n"
        "0.005,2\n"
    )
    table = read_columns(text)

    assert table.metadata == {
        "station": "IV.MODE",
        "Sampling rate": "200 Hz",
        "start": "2012-05-20T02:03:23",
    }
    assert table.names == ["time", "X"]


def test_read_columns_commented_names():
    """A commented line of names (numpy.savetxt) is used."""
    table = read_columns("# t acc\n0.0 1e-3\n0.5 2e-3\n")

    assert table.names == ["t", "acc"]
    np.testing.assert_array_equal(table["acc"], [1e-3, 2e-3])


def test_read_columns_without_names():
    """Without a line of names, names are None or given."""
    text = "1 2 3\n4 5 6\n"
    assert read_columns(text).names is None

    table = read_columns(text, names=["x", "y", "z"])
    np.testing.assert_array_equal(table["z"], [3.0, 6.0])
    with pytest.raises(ValueError, match="2 names for 3 columns"):
        read_columns(text, names=["x", "y"])


def test_read_columns_names_mismatch_warns():
    """A line of names with a different number of fields is ignored."""
    with pytest.warns(AsciiWarning, match="names ignored"):
        table = read_columns("a,b,c\n1,2\n3,4\n")
    assert table.names is None


def test_read_columns_dates():
    """Dates (also with a blank) are kept as text."""
    table = read_columns(
        "date\tvalue\n2026-11-05 10:00:00.000\t1.5\n"
        "2026-11-05T10:00:00.010Z\t1.6\n"
    )
    assert table["date"].dtype.kind == "U"
    assert list(table["date"]) == [
        "2026-11-05T10:00:00.000", "2026-11-05T10:00:00.010"
    ]


def test_read_columns_empty_and_invalid_values():
    """Empty values are NaN; other text in numbers gives a warning."""
    with pytest.warns(AsciiWarning, match="column b: 1 values"):
        table = read_columns("a,b\n1,2\n3,\n5,x\n")
    np.testing.assert_array_equal(table["a"], [1.0, 3.0, 5.0])
    assert np.isnan(table["b"][1]) and np.isnan(table["b"][2])


def test_read_columns_errors():
    """Missing tables and rows of a different length are rejected."""
    with pytest.raises(ValueError, match="no table"):
        read_columns("only: metadata\n")
    with pytest.raises(ValueError, match="line 3: 1 fields, expected 2"):
        read_columns("1,2\n3,4\n5\n")
    with pytest.raises(ValueError, match="decimal"):
        read_columns("1\n", decimal=";")


def test_read_columns_skipline_and_comments():
    """Skipped lines and comment lines inside the table."""
    table = read_columns(
        "garbage 1 2 3\n% a b\n1 2\n% note\n3 4\n", skipline=1,
        comment="%",
    )
    assert table.names == ["a", "b"]
    np.testing.assert_array_equal(table["a"], [1.0, 3.0])


def test_read_columns_large_table_fast_path():
    """Tables of numbers only are read in one pass."""
    data = np.arange(30000.0).reshape(-1, 3)
    text = "\n".join(",".join(f"{v:g}" for v in row) for row in data)
    table = read_columns(text + "\n")
    np.testing.assert_array_equal(np.column_stack(table.columns), data)


def test_write_columns_round_trip(tmp_path):
    """Written tables are read back with names and metadata."""
    path = tmp_path / "t.csv"
    write_columns(
        path, [[0.0, 0.5], np.array(["2026-01-01T00:00:00", "x"]),
               [1.25, -3.0]],
        names=["t", "label", "v"], metadata={"units": "M/S"},
    )
    assert path.read_text().splitlines()[:2] == ["# units: M/S", "t,label,v"]

    table = read_columns(path)
    assert table.metadata == {"units": "M/S"}
    np.testing.assert_array_equal(table["v"], [1.25, -3.0])
    with pytest.raises(FileExistsError):
        write_columns(path, [[1.0]])
    with pytest.raises(ValueError, match="same length"):
        write_columns(io.StringIO(), [[1.0], [1.0, 2.0]])
    with pytest.raises(ValueError, match="names"):
        write_columns(io.StringIO(), [[1.0]], names=["a", "b"])


def test_text_columns_access():
    """Columns by name or index."""
    table = TextColumns([[1.0], [2.0]], names=["a", "b"])
    assert table["b"][0] == 2.0 and table[0][0] == 1.0
    with pytest.raises(KeyError):
        table["c"]
    with pytest.raises(KeyError):
        TextColumns([[1.0]])["a"]
    with pytest.raises(ValueError):
        TextColumns([[1.0]], names=["a", "b"])
