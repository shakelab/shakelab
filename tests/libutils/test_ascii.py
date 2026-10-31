# ****************************************************************************
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
# ****************************************************************************
"""
Unit tests for :mod:`shakelab.libutils.ascii`.

The tests cover row and column editing, broadcasting and validation
of column data, reading files with comments, blank lines, skipped
lines, unnamed columns and per-column types, and the write/read round
trip.
"""

from __future__ import annotations

import pytest

from shakelab.libutils.ascii import AsciiTable


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
