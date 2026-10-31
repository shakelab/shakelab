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
Unit tests for :mod:`shakelab.libutils.json`.

The reader is checked against standard JSON documents, escapes,
numbers and the NaN/Infinity extension, and must reject malformed
input with a positioned error. The writer is checked for valid
escaping, float precision, NumPy support, layout options and the
round trip through the reader, plus the optional features (comments,
float formats, sorted keys). The standard ``json`` module is used only
as an independent reference.
"""

from __future__ import annotations

import io
import json as std_json
import math

import numpy as np
import pytest

from shakelab.libutils.json import (
    JsonParseError,
    read,
    read_str,
    write,
    write_str,
)


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


SAMPLE = """
{
  "type": "Feature",
  "geometry": {"type": "Point", "coordinates": [13.25, 46.06, -1e3]},
  "properties": {
    "name": "Udine",
    "count": 0,
    "ok": true,
    "missing": null,
    "list": [],
    "nested": {"a": [1, [2, {"b": false}]]}
  }
}
"""


def test_read_matches_standard_library():
    """A standard document parses as with the reference parser."""
    assert read_str(SAMPLE) == std_json.loads(SAMPLE)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("0", 0),
        ("-12", -12),
        ("3.5", 3.5),
        ("1e3", 1000.0),
        ("-2.5E-2", -0.025),
        ("12345678901234567890", 12345678901234567890),
    ],
)
def test_read_numbers(text, expected):
    """Integers stay int, other numbers become float."""
    value = read_str(text)
    assert value == expected
    assert type(value) is type(expected)


def test_read_nan_and_infinity():
    """NaN and infinities are accepted as an extension."""
    value = read_str("[NaN, Infinity, -Infinity]")
    assert math.isnan(value[0])
    assert value[1:] == [math.inf, -math.inf]


def test_read_string_escapes():
    """All JSON escapes are decoded, surrogate pairs included."""
    text = r'"q\" b\\ s\/ \b\f\n\r\t è 🌍"'
    assert read_str(text) == 'q" b\\ s/ \b\f\n\r\t è 🌍'


def test_read_unicode_text_and_duplicate_keys():
    """Non-ASCII text is kept and the last duplicate key wins."""
    assert read_str('{"città": "Gemona", "a": 1, "a": 2}') == {
        "città": "Gemona",
        "a": 2,
    }


@pytest.mark.parametrize(
    "text",
    [
        "",
        "[1, 2",
        "[1, 2,]",
        '{"a": 1,}',
        "{'a': 1}",
        '{"a" 1}',
        "01",
        "1.",
        ".5",
        "+1",
        "tru",
        '"open',
        '"tab\there"',
        r'"\x"',
        r'"\u12"',
        "[1] [2]",
    ],
)
def test_read_rejects_malformed(text):
    """Malformed documents raise JsonParseError."""
    with pytest.raises(JsonParseError):
        read_str(text)


def test_parse_error_position():
    """Errors report line and column."""
    with pytest.raises(JsonParseError) as info:
        read_str('{\n  "a": 1,\n  "b" 2\n}')

    err = info.value
    assert (err.lineno, err.colno) == (3, 7)
    assert isinstance(err, ValueError)


def test_read_from_path_and_stream(tmp_path):
    """read accepts paths and open text streams."""
    path = tmp_path / "data.json"
    path.write_text('{"a": [1, 2]}', encoding="utf-8")

    assert read(path) == {"a": [1, 2]}
    assert read(str(path)) == {"a": [1, 2]}
    assert read(io.StringIO("[true]")) == [True]


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


def test_write_escapes_strings():
    """Strings with quotes, backslashes and controls stay valid JSON."""
    data = {'k"\\': 'a "b" \\ \n\t\x01 è'}
    text = write_str(data)
    assert std_json.loads(text) == data
    assert "è" in text


def test_write_layout():
    """Objects and arrays are indented; numeric arrays stay inline."""
    data = {"a": [1, 2.5], "b": ["x", {"c": None}], "e": {}, "f": []}
    assert write_str(data) == (
        "{\n"
        '  "a": [1, 2.5],\n'
        '  "b": [\n'
        '    "x",\n'
        "    {\n"
        '      "c": null\n'
        "    }\n"
        "  ],\n"
        '  "e": {},\n'
        '  "f": []\n'
        "}"
    )


def test_write_single_line_and_no_inline():
    """indent=None writes one line; inline_numbers=False expands."""
    data = {"a": [1, [2, 3]], "b": True}
    assert write_str(data, indent=None) == '{"a": [1, [2, 3]], "b": true}'
    assert write_str([1, 2], indent=1, inline_numbers=False) == (
        "[\n 1,\n 2\n]"
    )


def test_write_multipolygon_coordinates():
    """Deeply nested coordinates are written point by point."""
    coords = [[[[12.0, 45.0], [13.0, 45.0], [12.0, 46.0]]]]
    text = write_str({"coordinates": coords})
    assert "[12.0, 45.0]" in text
    assert read_str(text) == {"coordinates": coords}


def test_write_float_precision():
    """Floats are written with fixed decimals, without negative zero."""
    text = write_str([1.23456, 2.0, -0.0001, 7], float_precision=3)
    assert text == "[1.235, 2.000, 0.000, 7]"


def test_write_nan_and_infinity():
    """Non-finite floats use the NaN/Infinity extension or raise."""
    assert write_str([math.nan, math.inf, -math.inf]) == (
        "[NaN, Infinity, -Infinity]"
    )
    with pytest.raises(ValueError):
        write_str([math.nan], allow_nan=False)


def test_write_numpy_tuples_and_keys():
    """NumPy objects, tuples and non-string keys are converted."""
    data = {
        "arr": np.array([[1.5, 2.0]]),
        "int": np.int64(3),
        "flag": np.bool_(True),
        "tup": (1, "a"),
        1: "one",
        None: "none",
    }
    assert std_json.loads(write_str(data)) == {
        "arr": [[1.5, 2.0]],
        "int": 3,
        "flag": True,
        "tup": [1, "a"],
        "1": "one",
        "null": "none",
    }


def test_write_rejects_unsupported_types():
    """Objects that cannot be represented raise TypeError."""
    with pytest.raises(TypeError):
        write_str({"a": object()})
    with pytest.raises(TypeError):
        write_str({(1, 2): "tuple key"})


def test_write_read_round_trip(tmp_path):
    """Files written by write are read back unchanged."""
    data = std_json.loads(SAMPLE)
    path = tmp_path / "out.json"
    write(data, path)

    assert path.read_text(encoding="utf-8").endswith("}\n")
    assert read(path) == data

    stream = io.StringIO()
    write(data, stream, indent=None)
    assert read_str(stream.getvalue()) == data


# ----------------------------------------------------------------------------
# Optional features
# ----------------------------------------------------------------------------


def test_read_comments():
    """Line comments are skipped only when enabled and outside strings."""
    text = (
        "# configuration\n"
        "{\n"
        '  "url": "http://host/#x",  // not a comment inside strings\n'
        "  # magnitude threshold\n"
        '  "mmin": 3.0\n'
        "} // end"
    )
    assert read_str(text, allow_comments=True) == {
        "url": "http://host/#x",
        "mmin": 3.0,
    }
    with pytest.raises(JsonParseError):
        read_str(text)


def test_read_comments_from_file(tmp_path):
    """read forwards allow_comments."""
    path = tmp_path / "config.json"
    path.write_text("[1, // one\n 2]", encoding="utf-8")
    assert read(path, allow_comments=True) == [1, 2]


def test_write_float_format_significant_digits():
    """Significant-digit formats keep small values and float type."""
    text = write_str(
        [1.2e-7, 2.0, 1234567.0, -0.0, 3],
        float_format=".6g",
    )
    assert text == "[1.2e-07, 2.0, 1.23457e+06, 0.0, 3]"

    values = read_str(text)
    assert all(type(v) is float for v in values[:4])


def test_write_float_options_validation():
    """Float options are exclusive and formats are validated."""
    with pytest.raises(ValueError):
        write_str([1.0], float_precision=2, float_format=".3g")
    with pytest.raises(ValueError):
        write_str([1.0], float_format="not a format")


def test_write_zero_decimals_keeps_float():
    """A float written with no decimals still reads back as float."""
    assert write_str([2.4], float_precision=0) == "[2.0]"


def test_write_sort_keys(tmp_path):
    """Keys are written in sorted order at every level."""
    data = {"b": 1, "a": {"d": 2, "c": 3}}
    assert write_str(data, indent=None, sort_keys=True) == (
        '{"a": {"c": 3, "d": 2}, "b": 1}'
    )
    assert write_str(data, indent=None) == (
        '{"b": 1, "a": {"d": 2, "c": 3}}'
    )

    path = tmp_path / "sorted.json"
    write(data, path, indent=None, sort_keys=True)
    assert path.read_text(encoding="utf-8") == (
        '{"a": {"c": 3, "d": 2}, "b": 1}\n'
    )
