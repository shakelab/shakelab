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
Minimal JSON reader and writer, independent of the standard ``json``
module.

The module provides a strict recursive-descent parser (:func:`read`,
:func:`read_str`) and a writer (:func:`write`, :func:`write_str`) with
controlled float precision and compact output of numeric arrays, which
keeps GeoJSON-like ``coordinates`` readable. Being self-contained, it
is the place where ShakeLab-specific JSON features can be added.

Conventions
-----------
- Input must be standard JSON (RFC 8259), with one extension: the
  tokens ``NaN``, ``Infinity`` and ``-Infinity`` are accepted, as in
  the output of most scientific tools. Line comments starting with
  ``//`` or ``#`` can be enabled with ``allow_comments``.
- Integers are returned as ``int``, other numbers as ``float``;
  duplicate object keys keep the last value.
- Strings are written as UTF-8 without escaping non-ASCII characters.
- Tuples and NumPy arrays are written as arrays, NumPy scalars as the
  corresponding Python scalars. Non-string keys (numbers, booleans,
  ``None``) are converted to strings.
- Arrays whose items are all numbers are written on a single line
  (e.g. ``[12.0, 45.0, 0.0]``) unless ``inline_numbers`` is False.
- Floats are written with full precision by default, with a fixed
  number of decimals (``float_precision``) or with any format
  specification (``float_format``, e.g. ``".6g"``); a float is always
  written so that it reads back as a float.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Optional, TextIO, Union

import numpy as np

__all__ = [
    "JsonParseError",
    "read",
    "read_str",
    "write",
    "write_str",
]

FileLike = Union[str, Path, TextIO]


class JsonParseError(ValueError):
    """
    Error raised for malformed JSON input.

    Attributes
    ----------
    pos : int
        Character offset of the error.
    lineno, colno : int
        Line and column of the error (1-based).
    """

    def __init__(self, msg: str, text: str, pos: int) -> None:
        self.pos = pos
        self.lineno = text.count("\n", 0, pos) + 1
        self.colno = pos - text.rfind("\n", 0, pos)
        super().__init__(
            f"{msg}: line {self.lineno} column {self.colno} (char {pos})"
        )


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

_WHITESPACE = " \t\n\r"
_NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?")
_ESCAPES = {
    '"': '"', "\\": "\\", "/": "/", "b": "\b",
    "f": "\f", "n": "\n", "r": "\r", "t": "\t",
}
_LITERALS = (
    ("true", True),
    ("false", False),
    ("null", None),
    ("NaN", math.nan),
    ("Infinity", math.inf),
    ("-Infinity", -math.inf),
)


class _Parser:
    """Recursive-descent JSON parser over a string."""

    def __init__(self, text: str, allow_comments: bool = False) -> None:
        self.text = text
        self.pos = 0
        self.allow_comments = allow_comments

    def error(self, msg: str, pos: int | None = None) -> JsonParseError:
        return JsonParseError(
            msg, self.text, self.pos if pos is None else pos
        )

    def skip_whitespace(self) -> None:
        """Skip blanks and, if enabled, line comments."""
        text, pos = self.text, self.pos
        while pos < len(text):
            if text[pos] in _WHITESPACE:
                pos += 1
            elif self.allow_comments and (
                text[pos] == "#" or text.startswith("//", pos)
            ):
                end = text.find("\n", pos)
                pos = len(text) if end < 0 else end + 1
            else:
                break
        self.pos = pos

    def parse(self) -> Any:
        value = self.parse_value()
        self.skip_whitespace()
        if self.pos != len(self.text):
            raise self.error("Extra data")
        return value

    def parse_value(self) -> Any:
        self.skip_whitespace()
        if self.pos >= len(self.text):
            raise self.error("Unexpected end of input")

        ch = self.text[self.pos]
        if ch == "{":
            return self.parse_object()
        if ch == "[":
            return self.parse_array()
        if ch == '"':
            return self.parse_string()

        for word, value in _LITERALS:
            if self.text.startswith(word, self.pos):
                self.pos += len(word)
                return value

        if ch == "-" or ch.isdigit():
            return self.parse_number()

        raise self.error(f"Unexpected character {ch!r}")

    def parse_number(self) -> int | float:
        match = _NUMBER.match(self.text, self.pos)
        if match is None:
            raise self.error("Invalid number")
        self.pos = match.end()
        token = match.group(0)
        if match.group(1) or match.group(2):
            return float(token)
        return int(token)

    def parse_string(self) -> str:
        text = self.text
        start = self.pos
        pos = start + 1
        parts: list[str] = []

        while True:
            chunk = pos
            while pos < len(text) and text[pos] not in '"\\' \
                    and text[pos] >= " ":
                pos += 1
            parts.append(text[chunk:pos])

            if pos >= len(text):
                raise self.error("Unterminated string", start)

            ch = text[pos]
            if ch == '"':
                self.pos = pos + 1
                return "".join(parts)

            if ch != "\\":
                raise self.error("Invalid control character", pos)

            esc = text[pos + 1:pos + 2]
            if esc in _ESCAPES:
                parts.append(_ESCAPES[esc])
                pos += 2
            elif esc == "u":
                char, pos = self.parse_unicode(pos)
                parts.append(char)
            else:
                raise self.error("Invalid escape", pos)

    def parse_unicode(self, pos: int) -> tuple[str, int]:
        """Decode a ``\\uXXXX`` escape, joining surrogate pairs."""
        code = self.hex4(pos)
        pos += 6
        if 0xD800 <= code <= 0xDBFF \
                and self.text.startswith("\\u", pos):
            low = self.hex4(pos)
            if 0xDC00 <= low <= 0xDFFF:
                code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00)
                pos += 6
        return chr(code), pos

    def hex4(self, pos: int) -> int:
        digits = self.text[pos + 2:pos + 6]
        if len(digits) != 4 or not all(
            c in "0123456789abcdefABCDEF" for c in digits
        ):
            raise self.error("Invalid \\u escape", pos)
        return int(digits, 16)

    def parse_array(self) -> list:
        self.pos += 1
        items: list = []
        self.skip_whitespace()
        if self.text.startswith("]", self.pos):
            self.pos += 1
            return items

        while True:
            items.append(self.parse_value())
            self.skip_whitespace()
            ch = self.text[self.pos:self.pos + 1]
            self.pos += 1
            if ch == "]":
                return items
            if ch != ",":
                raise self.error("Expected ',' or ']'", self.pos - 1)

    def parse_object(self) -> dict:
        self.pos += 1
        obj: dict = {}
        self.skip_whitespace()
        if self.text.startswith("}", self.pos):
            self.pos += 1
            return obj

        while True:
            self.skip_whitespace()
            if not self.text.startswith('"', self.pos):
                raise self.error("Expected property name")
            key = self.parse_string()

            self.skip_whitespace()
            if not self.text.startswith(":", self.pos):
                raise self.error("Expected ':'")
            self.pos += 1

            obj[key] = self.parse_value()

            self.skip_whitespace()
            ch = self.text[self.pos:self.pos + 1]
            self.pos += 1
            if ch == "}":
                return obj
            if ch != ",":
                raise self.error("Expected ',' or '}'", self.pos - 1)


def read_str(text: str, allow_comments: bool = False) -> Any:
    """
    Parse JSON content from a string.

    Parameters
    ----------
    text
        JSON document.
    allow_comments
        Ignore line comments starting with ``//`` or ``#`` outside
        strings.

    Returns
    -------
    Any
        Parsed Python object.

    Raises
    ------
    JsonParseError
        If the document is malformed.
    """
    return _Parser(text, allow_comments).parse()


def read(file: FileLike, allow_comments: bool = False) -> Any:
    """
    Parse JSON content from a file.

    Parameters
    ----------
    file
        Path to a UTF-8 file, or an open readable text stream.
    allow_comments
        Ignore line comments starting with ``//`` or ``#`` outside
        strings.

    Returns
    -------
    Any
        Parsed Python object.

    Raises
    ------
    JsonParseError
        If the document is malformed.
    """
    if isinstance(file, (str, Path)):
        with open(file, "r", encoding="utf-8") as fobj:
            return read_str(fobj.read(), allow_comments)
    return read_str(file.read(), allow_comments)


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------

_STRING_ESCAPES = {
    '"': '\\"', "\\": "\\\\", "\n": "\\n", "\r": "\\r",
    "\t": "\\t", "\b": "\\b", "\f": "\\f",
}


class _Writer:
    """Recursive JSON encoder with formatting options."""

    def __init__(
        self,
        indent: Optional[int],
        float_precision: Optional[int],
        float_format: Optional[str],
        inline_numbers: bool,
        allow_nan: bool,
        sort_keys: bool,
    ) -> None:
        if float_precision is not None and float_format is not None:
            raise ValueError(
                "float_precision and float_format are mutually exclusive"
            )
        if float_precision is not None:
            float_format = f".{int(float_precision)}f"
        if float_format is not None:
            format(1.0, float_format)

        self.indent = indent if indent else None
        self.float_format = float_format
        self.inline_numbers = inline_numbers
        self.allow_nan = allow_nan
        self.sort_keys = sort_keys

    def encode(self, obj: Any, level: int = 0) -> str:
        obj = _to_python(obj)

        if obj is None:
            return "null"
        if obj is True:
            return "true"
        if obj is False:
            return "false"
        if isinstance(obj, int):
            return str(int(obj))
        if isinstance(obj, float):
            return self.encode_float(obj)
        if isinstance(obj, str):
            return _encode_string(obj)
        if isinstance(obj, dict):
            return self.encode_object(obj, level)
        if isinstance(obj, (list, tuple)):
            return self.encode_array(obj, level)

        raise TypeError(
            f"Object of type {type(obj).__name__} is not JSON serializable"
        )

    def encode_float(self, x: float) -> str:
        if math.isnan(x) or math.isinf(x):
            if not self.allow_nan:
                raise ValueError(f"Out of range float value: {x!r}")
            if math.isnan(x):
                return "NaN"
            return "Infinity" if x > 0 else "-Infinity"

        if self.float_format is None:
            return repr(float(x))

        text = format(float(x), self.float_format)
        if float(text) == 0.0:
            text = format(0.0, self.float_format)
        if not any(c in text for c in ".eE"):
            text += ".0"
        return text

    def encode_object(self, obj: dict, level: int) -> str:
        if not obj:
            return "{}"
        pairs = [(_key_to_str(k), v) for k, v in obj.items()]
        if self.sort_keys:
            pairs.sort(key=lambda kv: kv[0])
        items = [
            f"{_encode_string(k)}: {self.encode(v, level + 1)}"
            for k, v in pairs
        ]
        return self.wrap("{", items, "}", level)

    def encode_array(self, obj: list | tuple, level: int) -> str:
        if not obj:
            return "[]"
        items = [self.encode(v, level + 1) for v in obj]
        if self.inline_numbers and all(_is_number(v) for v in obj):
            return "[" + ", ".join(items) + "]"
        return self.wrap("[", items, "]", level)

    def wrap(self, open_: str, items: list[str], close: str,
             level: int) -> str:
        if self.indent is None:
            return open_ + ", ".join(items) + close
        inner = " " * (self.indent * (level + 1))
        outer = " " * (self.indent * level)
        body = ",\n".join(inner + item for item in items)
        return f"{open_}\n{body}\n{outer}{close}"


def _to_python(obj: Any) -> Any:
    """Convert NumPy arrays and scalars to Python objects."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _is_number(obj: Any) -> bool:
    obj = _to_python(obj)
    return isinstance(obj, (int, float)) and not isinstance(obj, bool)


def _key_to_str(key: Any) -> str:
    key = _to_python(key)
    if isinstance(key, str):
        return key
    if key is None:
        return "null"
    if isinstance(key, bool):
        return "true" if key else "false"
    if isinstance(key, (int, float)):
        return repr(key)
    raise TypeError(
        f"Keys must be str, int, float, bool or None, "
        f"not {type(key).__name__}"
    )


def _encode_string(text: str) -> str:
    out = []
    for ch in text:
        if ch in _STRING_ESCAPES:
            out.append(_STRING_ESCAPES[ch])
        elif ch < " ":
            out.append(f"\\u{ord(ch):04x}")
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def write_str(
    data: Any,
    indent: Optional[int] = 2,
    float_precision: Optional[int] = None,
    inline_numbers: bool = True,
    allow_nan: bool = True,
    float_format: Optional[str] = None,
    sort_keys: bool = False,
) -> str:
    """
    Serialize a Python object to a JSON string.

    Parameters
    ----------
    data
        Object to serialize (dicts, lists, tuples, strings, numbers,
        booleans, None, NumPy arrays and scalars).
    indent
        Spaces per indentation level; None or 0 writes a single line.
    float_precision
        If set, floats are written with this number of decimals.
    inline_numbers
        Write arrays of numbers on a single line.
    allow_nan
        Write NaN and infinities as ``NaN``, ``Infinity`` and
        ``-Infinity``; if False they raise ValueError.
    float_format
        Format specification for floats (e.g. ``".6g"`` for six
        significant digits); exclusive with ``float_precision``.
    sort_keys
        Write object keys in sorted order.

    Returns
    -------
    str
        JSON document.

    Raises
    ------
    TypeError
        If an object cannot be serialized.
    ValueError
        If both float options are given or the format is invalid.
    """
    writer = _Writer(
        indent,
        float_precision,
        float_format,
        inline_numbers,
        allow_nan,
        sort_keys,
    )
    return writer.encode(data)


def write(
    data: Any,
    file: FileLike,
    indent: Optional[int] = 2,
    float_precision: Optional[int] = None,
    inline_numbers: bool = True,
    allow_nan: bool = True,
    float_format: Optional[str] = None,
    sort_keys: bool = False,
) -> None:
    """
    Write a Python object to a file in JSON format.

    Parameters
    ----------
    data
        Object to serialize (see :func:`write_str`).
    file
        Output path, or an open writable text stream.
    indent
        Spaces per indentation level; None or 0 writes a single line.
    float_precision
        If set, floats are written with this number of decimals.
    inline_numbers
        Write arrays of numbers on a single line.
    allow_nan
        Write NaN and infinities instead of raising ValueError.
    float_format
        Format specification for floats; exclusive with
        ``float_precision``.
    sort_keys
        Write object keys in sorted order.
    """
    text = write_str(
        data,
        indent=indent,
        float_precision=float_precision,
        inline_numbers=inline_numbers,
        allow_nan=allow_nan,
        float_format=float_format,
        sort_keys=sort_keys,
    ) + "\n"

    if isinstance(file, (str, Path)):
        with open(file, "w", encoding="utf-8") as fobj:
            fobj.write(text)
    else:
        file.write(text)
