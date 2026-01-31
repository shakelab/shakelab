# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
Minimal JSON reader/writer utilities using only standard libraries.

This module provides:

* a lightweight, recursive JSON parser (:func:`read`,
  :func:`read_str`), and
* a small JSON writer (:func:`write`) with support for controlled
  float precision and pretty-printing of GeoJSON-like ``coordinates``
  (including 2D and 3D points).

The parser is intentionally minimal and non-validating: it is intended
for well-formed JSON produced by ShakeLab tools. Strings are written
as UTF-8 without ASCII-escaping non-ASCII characters.
"""

from typing import Any, Optional, TextIO, Union


# ============================================================================
# JSON reading utilities
# ============================================================================


def _parse_json_text(text: str) -> Any:
    """
    Parse JSON content from a string and return the corresponding
    Python object.

    This is the core parser used by :func:`read` and :func:`read_str`.
    It implements a minimal subset of JSON sufficient for ShakeLab
    internal use.
    """

    def next_non_whitespace(it: "CharIterator") -> str:
        while True:
            ch = next(it)
            if ch not in " \t\n\r":
                return ch

    def parse_string(it: "CharIterator") -> str:
        result = ""
        while True:
            ch = next(it)
            if ch == '"':
                return result
            if ch == "\\":
                esc = next(it)
                if esc == '"':
                    result += '"'
                elif esc == "\\":
                    result += "\\"
                elif esc == "n":
                    result += "\n"
                elif esc == "t":
                    result += "\t"
                else:
                    raise ValueError(f"Unsupported escape: \\{esc}")
            else:
                result += ch

    def parse_number(it: "CharIterator", first: str) -> Union[int, float]:
        num = first
        while True:
            try:
                ch = peek(it)
                if ch.isdigit() or ch in ".eE+-":
                    num += next(it)
                else:
                    break
            except StopIteration:
                break
        if "." in num or "e" in num or "E" in num:
            return float(num)
        return int(num)

    def parse_array(it: "CharIterator") -> list:
        arr = []
        ch = next_non_whitespace(it)
        if ch == "]":
            return arr
        while True:
            putback(it, ch)
            arr.append(parse_value(it))
            ch = next_non_whitespace(it)
            if ch == "]":
                return arr
            if ch != ",":
                raise ValueError("Expected ',' or ']' in array")
            ch = next_non_whitespace(it)

    def parse_object(it: "CharIterator") -> dict:
        obj = {}
        ch = next_non_whitespace(it)
        if ch == "}":
            return obj
        while True:
            if ch != '"':
                raise ValueError("Expected key string")
            key = parse_string(it)
            ch = next_non_whitespace(it)
            if ch != ":":
                raise ValueError("Expected ':' after key")
            val = parse_value(it)
            obj[key] = val
            ch = next_non_whitespace(it)
            if ch == "}":
                return obj
            if ch != ",":
                raise ValueError("Expected ',' or '}' in object")
            ch = next_non_whitespace(it)

    def expect(it: "CharIterator", string: str) -> None:
        for c in string:
            if next(it) != c:
                raise ValueError(f"Expected {string}")

    def peek(it: "CharIterator") -> str:
        ch = next(it)
        putback(it, ch)
        return ch

    def putback(it: "CharIterator", ch: str) -> None:
        it.stack.append(ch)

    def parse_value(it: "CharIterator") -> Any:
        ch = next_non_whitespace(it)
        if ch == '"':
            return parse_string(it)
        if ch == "{":
            return parse_object(it)
        if ch == "[":
            return parse_array(it)
        if ch == "t":
            expect(it, "rue")
            return True
        if ch == "f":
            expect(it, "alse")
            return False
        if ch == "n":
            expect(it, "ull")
            return None
        if ch == "-" or ch.isdigit():
            return parse_number(it, ch)
        raise ValueError(f"Unexpected character: {ch}")

    class CharIterator:
        """Simple iterator over a string with a push-back stack."""

        def __init__(self, text: str) -> None:
            self.text = text
            self.index = 0
            self.stack: list[str] = []

        def __iter__(self) -> "CharIterator":
            return self

        def __next__(self) -> str:
            if self.stack:
                return self.stack.pop()
            if self.index >= len(self.text):
                raise StopIteration
            ch = self.text[self.index]
            self.index += 1
            return ch

    it = CharIterator(text)
    return parse_value(it)


def read(file: Union[str, TextIO]) -> Any:
    """
    Read JSON content from a file or filename.

    This uses a minimal recursive parser implemented with only
    standard libraries.

    Parameters
    ----------
    file : str or file-like
        Path to a file, or an open readable text stream.

    Returns
    -------
    Any
        Parsed JSON as native Python object.
    """
    close_file = False
    if isinstance(file, str):
        file = open(file, "r", encoding="utf-8")
        close_file = True

    try:
        content = file.read()
        return _parse_json_text(content)
    finally:
        if close_file:
            file.close()


def read_str(text: str) -> Any:
    """
    Parse JSON content from a string and return the Python object.

    This is a convenience wrapper around the internal parser used by
    :func:`read`.
    """
    return _parse_json_text(text)


# ============================================================================
# JSON writing utilities
# ============================================================================


def write(
    data: Any,
    file: Union[str, TextIO],
    indent: int = 2,
    float_precision: Optional[int] = None,
) -> None:
    """
    Write a Python object to a file in JSON format using only standard
    libraries.

    Parameters
    ----------
    data : Any
        The JSON-serializable object to serialize.
    file : str or file-like
        Path to the output file or an open writable text stream.
    indent : int, optional
        Number of spaces for indentation (default is 2).
    float_precision : int, optional
        If set, float values will be rounded to the given number of
        decimal places.
    """

    def _format_float(x: float) -> str:
        if float_precision is not None:
            return f"{round(x, float_precision):.{float_precision}f}"
        return str(x)

    def _dump(obj: Any, level: int = 0, parent_key: str = "") -> str:
        sp = " " * (level * indent)

        if isinstance(obj, dict):
            if not obj:
                return "{}"
            items = []
            for k, v in obj.items():
                key = f'"{k}"'
                val = _dump(v, level + 1, k)
                items.append(f'{sp}{" " * indent}{key}: {val}')
            return "{\n" + ",\n".join(items) + "\n" + sp + "}"

        if isinstance(obj, list):
            if not obj:
                return "[]"
            if parent_key == "coordinates":
                return _format_coordinates(obj, level)
            items = [_dump(i, level + 1, parent_key) for i in obj]
            inner = ",\n".join(
                f'{sp}{" " * indent}{i}' for i in items
            )
            return "[\n" + inner + "\n" + sp + "]"

        if isinstance(obj, str):
            return '"' + obj.replace('"', '\\"') + '"'

        if isinstance(obj, bool):
            return "true" if obj else "false"

        if obj is None:
            return "null"

        if isinstance(obj, float):
            return _format_float(obj)

        return str(obj)

    def _format_coordinates(obj: Any, level: int) -> str:
        """
        Pretty-print GeoJSON-like coordinates.

        Supports the most common 2D/3D patterns:

        * [x, y] or [x, y, z]
        * [[x, y], [x, y], ...]
        * [[x, y, z], [x, y, z], ...]
        * polygon-like [[[x, y], ...], [...]]
        """
        sp0 = " " * (level * indent)
        sp1 = " " * ((level + 1) * indent)
        sp2 = " " * ((level + 2) * indent)

        def fmtpt(pt: list[float]) -> str:
            coords = ", ".join(_format_float(x) for x in pt)
            return f"[{coords}]"

        # Single 2D/3D coordinate: [x, y] or [x, y, z]
        if (
            isinstance(obj, list)
            and len(obj) in (2, 3)
            and all(isinstance(x, float) for x in obj)
        ):
            coords = ", ".join(_format_float(x) for x in obj)
            return f"[{coords}]"

        # Simple LineString-like: [[x, y], [x, y], ...] or 3D
        if (
            isinstance(obj, list)
            and obj
            and all(
                isinstance(pt, list)
                and len(pt) in (2, 3)
                and all(isinstance(x, float) for x in pt)
                for pt in obj
            )
        ):
            lines = [sp0 + "[", sp1 + "["]
            for pt in obj:
                lines.append(f"{sp2}{fmtpt(pt)},")
            lines[-1] = lines[-1].rstrip(",")
            lines.append(sp1 + "]")
            lines.append(sp0 + "]")
            return "\n".join(lines)

        # Polygon-like: [[[x, y], ...], [...], ...]
        if isinstance(obj, list) and obj and all(
            isinstance(ring, list) for ring in obj
        ):
            lines = [sp0 + "["]
            for ring in obj:
                lines.append(sp1 + "[")
                for pt in ring:
                    lines.append(f"{sp2}{fmtpt(pt)},")
                lines[-1] = lines[-1].rstrip(",")
                lines.append(sp1 + "]")
            lines.append(sp0 + "]")
            return "\n".join(lines)

        # Fallback: use generic list formatting
        return _dump(obj, level)

    close_file = False
    if isinstance(file, str):
        file = open(file, "w", encoding="utf-8")
        close_file = True

    try:
        file.write(_dump(data, level=0))
        file.write("\n")
    finally:
        if close_file:
            file.close()
