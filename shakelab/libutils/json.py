# ****************************************************************************
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
This module provides minimal read and write utilities for JSON-formatted data,
implemented using only Python standard libraries. It includes a custom JSON
writer with support for controlled float precision and pretty-printing of
GeoJSON coordinates, and a lightweight recursive JSON parser.
These tools are particularly useful with ShakeLab's internal data models.
"""

from typing import Any, Union, TextIO, Any


def read(file: Union[str, TextIO]) -> Any:
    """
    Minimal JSON parser using only only standard libraries.

    Parameters
    ----------
    file : str or file-like

    Returns
    -------
    Any
        Parsed JSON as native Python object.
    """

    def parse_value(it):
        ch = next_non_whitespace(it)
        if ch == '"':
            return parse_string(it)
        elif ch == '{':
            return parse_object(it)
        elif ch == '[':
            return parse_array(it)
        elif ch == 't':
            expect(it, 'rue')
            return True
        elif ch == 'f':
            expect(it, 'alse')
            return False
        elif ch == 'n':
            expect(it, 'ull')
            return None
        elif ch == '-' or ch.isdigit():
            return parse_number(it, ch)
        else:
            raise ValueError(f"Unexpected character: {ch}")

    def next_non_whitespace(it):
        while True:
            ch = next(it)
            if ch not in ' \t\n\r':
                return ch

    def parse_string(it):
        result = ''
        while True:
            ch = next(it)
            if ch == '"':
                return result
            elif ch == '\\':
                esc = next(it)
                if esc == '"':
                    result += '"'
                elif esc == '\\':
                    result += '\\'
                elif esc == 'n':
                    result += '\n'
                elif esc == 't':
                    result += '\t'
                else:
                    raise ValueError(f"Unsupported escape: \\{esc}")
            else:
                result += ch

    def parse_number(it, first):
        num = first
        while True:
            try:
                ch = peek(it)
                if ch.isdigit() or ch in '.eE+-':
                    num += next(it)
                else:
                    break
            except StopIteration:
                break
        if '.' in num or 'e' in num or 'E' in num:
            return float(num)
        return int(num)

    def parse_array(it):
        arr = []
        ch = next_non_whitespace(it)
        if ch == ']':
            return arr
        while True:
            putback(it, ch)
            arr.append(parse_value(it))
            ch = next_non_whitespace(it)
            if ch == ']':
                return arr
            elif ch != ',':
                raise ValueError("Expected ',' or ']' in array")
            ch = next_non_whitespace(it)

    def parse_object(it):
        obj = {}
        ch = next_non_whitespace(it)
        if ch == '}':
            return obj
        while True:
            if ch != '"':
                raise ValueError("Expected key string")
            key = parse_string(it)
            ch = next_non_whitespace(it)
            if ch != ':':
                raise ValueError("Expected ':' after key")
            val = parse_value(it)
            obj[key] = val
            ch = next_non_whitespace(it)
            if ch == '}':
                return obj
            elif ch != ',':
                raise ValueError("Expected ',' or '}' in object")
            ch = next_non_whitespace(it)

    def expect(it, string):
        for c in string:
            if next(it) != c:
                raise ValueError(f"Expected {string}")

    def peek(it):
        ch = next(it)
        putback(it, ch)
        return ch

    def putback(it, ch):
        it.stack.append(ch)

    class CharIterator:
        def __init__(self, text):
            self.text = text
            self.index = 0
            self.stack = []

        def __iter__(self):
            return self

        def __next__(self):
            if self.stack:
                return self.stack.pop()
            if self.index >= len(self.text):
                raise StopIteration
            ch = self.text[self.index]
            self.index += 1
            return ch

    close_file = False
    if isinstance(file, str):
        file = open(file, 'r', encoding='utf-8')
        close_file = True

    try:
        content = file.read()
        it = CharIterator(content)
        return parse_value(it)
    finally:
        if close_file:
            file.close()


def write(data: Any,
               file: Union[str, TextIO],
               indent: int = 2,
               float_precision: int = None) -> None:
    """
    Write a dictionary or list to a file in JSON format using only
    standard libraries.

    Parameters
    ----------
    data : Any
        The dictionary or list to serialize.
    file : str or file-like
        Path to the output file or an open writable stream.
    indent : int, optional
        Number of spaces for indentation (default is 2).
    float_precision : int, optional
        If set, float values will be rounded to the given number of
        decimal places.
    """

    def _format_float(x: float) -> str:
        if float_precision is not None:
            return f'{round(x, float_precision):.{float_precision}f}'
        return str(x)

    def _dump(obj: Any, level: int = 0, parent_key: str = '') -> str:
        sp = ' ' * (level * indent)

        if isinstance(obj, dict):
            if not obj:
                return '{}'
            items = []
            for k, v in obj.items():
                key = f'"{k}"'
                val = _dump(v, level + 1, k)
                items.append(f'{sp}{" " * indent}{key}: {val}')
            return '{\n' + ',\n'.join(items) + '\n' + sp + '}'

        elif isinstance(obj, list):
            if not obj:
                return '[]'
            if parent_key == 'coordinates':
                return _format_coordinates(obj, level)
            items = [_dump(i, level + 1, parent_key) for i in obj]
            return '[\n' + ',\n'.join(
                f'{sp}{" " * indent}{i}' for i in items
            ) + '\n' + sp + ']'

        elif isinstance(obj, str):
            return '"' + obj.replace('"', '\\"') + '"'

        elif isinstance(obj, bool):
            return 'true' if obj else 'false'

        elif obj is None:
            return 'null'

        elif isinstance(obj, float):
            return _format_float(obj)

        else:
            return str(obj)

    def _format_coordinates(obj: Any, level: int) -> str:
        sp0 = ' ' * (level * indent)
        sp1 = ' ' * ((level + 1) * indent)
        sp2 = ' ' * ((level + 2) * indent)
        sp3 = ' ' * ((level + 3) * indent)

        def fmtpt(pt):
            return f'[{_format_float(pt[0])}, {_format_float(pt[1])}]'

        if isinstance(obj, list) and all(isinstance(x, float) for x in obj):
            return f'[{_format_float(obj[0])}, {_format_float(obj[1])}]'

        if all(isinstance(pt, list) and len(pt) == 2 and
               all(isinstance(x, float) for x in pt) for pt in obj):
            lines = [sp0 + '[', sp1 + '[']
            for pt in obj:
                lines.append(f'{sp2}{fmtpt(pt)},')
            lines[-1] = lines[-1].rstrip(',')
            lines.append(sp1 + ']', sp0 + ']')
            return '\n'.join(lines)

        if all(isinstance(ring, list) for ring in obj):
            lines = [sp0 + '[']
            for ring in obj:
                lines.append(sp1 + '[')
                for pt in ring:
                    lines.append(f'{sp2}{fmtpt(pt)},')
                lines[-1] = lines[-1].rstrip(',')
                lines.append(sp1 + ']')
            lines.append(sp0 + ']')
            return '\n'.join(lines)

        raise ValueError('Unsupported coordinates structure')

    close_file = False
    if isinstance(file, str):
        file = open(file, 'w', encoding='utf-8')
        close_file = True

    try:
        file.write(_dump(data, level=0))
        file.write('\n')
    finally:
        if close_file:
            file.close()
