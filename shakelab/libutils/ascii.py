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
Tabular data from and to delimited ASCII files.

:class:`AsciiTable` stores a table as an ordered list of column names
(``header``) and a list of rows (``db``), each row being a dict that
maps column names to values.

:func:`read_columns` reads numeric tables of any layout, such as the
files written by simple user programs or exported from spreadsheets,
into NumPy columns (:class:`TextColumns`): it recognizes the
delimiter, the decimal comma, the line of column names and
``key: value`` (or ``key = value``) lines before the table, commented
or not. :func:`write_columns` writes them back.

:func:`read_text` and :func:`write_text` read the text of a path, a
content or a file object, and write text to a path or a file object;
they are shared by the text formats of ShakeLab.

Conventions
-----------
- Lines starting with the comment marker (leading blanks ignored) and
  blank lines are skipped when reading.
- Fields are stripped of surrounding blanks; empty fields and null
  markers (see :func:`shakelab.libutils.utils.is_empty`) are read as
  the ``empty`` value.
- Fields may be enclosed in quotes (``"`` by default) to contain the
  delimiter or surrounding blanks; a doubled quote inside a quoted
  field stands for a literal quote. Quoted fields are kept verbatim.
  When writing, only the fields that need it are quoted.
- Columns with an empty name in the header are ignored when reading.
- Errors are reported with exceptions, anomalies that do not prevent
  reading with warnings; nothing is printed.
"""

from __future__ import annotations

import os
import re
import warnings
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping, Sequence, Union

import numpy as np

from shakelab.libutils.utils import cast_value, is_empty

__all__ = [
    "AsciiTable",
    "AsciiWarning",
    "TextColumns",
    "read_columns",
    "read_text",
    "write_columns",
    "write_text",
]

PathLike = Union[str, Path]
DType = Union[Callable, Sequence[Callable]]


class AsciiWarning(UserWarning):
    """Warning about values or lines of a text table."""


class AsciiTable:
    """
    In-memory table with named columns.

    Parameters
    ----------
    header
        Column names. The list is copied.

    Attributes
    ----------
    header : list of str
        Column names, in order.
    db : list of dict
        Table rows, each mapping column names to values.
    """

    def __init__(self, header: Sequence[str] | None = None) -> None:
        self.header: list[str] = list(header) if header else []
        self.db: list[dict[str, Any]] = []

    # ------------------------------------------------------------------------
    # Container protocol
    # ------------------------------------------------------------------------

    def __getitem__(self, key: str) -> list[Any]:
        return self.extract(key)

    def __setitem__(self, key: str, value: Any) -> None:
        values = self._expand(value)

        if key not in self.header:
            self.header.append(key)

        for row, val in zip(self.db, values):
            row[key] = val

    def __iter__(self) -> Iterator[dict[str, Any]]:
        return iter(self.db)

    def __len__(self) -> int:
        return len(self.db)

    def __repr__(self) -> str:
        return (
            f"AsciiTable(rows={len(self.db)}, "
            f"header={self.header!r})"
        )

    @property
    def size(self) -> tuple[int, int]:
        """Number of rows and number of columns."""
        return len(self.db), len(self.header)

    # ------------------------------------------------------------------------
    # Rows and columns
    # ------------------------------------------------------------------------

    def add_element(self, data: Sequence[Any]) -> None:
        """
        Append a row given as values in header order.

        Parameters
        ----------
        data
            Row values, one per column.

        Raises
        ------
        ValueError
            If the number of values differs from the number of columns.
        """
        if len(data) != len(self.header):
            raise ValueError(
                f"row has {len(data)} values, expected "
                f"{len(self.header)}"
            )
        self.db.append(dict(zip(self.header, data)))

    def add_row(self, row: Mapping[str, Any]) -> None:
        """
        Append a row given as a mapping of column names to values.

        Parameters
        ----------
        row
            Row values by column name. Missing columns are set to
            ``None``.

        Raises
        ------
        KeyError
            If ``row`` contains a name that is not in the header.
        """
        unknown = [key for key in row if key not in self.header]
        if unknown:
            raise KeyError(f"keys not in header: {unknown!r}")
        self.db.append({key: row.get(key) for key in self.header})

    def add_key(
        self,
        key: str,
        index: int = -1,
        data: Any = None,
    ) -> None:
        """
        Add a column at a given position.

        Parameters
        ----------
        key
            New column name.
        index
            Position in the header; ``-1`` appends the column.
        data
            Column values: a scalar (broadcast to all rows) or a
            sequence with one value per row.

        Raises
        ------
        ValueError
            If the column already exists or ``data`` has the wrong
            length.
        """
        if key in self.header:
            raise ValueError(f"key already in header: {key!r}")

        values = self._expand(data)

        if index == -1:
            index = len(self.header)
        self.header.insert(index, key)

        for row, val in zip(self.db, values):
            row[key] = val

    def remove_key(self, key: str) -> None:
        """
        Remove a column.

        Parameters
        ----------
        key
            Column name.
        """
        self._require_key(key)
        self.header.remove(key)
        for row in self.db:
            row.pop(key, None)

    def rename_key(self, old_key: str, new_key: str) -> None:
        """
        Rename a column, keeping its position.

        Parameters
        ----------
        old_key
            Current column name.
        new_key
            New column name.
        """
        self._require_key(old_key)
        if new_key in self.header:
            raise ValueError(f"key already in header: {new_key!r}")

        self.header[self.header.index(old_key)] = new_key
        for row in self.db:
            row[new_key] = row.pop(old_key, None)

    def recast(self, key: str, dtype: Callable) -> None:
        """
        Convert the values of a column in place.

        Parameters
        ----------
        key
            Column name.
        dtype
            Target type; empty values become ``None``.
        """
        self._require_key(key)
        for row in self.db:
            row[key] = cast_value(row.get(key), dtype)

    def extract(self, key: str, dtype: Callable | None = None) -> list[Any]:
        """
        Return the values of a column.

        Parameters
        ----------
        key
            Column name.
        dtype
            Optional target type; empty values become ``None``.

        Returns
        -------
        list
            Column values, in row order.
        """
        self._require_key(key)
        values = [row.get(key) for row in self.db]
        if dtype is None:
            return values
        return [cast_value(v, dtype) for v in values]

    def to_array(
        self,
        key: str,
        dtype: Callable = float,
        empty: Any = np.nan,
    ) -> np.ndarray:
        """
        Return the values of a column as a NumPy array.

        Parameters
        ----------
        key
            Column name.
        dtype
            Type applied to each value.
        empty
            Value used for empty entries (NaN by default).

        Returns
        -------
        numpy.ndarray
            Column values, in row order.
        """
        self._require_key(key)
        return np.array(
            [cast_value(row.get(key), dtype, empty) for row in self.db]
        )

    def unique(self, key: str) -> list[Any]:
        """
        Return the distinct values of a column.

        Parameters
        ----------
        key
            Column name.

        Returns
        -------
        list
            Distinct values in order of first appearance.
        """
        self._require_key(key)
        return list(dict.fromkeys(row.get(key) for row in self.db))

    # ------------------------------------------------------------------------
    # Selection and ordering
    # ------------------------------------------------------------------------

    def filter(
        self,
        predicate: Callable[[dict[str, Any]], bool],
    ) -> "AsciiTable":
        """
        Return a new table with the rows that satisfy a condition.

        Parameters
        ----------
        predicate
            Function called with each row (a dict) and returning True
            for the rows to keep.

        Returns
        -------
        AsciiTable
            New table with the same header and copies of the rows.

        Examples
        --------
        >>> strong = tab.filter(lambda r: r["Magnitude"] >= 4.0)
        """
        out = AsciiTable(self.header)
        out.db = [dict(row) for row in self.db if predicate(row)]
        return out

    def select(self, keys: Sequence[str]) -> "AsciiTable":
        """
        Return a new table with a subset of the columns.

        Parameters
        ----------
        keys
            Column names to keep, in the desired order.

        Returns
        -------
        AsciiTable
            New table with copies of the selected values.
        """
        keys = list(keys)
        for key in keys:
            self._require_key(key)
        out = AsciiTable(keys)
        out.db = [{key: row.get(key) for key in keys} for row in self.db]
        return out

    def sort(self, key: str, reverse: bool = False) -> None:
        """
        Sort the rows in place by the values of a column.

        Empty values are placed last, whatever the direction.

        Parameters
        ----------
        key
            Column name.
        reverse
            Sort in descending order.
        """
        self._require_key(key)
        full = [row for row in self.db if not is_empty(row.get(key))]
        void = [row for row in self.db if is_empty(row.get(key))]
        full.sort(key=lambda row: row[key], reverse=reverse)
        self.db = full + void

    def concatenate(self, table: "AsciiTable") -> None:
        """
        Append the rows of another table with the same header.

        Parameters
        ----------
        table
            Table to append; its rows are copied.

        Raises
        ------
        ValueError
            If the headers differ.
        """
        if self.header != table.header:
            raise ValueError("headers do not match")
        self.db.extend(dict(row) for row in table.db)

    # ------------------------------------------------------------------------
    # File input/output
    # ------------------------------------------------------------------------

    def read(
        self,
        ascii_file: PathLike,
        header: Sequence[str] | None = None,
        dtype: DType = str,
        delimiter: str = ",",
        skipline: int = 0,
        comment: str = "#",
        empty: Any = None,
        quotechar: str | None = '"',
    ) -> None:
        """
        Load the table from a delimited ASCII file.

        Any existing content is replaced.

        Parameters
        ----------
        ascii_file
            Input file path.
        header
            Column names. If None, the first data line is used.
        dtype
            Type applied to all columns, or one type per column of the
            file (including unnamed ones).
        delimiter
            Field separator.
        skipline
            Number of initial lines to skip before anything else.
        comment
            Comment marker.
        empty
            Value assigned to empty fields.
        quotechar
            Quote character; None disables quoting.

        Raises
        ------
        ValueError
            If a line has fewer fields than the header, a quoted field
            is not terminated or ``dtype`` has the wrong length.
        """
        with open(ascii_file, "r", encoding="utf-8") as fobj:
            lines = fobj.read().splitlines()[skipline:]

        rows = [
            (n, line) for n, line in enumerate(lines, start=skipline + 1)
            if line.strip() and not line.lstrip().startswith(comment)
        ]

        if header is None:
            if not rows:
                raise ValueError(f"no header found in {ascii_file}")
            n, line = rows.pop(0)
            header = _split(line, delimiter, quotechar, n)

        header = list(header)
        if isinstance(dtype, (list, tuple)):
            if len(dtype) != len(header):
                raise ValueError("dtype and header lengths differ")
            dtypes = list(dtype)
        else:
            dtypes = [dtype] * len(header)

        columns = [
            (i, name, dtypes[i]) for i, name in enumerate(header)
            if not is_empty(name)
        ]

        self.header = [name for _, name, _ in columns]
        self.db = []

        for n, line in rows:
            fields = _split(line, delimiter, quotechar, n)
            if len(fields) < len(header):
                raise ValueError(
                    f"line {n}: {len(fields)} fields, expected "
                    f"{len(header)}"
                )
            self.db.append({
                name: cast_value(fields[i], dtp, empty)
                for i, name, dtp in columns
            })

    def write(
        self,
        ascii_file: PathLike,
        header: bool = True,
        delimiter: str = ",",
        quotechar: str | None = '"',
    ) -> None:
        """
        Save the table to a delimited ASCII file.

        Parameters
        ----------
        ascii_file
            Output file path.
        header
            Whether to write the header line.
        delimiter
            Field separator.
        quotechar
            Quote character used for fields that contain the
            delimiter, the quote itself or surrounding blanks; None
            writes fields verbatim.
        """
        def fmt(value: Any) -> str:
            return _quote(cast_value(value, str, ""), delimiter, quotechar)

        lines = []
        if header:
            lines.append(delimiter.join(fmt(key) for key in self.header))

        for row in self.db:
            lines.append(delimiter.join(
                fmt(row.get(key)) for key in self.header
            ))

        with open(ascii_file, "w", encoding="utf-8") as fobj:
            for line in lines:
                fobj.write(line + "\n")

    # ------------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------------

    def _require_key(self, key: str) -> None:
        if key not in self.header:
            raise KeyError(key)

    def _expand(self, value: Any) -> list[Any]:
        """Broadcast a scalar or validate a per-row sequence."""
        if value is None or isinstance(value, (str, bytes)) \
                or not hasattr(value, "__len__"):
            return [value] * len(self.db)

        if len(value) != len(self.db):
            raise ValueError(
                f"got {len(value)} values for {len(self.db)} rows"
            )
        return list(value)


# ----------------------------------------------------------------------------
# Text sources
# ----------------------------------------------------------------------------


def read_text(source: Any) -> str:
    """
    Return the text of a path, a content or a file object.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        File path, content (a string with line breaks, or bytes in
        UTF-8), or a file opened in text or binary mode.

    Raises
    ------
    TypeError
        If ``source`` is none of these.
    """
    if isinstance(source, (bytes, bytearray, memoryview)):
        return bytes(source).decode("utf-8", errors="replace")
    if isinstance(source, str) and ("\n" in source or "\r" in source):
        return source
    if isinstance(source, (str, os.PathLike)):
        with open(source, "r", encoding="utf-8", errors="replace") as fid:
            return fid.read()
    if hasattr(source, "read"):
        content = source.read()
        if isinstance(content, (bytes, bytearray)):
            content = content.decode("utf-8", errors="replace")
        return content
    raise TypeError("source must be a path, its content or a file object")


def write_text(target: Any, content: str, owrite: bool = False) -> None:
    """
    Write text to a path or a text file object.

    Parameters
    ----------
    target : str, path-like or text file object
        Output file, or a file opened in text mode.
    content
        Text to write.
    owrite
        If True, overwrite an existing file.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    TypeError
        If ``target`` is neither a path nor a file object.
    """
    if isinstance(target, (str, os.PathLike)):
        if os.path.exists(target) and not owrite:
            raise FileExistsError(f"File exists: {target}")
        with open(target, "w", encoding="utf-8") as fid:
            fid.write(content)
    elif hasattr(target, "write"):
        target.write(content)
    else:
        raise TypeError("target must be a path or a text file object")


# ----------------------------------------------------------------------------
# Numeric columns
# ----------------------------------------------------------------------------

# Delimiters tried, in order of preference; None stands for blanks.
DELIMITERS = ("\t", ";", "|", ",", None)

_KEY_VALUE = re.compile(r"^([A-Za-z][^:=]*?)\s*[:=]\s*(.*)$")
_NUMBER = re.compile(r"^[+-]?(\d+([.,]\d*)?|[.,]\d+)([eE][+-]?\d+)?$")
_DECIMAL_COMMA = re.compile(r"^[+-]?\d*,\d+([eE][+-]?\d+)?$")
_DATE = re.compile(
    r"^\d{4}-\d{2}-\d{2}([T ]\d{2}:\d{2}(:\d{2}(\.\d+)?)?)?Z?$"
)
_DATE_SPACE = re.compile(r"(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2})")
_SPECIAL = {"nan", "inf", "-inf", "+inf", "infinity", "-infinity"}


class TextColumns:
    """
    Columns of a text table.

    Attributes
    ----------
    names : list of str or None
        Column names, None if the file has no line of names and none
        were given.
    columns : list of numpy.ndarray
        Columns in file order: float64 for numbers (NaN for empty or
        invalid values), strings for dates (``YYYY-MM-DDTHH:MM:SS``,
        a space in place of ``T`` replaced) and other text.
    metadata : dict
        ``key: value`` pairs found before the table, keys and values
        as written (stripped).
    delimiter : str or None
        Field delimiter, None for blanks.
    decimal : str
        Decimal separator, ``"."`` or ``","``.
    """

    def __init__(
        self,
        columns: Sequence[Any],
        names: Sequence[str] | None = None,
        metadata: Mapping[str, str] | None = None,
        delimiter: str | None = None,
        decimal: str = ".",
    ) -> None:
        self.columns: list[np.ndarray] = [np.asarray(c) for c in columns]
        if names is not None and len(names) != len(self.columns):
            raise ValueError(
                f"{len(names)} names for {len(self.columns)} columns"
            )
        self.names: list[str] | None = (
            None if names is None else list(names)
        )
        self.metadata: dict[str, str] = dict(metadata or {})
        self.delimiter = delimiter
        self.decimal = decimal

    def __repr__(self) -> str:
        return (
            f"TextColumns(names={self.names!r}, nrows={self.nrows}, "
            f"metadata={len(self.metadata)} keys)"
        )

    def __len__(self) -> int:
        return len(self.columns)

    def __getitem__(self, key: int | str) -> np.ndarray:
        if isinstance(key, str):
            if self.names is None or key not in self.names:
                raise KeyError(key)
            key = self.names.index(key)
        return self.columns[key]

    @property
    def nrows(self) -> int:
        """Number of rows."""
        return len(self.columns[0]) if self.columns else 0


def read_columns(
    source: Any,
    names: Sequence[str] | None = None,
    delimiter: str | None = None,
    decimal: str | None = None,
    comment: str = "#",
    skipline: int = 0,
) -> TextColumns:
    """
    Read a numeric table of any layout.

    Before the table, blank lines are skipped, ``key: value`` and
    ``key = value`` lines (commented or not) become metadata and the
    last other line with as many fields as the table gives the column
    names (also when commented, as written by ``numpy.savetxt``). The
    table starts at the first line made of numbers or dates.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        File path, content (a string with line breaks, or bytes), or an
        open file.
    names
        Column names; they replace the line of names of the file.
    delimiter
        Field delimiter; recognized if None (tabulation, semicolon,
        vertical bar, comma or blanks).
    decimal
        Decimal separator, ``"."`` or ``","``; recognized if None (a
        comma when it is not the delimiter and numbers use it).
    comment
        Comment marker.
    skipline
        Number of initial lines to skip before anything else.

    Returns
    -------
    TextColumns

    Raises
    ------
    ValueError
        If no table is found, a row has a different number of fields,
        or ``names`` does not match the number of columns.

    Warns
    -----
    AsciiWarning
        Values that are not numbers in a numeric column (read as NaN),
        and a line of names that does not match the table.
    """
    if decimal not in (None, ".", ","):
        raise ValueError(f"decimal must be '.' or ',': {decimal!r}")
    lines = read_text(source).lstrip("﻿").splitlines()[skipline:]

    metadata: dict[str, str] = {}
    candidates: list[tuple[str, bool]] = []
    start = None
    for k, line in enumerate(lines):
        text = line.strip()
        if not text:
            continue
        commented = bool(comment) and text.startswith(comment)
        if commented:
            text = text[len(comment):].strip()
        elif _is_row(text):
            start = k
            break
        match = _KEY_VALUE.match(text)
        if match:
            metadata[match.group(1).strip()] = match.group(2).strip()
        elif text:
            candidates.append((text, commented))
    if start is None:
        raise ValueError("no table of numbers found")

    rows = [
        (skipline + k + 1, line.strip())
        for k, line in enumerate(lines[start:], start)
        if line.strip()
        and not (comment and line.strip().startswith(comment))
    ]
    sample = [_DATE_SPACE.sub(r"\1T\2", text) for _, text in rows[:20]]
    if delimiter is None:
        delimiter = _guess_delimiter(sample)
    sample_fields = [_fields(text, delimiter) for text in sample]
    ncol = len(sample_fields[0])
    if decimal is None:
        decimal = "," if delimiter != "," and any(
            _DECIMAL_COMMA.match(v) for row in sample_fields for v in row
        ) else "."

    if names is None:
        names = _guess_names(candidates, delimiter, ncol)
    elif len(names) != ncol:
        raise ValueError(f"{len(names)} names for {ncol} columns")

    columns = _numeric_columns(
        [text for _, text in rows], delimiter, decimal, ncol
    )
    if columns is None:
        fields = [
            _fields(_DATE_SPACE.sub(r"\1T\2", text), delimiter)
            for _, text in rows
        ]
        for (n, _), values in zip(rows, fields):
            if len(values) != ncol:
                raise ValueError(
                    f"line {n}: {len(values)} fields, expected {ncol}"
                )
        columns = [
            _column([row[j] for row in fields], decimal,
                    names[j] if names else j + 1)
            for j in range(ncol)
        ]
    return TextColumns(columns, names, metadata, delimiter, decimal)


def write_columns(
    target: Any,
    columns: Sequence[Any],
    names: Sequence[str] | None = None,
    metadata: Mapping[str, Any] | None = None,
    delimiter: str = ",",
    fmt: str = "%.9g",
    comment: str = "#",
    owrite: bool = False,
) -> None:
    """
    Write columns to a text table readable by :func:`read_columns`.

    Parameters
    ----------
    target : str, path-like or text file object
        Output file, or a file opened in text mode.
    columns
        Columns of equal length; numbers are written with ``fmt``,
        other values as text.
    names
        Column names, written as the first line of the table.
    metadata
        Pairs written as ``# key: value`` lines before the table.
    delimiter
        Field delimiter.
    fmt
        ``%`` format of the numbers.
    comment
        Comment marker of the metadata lines.
    owrite
        If True, overwrite an existing file.

    Raises
    ------
    ValueError
        If the columns have different lengths or the names do not
        match them.
    """
    columns = [np.asarray(c) for c in columns]
    if len({c.size for c in columns}) > 1:
        raise ValueError("columns must have the same length")
    if names is not None and len(names) != len(columns):
        raise ValueError(f"{len(names)} names for {len(columns)} columns")

    lines = [
        f"{comment} {key}: {value}" for key, value in (metadata or {}).items()
    ]
    if names is not None:
        lines.append(delimiter.join(str(name) for name in names))
    texts = [
        [fmt % v for v in c] if np.issubdtype(c.dtype, np.number)
        else [str(v) for v in c]
        for c in columns
    ]
    lines.extend(delimiter.join(row) for row in zip(*texts))
    write_text(target, "\n".join(lines) + "\n", owrite)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _split(
    line: str,
    delimiter: str,
    quotechar: str | None = '"',
    lineno: int | None = None,
) -> list[str]:
    """
    Split a line into fields, honouring quoted fields.

    Unquoted fields are stripped of surrounding blanks; quoted fields
    are returned verbatim, with doubled quotes collapsed.
    """
    if not quotechar or quotechar not in line:
        return [field.strip() for field in line.split(delimiter)]

    where = f"line {lineno}: " if lineno is not None else ""
    fields: list[str] = []
    buf: list[str] = []
    quoted = in_quotes = closed = False
    i, n = 0, len(line)

    while i < n:
        ch = line[i]

        if in_quotes:
            if ch == quotechar:
                if line.startswith(quotechar, i + 1):
                    buf.append(ch)
                    i += 2
                    continue
                in_quotes = False
                closed = True
            else:
                buf.append(ch)
            i += 1
            continue

        if line.startswith(delimiter, i):
            text = "".join(buf)
            fields.append(text if quoted else text.strip())
            buf = []
            quoted = closed = False
            i += len(delimiter)
            continue

        if closed:
            if not ch.isspace():
                raise ValueError(
                    f"{where}unexpected character after quoted field"
                )
        elif ch == quotechar and not "".join(buf).strip():
            buf = []
            quoted = in_quotes = True
        else:
            buf.append(ch)
        i += 1

    if in_quotes:
        raise ValueError(f"{where}unterminated quoted field")

    text = "".join(buf)
    fields.append(text if quoted else text.strip())
    return fields


def _quote(text: str, delimiter: str, quotechar: str | None) -> str:
    """Quote a field for writing only when needed."""
    if not quotechar:
        return text

    needs = (
        delimiter in text
        or quotechar in text
        or "\n" in text
        or text != text.strip()
    )
    if not needs:
        return text

    doubled = text.replace(quotechar, quotechar * 2)
    return f"{quotechar}{doubled}{quotechar}"


def _is_value(text: str) -> bool:
    """Whether a field is a number, a date or a null marker."""
    return bool(
        _NUMBER.match(text) or _DATE.match(text)
        or text.lower() in _SPECIAL or is_empty(text)
    )


def _is_row(text: str) -> bool:
    """Whether a line is a row of the table (numbers or dates)."""
    text = _DATE_SPACE.sub(r"\1T\2", text)
    fields = [f for f in re.split(r"[\t;|,\s]+", text) if f]
    return bool(fields) and all(
        _NUMBER.match(f) or _DATE.match(f) or f.lower() in _SPECIAL
        for f in fields
    )


def _fields(text: str, delimiter: str | None) -> list[str]:
    """Fields of a row."""
    if delimiter is None:
        return text.split()
    return [f.strip() for f in text.split(delimiter)]


def _guess_delimiter(rows: list[str]) -> str | None:
    """
    First delimiter splitting the first row into more than one valid
    field; None (blanks) otherwise. Rows with a different number of
    fields are reported by the caller.
    """
    for delimiter in DELIMITERS:
        if delimiter is None or delimiter not in rows[0]:
            continue
        fields = _fields(rows[0], delimiter)
        if len(fields) > 1 and all(_is_value(v) for v in fields):
            return delimiter
    return None


def _numeric_columns(
    rows: list[str],
    delimiter: str | None,
    decimal: str,
    ncol: int,
) -> list[np.ndarray] | None:
    """Columns of a table made only of numbers, None otherwise."""
    if decimal == ",":
        rows = [row.replace(",", ".") for row in rows]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = np.loadtxt(
                rows, delimiter=delimiter, dtype=np.float64, ndmin=2,
                comments=None,
            )
    except ValueError:
        return None
    if table.shape != (len(rows), ncol):
        return None
    return [table[:, j].copy() for j in range(ncol)]


def _guess_names(
    candidates: list[tuple[str, bool]],
    delimiter: str | None,
    ncol: int,
) -> list[str] | None:
    """Column names from the last fitting line before the table."""
    for text, commented in reversed(candidates):
        names = [
            name.strip().strip('"').strip()
            for name in _fields(text, delimiter)
        ]
        if len(names) == ncol:
            return names
        if not commented:
            warnings.warn(
                f"line of names {text!r} has {len(names)} fields, the "
                f"table {ncol}: names ignored",
                AsciiWarning,
                stacklevel=3,
            )
            return None
    return None


def _column(values: list[str], decimal: str, name: Any) -> np.ndarray:
    """Column of floats, or of strings if it holds dates or text."""
    texts = [v.replace(",", ".") if decimal == "," else v for v in values]
    try:
        return np.array(
            [np.nan if is_empty(v) else v for v in texts], dtype=np.float64
        )
    except ValueError:
        pass
    filled = [v for v in values if not is_empty(v)]
    if filled and all(_DATE.match(v) for v in filled):
        return np.array([v.rstrip("Z") for v in values])
    numbers = []
    bad = 0
    for v in texts:
        try:
            numbers.append(float(v) if not is_empty(v) else np.nan)
        except ValueError:
            numbers.append(np.nan)
            bad += 1
    if bad == len(values):
        return np.array(values)
    warnings.warn(
        f"column {name}: {bad} values that are not numbers read as NaN",
        AsciiWarning,
        stacklevel=3,
    )
    return np.array(numbers, dtype=np.float64)
