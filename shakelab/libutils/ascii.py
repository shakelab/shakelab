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
Tabular data from and to delimited ASCII files.

:class:`AsciiTable` stores a table as an ordered list of column names
(``header``) and a list of rows (``db``), each row being a dict that
maps column names to values.

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
- Errors are reported with exceptions; nothing is printed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterator, Mapping, Sequence, Union

import numpy as np

from shakelab.libutils.utils import cast_value, is_empty

__all__ = ["AsciiTable"]

PathLike = Union[str, Path]
DType = Union[Callable, Sequence[Callable]]


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
