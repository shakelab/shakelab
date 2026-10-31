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
- Columns with an empty name in the header are ignored when reading.
- Errors are reported with exceptions; nothing is printed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterator, Sequence, Union

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

        Raises
        ------
        ValueError
            If a line has fewer fields than the header or ``dtype`` has
            the wrong length.
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
            header = _split(rows.pop(0)[1], delimiter)

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
            fields = _split(line, delimiter)
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
        """
        lines = []
        if header:
            lines.append(delimiter.join(self.header))

        for row in self.db:
            lines.append(delimiter.join(
                cast_value(row.get(key), str, "") for key in self.header
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


def _split(line: str, delimiter: str) -> list[str]:
    """Split a line into stripped fields."""
    return [field.strip() for field in line.split(delimiter)]
