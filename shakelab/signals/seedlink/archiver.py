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
MiniSEED archive of the records received by a SeedLink client.

An :class:`Archive` appends the MiniSEED records, as received, to files
whose paths follow a template:

- the fields ``{net}``, ``{sta}``, ``{loc}`` and ``{cha}`` separate the
  files by network, station or channel (a template without fields
  writes a single file);
- ``rotate`` starts a new file at regular intervals, aligned with the
  clock; ``{time}`` (with a ``strftime`` format, e.g.
  ``{time:%Y%m%dT%H%M}``) is the start of the interval, ``{year}`` and
  ``{jday}`` (day of the year, three digits) its date;
- ``layout="sds"`` writes the SeisComP Data Structure: one file per
  channel and day in ``YEAR/NET/STA/CHA.D/NET.STA.LOC.CHA.D.YEAR.DAY``.

Each record goes to the file of its start time, unsplit. Files are
appended to, so that an interrupted archive can be resumed, and every
record is flushed to the operating system as soon as it is written.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import string
from collections import OrderedDict
from numbers import Real
from pathlib import Path
from typing import Any, BinaryIO, Optional

__all__ = ["SDS_TEMPLATE", "Archive"]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

# Path of a file of the SeisComP Data Structure, below the archive root.
SDS_TEMPLATE = (
    "{year}/{net}/{sta}/{cha}.D/{net}.{sta}.{loc}.{cha}.D.{year}.{jday}"
)

# Files kept open at the same time; the least recently used is closed.
MAX_OPEN_FILES = 64

FIELDS = ("net", "sta", "loc", "cha", "time", "year", "jday")

_UNITS = {"s": 1, "min": 60, "h": 3600, "d": 86400}
_ROTATE = re.compile(r"^\s*(\d+(?:\.\d*)?)\s*(s|min|h|d)\s*$")
_NS = 10 ** 9
_EPOCH = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)


# ----------------------------------------------------------------------------
# Archive
# ----------------------------------------------------------------------------


class Archive:
    """
    Files of MiniSEED records, organized by a path template.

    Parameters
    ----------
    path : str or path-like
        Path template (see the module notes), or the root directory
        with ``layout="sds"``.
    rotate : float or str, optional
        Duration of each file: seconds, or a string with a unit
        (``"30s"``, ``"10min"``, ``"1h"``, ``"1d"``). Intervals are
        aligned with the clock (multiples of the duration from the Unix
        epoch, hence from midnight UTC for divisors of a day). Default:
        no rotation.
    layout : {"sds"}, optional
        Predefined layout; ``"sds"`` implies daily files.

    Raises
    ------
    TypeError
        If an argument has an unsupported type.
    ValueError
        If the template has unknown fields, uses ``{time}`` without
        ``rotate``, or if ``rotate`` is invalid or combined with
        ``layout="sds"``.

    Examples
    --------
    >>> Archive("data/{net}.{sta}.{loc}.{cha}.{time:%Y%m%dT%H}.mseed",
    ...         rotate="1h")
    >>> Archive("archive", layout="sds")
    """

    def __init__(
        self,
        path: Any,
        rotate: Any = None,
        layout: Optional[str] = None,
    ) -> None:
        if not isinstance(path, (str, os.PathLike)):
            raise TypeError("path must be a string or a path")
        if layout not in (None, "sds"):
            raise ValueError("layout must be None or 'sds'")

        if layout == "sds":
            if rotate is not None:
                raise ValueError("an SDS archive has daily files: no rotate")
            template = os.path.join(os.fspath(path), SDS_TEMPLATE)
        else:
            template = os.fspath(path)

        self.template = template
        self.layout = layout
        self.rotate_ns = _rotate_ns(rotate)
        _check_template(template, self.rotate_ns is not None)

        self._files: OrderedDict = OrderedDict()

    def __repr__(self) -> str:
        return f"Archive({self.template!r})"

    def path(self, sid: str, start_ns: int) -> Path:
        """
        Path of the file of a record.

        Parameters
        ----------
        sid : str
            Source identifier ``NET.STA.LOC.CHA``.
        start_ns : int
            Start time of the record, in nanoseconds from the Unix
            epoch.

        Returns
        -------
        pathlib.Path
        """
        network, station, location, channel = sid.split(".")
        if self.rotate_ns is not None:
            start_ns -= start_ns % self.rotate_ns
        time = _EPOCH + _dt.timedelta(microseconds=start_ns // 1000)
        return Path(self.template.format(
            net=network,
            sta=station,
            loc=location,
            cha=channel,
            time=time,
            year=f"{time.year:04d}",
            jday=f"{time.timetuple().tm_yday:03d}",
        ))

    def write(self, payload: bytes, sid: str, start_ns: int) -> Path:
        """
        Append a record to its file.

        Parameters
        ----------
        payload : bytes
            MiniSEED record, as received.
        sid : str
            Source identifier of the record.
        start_ns : int
            Start time of the record (ns from the Unix epoch).

        Returns
        -------
        pathlib.Path
            Path of the file written.
        """
        path = self.path(sid, start_ns)
        handle = self._open(path)
        handle.write(payload)
        handle.flush()
        return path

    def close(self) -> None:
        """Close the open files."""
        while self._files:
            _, handle = self._files.popitem(last=False)
            handle.close()

    def _open(self, path: Path) -> BinaryIO:
        """File of a path, opened if needed."""
        handle = self._files.get(path)
        if handle is not None:
            self._files.move_to_end(path)
            return handle
        if len(self._files) >= MAX_OPEN_FILES:
            _, oldest = self._files.popitem(last=False)
            oldest.close()
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(path, "ab")
        self._files[path] = handle
        return handle


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _rotate_ns(rotate: Any) -> Optional[int]:
    """Rotation interval in nanoseconds (None without rotation)."""
    if rotate is None:
        return None
    if isinstance(rotate, str):
        match = _ROTATE.match(rotate)
        if match is None:
            raise ValueError(
                f"invalid rotate {rotate!r}: use seconds or a number with "
                "a unit among s, min, h, d (e.g. '10min')"
            )
        seconds = float(match.group(1)) * _UNITS[match.group(2)]
    elif isinstance(rotate, bool) or not isinstance(rotate, Real):
        raise TypeError("rotate must be a number of seconds or a string")
    else:
        seconds = float(rotate)
    if not seconds > 0 or seconds != seconds or seconds == float("inf"):
        raise ValueError("rotate must be positive and finite")
    return round(seconds * _NS)


def _check_template(template: str, rotate: bool) -> None:
    """Validate the fields of a path template."""
    try:
        fields = [
            name for _, name, _, _ in string.Formatter().parse(template)
            if name is not None
        ]
    except ValueError as error:
        raise ValueError(f"invalid path template: {error}") from None

    for name in fields:
        if name not in FIELDS:
            raise ValueError(
                f"unknown field {{{name}}} in the path template: use "
                + ", ".join(f"{{{field}}}" for field in FIELDS)
            )
    if "time" in fields and not rotate:
        raise ValueError(
            "{time} needs rotate (otherwise every record would have its "
            "own file)"
        )
