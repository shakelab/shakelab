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
Waveforms from and to text tables of any layout.

Data are often exchanged as text or CSV files written by simple
programs or spreadsheets, without a standard format. :func:`asciiread`
reads them with :func:`shakelab.libutils.ascii.read_columns`, which
recognizes the delimiter, the decimal comma, the line of column names
and the ``key: value`` lines before the table, and builds one record
per data column.

Parameters
----------
The parameters of the records can be given as arguments or as
``key: value`` lines in the file, with the same names and values:

==========  ==========================================================
Name        Accepted keys in the file (case, blanks and units ignored)
==========  ==========================================================
rate        ``rate``, ``sampling_rate``, ``sample_rate``, ``sps``,
            ``fs``, ``sampling_frequency`` (Hz; ``100 Hz`` accepted)
delta       ``delta``, ``dt``, ``sampling_interval``, ``time_step``
            (s; ``0.01 s`` accepted)
start       ``start``, ``starttime``, ``start_time``, ``t0``
            (ISO 8601, UTC)
sid         ``sid``, ``source_id`` (``NET.STA.LOC.CHA``)
network     ``network``, ``net``
station     ``station``, ``sta``
location    ``location``, ``loc``
units       ``units``, ``unit`` (one value, or one per column)
==========  ==========================================================

Arguments take precedence over the file, and the file over what is
derived from the data (the sampling interval and the start time from a
column of times). Only the samples and the sampling interval are
required; missing start time and units are left None, missing codes
of the source identifiers empty.

Columns
-------
- The column of times is the one named ``time``, ``t``, ``date``,
  ``datetime``, ``timestamp``, ``seconds`` or ``tempo`` (see
  :data:`TIME_COLUMNS`; or given with
  ``time_column``), or the first column if it holds dates. Dates give
  the start time; numbers are seconds from ``start``. Without a
  column of times, the times follow from ``start`` and the sampling
  interval.
- Without a line of names, the names are given with ``columns``; a
  name None skips the column.
- The source identifier of a column is ``sid`` formatted with
  ``{column}`` (e.g. ``"IV.MODE..{column}"``), one of a list of
  ``sid``, the column name if it is an identifier, or otherwise
  ``network.station.location.<column name>``, with empty codes when
  they are unknown (e.g. ``...HNE``, or ``...1`` for the first column
  of a table without names).
- The column name and the metadata of the file are kept in
  ``record.head.extra["ascii"]``.

Writing
-------
:func:`asciiwrite` writes records with the same start time, sampling
and length as columns of one table, preceded by their parameters as
``# key: value`` lines, so that :func:`asciiread` reads them back.
"""

from __future__ import annotations

import re
import warnings
from typing import TYPE_CHECKING, Any, Optional, Sequence

import numpy as np

from shakelab.libutils.ascii import read_columns, write_columns
from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import Record, StreamCollection

__all__ = [
    "PARAMETER_KEYS",
    "TIME_COLUMNS",
    "AsciiSignalWarning",
    "asciiread",
    "asciiwrite",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

#: Keys accepted in the file for each parameter (normalized: lower
#: case, words joined by underscores, units in brackets removed).
PARAMETER_KEYS = {
    "rate": (
        "rate", "sampling_rate", "sample_rate", "samplerate", "sps", "fs",
        "sampling_frequency", "sample_frequency", "sampling_rate_hz",
    ),
    "delta": (
        "delta", "dt", "sampling_interval", "sample_interval",
        "sampling_period", "time_step", "timestep", "sampling_interval_s",
    ),
    "start": (
        "start", "starttime", "start_time", "t0", "first_sample_time",
    ),
    "sid": ("sid", "source_id", "source_identifier"),
    "network": ("network", "net", "network_code"),
    "station": ("station", "sta", "station_code"),
    "location": ("location", "loc", "location_code"),
    "units": ("units", "unit"),
}

#: Names of the column of times (normalized).
TIME_COLUMNS = (
    "time", "t", "times", "date", "datetime", "date_time", "timestamp",
    "time_s", "seconds", "tempo", "temps", "zeit", "tiempo",
)

# Relative tolerance on the sampling interval.
_TOLERANCE = 1e-3

_LEADING_NUMBER = re.compile(r"^\s*([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)")


class AsciiSignalWarning(UserWarning):
    """Warning about inconsistent sampling or parameters of a table."""


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def asciiread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
    columns: Optional[Sequence[Optional[str]]] = None,
    time_column: Optional[int | str] = None,
    sid: Optional[str | Sequence[Optional[str]]] = None,
    network: Optional[str] = None,
    station: Optional[str] = None,
    location: Optional[str] = None,
    rate: Optional[float | str] = None,
    delta: Optional[float | str] = None,
    start: Any = None,
    units: Optional[str | Sequence[Optional[str]]] = None,
    delimiter: Optional[str] = None,
    decimal: Optional[str] = None,
    comment: str = "#",
    skipline: int = 0,
) -> StreamCollection:
    """
    Read the waveforms of a text table.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        Text file, its content (a string with line breaks, or bytes),
        or an open file.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.
    columns : sequence of str or None, optional
        Column names, in file order, replacing the line of names of the
        file; a name None skips the column.
    time_column : int or str, optional
        Index or name of the column of times; recognized if None.
    sid : str or sequence, optional
        Source identifier: a pattern with ``{column}``, one identifier
        (single data column) or one per data column. If None, built
        from the codes and the column names.
    network, station, location : str, optional
        Codes of the source identifiers built from the column names.
    rate : float or str, optional
        Sampling rate in Hz.
    delta : float or str, optional
        Sampling interval in s.
    start : Date or str, optional
        Time of the first sample (or of time zero of a column of
        seconds), UTC.
    units : str or sequence, optional
        Units of all data columns, or one per data column.
    delimiter, decimal, comment, skipline : optional
        Table layout, recognized if None (see
        :func:`shakelab.libutils.ascii.read_columns`).

    Returns
    -------
    StreamCollection
        One record per data column.

    Raises
    ------
    ValueError
        If the file has no numeric data column, the sampling interval
        cannot be determined, or the parameters are contradictory.
    """
    from shakelab.signals.base import Record, StreamCollection

    table = read_columns(
        source, delimiter=delimiter, decimal=decimal, comment=comment,
        skipline=skipline,
        names=None if columns is None else [
            "" if name is None else str(name) for name in columns
        ],
    )
    names = table.names or [str(k + 1) for k in range(len(table))]
    skipped = {
        k for k, name in enumerate(names)
        if columns is not None and columns[k] is None
    }
    params = _parameters(
        table.metadata,
        dict(sid=sid, network=network, station=station, location=location,
             rate=rate, delta=delta, start=start, units=units),
    )

    tcol = _time_column(table, names, time_column, skipped)
    data_cols = [
        k for k in range(len(table))
        if k != tcol and k not in skipped
        and table.columns[k].dtype.kind == "f"
    ]
    if not data_cols:
        raise ValueError("no numeric data column found")

    step, first = _sampling(table, tcol, params)
    sids = _sids(params, [names[k] for k in data_cols])
    units_list = _per_column(params.get("units"), len(data_cols), "units")

    if stream_collection is None:
        stream_collection = StreamCollection()
    for k, code, unit in zip(data_cols, sids, units_list):
        record = Record(
            data=table.columns[k].copy(),
            delta=step,
            time=first,
            sid=code,
            units=unit,
        )
        extra = {"column": names[k]}
        if table.metadata:
            extra["metadata"] = dict(table.metadata)
        record.head.extra["ascii"] = extra
        stream_collection.append(record)
    return stream_collection


def asciiwrite(
    target: Any,
    stream_collection: Any,
    time: Optional[str] = "relative",
    delimiter: str = ",",
    fmt: str = "%.9g",
    owrite: bool = False,
) -> None:
    """
    Write waveforms as columns of a text table.

    The records must have the same start time, sampling interval and
    number of samples. The parameters are written as ``# key: value``
    lines and the source identifiers (or the column names of records
    read from tables) as column names.

    Parameters
    ----------
    target : str, path-like or text file object
        Output file, or a file opened in text mode.
    stream_collection : StreamCollection, Stream or Record
        Waveforms to write.
    time : {"relative", "absolute", None}, optional
        Column of times: seconds from the start time (default), dates
        (ISO 8601 to the microsecond), or none.
    delimiter : str, optional
        Field delimiter. Default is ``","``.
    fmt : str, optional
        ``%`` format of the samples. Default is ``"%.9g"``.
    owrite : bool, optional
        If True, overwrite an existing file.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    ValueError
        If the records cannot form one table or ``time`` is invalid.
    """
    if time not in ("relative", "absolute", None):
        raise ValueError(f"time must be 'relative', 'absolute' or None: "
                         f"{time!r}")
    records = _records(stream_collection)
    if not records:
        raise ValueError("no records to write")
    first = records[0]
    for rec in records[1:]:
        if (
            rec.nsamp != first.nsamp
            or not np.isclose(rec.delta, first.delta, rtol=1e-9, atol=0)
            or (rec.head.time is None) != (first.head.time is None)
            or (rec.head.time is not None
                and rec.head.time.to_unix_ns()
                != first.head.time.to_unix_ns())
        ):
            raise ValueError(
                "the records must have the same start time, sampling "
                "interval and number of samples; write them to separate "
                "files"
            )
    if time == "absolute" and first.head.time is None:
        raise ValueError("absolute times need the start time")

    metadata = {"delta": repr(float(first.delta))}
    if first.head.time is not None:
        metadata["start"] = first.head.time.iso8601
    unit_set = {rec.head.units for rec in records}
    if unit_set != {None}:
        metadata["units"] = (
            str(unit_set.pop()) if len(unit_set) == 1
            else ", ".join(str(rec.head.units or "") for rec in records)
        )

    names = [
        rec.head.sid or rec.head.extra.get("ascii", {}).get("column")
        or f"{k + 1}"
        for k, rec in enumerate(records)
    ]
    columns = [np.asarray(rec.data) for rec in records]
    if time == "relative":
        names.insert(0, "time")
        columns.insert(0, np.array(
            [f"{k * first.delta:.9g}" for k in range(first.nsamp)]
        ))
    elif time == "absolute":
        start = first.head.time.to_unix_ns()
        names.insert(0, "time")
        columns.insert(0, np.array([
            _format_date(start + round(k * first.delta * 1e9))
            for k in range(first.nsamp)
        ]))
    write_columns(
        target, columns, names=names, metadata=metadata,
        delimiter=delimiter, fmt=fmt, owrite=owrite,
    )


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _normalize_key(key: str) -> str:
    """Key in lower case, units in brackets removed, words joined."""
    key = re.sub(r"[\[(].*?[\])]", " ", key).lower()
    return "_".join(re.findall(r"[a-z0-9]+", key))


def _parameters(metadata: dict, arguments: dict) -> dict:
    """Parameters from the arguments, then from the file."""
    found = {}
    for key, value in metadata.items():
        name = _normalize_key(key)
        for param, keys in PARAMETER_KEYS.items():
            if name in keys and param not in found and value != "":
                found[param] = value
    for param, value in arguments.items():
        if value is not None:
            found[param] = value
    return found


def _number(value: Any, name: str) -> float:
    """Positive number of a parameter (a leading number in a string)."""
    if isinstance(value, str):
        match = _LEADING_NUMBER.match(value.replace(",", "."))
        if match is None:
            raise ValueError(f"invalid {name}: {value!r}")
        value = match.group(1)
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"invalid {name}: {value!r}") from None
    if not np.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be positive: {value!r}")
    return number


def _date(value: Any) -> Date:
    """Date of a parameter (Date, ISO string with T or blank)."""
    if isinstance(value, Date):
        return value
    text = str(value).strip().replace(" ", "T", 1)
    try:
        return Date(text)
    except Exception:
        raise ValueError(f"invalid start time: {value!r}") from None


def _time_column(
    table: Any, names: list, time_column: Any, skipped: set
) -> Optional[int]:
    """Index of the column of times, None if there is none."""
    if time_column is not None:
        if isinstance(time_column, str):
            if time_column not in names:
                raise ValueError(f"time column {time_column!r} not found")
            return names.index(time_column)
        if not 0 <= time_column < len(table):
            raise ValueError(f"time column {time_column} out of range")
        return int(time_column)
    for k, name in enumerate(names):
        if k not in skipped and _normalize_key(name) in TIME_COLUMNS:
            return k
    if table.columns and table.columns[0].dtype.kind == "U" \
            and 0 not in skipped:
        return 0
    return None


def _sampling(table: Any, tcol: Optional[int], params: dict) -> tuple:
    """Sampling interval and time of the first sample."""
    step = None
    if "rate" in params:
        step = 1.0 / _number(params["rate"], "sampling rate")
    if "delta" in params:
        given = _number(params["delta"], "sampling interval")
        if step is not None and not np.isclose(
            step, given, rtol=_TOLERANCE, atol=0
        ):
            raise ValueError(
                f"sampling rate {params['rate']!r} and interval "
                f"{params['delta']!r} disagree"
            )
        step = given if step is None else step
    start = _date(params["start"]) if "start" in params else None

    times = None
    if tcol is not None:
        column = table.columns[tcol]
        if column.dtype.kind == "U":
            try:
                ns = column.astype("datetime64[ns]").astype(np.int64)
            except ValueError:
                raise ValueError("invalid dates in the column of times") \
                    from None
            times = (ns - ns[0]) / 1e9
            if "start" not in params:
                start = Date.from_unix_ns(int(ns[0]))
        else:
            times = column.astype(np.float64)
            if start is not None:
                start = Date.from_unix_ns(
                    start.to_unix_ns() + round(times[0] * 1e9)
                )
            elif times[0] != 0:
                warnings.warn(
                    f"times start at {times[0]} s and no start time is "
                    "given: the start time is None",
                    AsciiSignalWarning,
                    stacklevel=4,
                )

    if times is not None and times.size > 1:
        steps = np.diff(times)
        derived = float(np.median(steps))
        if derived <= 0:
            raise ValueError("the times of the column do not increase")
        if np.max(np.abs(steps - derived)) > _TOLERANCE * derived:
            warnings.warn(
                "irregular sampling in the column of times: the median "
                f"interval {derived} s is used",
                AsciiSignalWarning,
                stacklevel=4,
            )
        if step is None:
            step = derived
        elif not np.isclose(step, derived, rtol=_TOLERANCE, atol=0):
            warnings.warn(
                f"sampling interval {step} s given, {derived} s in the "
                "column of times: the given one is used",
                AsciiSignalWarning,
                stacklevel=4,
            )
    if step is None:
        raise ValueError(
            "the sampling interval is unknown: give rate or delta (as "
            "arguments or in the file) or a column of times"
        )
    return step, start


def _sids(params: dict, names: list) -> list:
    """Source identifier of each data column."""
    sid = params.get("sid")
    if isinstance(sid, str):
        if "{column}" in sid:
            return [sid.format(column=name) for name in names]
        if len(names) == 1:
            return [sid]
        raise ValueError(
            "one source identifier for several columns: use a pattern "
            "with {column} or a list"
        )
    if sid is not None:
        return _per_column(list(sid), len(names), "sid")
    net = params.get("network") or ""
    sta = params.get("station") or ""
    loc = params.get("location") or ""
    if "." in str(sta) and not net:
        net, sta = str(sta).split(".", 1)
    return [
        name if name.count(".") == 3 else f"{net}.{sta}.{loc}.{name}"
        for name in names
    ]


def _per_column(value: Any, n: int, name: str) -> list:
    """One value per data column."""
    if value is None or isinstance(value, str):
        if isinstance(value, str) and n > 1 and "," in value:
            value = [v.strip() or None for v in value.split(",")]
        else:
            return [value] * n
    value = list(value)
    if len(value) != n:
        raise ValueError(f"{len(value)} {name} values for {n} data columns")
    return value


def _records(stream_collection: Any) -> list:
    """Records of a StreamCollection, Stream or Record."""
    from shakelab.signals.base import Record, Stream, StreamCollection

    if isinstance(stream_collection, Record):
        return [stream_collection]
    if isinstance(stream_collection, Stream):
        return list(stream_collection)
    if isinstance(stream_collection, StreamCollection):
        return [rec for stream in stream_collection for rec in stream]
    raise TypeError(
        "expected a StreamCollection, a Stream or a Record, got "
        f"{type(stream_collection).__name__}"
    )


def _format_date(unix_ns: int) -> str:
    """ISO date of a time, to the microsecond."""
    return Date.from_unix_ns(unix_ns).iso8601
