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
Simple ASCII time series reader and writer.

Simple ASCII is the text format of the EarthScope (formerly IRIS) data
services and of ``mseed2ascii`` ("Simple ASCII time series format",
2010.6.18). Each contiguous segment is a header line::

    TIMESERIES NL_HGN_00_BHZ_R, 12 samples, 40 sps,
    2003-05-29T02:13:22.043400, SLIST, INTEGER, Counts

(on one line) followed by its samples, either as a list (``SLIST``,
one or more values per line, left to right) or as time-sample pairs
(``TSPAIR``, one sample per line). A file may hold any number of
segments, of one or more channels. Sample types are ``INTEGER``,
``FLOAT`` and ``ASCII`` (text, such as station logs).

Reading
-------
- Each segment becomes a record; the source name ``Net_Sta_Loc_Chan``
  (with an optional quality code) becomes the source identifier, the
  units are declared as written and the quality code is kept in
  ``record.head.extra["simpleascii"]``.
- Integers are read as int32 when they fit, floating-point values as
  float64.
- A number of samples different from the header, or a first
  time-sample pair that does not match the start time, gives a
  warning.
- Text segments are read by :class:`SimpleAscii` only.

Writing
-------
- One segment per record, with the number formats of ``mseed2ascii``
  (``%d`` for integers, ``%.10g`` for float64, ``%.8g`` for float32);
  the sampling rate is written exactly (``mseed2ascii`` writes six
  significant digits, too few for rates such as 1/60 Hz).
- Times are written to the microsecond, as required by the format.
- Records without units are written as ``Unknown``.
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from shakelab.libutils.ascii import read_text, write_text
from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import Record, StreamCollection

__all__ = [
    "HEADER_FIELDS",
    "SimpleAscii",
    "SimpleAsciiWarning",
    "simpleasciiread",
    "simpleasciiwrite",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

#: Header fields, in the order of the header line.
HEADER_FIELDS = (
    "SOURCE", "SAMPLES", "SPS", "TIME", "FORMAT", "TYPE", "UNITS",
)

KEYWORD = "TIMESERIES"
STYLES = ("SLIST", "TSPAIR")
TYPES = ("INTEGER", "FLOAT", "ASCII")
DEFAULT_UNITS = "Unknown"

_INT32 = np.iinfo(np.int32)


class SimpleAsciiWarning(UserWarning):
    """Warning about inconsistent Simple ASCII segments."""


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def simpleasciiread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
) -> StreamCollection:
    """
    Read Simple ASCII waveforms into a StreamCollection.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        Simple ASCII file, its content (a string with line breaks, or
        bytes), or an open file.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.

    Returns
    -------
    StreamCollection
        Collection with one record per segment.

    Raises
    ------
    ValueError
        If the file is invalid or holds a text segment.
    """
    from shakelab.signals.base import StreamCollection

    if stream_collection is None:
        stream_collection = StreamCollection()
    for segment in SimpleAscii.read(source):
        stream_collection.append(segment.to_shakelab())
    return stream_collection


def simpleasciiwrite(
    target: Any,
    stream_collection: Any,
    style: str = "slist",
    columns: int = 1,
    owrite: bool = False,
) -> None:
    """
    Write waveforms to a Simple ASCII file.

    Parameters
    ----------
    target : str, path-like or text file object
        Output file, or a file opened in text mode.
    stream_collection : StreamCollection, Stream or Record
        Waveforms to write, one segment per record.
    style : {"slist", "tspair"}, optional
        Sample list or time-sample pairs. Default is ``"slist"``.
    columns : int, optional
        Values per line of a sample list. Default is 1.
    owrite : bool, optional
        If True, overwrite an existing file.

    Raises
    ------
    TypeError
        If ``stream_collection`` is not a waveform container.
    FileExistsError
        If the file exists and ``owrite`` is False.
    ValueError
        If the style, the number of columns or a record is invalid.
    """
    from shakelab.signals.base import Record, Stream, StreamCollection

    if not isinstance(stream_collection, (Record, Stream, StreamCollection)):
        raise TypeError(
            "expected a StreamCollection, a Stream or a Record, got "
            f"{type(stream_collection).__name__}"
        )
    content = "".join(
        SimpleAscii.from_shakelab(record).to_text(style, columns)
        for record in stream_collection.records()
    )
    write_text(target, content, owrite)


# ----------------------------------------------------------------------------
# Simple ASCII container
# ----------------------------------------------------------------------------


class SimpleAscii:
    """
    One segment of a Simple ASCII file.

    Parameters
    ----------
    header : dict, optional
        Header fields by name (see :data:`HEADER_FIELDS`).
    data : array_like or str, optional
        Samples, or the text of an ``ASCII`` segment (one sample per
        character).

    Attributes
    ----------
    header : dict
        ``SOURCE`` (str), ``SAMPLES`` (int), ``SPS`` (float), ``TIME``
        (str, ISO 8601 without zone), ``FORMAT`` (``SLIST`` or
        ``TSPAIR``), ``TYPE`` (``INTEGER``, ``FLOAT`` or ``ASCII``)
        and ``UNITS`` (str).
    data : numpy.ndarray or str
        Samples (int32, int64 or float64), or text.
    """

    def __init__(
        self,
        header: Optional[dict] = None,
        data: Any = None,
    ) -> None:
        if header is not None and not isinstance(header, dict):
            raise TypeError("header must be a dict")
        self.header: dict = dict(header or {})
        if isinstance(data, str):
            self.data: Any = data
        else:
            self.data = np.asarray([] if data is None else data)

    def __repr__(self) -> str:
        return (
            f"SimpleAscii(sid={self.sid!r}, "
            f"type={self.header.get('TYPE')!r}, nsamp={len(self.data)})"
        )

    # -- properties -----------------------------------------------------

    @property
    def sid(self) -> Optional[str]:
        """Source identifier ``NET.STA.LOC.CHA``, None if invalid."""
        codes = _source_codes(self.header.get("SOURCE"))
        return None if codes is None else ".".join(codes[:4])

    @property
    def quality(self) -> Optional[str]:
        """Quality code of the source name, None if absent."""
        codes = _source_codes(self.header.get("SOURCE"))
        return None if codes is None else codes[4]

    @property
    def is_timeseries(self) -> bool:
        """Whether the segment holds samples rather than text."""
        return self.header.get("TYPE") != "ASCII"

    @property
    def start_ns(self) -> Optional[int]:
        """Time of the first sample, Unix ns (UTC); None if not given."""
        value = self.header.get("TIME")
        if not value:
            return None
        return _parse_time(value)

    # -- conversion -----------------------------------------------------

    @classmethod
    def read(cls, source: Any) -> list:
        """
        Read the segments of a Simple ASCII file.

        Parameters
        ----------
        source : str, path-like, bytes or file object
            Simple ASCII file, its content (a string with line breaks,
            or bytes), or an open file.

        Returns
        -------
        list of SimpleAscii

        Raises
        ------
        ValueError
            If the file is invalid.
        """
        return cls.from_text(read_text(source))

    @classmethod
    def from_text(cls, text: str) -> list:
        """
        Parse the segments of a Simple ASCII text.

        Raises
        ------
        ValueError
            If the text is not Simple ASCII or a segment is invalid.
        """
        lines = text.splitlines()
        starts = [
            k for k, line in enumerate(lines)
            if line.startswith(KEYWORD + " ")
        ]
        if not starts or any(line.strip() for line in lines[:starts[0]]):
            raise ValueError(
                f"not a Simple ASCII file: no {KEYWORD} header line"
            )
        starts.append(len(lines))
        return [
            cls._parse_segment(lines[a], lines[a + 1:b], a + 1)
            for a, b in zip(starts[:-1], starts[1:])
        ]

    def to_text(self, style: Optional[str] = None, columns: int = 1) -> str:
        """
        Return the text of the segment.

        Parameters
        ----------
        style : {"slist", "tspair"}, optional
            Sample list or time-sample pairs; the ``FORMAT`` of the
            header if None.
        columns : int, optional
            Values per line of a sample list.

        Raises
        ------
        ValueError
            If the style or the number of columns is invalid, or a
            header field is missing.
        """
        style = (style or self.header.get("FORMAT") or "SLIST").upper()
        if style not in STYLES:
            raise ValueError(f"style must be 'slist' or 'tspair': {style!r}")
        if isinstance(columns, bool) or not isinstance(columns, int):
            raise TypeError("columns must be an integer")
        if columns < 1:
            raise ValueError("columns must be positive")
        h = self.header
        for name in ("SOURCE", "SPS", "TIME", "TYPE", "UNITS"):
            if h.get(name) in (None, ""):
                raise ValueError(f"header field {name} is missing")
        if " " in h["UNITS"] or " " in h["SOURCE"]:
            raise ValueError("source name and units cannot contain spaces")

        line = (
            f"{KEYWORD} {h['SOURCE']}, {len(self.data)} samples, "
            f"{_format_rate(h['SPS'])} sps, {h['TIME']}, {style}, "
            f"{h['TYPE']}, {h['UNITS']}"
        )
        if not self.is_timeseries:
            return f"{line}\n{self.data}\n"

        fmt = _sample_format(self.data)
        values = [fmt % v for v in self.data]
        if style == "TSPAIR":
            start = self.start_ns
            step = 1e9 / h["SPS"]
            body = [
                f"{_format_time(start + round(k * step))}  {v}"
                for k, v in enumerate(values)
            ]
        else:
            body = [
                "".join(
                    v.ljust(10) + "  " if j < columns - 1 else v
                    for j, v in enumerate(values[i:i + columns])
                )
                for i in range(0, len(values), columns)
            ]
        return "\n".join([line] + body) + "\n"

    def write(
        self,
        target: Any,
        style: Optional[str] = None,
        columns: int = 1,
        owrite: bool = False,
    ) -> None:
        """
        Write the segment to a Simple ASCII file.

        Parameters
        ----------
        target : str, path-like or text file object
            Output file, or a file opened in text mode.
        style : {"slist", "tspair"}, optional
            Sample list or time-sample pairs; the ``FORMAT`` of the
            header if None.
        columns : int, optional
            Values per line of a sample list.
        owrite : bool, optional
            If True, overwrite an existing file.

        Raises
        ------
        FileExistsError
            If the file exists and ``owrite`` is False.
        """
        write_text(target, self.to_text(style, columns), owrite)

    def to_shakelab(self) -> Record:
        """
        Convert the segment to a ShakeLab Record.

        Raises
        ------
        ValueError
            If the segment holds text or its sampling rate is not
            positive.
        """
        from shakelab.signals.base import Record

        if not self.is_timeseries:
            raise ValueError(
                "the segment holds text (type ASCII), not a waveform; "
                "read it with SimpleAscii"
            )
        rate = self.header.get("SPS")
        if not rate or rate <= 0:
            raise ValueError(f"invalid sampling rate: {rate!r}")
        start = self.start_ns
        record = Record(
            data=self.data.copy(),
            delta=1.0 / rate,
            time=None if start is None else Date.from_unix_ns(start),
            sid=self.sid,
            units=self.header.get("UNITS"),
        )
        if self.quality:
            record.head.extra["simpleascii"] = {"QUALITY": self.quality}
        return record

    @classmethod
    def from_shakelab(cls, record: Record) -> SimpleAscii:
        """
        Build the Simple ASCII segment of a ShakeLab Record.

        The quality code is taken from
        ``record.head.extra["simpleascii"]``, if present.

        Raises
        ------
        ValueError
            If the record has no source identifier, start time or
            sampling interval.
        """
        head = record.head
        if not head.sid or head.time is None or head.delta is None:
            raise ValueError(
                "a Simple ASCII segment needs the source identifier, "
                "start time and sampling interval of the record"
            )
        quality = head.extra.get("simpleascii", {}).get("QUALITY")
        source = "_".join(head.sid.split(".") + ([quality] if quality
                                                  else []))
        data = np.asarray(record.data)
        if np.iscomplexobj(data):
            raise ValueError("complex data cannot be written")
        header = {
            "SOURCE": source,
            "SAMPLES": data.size,
            "SPS": 1.0 / head.delta,
            "TIME": _format_time(head.time.to_unix_ns()),
            "FORMAT": "SLIST",
            "TYPE": (
                "INTEGER" if np.issubdtype(data.dtype, np.integer)
                else "FLOAT"
            ),
            "UNITS": head.units or DEFAULT_UNITS,
        }
        return cls(header, data)

    def info(self) -> None:
        """Print the header fields."""
        print(repr(self))
        for name in HEADER_FIELDS:
            print(f"{name:>8} = {self.header.get(name)}")

    # -- helpers --------------------------------------------------------

    @classmethod
    def _parse_segment(cls, line: str, body: list, n: int) -> SimpleAscii:
        """Segment of a header line (line ``n``) and its data lines."""
        fields = [f.strip() for f in line[len(KEYWORD):].split(",")]
        if len(fields) != len(HEADER_FIELDS):
            raise ValueError(
                f"line {n}: {len(fields)} header fields, expected "
                f"{len(HEADER_FIELDS)}"
            )
        source, samples, sps, time, style, kind, units = fields
        try:
            header = {
                "SOURCE": source,
                "SAMPLES": int(samples.split()[0]),
                "SPS": float(sps.split()[0]),
                "TIME": time,
                "FORMAT": style.upper(),
                "TYPE": kind.upper(),
                "UNITS": units,
            }
            _parse_time(time)
        except (ValueError, IndexError):
            raise ValueError(f"line {n}: invalid header: {line!r}") from None
        if header["FORMAT"] not in STYLES or header["TYPE"] not in TYPES:
            raise ValueError(
                f"line {n}: unknown format or type: {style}, {kind}"
            )

        segment = cls(header)
        if header["TYPE"] == "ASCII":
            segment.data = "\n".join(body).rstrip("\n")
            return segment

        tokens = " ".join(body).split()
        if header["FORMAT"] == "TSPAIR":
            if len(tokens) % 2:
                raise ValueError(f"line {n}: incomplete time-sample pair")
            times, tokens = tokens[0::2], tokens[1::2]
        try:
            if header["TYPE"] == "INTEGER":
                data = np.array(tokens, dtype=np.int64)
                if data.size and (
                    data.min() >= _INT32.min and data.max() <= _INT32.max
                ):
                    data = data.astype(np.int32)
            else:
                data = np.array(tokens, dtype=np.float64)
        except ValueError as err:
            raise ValueError(f"line {n}: invalid sample: {err}") from None
        segment.data = data

        name = source or f"line {n}"
        if data.size != header["SAMPLES"]:
            warnings.warn(
                f"{name}: {header['SAMPLES']} samples declared, "
                f"{data.size} read",
                SimpleAsciiWarning,
                stacklevel=5,
            )
        if header["FORMAT"] == "TSPAIR" and times:
            offset = abs(_parse_time(times[0]) - _parse_time(time))
            if header["SPS"] > 0 and offset > 0.5e9 / header["SPS"]:
                warnings.warn(
                    f"{name}: first time-sample pair at {times[0]}, "
                    f"start time {time}",
                    SimpleAsciiWarning,
                    stacklevel=5,
                )
        return segment


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _source_codes(source: Optional[str]) -> Optional[tuple]:
    """Codes (net, sta, loc, cha, quality) of a source name."""
    if not source:
        return None
    codes = source.split("_")
    if len(codes) == 4:
        return (*codes, None)
    if len(codes) == 5:
        return (*codes[:4], codes[4] or None)
    return None


def _parse_time(value: str) -> int:
    """Unix ns of an ISO time ``YYYY-MM-DDTHH:MM:SS[.ffffff]`` (UTC)."""
    text = value.strip().rstrip("Z")
    seconds, _, fraction = text.partition(".")
    try:
        moment = datetime.strptime(seconds, "%Y-%m-%dT%H:%M:%S").replace(
            tzinfo=timezone.utc
        )
        ns = int(round(float("0." + (fraction or "0")) * 1e9))
    except ValueError:
        raise ValueError(f"invalid time: {value!r}") from None
    return int(moment.timestamp()) * 10**9 + ns


def _format_time(unix_ns: int) -> str:
    """ISO text of a time, rounded to the microsecond."""
    us = round(unix_ns / 1e3)
    moment = datetime.fromtimestamp(us // 10**6, tz=timezone.utc)
    return moment.strftime("%Y-%m-%dT%H:%M:%S") + f".{us % 10**6:06d}"


def _format_rate(rate: float) -> str:
    """Sampling rate: ten significant digits, more if they are needed."""
    text = f"{rate:.10g}"
    return text if float(text) == rate else repr(float(rate))


def _sample_format(data: np.ndarray) -> str:
    """``%`` format of the samples (as ``mseed2ascii``)."""
    if np.issubdtype(data.dtype, np.integer):
        return "%d"
    if data.dtype == np.float32:
        return "%.8g"
    return "%.10g"
