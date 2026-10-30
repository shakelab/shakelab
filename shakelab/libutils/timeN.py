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
Time and date utilities.

This module provides :class:`Date`, an absolute UTC instant, together
with calendar helpers, ISO 8601 parsing and formatting, and conversions
between calendar components and seconds from the internal epoch.

Design choices
--------------
- The implementation intentionally does not use the standard
  ``datetime`` module.
- A :class:`Date` stores a single integer: nanoseconds from the internal
  epoch. Comparisons, hashing and arithmetic are exact integer
  operations; calendar fields are derived on demand and cached.
- Conversions of real numbers and strings use :class:`decimal.Decimal`
  with a private context, independent of the global decimal context.
- ISO 8601 strings with an explicit offset are accepted, but they are
  converted to UTC; the offset itself is not retained.

Conventions
-----------
- The internal epoch is ``0001-01-01T00:00:00Z``. This is not Unix time;
  use :meth:`Date.from_unix`, :meth:`Date.to_unix` and their
  nanosecond variants to exchange Unix timestamps.
- Dates use the proleptic Gregorian calendar with astronomical year
  numbering, as in ISO 8601: year ``0`` is 1 BC, year ``-1`` is 2 BC.
  Instants before the internal epoch have negative seconds.
- ISO 8601 years outside ``0000``-``9999`` carry an explicit sign and
  at least four digits (e.g. ``-0459``, ``+10000``).
- Leap seconds are not represented: seconds lie in ``[0, 60)``.
- The resolution of a :class:`Date` is one nanosecond. Finer inputs
  are rounded to the nearest nanosecond, ties to even.
- Time shifts are real numbers expressed in seconds: Python ``int`` and
  ``float``, :class:`decimal.Decimal`, or NumPy integer and floating
  scalars. Booleans and non-finite values are rejected.
- A bare number used as a date (construction, comparison, equality) is
  always interpreted as seconds from the internal epoch, consistently
  with :meth:`Date.to_seconds`.
- The timing error of a date is metadata: it is preserved by copies and
  time shifts, and ignored by comparisons and hashing.
"""

from __future__ import annotations

from decimal import (
    ROUND_HALF_EVEN,
    Context,
    Decimal,
    InvalidOperation,
    localcontext,
)
import math
from numbers import Integral, Real
import re
import time as systime
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union


__all__ = [
    "Date",
    "MINUTE_SECONDS",
    "HOUR_SECONDS",
    "DAY_SECONDS",
    "YEAR_DAYS",
    "UNIX_EPOCH_SECONDS",
    "leap_check",
    "leap_num",
    "days_in_month",
    "year_days",
    "date_to_second",
    "second_to_date",
    "from_ordinal_day",
    "to_ordinal_day",
    "read_iso8601_date",
    "write_iso8601_date",
    "utc_now",
    "to_decimal",
]


MINUTE_SECONDS = Decimal("60")
HOUR_SECONDS = Decimal("3600")
DAY_SECONDS = Decimal("86400")
YEAR_DAYS = 365

# Seconds from the internal epoch to 1970-01-01T00:00:00Z.
UNIX_EPOCH_SECONDS = Decimal("62135596800")

_NS_PER_SECOND = 1_000_000_000
_NS_PER_MINUTE = 60 * _NS_PER_SECOND
_NS_PER_HOUR = 60 * _NS_PER_MINUTE
_NS_PER_DAY = 24 * _NS_PER_HOUR
_UNIX_EPOCH_NS = 62_135_596_800 * _NS_PER_SECOND
_NS_DECIMAL = Decimal(_NS_PER_SECOND)

# Private context: results do not depend on the global decimal context.
_DECIMAL_CONTEXT = Context(prec=100, rounding=ROUND_HALF_EVEN)

# Sentinel for optional keyword arguments of Date.replace().
_UNSET = object()


ISO8601_PATTERN = re.compile(
    r"^"
    r"(?P<year>[+-][0-9]{4,}|[0-9]{4})"
    r"(?:-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})"
    r"|-?(?P<ordinal>[0-9]{3}))"
    r"T"
    r"(?P<hour>[0-9]{2}):"
    r"(?P<minute>[0-9]{2}):"
    r"(?P<second>[0-9]{2}(?:\.[0-9]+)?)"
    r"(?P<timezone>Z|[+-][0-9]{2}(?::?[0-9]{2})?"
    r"(?::?[0-9]{2})?)?"
    r"$"
)


TimeShift = Union[int, float, Decimal]
DateComponents = List[Union[int, float]]
DateFields = Tuple[int, int, int, int, int, int, int]


# ---------------------------------------------------------------------------
# Date class
# ---------------------------------------------------------------------------

class Date:
    """
    Absolute UTC date and time.

    Parameters
    ----------
    date
        Input date. Supported formats are:

        - an ISO 8601 string, or ``"now"`` for the current UTC time;
        - seconds from the internal epoch (see module conventions);
        - a calendar sequence ``[year, month, day, hour, minute,
          second]`` or an ordinal sequence ``[year, ordinal_day, hour,
          minute, second]``;
        - a dictionary with ``year`` and either ``month`` and ``day`` or
          one of ``ordinal_day``, ``day_of_year``, ``year_day``; time
          fields default to zero;
        - another :class:`Date`.

        If omitted, the date is left unset.
    error
        Timing error in seconds (finite, non-negative). If omitted, it
        is inherited when ``date`` is a :class:`Date`, and is zero
        otherwise.
    timezone
        Only ``"UTC"`` is supported.

    Raises
    ------
    TypeError
        If the input format or the error type is not supported.
    ValueError
        If the input does not describe a valid date, or the error is
        negative or not finite.
    NotImplementedError
        If a time zone other than UTC is requested.

    Notes
    -----
    Internally a date is an integer number of nanoseconds from the
    internal epoch. The calendar fields (:attr:`year`, :attr:`month`,
    :attr:`day`, :attr:`hour`, :attr:`minute`, :attr:`second`,
    :attr:`nanosecond`) are derived from it.

    Assigning a field, e.g. ``date.year = 2027``, validates the
    resulting date and updates it in place; use :meth:`replace` to
    change several fields at once, since intermediate dates (such as
    February 31) are rejected.

    Bare numbers compare with dates as seconds from the internal epoch,
    exactly (without rounding to nanoseconds). Equality with a
    non-finite number is False; ordering against a non-finite number
    raises ``ValueError``.

    An unset date can be filled later, for example with
    :meth:`strparse` or :meth:`set_date`. Operations that need an
    instant (comparisons, arithmetic, hashing, field assignment,
    conversion to seconds or strings) raise ``ValueError`` on an unset
    date, and its fields read as ``None``.

    Dates are hashable and dates that compare equal (also with bare
    numbers) have the same hash. Do not modify a date while it is used
    as a dictionary key or set member.
    """

    # Make NumPy scalars defer binary operators to Date (NEP 13), so that
    # expressions such as ``np.int64(5) + date`` call Date.__radd__.
    __array_ufunc__ = None

    def __init__(
        self,
        date: Any = None,
        error: Optional[float] = None,
        timezone: str = "UTC",
    ) -> None:
        self._ns: Optional[int] = None
        self._fields: Optional[DateFields] = None

        if error is None:
            error = date.error if isinstance(date, Date) else 0.0

        self.error = _validate_error(error)

        if timezone != "UTC":
            raise NotImplementedError("Only UTC timezone is supported.")

        if date is None:
            return

        if isinstance(date, str) and date == "now":
            self._set_ns(systime.time_ns() + _UNIX_EPOCH_NS)
            return

        self.set_date(date)

    def __repr__(self) -> str:
        if not self.is_set:
            return "Date(None)"

        return f"Date({self.iso8601!r})"

    def __str__(self) -> str:
        if not self.is_set:
            return "Date(None)"

        return self.iso8601

    def __eq__(self, value: object) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        if not isinstance(value, Date) and not _is_finite_number(value):
            return False

        return self._compare(value) == 0

    def __hash__(self) -> int:
        ns = self._require_ns()
        seconds, fraction = divmod(ns, _NS_PER_SECOND)

        # Consistent with the hash of equal numbers (int, float, Decimal).
        if fraction == 0:
            return hash(seconds)

        return hash(_ns_to_decimal(ns))

    def __lt__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self._compare(value) < 0

    def __le__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self._compare(value) <= 0

    def __gt__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self._compare(value) > 0

    def __ge__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self._compare(value) >= 0

    def __add__(self, value: TimeShift) -> Date:
        if not _is_time_shift(value):
            return NotImplemented

        ns = self._require_ns() + _seconds_to_ns(value)

        return Date._from_ns(ns, self.error)

    __radd__ = __add__

    def __sub__(self, value: Union[Date, TimeShift]) -> Union[Date, float]:
        if isinstance(value, Date):
            difference = self._require_ns() - value._require_ns()
            return difference / _NS_PER_SECOND

        if _is_time_shift(value):
            ns = self._require_ns() - _seconds_to_ns(value)
            return Date._from_ns(ns, self.error)

        return NotImplemented

    # -----------------------------------------------------------------------
    # Alternative constructors
    # -----------------------------------------------------------------------

    @classmethod
    def now(cls) -> Date:
        """
        Return the current UTC time.

        Returns
        -------
        Date
            Current UTC date, with the nanosecond resolution of
            :func:`time.time_ns`.
        """
        return cls._from_ns(systime.time_ns() + _UNIX_EPOCH_NS)

    @classmethod
    def from_iso8601(cls, date_string: str) -> Date:
        """
        Create a Date from an ISO 8601 string.

        Parameters
        ----------
        date_string
            ISO 8601 date string.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls(date_string)

    @classmethod
    def from_dict(cls, date_dict: Dict[str, Any]) -> Date:
        """
        Create a Date from a dictionary.

        Parameters
        ----------
        date_dict
            Dictionary with calendar fields or ordinal day.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls(date_dict)

    @classmethod
    def from_ordinal(
        cls,
        year: int,
        ordinal_day: int,
        hour: int = 0,
        minute: int = 0,
        second: float = 0.0,
    ) -> Date:
        """
        Create a Date from year and ordinal day.

        Parameters
        ----------
        year
            Year.
        ordinal_day
            Day of year.
        hour
            Hour.
        minute
            Minute.
        second
            Second.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls([year, ordinal_day, hour, minute, second])

    @classmethod
    def from_unix(cls, seconds: TimeShift) -> Date:
        """
        Create a Date from a Unix timestamp.

        Parameters
        ----------
        seconds
            Seconds from 1970-01-01T00:00:00Z. Negative values describe
            instants before 1970.

        Returns
        -------
        Date
            Corresponding UTC date.
        """
        return cls._from_ns(_seconds_to_ns(seconds) + _UNIX_EPOCH_NS)

    @classmethod
    def from_unix_ns(cls, nanoseconds: int) -> Date:
        """
        Create a Date from a Unix timestamp in integer nanoseconds.

        Parameters
        ----------
        nanoseconds
            Nanoseconds from 1970-01-01T00:00:00Z, as returned by
            :func:`time.time_ns`.

        Returns
        -------
        Date
            Corresponding UTC date.

        Raises
        ------
        TypeError
            If ``nanoseconds`` is not an integer.
        """
        if (
            isinstance(nanoseconds, bool)
            or not isinstance(nanoseconds, Integral)
        ):
            raise TypeError("Nanoseconds must be an integer.")

        return cls._from_ns(int(nanoseconds) + _UNIX_EPOCH_NS)

    @classmethod
    def _from_ns(cls, ns: int, error: float = 0.0) -> Date:
        """
        Create a Date from nanoseconds from the internal epoch.
        """
        date = cls(error=error)
        date._set_ns(ns)

        return date

    # -----------------------------------------------------------------------
    # Calendar fields
    # -----------------------------------------------------------------------

    @property
    def is_set(self) -> bool:
        """
        Whether the date holds an instant.
        """
        return self._ns is not None

    @property
    def year(self) -> Optional[int]:
        """
        Year, with astronomical numbering (year 0 is 1 BC).
        """
        return self._field(0)

    @year.setter
    def year(self, value: int) -> None:
        self._set_ns(self._ns_with(year=value))

    @property
    def month(self) -> Optional[int]:
        """
        Month, from 1 to 12.
        """
        return self._field(1)

    @month.setter
    def month(self, value: int) -> None:
        self._set_ns(self._ns_with(month=value))

    @property
    def day(self) -> Optional[int]:
        """
        Day of month.
        """
        return self._field(2)

    @day.setter
    def day(self, value: int) -> None:
        self._set_ns(self._ns_with(day=value))

    @property
    def hour(self) -> Optional[int]:
        """
        Hour, from 0 to 23.
        """
        return self._field(3)

    @hour.setter
    def hour(self, value: int) -> None:
        self._set_ns(self._ns_with(hour=value))

    @property
    def minute(self) -> Optional[int]:
        """
        Minute, from 0 to 59.
        """
        return self._field(4)

    @minute.setter
    def minute(self, value: int) -> None:
        self._set_ns(self._ns_with(minute=value))

    @property
    def second(self) -> Optional[float]:
        """
        Second of minute, including its fractional part, as a float.

        See :attr:`nanosecond` for the exact fractional part.
        """
        fields = self._get_fields()

        if fields is None:
            return None

        return (fields[5] * _NS_PER_SECOND + fields[6]) / _NS_PER_SECOND

    @second.setter
    def second(self, value: float) -> None:
        self._set_ns(self._ns_with(second=value))

    @property
    def nanosecond(self) -> Optional[int]:
        """
        Exact fractional part of the second, in nanoseconds (read-only).
        """
        return self._field(6)

    def replace(
        self,
        year: Any = _UNSET,
        month: Any = _UNSET,
        day: Any = _UNSET,
        hour: Any = _UNSET,
        minute: Any = _UNSET,
        second: Any = _UNSET,
        error: Any = _UNSET,
    ) -> Date:
        """
        Return a new date with some fields replaced.

        Parameters
        ----------
        year, month, day, hour, minute, second
            New values of the calendar fields. Omitted fields are kept;
            when ``second`` is omitted its fractional part is kept too.
        error
            New timing error. By default the current one is kept.

        Returns
        -------
        Date
            New date. The original date is not modified.

        Raises
        ------
        ValueError
            If the date is unset or the resulting date is invalid.
        """
        ns = self._ns_with(
            year=year,
            month=month,
            day=day,
            hour=hour,
            minute=minute,
            second=second,
        )

        return Date._from_ns(ns, self.error if error is _UNSET else error)

    def copy(self) -> Date:
        """
        Return an independent copy of the date.

        Returns
        -------
        Date
            Copied date, including the timing error. Copying an unset
            date returns an unset date.
        """
        if not self.is_set:
            return Date(error=self.error)

        return Date._from_ns(self._ns, self.error)

    @property
    def iso8601(self) -> str:
        """
        Date as an ISO 8601 UTC string with microsecond precision.
        """
        return self.to_iso8601()

    @property
    def ordinal_day(self) -> int:
        """
        Day of year.

        Notes
        -----
        This is the ordinal day within the year, not the astronomical
        Julian Date.
        """
        self._require_set()

        return to_ordinal_day(self.year, self.month, self.day)

    @ordinal_day.setter
    def ordinal_day(self, value: int) -> None:
        self._require_set()

        month, day = from_ordinal_day(self.year, value)
        self._set_ns(self._ns_with(month=month, day=day))

    @property
    def day_of_year(self) -> int:
        """
        Day of year. Alias of :attr:`ordinal_day`.
        """
        return self.ordinal_day

    @day_of_year.setter
    def day_of_year(self, value: int) -> None:
        self.ordinal_day = value

    @property
    def year_day(self) -> int:
        """
        Day of year. Alias of :attr:`ordinal_day`.
        """
        return self.ordinal_day

    @year_day.setter
    def year_day(self, value: int) -> None:
        self.ordinal_day = value

    @property
    def seconds(self) -> float:
        """
        Seconds from the internal epoch.
        """
        return self.to_seconds()

    @seconds.setter
    def seconds(self, value: TimeShift) -> None:
        self.from_seconds(value)

    # -----------------------------------------------------------------------
    # Setting and representing the date
    # -----------------------------------------------------------------------

    def set_date(self, date: Any) -> None:
        """
        Set the date from one of the supported input formats.

        Parameters
        ----------
        date
            Input date, in any format accepted by :class:`Date` except
            ``"now"``. The timing error is not changed.

        Raises
        ------
        TypeError
            If the input format is not supported.
        ValueError
            If the input does not describe a valid date, or if it is an
            unset Date.
        """
        if isinstance(date, Date):
            ns = date._require_ns()

        elif isinstance(date, (list, tuple)):
            ns = _sequence_to_ns(date)

        elif isinstance(date, dict):
            ns = _dict_to_ns(date)

        elif isinstance(date, str):
            ns, _ = _parse_iso8601(date)

        elif _is_time_shift(date):
            ns = _seconds_to_ns(date)

        else:
            raise TypeError("Unsupported date format.")

        self._set_ns(ns)

    def get_date(
        self,
        dtype: str = "list",
        format: str = "calendar",
    ) -> Union[List[Any], Dict[str, Any], str]:
        """
        Return the date in the requested representation.

        Parameters
        ----------
        dtype
            Output type: ``"list"``, ``"dict"``, ``"string"`` or
            ``"iso8601"`` (short forms ``"l"``, ``"d"``, ``"s"``,
            ``"str"`` are also accepted).
        format
            Date format: ``"calendar"`` or ``"ordinal"``.

        Returns
        -------
        list, dict or str
            Date representation.

        Raises
        ------
        ValueError
            If an option is not supported or the date is unset.
        """
        if format not in ("calendar", "ordinal"):
            raise ValueError("Unsupported date format.")

        self._require_set()

        if dtype in ("l", "list"):
            if format == "calendar":
                return [
                    self.year,
                    self.month,
                    self.day,
                    self.hour,
                    self.minute,
                    self.second,
                ]

            return [
                self.year,
                self.ordinal_day,
                self.hour,
                self.minute,
                self.second,
            ]

        if dtype in ("d", "dict"):
            return self.to_dict(format=format)

        if dtype in ("s", "str", "string", "iso8601"):
            return self.to_iso8601(format=format)

        raise ValueError("Unsupported output type.")

    def to_dict(self, format: str = "calendar") -> Dict[str, Any]:
        """
        Return the date as a dictionary.

        Parameters
        ----------
        format
            Either ``"calendar"`` or ``"ordinal"``.

        Returns
        -------
        dict
            Date dictionary.
        """
        self._require_set()

        if format == "calendar":
            return {
                "year": self.year,
                "month": self.month,
                "day": self.day,
                "hour": self.hour,
                "minute": self.minute,
                "second": self.second,
            }

        if format == "ordinal":
            return {
                "year": self.year,
                "ordinal_day": self.ordinal_day,
                "hour": self.hour,
                "minute": self.minute,
                "second": self.second,
            }

        raise ValueError("Unsupported date format.")

    def to_ordinal(self) -> List[Any]:
        """
        Return the date as ``[year, ordinal_day, hour, minute, second]``.
        """
        return self.get_date(dtype="list", format="ordinal")

    def to_iso8601(
        self,
        precision: int = 6,
        format: str = "calendar",
        timezone: str = "Z",
    ) -> str:
        """
        Return an ISO 8601 string.

        Parameters
        ----------
        precision
            Number of decimal digits for seconds. With ``0`` the
            fractional part is omitted; digits beyond the nanosecond
            are zeros.
        format
            Either ``"calendar"`` or ``"ordinal"``.
        timezone
            ``"Z"`` for UTC, or a numeric offset such as ``"+01:00"``,
            ``"-0430"`` or ``"+02"``. With an offset, the written time
            is the local time at that offset.

        Returns
        -------
        str
            ISO 8601 string.

        Notes
        -----
        Seconds are rounded exactly to ``precision`` digits (ties to
        even); the carry propagates to the following fields.
        """
        return _format_iso8601(
            self._require_ns(),
            precision,
            format,
            timezone,
        )

    def strformat(self, fmt: str, precision: int = 4) -> str:
        """
        Format the date with a restricted set of placeholders.

        Parameters
        ----------
        fmt
            Format string. Supported placeholders are ``%Y``, ``%j``,
            ``%m``, ``%d``, ``%H``, ``%M`` and ``%S``.
        precision
            Number of decimal digits for ``%S``. With ``0`` the
            fractional part is omitted.

        Returns
        -------
        str
            Formatted date.

        Notes
        -----
        Seconds are rounded exactly to ``precision`` digits (ties to
        even) before formatting; the carry propagates to the following
        fields.
        """
        ns = self._require_ns()
        precision = _validate_precision(precision)

        year, month, day, hour, minute, second, fraction = _ns_to_fields(
            _round_ns(ns, precision)
        )

        replacements = {
            "%Y": _format_year(year),
            "%j": f"{to_ordinal_day(year, month, day):03d}",
            "%m": f"{month:02d}",
            "%d": f"{day:02d}",
            "%H": f"{hour:02d}",
            "%M": f"{minute:02d}",
            "%S": _format_seconds(second, fraction, precision),
        }

        for key, value in replacements.items():
            fmt = fmt.replace(key, value)

        return fmt

    def strparse(self, date_string: str, fmt: str) -> None:
        """
        Set the date by parsing a string with a restricted format.

        Parameters
        ----------
        date_string
            String to parse.
        fmt
            Format string with fixed-width fields. Supported
            placeholders are ``%Y``, ``%j``, ``%m``, ``%d``, ``%H``,
            ``%M`` and ``%S`` (seconds may include a fractional part).
            ``%Y`` has four digits and an optional sign.

        Raises
        ------
        ValueError
            If the format is invalid (missing ``%Y``, repeated
            placeholders, or ``%j`` mixed with ``%m``/``%d``), or the
            string does not match it or does not describe a valid date.
        """
        _validate_parse_format(fmt)
        pattern = _format_to_regex(fmt)
        match = re.match(pattern, date_string)

        if match is None:
            raise ValueError("Date string does not match format.")

        values = match.groupdict()
        year = int(values["Y"])
        hour = int(values.get("H") or 0)
        minute = int(values.get("M") or 0)
        second = values.get("S") or "0"

        if values.get("j") is not None:
            month, day = from_ordinal_day(year, values["j"])
        else:
            month = int(values.get("m") or 1)
            day = int(values.get("d") or 1)

        self._set_ns(_fields_to_ns(year, month, day, hour, minute, second))

    # -----------------------------------------------------------------------
    # Conversions to and from seconds
    # -----------------------------------------------------------------------

    def to_seconds(self, decimal: bool = False) -> Union[float, Decimal]:
        """
        Return seconds from the internal epoch.

        Parameters
        ----------
        decimal
            If True, return an exact :class:`decimal.Decimal`; otherwise
            return a float.

        Returns
        -------
        float or Decimal
            Seconds from ``0001-01-01T00:00:00Z``.

        Notes
        -----
        Present-day values are about ``6.4e10`` s, so the float result
        resolves only about 8 microseconds. Use ``decimal=True`` for
        exact computations, or :meth:`to_unix` for a float with finer
        resolution.
        """
        ns = self._require_ns()

        if decimal:
            return _ns_to_decimal(ns)

        return ns / _NS_PER_SECOND

    def to_unix(self, decimal: bool = False) -> Union[float, Decimal]:
        """
        Return the Unix timestamp of the date.

        Parameters
        ----------
        decimal
            If True, return an exact :class:`decimal.Decimal`; otherwise
            return a float.

        Returns
        -------
        float or Decimal
            Seconds from 1970-01-01T00:00:00Z.

        Notes
        -----
        Present-day float timestamps resolve about 0.24 microseconds.
        """
        ns = self.to_unix_ns()

        if decimal:
            return _ns_to_decimal(ns)

        return ns / _NS_PER_SECOND

    def to_unix_ns(self) -> int:
        """
        Return the Unix timestamp of the date in integer nanoseconds.

        Returns
        -------
        int
            Nanoseconds from 1970-01-01T00:00:00Z.
        """
        return self._require_ns() - _UNIX_EPOCH_NS

    def from_seconds(self, second: TimeShift) -> None:
        """
        Set the date from seconds from the internal epoch.

        Parameters
        ----------
        second
            Seconds from ``0001-01-01T00:00:00Z``; negative values
            describe earlier instants.
        """
        self._set_ns(_seconds_to_ns(second))

    def shift_time(self, time: TimeShift, units: str = "s") -> None:
        """
        Shift the date in place by a time interval.

        Parameters
        ----------
        time
            Time shift, in the selected units.
        units
            Time units: seconds (``"s"``, ``"second"``, ``"seconds"``),
            minutes (``"m"``, ``"minute"``, ``"minutes"``), hours
            (``"h"``, ``"hour"``, ``"hours"``) or days (``"d"``,
            ``"day"``, ``"days"``).

        Raises
        ------
        ValueError
            If the units are not supported or the date is unset.
        """
        ns = self._require_ns()
        value = to_decimal(time)

        if units in ("s", "second", "seconds"):
            factor = Decimal(1)
        elif units in ("m", "minute", "minutes"):
            factor = MINUTE_SECONDS
        elif units in ("h", "hour", "hours"):
            factor = HOUR_SECONDS
        elif units in ("d", "day", "days"):
            factor = DAY_SECONDS
        else:
            raise ValueError("Unsupported time unit.")

        with localcontext(_DECIMAL_CONTEXT):
            seconds = value * factor

        self._set_ns(ns + _seconds_to_ns(seconds))

    # -----------------------------------------------------------------------
    # Internal state
    # -----------------------------------------------------------------------

    def _require_set(self) -> None:
        """
        Raise ValueError if the date is unset.
        """
        if self._ns is None:
            raise ValueError("Date is not set.")

    def _require_ns(self) -> int:
        """
        Return the internal nanoseconds, or raise if the date is unset.
        """
        self._require_set()

        return self._ns

    def _set_ns(self, ns: int) -> None:
        """
        Set the internal nanoseconds and invalidate the field cache.
        """
        self._ns = int(ns)
        self._fields = None

    def _get_fields(self) -> Optional[DateFields]:
        """
        Return the cached calendar fields, or None for an unset date.
        """
        if self._ns is None:
            return None

        if self._fields is None:
            self._fields = _ns_to_fields(self._ns)

        return self._fields

    def _field(self, index: int) -> Optional[int]:
        """
        Return one calendar field, or None for an unset date.
        """
        fields = self._get_fields()

        if fields is None:
            return None

        return fields[index]

    def _ns_with(
        self,
        year: Any = _UNSET,
        month: Any = _UNSET,
        day: Any = _UNSET,
        hour: Any = _UNSET,
        minute: Any = _UNSET,
        second: Any = _UNSET,
    ) -> int:
        """
        Return the nanoseconds of the date with some fields replaced.
        """
        self._require_set()

        (
            old_year,
            old_month,
            old_day,
            old_hour,
            old_minute,
            old_second,
            old_fraction,
        ) = self._get_fields()

        if second is _UNSET:
            second = _ns_to_decimal(
                old_second * _NS_PER_SECOND + old_fraction
            )

        return _fields_to_ns(
            old_year if year is _UNSET else year,
            old_month if month is _UNSET else month,
            old_day if day is _UNSET else day,
            old_hour if hour is _UNSET else hour,
            old_minute if minute is _UNSET else minute,
            second,
        )

    def _compare(self, value: Union[Date, TimeShift]) -> int:
        """
        Compare with a date or a number; return -1, 0 or 1.

        Numbers are compared exactly, as seconds from the internal
        epoch.
        """
        ns = self._require_ns()

        if isinstance(value, Date):
            other = value._require_ns()
            return (ns > other) - (ns < other)

        seconds = to_decimal(value)

        with localcontext(_DECIMAL_CONTEXT):
            other = seconds * _NS_DECIMAL

        mine = Decimal(ns)

        return (mine > other) - (mine < other)


# ---------------------------------------------------------------------------
# Calendar helpers
# ---------------------------------------------------------------------------

def leap_check(year: int) -> bool:
    """
    Return True if ``year`` is a leap year.

    Gregorian rules are applied to every year (proleptic calendar), so
    year ``0`` (1 BC) is a leap year.
    """
    year = _to_int(year, "year")

    divisible_by_4 = year % 4 == 0
    not_century = year % 100 != 0
    divisible_by_400 = year % 400 == 0

    return divisible_by_4 and (not_century or divisible_by_400)


def leap_num(year: int) -> int:
    """
    Return the number of leap years from year 1 to ``year - 1``.

    For years before 1 the result is negative: its magnitude is the
    number of leap years from ``year`` to ``0``. With this convention
    ``(year - 1) * 365 + leap_num(year)`` is always the number of days
    from the internal epoch to January 1 of ``year``.
    """
    year = _to_int(year, "year")

    # Floor division keeps the formula valid for non-positive years.
    previous = year - 1
    return previous // 4 - previous // 100 + previous // 400


def days_in_month(year: int) -> List[int]:
    """
    Return the cumulative number of days at the end of each month.
    """
    if leap_check(year):
        return [31, 60, 91, 121, 152, 182,
                213, 244, 274, 305, 335, 366]

    return [31, 59, 90, 120, 151, 181,
            212, 243, 273, 304, 334, 365]


def year_days(year: int) -> int:
    """
    Return the number of days in ``year``.
    """
    if leap_check(year):
        return 366

    return 365


def from_ordinal_day(year: int, day: int) -> Tuple[int, int]:
    """
    Convert an ordinal day of year to ``(month, day)``.

    Raises
    ------
    ValueError
        If the ordinal day is outside the year.
    """
    year = _to_int(year, "year")
    day = _to_int(day, "ordinal day")

    max_day = year_days(year)

    if day < 1 or day > max_day:
        raise ValueError("Ordinal day is outside the valid range.")

    month_days = days_in_month(year)

    for index, cumulative_day in enumerate(month_days):
        if day <= cumulative_day:
            month = index + 1
            previous = 0 if index == 0 else month_days[index - 1]
            return month, day - previous

    raise ValueError("Invalid ordinal day.")


def to_ordinal_day(year: int, month: int, day: int) -> int:
    """
    Convert ``(month, day)`` to the ordinal day of year.
    """
    year, month, day, _, _, _ = _validate_components(
        year, month, day, 0, 0, 0.0,
    )

    if month == 1:
        return day

    return days_in_month(year)[month - 2] + day


# ---------------------------------------------------------------------------
# Conversions to and from seconds
# ---------------------------------------------------------------------------

def date_to_second(
    year: int = 1,
    month: int = 1,
    day: int = 1,
    hour: int = 0,
    minute: int = 0,
    second: float = 0.0,
    decimal: bool = False,
) -> Union[float, Decimal]:
    """
    Convert a UTC date to seconds from the internal epoch.

    Parameters
    ----------
    year, month, day, hour, minute, second
        Calendar components.
    decimal
        If True, return a :class:`decimal.Decimal`; otherwise a float.

    Returns
    -------
    float or Decimal
        Seconds from ``0001-01-01T00:00:00Z``, negative for earlier
        instants.
    """
    year, month, day, hour, minute, second = _validate_components(
        year, month, day, hour, minute, second,
    )

    days = _days_from_date(year, month, day)

    with localcontext(_DECIMAL_CONTEXT):
        total = Decimal(days) * DAY_SECONDS
        total += Decimal(hour) * HOUR_SECONDS
        total += Decimal(minute) * MINUTE_SECONDS
        total += to_decimal(second)

    if decimal:
        return total

    return float(total)


def second_to_date(second: TimeShift) -> DateComponents:
    """
    Convert seconds from the internal epoch to UTC components.

    Parameters
    ----------
    second
        Seconds from ``0001-01-01T00:00:00Z``; negative values describe
        earlier instants.

    Returns
    -------
    list
        ``[year, month, day, hour, minute, second]``.
    """
    second = to_decimal(second)

    with localcontext(_DECIMAL_CONTEXT):
        # Decimal floor division truncates towards zero: correct it so
        # that the seconds of the day are always in [0, 86400).
        total_days = int(second // DAY_SECONDS)
        day_seconds = second - Decimal(total_days) * DAY_SECONDS

        if day_seconds < 0:
            total_days -= 1
            day_seconds += DAY_SECONDS

        hour = int(day_seconds // HOUR_SECONDS)
        day_seconds -= Decimal(hour) * HOUR_SECONDS

        minute = int(day_seconds // MINUTE_SECONDS)
        day_seconds -= Decimal(minute) * MINUTE_SECONDS

    year, month, day = _date_from_days(total_days)

    return [year, month, day, hour, minute, float(day_seconds)]


def to_decimal(value: Union[Date, TimeShift]) -> Decimal:
    """
    Convert a date or a time value to :class:`decimal.Decimal`.

    Parameters
    ----------
    value
        A :class:`Date` (converted to seconds from the internal epoch)
        or a real number in seconds.

    Returns
    -------
    Decimal
        Converted value. Floats are converted through their shortest
        decimal representation, so ``0.1`` becomes ``Decimal("0.1")``.

    Raises
    ------
    TypeError
        If the value is not a Date or a real number, or is a boolean.
    ValueError
        If the value is not finite.
    """
    if isinstance(value, Date):
        return value.to_seconds(decimal=True)

    if isinstance(value, bool):
        raise TypeError("Value cannot be converted to Decimal.")

    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, Integral):
        result = Decimal(int(value))
    elif isinstance(value, Real):
        result = Decimal(str(float(value)))
    else:
        raise TypeError("Value cannot be converted to Decimal.")

    if not result.is_finite():
        raise ValueError("Time value must be finite.")

    return result


# ---------------------------------------------------------------------------
# ISO 8601 and current time
# ---------------------------------------------------------------------------

def read_iso8601_date(date_string: str) -> Tuple[Any, ...]:
    """
    Parse an ISO 8601 date string and return UTC components.

    Supported examples are:

    - ``2026-05-09T10:30:12Z``
    - ``2026-05-09T10:30:12.345Z``
    - ``2026-129T10:30:12Z``
    - ``2026-05-09T12:30:12+02:00``

    A missing time zone designator is interpreted as UTC.

    Returns
    -------
    tuple
        ``(year, month, day, hour, minute, second, offset)``. The first
        six fields are converted to UTC, with seconds rounded to the
        nanosecond. The offset is expressed in seconds and corresponds
        to the input offset from UTC.

    Raises
    ------
    ValueError
        If the string is not a supported ISO 8601 date.
    """
    ns, offset = _parse_iso8601(date_string)
    year, month, day, hour, minute, second, fraction = _ns_to_fields(ns)
    seconds = (second * _NS_PER_SECOND + fraction) / _NS_PER_SECOND

    return (year, month, day, hour, minute, seconds, float(offset))


def write_iso8601_date(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    precision: int = 6,
    timezone: str = "Z",
    format: str = "calendar",
) -> str:
    """
    Convert date components to an ISO 8601 string.

    Parameters
    ----------
    year, month, day, hour, minute, second
        Calendar components, in UTC.
    precision
        Number of decimal digits for seconds. With ``0`` the fractional
        part is omitted; digits beyond the nanosecond are zeros.
    timezone
        ``"Z"`` for UTC, or a numeric offset such as ``"+01:00"``,
        ``"-0430"`` or ``"+02"``. With an offset the components are
        converted to the local time at that offset, and the offset is
        written as ``±HH:MM`` (``±HH:MM:SS`` if it has seconds).
    format
        Either ``"calendar"`` or ``"ordinal"``.

    Returns
    -------
    str
        ISO 8601 string.

    Raises
    ------
    TypeError
        If ``timezone`` is not a string or ``precision`` is not an
        integer.
    ValueError
        If a component, the precision, the offset or the format is
        invalid.

    Notes
    -----
    Seconds are first rounded to the nanosecond, then exactly to
    ``precision`` digits (ties to even); the carry propagates to the
    following fields.

    Years outside ``0``-``9999`` are written with an explicit sign and at
    least four digits (ISO 8601 expanded representation).
    """
    ns = _fields_to_ns(year, month, day, hour, minute, second)

    return _format_iso8601(ns, precision, format, timezone)


def utc_now() -> str:
    """
    Return the current UTC time as an ISO 8601 string.

    Returns
    -------
    str
        Current UTC time with microsecond precision.
    """
    return Date.now().to_iso8601(precision=6)


# ---------------------------------------------------------------------------
# Private helpers: nanosecond core
# ---------------------------------------------------------------------------

def _seconds_to_ns(value: Any) -> int:
    """
    Convert seconds to integer nanoseconds, rounding ties to even.

    ``value`` may be a real number, a Decimal or a numeric string.
    """
    if isinstance(value, Integral) and not isinstance(value, bool):
        return int(value) * _NS_PER_SECOND

    if isinstance(value, str):
        try:
            seconds = Decimal(value)
        except InvalidOperation:
            raise ValueError("Invalid number of seconds.") from None

        if not seconds.is_finite():
            raise ValueError("Time value must be finite.")
    else:
        seconds = to_decimal(value)

    with localcontext(_DECIMAL_CONTEXT):
        nanoseconds = seconds * _NS_DECIMAL
        return int(nanoseconds.to_integral_value(rounding=ROUND_HALF_EVEN))


def _ns_to_decimal(ns: int) -> Decimal:
    """
    Convert integer nanoseconds to exact Decimal seconds.
    """
    with localcontext(_DECIMAL_CONTEXT):
        return Decimal(ns).scaleb(-9)


def _days_from_date(year: int, month: int, day: int) -> int:
    """
    Return the days from the internal epoch to a validated date.
    """
    month_days = [0] + days_in_month(year)[:-1]

    return _days_before_year(year) + month_days[month - 1] + day - 1


def _date_from_days(day_count: int) -> Tuple[int, int, int]:
    """
    Return ``(year, month, day)`` for a day counted from the epoch.
    """
    year = _year_from_day_count(day_count)
    month, day = from_ordinal_day(
        year,
        day_count - _days_before_year(year) + 1,
    )

    return year, month, day


def _fields_to_ns(
    year: Any,
    month: Any,
    day: Any,
    hour: Any,
    minute: Any,
    second: Any,
) -> int:
    """
    Validate calendar fields and convert them to nanoseconds.

    Seconds are rounded to the nanosecond; a value that rounds to 60 s
    carries into the next minute.
    """
    year, month, day, hour, minute, _ = _validate_components(
        year, month, day, hour, minute, second,
    )

    return (
        _days_from_date(year, month, day) * _NS_PER_DAY
        + hour * _NS_PER_HOUR
        + minute * _NS_PER_MINUTE
        + _seconds_to_ns(second)
    )


def _ns_to_fields(ns: int) -> DateFields:
    """
    Convert nanoseconds to ``(year, month, day, hour, minute, second,
    nanosecond)``, with integer seconds and nanoseconds.
    """
    days, remainder = divmod(ns, _NS_PER_DAY)
    year, month, day = _date_from_days(days)

    hour, remainder = divmod(remainder, _NS_PER_HOUR)
    minute, remainder = divmod(remainder, _NS_PER_MINUTE)
    second, fraction = divmod(remainder, _NS_PER_SECOND)

    return year, month, day, hour, minute, second, fraction


def _sequence_to_ns(values: Sequence[Any]) -> int:
    """
    Convert a calendar or ordinal sequence to nanoseconds.
    """
    if len(values) == 6:
        return _fields_to_ns(*values)

    if len(values) == 5:
        year = _to_int(values[0], "year")
        month, day = from_ordinal_day(year, values[1])
        return _fields_to_ns(
            year,
            month,
            day,
            values[2],
            values[3],
            values[4],
        )

    raise ValueError("Sequence must have 5 or 6 elements.")


def _dict_to_ns(date: Dict[str, Any]) -> int:
    """
    Convert a calendar or ordinal dictionary to nanoseconds.
    """
    year = _to_int(date["year"], "year")
    hour = date.get("hour", 0)
    minute = date.get("minute", 0)
    second = date.get("second", 0.0)

    if "month" in date and "day" in date:
        month = date["month"]
        day = date["day"]
    elif "ordinal_day" in date:
        month, day = from_ordinal_day(year, date["ordinal_day"])
    elif "day_of_year" in date:
        month, day = from_ordinal_day(year, date["day_of_year"])
    elif "year_day" in date:
        month, day = from_ordinal_day(year, date["year_day"])
    else:
        raise ValueError("Dictionary lacks calendar or ordinal fields.")

    return _fields_to_ns(year, month, day, hour, minute, second)


def _parse_iso8601(date_string: str) -> Tuple[int, Decimal]:
    """
    Parse an ISO 8601 string into UTC nanoseconds and the input offset.
    """
    match = ISO8601_PATTERN.match(date_string)

    if match is None:
        raise ValueError("Invalid ISO 8601 date string.")

    values = match.groupdict()
    year = int(values["year"])

    if values["ordinal"] is not None:
        month, day = from_ordinal_day(year, values["ordinal"])
    else:
        month = int(values["month"])
        day = int(values["day"])

    local_ns = _fields_to_ns(
        year,
        month,
        day,
        int(values["hour"]),
        int(values["minute"]),
        values["second"],
    )

    offset = _parse_timezone_offset(values["timezone"] or "Z")

    return local_ns - int(offset) * _NS_PER_SECOND, offset


def _format_iso8601(
    ns: int,
    precision: int,
    format: str,
    timezone: str,
) -> str:
    """
    Format nanoseconds from the internal epoch as ISO 8601.
    """
    precision = _validate_precision(precision)

    if not isinstance(timezone, str):
        raise TypeError("The timezone must be a string.")

    offset = _parse_timezone_offset(timezone)

    if timezone == "Z":
        designator = "Z"
    else:
        designator = _format_timezone_offset(offset)

    local_ns = _round_ns(ns + int(offset) * _NS_PER_SECOND, precision)
    year, month, day, hour, minute, second, fraction = _ns_to_fields(
        local_ns
    )

    if format == "calendar":
        date_part = f"{_format_year(year)}-{month:02d}-{day:02d}"
    elif format == "ordinal":
        ordinal = to_ordinal_day(year, month, day)
        date_part = f"{_format_year(year)}-{ordinal:03d}"
    else:
        raise ValueError("Unsupported ISO 8601 date format.")

    time_part = (
        f"{hour:02d}:{minute:02d}:"
        f"{_format_seconds(second, fraction, precision)}"
    )

    return f"{date_part}T{time_part}{designator}"


def _round_ns(ns: int, precision: int) -> int:
    """
    Round nanoseconds to ``precision`` decimal digits of a second.

    Ties are rounded to even. Precisions of 9 or more leave the value
    unchanged.
    """
    if precision >= 9:
        return ns

    unit = 10 ** (9 - precision)
    quotient, remainder = divmod(ns, unit)
    twice = 2 * remainder

    if twice > unit or (twice == unit and quotient % 2 == 1):
        quotient += 1

    return quotient * unit


def _format_seconds(second: int, fraction: int, precision: int) -> str:
    """
    Format integer seconds and nanoseconds with ``precision`` digits.

    The value must already be rounded with :func:`_round_ns`.
    """
    if precision == 0:
        return f"{second:02d}"

    if precision <= 9:
        digits = fraction // 10 ** (9 - precision)
        return f"{second:02d}.{digits:0{precision}d}"

    return f"{second:02d}.{fraction:09d}{'0' * (precision - 9)}"


# ---------------------------------------------------------------------------
# Private helpers: validation and parsing
# ---------------------------------------------------------------------------

def _is_time_shift(value: Any) -> bool:
    """
    Return True if ``value`` is a real number usable as a time shift.
    """
    if isinstance(value, bool):
        return False

    return isinstance(value, (Decimal, Real))


def _is_time_value(value: Any) -> bool:
    """
    Return True if ``value`` can be compared with a Date.
    """
    return isinstance(value, Date) or _is_time_shift(value)


def _is_finite_number(value: Union[Decimal, Real]) -> bool:
    """
    Return True if a real number or Decimal is finite.
    """
    if isinstance(value, Decimal):
        return value.is_finite()

    return math.isfinite(float(value))


def _to_int(value: Any, name: str) -> int:
    """
    Convert an integral value to int, rejecting silent truncation.

    Accepted values are integers (including NumPy integers),
    integral-valued floats and Decimals (e.g. ``2026.0``) and strings
    holding an integer (e.g. ``"2026"``).
    """
    if isinstance(value, bool):
        raise TypeError(f"The {name} must be an integer.")

    if isinstance(value, Integral):
        return int(value)

    if isinstance(value, (Decimal, Real)):
        if not _is_finite_number(value) or value != int(value):
            raise ValueError(f"The {name} must be an integer.")

        return int(value)

    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            raise ValueError(f"The {name} must be an integer.") from None

    raise TypeError(f"The {name} must be an integer.")


def _to_float(value: Any, name: str) -> float:
    """
    Convert a real value (or a numeric string) to a finite float.
    """
    if isinstance(value, bool):
        raise TypeError(f"The {name} must be a real number.")

    if isinstance(value, (Decimal, Real)):
        result = float(value)
    elif isinstance(value, str):
        try:
            result = float(value)
        except ValueError:
            raise ValueError(
                f"The {name} must be a real number."
            ) from None
    else:
        raise TypeError(f"The {name} must be a real number.")

    if not math.isfinite(result):
        raise ValueError(f"The {name} must be finite.")

    return result


def _validate_error(error: Any) -> float:
    """
    Validate a timing error (finite, non-negative real number).
    """
    if isinstance(error, bool) or not isinstance(error, (Decimal, Real)):
        raise TypeError("The timing error must be a real number.")

    error = float(error)

    if not math.isfinite(error) or error < 0:
        raise ValueError("The timing error must be finite and non-negative.")

    return error


def _validate_components(
    year: Any,
    month: Any,
    day: Any,
    hour: Any,
    minute: Any,
    second: Any,
) -> Tuple[int, int, int, int, int, float]:
    """
    Validate date components and return them as ints and a float.
    """
    year = _to_int(year, "year")
    month = _to_int(month, "month")
    day = _to_int(day, "day")
    hour = _to_int(hour, "hour")
    minute = _to_int(minute, "minute")
    second = _to_float(second, "second")

    if month < 1 or month > 12:
        raise ValueError("Month must be between 1 and 12.")

    month_days = days_in_month(year)
    previous = 0 if month == 1 else month_days[month - 2]
    max_day = month_days[month - 1] - previous

    if day < 1 or day > max_day:
        raise ValueError("Invalid day for the selected month.")

    if hour < 0 or hour > 23:
        raise ValueError("Hour must be between 0 and 23.")

    if minute < 0 or minute > 59:
        raise ValueError("Minute must be between 0 and 59.")

    if not 0 <= second < 60:
        raise ValueError("Second must be in the range [0, 60).")

    return year, month, day, hour, minute, second


def _validate_precision(precision: int) -> int:
    """
    Validate a number of decimal digits for seconds.
    """
    if isinstance(precision, bool) or not isinstance(precision, Integral):
        raise TypeError("The precision must be an integer.")

    precision = int(precision)

    if precision < 0:
        raise ValueError("The precision must be non-negative.")

    return precision


def _validate_parse_format(fmt: str) -> None:
    """
    Check that a strparse format is unambiguous.
    """
    if "%Y" not in fmt:
        raise ValueError("The format must contain %Y.")

    for key in ("%Y", "%j", "%m", "%d", "%H", "%M", "%S"):
        if fmt.count(key) > 1:
            raise ValueError(f"Placeholder {key} appears more than once.")

    if "%j" in fmt and ("%m" in fmt or "%d" in fmt):
        raise ValueError("Use either %j or %m and %d, not both.")


def _format_year(year: int) -> str:
    """
    Format a year as in ISO 8601, with a sign outside ``0``-``9999``.
    """
    if 0 <= year <= 9999:
        return f"{year:04d}"

    return f"{year:+05d}"


def _days_before_year(year: int) -> int:
    """
    Return the number of days before January 1 of ``year``.
    """
    return (year - 1) * YEAR_DAYS + leap_num(year)


def _year_from_day_count(day_count: int) -> int:
    """
    Return the year containing a day counted from the internal epoch.

    ``day_count`` may be negative for days before the epoch.
    """
    # A Gregorian cycle has 146097 days in 400 years: the estimate is
    # off by at most one year and is corrected below.
    year = (day_count * 400) // 146097 + 1

    while _days_before_year(year) > day_count:
        year -= 1

    while _days_before_year(year + 1) <= day_count:
        year += 1

    return year


def _parse_timezone_offset(timezone: Optional[str]) -> Decimal:
    """
    Parse an ISO 8601 time zone offset.

    Returns
    -------
    Decimal
        Offset from UTC in seconds.
    """
    if timezone in ("Z", None):
        return Decimal("0")

    sign = timezone[0]

    if sign not in ("+", "-"):
        raise ValueError("Invalid timezone offset.")

    value = timezone[1:].replace(":", "")

    if len(value) not in (2, 4, 6) or not value.isdigit():
        raise ValueError("Invalid timezone offset.")

    hour = int(value[0:2])
    minute = int(value[2:4] or 0)
    second = int(value[4:6] or 0)

    if hour > 23 or minute > 59 or second > 59:
        raise ValueError("Invalid timezone offset.")

    with localcontext(_DECIMAL_CONTEXT):
        offset = Decimal(hour) * HOUR_SECONDS
        offset += Decimal(minute) * MINUTE_SECONDS
        offset += Decimal(second)

        if sign == "-":
            offset = -offset

    return offset


def _format_timezone_offset(offset: Decimal) -> str:
    """
    Format an offset in seconds as ``±HH:MM`` or ``±HH:MM:SS``.
    """
    sign = "-" if offset < 0 else "+"
    total = int(abs(offset))
    hour, remainder = divmod(total, 3600)
    minute, second = divmod(remainder, 60)

    if second:
        return f"{sign}{hour:02d}:{minute:02d}:{second:02d}"

    return f"{sign}{hour:02d}:{minute:02d}"


def _format_to_regex(fmt: str) -> str:
    """
    Convert a restricted date format to a regular expression.
    """
    fields = {
        "%Y": r"(?P<Y>[+-]?[0-9]{4})",
        "%j": r"(?P<j>[0-9]{3})",
        "%m": r"(?P<m>[0-9]{2})",
        "%d": r"(?P<d>[0-9]{2})",
        "%H": r"(?P<H>[0-9]{2})",
        "%M": r"(?P<M>[0-9]{2})",
        "%S": r"(?P<S>[0-9]{2}(?:\.[0-9]+)?)",
    }

    pattern = re.escape(fmt)

    for key, value in fields.items():
        pattern = pattern.replace(re.escape(key), value)

    return f"^{pattern}$"
