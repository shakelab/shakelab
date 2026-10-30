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
- Comparisons and arithmetic are evaluated with
  :class:`decimal.Decimal`, so that time shifts are not affected by
  binary floating-point rounding of large epoch offsets.
- ISO 8601 strings with an explicit offset are accepted, but they are
  converted to UTC; the offset itself is not retained.

Conventions
-----------
- The internal epoch is ``0001-01-01T00:00:00Z``. This is not Unix time.
- Leap seconds are not represented: seconds lie in ``[0, 60)``.
- Time shifts are real numbers expressed in seconds: Python ``int`` and
  ``float``, :class:`decimal.Decimal`, or NumPy integer and floating
  scalars. Booleans and non-finite values are rejected.
"""

from __future__ import annotations

from decimal import Decimal
from numbers import Integral, Real
import re
import time as systime
from typing import Any, Dict, List, Sequence, Tuple, Union


__all__ = [
    "Date",
    "MINUTE_SECONDS",
    "HOUR_SECONDS",
    "DAY_SECONDS",
    "YEAR_DAYS",
    "MSEC",
    "HSEC",
    "DSEC",
    "YDAYS",
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


# Deprecated aliases kept temporarily for backward compatibility.
# Prefer MINUTE_SECONDS, HOUR_SECONDS, DAY_SECONDS and YEAR_DAYS.
MSEC = MINUTE_SECONDS
HSEC = HOUR_SECONDS
DSEC = DAY_SECONDS
YDAYS = YEAR_DAYS


ISO8601_PATTERN = re.compile(
    r"^"
    r"(?P<year>[0-9]{4})"
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
        Optional timing error, in seconds.
    timezone
        Only ``"UTC"`` is supported.

    Raises
    ------
    TypeError
        If the input format is not supported.
    ValueError
        If the input does not describe a valid date.
    NotImplementedError
        If a time zone other than UTC is requested.

    Notes
    -----
    An unset date can be filled later, for example with
    :meth:`strparse`. Operations that need an instant (comparisons,
    arithmetic, hashing, conversion to seconds or strings) raise
    ``ValueError`` on an unset date.

    Dates are hashable and dates that compare equal have the same hash.
    Date fields are mutable: do not modify a date while it is used as a
    dictionary key or set member.
    """

    # Make NumPy scalars defer binary operators to Date (NEP 13), so that
    # expressions such as ``np.int64(5) + date`` call Date.__radd__.
    __array_ufunc__ = None

    def __init__(
        self,
        date: Any = None,
        error: float = 0.0,
        timezone: str = "UTC",
    ) -> None:
        self.year = None
        self.month = None
        self.day = None
        self.hour = None
        self.minute = None
        self.second = None
        self.error = float(error)

        if timezone != "UTC":
            raise NotImplementedError("Only UTC timezone is supported.")

        if date is None:
            return

        if isinstance(date, str) and date == "now":
            self._set_from_sequence(_utc_now_components())
            self._selfcheck()
            return

        self.set_date(date)

    def __repr__(self) -> str:
        if not self.is_set:
            return "Date(None)"

        return self.iso8601

    def __eq__(self, value: object) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self.to_seconds(decimal=True) == to_decimal(value)

    def __hash__(self) -> int:
        return hash(self.to_seconds(decimal=True))

    def __lt__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self.to_seconds(decimal=True) < to_decimal(value)

    def __le__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self.to_seconds(decimal=True) <= to_decimal(value)

    def __gt__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self.to_seconds(decimal=True) > to_decimal(value)

    def __ge__(self, value: Any) -> bool:
        if not _is_time_value(value):
            return NotImplemented

        return self.to_seconds(decimal=True) >= to_decimal(value)

    def __add__(self, value: TimeShift) -> Date:
        if not _is_time_shift(value):
            return NotImplemented

        return Date(self.to_seconds(decimal=True) + to_decimal(value))

    __radd__ = __add__

    def __sub__(self, value: Union[Date, TimeShift]) -> Union[Date, float]:
        if isinstance(value, Date):
            seconds = self.to_seconds(decimal=True)
            seconds -= value.to_seconds(decimal=True)
            return float(seconds)

        if _is_time_shift(value):
            return Date(self.to_seconds(decimal=True) - to_decimal(value))

        return NotImplemented

    @classmethod
    def now(cls) -> Date:
        """
        Return the current UTC time.

        Returns
        -------
        Date
            Current UTC date, with the sub-second resolution of the
            system clock.
        """
        return cls(_utc_now_components())

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

    @property
    def is_set(self) -> bool:
        """
        Whether the date holds an instant.
        """
        return None not in (
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        )

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

        return Date(self.to_seconds(decimal=True), error=self.error)

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
        self.month = month
        self.day = day
        self._selfcheck()

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

    def set_date(self, date: Any) -> None:
        """
        Set the date from one of the supported input formats.

        Parameters
        ----------
        date
            Input date, in any format accepted by :class:`Date` except
            ``"now"``.

        Raises
        ------
        TypeError
            If the input format is not supported.
        ValueError
            If the input does not describe a valid date, or if it is an
            unset Date.
        """
        if isinstance(date, Date):
            date._require_set()
            self._set_from_sequence(date.get_date())

        elif isinstance(date, (list, tuple)):
            self._set_from_sequence(date)

        elif isinstance(date, dict):
            self._set_from_dict(date)

        elif isinstance(date, str):
            values = read_iso8601_date(date)
            self._set_from_sequence(values[:6])

        elif _is_time_shift(date):
            self.from_seconds(date)
            return

        else:
            raise TypeError("Unsupported date format.")

        self._selfcheck()

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
    ) -> str:
        """
        Return an ISO 8601 UTC string.

        Parameters
        ----------
        precision
            Number of decimal digits for seconds. With ``0`` the
            fractional part is omitted.
        format
            Either ``"calendar"`` or ``"ordinal"``.

        Returns
        -------
        str
            ISO 8601 UTC string.
        """
        self._require_set()

        return write_iso8601_date(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            precision=precision,
            format=format,
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
        Seconds are rounded to ``precision`` digits before formatting.
        When rounding reaches 60 seconds, the carry propagates to the
        minute and, if needed, to the following fields.
        """
        self._require_set()
        precision = _validate_precision(precision)

        year, month, day, hour, minute, second = _round_components(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            precision,
        )

        replacements = {
            "%Y": f"{year:04d}",
            "%j": f"{to_ordinal_day(year, month, day):03d}",
            "%m": f"{month:02d}",
            "%d": f"{day:02d}",
            "%H": f"{hour:02d}",
            "%M": f"{minute:02d}",
            "%S": _format_seconds(second, precision),
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

        Raises
        ------
        ValueError
            If the string does not match the format or does not describe
            a valid date.
        """
        pattern = _format_to_regex(fmt)
        match = re.match(pattern, date_string)

        if match is None:
            raise ValueError("Date string does not match format.")

        values = match.groupdict()
        year = int(values["Y"])
        hour = int(values.get("H") or 0)
        minute = int(values.get("M") or 0)
        second = _decimal_to_float(values.get("S") or "0")

        if values.get("j") is not None:
            month, day = from_ordinal_day(year, int(values["j"]))
        else:
            month = int(values.get("m") or 1)
            day = int(values.get("d") or 1)

        self._set_components(year, month, day, hour, minute, second)
        self._selfcheck()

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
        """
        self._require_set()

        return date_to_second(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            decimal=decimal,
        )

    def from_seconds(self, second: TimeShift) -> None:
        """
        Set the date from seconds from the internal epoch.

        Parameters
        ----------
        second
            Non-negative seconds from ``0001-01-01T00:00:00Z``.
        """
        values = second_to_date(second)
        self._set_from_sequence(values)
        self._selfcheck()

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
        self._require_set()
        value = to_decimal(time)

        if units in ("s", "second", "seconds"):
            seconds = value
        elif units in ("m", "minute", "minutes"):
            seconds = value * MINUTE_SECONDS
        elif units in ("h", "hour", "hours"):
            seconds = value * HOUR_SECONDS
        elif units in ("d", "day", "days"):
            seconds = value * DAY_SECONDS
        else:
            raise ValueError("Unsupported time unit.")

        self.from_seconds(self.to_seconds(decimal=True) + seconds)

    def _require_set(self) -> None:
        """
        Raise ValueError if the date is unset.
        """
        if not self.is_set:
            raise ValueError("Date is not set.")

    def _set_from_sequence(self, values: Sequence[Any]) -> None:
        """
        Set date fields from a calendar or ordinal sequence.
        """
        if len(values) == 6:
            self._set_components(*values)
            return

        if len(values) == 5:
            year = int(values[0])
            month, day = from_ordinal_day(year, int(values[1]))
            self._set_components(
                year,
                month,
                day,
                values[2],
                values[3],
                values[4],
            )
            return

        raise ValueError("Sequence must have 5 or 6 elements.")

    def _set_from_dict(self, date: Dict[str, Any]) -> None:
        """
        Set date fields from a dictionary.
        """
        year = int(date["year"])
        hour = int(date.get("hour", 0))
        minute = int(date.get("minute", 0))
        second = float(date.get("second", 0.0))

        if "month" in date and "day" in date:
            month = int(date["month"])
            day = int(date["day"])
        elif "ordinal_day" in date:
            month, day = from_ordinal_day(year, int(date["ordinal_day"]))
        elif "day_of_year" in date:
            month, day = from_ordinal_day(year, int(date["day_of_year"]))
        elif "year_day" in date:
            month, day = from_ordinal_day(year, int(date["year_day"]))
        else:
            raise ValueError("Dictionary lacks calendar or ordinal fields.")

        self._set_components(year, month, day, hour, minute, second)

    def _set_components(
        self,
        year: int,
        month: int,
        day: int,
        hour: int,
        minute: int,
        second: float,
    ) -> None:
        """
        Set raw date components without validation.
        """
        self.year = int(year)
        self.month = int(month)
        self.day = int(day)
        self.hour = int(hour)
        self.minute = int(minute)
        self.second = float(second)

    def _selfcheck(self) -> None:
        """
        Validate date fields.
        """
        _validate_components(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        )


# ---------------------------------------------------------------------------
# Calendar helpers
# ---------------------------------------------------------------------------

def leap_check(year: int) -> bool:
    """
    Return True if ``year`` is a leap year (Gregorian rules).
    """
    year = int(year)

    if year < 1:
        raise ValueError("Year must be greater than or equal to 1.")

    divisible_by_4 = year % 4 == 0
    not_century = year % 100 != 0
    divisible_by_400 = year % 400 == 0

    return divisible_by_4 and (not_century or divisible_by_400)


def leap_num(year: int) -> int:
    """
    Return the number of leap years before ``year``.
    """
    year = int(year)

    if year < 1:
        raise ValueError("Year must be greater than or equal to 1.")

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
    year = int(year)
    day = int(day)

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
    year = int(year)
    month = int(month)
    day = int(day)

    _validate_components(year, month, day, 0, 0, 0.0)

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
        Seconds from ``0001-01-01T00:00:00Z``.
    """
    _validate_components(year, month, day, hour, minute, second)

    month_days = [0] + days_in_month(year)[:-1]
    days = (year - 1) * YEAR_DAYS
    days += leap_num(year)
    days += month_days[month - 1]
    days += day - 1

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
        Non-negative seconds from ``0001-01-01T00:00:00Z``.

    Returns
    -------
    list
        ``[year, month, day, hour, minute, second]``.
    """
    second = to_decimal(second)

    if second < 0:
        raise ValueError("Seconds must be greater than or equal to zero.")

    total_days = int(second // DAY_SECONDS)
    day_seconds = second - Decimal(total_days) * DAY_SECONDS

    year = _year_from_day_count(total_days)
    days_before_year = _days_before_year(year)
    ordinal = total_days - days_before_year + 1

    month, day = from_ordinal_day(year, ordinal)

    hour = int(day_seconds // HOUR_SECONDS)
    day_seconds -= Decimal(hour) * HOUR_SECONDS

    minute = int(day_seconds // MINUTE_SECONDS)
    day_seconds -= Decimal(minute) * MINUTE_SECONDS

    return [
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        float(day_seconds),
    ]


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
        six fields are converted to UTC. The offset is expressed in
        seconds and corresponds to the input offset from UTC.

    Raises
    ------
    ValueError
        If the string is not a supported ISO 8601 date.
    """
    match = ISO8601_PATTERN.match(date_string)

    if match is None:
        raise ValueError("Invalid ISO 8601 date string.")

    values = match.groupdict()
    year = int(values["year"])

    if values["ordinal"] is not None:
        month, day = from_ordinal_day(year, int(values["ordinal"]))
    else:
        month = int(values["month"])
        day = int(values["day"])

    hour = int(values["hour"])
    minute = int(values["minute"])
    second = float(values["second"])
    timezone = values["timezone"] or "Z"

    _validate_components(year, month, day, hour, minute, second)

    offset = _parse_timezone_offset(timezone)
    absolute_seconds = date_to_second(
        year,
        month,
        day,
        hour,
        minute,
        second,
        decimal=True,
    )

    absolute_seconds -= offset
    utc_date = second_to_date(absolute_seconds)

    return tuple(utc_date) + (float(offset),)


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
    Convert date components to an ISO 8601 UTC string.

    Parameters
    ----------
    year, month, day, hour, minute, second
        Calendar components.
    precision
        Number of decimal digits for seconds. With ``0`` the fractional
        part is omitted.
    timezone
        Time zone designator appended to the string.
    format
        Either ``"calendar"`` or ``"ordinal"``.

    Returns
    -------
    str
        ISO 8601 string.

    Notes
    -----
    Seconds are rounded to ``precision`` digits. When rounding reaches
    60 seconds, the carry propagates to the following fields.
    """
    _validate_components(year, month, day, hour, minute, second)
    precision = _validate_precision(precision)

    year, month, day, hour, minute, second = _round_components(
        year,
        month,
        day,
        hour,
        minute,
        second,
        precision,
    )

    if format == "calendar":
        date_part = f"{year:04d}-{month:02d}-{day:02d}"
    elif format == "ordinal":
        ordinal = to_ordinal_day(year, month, day)
        date_part = f"{year:04d}-{ordinal:03d}"
    else:
        raise ValueError("Unsupported ISO 8601 date format.")

    time_part = (
        f"{hour:02d}:{minute:02d}:"
        f"{_format_seconds(second, precision)}"
    )

    return f"{date_part}T{time_part}{timezone}"


def utc_now() -> str:
    """
    Return the current UTC time as an ISO 8601 string.

    Returns
    -------
    str
        Current UTC time with microsecond precision.
    """
    return write_iso8601_date(*_utc_now_components(), precision=6)


# ---------------------------------------------------------------------------
# Private helpers
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


def _utc_now_components() -> DateComponents:
    """
    Return the current UTC time as calendar components.

    Seconds keep the nanosecond resolution of :func:`time.time_ns`.
    """
    seconds, nanoseconds = divmod(systime.time_ns(), 1_000_000_000)
    current = systime.gmtime(seconds)

    return [
        current.tm_year,
        current.tm_mon,
        current.tm_mday,
        current.tm_hour,
        current.tm_min,
        current.tm_sec + nanoseconds / 1e9,
    ]


def _validate_components(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
) -> None:
    """
    Validate raw date components.
    """
    year = int(year)
    month = int(month)
    day = int(day)
    hour = int(hour)
    minute = int(minute)
    second = float(second)

    if year < 1:
        raise ValueError("Year must be greater than or equal to 1.")

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

    if second < 0 or second >= 60:
        raise ValueError("Second must be in the range [0, 60).")


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


def _round_components(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    precision: int,
) -> Tuple[int, int, int, int, int, float]:
    """
    Round seconds to ``precision`` digits, carrying into the minute.

    When the rounded seconds reach 60, the date is advanced to the next
    minute (and, if needed, hour, day, month and year) and the seconds
    are set to zero.
    """
    rounded = round(float(second), precision)

    if rounded < 60:
        return year, month, day, hour, minute, rounded

    next_minute = date_to_second(
        year,
        month,
        day,
        hour,
        minute,
        0.0,
        decimal=True,
    )
    next_minute += MINUTE_SECONDS

    year, month, day, hour, minute, _ = second_to_date(next_minute)

    return year, month, day, hour, minute, 0.0


def _format_seconds(second: float, precision: int) -> str:
    """
    Format already rounded seconds with a two-digit integer part.
    """
    width = 2 if precision == 0 else precision + 3

    return f"{second:0{width}.{precision}f}"


def _days_before_year(year: int) -> int:
    """
    Return the number of days before January 1 of ``year``.
    """
    return (year - 1) * YEAR_DAYS + leap_num(year)


def _year_from_day_count(day_count: int) -> int:
    """
    Return the year from the number of elapsed days since the epoch.
    """
    low = 1
    high = max(2, int(day_count // 365) + 2)

    while _days_before_year(high) <= day_count:
        high *= 2

    while low < high:
        middle = (low + high) // 2

        if _days_before_year(middle + 1) <= day_count:
            low = middle + 1
        else:
            high = middle

    return low


def _parse_timezone_offset(timezone: str) -> Decimal:
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

    if len(value) not in (2, 4, 6):
        raise ValueError("Invalid timezone offset.")

    hour = int(value[0:2])
    minute = int(value[2:4] or 0)
    second = int(value[4:6] or 0)

    if hour > 23 or minute > 59 or second > 59:
        raise ValueError("Invalid timezone offset.")

    offset = Decimal(hour) * HOUR_SECONDS
    offset += Decimal(minute) * MINUTE_SECONDS
    offset += Decimal(second)

    if sign == "-":
        offset = -offset

    return offset


def _format_to_regex(fmt: str) -> str:
    """
    Convert a restricted date format to a regular expression.
    """
    fields = {
        "%Y": r"(?P<Y>[0-9]{4})",
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


def _decimal_to_float(value: str) -> float:
    """
    Convert a decimal string to float without binary parsing artefacts.
    """
    return float(Decimal(str(value)))
