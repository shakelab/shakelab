# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Time and date utilities.

The Date class stores an absolute UTC instant. ISO 8601 strings with an
explicit offset are accepted, but they are converted internally to UTC.

The internal epoch is 0001-01-01T00:00:00Z.

The implementation intentionally does not use the standard datetime module.
"""

from decimal import Decimal
import re
import time as systime


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


class Date(object):
    """
    Represent an absolute UTC date and time.

    Parameters
    ----------
    date : str, int, float, Decimal, list, tuple, dict, Date, optional
        Input date. Supported formats are ISO 8601 strings, seconds from
        the internal epoch, calendar tuples, ordinal tuples, dictionaries,
        or another Date object.
    error : float, optional
        Optional timing error, in seconds.
    timezone : str, optional
        Currently only ``"UTC"`` is supported.

    Notes
    -----
    The internal epoch is ``0001-01-01T00:00:00Z``. This is not Unix time.
    """

    def __init__(self, date=None, error=0.0, timezone="UTC"):
        self.year = None
        self.month = None
        self.day = None
        self.hour = None
        self.minute = None
        self.second = None
        self.error = float(error)

        if timezone != "UTC":
            raise NotImplementedError("Only UTC timezone is supported.")

        if date is not None:
            if isinstance(date, str) and date == "now":
                self.set_date(utc_now())
            else:
                self.set_date(date)

    def __repr__(self):
        return self.iso8601

    def __eq__(self, value):
        try:
            return self.to_seconds(decimal=True) == to_decimal(value)
        except TypeError:
            return False

    def __lt__(self, value):
        return self.to_seconds(decimal=True) < to_decimal(value)

    def __le__(self, value):
        return self.to_seconds(decimal=True) <= to_decimal(value)

    def __gt__(self, value):
        return self.to_seconds(decimal=True) > to_decimal(value)

    def __ge__(self, value):
        return self.to_seconds(decimal=True) >= to_decimal(value)

    def __add__(self, value):
        if isinstance(value, (int, float, Decimal)):
            seconds = self.to_seconds(decimal=True)
            seconds += to_decimal(value)
            return Date(seconds)

        raise TypeError("Date can only be added to a time shift.")

    def __sub__(self, value):
        if isinstance(value, (int, float, Decimal)):
            seconds = self.to_seconds(decimal=True)
            seconds -= to_decimal(value)
            return Date(seconds)

        if isinstance(value, Date):
            seconds = self.to_seconds(decimal=True)
            seconds -= value.to_seconds(decimal=True)
            return float(seconds)

        raise TypeError("Unsupported subtraction.")

    @classmethod
    def now(cls):
        """
        Return the current UTC time.

        Returns
        -------
        Date
            Current UTC date.
        """
        return cls(utc_now())

    @classmethod
    def from_iso8601(cls, date_string):
        """
        Create a Date from an ISO 8601 string.

        Parameters
        ----------
        date_string : str
            ISO 8601 date string.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls(date_string)

    @classmethod
    def from_dict(cls, date_dict):
        """
        Create a Date from a dictionary.

        Parameters
        ----------
        date_dict : dict
            Dictionary with calendar fields or ordinal day.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls(date_dict)

    @classmethod
    def from_ordinal(cls, year, ordinal_day, hour=0, minute=0,
                     second=0.0):
        """
        Create a Date from year and ordinal day.

        Parameters
        ----------
        year : int
            Year.
        ordinal_day : int
            Day of year.
        hour : int, optional
            Hour.
        minute : int, optional
            Minute.
        second : float, optional
            Second.

        Returns
        -------
        Date
            Parsed UTC date.
        """
        return cls([year, ordinal_day, hour, minute, second])

    def copy(self):
        """
        Return an independent copy of the date.

        Returns
        -------
        Date
            Copied date.
        """
        return Date(self.to_seconds(decimal=True), error=self.error)

    @property
    def iso8601(self):
        """
        Return the date as an ISO 8601 UTC string.
        """
        return self.to_iso8601()

    @property
    def ordinal_day(self):
        """
        Return the day of year.

        Notes
        -----
        This is the ordinal day within the year, not the astronomical
        Julian Date.
        """
        return to_ordinal_day(self.year, self.month, self.day)

    @ordinal_day.setter
    def ordinal_day(self, value):
        month, day = from_ordinal_day(self.year, value)
        self.month = month
        self.day = day
        self._selfcheck()

    @property
    def day_of_year(self):
        """
        Return the day of year.

        This is an alias of ``ordinal_day``.
        """
        return self.ordinal_day

    @day_of_year.setter
    def day_of_year(self, value):
        self.ordinal_day = value

    @property
    def year_day(self):
        """
        Return the day of year.

        This is an alias of ``ordinal_day``.
        """
        return self.ordinal_day

    @year_day.setter
    def year_day(self, value):
        self.ordinal_day = value

    @property
    def seconds(self):
        """
        Return seconds from the internal epoch.
        """
        return self.to_seconds()

    @seconds.setter
    def seconds(self, value):
        self.from_seconds(value)

    def set_date(self, date):
        """
        Set the date from one of the supported input formats.
        """
        if isinstance(date, Date):
            self._set_from_sequence(date.get_date())

        elif isinstance(date, (list, tuple)):
            self._set_from_sequence(date)

        elif isinstance(date, dict):
            self._set_from_dict(date)

        elif isinstance(date, str):
            values = read_iso8601_date(date)
            self._set_from_sequence(values[:6])

        elif isinstance(date, (int, float, Decimal)):
            self.from_seconds(date)
            return

        else:
            raise TypeError("Unsupported date format.")

        self._selfcheck()

    def get_date(self, dtype="list", format="calendar"):
        """
        Return the date in the requested format.

        Parameters
        ----------
        dtype : str, optional
            Output data type. Supported values are ``"list"``,
            ``"dict"``, ``"string"``, and ``"iso8601"``.
        format : str, optional
            Date format. Supported values are ``"calendar"`` and
            ``"ordinal"``.

        Returns
        -------
        list, dict or str
            Date representation.
        """
        if format not in ("calendar", "ordinal"):
            raise ValueError("Unsupported date format.")

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
            if format == "calendar":
                return self.to_iso8601()

            return self.to_iso8601(format="ordinal")

        raise ValueError("Unsupported output type.")

    def to_dict(self, format="calendar"):
        """
        Return the date as a dictionary.

        Parameters
        ----------
        format : str, optional
            Either ``"calendar"`` or ``"ordinal"``.

        Returns
        -------
        dict
            Date dictionary.
        """
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

    def to_ordinal(self):
        """
        Return the date as an ordinal-date list.

        Returns
        -------
        list
            Date as [year, ordinal_day, hour, minute, second].
        """
        return self.get_date(dtype="list", format="ordinal")

    def to_iso8601(self, precision=6, format="calendar"):
        """
        Return an ISO 8601 UTC string.

        Parameters
        ----------
        precision : int, optional
            Number of decimal digits for seconds.
        format : str, optional
            Either ``"calendar"`` or ``"ordinal"``.

        Returns
        -------
        str
            ISO 8601 UTC string.
        """
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

    def strformat(self, fmt, precision=4):
        """
        Format the date using a simple custom formatter.

        Supported placeholders are ``%Y``, ``%j``, ``%m``, ``%d``,
        ``%H``, ``%M`` and ``%S``.
        """
        second = round(float(self.second), precision)
        second_int = int(second)
        second_frac = second - second_int
        frac_value = round(second_frac * 10 ** precision)

        if frac_value == 10 ** precision:
            second_int += 1
            frac_value = 0

        second_string = f"{second_int:02d}.{frac_value:0{precision}d}"

        replacements = {
            "%Y": f"{self.year:04d}",
            "%j": f"{self.ordinal_day:03d}",
            "%m": f"{self.month:02d}",
            "%d": f"{self.day:02d}",
            "%H": f"{self.hour:02d}",
            "%M": f"{self.minute:02d}",
            "%S": second_string,
        }

        for key, value in replacements.items():
            fmt = fmt.replace(key, value)

        return fmt

    def strparse(self, date_string, fmt):
        """
        Parse a date string using a restricted custom format.

        This method supports fixed-width fields only. The supported
        placeholders are ``%Y``, ``%j``, ``%m``, ``%d``, ``%H``,
        ``%M`` and ``%S``.
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

    def to_seconds(self, decimal=False):
        """
        Return seconds from the internal epoch.
        """
        return date_to_second(
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            decimal=decimal,
        )

    def from_seconds(self, second):
        """
        Set the date from seconds from the internal epoch.
        """
        values = second_to_date(second)
        self._set_from_sequence(values)
        self._selfcheck()

    def shift_time(self, time, units="s"):
        """
        Shift the date by a time interval.

        Parameters
        ----------
        time : int, float or Decimal
            Time shift.
        units : str, optional
            Time units. Supported values are seconds, minutes, hours,
            and days.
        """
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

    def _set_from_sequence(self, values):
        """
        Set date fields from a calendar or ordinal sequence.
        """
        if len(values) == 6:
            year = int(values[0])
            month = int(values[1])
            day = int(values[2])
            hour = int(values[3])
            minute = int(values[4])
            second = float(values[5])
            self._set_components(year, month, day, hour, minute, second)
            return

        if len(values) == 5:
            year = int(values[0])
            ordinal = int(values[1])
            hour = int(values[2])
            minute = int(values[3])
            second = float(values[4])
            month, day = from_ordinal_day(year, ordinal)
            self._set_components(year, month, day, hour, minute, second)
            return

        raise ValueError("Sequence must have 5 or 6 elements.")

    def _set_from_dict(self, date):
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

    def _set_components(self, year, month, day, hour, minute, second):
        """
        Set raw date components.
        """
        self.year = int(year)
        self.month = int(month)
        self.day = int(day)
        self.hour = int(hour)
        self.minute = int(minute)
        self.second = float(second)

    def _selfcheck(self):
        """
        Validate date fields.
        """
        if self.year < 1:
            raise ValueError("Year must be greater than or equal to 1.")

        if self.month < 1 or self.month > 12:
            raise ValueError("Month must be between 1 and 12.")

        month_days = days_in_month(self.year)
        previous = 0 if self.month == 1 else month_days[self.month - 2]
        max_day = month_days[self.month - 1] - previous

        if self.day < 1 or self.day > max_day:
            raise ValueError("Invalid day for the selected month.")

        if self.hour < 0 or self.hour > 23:
            raise ValueError("Hour must be between 0 and 23.")

        if self.minute < 0 or self.minute > 59:
            raise ValueError("Minute must be between 0 and 59.")

        if self.second < 0 or self.second >= 60:
            raise ValueError("Second must be in the range [0, 60).")


def leap_check(year):
    """
    Return True if year is a leap year.
    """
    year = int(year)

    if year < 1:
        raise ValueError("Year must be greater than or equal to 1.")

    divisible_by_4 = year % 4 == 0
    not_century = year % 100 != 0
    divisible_by_400 = year % 400 == 0

    return divisible_by_4 and (not_century or divisible_by_400)


def leap_num(year):
    """
    Return the number of leap years before the given year.
    """
    year = int(year)

    if year < 1:
        raise ValueError("Year must be greater than or equal to 1.")

    previous = year - 1
    return previous // 4 - previous // 100 + previous // 400


def days_in_month(year):
    """
    Return cumulative days at the end of each month.
    """
    if leap_check(year):
        return [31, 60, 91, 121, 152, 182,
                213, 244, 274, 305, 335, 366]

    return [31, 59, 90, 120, 151, 181,
            212, 243, 273, 304, 334, 365]


def year_days(year):
    """
    Return the number of days in a year.
    """
    if leap_check(year):
        return 366

    return 365


def date_to_second(year=1, month=1, day=1, hour=0, minute=0,
                   second=0.0, decimal=False):
    """
    Convert a UTC date to seconds from the internal epoch.
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


def second_to_date(second):
    """
    Convert seconds from the internal epoch to a UTC date.
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


def from_ordinal_day(year, day):
    """
    Convert ordinal day of year to month and day.
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


def to_ordinal_day(year, month, day):
    """
    Convert month and day to ordinal day of year.
    """
    year = int(year)
    month = int(month)
    day = int(day)

    _validate_components(year, month, day, 0, 0, 0.0)

    if month == 1:
        return day

    return days_in_month(year)[month - 2] + day


def read_iso8601_date(date_string):
    """
    Parse an ISO 8601 date string and return UTC components.

    Supported examples are:

    - ``2026-05-09T10:30:12Z``
    - ``2026-05-09T10:30:12.345Z``
    - ``2026-129T10:30:12Z``
    - ``2026-05-09T12:30:12+02:00``

    Returns
    -------
    tuple
        ``(year, month, day, hour, minute, second, offset)``. The first
        six fields are converted to UTC. The offset is expressed in
        seconds and corresponds to the input offset from UTC.
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


def write_iso8601_date(year, month, day, hour, minute, second,
                       precision=6, timezone="Z", format="calendar"):
    """
    Convert date components to an ISO 8601 UTC string.
    """
    _validate_components(year, month, day, hour, minute, second)

    second_value = round(float(second), precision)

    if second_value >= 60:
        date = Date([year, month, day, hour, minute, 0.0])
        date.shift_time(1, "minute")
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        second_value = 0.0

    if format == "calendar":
        date_part = f"{year:04d}-{month:02d}-{day:02d}"
    elif format == "ordinal":
        ordinal = to_ordinal_day(year, month, day)
        date_part = f"{year:04d}-{ordinal:03d}"
    else:
        raise ValueError("Unsupported ISO 8601 date format.")

    width = 3 + precision
    time_part = (
        f"{hour:02d}:{minute:02d}:"
        f"{second_value:0{width}.{precision}f}"
    )

    return f"{date_part}T{time_part}{timezone}"


def utc_now():
    """
    Return the current UTC time as an ISO 8601 string.
    """
    current = systime.gmtime(systime.time())

    return write_iso8601_date(
        current.tm_year,
        current.tm_mon,
        current.tm_mday,
        current.tm_hour,
        current.tm_min,
        float(current.tm_sec),
    )


def to_decimal(value):
    """
    Convert a numeric time value to Decimal.
    """
    if isinstance(value, Date):
        return value.to_seconds(decimal=True)

    if isinstance(value, Decimal):
        return value

    if isinstance(value, (int, float)):
        return Decimal(str(value))

    raise TypeError("Value cannot be converted to Decimal.")


def _validate_components(year, month, day, hour, minute, second):
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


def _days_before_year(year):
    """
    Return the number of days before January 1 of year.
    """
    return (year - 1) * YEAR_DAYS + leap_num(year)


def _year_from_day_count(day_count):
    """
    Return year from the number of elapsed days since the epoch.
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


def _parse_timezone_offset(timezone):
    """
    Parse an ISO 8601 timezone offset.

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


def _format_to_regex(fmt):
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


def _decimal_to_float(value):
    """
    Convert a decimal-like string to float.
    """
    return float(Decimal(str(value)))