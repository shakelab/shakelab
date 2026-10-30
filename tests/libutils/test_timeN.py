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
Unit tests for :mod:`shakelab.libutils.timeN`.

The tests cover calendar helpers, conversions between calendar
components and seconds from the internal epoch, ISO 8601 parsing and
formatting, custom formats, comparisons, arithmetic, hashing and the
handling of unset dates.

Conventions
-----------
- ``REFERENCE`` (2026-10-30T15:25:07Z, ordinal day 303) is the common
  reference instant.
- Tests on the current time freeze the system clock with the
  ``fixed_clock`` fixture.
"""

from __future__ import annotations

from decimal import Decimal

import numpy as np
import pytest

from shakelab.libutils import timeN
from shakelab.libutils.timeN import (
    Date,
    date_to_second,
    days_in_month,
    from_ordinal_day,
    leap_check,
    leap_num,
    read_iso8601_date,
    second_to_date,
    to_decimal,
    to_ordinal_day,
    utc_now,
    write_iso8601_date,
    year_days,
)


UNIX_EPOCH_SECONDS = 62135596800
FIXED_TIME_NS = 1_700_000_000_123_456_789
REFERENCE = "2026-10-30T15:25:07Z"
REFERENCE_SECONDS = date_to_second(2026, 10, 30, 15, 25, 7, decimal=True)


@pytest.fixture
def fixed_clock(monkeypatch):
    """Freeze the system clock at 2023-11-14T22:13:20.123456789Z."""
    monkeypatch.setattr(timeN.systime, "time_ns", lambda: FIXED_TIME_NS)


# ---------------------------------------------------------------------------
# Calendar helpers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("year, expected", [
    (1, False),
    (4, True),
    (1900, False),
    (2000, True),
    (2024, True),
    (2026, False),
])
def test_leap_check(year, expected):
    """Gregorian leap-year rules are applied."""
    assert leap_check(year) is expected


@pytest.mark.parametrize("year, expected", [(1, 0), (5, 1), (2001, 485)])
def test_leap_num_counts_previous_leap_years(year, expected):
    """Leap years strictly before the given year are counted."""
    assert leap_num(year) == expected


@pytest.mark.parametrize("function", [leap_check, leap_num])
def test_calendar_helpers_reject_year_zero(function):
    """Year zero is outside the supported calendar."""
    with pytest.raises(ValueError, match="Year"):
        function(0)


def test_days_in_month_is_cumulative():
    """Month boundaries are cumulative day counts."""
    assert days_in_month(2026)[-1] == year_days(2026) == 365
    assert days_in_month(2024)[1] == 60
    assert year_days(2024) == 366


@pytest.mark.parametrize("year", [2024, 2026])
def test_ordinal_day_round_trip(year):
    """Every ordinal day survives a round trip through (month, day)."""
    for ordinal in range(1, year_days(year) + 1):
        month, day = from_ordinal_day(year, ordinal)
        assert to_ordinal_day(year, month, day) == ordinal


@pytest.mark.parametrize("year, ordinal, expected", [
    (2024, 60, (2, 29)),
    (2026, 60, (3, 1)),
    (2026, 365, (12, 31)),
])
def test_from_ordinal_day(year, ordinal, expected):
    """Ordinal days are converted to the expected (month, day)."""
    assert from_ordinal_day(year, ordinal) == expected


@pytest.mark.parametrize("ordinal", [0, 366])
def test_from_ordinal_day_rejects_out_of_range(ordinal):
    """Ordinal days outside the year are rejected."""
    with pytest.raises(ValueError, match="Ordinal day"):
        from_ordinal_day(2026, ordinal)


# ---------------------------------------------------------------------------
# Seconds from the internal epoch
# ---------------------------------------------------------------------------

def test_epoch_origin():
    """The internal epoch is 0001-01-01T00:00:00Z."""
    assert date_to_second(1, 1, 1) == 0.0
    assert second_to_date(0) == [1, 1, 1, 0, 0, 0.0]


def test_unix_epoch_offset():
    """The Unix epoch lies 62135596800 s after the internal epoch."""
    seconds = date_to_second(1970, 1, 1, decimal=True)
    assert seconds == Decimal(UNIX_EPOCH_SECONDS)


@pytest.mark.parametrize("components", [
    [1, 1, 1, 0, 0, 0.0],
    [2000, 12, 31, 23, 59, 59.5],
    [2024, 2, 29, 12, 0, 0.25],
    [2026, 10, 30, 15, 25, 7.125],
])
def test_seconds_round_trip(components):
    """Calendar components survive a round trip through seconds."""
    seconds = date_to_second(*components, decimal=True)
    assert second_to_date(seconds) == components


def test_second_to_date_rejects_negative_values():
    """Instants before the internal epoch are rejected."""
    with pytest.raises(ValueError, match="greater than or equal to zero"):
        second_to_date(-1)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value, expected", [
    ("2026-10-30T15:25:07Z", [2026, 10, 30, 15, 25, 7.0]),
    ("2026-10-30T15:25:07.125Z", [2026, 10, 30, 15, 25, 7.125]),
    ("2026-10-30T15:25:07", [2026, 10, 30, 15, 25, 7.0]),
    ("2026-303T15:25:07Z", [2026, 10, 30, 15, 25, 7.0]),
    ("2026303T15:25:07Z", [2026, 10, 30, 15, 25, 7.0]),
    ("2026-10-30T17:25:07+02:00", [2026, 10, 30, 15, 25, 7.0]),
    ("2026-10-30T10:55:07-0430", [2026, 10, 30, 15, 25, 7.0]),
    ("2026-01-01T01:00:00+02:00", [2025, 12, 31, 23, 0, 0.0]),
])
def test_date_from_iso8601(value, expected):
    """ISO 8601 strings are parsed and converted to UTC."""
    assert Date(value).get_date() == expected


def test_read_iso8601_date_returns_utc_components_and_offset():
    """The parser returns UTC components and the input offset."""
    assert read_iso8601_date("2026-10-30T17:25:07+02:00") == (
        2026, 10, 30, 15, 25, 7.0, 7200.0,
    )


@pytest.mark.parametrize("value", [
    "2026-02-29T00:00:00Z",
    "2026-13-01T00:00:00Z",
    "2026-366T00:00:00Z",
    "2026-01-01T24:00:00Z",
    "2026-01-01T00:60:00Z",
    "2026-01-01T00:00:60Z",
    "2026-01-01 00:00:00Z",
    "2026-01-01",
    "not a date",
])
def test_invalid_iso8601_strings_are_rejected(value):
    """Malformed or out-of-range ISO 8601 strings are rejected."""
    with pytest.raises(ValueError):
        Date(value)


@pytest.mark.parametrize("value", [
    [2026, 10, 30, 15, 25, 7.0],
    (2026, 303, 15, 25, 7.0),
    {
        "year": 2026, "month": 10, "day": 30,
        "hour": 15, "minute": 25, "second": 7.0,
    },
    {"year": 2026, "ordinal_day": 303, "hour": 15, "minute": 25,
     "second": 7.0},
    {"year": 2026, "day_of_year": 303, "hour": 15, "minute": 25,
     "second": 7.0},
    {"year": 2026, "year_day": 303, "hour": 15, "minute": 25,
     "second": 7.0},
])
def test_date_from_components(value):
    """Calendar and ordinal sequences and dictionaries are accepted."""
    assert Date(value) == Date(REFERENCE)


@pytest.mark.parametrize("value", [
    REFERENCE_SECONDS,
    float(REFERENCE_SECONDS),
    int(REFERENCE_SECONDS),
])
def test_date_from_seconds(value):
    """Seconds from the internal epoch are accepted."""
    assert Date(value) == Date(REFERENCE)


def test_dictionary_without_day_fields_is_rejected():
    """A dictionary needs calendar or ordinal day fields."""
    with pytest.raises(ValueError, match="calendar or ordinal"):
        Date({"year": 2026, "hour": 1})


@pytest.mark.parametrize("value", [object(), 1j, True])
def test_unsupported_input_types_are_rejected(value):
    """Unsupported input types raise TypeError."""
    with pytest.raises(TypeError):
        Date(value)


def test_only_utc_timezone_is_supported():
    """Time zones other than UTC are not implemented."""
    with pytest.raises(NotImplementedError):
        Date(REFERENCE, timezone="CET")


def test_year_zero_is_rejected():
    """Year zero is rejected on construction."""
    with pytest.raises(ValueError, match="Year"):
        Date([0, 1, 1, 0, 0, 0.0])


def test_date_constructor_copies_another_date():
    """Building a Date from another Date creates an independent copy."""
    original = Date(REFERENCE)
    copied = Date(original)
    original.shift_time(1, "h")
    assert copied == Date(REFERENCE)


def test_copy_preserves_value_and_error():
    """copy() returns an independent date with the same error."""
    original = Date("2026-10-30T15:25:07.5Z", error=0.01)
    copied = original.copy()
    assert copied == original
    assert copied is not original
    assert copied.error == 0.01


# ---------------------------------------------------------------------------
# Current time
# ---------------------------------------------------------------------------

def test_now_keeps_sub_second_resolution(fixed_clock):
    """Date.now() keeps the sub-second part of the system clock."""
    now = Date.now()
    assert now.get_date()[:5] == [2023, 11, 14, 22, 13]
    assert now.second == pytest.approx(20.123456789)


def test_now_string_matches_now_classmethod(fixed_clock):
    """Date("now") and Date.now() agree."""
    assert Date("now") == Date.now()


def test_utc_now_reports_microseconds(fixed_clock):
    """utc_now() reports microseconds."""
    assert utc_now() == "2023-11-14T22:13:20.123457Z"


# ---------------------------------------------------------------------------
# Unix time
# ---------------------------------------------------------------------------

def test_unix_epoch_constant():
    """The exported Unix epoch matches the calendar computation."""
    assert timeN.UNIX_EPOCH_SECONDS == Decimal(UNIX_EPOCH_SECONDS)
    assert Date.from_unix(0) == Date("1970-01-01T00:00:00Z")


@pytest.mark.parametrize("seconds, expected", [
    (1793373907, REFERENCE),
    (1700000000.5, "2023-11-14T22:13:20.5Z"),
    (Decimal("1700000000.000000001"), "2023-11-14T22:13:20.000000001Z"),
    (np.int64(1700000000), "2023-11-14T22:13:20Z"),
    (-1, "1969-12-31T23:59:59Z"),
])
def test_from_unix(seconds, expected):
    """Unix timestamps are converted to UTC dates."""
    assert Date.from_unix(seconds) == Date(expected)


def test_to_unix():
    """to_unix() returns float or exact Decimal Unix seconds."""
    date = Date(REFERENCE)
    assert date.to_unix() == 1793373907.0
    assert date.to_unix(decimal=True) == Decimal(1793373907)


def test_to_unix_keeps_microseconds():
    """Float Unix timestamps resolve microseconds at present dates."""
    date = Date("2026-10-30T15:25:07.000001Z")
    assert date.to_unix() == pytest.approx(1793373907.000001, abs=5e-7)


def test_unix_ns_round_trip():
    """Integer nanoseconds survive a round trip through Date."""
    assert Date.from_unix_ns(FIXED_TIME_NS).to_unix_ns() == FIXED_TIME_NS


def test_to_unix_ns_rounds_to_nearest_nanosecond():
    """Sub-nanosecond parts are rounded to the nearest nanosecond."""
    assert Date.from_unix(Decimal("1.0000000006")).to_unix_ns() == (
        1_000_000_001
    )
    assert Date.from_unix(Decimal("1.0000000004")).to_unix_ns() == (
        1_000_000_000
    )


@pytest.mark.parametrize("value", [1.5, True, "1"])
def test_from_unix_ns_rejects_non_integers(value):
    """Only integers are accepted as nanoseconds."""
    with pytest.raises(TypeError, match="integer"):
        Date.from_unix_ns(value)


def test_now_is_exact_in_nanoseconds(fixed_clock):
    """Date.now() reproduces the system clock to the nanosecond."""
    assert Date.now().to_unix_ns() == FIXED_TIME_NS


# ---------------------------------------------------------------------------
# Formatting and parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("precision, expected", [
    (6, "2026-10-30T15:25:07.125000Z"),
    (3, "2026-10-30T15:25:07.125Z"),
    (1, "2026-10-30T15:25:07.1Z"),
    (0, "2026-10-30T15:25:07Z"),
])
def test_to_iso8601_precision(precision, expected):
    """The number of decimal digits follows the precision."""
    date = Date("2026-10-30T15:25:07.125Z")
    assert date.to_iso8601(precision) == expected


def test_iso8601_property_repr_and_str():
    """iso8601, repr() and str() share the ISO 8601 representation."""
    date = Date(REFERENCE)
    assert date.iso8601 == "2026-10-30T15:25:07.000000Z"
    assert repr(date) == date.iso8601
    assert str(date) == date.iso8601


def test_to_iso8601_ordinal_format():
    """The ordinal ISO 8601 format uses the day of year."""
    date = Date(REFERENCE)
    assert date.to_iso8601(0, format="ordinal") == "2026-303T15:25:07Z"


@pytest.mark.parametrize("precision, expected", [
    (6, "2027-01-01T00:00:00.000000Z"),
    (0, "2027-01-01T00:00:00Z"),
])
def test_iso8601_rounding_carries_into_next_minute(precision, expected):
    """Rounding to 60 s carries into the next minute and year."""
    date = Date([2026, 12, 31, 23, 59, 59.9999996])
    assert date.to_iso8601(precision) == expected


def test_write_iso8601_date_rejects_negative_precision():
    """A negative precision is rejected."""
    with pytest.raises(ValueError, match="precision"):
        write_iso8601_date(2026, 10, 30, 15, 25, 7.0, precision=-1)


@pytest.mark.parametrize("fmt, precision, expected", [
    ("%Y-%m-%dT%H:%M:%S", 4, "2026-10-30T15:25:07.1250"),
    ("%Y%m%d_%H%M%S", 0, "20261030_152507"),
    ("%Y.%j", 4, "2026.303"),
])
def test_strformat(fmt, precision, expected):
    """Custom formats are rendered field by field."""
    date = Date("2026-10-30T15:25:07.125Z")
    assert date.strformat(fmt, precision=precision) == expected


def test_strformat_rounding_carries_into_next_minute():
    """strformat() carries rounding to 60 s into the next minute."""
    date = Date([2026, 12, 31, 23, 59, 59.99996])
    expected = "2027-01-01T00:00:00.0000"
    assert date.strformat("%Y-%m-%dT%H:%M:%S") == expected


@pytest.mark.parametrize("text, fmt", [
    ("20261030_152507", "%Y%m%d_%H%M%S"),
    ("2026-303 15:25:07", "%Y-%j %H:%M:%S"),
    ("2026/10/30 15:25:07.0", "%Y/%m/%d %H:%M:%S"),
])
def test_strparse(text, fmt):
    """Fixed-width custom formats are parsed."""
    date = Date()
    date.strparse(text, fmt)
    assert date == Date(REFERENCE)


def test_strparse_rejects_mismatching_strings():
    """A string that does not match the format is rejected."""
    with pytest.raises(ValueError, match="does not match"):
        Date().strparse("2026-10-30", "%Y%m%d")


def test_get_date_representations():
    """get_date(), to_dict() and to_ordinal() agree."""
    date = Date(REFERENCE)
    assert date.get_date("list", "ordinal") == [2026, 303, 15, 25, 7.0]
    assert date.to_dict()["month"] == 10
    assert date.to_dict("ordinal")["ordinal_day"] == 303
    assert date.to_ordinal() == [2026, 303, 15, 25, 7.0]
    assert date.get_date("iso8601") == date.iso8601


@pytest.mark.parametrize("dtype, fmt", [
    ("set", "calendar"),
    ("list", "julian"),
])
def test_get_date_rejects_unknown_options(dtype, fmt):
    """Unknown output types and date formats are rejected."""
    with pytest.raises(ValueError):
        Date(REFERENCE).get_date(dtype, fmt)


def test_ordinal_day_setter_changes_month_and_day():
    """Setting the ordinal day updates month and day."""
    date = Date("2024-01-01T00:00:00Z")
    date.ordinal_day = 60
    assert (date.month, date.day) == (2, 29)
    assert date.day_of_year == date.year_day == 60


# ---------------------------------------------------------------------------
# Comparisons and arithmetic
# ---------------------------------------------------------------------------

def test_comparisons():
    """Dates are ordered by instant, independently of the offset."""
    early = Date(REFERENCE)
    late = Date("2026-10-30T15:25:08Z")
    assert early < late
    assert early <= late
    assert late > early
    assert late >= early
    assert early != late
    assert early == Date("2026-10-30T17:25:07+02:00")


def test_equality_with_seconds():
    """A date equals its seconds from the internal epoch."""
    date = Date(REFERENCE)
    assert date == date.to_seconds(decimal=True)


def test_equality_with_unsupported_types_is_false():
    """Equality with unsupported types is False, not an error."""
    assert Date(REFERENCE) != REFERENCE


def test_ordering_with_unsupported_types_raises():
    """Ordering against unsupported types raises TypeError."""
    with pytest.raises(TypeError):
        Date(REFERENCE) < "2026-10-30T15:25:08Z"


@pytest.mark.parametrize("shift, expected", [
    (1.5, "2026-10-30T15:25:08.500000Z"),
    (Decimal("0.000001"), "2026-10-30T15:25:07.000001Z"),
    (-86400, "2026-10-29T15:25:07.000000Z"),
    (86400 * 62, "2026-12-31T15:25:07.000000Z"),
])
def test_addition_of_time_shift(shift, expected):
    """Adding a time shift moves the instant."""
    assert (Date(REFERENCE) + shift).iso8601 == expected


def test_addition_crosses_leap_day():
    """Arithmetic crosses the leap day correctly."""
    date = Date("2024-02-28T23:59:59Z") + 1
    assert date.get_date() == [2024, 2, 29, 0, 0, 0.0]


def test_subtraction():
    """Date minus Date is a float; Date minus a shift is a Date."""
    late = Date("2026-10-30T15:25:07.25Z")
    early = Date("2026-10-30T15:24:07Z")
    assert late - early == pytest.approx(60.25)
    assert late - 60.25 == early


def test_addition_does_not_modify_operand():
    """Arithmetic returns a new date."""
    date = Date(REFERENCE)
    shifted = date + 10
    assert shifted != date
    assert date == Date(REFERENCE)


@pytest.mark.parametrize("operand", ["1", None, Date(REFERENCE)])
def test_addition_rejects_unsupported_operands(operand):
    """Adding unsupported operands raises TypeError."""
    with pytest.raises(TypeError):
        Date(REFERENCE) + operand


@pytest.mark.parametrize("units, expected", [
    ("s", "2026-10-30T15:25:37Z"),
    ("minutes", "2026-10-30T15:55:07Z"),
    ("h", "2026-10-31T21:25:07Z"),
    ("days", "2026-11-29T15:25:07Z"),
])
def test_shift_time_units(units, expected):
    """shift_time() supports seconds, minutes, hours and days."""
    date = Date(REFERENCE)
    date.shift_time(30, units)
    assert date == Date(expected)


def test_shift_time_rejects_unknown_units():
    """Unknown time units are rejected."""
    with pytest.raises(ValueError, match="time unit"):
        Date(REFERENCE).shift_time(1, "weeks")


# ---------------------------------------------------------------------------
# NumPy scalars and invalid time values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("shift", [
    np.int64(5),
    np.int32(5),
    np.float64(5.0),
    np.float32(5.0),
])
def test_numpy_scalars_are_time_shifts(shift):
    """NumPy scalars work as time shifts on both sides."""
    date = Date(REFERENCE)
    expected = Date("2026-10-30T15:25:12Z")
    assert date + shift == expected
    assert shift + date == expected
    assert expected - shift == date


@pytest.mark.parametrize("value", [
    np.int64(int(REFERENCE_SECONDS)),
    np.float64(float(REFERENCE_SECONDS)),
])
def test_numpy_scalars_are_accepted_as_seconds(value):
    """NumPy scalars are accepted as seconds from the epoch."""
    assert Date(value) == Date(REFERENCE)


def test_boolean_time_shift_is_rejected():
    """Booleans are not time shifts."""
    with pytest.raises(TypeError):
        Date(REFERENCE) + True


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_non_finite_time_shift_is_rejected(value):
    """Non-finite time shifts are rejected."""
    with pytest.raises(ValueError, match="finite"):
        Date(REFERENCE) + value


def test_to_decimal_conversions():
    """to_decimal() converts dates and real numbers exactly."""
    assert to_decimal(np.int64(3)) == Decimal(3)
    assert to_decimal(0.1) == Decimal("0.1")
    assert to_decimal(Decimal("0.25")) == Decimal("0.25")
    assert to_decimal(Date(REFERENCE)) == REFERENCE_SECONDS


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

def test_equal_dates_have_equal_hashes():
    """Equal dates have equal hashes."""
    utc = Date(REFERENCE)
    offset = Date("2026-10-30T17:25:07+02:00")
    assert hash(utc) == hash(offset)
    assert len({utc, offset}) == 1


def test_dates_can_be_dictionary_keys():
    """Dates can be used as dictionary keys."""
    index = {Date(REFERENCE): "event"}
    assert index[Date("2026-303T15:25:07Z")] == "event"


# ---------------------------------------------------------------------------
# Unset dates
# ---------------------------------------------------------------------------

def test_unset_date_has_explicit_representation():
    """An unset date is reported as Date(None)."""
    date = Date()
    assert not date.is_set
    assert repr(date) == "Date(None)"


@pytest.mark.parametrize("operation", [
    lambda date: date.to_seconds(),
    lambda date: date.to_iso8601(),
    lambda date: date.strformat("%Y"),
    lambda date: date.shift_time(1),
    lambda date: date + 1,
    lambda date: date < Date(REFERENCE),
    lambda date: hash(date),
    lambda date: Date(date),
])
def test_unset_date_operations_raise_clear_error(operation):
    """Operations on an unset date raise a clear ValueError."""
    with pytest.raises(ValueError, match="not set"):
        operation(Date())


def test_unset_date_copy_is_unset():
    """Copying an unset date returns an unset date."""
    copied = Date(error=0.5).copy()
    assert not copied.is_set
    assert copied.error == 0.5
