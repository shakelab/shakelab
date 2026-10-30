"""Unit tests for :mod:`shakelab.libutils.timeN`.

The tests cover calendar helpers, conversions between calendar
components and seconds from the internal epoch, ISO 8601 parsing and
formatting, custom formats, comparisons, arithmetic, hashing and the
handling of unset dates.
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
    assert leap_check(year) is expected


@pytest.mark.parametrize("year, expected", [(1, 0), (5, 1), (2001, 485)])
def test_leap_num_counts_previous_leap_years(year, expected):
    assert leap_num(year) == expected


@pytest.mark.parametrize("function", [leap_check, leap_num])
def test_calendar_helpers_reject_year_zero(function):
    with pytest.raises(ValueError, match="Year"):
        function(0)


def test_days_in_month_is_cumulative():
    assert days_in_month(2026)[-1] == year_days(2026) == 365
    assert days_in_month(2024)[1] == 60
    assert year_days(2024) == 366


@pytest.mark.parametrize("year", [2024, 2026])
def test_ordinal_day_round_trip(year):
    for ordinal in range(1, year_days(year) + 1):
        month, day = from_ordinal_day(year, ordinal)
        assert to_ordinal_day(year, month, day) == ordinal


@pytest.mark.parametrize("year, ordinal, expected", [
    (2024, 60, (2, 29)),
    (2026, 60, (3, 1)),
    (2026, 365, (12, 31)),
])
def test_from_ordinal_day(year, ordinal, expected):
    assert from_ordinal_day(year, ordinal) == expected


@pytest.mark.parametrize("ordinal", [0, 366])
def test_from_ordinal_day_rejects_out_of_range(ordinal):
    with pytest.raises(ValueError, match="Ordinal day"):
        from_ordinal_day(2026, ordinal)


# ---------------------------------------------------------------------------
# Seconds from the internal epoch
# ---------------------------------------------------------------------------

def test_epoch_origin():
    assert date_to_second(1, 1, 1) == 0.0
    assert second_to_date(0) == [1, 1, 1, 0, 0, 0.0]


def test_unix_epoch_offset():
    seconds = date_to_second(1970, 1, 1, decimal=True)
    assert seconds == Decimal(UNIX_EPOCH_SECONDS)


@pytest.mark.parametrize("components", [
    [1, 1, 1, 0, 0, 0.0],
    [2000, 12, 31, 23, 59, 59.5],
    [2024, 2, 29, 12, 0, 0.25],
    [2026, 10, 30, 15, 25, 7.125],
])
def test_seconds_round_trip(components):
    seconds = date_to_second(*components, decimal=True)
    assert second_to_date(seconds) == components


def test_second_to_date_rejects_negative_values():
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
    assert Date(value).get_date() == expected


def test_read_iso8601_date_returns_utc_components_and_offset():
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
    assert Date(value) == Date(REFERENCE)


@pytest.mark.parametrize("value", [
    REFERENCE_SECONDS,
    float(REFERENCE_SECONDS),
    int(REFERENCE_SECONDS),
])
def test_date_from_seconds(value):
    assert Date(value) == Date(REFERENCE)


def test_dictionary_without_day_fields_is_rejected():
    with pytest.raises(ValueError, match="calendar or ordinal"):
        Date({"year": 2026, "hour": 1})


@pytest.mark.parametrize("value", [object(), 1j, True])
def test_unsupported_input_types_are_rejected(value):
    with pytest.raises(TypeError):
        Date(value)


def test_only_utc_timezone_is_supported():
    with pytest.raises(NotImplementedError):
        Date(REFERENCE, timezone="CET")


def test_year_zero_is_rejected():
    with pytest.raises(ValueError, match="Year"):
        Date([0, 1, 1, 0, 0, 0.0])


def test_date_constructor_copies_another_date():
    original = Date(REFERENCE)
    copied = Date(original)
    original.shift_time(1, "h")
    assert copied == Date(REFERENCE)


def test_copy_preserves_value_and_error():
    original = Date("2026-10-30T15:25:07.5Z", error=0.01)
    copied = original.copy()
    assert copied == original
    assert copied is not original
    assert copied.error == 0.01


# ---------------------------------------------------------------------------
# Current time
# ---------------------------------------------------------------------------

def test_now_keeps_sub_second_resolution(fixed_clock):
    now = Date.now()
    assert now.get_date()[:5] == [2023, 11, 14, 22, 13]
    assert now.second == pytest.approx(20.123456789)


def test_now_string_matches_now_classmethod(fixed_clock):
    assert Date("now") == Date.now()


def test_utc_now_reports_microseconds(fixed_clock):
    assert utc_now() == "2023-11-14T22:13:20.123457Z"


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
    date = Date("2026-10-30T15:25:07.125Z")
    assert date.to_iso8601(precision) == expected


def test_iso8601_property_repr_and_str():
    date = Date(REFERENCE)
    assert date.iso8601 == "2026-10-30T15:25:07.000000Z"
    assert repr(date) == date.iso8601
    assert str(date) == date.iso8601


def test_to_iso8601_ordinal_format():
    date = Date(REFERENCE)
    assert date.to_iso8601(0, format="ordinal") == "2026-303T15:25:07Z"


@pytest.mark.parametrize("precision, expected", [
    (6, "2027-01-01T00:00:00.000000Z"),
    (0, "2027-01-01T00:00:00Z"),
])
def test_iso8601_rounding_carries_into_next_minute(precision, expected):
    date = Date([2026, 12, 31, 23, 59, 59.9999996])
    assert date.to_iso8601(precision) == expected


def test_write_iso8601_date_rejects_negative_precision():
    with pytest.raises(ValueError, match="precision"):
        write_iso8601_date(2026, 10, 30, 15, 25, 7.0, precision=-1)


@pytest.mark.parametrize("fmt, precision, expected", [
    ("%Y-%m-%dT%H:%M:%S", 4, "2026-10-30T15:25:07.1250"),
    ("%Y%m%d_%H%M%S", 0, "20261030_152507"),
    ("%Y.%j", 4, "2026.303"),
])
def test_strformat(fmt, precision, expected):
    date = Date("2026-10-30T15:25:07.125Z")
    assert date.strformat(fmt, precision=precision) == expected


def test_strformat_rounding_carries_into_next_minute():
    date = Date([2026, 12, 31, 23, 59, 59.99996])
    expected = "2027-01-01T00:00:00.0000"
    assert date.strformat("%Y-%m-%dT%H:%M:%S") == expected


@pytest.mark.parametrize("text, fmt", [
    ("20261030_152507", "%Y%m%d_%H%M%S"),
    ("2026-303 15:25:07", "%Y-%j %H:%M:%S"),
    ("2026/10/30 15:25:07.0", "%Y/%m/%d %H:%M:%S"),
])
def test_strparse(text, fmt):
    date = Date()
    date.strparse(text, fmt)
    assert date == Date(REFERENCE)


def test_strparse_rejects_mismatching_strings():
    with pytest.raises(ValueError, match="does not match"):
        Date().strparse("2026-10-30", "%Y%m%d")


def test_get_date_representations():
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
    with pytest.raises(ValueError):
        Date(REFERENCE).get_date(dtype, fmt)


def test_ordinal_day_setter_changes_month_and_day():
    date = Date("2024-01-01T00:00:00Z")
    date.ordinal_day = 60
    assert (date.month, date.day) == (2, 29)
    assert date.day_of_year == date.year_day == 60


# ---------------------------------------------------------------------------
# Comparisons and arithmetic
# ---------------------------------------------------------------------------

def test_comparisons():
    early = Date(REFERENCE)
    late = Date("2026-10-30T15:25:08Z")
    assert early < late
    assert early <= late
    assert late > early
    assert late >= early
    assert early != late
    assert early == Date("2026-10-30T17:25:07+02:00")


def test_equality_with_seconds():
    date = Date(REFERENCE)
    assert date == date.to_seconds(decimal=True)


def test_equality_with_unsupported_types_is_false():
    assert Date(REFERENCE) != REFERENCE


def test_ordering_with_unsupported_types_raises():
    with pytest.raises(TypeError):
        Date(REFERENCE) < "2026-10-30T15:25:08Z"


@pytest.mark.parametrize("shift, expected", [
    (1.5, "2026-10-30T15:25:08.500000Z"),
    (Decimal("0.000001"), "2026-10-30T15:25:07.000001Z"),
    (-86400, "2026-10-29T15:25:07.000000Z"),
    (86400 * 62, "2026-12-31T15:25:07.000000Z"),
])
def test_addition_of_time_shift(shift, expected):
    assert (Date(REFERENCE) + shift).iso8601 == expected


def test_addition_crosses_leap_day():
    date = Date("2024-02-28T23:59:59Z") + 1
    assert date.get_date() == [2024, 2, 29, 0, 0, 0.0]


def test_subtraction():
    late = Date("2026-10-30T15:25:07.25Z")
    early = Date("2026-10-30T15:24:07Z")
    assert late - early == pytest.approx(60.25)
    assert late - 60.25 == early


def test_addition_does_not_modify_operand():
    date = Date(REFERENCE)
    shifted = date + 10
    assert shifted != date
    assert date == Date(REFERENCE)


@pytest.mark.parametrize("operand", ["1", None, Date(REFERENCE)])
def test_addition_rejects_unsupported_operands(operand):
    with pytest.raises(TypeError):
        Date(REFERENCE) + operand


@pytest.mark.parametrize("units, expected", [
    ("s", "2026-10-30T15:25:37Z"),
    ("minutes", "2026-10-30T15:55:07Z"),
    ("h", "2026-10-31T21:25:07Z"),
    ("days", "2026-11-29T15:25:07Z"),
])
def test_shift_time_units(units, expected):
    date = Date(REFERENCE)
    date.shift_time(30, units)
    assert date == Date(expected)


def test_shift_time_rejects_unknown_units():
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
    assert Date(value) == Date(REFERENCE)


def test_boolean_time_shift_is_rejected():
    with pytest.raises(TypeError):
        Date(REFERENCE) + True


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_non_finite_time_shift_is_rejected(value):
    with pytest.raises(ValueError, match="finite"):
        Date(REFERENCE) + value


def test_to_decimal_conversions():
    assert to_decimal(np.int64(3)) == Decimal(3)
    assert to_decimal(0.1) == Decimal("0.1")
    assert to_decimal(Decimal("0.25")) == Decimal("0.25")
    assert to_decimal(Date(REFERENCE)) == REFERENCE_SECONDS


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

def test_equal_dates_have_equal_hashes():
    utc = Date(REFERENCE)
    offset = Date("2026-10-30T17:25:07+02:00")
    assert hash(utc) == hash(offset)
    assert len({utc, offset}) == 1


def test_dates_can_be_dictionary_keys():
    index = {Date(REFERENCE): "event"}
    assert index[Date("2026-303T15:25:07Z")] == "event"


# ---------------------------------------------------------------------------
# Unset dates
# ---------------------------------------------------------------------------

def test_unset_date_has_explicit_representation():
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
    with pytest.raises(ValueError, match="not set"):
        operation(Date())


def test_unset_date_copy_is_unset():
    copied = Date(error=0.5).copy()
    assert not copied.is_set
    assert copied.error == 0.5
