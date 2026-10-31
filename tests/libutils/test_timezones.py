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
Unit tests for :mod:`shakelab.libutils.timezones`.

The tests check the consistency of the time zone table, the offset
formatting, standard and daylight-saving offsets for table codes and
IANA names, the fallback when the time zone database is missing, the
integration with :class:`shakelab.libutils.timeN.Date` and the search
helper. Tests that need the system time zone database are skipped when
it is not available.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import pytest

from shakelab.libutils import timezones
from shakelab.libutils.timeN import Date
from shakelab.libutils.timezones import (
    TIMEZONE,
    find_timezone,
    format_offset,
    offset_string,
    utc_offset,
)


def _has_tzdata() -> bool:
    try:
        ZoneInfo("Europe/Rome")
    except ZoneInfoNotFoundError:
        return False
    return True


needs_tzdata = pytest.mark.skipif(
    not _has_tzdata(), reason="time zone database not available"
)

WINTER = Date("2026-01-15T12:00:00Z")
SUMMER = Date("2026-07-15T12:00:00Z")


# ----------------------------------------------------------------------------
# Table
# ----------------------------------------------------------------------------


def test_table_entries_are_consistent():
    """Every entry has a label matching its offset and an IANA name."""
    pattern = re.compile(r"^\(UTC(?:[+-]\d{2}:\d{2})?\) \S")

    for code, entry in TIMEZONE.items():
        assert set(entry) == {"label", "delay", "iana"}, code
        assert isinstance(entry["delay"], float), code
        assert pattern.match(entry["label"]), code

        if entry["delay"] == 0.0:
            assert entry["label"].startswith("(UTC) "), code
        else:
            assert format_offset(entry["delay"]) in entry["label"], code


@needs_tzdata
def test_table_iana_names_exist():
    """All IANA names are known to the time zone database."""
    for code, entry in TIMEZONE.items():
        ZoneInfo(entry["iana"])


@needs_tzdata
@pytest.mark.parametrize(
    "code, delay",
    [
        ("WET", 1.0), ("GMT", 0.0), ("IST", 5.5), ("NPT", 5.75),
        ("NST", -3.5), ("VST", -4.0), ("SMST", 13.0), ("JST", 3.0),
        ("TST", 9.0), ("ET", -5.0),
    ],
)
def test_standard_offsets_match_database(code, delay):
    """Standard offsets agree with the database in local winter."""
    assert TIMEZONE[code]["delay"] == delay
    assert utc_offset(code, WINTER) == delay


# ----------------------------------------------------------------------------
# Offset formatting
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "hours, text",
    [
        (0.0, "+00:00"), (1.0, "+01:00"), (5.75, "+05:45"),
        (-3.5, "-03:30"), (-12, "-12:00"), (13.0, "+13:00"),
        (0.01, "+00:01"),
    ],
)
def test_format_offset(hours, text):
    """Offsets are written as signed hours and minutes."""
    assert format_offset(hours) == text


def test_format_offset_out_of_range():
    """Offsets of a day or more are rejected."""
    with pytest.raises(ValueError):
        format_offset(24.0)


# ----------------------------------------------------------------------------
# Offsets
# ----------------------------------------------------------------------------


def test_table_offset_without_date_is_standard():
    """Without a date a table zone returns its standard offset."""
    assert utc_offset("WET") == 1.0
    assert offset_string("NST") == "-03:30"


@needs_tzdata
def test_daylight_saving_time():
    """With a date the offset includes daylight saving time."""
    assert utc_offset("WET", WINTER) == 1.0
    assert utc_offset("WET", SUMMER) == 2.0
    assert utc_offset("AEST", WINTER) == 11.0
    assert utc_offset("AEST", SUMMER) == 10.0
    assert utc_offset("IST", SUMMER) == 5.5


@needs_tzdata
def test_iana_names_and_datetimes():
    """IANA names and datetime inputs are accepted."""
    assert utc_offset("Europe/Rome", SUMMER) == 2.0
    assert utc_offset("Europe/Rome", datetime(2026, 1, 15, 12)) == 1.0

    aware = datetime(2026, 7, 15, 12, tzinfo=timezone.utc)
    assert utc_offset("America/New_York", aware) == -4.0
    assert isinstance(utc_offset("Europe/Rome"), float)


@needs_tzdata
def test_dst_transition_instant():
    """The offset changes exactly at the transition instant (UTC)."""
    before = Date("2026-03-29T00:59:59Z")
    after = Date("2026-03-29T01:00:00Z")
    assert utc_offset("Europe/Rome", before) == 1.0
    assert utc_offset("Europe/Rome", after) == 2.0


@needs_tzdata
def test_integration_with_date():
    """offset_string feeds Date.to_iso8601 directly."""
    offset = offset_string("Europe/Rome", SUMMER)
    assert offset == "+02:00"
    assert SUMMER.to_iso8601(precision=0, timezone=offset) == (
        "2026-07-15T14:00:00+02:00"
    )


def test_unknown_zone_and_bad_types():
    """Unknown zones raise ValueError, bad types TypeError."""
    with pytest.raises(ValueError):
        utc_offset("Mars/Olympus_Mons", SUMMER)
    with pytest.raises(ValueError):
        utc_offset("../etc", SUMMER)
    with pytest.raises(TypeError):
        utc_offset(1.0)
    with pytest.raises(TypeError):
        utc_offset("WET", "2026-07-15")


def test_fallback_without_database(monkeypatch):
    """Without the database, table zones use the standard offset."""
    def missing(name):
        raise ZoneInfoNotFoundError(name)

    monkeypatch.setattr(timezones, "ZoneInfo", missing)

    with pytest.warns(RuntimeWarning):
        assert utc_offset("WET", SUMMER) == 1.0

    with pytest.raises(ValueError):
        utc_offset("Europe/Rome", SUMMER)


# ----------------------------------------------------------------------------
# Search
# ----------------------------------------------------------------------------


def test_find_timezone():
    """Zones are found by city, code or IANA name, ignoring case."""
    assert find_timezone("rome") == ["WET"]
    assert find_timezone("Europe/Berlin") == ["WET"]
    assert "NZST" in find_timezone("auckland")
    assert find_timezone("atlantis") == []
    assert find_timezone("  ") == []
