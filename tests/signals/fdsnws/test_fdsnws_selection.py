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
Tests of :mod:`shakelab.signals.fdsnws.selection`.

Patterns of :mod:`shakelab.signals.codes` translated into the codes,
GET parameters and POST lines of the FDSN services.
"""

from __future__ import annotations

import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.codes import parse_sid_pattern
from shakelab.signals.fdsnws import selection
from shakelab.signals.fdsnws.selection import (
    get_params,
    parse_selection,
    post_body,
    server_codes,
)

T0 = "2012-05-20T02:00:00"
T1 = "2012-05-20T03:00:00"


def _codes(pattern: str) -> tuple:
    return server_codes(parse_sid_pattern(pattern)[0])


# ----------------------------------------------------------------------------
# Selections
# ----------------------------------------------------------------------------


def test_parse_single_pattern():
    """A pattern with the window of the arguments."""
    (request,) = parse_selection("IV.MODE", T0, T1)
    assert str(request.pattern) == "IV.MODE.*.*"
    assert request.starttime == Date(T0)
    assert request.endtime == Date(T1)


def test_parse_braces_and_lists():
    """Brace alternatives become requests, without repetitions."""
    requests = parse_selection(["IV.{MODE,CAVE}..HNZ", "IV.MODE..HNZ"])
    assert [str(r.pattern) for r in requests] == [
        "IV.MODE..HNZ", "IV.CAVE..HNZ"
    ]
    assert requests[0].starttime is None


def test_parse_windows_of_items():
    """Items with their own windows; None takes the arguments."""
    requests = parse_selection(
        [("IV.MODE", "2012-05-20T02:30:00", None), "OX.*"], T0, T1
    )
    assert requests[0].starttime == Date("2012-05-20T02:30:00")
    assert requests[0].endtime == Date(T1)
    assert requests[1].starttime == Date(T0)


def test_tuple_of_patterns():
    """A top-level tuple is a sequence of patterns."""
    requests = parse_selection(("IV.A", "IV.B", "IV.C"))
    assert len(requests) == 3


@pytest.mark.parametrize(
    "select, error",
    [
        (None, TypeError),
        ([], ValueError),
        ([1], TypeError),
        ([("IV.MODE", T0)], TypeError),
        ("IV..MODE", ValueError),
        ([("IV.MODE", T1, T0)], ValueError),
    ],
)
def test_parse_errors(select, error):
    """Invalid selections are rejected."""
    with pytest.raises(error):
        parse_selection(select)


def test_time_errors():
    """Times are Dates or ISO strings."""
    with pytest.raises(TypeError):
        parse_selection("IV", 0.0)
    with pytest.raises(ValueError):
        parse_selection("IV", "yesterday")


def test_matches():
    """The original pattern decides, not the server codes."""
    (request,) = parse_selection("IV.MODE..HN[ZN]")
    assert request.matches("IV.MODE..HNZ")
    assert not request.matches("IV.MODE..HNE")
    assert not request.matches("IV.MODE.00.HNZ")


# ----------------------------------------------------------------------------
# Server codes
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "pattern, expected",
    [
        ("IV.MODE..HNZ", ("IV", "MODE", "--", "HNZ")),
        ("IV.MODE", ("IV", "MODE", "*", "*")),
        ("IV.MODE.*.HH[ZNE]", ("IV", "MODE", "*", "HHZ,HHN,HHE")),
        ("IV.MODE.0[0-2].HHZ", ("IV", "MODE", "00,01,02", "HHZ")),
        ("IV.MODE.*.HH[!Z]", ("IV", "MODE", "*", "HH?")),
        ("IV.M?DE.*.?H*", ("IV", "M?DE", "*", "?H*")),
    ],
)
def test_server_codes(pattern, expected):
    """Empty locations, sets, ranges and negated sets."""
    assert _codes(pattern) == expected


def test_large_sets_become_jokers():
    """Too many combinations: the sets become ``?``."""
    codes = _codes("IV.[A-Z][A-Z]")
    assert codes[1] == "??"


# ----------------------------------------------------------------------------
# GET and POST
# ----------------------------------------------------------------------------


def test_get_params():
    """One request: GET parameters; more: None."""
    requests = parse_selection("IV.MODE..HN[ZN]", T0, T1)
    params = get_params(requests)
    assert params == {
        "network": "IV",
        "station": "MODE",
        "location": "--",
        "channel": "HNZ,HNN",
        "starttime": Date(T0),
        "endtime": Date(T1),
    }
    assert get_params(parse_selection(["IV.A", "IV.B"])) is None


def test_post_body():
    """Parameter lines, then one line per combination of codes."""
    requests = parse_selection(
        ["IV.MODE..HN[ZN]", ("OX.*", None, None)], T0, T1
    )
    body = post_body(requests, {"quality": "B", "minimumlength": None})
    assert body.splitlines() == [
        "quality=B",
        "IV MODE -- HNZ 2012-05-20T02:00:00.000000 "
        "2012-05-20T03:00:00.000000",
        "IV MODE -- HNN 2012-05-20T02:00:00.000000 "
        "2012-05-20T03:00:00.000000",
        "OX * * * 2012-05-20T02:00:00.000000 2012-05-20T03:00:00.000000",
    ]


def test_post_body_without_times():
    """Missing times are written as ``*``."""
    body = post_body(parse_selection(["IV.A", "IV.B"]))
    assert body.splitlines() == ["IV A * * * *", "IV B * * * *"]


def test_to_date():
    """Dates are kept, strings parsed, None kept."""
    date = Date(T0)
    assert selection.to_date(date) is date
    assert selection.to_date(T0) == date
    assert selection.to_date(None) is None
