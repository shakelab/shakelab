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
Tests of :mod:`shakelab.signals.codes`.
"""

from __future__ import annotations

import pytest

from shakelab.signals.codes import (
    SidPattern,
    expand_braces,
    match_code,
    match_sid,
    parse_sid_pattern,
)


@pytest.mark.parametrize("pattern, expected", [
    ("OX.CAE..HHZ", ["OX.CAE..HHZ"]),
    ("OX.{CAE,BAD}", ["OX.CAE", "OX.BAD"]),
    ("HH{Z,N}", ["HHZ", "HHN"]),
    ("{A,B}.{1,2}", ["A.1", "A.2", "B.1", "B.2"]),
    ("X{a,{b,c}d}", ["Xa", "Xbd", "Xcd"]),
    ("{A,A,B}", ["A", "B"]),
    ("{,00}", ["", "00"]),
])
def test_expand_braces(pattern, expected):
    assert expand_braces(pattern) == expected


@pytest.mark.parametrize("pattern", ["{A,B", "A}B", "A}{B", "{{A}"])
def test_unbalanced_braces(pattern):
    with pytest.raises(ValueError):
        expand_braces(pattern)


def test_parse_sid_pattern():
    assert parse_sid_pattern("IV.MODE") == [
        SidPattern("IV", "MODE", "*", "*")
    ]
    assert parse_sid_pattern(" IV.MODE..HN? ") == [
        SidPattern("IV", "MODE", "", "HN?")
    ]
    assert parse_sid_pattern("OX.{CAE,BAD}.*.HHZ") == [
        SidPattern("OX", "CAE", "*", "HHZ"),
        SidPattern("OX", "BAD", "*", "HHZ"),
    ]
    assert parse_sid_pattern("IV", parts=2) == [("IV", "*")]
    assert str(parse_sid_pattern("IV")[0]) == "IV.*.*.*"


@pytest.mark.parametrize("pattern", [
    "", "  ", "A.B.C.D.E", "IV..*.HHZ", "IV.MO DE", "IV.MODE.00.",
    "IV.{MODE,}",
])
def test_parse_invalid(pattern):
    with pytest.raises(ValueError):
        parse_sid_pattern(pattern)


def test_parse_types():
    with pytest.raises(TypeError):
        parse_sid_pattern(None)
    with pytest.raises(ValueError):
        parse_sid_pattern("IV", parts=3)


@pytest.mark.parametrize("sid, pattern, expected", [
    ("OX.CAE..HHZ", "OX.{CAE,BAD}", True),
    ("OX.BAD..HHN", "OX.{CAE,BAD}.*.HH[ZN]", True),
    ("OX.BAD..HHE", "OX.{CAE,BAD}.*.HH[ZN]", False),
    ("OX.CAE..HHZ", "*Z", False),
    ("OX.CAE..HHZ", "*.*.*.*Z", True),
    ("OX.CAE..HHZ", "OX.CAE..HHZ", True),
    ("OX.CAE.00.HHZ", "OX.CAE..HHZ", False),
    ("OX.CAE..HHZ", "OX.C*", True),
    ("OX.CAE..HHZ", "O*.*", True),
    (None, "*", False),
    ("OX.CAE", "*", False),
])
def test_match_sid(sid, pattern, expected):
    assert match_sid(sid, pattern) is expected


def test_match_code():
    assert match_code("IV.MODE", "IV")
    assert match_code("IV.MODE", "{OX,IV}.MO*")
    assert not match_code("IV.MODE", "OX")
    assert not match_code(None, "*")
    assert not match_code("IV.MODE..HNZ", "*")


def test_metadata_select_code_by_code():
    from shakelab.signals.metadata import ChannelEpoch, MetadataCollection

    collection = MetadataCollection()
    for sid in ["IV.MODE..HNZ", "IV.MODE..HNE", "OX.CAE..HHZ"]:
        collection.append(ChannelEpoch(sid=sid))
    assert [e.sid for e in collection.select("IV")] == [
        "IV.MODE..HNZ", "IV.MODE..HNE",
    ]
    assert [e.sid for e in collection.select("{IV,OX}.*.*.H?Z")] == [
        "IV.MODE..HNZ", "OX.CAE..HHZ",
    ]
    assert len(collection.select("*Z")) == 0
