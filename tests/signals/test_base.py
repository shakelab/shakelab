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
Unit tests for :mod:`shakelab.signals.base`.

The tests are organized by class (``Header``, ``Record``, ``Stream``,
``StreamCollection``). They fix the behavior of the public API used in
the ShakeLab demos and check the corrected defects. Real data come from
``tests/signals/data`` (station IV.MODE, Emilia earthquake of 20 May
2012; see the README in that folder).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Header, Record

DATA = Path(__file__).parent / "data"


# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------


def test_header_initialization():
    """All constructor arguments are stored and converted."""
    header = Header(
        delta=0.01,
        time="2026-01-01T00:00:00",
        location=(13.0, 46.0),
        sid="IV.TEST..HHZ",
        eid="event_001",
        units="m/s",
    )

    assert header.delta == 0.01
    assert header.rate == 100.0
    assert header.sid == "IV.TEST..HHZ"
    assert header.eid == "event_001"
    assert header.units == "m/s"
    assert header.time == Date("2026-01-01T00:00:00")
    assert header.location == WgsPoint(longitude=13.0, latitude=46.0)
    assert header.response is None
    assert header.meta == {}


def test_header_defaults():
    """An empty header has no sampling, time or location."""
    header = Header()

    assert header.delta is None
    assert header.rate is None
    assert header.time is None
    assert header.location is None
    assert header.nsamp is None


def test_header_rate_and_delta_are_synchronized():
    """Setting the rate updates the interval and vice versa."""
    header = Header()

    header.rate = 50.0
    assert header.delta == 0.02

    header.delta = 0.005
    assert header.rate == 200.0

    header.delta = None
    assert header.rate is None


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
@pytest.mark.parametrize("name", ["delta", "rate"])
def test_header_invalid_sampling_values(name, value):
    """Sampling values must be finite and strictly positive."""
    with pytest.raises(ValueError):
        setattr(Header(), name, value)


@pytest.mark.parametrize("value", ["0.01", True, [0.01]])
def test_header_sampling_type_errors(value):
    """Sampling values must be real numbers, booleans excluded."""
    with pytest.raises(TypeError):
        Header(delta=value)


def test_header_accepts_numpy_sampling_values():
    """NumPy scalars are accepted and stored as float."""
    header = Header()
    header.rate = np.int64(100)

    assert header.delta == 0.01
    assert type(header.rate) is float


@pytest.mark.parametrize(
    "value",
    [
        "2012-05-20T02:03:23Z",
        Date("2012-05-20T02:03:23Z"),
        [2012, 5, 20, 2, 3, 23],
        Date("2012-05-20T02:03:23Z").to_seconds(),
    ],
)
def test_header_time_accepts_date_formats(value):
    """The start time accepts every format accepted by Date."""
    header = Header(time=value)

    assert header.time == Date("2012-05-20T02:03:23Z")


def test_header_time_is_copied():
    """Changing an assigned Date afterwards does not move the header."""
    date = Date("2012-05-20T02:03:23Z")
    header = Header(time=date)

    date.year = 2000

    assert header.time == Date("2012-05-20T02:03:23Z")
    assert header.time is not date


def test_header_invalid_time():
    """Unsupported or malformed times are rejected."""
    with pytest.raises(TypeError):
        Header(time=object())
    with pytest.raises(ValueError):
        Header(time="not a date")


@pytest.mark.parametrize(
    "value",
    [
        (13.0, 46.0, 250.0),
        [13.0, 46.0, 250.0],
        np.array([13.0, 46.0, 250.0]),
    ],
)
def test_header_location_sequences(value):
    """Locations are (longitude, latitude[, elevation]) sequences."""
    header = Header(location=value)

    assert header.location == WgsPoint(
        longitude=13.0, latitude=46.0, elevation=250.0
    )


def test_header_location_point_is_kept():
    """A WgsPoint is stored as is."""
    point = WgsPoint(longitude=13.0, latitude=46.0)

    assert Header(location=point).location is point


def test_header_invalid_location():
    """Wrong lengths, strings and invalid latitudes are rejected."""
    with pytest.raises(ValueError):
        Header(location=(13.0,))
    with pytest.raises(TypeError):
        Header(location="13.0, 46.0")
    with pytest.raises(TypeError):
        Header(location=13.0)
    with pytest.raises(ValueError):
        Header(location=(13.0, 95.0))


def test_header_nsamp_from_parent():
    """nsamp is read from the parent record."""
    record = Record(data=np.zeros(25), delta=0.01)

    assert record.head.nsamp == 25


def test_header_copy_is_independent():
    """A copy shares no mutable state and drops the parent."""
    parent = object()
    header = Header(
        delta=0.01,
        time="2026-01-01T00:00:00",
        sid="IV.TEST..HHZ",
        parent=parent,
    )
    header.meta["quality"] = "D"

    copied = header.copy()

    assert copied is not header
    assert copied._parent is None
    assert copied.delta == header.delta
    assert copied.sid == header.sid
    assert copied.time == header.time
    assert copied.time is not header.time
    assert copied.meta == header.meta
    assert copied.meta is not header.meta

    other = object()
    assert header.copy(parent=other)._parent is other


def test_header_text_representations(capsys):
    """repr is compact; str and info list all fields."""
    header = Header(delta=0.01, sid="IV.TEST..HHZ")

    assert repr(header).startswith("Header(sid='IV.TEST..HHZ'")
    assert "sampling rate: 100.0" in str(header)

    header.info()
    assert "source identifier: IV.TEST..HHZ" in capsys.readouterr().out
