# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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
Module for basic waveform analysis - TEST
"""
import numpy as np
import pytest

from shakelab.signals.base import Header, Record
from shakelab.libutils.time import Date
from shakelab.libutils.geodetic import WgsPoint


def test_header_initialization():
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
    assert isinstance(header.time, Date)
    assert isinstance(header.location, WgsPoint)


def test_header_rate_sync():
    header = Header()

    header.rate = 50.0

    assert header.rate == 50.0
    assert header.delta == 0.02


def test_header_invalid_sampling_values():
    header = Header()

    for value in (0, -1, np.nan, np.inf):
        try:
            header.delta = value
        except ValueError:
            pass
        else:
            raise AssertionError(
                "Invalid delta accepted: {0}".format(value)
            )

def test_header_copy_does_not_copy_parent():
    parent = object()
    header = Header(
        delta=0.01,
        sid="IV.TEST..HHZ",
        parent=parent,
    )
    header.meta["quality"] = "D"

    copied = header.copy()

    assert copied is not header
    assert copied._parent is None
    assert copied.meta == header.meta
    assert copied.meta is not header.meta


