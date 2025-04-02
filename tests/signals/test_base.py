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
    header = Header()
    assert header._rate is None
    assert header._delta is None
    assert header.sid is None
    assert header.eid is None
    assert isinstance(header.time, Date)
    assert isinstance(header.location, WgsPoint)
    assert header.units is None
    assert header.response is None
    assert header.meta == {}
    assert header._parent is None

def test_delta_rate_setters():
    header = Header()

    header.delta = 0.1
    assert header._delta == 0.1
    assert header._rate == 10.0

    header.rate = 20.0
    assert header._rate == 20.0
    assert header._delta == 0.05

def test_copy_method():
    header = Header()
    header.sid = '123'
    header.eid = '456'
    header.delta = 0.1
    header.units = 'm/s'
    header.meta = {'key': 'value'}

    header_copy = header.copy()
    assert header_copy.sid == header.sid
    assert header_copy.eid == header.eid
    assert header_copy.delta == header.delta
    assert header_copy.rate == header.rate
    assert header_copy.units == header.units
    assert header_copy.meta == header.meta
    assert header_copy is not header
    assert header_copy.meta is not header.meta  # Ensures a deep copy


@pytest.fixture
def sample_data():
    return np.array([1, 2, 3, 4, 5])

@pytest.fixture
def record(sample_data):
    return Record(data=sample_data, delta=0.1, time=0)

def test_record_initialization(record, sample_data):
    assert np.array_equal(record.data, sample_data)
    assert record.delta == 0.1
    assert record.time == 0