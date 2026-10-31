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
Unit tests for :mod:`shakelab.libutils.constants`.

The tests check the values, types and consistency of the exported
constants.
"""

from __future__ import annotations

import math

import pytest

from shakelab.libutils import constants
from shakelab.libutils.constants import (
    GRAVITY,
    GRAVITY_HIGH,
    GRAVITY_LOW,
    PI,
)


# ----------------------------------------------------------------------------
# Values
# ----------------------------------------------------------------------------


def test_standard_gravity_value():
    """GRAVITY is the exact standard acceleration of gravity."""
    assert GRAVITY == 9.80665


def test_gravity_bounds_enclose_standard_value():
    """The surface gravity range must contain the standard value."""
    assert GRAVITY_LOW < GRAVITY < GRAVITY_HIGH


def test_pi_matches_math_module():
    """PI must be the double-precision value of pi."""
    assert PI == math.pi


# ----------------------------------------------------------------------------
# Public interface
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("name", constants.__all__)
def test_exported_constants_are_floats(name):
    """Every exported constant is a plain float."""
    assert type(getattr(constants, name)) is float
