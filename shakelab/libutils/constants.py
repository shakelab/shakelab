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
Physical and mathematical constants used across ShakeLab.

Conventions
-----------
- Values are in SI units unless stated otherwise.
- Geodetic constants (Earth radii, WGS84 ellipsoid parameters) live in
  :mod:`shakelab.libutils.geodeticN.utilities` and are not duplicated
  here.
"""

from __future__ import annotations

import math

__all__ = [
    "GRAVITY",
    "GRAVITY_LOW",
    "GRAVITY_HIGH",
    "PI",
]


# ----------------------------------------------------------------------------
# Gravity
# ----------------------------------------------------------------------------

GRAVITY: float = 9.80665
"""Standard acceleration of gravity, g0 (m/s^2), exact by definition."""

GRAVITY_LOW: float = 9.764
"""Approximate lowest gravity at the Earth's surface (m/s^2)."""

GRAVITY_HIGH: float = 9.834
"""Approximate highest gravity at the Earth's surface (m/s^2)."""


# ----------------------------------------------------------------------------
# Mathematical constants
# ----------------------------------------------------------------------------

PI: float = math.pi
"""Ratio of a circle's circumference to its diameter."""
