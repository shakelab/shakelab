#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Unit tests for :mod:`shakelab.libutils.geodeticN.utilities`.

These tests validate the core geodetic helpers implemented in the
``utilities`` module: angle conversions, longitude normalization,
DMS parsing, spherical distances and basic ellipsoid helpers.

Run with, e.g.:

    python -m unittest test_utilities
or:

    python test_utilities.py
"""

from __future__ import annotations

import unittest
from math import isclose, pi

from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    EQUATORIAL_EARTH_RADIUS_M,
    POLAR_EARTH_RADIUS_M,
    deg_to_rad,
    rad_to_deg,
    normalize_longitude,
    dms2dec,
    great_circle_distance,
    chord_distance,
    geocentric_radius,
    forward_azimuth,
)


class TestGeodeticUtilities(unittest.TestCase):
    """
    Basic validation tests for the geodeticN utilities.

    Tests are intentionally simple and deterministic, so they can be
    used as smoke tests whenever refactoring the geodetic code.
    """

    # ------------------------------------------------------------------
    # Angle helpers
    # ------------------------------------------------------------------

    def test_angle_conversions(self) -> None:
        """Test deg_to_rad and rad_to_deg for scalars and arrays."""
        rad = float(deg_to_rad(180.0))
        self.assertTrue(
            isclose(rad, pi, rel_tol=1e-12),
            msg="deg_to_rad(180) failed",
        )

        deg = float(rad_to_deg(pi))
        self.assertTrue(
            isclose(deg, 180.0, rel_tol=1e-12),
            msg="rad_to_deg(pi) failed",
        )

        arr = deg_to_rad([0.0, 90.0])
        self.assertEqual(arr.shape, (2,), msg="deg_to_rad array shape")

        self.assertTrue(
            isclose(float(arr[1]), pi / 2.0, rel_tol=1e-12),
            msg="deg_to_rad array value failed",
        )

    def test_normalize_longitude(self) -> None:
        """Check longitude normalization into [-180, 180)."""
        self.assertTrue(
            isclose(
                float(normalize_longitude(190.0)),
                -170.0,
                rel_tol=1e-12,
            ),
            msg="normalize_longitude(190) failed",
        )

        self.assertTrue(
            isclose(
                float(normalize_longitude(-190.0)),
                170.0,
                rel_tol=1e-12,
            ),
            msg="normalize_longitude(-190) failed",
        )

        arr = normalize_longitude([0.0, 360.0, -540.0])
        expected = [0.0, 0.0, -180.0]
        for v, e in zip(arr, expected):
            self.assertTrue(
                isclose(float(v), e, rel_tol=1e-12),
                msg="normalize_longitude array failed",
            )

    def test_dms2dec(self) -> None:
        """Validate DMS string parsing to decimal degrees."""
        self.assertTrue(
            isclose(dms2dec("12"), 12.0, rel_tol=1e-12),
            msg="dms2dec simple degrees failed",
        )

        self.assertTrue(
            isclose(dms2dec("12 30"), 12.5, rel_tol=1e-12),
            msg="dms2dec degrees+minutes failed",
        )

        self.assertTrue(
            isclose(dms2dec("12 30 0"), 12.5, rel_tol=1e-12),
            msg="dms2dec degrees+minutes+seconds failed",
        )

        self.assertTrue(
            isclose(dms2dec("12°30'00\"N"), 12.5, rel_tol=1e-12),
            msg="dms2dec DMS north failed",
        )

        self.assertTrue(
            isclose(dms2dec("12°30'S"), -12.5, rel_tol=1e-12),
            msg="dms2dec DMS south failed",
        )

        # Hemisphere overrides numeric sign
        self.assertTrue(
            isclose(dms2dec("-12°30'N"), 12.5, rel_tol=1e-12),
            msg="dms2dec hemisphere overrides sign failed",
        )

    # ------------------------------------------------------------------
    # Distances on the sphere
    # ------------------------------------------------------------------

    def test_great_circle_distance_equator_90deg(self) -> None:
        """Quarter-circumference along equator must match expectation."""
        lon1, lat1 = 0.0, 0.0
        lon2, lat2 = 90.0, 0.0

        d = great_circle_distance(lon1, lat1, lon2, lat2)
        quarter_circ = 0.5 * pi * MEAN_EARTH_RADIUS_M

        self.assertTrue(
            isclose(d, quarter_circ, rel_tol=1e-6),
            msg="great_circle_distance equator 90° failed",
        )

    def test_chord_vs_arc_distance(self) -> None:
        """Chord distance must be shorter than the corresponding arc."""
        lon1, lat1 = 0.0, 0.0
        lon2, lat2 = 60.0, 30.0

        d_arc = great_circle_distance(lon1, lat1, lon2, lat2)
        d_chord = chord_distance(lon1, lat1, lon2, lat2)

        self.assertLess(d_chord, d_arc, msg="Chord must be shorter")

        ratio = d_chord / d_arc
        self.assertTrue(
            0.9 < ratio < 1.0,
            msg="Chord/arc ratio out of expected bounds",
        )

    # ------------------------------------------------------------------
    # Ellipsoid-related helpers
    # ------------------------------------------------------------------

    def test_geocentric_radius(self) -> None:
        """Check geocentric_radius at key latitudes and monotonicity."""
        r_eq = geocentric_radius(0.0)
        r_pole = geocentric_radius(90.0)
        r_mid = geocentric_radius(45.0)

        self.assertTrue(
            isclose(
                float(r_eq),
                EQUATORIAL_EARTH_RADIUS_M,
                rel_tol=1e-6,
            ),
            msg="geocentric_radius at equator failed",
        )

        self.assertTrue(
            isclose(
                float(r_pole),
                POLAR_EARTH_RADIUS_M,
                rel_tol=1e-6,
            ),
            msg="geocentric_radius at pole failed",
        )

        self.assertTrue(
            EQUATORIAL_EARTH_RADIUS_M > float(r_mid) > POLAR_EARTH_RADIUS_M,
            msg="geocentric_radius monotonicity failed",
        )

    def test_forward_azimuth(self) -> None:
        """Validate forward_azimuth for simple great-circle paths."""
        # Eastward along equator
        az = forward_azimuth(0.0, 0.0, 90.0, 0.0)
        self.assertTrue(
            isclose(az, 90.0, rel_tol=1e-12),
            msg="forward_azimuth eastward failed",
        )

        # Northward along meridian
        az = forward_azimuth(0.0, 0.0, 0.0, 10.0)
        self.assertTrue(
            isclose(az, 0.0, rel_tol=1e-12),
            msg="forward_azimuth northward failed",
        )

        # Southward along meridian
        az = forward_azimuth(0.0, 0.0, 0.0, -10.0)
        self.assertTrue(
            isclose(az, 180.0, rel_tol=1e-12),
            msg="forward_azimuth southward failed",
        )


if __name__ == "__main__":
    unittest.main()
