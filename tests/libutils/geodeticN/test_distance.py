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
Unit tests for :mod:`shakelab.libutils.geodeticN.distance`.

These tests validate the seismological distance helpers implemented
in ``distance.py``: epicentral and hypocentral distances, and the
two tunnel-distance variants.

Run with, e.g.:

    python -m unittest test_distance
or:

    python test_distance.py
"""

from __future__ import annotations

import unittest
from math import isclose

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.distance import (
    epicentral_distance,
    hypocentral_distance,
    tunnel_distance_sphere,
    tunnel_distance_ellipsoid,
)
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
)


class TestGeodeticDistance(unittest.TestCase):
    """
    Validation tests for the distance helpers.

    The focus here is on consistency between the different distance
    measures and on the correct use of WgsPoint attributes
    (longitude, latitude, elevation).
    """

    def test_epicentral_distance_matches_great_circle(self) -> None:
        """
        Epicentral distance must match the great-circle distance for
        WgsPoint inputs.

        This codifies the convention that distance functions operate
        on objects exposing ``longitude`` and ``latitude``.
        """
        hypo = WgsPoint(longitude=12.0, latitude=45.0)
        sta = WgsPoint(longitude=13.0, latitude=45.0)

        d_epi = epicentral_distance(hypo, sta)
        d_gc = great_circle_distance(
            hypo.longitude,
            hypo.latitude,
            sta.longitude,
            sta.latitude,
            radius_m=MEAN_EARTH_RADIUS_M,
        )

        self.assertTrue(
            isclose(d_epi, d_gc, rel_tol=1e-12),
            msg="epicentral_distance mismatch with great_circle_distance",
        )

    def test_hypocentral_distance_greater_than_epicentral(self) -> None:
        """
        Hypocentral distance must be larger than epicentral distance
        when the hypocenter is below the surface.
        """
        hypo = WgsPoint(longitude=12.0, latitude=45.0, elevation=-10000.0)
        sta = WgsPoint(longitude=13.0, latitude=45.0, elevation=0.0)

        d_epi = epicentral_distance(hypo, sta)
        d_hyp = hypocentral_distance(hypo, sta, approx="ellipsoid")

        self.assertGreater(
            d_hyp,
            d_epi,
            msg="hypocentral_distance should exceed epicentral_distance",
        )

    def test_tunnel_distance_sphere_and_ellipsoid_consistency(self) -> None:
        """
        Tunnel distances computed on the sphere and ellipsoid should
        be close for shallow depths and moderate offsets.
        """
        p1 = WgsPoint(longitude=12.0, latitude=45.0, elevation=-5000.0)
        p2 = WgsPoint(longitude=12.5, latitude=45.2, elevation=0.0)

        d_sph = tunnel_distance_sphere(p1, p2)
        d_ell = tunnel_distance_ellipsoid(p1, p2)

        # Distances must be positive and reasonably close
        self.assertGreater(d_sph, 0.0, msg="sphere tunnel distance <= 0")
        self.assertGreater(d_ell, 0.0, msg="ellipsoid tunnel distance <= 0")

        rel_diff = abs(d_sph - d_ell) / d_ell
        self.assertLess(
            rel_diff,
            0.01,
            msg="sphere vs ellipsoid tunnel distance differ > 1%",
        )

    def test_hypocentral_invalid_approx_raises(self) -> None:
        """hypocentral_distance must raise on invalid approx option."""
        hypo = WgsPoint(longitude=12.0, latitude=45.0, elevation=-5000.0)
        sta = WgsPoint(longitude=13.0, latitude=45.0, elevation=0.0)

        with self.assertRaises(ValueError):
            hypocentral_distance(hypo, sta, approx="invalid-model")


if __name__ == "__main__":
    unittest.main()
