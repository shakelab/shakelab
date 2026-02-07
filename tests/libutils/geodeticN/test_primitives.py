#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
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
Unit tests for :mod:`shakelab.libutils.geodeticN.primitives`.

These tests validate the basic geometric primitives:

* :class:`WgsPoint`
* :class:`WgsPolygon`
* :class:`WgsMesh`

Run with, e.g.:

    python -m unittest test_primitives
or:

    python test_primitives.py
"""

from __future__ import annotations

import unittest
from math import isclose

import numpy as np

from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
)
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
)


class TestWgsPoint(unittest.TestCase):
    """
    Basic tests for the WgsPoint primitive.
    """

    def test_basic_properties(self) -> None:
        """Check attribute storage and to_tuple()."""
        p = WgsPoint(longitude=12.0, latitude=45.0, elevation=100.0)

        self.assertTrue(
            isclose(p.longitude, 12.0, rel_tol=0.0),
            msg="WgsPoint.longitude value failed",
        )
        self.assertTrue(
            isclose(p.latitude, 45.0, rel_tol=0.0),
            msg="WgsPoint.latitude value failed",
        )
        self.assertTrue(
            isclose(p.elevation, 100.0, rel_tol=0.0),
            msg="WgsPoint.elevation value failed",
        )

        lon, lat = p.to_tuple()
        self.assertTrue(
            isclose(lon, 12.0, rel_tol=0.0)
            and isclose(lat, 45.0, rel_tol=0.0),
            msg="WgsPoint.to_tuple failed",
        )

    def test_epicentral_distance_to(self) -> None:
        """
        epicentral_distance_to must match great_circle_distance.
        """
        p = WgsPoint(longitude=12.0, latitude=45.0, elevation=0.0)
        q = WgsPoint(longitude=13.0, latitude=45.0, elevation=0.0)

        d_epi = p.epicentral_distance_to(q)
        d_ref = great_circle_distance(
            12.0,
            45.0,
            13.0,
            45.0,
            radius_m=MEAN_EARTH_RADIUS_M,
        )

        self.assertTrue(
            isclose(d_epi, d_ref, rel_tol=1e-9),
            msg=(
                "epicentral_distance_to mismatch with "
                "great_circle_distance"
            ),
        )

    def test_hypocentral_distance_and_sub(self) -> None:
        """
        __sub__ must be consistent with hypocentral_distance_to().
        """
        p = WgsPoint(longitude=12.0, latitude=45.0, elevation=-5000.0)
        q = WgsPoint(longitude=12.1, latitude=45.1, elevation=0.0)

        d_hyp = p.hypocentral_distance_to(q)
        d_op = p - q

        self.assertTrue(
            isclose(d_op, d_hyp, rel_tol=1e-12),
            msg="WgsPoint.__sub__ mismatch with hypocentral_distance_to",
        )

    def test_hypocentral_distance_approx_models(self) -> None:
        """
        Sphere and ellipsoid approximations should be close for
        moderate offsets and shallow depths.
        """
        p = WgsPoint(longitude=12.0, latitude=45.0, elevation=-5000.0)
        q = WgsPoint(longitude=12.2, latitude=45.2, elevation=0.0)

        d_sphere = p.hypocentral_distance_to(q, approx="sphere")
        d_ell = p.hypocentral_distance_to(q, approx="ellipsoid")

        self.assertGreater(
            d_sphere,
            0.0,
            msg="sphere hypocentral distance should be positive",
        )
        self.assertGreater(
            d_ell,
            0.0,
            msg="ellipsoid hypocentral distance should be positive",
        )

        rel_diff = abs(d_sphere - d_ell) / d_ell
        self.assertLess(
            rel_diff,
            0.01,
            msg="sphere vs ellipsoid hypocentral distance differ > 1%",
        )


class TestWgsPolygon(unittest.TestCase):
    """
    Tests for WgsPolygon: bbox, centroid, area and point-in-polygon.
    """

    def test_basic_operations(self) -> None:
        """Basic operations on a simple square polygon."""
        coords = [
            (-1.0, -1.0),
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
        ]
        poly = WgsPolygon.from_lonlat(coords)

        self.assertFalse(poly.is_empty(), msg="Polygon should not be empty")
        self.assertEqual(
            poly.n_parts(),
            1,
            msg="Polygon should have one part",
        )

        lon_min, lon_max, lat_min, lat_max = poly.bounding_box()
        self.assertTrue(
            isclose(lon_min, -1.0, rel_tol=1e-12)
            and isclose(lon_max, 1.0, rel_tol=1e-12)
            and isclose(lat_min, -1.0, rel_tol=1e-12)
            and isclose(lat_max, 1.0, rel_tol=1e-12),
            msg="WgsPolygon.bounding_box failed",
        )

        c = poly.centroid()
        self.assertTrue(
            isclose(c.longitude, 0.0, rel_tol=1e-12)
            and isclose(c.latitude, 0.0, rel_tol=1e-12),
            msg="WgsPolygon.centroid failed",
        )

        area = poly.area()
        self.assertGreater(
            area,
            0.0,
            msg="WgsPolygon.area should be positive",
        )

        # Rough sanity check: area of a 2°x2° box near equator.
        side = MEAN_EARTH_RADIUS_M * (2.0 * np.pi / 180.0)
        approx_area = side * side
        ratio = area / approx_area
        self.assertTrue(
            0.5 < ratio < 2.0,
            msg="WgsPolygon.area out of rough expected bounds",
        )

        inside = poly.contains_point(
            WgsPoint(longitude=0.0, latitude=0.0),
        )
        outside = poly.contains_point(
            WgsPoint(longitude=10.0, latitude=10.0),
        )

        self.assertTrue(
            inside,
            msg="Point (0,0) should be inside polygon",
        )
        self.assertFalse(
            outside,
            msg="Point (10,10) should be outside polygon",
        )


class TestWgsMesh(unittest.TestCase):
    """
    Tests for WgsMesh container.
    """

    def test_basic_operations(self) -> None:
        """Check construction, length and bounding box."""
        coords = [
            (0.0, 0.0),
            (1.0, 2.0),
            (-1.0, -2.0),
        ]
        mesh = WgsMesh.from_lonlat(coords)

        self.assertEqual(
            len(mesh),
            3,
            msg="WgsMesh length failed",
        )

        lon_min, lon_max, lat_min, lat_max = mesh.bounding_box()
        self.assertTrue(
            isclose(lon_min, -1.0, rel_tol=1e-12)
            and isclose(lon_max, 1.0, rel_tol=1e-12)
            and isclose(lat_min, -2.0, rel_tol=1e-12)
            and isclose(lat_max, 2.0, rel_tol=1e-12),
            msg="WgsMesh.bounding_box failed",
        )


if __name__ == "__main__":
    unittest.main()
