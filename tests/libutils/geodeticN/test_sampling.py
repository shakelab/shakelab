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
Unit tests for :mod:`shakelab.libutils.geodeticN.sampling`.

These tests validate the spatial sampling utilities implemented in
``sampling.py``:

* :func:`sample_uniform_in_polygon_xy`
* :func:`sample_uniform_in_polygon_wgs84`

Run with, e.g.:

    python -m unittest test_sampling
or:

    python test_sampling.py
"""

from __future__ import annotations

import unittest
from math import isclose

import numpy as np

from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
    _ray_casting_contains,
)
from shakelab.libutils.geodeticN.sampling import (
    sample_uniform_in_polygon_xy,
    sample_uniform_in_polygon_wgs84,
)
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
)


class TestSamplingXY(unittest.TestCase):
    """
    Tests for uniform sampling in planar polygons.
    """

    def test_sample_xy_invalid_args(self) -> None:
        """n_samples and vertices_xy validation."""
        vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_xy(vertices, 0)

        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_xy(vertices, -10)

        # Wrong shape
        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_xy(
                np.array([1.0, 2.0]), 10
            )

        # Too few vertices
        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_xy(
                np.array([[0.0, 0.0], [1.0, 0.0]]),
                10,
            )

    def test_sample_xy_triangle_inside(self) -> None:
        """
        All returned points must lie inside the input triangle.

        The check uses the same ray-casting predicate used internally
        by the sampling functions.
        """
        np.random.seed(0)

        vertices = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ],
            dtype=float,
        )
        n_samples = 200

        pts = sample_uniform_in_polygon_xy(vertices, n_samples)
        self.assertEqual(
            pts.shape,
            (n_samples, 2),
            msg="Unexpected shape of sampled points",
        )

        x = vertices[:, 0]
        y = vertices[:, 1]

        for px, py in pts:
            inside = _ray_casting_contains(
                x,
                y,
                float(px),
                float(py),
            )
            self.assertTrue(
                inside,
                msg="Sampled point lies outside planar triangle",
            )

        # Rough sanity check on centroid location: for a right triangle
        # with vertices (0,0), (1,0), (0,1), centroid is at (1/3, 1/3).
        mean_x = float(np.mean(pts[:, 0]))
        mean_y = float(np.mean(pts[:, 1]))

        self.assertTrue(
            isclose(mean_x, 1.0 / 3.0, abs_tol=0.1)
            and isclose(mean_y, 1.0 / 3.0, abs_tol=0.1),
            msg="Sampled centroid far from expected triangle centroid",
        )


class TestSamplingWGS84(unittest.TestCase):
    """
    Tests for uniform sampling in WGS84 polygons.
    """

    def test_sample_wgs84_invalid_args(self) -> None:
        """n_samples validation and empty polygon handling."""
        square = WgsPolygon.from_lonlat(
            [
                (-1.0, -1.0),
                (1.0, -1.0),
                (1.0, 1.0),
                (-1.0, 1.0),
            ]
        )

        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_wgs84(square, 0)

        with self.assertRaises(ValueError):
            sample_uniform_in_polygon_wgs84(square, -5)

        # Empty polygon returns an empty mesh.
        empty_poly = WgsPolygon(parts=[])
        mesh = sample_uniform_in_polygon_wgs84(empty_poly, 10)
        self.assertIsInstance(mesh, WgsMesh)
        self.assertEqual(len(mesh), 0)

    def test_sample_wgs84_basic_square(self) -> None:
        """
        Basic sampling inside a simple square near the equator.

        Checks:
        * the number of samples matches the request;
        * points remain within a loose bounding box around the polygon;
        * the sample centroid is close to the polygon centroid.
        """
        np.random.seed(0)

        coords = [
            (-1.0, -1.0),
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
        ]
        polygon = WgsPolygon.from_lonlat(coords)
        n_samples = 300

        mesh = sample_uniform_in_polygon_wgs84(polygon, n_samples)

        self.assertIsInstance(mesh, WgsMesh)
        self.assertEqual(
            len(mesh),
            n_samples,
            msg="Sampled mesh size mismatch",
        )

        # Bounding box of the original polygon.
        lon_min, lon_max, lat_min, lat_max = polygon.bounding_box()

        # All sampled points should remain reasonably close to the
        # original bbox (allowing for small numerical effects).
        margin = 0.2  # degrees

        for p in mesh:
            self.assertGreaterEqual(
                p.longitude,
                lon_min - margin,
                msg="Sampled longitude smaller than polygon bbox",
            )
            self.assertLessEqual(
                p.longitude,
                lon_max + margin,
                msg="Sampled longitude greater than polygon bbox",
            )
            self.assertGreaterEqual(
                p.latitude,
                lat_min - margin,
                msg="Sampled latitude smaller than polygon bbox",
            )
            self.assertLessEqual(
                p.latitude,
                lat_max + margin,
                msg="Sampled latitude greater than polygon bbox",
            )

        # Centroid sanity check. For the square [-1,1]x[-1,1], the
        # centroid is at (0,0).
        lons = np.array([p.longitude for p in mesh])
        lats = np.array([p.latitude for p in mesh])

        mean_lon = float(lons.mean())
        mean_lat = float(lats.mean())

        self.assertTrue(
            abs(mean_lon) < 0.2 and abs(mean_lat) < 0.2,
            msg="Sampled centroid far from expected polygon centroid",
        )

    def test_sample_wgs84_multipolygon(self) -> None:
        """
        Sampling on a simple multipolygon: ensure non-empty result and
        that sampled points roughly cover both parts.
        """
        np.random.seed(1)

        # Two small disjoint squares.
        poly1 = [
            (-2.0, -1.0),
            (-1.0, -1.0),
            (-1.0, 0.0),
            (-2.0, 0.0),
        ]
        poly2 = [
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
        ]

        parts = [poly1, poly2]
        polygon = WgsPolygon.from_parts(parts)
        n_samples = 400

        mesh = sample_uniform_in_polygon_wgs84(polygon, n_samples)
        self.assertIsInstance(mesh, WgsMesh)
        self.assertEqual(len(mesh), n_samples)

        # Count how many points fall "near" each part by simple bbox
        # tests. This is just a sanity check that both parts receive
        # a non-negligible number of samples.
        n1 = 0
        n2 = 0
        for p in mesh:
            if -2.1 <= p.longitude <= -0.9 and -1.1 <= p.latitude <= 0.1:
                n1 += 1
            if 0.9 <= p.longitude <= 2.1 and -0.1 <= p.latitude <= 1.1:
                n2 += 1

        self.assertGreater(
            n1,
            0,
            msg="No samples found near first polygon part",
        )
        self.assertGreater(
            n2,
            0,
            msg="No samples found near second polygon part",
        )

        # The two sets should not be extremely unbalanced given similar
        # areas. Allow for randomness with a generous factor.
        ratio = n1 / max(1, n2)
        self.assertTrue(
            0.25 < ratio < 4.0,
            msg="Samples are too unevenly distributed across parts",
        )


if __name__ == "__main__":
    unittest.main()
