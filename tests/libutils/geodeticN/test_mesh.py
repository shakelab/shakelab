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
Unit tests for :mod:`shakelab.libutils.geodeticN.mesh`.

These tests validate the mesh generation utilities implemented in
``mesh.py``:

* :func:`spherical_mesh_global`
* :func:`spherical_cap_mesh_fibonacci`
* :func:`spherical_cap_mesh_fibonacci_bbox`
* :func:`mesh_cartesian`

Run with, e.g.:

    python -m unittest test_mesh
or:

    python test_mesh.py
"""

from __future__ import annotations

import unittest
from math import isclose

from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsMesh,
)
from shakelab.libutils.geodeticN.mesh import (
    mesh_fibonacci_global,
    mesh_fibonacci_spherical_cap,
    mesh_fibonacci_bbox,
    mesh_cartesian,
    mesh_sinusoidal,
)
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
)


class TestSphericalMeshGlobal(unittest.TestCase):
    """
    Tests for the global Fibonacci mesh generator.
    """

    def test_spacing_positive_required(self) -> None:
        """spacing_km must be strictly positive."""
        with self.assertRaises(ValueError):
            mesh_fibonacci_global(0.0)

        with self.assertRaises(ValueError):
            mesh_fibonacci_global(-10.0)

    def test_global_mesh_basic_properties(self) -> None:
        """
        Check that the global mesh is non-empty, has the expected
        number of nodes and lies within valid geographic bounds.
        """
        spacing_km = 200.0
        radius_m = MEAN_EARTH_RADIUS_M

        mesh = mesh_fibonacci_global(spacing_km=spacing_km, radius_m=radius_m)
        self.assertIsInstance(mesh, WgsMesh, msg="Returned type is not WgsMesh")

        n = len(mesh)
        self.assertGreater(n, 0, msg="Global mesh should not be empty")

        # Expected node count from the same analytical formula used
        # in spherical_mesh_global.
        a_target = spacing_km * spacing_km * 1e6
        a_sphere = 4.0 * 3.141592653589793 * (radius_m ** 2)
        n_expected = max(1, int(round(a_sphere / a_target)))

        self.assertEqual(
            n,
            n_expected,
            msg="Global mesh node count mismatch with analytical formula",
        )

        # All points must lie on the sphere and within valid lon/lat bounds.
        for p in mesh:
            self.assertIsInstance(p, WgsPoint)
            self.assertTrue(
                -180.0 <= p.longitude <= 180.0,
                msg="Longitude out of bounds in global mesh",
            )
            self.assertTrue(
                -90.0 <= p.latitude <= 90.0,
                msg="Latitude out of bounds in global mesh",
            )


class TestSphericalCapMesh(unittest.TestCase):
    """
    Tests for spherical_cap_mesh_fibonacci and the bbox variant.
    """

    def test_cap_mesh_basic_properties(self) -> None:
        """
        Generated nodes must lie within the requested cap radius
        (within a small tolerance).
        """
        center_lon = 13.0
        center_lat = 46.0
        max_dist_km = 100.0
        spacing_km = 20.0

        mesh = mesh_fibonacci_spherical_cap(
            center_lon_deg=center_lon,
            center_lat_deg=center_lat,
            max_distance_km=max_dist_km,
            spacing_km=spacing_km,
        )

        self.assertIsInstance(mesh, WgsMesh)
        self.assertGreater(len(mesh), 0, msg="Cap mesh should not be empty")

        for p in mesh:
            d_m = great_circle_distance(
                center_lon,
                center_lat,
                p.longitude,
                p.latitude,
                radius_m=MEAN_EARTH_RADIUS_M,
            )
            d_km = d_m / 1000.0
            self.assertLessEqual(
                d_km,
                max_dist_km * 1.001,
                msg="Node lies outside spherical cap radius",
            )

    def test_cap_mesh_invalid_parameters(self) -> None:
        """spacing_km and max_distance_km must be positive."""
        center_lon = 13.0
        center_lat = 46.0

        with self.assertRaises(ValueError):
            mesh_fibonacci_spherical_cap(
                center_lon_deg=center_lon,
                center_lat_deg=center_lat,
                max_distance_km=0.0,
                spacing_km=20.0,
            )

        with self.assertRaises(ValueError):
            mesh_fibonacci_spherical_cap(
                center_lon_deg=center_lon,
                center_lat_deg=center_lat,
                max_distance_km=100.0,
                spacing_km=0.0,
            )

    def test_cap_mesh_bbox_clipping(self) -> None:
        """
        The bbox variant must return nodes that all lie inside the
        provided bounding box.
        """
        lon_min, lon_max = 12.0, 14.0
        lat_min, lat_max = 45.0, 47.0

        spacing_km = 20.0

        mesh = mesh_fibonacci_bbox(
            lon_min=lon_min,
            lon_max=lon_max,
            lat_min=lat_min,
            lat_max=lat_max,
            spacing_km=spacing_km,
        )

        self.assertIsInstance(mesh, WgsMesh)
        self.assertGreater(len(mesh), 0, msg="BBox cap mesh should not be empty")

        for p in mesh:
            self.assertGreaterEqual(
                p.longitude,
                lon_min - 1e-8,
                msg="Longitude smaller than bbox minimum",
            )
            self.assertLessEqual(
                p.longitude,
                lon_max + 1e-8,
                msg="Longitude greater than bbox maximum",
            )
            self.assertGreaterEqual(
                p.latitude,
                lat_min - 1e-8,
                msg="Latitude smaller than bbox minimum",
            )
            self.assertLessEqual(
                p.latitude,
                lat_max + 1e-8,
                msg="Latitude greater than bbox maximum",
            )

    def test_cap_mesh_bbox_invalid_spacing(self) -> None:
        """spacing_km must be positive in the bbox wrapper."""
        with self.assertRaises(ValueError):
            mesh_fibonacci_bbox(
                lon_min=0.0,
                lon_max=1.0,
                lat_min=0.0,
                lat_max=1.0,
                spacing_km=0.0,
            )


class TestMeshCartesian(unittest.TestCase):
    """
    Tests for the regular Cartesian (lat/lon grid) mesh generator.
    """

    def test_mesh_cartesian_basic(self) -> None:
        """
        mesh_cartesian must generate a regular grid with the expected
        number of points.
        """
        lon_min, lon_max = -1.0, 1.0
        lat_min, lat_max = -1.0, 1.0
        lon_step = 1.0
        lat_step = 1.0

        mesh = mesh_cartesian(
            lon_min=lon_min,
            lon_max=lon_max,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_step_deg=lon_step,
            lat_step_deg=lat_step,
        )

        self.assertIsInstance(mesh, WgsMesh)

        # Expected number of nodes: 3 longitudes x 3 latitudes = 9
        self.assertEqual(
            len(mesh),
            9,
            msg="mesh_cartesian node count mismatch",
        )

        # Bounding box must match the requested limits.
        lon0, lon1, lat0, lat1 = mesh.bounding_box()
        self.assertTrue(
            isclose(lon0, lon_min, abs_tol=1e-12)
            and isclose(lon1, lon_max, abs_tol=1e-12)
            and isclose(lat0, lat_min, abs_tol=1e-12)
            and isclose(lat1, lat_max, abs_tol=1e-12),
            msg="mesh_cartesian bounding_box mismatch",
        )

    def test_mesh_cartesian_invalid_steps(self) -> None:
        """Step sizes must be positive."""
        with self.assertRaises(ValueError):
            mesh_cartesian(
                lon_min=0.0,
                lon_max=1.0,
                lat_min=0.0,
                lat_max=1.0,
                lon_step_deg=0.0,
                lat_step_deg=1.0,
            )

        with self.assertRaises(ValueError):
            mesh_cartesian(
                lon_min=0.0,
                lon_max=1.0,
                lat_min=0.0,
                lat_max=1.0,
                lon_step_deg=1.0,
                lat_step_deg=-1.0,
            )


class TestMeshSinusoidal(unittest.TestCase):
    """
    Tests for the equal-area sinusoidal mesh generator.
    """

    def test_sinusoidal_invalid_spacing(self) -> None:
        """spacing_km must be strictly positive."""
        with self.assertRaises(ValueError):
            mesh_sinusoidal(spacing_km=0.0)

        with self.assertRaises(ValueError):
            mesh_sinusoidal(spacing_km=-10.0)

    def test_sinusoidal_global_basic(self) -> None:
        """
        Global sinusoidal mesh must be non-empty and cover the full
        longitude/latitude extent within reasonable bounds.
        """
        spacing_km = 500.0

        mesh = mesh_sinusoidal(spacing_km=spacing_km)
        self.assertIsInstance(mesh, WgsMesh)
        self.assertGreater(
            len(mesh),
            0,
            msg="Global sinusoidal mesh should not be empty",
        )

        lon_min, lon_max, lat_min, lat_max = mesh.bounding_box()

        # Con tolleranze generose, dato che la mesh può non centrare
        # esattamente gli estremi ±180/±90.
        self.assertLessEqual(
            lon_min,
            -150.0,
            msg="Global sinusoidal mesh does not extend far enough west",
        )
        self.assertGreaterEqual(
            lon_max,
            150.0,
            msg="Global sinusoidal mesh does not extend far enough east",
        )
        self.assertLessEqual(
            lat_min,
            -80.0,
            msg="Global sinusoidal mesh does not extend far enough south",
        )
        self.assertGreaterEqual(
            lat_max,
            80.0,
            msg="Global sinusoidal mesh does not extend far enough north",
        )

    def test_sinusoidal_bbox_subset(self) -> None:
        """
        A mesh generated over a geographic subset must be non-empty,
        bounded within the requested limits and contain fewer nodes
        than the corresponding global mesh (for the same spacing).
        """
        spacing_km = 300.0

        # Global mesh for reference
        global_mesh = mesh_sinusoidal(spacing_km=spacing_km)

        lon_min, lon_max = 10.0, 20.0
        lat_min, lat_max = 40.0, 50.0

        bbox_mesh = mesh_sinusoidal(
            spacing_km=spacing_km,
            lon_min=lon_min,
            lon_max=lon_max,
            lat_min=lat_min,
            lat_max=lat_max,
        )

        self.assertIsInstance(bbox_mesh, WgsMesh)
        self.assertGreater(
            len(bbox_mesh),
            0,
            msg="Sinusoidal bbox mesh should not be empty",
        )

        # Tutti i punti devono stare dentro (con piccola tolleranza)
        for p in bbox_mesh:
            self.assertGreaterEqual(
                p.longitude,
                lon_min - 1e-6,
                msg="Longitude smaller than bbox minimum",
            )
            self.assertLessEqual(
                p.longitude,
                lon_max + 1e-6,
                msg="Longitude greater than bbox maximum",
            )
            self.assertGreaterEqual(
                p.latitude,
                lat_min - 1e-6,
                msg="Latitude smaller than bbox minimum",
            )
            self.assertLessEqual(
                p.latitude,
                lat_max + 1e-6,
                msg="Latitude greater than bbox maximum",
            )

        # Il sottoinsieme deve avere meno nodi del globale.
        self.assertLess(
            len(bbox_mesh),
            len(global_mesh),
            msg="BBox sinusoidal mesh should have fewer nodes than "
            "the global mesh for the same spacing.",
        )


if __name__ == "__main__":
    unittest.main()
