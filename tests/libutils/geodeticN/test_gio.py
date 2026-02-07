#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for shakelab.libutils.geodeticN.gio.

These tests validate the basic behaviour of the GeoJSON I/O helpers:

* read_geojson
* write_geojson

Run with:
    python -m unittest test_gio
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from typing import List

from shakelab.libutils.geodeticN import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
    PolygonKind,
)
from shakelab.libutils.geodeticN.gio import (
    read_geojson,
    write_geojson,
)


class TestGio(unittest.TestCase):
    """
    Test suite for GeoJSON I/O helpers.
    """

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _write_temp_geojson(self, data: dict) -> Path:
        """
        Write a GeoJSON dictionary to a temporary file and return path.
        """
        tmp = tempfile.NamedTemporaryFile(
            suffix=".geojson",
            delete=False,
        )
        path = Path(tmp.name)
        tmp.close()
        with path.open("w", encoding="utf-8") as fobj:
            json.dump(data, fobj)
        return path

    # ------------------------------------------------------------------
    # Reading tests
    # ------------------------------------------------------------------

    def test_read_point_feature(self) -> None:
        """
        Read a Feature with a Point geometry.
        """
        data = {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [12.0, 45.0],
            },
            "properties": {"name": "test"},
        }

        path = self._write_temp_geojson(data)
        try:
            objs = read_geojson(path)
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual(len(objs), 1)
        pt = objs[0]
        self.assertIsInstance(pt, WgsPoint)
        self.assertAlmostEqual(pt.longitude, 12.0)
        self.assertAlmostEqual(pt.latitude, 45.0)

    def test_read_feature_collection_mixed(self) -> None:
        """
        Read a FeatureCollection with different geometry types.
        """
        data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": [0.0, 0.0],
                    },
                    "properties": {},
                },
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [
                            [0.0, 0.0],
                            [1.0, 1.0],
                        ],
                    },
                    "properties": {},
                },
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [0.0, 0.0],
                            [1.0, 0.0],
                            [1.0, 1.0],
                            [0.0, 1.0],
                            [0.0, 0.0],
                        ]],
                    },
                    "properties": {},
                },
            ],
        }

        path = self._write_temp_geojson(data)
        try:
            objs = read_geojson(path)
        finally:
            path.unlink(missing_ok=True)

        # Expected: Point -> WgsPoint, LineString -> WgsMesh,
        # Polygon -> WgsPolygon
        self.assertEqual(len(objs), 3)
        self.assertIsInstance(objs[0], WgsPoint)
        self.assertIsInstance(objs[1], WgsMesh)
        self.assertIsInstance(objs[2], WgsPolygon)

    def test_read_multipolygon(self) -> None:
        """
        Read a MultiPolygon and check it is mapped to a WgsPolygon.
        """
        data = {
            "type": "Feature",
            "geometry": {
                "type": "MultiPolygon",
                "coordinates": [
                    [[
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [1.0, 1.0],
                        [0.0, 1.0],
                        [0.0, 0.0],
                    ]],
                    [[
                        [2.0, 2.0],
                        [3.0, 2.0],
                        [3.0, 3.0],
                        [2.0, 3.0],
                        [2.0, 2.0],
                    ]],
                ],
            },
            "properties": {},
        }

        path = self._write_temp_geojson(data)
        try:
            objs = read_geojson(path)
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual(len(objs), 1)
        poly = objs[0]
        self.assertIsInstance(poly, WgsPolygon)
        self.assertEqual(poly.kind, PolygonKind.MULTIPOLYGON)
        self.assertEqual(poly.n_parts(), 2)

    # ------------------------------------------------------------------
    # Writing tests
    # ------------------------------------------------------------------

    def test_write_point_and_polygon(self) -> None:
        """
        Write WgsPoint and WgsPolygon and validate GeoJSON structure.
        """
        pt = WgsPoint(longitude=12.0, latitude=45.0)
        poly_coords = [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        poly = WgsPolygon.from_lonlat(poly_coords)

        tmp = tempfile.NamedTemporaryFile(
            suffix=".geojson",
            delete=False,
        )
        path = Path(tmp.name)
        tmp.close()

        try:
            write_geojson(path, [pt, poly])
            with path.open("r", encoding="utf-8") as fobj:
                data = json.load(fobj)
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual(data.get("type"), "FeatureCollection")
        features: List[dict] = data.get("features", [])
        self.assertEqual(len(features), 2)

        geom0 = features[0].get("geometry", {})
        self.assertEqual(geom0.get("type"), "Point")
        self.assertEqual(geom0.get("coordinates"), [12.0, 45.0])

        geom1 = features[1].get("geometry", {})
        self.assertEqual(geom1.get("type"), "Polygon")
        coords = geom1.get("coordinates")
        self.assertIsInstance(coords, list)
        self.assertGreater(len(coords), 0)

    def test_polygon_ring_closure_on_write(self) -> None:
        """
        GeoJSON rings written for polygons must be explicitly closed.
        """
        coords = [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        poly = WgsPolygon.from_lonlat(coords)

        tmp = tempfile.NamedTemporaryFile(
            suffix=".geojson",
            delete=False,
        )
        path = Path(tmp.name)
        tmp.close()

        try:
            write_geojson(path, [poly])
            with path.open("r", encoding="utf-8") as fobj:
                data = json.load(fobj)
        finally:
            path.unlink(missing_ok=True)

        features = data.get("features", [])
        self.assertEqual(len(features), 1)

        geom = features[0].get("geometry", {})
        self.assertEqual(geom.get("type"), "Polygon")
        rings = geom.get("coordinates", [])
        self.assertGreater(len(rings), 0)

        ring = rings[0]
        self.assertGreaterEqual(len(ring), 5)
        self.assertEqual(ring[0], ring[-1])

    def test_write_empty_mesh_is_skipped(self) -> None:
        """
        Empty WgsMesh should not produce a geometry in the output.
        """
        mesh = WgsMesh(points=[])

        tmp = tempfile.NamedTemporaryFile(
            suffix=".geojson",
            delete=False,
        )
        path = Path(tmp.name)
        tmp.close()

        try:
            write_geojson(path, [mesh])
            with path.open("r", encoding="utf-8") as fobj:
                data = json.load(fobj)
        finally:
            path.unlink(missing_ok=True)

        features = data.get("features", [])
        self.assertEqual(len(features), 0)


if __name__ == "__main__":
    unittest.main()
