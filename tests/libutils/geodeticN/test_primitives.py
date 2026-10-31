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
Unit tests for :mod:`shakelab.libutils.geodeticN.primitives`.

The tests cover the validation and the methods of :class:`WgsPoint`,
the geometry of :class:`WgsPolygon` (including multipolygons, closed
rings, the antimeridian and the ``[0, 360)`` convention) and the
:class:`WgsMesh` container.

Conventions
-----------
- ``SQUARE`` is the 2° x 2° box centered at (0°, 0°).
"""

from __future__ import annotations

import dataclasses
from decimal import Decimal
import math

import numpy as np
import pytest

from shakelab.libutils.geodeticN.primitives import (
    PolygonKind,
    WgsMesh,
    WgsPoint,
    WgsPolygon,
)
from shakelab.libutils.geodeticN.transform import MetricPoint
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
)


SQUARE = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]


def _box_area(dlon_deg, lat1_deg, lat2_deg):
    """Return the exact spherical area of a longitude-latitude box."""
    return (
        MEAN_EARTH_RADIUS_M ** 2
        * math.radians(dlon_deg)
        * (math.sin(math.radians(lat2_deg)) - math.sin(math.radians(lat1_deg)))
    )


# ---------------------------------------------------------------------------
# WgsPoint
# ---------------------------------------------------------------------------

def test_point_attributes_and_tuple():
    """Coordinates are stored as floats and exported as (lon, lat)."""
    point = WgsPoint(longitude=12, latitude=np.float32(45.0), elevation=100)
    assert (point.longitude, point.latitude, point.elevation) == (
        12.0, 45.0, 100.0,
    )
    assert all(
        type(value) is float
        for value in (point.longitude, point.latitude, point.elevation)
    )
    assert point.to_tuple() == (12.0, 45.0)


@pytest.mark.parametrize("longitude", [350.0, -190.0, 720.0])
def test_point_keeps_longitude_convention(longitude):
    """Longitudes outside [-180, 180) are kept as given."""
    assert WgsPoint(longitude=longitude, latitude=0.0).longitude == longitude


def test_point_accepts_numeric_types():
    """Decimal and NumPy scalars are accepted."""
    point = WgsPoint(Decimal("12.5"), np.int64(45), np.float64(-10.0))
    assert point == WgsPoint(12.5, 45.0, -10.0)


@pytest.mark.parametrize("kwargs, error", [
    ({"longitude": 0.0, "latitude": 90.5}, ValueError),
    ({"longitude": 0.0, "latitude": -91.0}, ValueError),
    ({"longitude": math.nan, "latitude": 0.0}, ValueError),
    ({"longitude": 0.0, "latitude": 0.0, "elevation": math.inf},
     ValueError),
    ({"longitude": "12", "latitude": 0.0}, TypeError),
    ({"longitude": True, "latitude": 0.0}, TypeError),
    ({"longitude": None, "latitude": 0.0}, TypeError),
])
def test_point_validation(kwargs, error):
    """Invalid coordinates are rejected."""
    with pytest.raises(error):
        WgsPoint(**kwargs)


def test_point_is_immutable_and_hashable():
    """Points cannot be modified and can be used in sets."""
    point = WgsPoint(12.0, 45.0)
    with pytest.raises(dataclasses.FrozenInstanceError):
        point.latitude = 46.0
    assert len({point, WgsPoint(12.0, 45.0)}) == 1


def test_point_epicentral_distance():
    """epicentral_distance_to is the great-circle distance."""
    distance = WgsPoint(12.0, 45.0).epicentral_distance_to(
        WgsPoint(13.0, 45.0)
    )
    assert distance == pytest.approx(
        great_circle_distance(12.0, 45.0, 13.0, 45.0), rel=1e-12
    )


def test_point_hypocentral_distance_and_subtraction():
    """The subtraction operator is the hypocentral distance."""
    hypocentre = WgsPoint(12.0, 45.0, -5000.0)
    site = WgsPoint(12.1, 45.1, 0.0)
    assert hypocentre - site == pytest.approx(
        hypocentre.hypocentral_distance_to(site), rel=1e-14
    )
    sphere = hypocentre.hypocentral_distance_to(site, approx="sphere")
    ellipsoid = hypocentre.hypocentral_distance_to(site, approx="ellipsoid")
    assert abs(sphere - ellipsoid) / ellipsoid < 0.01


def test_point_subtraction_with_other_types():
    """Subtracting a non-point is not supported."""
    with pytest.raises(TypeError):
        WgsPoint(12.0, 45.0) - 1.0


def test_point_to_metric():
    """to_metric returns ECEF coordinates by default."""
    metric = WgsPoint(0.0, 0.0).to_metric()
    assert isinstance(metric, MetricPoint)
    assert metric.x == pytest.approx(6378137.0)


# ---------------------------------------------------------------------------
# WgsPolygon
# ---------------------------------------------------------------------------

def test_polygon_basic_geometry():
    """Bounding box, centroid, area and containment of a square."""
    polygon = WgsPolygon.from_lonlat(SQUARE)
    assert not polygon.is_empty()
    assert polygon.n_parts() == 1
    assert polygon.kind == PolygonKind.SIMPLE
    assert polygon.bounding_box() == (-1.0, 1.0, -1.0, 1.0)

    centroid = polygon.centroid()
    assert centroid.longitude == pytest.approx(0.0, abs=1e-12)
    assert centroid.latitude == pytest.approx(0.0, abs=1e-12)

    assert polygon.area() == pytest.approx(_box_area(2.0, -1.0, 1.0),
                                           rel=1e-3)
    assert polygon.contains_point(WgsPoint(0.0, 0.0))
    assert not polygon.contains_point(WgsPoint(10.0, 10.0))


def test_polygon_kind_is_inferred():
    """The kind follows the number of parts unless given."""
    square = [WgsPoint(*c) for c in SQUARE]
    shifted = [WgsPoint(c[0] + 5.0, c[1]) for c in SQUARE]
    assert WgsPolygon(parts=[square, shifted]).kind == PolygonKind.MULTIPOLYGON
    assert WgsPolygon(parts=[square]).kind == PolygonKind.SIMPLE
    assert WgsPolygon().kind == PolygonKind.SIMPLE
    single = WgsPolygon.from_parts([SQUARE], kind=PolygonKind.MULTIPOLYGON)
    assert single.kind == PolygonKind.MULTIPOLYGON


def test_polygon_rejects_inconsistent_input():
    """A SIMPLE kind with several parts and non-point vertices fail."""
    square = [WgsPoint(*c) for c in SQUARE]
    with pytest.raises(ValueError, match="SIMPLE"):
        WgsPolygon(parts=[square, square], kind=PolygonKind.SIMPLE)
    with pytest.raises(TypeError, match="WgsPoint"):
        WgsPolygon(parts=[[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]])


def test_polygon_closed_ring_is_equivalent():
    """Repeating the first vertex does not change area or centroid."""
    open_polygon = WgsPolygon.from_lonlat(SQUARE)
    closed_polygon = WgsPolygon.from_lonlat(SQUARE + [SQUARE[0]])
    assert closed_polygon.area() == pytest.approx(open_polygon.area())
    assert closed_polygon.centroid() == open_polygon.centroid()


def test_polygon_centroid_is_area_weighted():
    """The centroid of an L-shape is not the mean of its vertices."""
    shape = WgsPolygon.from_lonlat([
        (0.0, 0.0), (2.0, 0.0), (2.0, 0.2), (0.2, 0.2), (0.2, 2.0),
        (0.0, 2.0),
    ])
    centroid = shape.centroid()
    assert centroid.longitude == pytest.approx(0.5737, abs=1e-3)
    assert centroid.latitude == pytest.approx(0.5737, abs=1e-3)


def test_polygon_degenerate_centroid():
    """Zero-area polygons fall back to the mean of the vertices."""
    line = WgsPolygon.from_lonlat([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
    assert line.area() == 0.0
    assert line.centroid() == WgsPoint(1.0, 0.0)


def test_multipolygon_geometry():
    """Area adds up and containment works in every part."""
    polygon = WgsPolygon.from_parts([
        SQUARE,
        [(c[0] + 10.0, c[1]) for c in SQUARE],
    ])
    assert polygon.kind == PolygonKind.MULTIPOLYGON
    assert polygon.area() == pytest.approx(
        2.0 * WgsPolygon.from_lonlat(SQUARE).area(), rel=1e-12
    )
    assert polygon.contains_point(WgsPoint(0.0, 0.0))
    assert polygon.contains_point(WgsPoint(10.0, 0.0))
    assert not polygon.contains_point(WgsPoint(5.0, 0.0))
    assert polygon.centroid().longitude == pytest.approx(5.0, abs=1e-9)


def test_polygon_contains_points_vectorized():
    """contains_points returns one boolean per point."""
    polygon = WgsPolygon.from_lonlat(SQUARE)
    inside = polygon.contains_points([0.0, 0.5, 2.0], [0.0, -0.5, 0.0])
    np.testing.assert_array_equal(inside, [True, True, False])


def test_polygon_across_antimeridian():
    """A box from 170° to -170° behaves like a box around 0°."""
    crossing = WgsPolygon.from_lonlat([
        (170.0, -5.0), (-170.0, -5.0), (-170.0, 5.0), (170.0, 5.0),
    ])
    reference = WgsPolygon.from_lonlat([
        (-10.0, -5.0), (10.0, -5.0), (10.0, 5.0), (-10.0, 5.0),
    ])
    assert crossing.area() == pytest.approx(reference.area(), rel=1e-12)
    assert crossing.contains_point(WgsPoint(180.0, 0.0))
    assert crossing.contains_point(WgsPoint(-175.0, 0.0))
    assert not crossing.contains_point(WgsPoint(0.0, 0.0))
    assert crossing.centroid().longitude == pytest.approx(180.0, abs=1e-9)


def test_polygon_longitude_convention():
    """Vertices and points in [0, 360) give the same results."""
    mixed = WgsPolygon.from_lonlat([
        (350.0, -1.0), (10.0, -1.0), (10.0, 1.0), (350.0, 1.0),
    ])
    reference = WgsPolygon.from_lonlat([
        (-10.0, -1.0), (10.0, -1.0), (10.0, 1.0), (-10.0, 1.0),
    ])
    assert mixed.area() == pytest.approx(reference.area(), rel=1e-12)
    assert mixed.contains_point(WgsPoint(0.0, 0.0))
    assert mixed.contains_point(WgsPoint(360.0, 0.0))
    assert mixed.contains_point(WgsPoint(355.0, 0.0))


def test_polygon_to_arrays():
    """to_arrays returns the (lon, lat) arrays of each part."""
    (lon, lat), = WgsPolygon.from_lonlat(SQUARE).to_arrays()
    np.testing.assert_array_equal(lon, [-1.0, 1.0, 1.0, -1.0])
    np.testing.assert_array_equal(lat, [-1.0, -1.0, 1.0, 1.0])


def test_empty_polygon():
    """Empty polygons have zero area and contain nothing."""
    polygon = WgsPolygon()
    assert polygon.is_empty()
    assert polygon.area() == 0.0
    assert not polygon.contains_point(WgsPoint(0.0, 0.0))
    with pytest.raises(ValueError):
        polygon.bounding_box()
    with pytest.raises(ValueError):
        polygon.centroid()


# ---------------------------------------------------------------------------
# WgsMesh
# ---------------------------------------------------------------------------

def test_mesh_basic_operations():
    """Length, iteration, indexing and bounding box."""
    mesh = WgsMesh.from_lonlat([(0.0, 0.0), (1.0, 2.0), (-1.0, -2.0, 50.0)])
    assert len(mesh) == 3
    assert [p.longitude for p in mesh] == [0.0, 1.0, -1.0]
    assert mesh[2].elevation == 50.0
    assert mesh.bounding_box() == (-1.0, 1.0, -2.0, 2.0)


def test_mesh_arrays_round_trip():
    """to_arrays and from_arrays are inverse."""
    mesh = WgsMesh.from_arrays([0.0, 1.0], [10.0, 11.0], [5.0, -5.0])
    lon, lat, elev = mesh.to_arrays()
    np.testing.assert_array_equal(lon, [0.0, 1.0])
    np.testing.assert_array_equal(lat, [10.0, 11.0])
    np.testing.assert_array_equal(elev, [5.0, -5.0])
    assert WgsMesh.from_arrays(lon, lat).points[0].elevation == 0.0


def test_mesh_from_arrays_rejects_mismatched_lengths():
    """Coordinate arrays must have the same length."""
    with pytest.raises(ValueError, match="same length"):
        WgsMesh.from_arrays([0.0, 1.0], [0.0])


def test_empty_mesh():
    """An empty mesh has no bounding box and clips to an empty mesh."""
    mesh = WgsMesh()
    assert len(mesh) == 0
    with pytest.raises(ValueError):
        mesh.bounding_box()
    assert len(mesh.clip(WgsPolygon.from_lonlat(SQUARE))) == 0


def test_mesh_clip():
    """clip keeps the points inside the polygon, in order."""
    mesh = WgsMesh.from_lonlat([
        (0.0, 0.0), (5.0, 5.0), (0.5, -0.5), (-2.0, 0.0),
    ])
    clipped = mesh.clip(WgsPolygon.from_lonlat(SQUARE))
    assert [p.to_tuple() for p in clipped] == [(0.0, 0.0), (0.5, -0.5)]
    assert len(mesh) == 4


def test_mesh_clip_across_antimeridian():
    """Clipping works for polygons crossing the antimeridian."""
    polygon = WgsPolygon.from_lonlat([
        (170.0, -5.0), (-170.0, -5.0), (-170.0, 5.0), (170.0, 5.0),
    ])
    mesh = WgsMesh.from_lonlat([(175.0, 0.0), (185.0, 0.0), (0.0, 0.0)])
    clipped = mesh.clip(polygon)
    assert [p.longitude for p in clipped] == [175.0, 185.0]
