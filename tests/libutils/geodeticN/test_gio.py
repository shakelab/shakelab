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
Unit tests for :mod:`shakelab.libutils.geodeticN.gio`.

The tests cover GeoJSON reading (all supported geometry types, 3-D
positions, skipped content), GeoJSON writing (structure, ring closure,
elevations, multipolygons, round trips) and CSV input and output.

Conventions
-----------
- Files are written in the ``tmp_path`` directory provided by pytest.
"""

from __future__ import annotations

import json

import pytest

from shakelab.libutils.geodeticN import (
    PolygonKind,
    WgsMesh,
    WgsPoint,
    WgsPolygon,
)
from shakelab.libutils.geodeticN.gio import (
    read_csv,
    read_geojson,
    write_csv,
    write_geojson,
)


SQUARE = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def _feature(geometry):
    """Wrap a geometry in a GeoJSON Feature."""
    return {"type": "Feature", "geometry": geometry, "properties": {}}


def _read(tmp_path, data):
    """Write a GeoJSON dictionary and read it back."""
    path = tmp_path / "input.geojson"
    path.write_text(json.dumps(data), encoding="utf-8")
    return read_geojson(path)


def _write(tmp_path, objects):
    """Write objects to GeoJSON and return the parsed file."""
    path = tmp_path / "output.geojson"
    write_geojson(path, objects)
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Reading GeoJSON
# ---------------------------------------------------------------------------

def test_read_point_feature(tmp_path):
    """A Point feature becomes a WgsPoint."""
    (point,) = _read(tmp_path, _feature({
        "type": "Point", "coordinates": [12.0, 45.0],
    }))
    assert point == WgsPoint(12.0, 45.0)


def test_read_point_with_elevation(tmp_path):
    """A third coordinate is read as elevation."""
    (point,) = _read(tmp_path, {
        "type": "Point", "coordinates": [12.0, 45.0, -5000.0],
    })
    assert point == WgsPoint(12.0, 45.0, -5000.0)


def test_read_feature_collection_mixed(tmp_path):
    """Each geometry type maps to the expected class."""
    objects = _read(tmp_path, {
        "type": "FeatureCollection",
        "features": [
            _feature({"type": "Point", "coordinates": [0.0, 0.0]}),
            _feature({"type": "LineString",
                      "coordinates": [[0.0, 0.0], [1.0, 1.0]]}),
            _feature({"type": "MultiPoint",
                      "coordinates": [[0.0, 0.0], [2.0, 2.0, 10.0]]}),
            _feature({"type": "Polygon",
                      "coordinates": [[list(c) for c in SQUARE]]}),
        ],
    })
    assert [type(obj) for obj in objects] == [
        WgsPoint, WgsMesh, WgsMesh, WgsPolygon,
    ]
    assert objects[2][1].elevation == 10.0


def test_read_multilinestring(tmp_path):
    """Each line of a MultiLineString becomes a mesh."""
    objects = _read(tmp_path, {
        "type": "MultiLineString",
        "coordinates": [[[0.0, 0.0], [1.0, 1.0]], [[2.0, 2.0], [3.0, 3.0]]],
    })
    assert [len(mesh) for mesh in objects] == [2, 2]


def test_read_multipolygon(tmp_path):
    """A MultiPolygon becomes a WgsPolygon of kind MULTIPOLYGON."""
    shifted = [[c[0] + 2.0, c[1] + 2.0] for c in SQUARE]
    (polygon,) = _read(tmp_path, _feature({
        "type": "MultiPolygon",
        "coordinates": [[[list(c) for c in SQUARE]], [shifted]],
    }))
    assert polygon.kind == PolygonKind.MULTIPOLYGON
    assert polygon.n_parts() == 2


def test_read_single_part_multipolygon_keeps_kind(tmp_path):
    """A MultiPolygon with one polygon keeps the MULTIPOLYGON kind."""
    (polygon,) = _read(tmp_path, {
        "type": "MultiPolygon",
        "coordinates": [[[list(c) for c in SQUARE]]],
    })
    assert polygon.kind == PolygonKind.MULTIPOLYGON
    assert polygon.n_parts() == 1


def test_read_skips_unsupported_content(tmp_path):
    """Null geometries and unsupported types are skipped."""
    objects = _read(tmp_path, {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "geometry": None, "properties": {}},
            _feature({"type": "GeometryCollection", "geometries": []}),
            _feature({"type": "Point", "coordinates": [1.0, 2.0]}),
        ],
    })
    assert objects == [WgsPoint(1.0, 2.0)]


# ---------------------------------------------------------------------------
# Writing GeoJSON
# ---------------------------------------------------------------------------

def test_write_point_and_polygon(tmp_path):
    """Points and polygons are written as features with closed rings."""
    data = _write(tmp_path, [
        WgsPoint(12.0, 45.0),
        WgsPolygon.from_lonlat(SQUARE),
    ])
    assert data["type"] == "FeatureCollection"
    point, polygon = (feature["geometry"] for feature in data["features"])
    assert point == {"type": "Point", "coordinates": [12.0, 45.0]}
    assert polygon["type"] == "Polygon"
    ring = polygon["coordinates"][0]
    assert len(ring) == 5
    assert ring[0] == ring[-1]


def test_write_closed_ring_is_not_closed_twice(tmp_path):
    """A ring that already repeats its first vertex is kept as is."""
    data = _write(tmp_path, [WgsPolygon.from_lonlat(SQUARE + [SQUARE[0]])])
    assert len(data["features"][0]["geometry"]["coordinates"][0]) == 5


def test_write_elevation_only_when_present(tmp_path):
    """Positions are 3-D only for geometries with some elevation."""
    data = _write(tmp_path, [
        WgsPoint(12.0, 45.0),
        WgsPoint(12.0, 45.0, -5000.0),
        WgsMesh.from_lonlat([(0.0, 0.0), (1.0, 1.0, 50.0)]),
    ])
    flat, deep, mesh = (feature["geometry"] for feature in data["features"])
    assert flat["coordinates"] == [12.0, 45.0]
    assert deep["coordinates"] == [12.0, 45.0, -5000.0]
    assert mesh["coordinates"] == [[0.0, 0.0, 0.0], [1.0, 1.0, 50.0]]


def test_write_multipart_polygon_as_multipolygon(tmp_path):
    """A polygon with several parts is never truncated to one part."""
    polygon = WgsPolygon.from_parts([
        SQUARE, [(c[0] + 2.0, c[1]) for c in SQUARE],
    ])
    geometry = _write(tmp_path, [polygon])["features"][0]["geometry"]
    assert geometry["type"] == "MultiPolygon"
    assert len(geometry["coordinates"]) == 2


def test_write_meshes(tmp_path):
    """Empty meshes are skipped; single nodes become Points."""
    data = _write(tmp_path, [
        WgsMesh(),
        WgsMesh.from_lonlat([(1.0, 2.0)]),
        WgsMesh.from_lonlat([(1.0, 2.0), (3.0, 4.0)]),
    ])
    types = [feature["geometry"]["type"] for feature in data["features"]]
    assert types == ["Point", "MultiPoint"]


def test_geojson_round_trip(tmp_path):
    """Written geometries are read back unchanged."""
    polygon = WgsPolygon.from_parts([
        SQUARE, [(c[0] + 2.0, c[1], 10.0) for c in SQUARE],
    ])
    point = WgsPoint(12.0, 45.0, -5000.0)
    path = tmp_path / "round.geojson"
    write_geojson(path, [point, polygon])
    point2, polygon2 = read_geojson(path)
    assert point2 == point
    assert polygon2.kind == PolygonKind.MULTIPOLYGON
    assert polygon2.area() == pytest.approx(polygon.area(), rel=1e-12)
    assert polygon2.parts[1][0] == polygon.parts[1][0]


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def test_csv_round_trip(tmp_path):
    """write_csv and read_csv preserve coordinates and elevations."""
    mesh = WgsMesh.from_lonlat([(12.5, 45.25, 100.0), (-0.1, 1e-7)])
    path = tmp_path / "mesh.csv"
    write_csv(path, mesh)
    assert path.read_text().splitlines()[0] == "longitude,latitude,elevation"
    assert read_csv(path).points == mesh.points


def test_csv_columns_by_name(tmp_path):
    """Columns are found by name, in any order and case."""
    path = tmp_path / "legacy.csv"
    path.write_text("Latitude, Longitude\n45.0,12.0\n\n46.0,13.0\n")
    mesh = read_csv(path)
    assert [p.to_tuple() for p in mesh] == [(12.0, 45.0), (13.0, 46.0)]


def test_csv_custom_delimiter(tmp_path):
    """A custom delimiter can be used for reading and writing."""
    mesh = WgsMesh.from_lonlat([(1.0, 2.0, 3.0)])
    path = tmp_path / "mesh.csv"
    write_csv(path, mesh, delimiter=";")
    assert read_csv(path, delimiter=";").points == mesh.points


def test_csv_missing_column(tmp_path):
    """A file without a latitude column is rejected."""
    path = tmp_path / "bad.csv"
    path.write_text("lon,elev\n1.0,2.0\n")
    with pytest.raises(ValueError, match="latitude"):
        read_csv(path)
