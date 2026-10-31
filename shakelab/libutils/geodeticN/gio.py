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
Input and output of geodetic objects (GeoJSON and CSV).

GeoJSON geometries are converted to and from :class:`WgsPoint`,
:class:`WgsMesh` and :class:`WgsPolygon`; point meshes can also be read
from and written to CSV files.

Design choices
--------------
- Only the geometry is handled: CRS, styling and feature properties are
  intentionally ignored, so that the geometry classes stay independent
  of any specific data model.
- Polygon holes are ignored: only exterior rings are read.
- Elevations are read from the third coordinate when present and are
  written only for geometries with some non-zero elevation.

Conventions
-----------
- Coordinates are ``(lon, lat[, elevation])``, longitude and latitude
  in degrees and elevation in meters.
- CSV files have a header; columns are identified by name
  (``longitude``/``lon``, ``latitude``/``lat`` and optionally
  ``elevation``/``elev``), case-insensitively and in any order.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

from shakelab.libutils.geodeticN.primitives import (
    PolygonKind,
    WgsMesh,
    WgsPoint,
    WgsPolygon,
)


__all__ = [
    "read_geojson",
    "write_geojson",
    "read_csv",
    "write_csv",
]


GeoObject = Union[WgsPoint, WgsPolygon, WgsMesh]
PathLike = Union[str, Path]

_LONGITUDE_NAMES = ("longitude", "lon")
_LATITUDE_NAMES = ("latitude", "lat")
_ELEVATION_NAMES = ("elevation", "elev")


# ---------------------------------------------------------------------------
# GeoJSON
# ---------------------------------------------------------------------------

def read_geojson(path: PathLike) -> List[GeoObject]:
    """
    Read geometries from a GeoJSON file.

    The file may contain a Geometry, a Feature or a FeatureCollection.
    Supported geometries are mapped as follows:

    - Point → :class:`WgsPoint`;
    - MultiPoint, LineString → :class:`WgsMesh`;
    - MultiLineString → one :class:`WgsMesh` per line;
    - Polygon → :class:`WgsPolygon` (exterior ring);
    - MultiPolygon → :class:`WgsPolygon` of kind MULTIPOLYGON.

    Other geometry types are skipped.

    Parameters
    ----------
    path
        GeoJSON file path.

    Returns
    -------
    list
        Geometries in file order.
    """
    with Path(path).open("r", encoding="utf-8") as fobj:
        data = json.load(fobj)

    geometries: List[Dict[str, Any]] = []
    gtype = data.get("type")

    if gtype == "FeatureCollection":
        for feature in data.get("features", []):
            if feature.get("geometry") is not None:
                geometries.append(feature["geometry"])
    elif gtype == "Feature":
        if data.get("geometry") is not None:
            geometries.append(data["geometry"])
    else:
        geometries.append(data)

    result: List[GeoObject] = []

    for geometry in geometries:
        result.extend(_parse_geometry(geometry))

    return result


def write_geojson(path: PathLike, geometries: Iterable[GeoObject]) -> None:
    """
    Write geometries to a GeoJSON FeatureCollection.

    Objects are converted as follows:

    - :class:`WgsPoint` → Point;
    - :class:`WgsMesh` → Point (one node) or MultiPoint;
    - :class:`WgsPolygon` → Polygon, or MultiPolygon for kind
      MULTIPOLYGON or several parts.

    Empty objects are skipped; feature properties are empty.

    Parameters
    ----------
    path
        Output file path.
    geometries
        Objects to write.
    """
    features = []

    for obj in geometries:
        geometry = _geometry_to_geojson(obj)

        if geometry is not None:
            features.append({
                "type": "Feature",
                "geometry": geometry,
                "properties": {},
            })

    collection = {"type": "FeatureCollection", "features": features}

    with Path(path).open("w", encoding="utf-8") as fobj:
        json.dump(collection, fobj, indent=2)


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def read_csv(path: PathLike, delimiter: str = ",") -> WgsMesh:
    """
    Read a point mesh from a CSV file with a header.

    Parameters
    ----------
    path
        CSV file path.
    delimiter
        Field delimiter.

    Returns
    -------
    WgsMesh
        Points in file order. Missing elevations are zero.

    Raises
    ------
    ValueError
        If the longitude or latitude column is missing.
    """
    with Path(path).open("r", encoding="utf-8", newline="") as fobj:
        reader = csv.reader(fobj, delimiter=delimiter)
        header = [name.strip().lower() for name in next(reader, [])]

        lon_index = _column_index(header, _LONGITUDE_NAMES)
        lat_index = _column_index(header, _LATITUDE_NAMES)
        elev_index = _column_index(header, _ELEVATION_NAMES, required=False)

        points = []

        for row in reader:
            if not row or all(not field.strip() for field in row):
                continue

            elevation = 0.0

            if elev_index is not None and row[elev_index].strip():
                elevation = float(row[elev_index])

            points.append(WgsPoint(
                longitude=float(row[lon_index]),
                latitude=float(row[lat_index]),
                elevation=elevation,
            ))

    return WgsMesh(points=points)


def write_csv(
    path: PathLike,
    mesh: Union[WgsMesh, Sequence[WgsPoint]],
    delimiter: str = ",",
) -> None:
    """
    Write a point mesh to a CSV file.

    The columns are ``longitude``, ``latitude`` and ``elevation``.

    Parameters
    ----------
    path
        Output file path.
    mesh
        Mesh or sequence of points.
    delimiter
        Field delimiter.
    """
    with Path(path).open("w", encoding="utf-8", newline="") as fobj:
        writer = csv.writer(fobj, delimiter=delimiter)
        writer.writerow(["longitude", "latitude", "elevation"])

        for point in mesh:
            writer.writerow([
                repr(point.longitude),
                repr(point.latitude),
                repr(point.elevation),
            ])


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _parse_geometry(geometry: Dict[str, Any]) -> List[GeoObject]:
    """
    Convert one GeoJSON geometry to geodetic objects.
    """
    gtype = geometry.get("type")
    coords = geometry.get("coordinates")

    if gtype is None or coords is None:
        return []

    if gtype == "Point":
        return [_point(coords)]

    if gtype in ("MultiPoint", "LineString"):
        return [WgsMesh(points=[_point(c) for c in coords])]

    if gtype == "MultiLineString":
        return [WgsMesh(points=[_point(c) for c in line]) for line in coords]

    if gtype == "Polygon":
        if not coords:
            return []

        return [WgsPolygon(parts=[[_point(c) for c in coords[0]]])]

    if gtype == "MultiPolygon":
        parts = [
            [_point(c) for c in polygon[0]]
            for polygon in coords
            if polygon
        ]

        if not parts:
            return []

        return [WgsPolygon(parts=parts, kind=PolygonKind.MULTIPOLYGON)]

    return []


def _point(coord: Sequence[float]) -> WgsPoint:
    """
    Convert a GeoJSON position to a :class:`WgsPoint`.
    """
    if len(coord) < 2:
        raise ValueError("GeoJSON positions need at least two values.")

    elevation = float(coord[2]) if len(coord) > 2 else 0.0

    return WgsPoint(
        longitude=float(coord[0]),
        latitude=float(coord[1]),
        elevation=elevation,
    )


def _geometry_to_geojson(obj: GeoObject) -> Optional[Dict[str, Any]]:
    """
    Convert a geodetic object to a GeoJSON geometry, or None.
    """
    if isinstance(obj, WgsPoint):
        return {"type": "Point", "coordinates": _positions([obj])[0]}

    if isinstance(obj, WgsMesh):
        if len(obj) == 0:
            return None

        positions = _positions(obj.points)

        if len(positions) == 1:
            return {"type": "Point", "coordinates": positions[0]}

        return {"type": "MultiPoint", "coordinates": positions}

    if isinstance(obj, WgsPolygon):
        if obj.is_empty():
            return None

        with_z = any(p.elevation != 0.0 for p in obj.iter_vertices())
        rings = [_ring(part, with_z) for part in obj.parts]

        if obj.kind == PolygonKind.MULTIPOLYGON or len(rings) > 1:
            return {
                "type": "MultiPolygon",
                "coordinates": [[ring] for ring in rings],
            }

        return {"type": "Polygon", "coordinates": [rings[0]]}

    return None


def _positions(
    points: Sequence[WgsPoint],
    with_z: Optional[bool] = None,
) -> List[List[float]]:
    """
    Return GeoJSON positions, 3-D if any elevation is non-zero.
    """
    if with_z is None:
        with_z = any(p.elevation != 0.0 for p in points)

    if with_z:
        return [[p.longitude, p.latitude, p.elevation] for p in points]

    return [[p.longitude, p.latitude] for p in points]


def _ring(points: Sequence[WgsPoint], with_z: bool) -> List[List[float]]:
    """
    Return a closed GeoJSON linear ring.
    """
    positions = _positions(points, with_z)

    if positions and positions[0] != positions[-1]:
        positions.append(list(positions[0]))

    return positions


def _column_index(
    header: Sequence[str],
    names: Sequence[str],
    required: bool = True,
) -> Optional[int]:
    """
    Return the index of the first header field matching ``names``.
    """
    for name in names:
        if name in header:
            return header.index(name)

    if required:
        raise ValueError(f"CSV header lacks a '{names[0]}' column.")

    return None
