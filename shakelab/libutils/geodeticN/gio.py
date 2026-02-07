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
Input/output utilities for geodetic objects in ShakeLab.

This module provides simple helpers to read and write GeoJSON files
using the WgsPoint, WgsPolygon and WgsMesh classes defined in the
:mod:`geodeticN.geometry` module.

The focus is on geometry only: CRS handling, styling and extra fields
are intentionally ignored or kept minimal. All coordinates are
interpreted as (longitude, latitude) in decimal degrees.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List, Sequence, Union

from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
    PolygonKind,
)


GeoObject = Union[WgsPoint, WgsPolygon, WgsMesh]


# ============================================================================
# Reading GeoJSON
# ============================================================================


def read_geojson(path: Union[str, Path]) -> List[GeoObject]:
    """
    Read geometries from a GeoJSON file.

    The file may contain a Geometry, a Feature or a
    FeatureCollection. Geometries supported are Point, MultiPoint,
    LineString, MultiLineString, Polygon and MultiPolygon. Unsupported
    types are silently skipped.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the GeoJSON file.

    Returns
    -------
    list of GeoObject
        List of WgsPoint, WgsMesh or WgsPolygon instances parsed from
        the file.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as fobj:
        data = json.load(fobj)

    result: List[GeoObject] = []

    gtype = data.get("type")
    if gtype == "FeatureCollection":
        features = data.get("features", [])
        for feat in features:
            geom = feat.get("geometry")
            if geom is None:
                continue
            result.extend(_parse_geometry(geom))
    elif gtype == "Feature":
        geom = data.get("geometry")
        if geom is not None:
            result.extend(_parse_geometry(geom))
    else:
        # Assume a raw Geometry object.
        result.extend(_parse_geometry(data))

    return result


def _parse_geometry(geom: dict) -> List[GeoObject]:
    """
    Internal helper to parse a single GeoJSON geometry object.

    Parameters
    ----------
    geom : dict
        GeoJSON geometry dictionary.

    Returns
    -------
    list of GeoObject
        Parsed objects (0 or more). Multi-geometries may expand into
        several objects.
    """
    gtype = geom.get("type")
    coords = geom.get("coordinates")

    if gtype is None or coords is None:
        return []

    out: List[GeoObject] = []

    if gtype == "Point":
        lon, lat = coords
        out.append(WgsPoint(longitude=float(lon), latitude=float(lat)))

    elif gtype == "MultiPoint":
        pts = [(float(c[0]), float(c[1])) for c in coords]
        out.append(WgsMesh.from_lonlat(pts))

    elif gtype == "LineString":
        pts = [(float(c[0]), float(c[1])) for c in coords]
        out.append(WgsMesh.from_lonlat(pts))

    elif gtype == "MultiLineString":
        # Each line string is returned as a separate mesh.
        for line in coords:
            pts = [(float(c[0]), float(c[1])) for c in line]
            out.append(WgsMesh.from_lonlat(pts))

    elif gtype == "Polygon":
        # GeoJSON: list of linear rings. We only use the exterior ring
        # for now and ignore possible holes.
        if not coords:
            return []
        exterior = coords[0]
        pts = [(float(c[0]), float(c[1])) for c in exterior]
        out.append(WgsPolygon.from_lonlat(pts))

    elif gtype == "MultiPolygon":
        # Each polygon has its own list of rings. We consider only the
        # exterior ring of each polygon and represent them as parts of
        # a single WgsPolygon (kind=MULTIPOLYGON).
        parts: List[Sequence[tuple[float, float]]] = []
        for poly in coords:
            if not poly:
                continue
            exterior = poly[0]
            pts = [(float(c[0]), float(c[1])) for c in exterior]
            parts.append(pts)

        if parts:
            out.append(WgsPolygon.from_parts(parts))

    # Other geometry types (e.g. GeometryCollection) are ignored.
    return out


# ============================================================================
# Writing GeoJSON
# ============================================================================


def write_geojson(
    path: Union[str, Path],
    geometries: Iterable[GeoObject],
) -> None:
    """
    Write a list of geodetic objects to a GeoJSON file.

    The output is a FeatureCollection where each object is converted
    to a corresponding geometry:

    * WgsPoint      -> Point
    * WgsMesh       -> Point or MultiPoint (depending on size)
    * WgsPolygon    -> Polygon or MultiPolygon (depending on kind)

    Properties are left empty (``{}``) for all features.

    Parameters
    ----------
    path : str or pathlib.Path
        Output file path.
    geometries : iterable of GeoObject
        Objects to be written.
    """
    path = Path(path)

    features = []
    for obj in geometries:
        geom = _geometry_to_geojson(obj)
        if geom is None:
            continue
        feature = {
            "type": "Feature",
            "geometry": geom,
            "properties": {},
        }
        features.append(feature)

    collection = {
        "type": "FeatureCollection",
        "features": features,
    }

    with path.open("w", encoding="utf-8") as fobj:
        json.dump(collection, fobj, indent=2)


def _geometry_to_geojson(obj: GeoObject) -> dict | None:
    """
    Internal helper to convert a geodetic object to GeoJSON geometry.

    Parameters
    ----------
    obj : GeoObject
        WgsPoint, WgsPolygon or WgsMesh instance.

    Returns
    -------
    dict or None
        GeoJSON geometry dictionary, or None if the object cannot be
        represented.
    """
    if isinstance(obj, WgsPoint):
        return {
            "type": "Point",
            "coordinates": [obj.longitude, obj.latitude],
        }

    if isinstance(obj, WgsMesh):
        if len(obj) == 0:
            return None
        if len(obj) == 1:
            p = obj.points[0]
            return {
                "type": "Point",
                "coordinates": [p.longitude, p.latitude],
            }
        coords = [[p.longitude, p.latitude] for p in obj.points]
        return {
            "type": "MultiPoint",
            "coordinates": coords,
        }

    if isinstance(obj, WgsPolygon):
        if obj.is_empty():
            return None

        if obj.kind == PolygonKind.SIMPLE:
            ring = _ring_to_geojson(obj.parts[0])
            return {
                "type": "Polygon",
                "coordinates": [ring],
            }

        # Multipolygon: one ring per part, each as an exterior ring.
        polys = []
        for part in obj.parts:
            ring = _ring_to_geojson(part)
            polys.append([ring])

        return {
            "type": "MultiPolygon",
            "coordinates": polys,
        }

    return None


def _ring_to_geojson(points: Sequence[WgsPoint]) -> List[List[float]]:
    """
    Convert a polygon ring to a GeoJSON linear ring.

    A closing vertex is added if needed.

    Parameters
    ----------
    points : sequence of WgsPoint
        Ring vertices.

    Returns
    -------
    list of [lon, lat]
        GeoJSON coordinate list.
    """
    if not points:
        return []

    coords = [[p.longitude, p.latitude] for p in points]

    # Ensure ring closure (first == last).
    if coords[0] != coords[-1]:
        coords.append(coords[0])

    return coords
