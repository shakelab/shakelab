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
Geometric primitives on the WGS84 ellipsoid.

This module defines light-weight containers for points, polygons and
unstructured point meshes in geographic coordinates, with basic
operations (distances, bounding box, centroid, area, point-in-polygon
and clipping).

Design choices
--------------
- :class:`WgsPoint` is immutable and validated: coordinates are finite
  floats and the latitude lies in ``[-90, 90]``. The longitude is kept
  as given, in any convention.
- Polygon area, centroid and point-in-polygon tests are computed in the
  equal-area sinusoidal projection, centered on the mean longitude of
  each part (computed relative to its first vertex). This makes them
  independent of the longitude convention, valid across the
  antimeridian for parts spanning less than 180 degrees of longitude,
  and symmetric for symmetric polygons.
- Polygon edges are straight lines in the projected plane, which is an
  approximation of geodesics adequate for regional polygons.

Conventions
-----------
- Geographic coordinates are ``(longitude, latitude)`` in decimal
  degrees; longitudes may use ``[-180, 180)`` or ``[0, 360)``.
- Elevation is in meters, positive upward (a hypocentre at 10 km depth
  has ``elevation = -10000``).
- Distances are in meters, areas in square meters.
- A polygon ring may or may not repeat its first vertex at the end.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
from enum import IntEnum
import math
from numbers import Real
from typing import Any, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from shakelab.libutils.geodeticN.utilities import (
    relative_longitude,
    sinusoidal_to_wgs84,
    wgs84_to_sinusoidal,
)


__all__ = [
    "WgsPoint",
    "PolygonKind",
    "WgsPolygon",
    "WgsMesh",
]


# ---------------------------------------------------------------------------
# Point
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WgsPoint:
    """
    Geographic point on the WGS84 ellipsoid.

    Parameters
    ----------
    longitude
        Longitude in degrees, in any convention.
    latitude
        Latitude in degrees, in ``[-90, 90]``.
    elevation
        Elevation above the ellipsoid in meters, positive upward.

    Raises
    ------
    TypeError
        If a coordinate is not a real number.
    ValueError
        If a coordinate is not finite or the latitude is out of range.
    """

    longitude: float
    latitude: float
    elevation: float = 0.0

    def __post_init__(self) -> None:
        longitude = _finite_float(self.longitude, "longitude")
        latitude = _finite_float(self.latitude, "latitude")
        elevation = _finite_float(self.elevation, "elevation")

        if not -90.0 <= latitude <= 90.0:
            raise ValueError("The latitude must be in [-90, 90].")

        object.__setattr__(self, "longitude", longitude)
        object.__setattr__(self, "latitude", latitude)
        object.__setattr__(self, "elevation", elevation)

    def to_tuple(self) -> Tuple[float, float]:
        """
        Return the point as a ``(lon, lat)`` tuple.
        """
        return self.longitude, self.latitude

    def to_metric(self, frame: Any = None) -> Any:
        """
        Convert the point to Cartesian coordinates.

        Parameters
        ----------
        frame
            :class:`~shakelab.libutils.geodeticN.transform.MetricFrame`
            of the output, or ``None`` for absolute ECEF.

        Returns
        -------
        MetricPoint
            Cartesian coordinates in meters.
        """
        from shakelab.libutils.geodeticN.transform import wgs_to_metric

        return wgs_to_metric(self, frame=frame)

    def epicentral_distance_to(self, other: WgsPoint) -> float:
        """
        Return the great-circle (epicentral) distance in meters.

        See :func:`~shakelab.libutils.geodeticN.distance.epicentral_distance`.
        """
        from shakelab.libutils.geodeticN.distance import (
            epicentral_distance,
        )

        return epicentral_distance(self, other)

    def hypocentral_distance_to(
        self,
        other: WgsPoint,
        approx: str = "ellipsoid",
    ) -> float:
        """
        Return the 3-D (hypocentral) distance in meters.

        Parameters
        ----------
        other
            Target point.
        approx
            Earth model: ``"ellipsoid"`` (default) or ``"sphere"``.

        See :func:`~shakelab.libutils.geodeticN.distance.hypocentral_distance`.
        """
        from shakelab.libutils.geodeticN.distance import (
            hypocentral_distance,
        )

        return hypocentral_distance(self, other, approx=approx)

    def __sub__(self, other: WgsPoint) -> float:
        """
        Return the hypocentral distance to ``other`` (ellipsoid model).
        """
        if not isinstance(other, WgsPoint):
            return NotImplemented

        return self.hypocentral_distance_to(other)


# ---------------------------------------------------------------------------
# Polygon
# ---------------------------------------------------------------------------

class PolygonKind(IntEnum):
    """
    Polygon type: a single polygon or a set of disjoint polygons.
    """

    SIMPLE = 0
    MULTIPOLYGON = 1


@dataclass
class WgsPolygon:
    """
    Polygon or multipolygon in geographic coordinates.

    The geometry is a list of parts; each part is the exterior ring of a
    simple polygon. Holes are not represented.

    Parameters
    ----------
    parts
        List of parts, each a list of :class:`WgsPoint`.
    kind
        :attr:`PolygonKind.SIMPLE` or :attr:`PolygonKind.MULTIPOLYGON`.
        If omitted, it is inferred from the number of parts.

    Raises
    ------
    TypeError
        If a vertex is not a :class:`WgsPoint`.
    ValueError
        If ``kind`` is SIMPLE but there is more than one part.
    """

    parts: List[List[WgsPoint]] = field(default_factory=list)
    kind: Optional[PolygonKind] = None

    def __post_init__(self) -> None:
        self.parts = [list(part) for part in self.parts]

        for part in self.parts:
            for vertex in part:
                if not isinstance(vertex, WgsPoint):
                    raise TypeError("Polygon vertices must be WgsPoint.")

        if self.kind is None:
            self.kind = (
                PolygonKind.MULTIPOLYGON
                if len(self.parts) > 1
                else PolygonKind.SIMPLE
            )
        else:
            self.kind = PolygonKind(self.kind)

        if self.kind == PolygonKind.SIMPLE and len(self.parts) > 1:
            raise ValueError("A SIMPLE polygon cannot have several parts.")

    # -----------------------------------------------------------------------
    # Constructors
    # -----------------------------------------------------------------------

    @classmethod
    def from_lonlat(cls, coords: Sequence[Sequence[float]]) -> WgsPolygon:
        """
        Build a simple polygon from ``(lon, lat[, elevation])`` vertices.
        """
        return cls(parts=[_points_from_coords(coords)])

    @classmethod
    def from_parts(
        cls,
        parts: Sequence[Sequence[Sequence[float]]],
        kind: Optional[PolygonKind] = None,
    ) -> WgsPolygon:
        """
        Build a polygon or multipolygon from several parts.

        Parameters
        ----------
        parts
            Each element is a sequence of ``(lon, lat[, elevation])``
            vertices.
        kind
            Polygon kind; inferred from the number of parts if omitted.
        """
        return cls(
            parts=[_points_from_coords(part) for part in parts],
            kind=kind,
        )

    # -----------------------------------------------------------------------
    # Basic properties
    # -----------------------------------------------------------------------

    def is_empty(self) -> bool:
        """
        Return True if the polygon has no parts.
        """
        return len(self.parts) == 0

    def n_parts(self) -> int:
        """
        Return the number of parts.
        """
        return len(self.parts)

    def iter_vertices(self) -> Iterator[WgsPoint]:
        """
        Iterate over the vertices of all parts, in order.
        """
        for part in self.parts:
            yield from part

    def to_arrays(self) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Return the ``(lon, lat)`` arrays of each part.
        """
        return [_lonlat_arrays(part) for part in self.parts]

    # -----------------------------------------------------------------------
    # Geometry
    # -----------------------------------------------------------------------

    def bounding_box(self) -> Tuple[float, float, float, float]:
        """
        Return ``(lon_min, lon_max, lat_min, lat_max)`` in degrees.

        Longitude bounds are the plain minimum and maximum of the vertex
        longitudes, so they are not meaningful for polygons crossing the
        antimeridian in the ``[-180, 180)`` convention.

        Raises
        ------
        ValueError
            If the polygon is empty.
        """
        if self.is_empty():
            raise ValueError("Cannot compute bbox of empty polygon.")

        return _bounding_box(list(self.iter_vertices()))

    def centroid(self) -> WgsPoint:
        """
        Return the area-weighted centroid of the polygon.

        The centroid is computed in the sinusoidal projection centered
        on the mean longitude of the vertices, and returned in the
        longitude convention of the first vertex. For degenerate
        polygons (zero area) the mean of the vertices is returned.

        Raises
        ------
        ValueError
            If the polygon is empty.
        """
        if self.is_empty():
            raise ValueError("Cannot compute centroid of empty polygon.")

        vertices = [
            vertex for part in self.parts for vertex in _open_ring(part)
        ]
        lon0 = _reference_longitude(vertices)
        total_area = 0.0
        moment_x = 0.0
        moment_y = 0.0

        for part in self.parts:
            ring = _open_ring(part)

            if len(ring) < 3:
                continue

            lon, lat = _lonlat_arrays(ring)
            x, y = wgs84_to_sinusoidal(lon, lat, lon0_deg=lon0)
            area, cx, cy = _planar_centroid(x, y)

            total_area += area
            moment_x += area * cx
            moment_y += area * cy

        if total_area == 0.0:
            return WgsPoint(
                longitude=lon0,
                latitude=float(np.mean([p.latitude for p in vertices])),
            )

        lon, lat = sinusoidal_to_wgs84(
            moment_x / total_area,
            moment_y / total_area,
            lon0_deg=lon0,
        )

        return WgsPoint(longitude=float(lon), latitude=float(lat))

    def area(self) -> float:
        """
        Return the polygon area in square meters.

        Each part is projected with the sinusoidal (equal-area)
        projection centered on its mean longitude, and the planar areas
        are summed.
        """
        total = 0.0

        for part in self.parts:
            ring = _open_ring(part)

            if len(ring) < 3:
                continue

            lon, lat = _lonlat_arrays(ring)
            x, y = wgs84_to_sinusoidal(
                lon, lat, lon0_deg=_reference_longitude(ring),
            )
            total += abs(_signed_area(x, y))

        return float(total)

    def contains_point(self, point: WgsPoint) -> bool:
        """
        Return True if ``point`` lies inside the polygon.

        For a multipolygon, the point is inside if it lies in any part.
        """
        inside = self.contains_points([point.longitude], [point.latitude])

        return bool(inside[0])

    def contains_points(
        self,
        longitudes: Sequence[float],
        latitudes: Sequence[float],
    ) -> np.ndarray:
        """
        Test many points at once.

        Parameters
        ----------
        longitudes, latitudes
            Coordinates of the points in degrees.

        Returns
        -------
        numpy.ndarray
            Boolean array, True for points inside the polygon.
        """
        lon = np.asarray(longitudes, dtype=float)
        lat = np.asarray(latitudes, dtype=float)
        inside = np.zeros(lon.shape, dtype=bool)

        for part in self.parts:
            ring = _open_ring(part)

            if len(ring) < 3:
                continue

            lon0 = _reference_longitude(ring)
            ring_lon, ring_lat = _lonlat_arrays(ring)
            x, y = wgs84_to_sinusoidal(ring_lon, ring_lat, lon0_deg=lon0)
            px, py = wgs84_to_sinusoidal(lon, lat, lon0_deg=lon0)

            inside |= _points_in_ring(x, y, px, py)

        return inside


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

@dataclass
class WgsMesh:
    """
    Unstructured collection of geographic points.

    Parameters
    ----------
    points
        List of :class:`WgsPoint` nodes.
    """

    points: List[WgsPoint] = field(default_factory=list)

    @classmethod
    def from_lonlat(cls, coords: Sequence[Sequence[float]]) -> WgsMesh:
        """
        Build a mesh from ``(lon, lat[, elevation])`` tuples.
        """
        return cls(points=_points_from_coords(coords))

    @classmethod
    def from_arrays(
        cls,
        longitudes: Sequence[float],
        latitudes: Sequence[float],
        elevations: Optional[Sequence[float]] = None,
    ) -> WgsMesh:
        """
        Build a mesh from coordinate arrays.

        Parameters
        ----------
        longitudes, latitudes
            Coordinates in degrees.
        elevations
            Elevations in meters; zero if omitted.

        Raises
        ------
        ValueError
            If the arrays have different lengths.
        """
        lon = np.asarray(longitudes, dtype=float).ravel()
        lat = np.asarray(latitudes, dtype=float).ravel()

        if elevations is None:
            elev = np.zeros_like(lon)
        else:
            elev = np.asarray(elevations, dtype=float).ravel()

        if not lon.size == lat.size == elev.size:
            raise ValueError("Coordinate arrays must have the same length.")

        return cls(points=[
            WgsPoint(longitude=x, latitude=y, elevation=z)
            for x, y, z in zip(lon, lat, elev)
        ])

    def __len__(self) -> int:
        """
        Return the number of points.
        """
        return len(self.points)

    def __iter__(self) -> Iterator[WgsPoint]:
        """
        Iterate over the points.
        """
        return iter(self.points)

    def __getitem__(self, index: int) -> WgsPoint:
        """
        Return the point at ``index``.
        """
        return self.points[index]

    def to_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return ``(longitudes, latitudes, elevations)`` as arrays.
        """
        lon, lat = _lonlat_arrays(self.points)
        elev = np.array([p.elevation for p in self.points], dtype=float)

        return lon, lat, elev

    def bounding_box(self) -> Tuple[float, float, float, float]:
        """
        Return ``(lon_min, lon_max, lat_min, lat_max)`` in degrees.

        Raises
        ------
        ValueError
            If the mesh is empty.
        """
        if not self.points:
            raise ValueError("Cannot compute bbox of empty mesh.")

        return _bounding_box(self.points)

    def clip(self, polygon: WgsPolygon) -> WgsMesh:
        """
        Return a new mesh with the points inside ``polygon``.

        Parameters
        ----------
        polygon
            Clipping polygon or multipolygon.

        Returns
        -------
        WgsMesh
            Points of this mesh inside the polygon, in the same order.
        """
        if not self.points:
            return WgsMesh()

        lon, lat, _ = self.to_arrays()
        inside = polygon.contains_points(lon, lat)

        return WgsMesh(points=[
            point for point, keep in zip(self.points, inside) if keep
        ])


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _finite_float(value: Any, name: str) -> float:
    """
    Convert a real number to a finite float.
    """
    if isinstance(value, bool) or not isinstance(value, (Real, Decimal)):
        raise TypeError(f"The {name} must be a real number.")

    result = float(value)

    if not math.isfinite(result):
        raise ValueError(f"The {name} must be finite.")

    return result


def _points_from_coords(coords: Sequence[Sequence[float]]) -> List[WgsPoint]:
    """
    Convert ``(lon, lat[, elevation])`` sequences to points.
    """
    points = []

    for coord in coords:
        if len(coord) not in (2, 3):
            raise ValueError(
                "Coordinates must be (lon, lat) or (lon, lat, elevation)."
            )

        points.append(WgsPoint(*coord))

    return points


def _lonlat_arrays(points: Sequence[WgsPoint]) -> Tuple[np.ndarray, ...]:
    """
    Return the longitude and latitude arrays of a list of points.
    """
    lon = np.array([p.longitude for p in points], dtype=float)
    lat = np.array([p.latitude for p in points], dtype=float)

    return lon, lat


def _bounding_box(
    points: Sequence[WgsPoint],
) -> Tuple[float, float, float, float]:
    """
    Return the plain bounding box of a list of points.
    """
    lon, lat = _lonlat_arrays(points)

    return (
        float(lon.min()),
        float(lon.max()),
        float(lat.min()),
        float(lat.max()),
    )


def _open_ring(part: Sequence[WgsPoint]) -> List[WgsPoint]:
    """
    Drop the closing vertex of a ring, if it repeats the first one.
    """
    ring = list(part)

    if len(ring) > 1 and ring[0] == ring[-1]:
        ring = ring[:-1]

    return ring


def _reference_longitude(points: Sequence[WgsPoint]) -> float:
    """
    Return the mean longitude of points, relative to the first one.

    The result is in the convention of the first point and is not
    affected by the antimeridian, for points spanning less than 180
    degrees of longitude.
    """
    lon0 = points[0].longitude
    dlon = relative_longitude([p.longitude for p in points], lon0)

    return lon0 + float(np.mean(dlon))


def _signed_area(x: np.ndarray, y: np.ndarray) -> float:
    """
    Return the signed planar area of a ring (shoelace formula).
    """
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _planar_centroid(
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[float, float, float]:
    """
    Return ``(area, cx, cy)`` of a planar ring, area non-negative.
    """
    x_next = np.roll(x, -1)
    y_next = np.roll(y, -1)
    cross = x * y_next - x_next * y
    signed = 0.5 * float(np.sum(cross))

    if signed == 0.0:
        return 0.0, float(np.mean(x)), float(np.mean(y))

    cx = float(np.sum((x + x_next) * cross)) / (6.0 * signed)
    cy = float(np.sum((y + y_next) * cross)) / (6.0 * signed)

    return abs(signed), cx, cy


def _points_in_ring(
    x: np.ndarray,
    y: np.ndarray,
    px: np.ndarray,
    py: np.ndarray,
) -> np.ndarray:
    """
    Vectorized ray-casting test of points against a planar ring.

    Parameters
    ----------
    x, y
        Ring vertices (open or closed).
    px, py
        Test points.

    Returns
    -------
    numpy.ndarray
        Boolean array, True for points inside the ring.
    """
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    inside = np.zeros(px.shape, dtype=bool)

    j = len(x) - 1

    for i in range(len(x)):
        xi, yi = x[i], y[i]
        xj, yj = x[j], y[j]

        crosses = (yi > py) != (yj > py)

        with np.errstate(divide="ignore", invalid="ignore"):
            x_cross = (xj - xi) * (py - yi) / (yj - yi) + xi

        inside ^= crosses & (px < x_cross)
        j = i

    return inside


def _ray_casting_contains(
    x: np.ndarray,
    y: np.ndarray,
    x0: float,
    y0: float,
) -> bool:
    """
    Ray-casting test of a single point against a planar ring.
    """
    return bool(_points_in_ring(x, y, np.array([x0]), np.array([y0]))[0])
