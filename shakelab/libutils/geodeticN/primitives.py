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
Geometric primitives for geodetic operations in ShakeLab.

This module defines WGS84-based geometric objects used across the
library. It provides light-weight containers for points, polygons,
multipolygons and unstructured meshes, together with a few basic
operations (distance, area, bounding boxes).

Conventions
----------
* Geographic coordinates are expressed as (longitude, latitude) in
  decimal degrees, unless explicitly stated otherwise.
* Distances are expressed in meters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Iterable, List, Sequence, Tuple

import numpy as np

from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
    wgs84_to_sinusoidal,
)


# ============================================================================
# Basic types
# ============================================================================


@dataclass(frozen=True)
class WgsPoint:
    """
    Simple WGS84 point.

    Parameters
    ----------
    longitude : float
        Longitude in degrees.
    latitude : float
        Latitude in degrees.
    elevation : float, optional
        Elevation above the reference surface in meters. The default
        is 0.0. The elevation is not used in most geometric
        operations, which are based on longitude and latitude only.
    """

    longitude: float
    latitude: float
    elevation: float = 0.0

    def to_tuple(self) -> Tuple[float, float]:
        """
        Return the point as a (lon, lat) tuple.
        """
        return self.longitude, self.latitude

    def to_metric(self, frame=None) -> "MetricPoint":
        """
        Convert this WgsPoint to a MetricPoint.

        Parameters
        ----------
        frame : MetricFrame or None, optional
            Output frame definition.

            If None, the returned coordinates are expressed in absolute
            ECEF (Earth-Centered, Earth-Fixed).

            If a MetricFrame is provided, the returned coordinates are
            expressed in the local Cartesian system defined by the frame
            (translated ECEF or ENU, depending on frame.orientation).

        Returns
        -------
        MetricPoint
            Cartesian coordinates in meters.
        """
        # Local import to avoid circular dependencies (transform imports
        # WgsPoint from this module).
        from shakelab.libutils.geodeticN.transform import wgs_to_metric

        return wgs_to_metric(self, frame=frame)

    def epicentral_distance_to(self, other: "WgsPoint") -> float:
        """
        Compute great-circle (epicentral) distance to another WgsPoint.

        This is a convenience wrapper around the function
        :func:`epicentral_distance` defined in the module
        ``distance.py``.

        Parameters
        ----------
        other : WgsPoint
            Target point.

        Returns
        -------
        float
            Great-circle distance in meters.
        """
        from shakelab.libutils.geodeticN.distance import (
            epicentral_distance,
        )
        return epicentral_distance(self, other)

    def hypocentral_distance_to(
        self,
        other: "WgsPoint",
        approx: str = "ellipsoid",
    ) -> float:
        """
        Compute hypocentral (3-D straight-line) distance to another WgsPoint.

        This is a wrapper around :func:`hypocentral_distance` defined
        in ``distance.py``. The distance is evaluated either using the
        accurate ECEF/WGS84 ellipsoid or a spherical Earth model
        depending on the ``approx`` parameter.

        Parameters
        ----------
        other : WgsPoint
            Target point.
        approx : {"ellipsoid", "sphere"}, optional
            Approximation model. Default is "ellipsoid".

        Returns
        -------
        float
            Hypocentral (tunnel) distance in meters.

        Raises
        ------
        ValueError
            If ``approx`` is not one of the supported options.
        """
        from shakelab.libutils.geodeticN.distance import (
            hypocentral_distance,
        )
        return hypocentral_distance(self, other, approx=approx)

    def __sub__(self, other: "WgsPoint") -> float:
        """
        Return hypocentral (3-D tunnel) distance to another WgsPoint.

        This operator overload is provided as syntactic sugar so that
        the expression ``p1 - p2`` returns the hypocentral distance
        between two geographic points. The calculation uses the same
        approximation model adopted by :meth:`hypocentral_distance_to`,
        which defaults to the WGS84 ellipsoid.

        Parameters
        ----------
        other : WgsPoint
            Target point.

        Returns
        -------
        float
            Hypocentral (straight-line) distance in meters.

        Raises
        ------
        ValueError
            If the approximation model is not recognized.
        """
        if not isinstance(other, WgsPoint):
            return NotImplemented
        return self.hypocentral_distance_to(other)


class PolygonKind(IntEnum):
    """
    Polygon type indicator.

    Attributes
    ----------
    SIMPLE : int
        Single simple polygon (one exterior ring).
    MULTIPOLYGON : int
        MultiPolygon: list of disjoint simple polygons.
    """

    SIMPLE = 0
    MULTIPOLYGON = 1


@dataclass
class WgsPolygon:
    """
    WGS84 polygon or multipolygon.

    The geometry is stored as a list of parts, where each part is a
    sequence of vertices forming a simple polygon (exterior ring).
    Holes are not explicitly represented at the moment. A MultiPolygon
    is represented as multiple parts.

    Parameters
    ----------
    parts : list of list of WgsPoint
        List of polygon parts. For a simple polygon the list contains
        a single part.
    kind : PolygonKind, optional
        Polygon type: :class:`PolygonKind.SIMPLE` or
        :class:`PolygonKind.MULTIPOLYGON`. If not provided, the
        kind is inferred from the number of parts.
    """

    parts: List[List[WgsPoint]] = field(default_factory=list)
    kind: PolygonKind = PolygonKind.SIMPLE

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_lonlat(
        cls,
        coords: Sequence[Tuple[float, float]],
    ) -> "WgsPolygon":
        """
        Build a simple polygon from a sequence of (lon, lat) pairs.

        Parameters
        ----------
        coords : sequence of (float, float)
            Sequence of (lon, lat) vertices.

        Returns
        -------
        WgsPolygon
            Polygon with a single part.
        """
        part = [WgsPoint(longitude=c[0], latitude=c[1]) for c in coords]
        return cls(parts=[part], kind=PolygonKind.SIMPLE)

    @classmethod
    def from_parts(
        cls,
        parts: Sequence[Sequence[Tuple[float, float]]],
    ) -> "WgsPolygon":
        """
        Build a polygon or multipolygon from multiple parts.

        Parameters
        ----------
        parts : sequence of sequences
            Each element is a sequence of (lon, lat) pairs forming a
            simple polygon.

        Returns
        -------
        WgsPolygon
            Polygon with one or more parts. The kind is set to
            SIMPLE if there is a single part, MULTIPOLYGON otherwise.
        """
        poly_parts: List[List[WgsPoint]] = []
        for seq in parts:
            poly_parts.append(
                [WgsPoint(longitude=c[0], latitude=c[1]) for c in seq]
            )

        kind = (
            PolygonKind.SIMPLE
            if len(poly_parts) == 1
            else PolygonKind.MULTIPOLYGON
        )
        return cls(parts=poly_parts, kind=kind)

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------

    def is_empty(self) -> bool:
        """
        Return True if the polygon has no parts.
        """
        return len(self.parts) == 0

    def n_parts(self) -> int:
        """
        Return the number of polygon parts.
        """
        return len(self.parts)

    def iter_vertices(self) -> Iterable[WgsPoint]:
        """
        Iterate over all vertices of all parts.

        Yields
        ------
        WgsPoint
            Vertices in the order they appear in the parts.
        """
        for part in self.parts:
            for pt in part:
                yield pt

    # ------------------------------------------------------------------
    # Bounding box and centroid
    # ------------------------------------------------------------------

    def bounding_box(self) -> Tuple[float, float, float, float]:
        """
        Compute polygon bounding box.

        Returns
        -------
        lon_min, lon_max, lat_min, lat_max : float
            Bounding box in degrees.
        """
        if self.is_empty():
            raise ValueError("Cannot compute bbox of empty polygon.")

        lons = [p.longitude for p in self.iter_vertices()]
        lats = [p.latitude for p in self.iter_vertices()]

        lon_min = float(min(lons))
        lon_max = float(max(lons))
        lat_min = float(min(lats))
        lat_max = float(max(lats))

        return lon_min, lon_max, lat_min, lat_max

    def centroid(self) -> WgsPoint:
        """
        Compute an approximate centroid of the polygon.

        The centroid is computed as the arithmetic mean of vertex
        coordinates. This is a simple approximation and not the
        exact centroid in spherical geometry.

        Returns
        -------
        WgsPoint
            Approximate centroid.
        """
        if self.is_empty():
            raise ValueError("Cannot compute centroid of empty polygon.")

        lons = [p.longitude for p in self.iter_vertices()]
        lats = [p.latitude for p in self.iter_vertices()]

        lon_c = float(np.mean(lons))
        lat_c = float(np.mean(lats))

        return WgsPoint(longitude=lon_c, latitude=lat_c)

    # ------------------------------------------------------------------
    # Area and point-in-polygon
    # ------------------------------------------------------------------

    def area(self) -> float:
        """
        Compute polygon area in square meters.

        The area is computed by projecting each part onto an equal-area
        sinusoidal projection and then summing the planar polygon areas.

        Returns
        -------
        float
            Total area in square meters.
        """
        if self.is_empty():
            return 0.0

        total_area = 0.0

        for part in self.parts:
            if len(part) < 3:
                continue

            lon = np.array([p.longitude for p in part])
            lat = np.array([p.latitude for p in part])
            x, y = wgs84_to_sinusoidal(lon, lat)

            # Shoelace formula for polygon area in the projected plane.
            x_shift = np.roll(x, -1)
            y_shift = np.roll(y, -1)
            area_part = 0.5 * np.sum(x * y_shift - x_shift * y)

            total_area += abs(area_part)

        return float(total_area)

    def contains_point(self, point: WgsPoint) -> bool:
        """
        Test whether a point lies inside the polygon.

        The test is performed in the sinusoidal projection using a
        standard ray-casting algorithm. For a multipolygon, the point
        is considered inside if it lies in at least one part.

        Parameters
        ----------
        point : WgsPoint
            Test point.

        Returns
        -------
        bool
            True if the point is inside the polygon, False otherwise.
        """
        if self.is_empty():
            return False

        # Project the point once; reuse in each part.
        x_pt, y_pt = wgs84_to_sinusoidal(
            np.array([point.longitude]),
            np.array([point.latitude]),
        )
        x0 = float(x_pt[0])
        y0 = float(y_pt[0])

        inside_any = False

        for part in self.parts:
            if len(part) < 3:
                continue

            lon = np.array([p.longitude for p in part])
            lat = np.array([p.latitude for p in part])
            x, y = wgs84_to_sinusoidal(lon, lat)

            inside = _ray_casting_contains(x, y, x0, y0)
            if inside:
                inside_any = True
                break

        return inside_any


def _ray_casting_contains(
    x: np.ndarray,
    y: np.ndarray,
    x0: float,
    y0: float,
) -> bool:
    """
    Ray-casting point-in-polygon test on planar coordinates.

    Parameters
    ----------
    x, y : numpy.ndarray
        Polygon vertices in the plane.
    x0, y0 : float
        Test point in the same plane.

    Returns
    -------
    bool
        True if the point is inside the polygon, False otherwise.
    """
    n = len(x)
    inside = False

    j = n - 1
    for i in range(n):
        xi = x[i]
        yi = y[i]
        xj = x[j]
        yj = y[j]

        intersects = (
            (yi > y0) != (yj > y0)
            and x0
            < (xj - xi) * (y0 - yi) / (yj - yi + 1e-15) + xi
        )
        if intersects:
            inside = not inside
        j = i

    return inside


# ============================================================================
# Mesh of points
# ============================================================================


@dataclass
class WgsMesh:
    """
    Collection of WGS84 points.

    This class represents an unstructured mesh of points in geographic
    coordinates. It is primarily used as a container for grids and
    sampling locations.

    Parameters
    ----------
    points : list of WgsPoint
        List of mesh nodes.
    """

    points: List[WgsPoint] = field(default_factory=list)

    @classmethod
    def from_lonlat(
        cls,
        coords: Sequence[Tuple[float, float]],
    ) -> "WgsMesh":
        """
        Build a mesh from a list of (lon, lat) pairs.

        Parameters
        ----------
        coords : sequence of (float, float)
            Vertex coordinates.

        Returns
        -------
        WgsMesh
            Mesh instance.
        """
        pts = [WgsPoint(longitude=c[0], latitude=c[1]) for c in coords]
        return cls(points=pts)

    def __len__(self) -> int:
        """
        Return the number of points in the mesh.
        """
        return len(self.points)

    def __iter__(self):
        """
        Iterate over mesh points.
        """
        return iter(self.points)

    def bounding_box(self) -> Tuple[float, float, float, float]:
        """
        Compute mesh bounding box.

        Returns
        -------
        lon_min, lon_max, lat_min, lat_max : float
            Bounding box in degrees.

        Raises
        ------
        ValueError
            If the mesh is empty.
        """
        if not self.points:
            raise ValueError("Cannot compute bbox of empty mesh.")

        lons = [p.longitude for p in self.points]
        lats = [p.latitude for p in self.points]

        lon_min = float(min(lons))
        lon_max = float(max(lons))
        lat_min = float(min(lats))
        lat_max = float(max(lats))

        return lon_min, lon_max, lat_min, lat_max
