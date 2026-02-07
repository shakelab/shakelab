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
Spatial sampling utilities for ShakeLab.

This module provides functions to draw random samples from geometric
objects. The focus is on uniform sampling within polygons, either in
planar coordinates or in geographic (WGS84) coordinates.

Conventions
----------
* Geographic coordinates are expressed as (longitude, latitude) in
  decimal degrees.
* Distances are expressed in meters.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np
from scipy.spatial import Delaunay

from shakelab.libutils.geodeticN.utilities import (
    wgs84_to_sinusoidal,
    sinusoidal_to_wgs84,
)
from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
    _ray_casting_contains,
)


# ============================================================================
# Low-level helpers
# ============================================================================


def _sample_point_in_triangle(
    v1: np.ndarray,
    v2: np.ndarray,
    v3: np.ndarray,
) -> np.ndarray:
    """
    Uniformly sample a point inside a triangle.

    Parameters
    ----------
    v1, v2, v3 : numpy.ndarray
        Vertices of the triangle as 1D arrays of shape (2,).

    Returns
    -------
    numpy.ndarray
        Sampled point of shape (2,).
    """
    r1 = np.sqrt(np.random.random())
    r2 = np.random.random()
    return (1.0 - r1) * v1 + (r1 * (1.0 - r2)) * v2 + (r1 * r2) * v3


def _triangle_areas(triangles: np.ndarray) -> np.ndarray:
    """
    Compute areas of triangles in the plane.

    Parameters
    ----------
    triangles : numpy.ndarray
        Array of shape (n_tri, 3, 2) containing triangle vertices.

    Returns
    -------
    numpy.ndarray
        Areas of shape (n_tri,). Values are always non-negative.
    """
    v1 = triangles[:, 0, :]
    v2 = triangles[:, 1, :]
    v3 = triangles[:, 2, :]
    area = 0.5 * np.abs(np.cross(v2 - v1, v3 - v1))
    return area


# ============================================================================
# Sampling in planar polygons
# ============================================================================


def sample_uniform_in_polygon_xy(
    vertices_xy: np.ndarray,
    n_samples: int,
) -> np.ndarray:
    """
    Uniformly sample points inside a planar polygon.

    The polygon is defined by its vertices in 2D Cartesian coordinates.
    The algorithm uses a Delaunay triangulation of the vertices as a
    proposal distribution and rejection sampling to ensure that all
    points lie inside the polygon.

    Parameters
    ----------
    vertices_xy : numpy.ndarray
        Polygon vertices of shape (n_vertices, 2). The polygon is
        assumed to be simple and non self-intersecting.
    n_samples : int
        Number of random points to generate.

    Returns
    -------
    numpy.ndarray
        Sampled points of shape (n_samples, 2).
    """
    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")

    vertices_xy = np.asarray(vertices_xy, dtype=float)
    if vertices_xy.ndim != 2 or vertices_xy.shape[1] != 2:
        raise ValueError("vertices_xy must have shape (n, 2).")

    if vertices_xy.shape[0] < 3:
        raise ValueError("Polygon must have at least three vertices.")

    delaunay = Delaunay(vertices_xy)
    triangles = vertices_xy[delaunay.simplices]
    areas = _triangle_areas(triangles)

    if not np.all(areas > 0.0):
        raise ValueError("Degenerate triangles found in triangulation.")

    cum_areas = np.cumsum(areas)
    total_area = cum_areas[-1]

    samples: List[np.ndarray] = []

    while len(samples) < n_samples:
        rnd = np.random.random() * total_area
        idx = int(np.searchsorted(cum_areas, rnd))
        tri = triangles[idx]

        p = _sample_point_in_triangle(tri[0], tri[1], tri[2])

        if _ray_casting_contains(
            vertices_xy[:, 0],
            vertices_xy[:, 1],
            float(p[0]),
            float(p[1]),
        ):
            samples.append(p)

    return np.vstack(samples)


# ============================================================================
# Sampling in WGS84 polygons
# ============================================================================


def sample_uniform_in_polygon_wgs84(
    polygon: WgsPolygon,
    n_samples: int,
) -> WgsMesh:
    """
    Uniformly sample points inside a WGS84 polygon or multipolygon.

    The polygon is projected to an equal-area sinusoidal projection.
    Sampling is performed in the projected plane using
    :func:`sample_uniform_in_polygon_xy`. Samples are then mapped
    back to geographic coordinates.

    Parameters
    ----------
    polygon : WgsPolygon
        Target polygon or multipolygon.
    n_samples : int
        Number of random points to generate.

    Returns
    -------
    WgsMesh
        Mesh containing the sampled WgsPoint nodes.
    """
    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")

    if polygon.is_empty():
        return WgsMesh(points=[])

    # Project each part to sinusoidal and compute area.
    part_xy: List[Tuple[np.ndarray, np.ndarray]] = []
    part_area: List[float] = []

    for part in polygon.parts:
        if len(part) < 3:
            continue

        lon = np.array([p.longitude for p in part])
        lat = np.array([p.latitude for p in part])
        x, y = wgs84_to_sinusoidal(lon, lat)

        part_xy.append((x, y))

        x_shift = np.roll(x, -1)
        y_shift = np.roll(y, -1)
        area = 0.5 * np.sum(x * y_shift - x_shift * y)
        part_area.append(abs(float(area)))

    if not part_area:
        return WgsMesh(points=[])

    area_arr = np.array(part_area)
    area_sum = float(area_arr.sum())
    if area_sum <= 0.0:
        return WgsMesh(points=[])

    # Allocate samples to each part using a multinomial draw.
    probs = area_arr / area_sum
    counts = np.random.multinomial(n_samples, probs)

    sampled_lon: List[float] = []
    sampled_lat: List[float] = []

    for (x, y), count in zip(part_xy, counts):
        if count <= 0:
            continue

        verts = np.column_stack((x, y))
        pts_xy = sample_uniform_in_polygon_xy(verts, count)

        lon_deg, lat_deg = sinusoidal_to_wgs84(
            pts_xy[:, 0],
            pts_xy[:, 1],
        )

        sampled_lon.extend(float(v) for v in lon_deg)
        sampled_lat.extend(float(v) for v in lat_deg)

    coords = list(zip(sampled_lon, sampled_lat))
    return WgsMesh.from_lonlat(coords)
