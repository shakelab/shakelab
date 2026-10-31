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
Generation of point meshes in geographic coordinates.

The module provides generators of :class:`WgsMesh` nodes:

- :func:`mesh_fibonacci_global`: quasi equal-area Fibonacci mesh on the
  whole sphere;
- :func:`mesh_fibonacci_spherical_cap`: Fibonacci mesh on a spherical
  cap;
- :func:`mesh_fibonacci_bbox`: Fibonacci cap mesh clipped to a box;
- :func:`mesh_cartesian`: regular longitude-latitude grid;
- :func:`mesh_sinusoidal`: equal-area grid in the sinusoidal projection,
  aligned on a global grid centered at (0°, 0°).

A mesh restricted to a polygon is obtained by clipping, e.g.
``mesh_sinusoidal(...).clip(polygon)``.

Design choices
--------------
- Spacings are expressed in kilometers, as the typical distance between
  neighbouring nodes.
- Generators restricted to a longitude range compute only the nodes in
  that range; the global sinusoidal grid alignment does not depend on
  the requested window.

Conventions
-----------
- Coordinates are ``(lon, lat)`` in degrees.
- A longitude range ``[lon_min, lon_max]`` with ``lon_max < lon_min``
  crosses the antimeridian (e.g. ``170`` to ``-170``). Longitudes of
  nodes generated in a range are returned in
  ``[lon_min, lon_min + span]``, i.e. in the convention of ``lon_min``.
"""

from __future__ import annotations

import math
from typing import List, Optional, Tuple

import numpy as np

from shakelab.libutils.geodeticN.primitives import WgsMesh, WgsPoint
from shakelab.libutils.geodeticN.utilities import (
    DEG2RAD,
    MEAN_EARTH_RADIUS_M,
    RAD2DEG,
    great_circle_distance,
    normalize_longitude,
)


__all__ = [
    "mesh_fibonacci_global",
    "mesh_fibonacci_spherical_cap",
    "mesh_fibonacci_bbox",
    "mesh_cartesian",
    "mesh_sinusoidal",
]


GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))

# Tolerance on grid indices, in units of grid steps.
_INDEX_TOLERANCE = 1e-9


# ---------------------------------------------------------------------------
# Fibonacci meshes
# ---------------------------------------------------------------------------

def mesh_fibonacci_global(
    spacing_km: float,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a global Fibonacci (golden-spiral) mesh.

    Nodes are distributed quasi-uniformly in area over the whole sphere.
    For regional grids use :func:`mesh_fibonacci_spherical_cap`.

    Parameters
    ----------
    spacing_km
        Target spacing in kilometers; the mean area per node is about
        ``spacing_km ** 2``.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Mesh nodes, with longitudes in ``[-180, 180)``.

    Raises
    ------
    ValueError
        If ``spacing_km`` is not positive.
    """
    _require_positive(spacing_km, "spacing_km")

    area_node = (spacing_km * 1e3) ** 2
    area_sphere = 4.0 * math.pi * radius_m ** 2
    n_points = max(1, int(round(area_sphere / area_node)))

    index = np.arange(n_points, dtype=float)
    lat = np.arcsin(1.0 - 2.0 * (index + 0.5) / n_points) * RAD2DEG
    lon = normalize_longitude(GOLDEN_ANGLE * index * RAD2DEG)

    return WgsMesh.from_arrays(np.atleast_1d(lon), lat)


def mesh_fibonacci_spherical_cap(
    center_lon_deg: float,
    center_lat_deg: float,
    max_distance_km: float,
    spacing_km: float,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a quasi equal-area Fibonacci mesh on a spherical cap.

    Parameters
    ----------
    center_lon_deg, center_lat_deg
        Cap center in degrees.
    max_distance_km
        Cap radius (great-circle distance) in kilometers.
    spacing_km
        Target spacing in kilometers.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Mesh nodes. Longitudes are ``center_lon_deg`` plus an offset in
        ``[-180, 180]``, i.e. in the convention of the center.

    Raises
    ------
    ValueError
        If ``spacing_km`` or ``max_distance_km`` is not positive.
    """
    _require_positive(spacing_km, "spacing_km")
    _require_positive(max_distance_km, "max_distance_km")

    alpha = min(max_distance_km * 1e3 / radius_m, math.pi)

    area_cap = 2.0 * math.pi * radius_m ** 2 * (1.0 - math.cos(alpha))
    area_node = (spacing_km * 1e3) ** 2
    n_points = max(1, int(round(area_cap / area_node)))

    lat0 = center_lat_deg * DEG2RAD
    sin_lat0 = math.sin(lat0)
    cos_lat0 = math.cos(lat0)

    index = np.arange(n_points, dtype=float)
    z = 1.0 - (index + 0.5) / n_points * (1.0 - math.cos(alpha))
    theta = np.arccos(z)
    phi = GOLDEN_ANGLE * index

    sin_lat = sin_lat0 * np.cos(theta) + cos_lat0 * np.sin(theta) * np.cos(phi)
    lat = np.arcsin(np.clip(sin_lat, -1.0, 1.0))

    dlon = np.arctan2(
        np.sin(phi) * np.sin(theta) * cos_lat0,
        np.cos(theta) - sin_lat0 * sin_lat,
    )

    return WgsMesh.from_arrays(
        center_lon_deg + dlon * RAD2DEG,
        lat * RAD2DEG,
    )


def mesh_fibonacci_bbox(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    spacing_km: float,
    center_lon_deg: Optional[float] = None,
    center_lat_deg: Optional[float] = None,
    max_distance_km: Optional[float] = None,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a Fibonacci cap mesh clipped to a bounding box.

    Parameters
    ----------
    lon_min, lon_max
        Longitude range in degrees (see the module conventions for
        ranges crossing the antimeridian).
    lat_min, lat_max
        Latitude range in degrees.
    spacing_km
        Target spacing in kilometers.
    center_lon_deg, center_lat_deg
        Cap center; the box center if omitted.
    max_distance_km
        Cap radius; if omitted, the largest distance from the center to
        the box corners and edge midpoints.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Nodes inside the box, with longitudes in the convention of
        ``lon_min``.

    Raises
    ------
    ValueError
        If ``spacing_km`` is not positive or the latitude range is
        invalid.
    """
    _require_positive(spacing_km, "spacing_km")
    _validate_lat_range(lat_min, lat_max)

    span = _lon_span(lon_min, lon_max)

    if center_lon_deg is None:
        center_lon_deg = lon_min + 0.5 * span

    if center_lat_deg is None:
        center_lat_deg = 0.5 * (lat_min + lat_max)

    if max_distance_km is None:
        lons = [lon_min, lon_min + 0.5 * span, lon_min + span]
        lats = [lat_min, center_lat_deg, lat_max]
        max_distance_km = max(
            great_circle_distance(center_lon_deg, center_lat_deg, x, y,
                                  radius_m=radius_m)
            for x in lons for y in lats
        ) / 1e3
        max_distance_km = max(max_distance_km, spacing_km)

    cap = mesh_fibonacci_spherical_cap(
        center_lon_deg=center_lon_deg,
        center_lat_deg=center_lat_deg,
        max_distance_km=max_distance_km,
        spacing_km=spacing_km,
        radius_m=radius_m,
    )

    lon, lat, _ = cap.to_arrays()
    offset = _lon_offset(lon, lon_min)
    keep = (offset <= span) & (lat >= lat_min) & (lat <= lat_max)

    return WgsMesh.from_arrays(lon_min + offset[keep], lat[keep])


# ---------------------------------------------------------------------------
# Regular grids
# ---------------------------------------------------------------------------

def mesh_cartesian(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    lon_step_deg: float,
    lat_step_deg: float,
) -> WgsMesh:
    """
    Generate a regular longitude-latitude grid.

    The grid does not preserve area and is meant for quick sampling,
    debugging or plotting. Nodes start at ``(lon_min, lat_min)`` and
    include the upper bounds when they fall on the grid.

    Parameters
    ----------
    lon_min, lon_max
        Longitude range in degrees (see the module conventions).
    lat_min, lat_max
        Latitude range in degrees.
    lon_step_deg, lat_step_deg
        Grid steps in degrees.

    Returns
    -------
    WgsMesh
        Grid nodes, row by row from south to north.

    Raises
    ------
    ValueError
        If a step is not positive or the latitude range is invalid.
    """
    if lon_step_deg <= 0.0 or lat_step_deg <= 0.0:
        raise ValueError("Step sizes must be positive.")

    _validate_lat_range(lat_min, lat_max)

    lons = lon_min + lon_step_deg * np.arange(
        _grid_count(_lon_span(lon_min, lon_max), lon_step_deg)
    )
    lats = lat_min + lat_step_deg * np.arange(
        _grid_count(lat_max - lat_min, lat_step_deg)
    )

    grid_lon, grid_lat = np.meshgrid(lons, lats)

    return WgsMesh.from_arrays(grid_lon.ravel(), grid_lat.ravel())


def mesh_sinusoidal(
    spacing_km: float,
    lon_min: Optional[float] = None,
    lon_max: Optional[float] = None,
    lat_min: Optional[float] = None,
    lat_max: Optional[float] = None,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate an equal-area grid in the sinusoidal projection.

    The grid has constant spacing in the sinusoidal (x, y) plane and is
    aligned on a global grid centered at (0°, 0°), so that meshes over
    different windows share their nodes. Only the nodes inside the
    requested window are computed.

    Parameters
    ----------
    spacing_km
        Grid spacing in kilometers; the area per node is about
        ``spacing_km ** 2``.
    lon_min, lon_max
        Longitude range in degrees; the whole globe if omitted (see the
        module conventions).
    lat_min, lat_max
        Latitude range in degrees; ``[-90, 90]`` if omitted.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Grid nodes, row by row from south to north, with longitudes in
        the convention of ``lon_min``.

    Raises
    ------
    ValueError
        If ``spacing_km`` is not positive or the latitude range is
        invalid.
    """
    _require_positive(spacing_km, "spacing_km")

    lon_min = -180.0 if lon_min is None else lon_min
    lon_max = 180.0 if lon_max is None else lon_max
    lat_min = -90.0 if lat_min is None else lat_min
    lat_max = 90.0 if lat_max is None else lat_max

    _validate_lat_range(lat_min, lat_max)

    step = spacing_km * 1e3
    span = _lon_span(lon_min, lon_max)
    windows = _lon_windows(lon_min, span)

    row_first = math.ceil(radius_m * lat_min * DEG2RAD / step
                          - _INDEX_TOLERANCE)
    row_last = math.floor(radius_m * lat_max * DEG2RAD / step
                          + _INDEX_TOLERANCE)

    lons: List[np.ndarray] = []
    lats: List[np.ndarray] = []

    for row in range(row_first, row_last + 1):
        lat_rad = row * step / radius_m
        scale = radius_m * max(math.cos(lat_rad), 1e-12)

        for start, end in windows:
            col_first = math.ceil(scale * start * DEG2RAD / step
                                  - _INDEX_TOLERANCE)
            col_last = math.floor(scale * end * DEG2RAD / step
                                  + _INDEX_TOLERANCE)

            if col_last < col_first:
                continue

            cols = np.arange(col_first, col_last + 1, dtype=float)
            lon = cols * step / scale * RAD2DEG
            lons.append(lon_min + _lon_offset(lon, lon_min))
            lats.append(np.full(cols.size, lat_rad * RAD2DEG))

    if not lons:
        return WgsMesh()

    return WgsMesh.from_arrays(np.concatenate(lons), np.concatenate(lats))


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _require_positive(value: float, name: str) -> None:
    """
    Raise ValueError if ``value`` is not strictly positive.
    """
    if not value > 0.0:
        raise ValueError(f"{name} must be positive.")


def _validate_lat_range(lat_min: float, lat_max: float) -> None:
    """
    Check that ``-90 <= lat_min <= lat_max <= 90``.
    """
    if not -90.0 <= lat_min <= lat_max <= 90.0:
        raise ValueError(
            "Latitude bounds must satisfy -90 <= lat_min <= lat_max <= 90."
        )


def _lon_span(lon_min: float, lon_max: float) -> float:
    """
    Return the eastward extent from ``lon_min`` to ``lon_max``.

    A range of 360 degrees or more is the whole globe; ``lon_max <
    lon_min`` means a range crossing the antimeridian.
    """
    if lon_max - lon_min >= 360.0:
        return 360.0

    return (lon_max - lon_min) % 360.0


def _lon_offset(lon: np.ndarray, lon_min: float) -> np.ndarray:
    """
    Return the eastward offset of ``lon`` from ``lon_min`` in [0, 360).

    Offsets within a tiny tolerance of 360 are mapped to 0, so that
    nodes on ``lon_min`` are not lost to rounding.
    """
    offset = (np.asarray(lon, dtype=float) - lon_min) % 360.0

    return np.where(offset > 360.0 - 1e-9, 0.0, offset)


def _lon_windows(lon_min: float, span: float) -> List[Tuple[float, float]]:
    """
    Split ``[lon_min, lon_min + span]`` into intervals within [-180, 180].

    The whole globe is a single interval, so that no meridian is
    generated twice.
    """
    if span >= 360.0:
        return [(-180.0, 180.0)]

    start = normalize_longitude(lon_min)
    end = start + span

    if end <= 180.0:
        return [(start, end)]

    return [(start, 180.0), (-180.0, end - 360.0)]


def _grid_count(extent: float, step: float) -> int:
    """
    Return the number of grid nodes from 0 to ``extent`` with ``step``.
    """
    return int(math.floor(extent / step + _INDEX_TOLERANCE)) + 1
