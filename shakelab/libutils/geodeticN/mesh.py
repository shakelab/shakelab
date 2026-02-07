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
Mesh generation utilities for ShakeLab.

This module provides several algorithms to generate WGS84 spatial
meshes, including:
* Fibonacci global mesh (quasi-equal-area on full sphere)
* Fibonacci spherical cap mesh (efficient local equal-area mesh)
* Spherical cap mesh clipped to a bounding box
* Regular Cartesian mesh in latitude and longitude

All functions return mesh nodes as lists of WgsPoint or WgsMesh
containers, without any projection or I/O. Coordinates are always
expressed as (lon, lat) in degrees.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    DEG2RAD,
    RAD2DEG,
    sinusoidal_to_wgs84,
)
from shakelab.libutils.geodeticN.primitives import (
    WgsPoint,
    WgsMesh,
)


# ============================================================================
# Global Fibonacci mesh
# ============================================================================


def mesh_fibonacci_global(
    spacing_km: float,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a global Fibonacci (golden-spiral) mesh on the sphere.

    Points are distributed quasi-uniformly in area over the entire
    sphere. This function is intended for global or continental
    grids; for local regions, use spherical_cap_mesh_fibonacci.

    Parameters
    ----------
    spacing_km : float
        Target spacing between nodes in kilometers. The average
        cell area is approximately (spacing_km**2).
    radius_m : float, optional
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Mesh containing the generated nodes.
    """
    if spacing_km <= 0.0:
        raise ValueError("spacing_km must be positive.")

    a_target = spacing_km * spacing_km * 1e6
    a_sphere = 4.0 * np.pi * (radius_m ** 2)
    n_pts = max(1, int(round(a_sphere / a_target)))

    golden_angle = np.pi * (3.0 - np.sqrt(5.0))

    pts: List[WgsPoint] = []
    for i in range(n_pts):
        z = 1.0 - 2.0 * (i + 0.5) / n_pts
        lat = np.arcsin(z)
        lon = golden_angle * i

        lon_deg = float((lon * RAD2DEG + 180.0) % 360.0 - 180.0)
        lat_deg = float(lat * RAD2DEG)

        pts.append(WgsPoint(lon_deg, lat_deg))

    return WgsMesh(points=pts)


# ============================================================================
# Local Fibonacci cap mesh (Equal-area on a spherical cap)
# ============================================================================


def mesh_fibonacci_spherical_cap(
    center_lon_deg: float,
    center_lat_deg: float,
    max_distance_km: float,
    spacing_km: float,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a quasi-equal-area Fibonacci mesh on a spherical cap.

    The cap is centered at (center_lon_deg, center_lat_deg), with a
    radius max_distance_km (great-circle distance). Points are
    distributed using a Fibonacci spiral restricted to the cap,
    without generating the whole sphere.

    Parameters
    ----------
    center_lon_deg : float
        Longitude of cap center in degrees.
    center_lat_deg : float
        Latitude of cap center in degrees.
    max_distance_km : float
        Cap radius in kilometers.
    spacing_km : float
        Desired linear spacing in kilometers.
    radius_m : float, optional
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Mesh containing the generated nodes.
    """
    if spacing_km <= 0.0:
        raise ValueError("spacing_km must be positive.")
    if max_distance_km <= 0.0:
        raise ValueError("max_distance_km must be positive.")

    alpha = max_distance_km * 1000.0 / radius_m

    a_cap = 2.0 * np.pi * (radius_m ** 2) * (1.0 - np.cos(alpha))
    a_node = spacing_km * spacing_km * 1e6
    n_pts = max(1, int(round(a_cap / a_node)))

    lat0 = center_lat_deg * DEG2RAD
    lon0 = center_lon_deg * DEG2RAD
    sin_lat0 = np.sin(lat0)
    cos_lat0 = np.cos(lat0)

    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    z_min = np.cos(alpha)

    pts: List[WgsPoint] = []

    for i in range(n_pts):
        u = (i + 0.5) / n_pts
        z = 1.0 - u * (1.0 - z_min)
        theta = np.arccos(z)
        phi = golden_angle * i

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        sin_lat = (
            sin_lat0 * cos_theta +
            cos_lat0 * sin_theta * cos_phi
        )
        lat = np.arcsin(sin_lat)

        y = sin_phi * sin_theta * cos_lat0
        x = cos_theta - sin_lat0 * sin_lat
        lon = lon0 + np.arctan2(y, x)

        pts.append(WgsPoint(float(lon * RAD2DEG),
                            float(lat * RAD2DEG)))

    return WgsMesh(points=pts)


# ============================================================================
# Fibonacci cap mesh with bounding-box clipping
# ============================================================================


def mesh_fibonacci_bbox(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    spacing_km: float,
    center_lon_deg: float | None = None,
    center_lat_deg: float | None = None,
    max_distance_km: float | None = None,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate a Fibonacci cap mesh and clip nodes to a bounding box.

    The center and cap radius can be automatically inferred if not
    provided.

    Parameters
    ----------
    lon_min, lon_max : float
        Longitude bounds in degrees.
    lat_min, lat_max : float
        Latitude bounds in degrees.
    spacing_km : float
        Target spacing in kilometers.
    center_lon_deg, center_lat_deg : float, optional
        Center of the cap. If None, the box center is used.
    max_distance_km : float, optional
        Cap radius. If None, it is computed as the greatest
        distance from the center to the four corners of the box.
    radius_m : float, optional
        Sphere radius in meters.

    Returns
    -------
    WgsMesh
        Mesh containing the clipped nodes.
    """
    if spacing_km <= 0.0:
        raise ValueError("spacing_km must be positive.")

    if center_lon_deg is None:
        center_lon_deg = lon_min + 0.5 * (lon_max - lon_min)

    if center_lat_deg is None:
        center_lat_deg = lat_min + 0.5 * (lat_max - lat_min)

    if max_distance_km is None:
        corners = [
            (lon_min, lat_min),
            (lon_min, lat_max),
            (lon_max, lat_min),
            (lon_max, lat_max),
        ]
        dist_max = 0.0
        for lon_c, lat_c in corners:
            d = _gc_km(center_lon_deg, center_lat_deg, lon_c, lat_c)
            if d > dist_max:
                dist_max = d
        max_distance_km = dist_max

    raw_mesh = mesh_fibonacci_spherical_cap(
        center_lon_deg=center_lon_deg,
        center_lat_deg=center_lat_deg,
        max_distance_km=max_distance_km,
        spacing_km=spacing_km,
        radius_m=radius_m,
    )

    pts = []
    for p in raw_mesh:
        if lon_min <= p.longitude <= lon_max:
            if lat_min <= p.latitude <= lat_max:
                pts.append(p)

    return WgsMesh(points=pts)


def _gc_km(
    lon1: float,
    lat1: float,
    lon2: float,
    lat2: float,
) -> float:
    """
    Internal helper: great-circle distance in kilometers.
    """
    lon1r = lon1 * DEG2RAD
    lat1r = lat1 * DEG2RAD
    lon2r = lon2 * DEG2RAD
    lat2r = lat2 * DEG2RAD

    dlon = lon2r - lon1r
    dlat = lat2r - lat1r

    a = (
        np.sin(dlat * 0.5) ** 2
        + np.cos(lat1r) * np.cos(lat2r) *
        np.sin(dlon * 0.5) ** 2
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return float((MEAN_EARTH_RADIUS_M * c) / 1000.0)


# ============================================================================
# Regular Cartesian mesh (lat/lon grid)
# ============================================================================


def mesh_cartesian(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    lon_step_deg: float,
    lat_step_deg: float,
) -> WgsMesh:
    """
    Generate a regular mesh in geographic coordinates.

    This is a simple lat/lon grid. It does not preserve area and is
    intended mainly for quick sampling in debugging or plotting.

    Parameters
    ----------
    lon_min, lon_max : float
        Longitude bounds in degrees.
    lat_min, lat_max : float
        Latitude bounds in degrees.
    lon_step_deg, lat_step_deg : float
        Step sizes in degrees.

    Returns
    -------
    WgsMesh
        Regular mesh of points.
    """
    if lon_step_deg <= 0.0 or lat_step_deg <= 0.0:
        raise ValueError("Step sizes must be positive.")

    lons = np.arange(lon_min, lon_max + 1e-12, lon_step_deg)
    lats = np.arange(lat_min, lat_max + 1e-12, lat_step_deg)

    pts = []
    for lat in lats:
        for lon in lons:
            pts.append(WgsPoint(float(lon), float(lat)))

    return WgsMesh(points=pts)


# ============================================================================
# Equal-area sinusoidal mesh (centered at 0°, 0°)
# ============================================================================


def mesh_sinusoidal(
    spacing_km: float,
    lon_min: Optional[float] = None,
    lon_max: Optional[float] = None,
    lat_min: Optional[float] = None,
    lat_max: Optional[float] = None,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> WgsMesh:
    """
    Generate an approximately equal-area mesh using the sinusoidal
    projection, with the grid centered at (0°, 0°).

    The algorithm builds a regular grid in the sinusoidal (x, y)
    plane with constant spacing in meters, symmetric around the
    projection origin (lon=0°, lat=0°). The geographic bounds, if
    provided, are used only to *clip* the resulting nodes.

    If geographic bounds are not provided, the clipping window
    defaults to the whole globe (lon in [-180, 180], lat in
    [-90, 90]).

    Parameters
    ----------
    spacing_km : float
        Target linear spacing between nodes in kilometers. The
        approximate cell area is (spacing_km**2) km².
    lon_min, lon_max : float, optional
        Longitude bounds in degrees for clipping. If None, default
        to -180 and 180 degrees, respectively.
    lat_min, lat_max : float, optional
        Latitude bounds in degrees for clipping. If None, default
        to -90 and 90 degrees, respectively.
    radius_m : float, optional
        Sphere radius in meters used in the sinusoidal projection.
        Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    WgsMesh
        Mesh containing the generated nodes.

    Raises
    ------
    ValueError
        If spacing_km is not positive or if the geographic bounds
        are invalid.
    """
    if spacing_km <= 0.0:
        raise ValueError("spacing_km must be positive.")

    if lon_min is None:
        lon_min = -180.0
    if lon_max is None:
        lon_max = 180.0
    if lat_min is None:
        lat_min = -90.0
    if lat_max is None:
        lat_max = 90.0

    if lon_max <= lon_min:
        raise ValueError("lon_max must be greater than lon_min.")
    if lat_max <= lat_min:
        raise ValueError("lat_max must be greater than lat_min.")

    dx = spacing_km * 1000.0
    dy = spacing_km * 1000.0

    # Determine symmetric latitude extent (in degrees) around 0.
    max_lat_abs_deg = max(abs(lat_min), abs(lat_max))
    if max_lat_abs_deg <= 0.0:
        max_lat_abs_deg = spacing_km / 111.0

    lat_max_rad = max_lat_abs_deg * DEG2RAD
    y_max = radius_m * lat_max_rad

    # Symmetric y grid: from -y_max to +y_max.
    n_y = max(1, int(np.floor(y_max / dy)))
    ys = dy * np.arange(-n_y, n_y + 1, dtype=float)

    # Determine symmetric longitude extent (in degrees) around 0.
    max_lon_abs_deg = max(abs(lon_min), abs(lon_max))
    if max_lon_abs_deg <= 0.0:
        max_lon_abs_deg = spacing_km / 111.0

    lon_max_rad_abs = max_lon_abs_deg * DEG2RAD

    pts: List[WgsPoint] = []

    for y in ys:
        lat_rad = y / radius_m

        # Avoid numerical issues extremely close to the poles.
        cos_lat = np.cos(lat_rad)
        if cos_lat <= 0.0:
            cos_lat = max(cos_lat, 1e-12)

        x_max = radius_m * lon_max_rad_abs * cos_lat
        if x_max <= 0.0:
            continue

        # Symmetric x grid: from -x_max to +x_max.
        n_x = max(1, int(np.floor(x_max / dx)))
        xs = dx * np.arange(-n_x, n_x + 1, dtype=float)

        # Ensure we stay within [-x_max, x_max].
        xs = np.clip(xs, -x_max, x_max)

        ys_row = np.full_like(xs, y, dtype=float)

        lon_arr, lat_arr = sinusoidal_to_wgs84(
            xs,
            ys_row,
            radius_m=radius_m,
        )

        for lon_deg, lat_deg in zip(lon_arr, lat_arr):
            if lon_deg < lon_min - 1e-8 or lon_deg > lon_max + 1e-8:
                continue
            if lat_deg < lat_min - 1e-8 or lat_deg > lat_max + 1e-8:
                continue
            pts.append(
                WgsPoint(
                    longitude=float(lon_deg),
                    latitude=float(lat_deg),
                )
            )

    return WgsMesh(points=pts)
