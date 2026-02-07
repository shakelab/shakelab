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
Core geodetic utilities for ShakeLab.

This module contains low–level geodetic constants and functions that do
not depend on higher–level geometry or mesh classes. It is intended to
be small, stable and reusable across the library.

Conventions
----------
* Geographic coordinates are always given as (longitude, latitude) in
  decimal degrees, unless explicitly stated otherwise.
* Distances are returned in meters by default.
* Angles in the internal computations are expressed in radians.
"""

from __future__ import annotations

from typing import Tuple, Union

import re
import numpy as np

Number = Union[float, int]
ArrayLike = Union[Number, np.ndarray]


# -----------------------------------------------------------------------------
# Geodetic constants
# -----------------------------------------------------------------------------


#: Mean Earth radius in meters (IUGG 2015 value, close to WGS84 mean).
MEAN_EARTH_RADIUS_M: float = 6_371_008.8

#: WGS84 equatorial (semi–major) radius in meters.
EQUATORIAL_EARTH_RADIUS_M: float = 6_378_137.0

#: WGS84 polar (semi–minor) radius in meters.
POLAR_EARTH_RADIUS_M: float = 6_356_752.3

#: Degrees–to–radians conversion factor.
DEG2RAD: float = np.pi / 180.0

#: Radians–to–degrees conversion factor.
RAD2DEG: float = 180.0 / np.pi


# -----------------------------------------------------------------------------
# Basic angle helpers
# -----------------------------------------------------------------------------


def deg_to_rad(angle_deg: ArrayLike) -> ArrayLike:
    """
    Convert degrees to radians.

    Parameters
    ----------
    angle_deg : float or array_like
        Angle in degrees.

    Returns
    -------
    float or numpy.ndarray
        Angle in radians with the same shape as the input.
    """
    return np.asarray(angle_deg) * DEG2RAD


def rad_to_deg(angle_rad: ArrayLike) -> ArrayLike:
    """
    Convert radians to degrees.

    Parameters
    ----------
    angle_rad : float or array_like
        Angle in radians.

    Returns
    -------
    float or numpy.ndarray
        Angle in degrees with the same shape as the input.
    """
    return np.asarray(angle_rad) * RAD2DEG


def normalize_longitude(lon_deg: ArrayLike) -> ArrayLike:
    """
    Normalize longitude to the range [-180, 180).

    Parameters
    ----------
    lon_deg : float or array_like
        Longitude in degrees.

    Returns
    -------
    float or numpy.ndarray
        Normalized longitude(s) in degrees.
    """
    lon = np.asarray(lon_deg)
    lon_norm = (lon + 180.0) % 360.0 - 180.0
    return lon_norm


def dms2dec(value: ArrayLike) -> ArrayLike:
    """
    Convert degrees-minutes-seconds notation to decimal degrees.

    This helper accepts both plain decimal strings (e.g. ``"12.5"``)
    and DMS-like strings with arbitrary separators, such as::

        "12°34'56.7\""
        "12 34 56.7N"
        "-12:34:56.7"
        "12°34'S"

    If a hemisphere letter is present (N, S, E, W), it controls the
    sign of the result and overrides any explicit sign on the degrees.
    For numeric inputs, the value is returned unchanged as float.

    Parameters
    ----------
    value : float, int or str
        Angle in decimal degrees or DMS notation.

    Returns
    -------
    float or numpy.ndarray
        Angle in decimal degrees.
    """
    # If already numeric, just cast to float.
    if not isinstance(value, str):
        return float(value)

    s = value.strip().upper().replace(",", ".")
    if not s:
        raise ValueError("Empty DMS string.")

    # Detect optional hemisphere and its sign.
    hemi_sign = 1.0
    has_hemi = False
    if s[-1] in ("N", "E", "S", "W"):
        hemi = s[-1]
        has_hemi = True
        s = s[:-1].strip()
        if hemi in ("S", "W"):
            hemi_sign = -1.0

    # Split on any non-numeric separator.
    parts = [p for p in re.split(r"[^\d.+-]+", s) if p]
    if not parts:
        raise ValueError(f"Cannot parse DMS string: {value!r}")

    nums = [float(p) for p in parts]

    # If a hemisphere is present, ignore the numeric sign: the
    # hemisphere fully controls the final sign.
    if has_hemi:
        nums = [abs(x) for x in nums]

    if len(nums) == 1:
        deg = nums[0]
        dec = deg
    elif len(nums) == 2:
        deg, minute = nums
        dec = deg + minute / 60.0
    else:
        deg, minute, sec = nums[:3]
        dec = deg + minute / 60.0 + sec / 3600.0

    return float(dec * hemi_sign)


# -----------------------------------------------------------------------------
# Geometry calculators
# -----------------------------------------------------------------------------


def great_circle_distance(
    lon1_deg: Number,
    lat1_deg: Number,
    lon2_deg: Number,
    lat2_deg: Number,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute great–circle distance between two points on a sphere.

    The calculation uses the Haversine formula and assumes a spherical
    Earth with radius ``radius_m``.

    Parameters
    ----------
    lon1_deg, lat1_deg : float
        Longitude and latitude of the first point in degrees.
    lon2_deg, lat2_deg : float
        Longitude and latitude of the second point in degrees.
    radius_m : float, optional
        Sphere radius in meters. Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    float
        Great–circle distance between the two points in meters.
    """
    lon1 = lon1_deg * DEG2RAD
    lat1 = lat1_deg * DEG2RAD
    lon2 = lon2_deg * DEG2RAD
    lat2 = lat2_deg * DEG2RAD

    dlon = 0.5 * (lon2 - lon1)
    dlat = 0.5 * (lat2 - lat1)

    sin_dlat = np.sin(dlat)
    sin_dlon = np.sin(dlon)

    a = (
        sin_dlat * sin_dlat
        + sin_dlon * sin_dlon * np.cos(lat1) * np.cos(lat2)
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    return float(radius_m * c)


def chord_distance(
    lon1_deg: Number,
    lat1_deg: Number,
    lon2_deg: Number,
    lat2_deg: Number,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute chord (tunnel) distance between two points on a sphere.

    The chord distance is the straight–line distance through the
    Earth between the two points, assuming a spherical Earth with
    radius ``radius_m``.

    Parameters
    ----------
    lon1_deg, lat1_deg : float
        Longitude and latitude of the first point in degrees.
    lon2_deg, lat2_deg : float
        Longitude and latitude of the second point in degrees.
    radius_m : float, optional
        Sphere radius in meters. Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    float
        Straight–line distance between the two points in meters.
    """
    lon1 = lon1_deg * DEG2RAD
    lat1 = lat1_deg * DEG2RAD
    lon2 = lon2_deg * DEG2RAD
    lat2 = lat2_deg * DEG2RAD

    # Unit vectors for the two points
    x1 = np.cos(lat1) * np.cos(lon1)
    y1 = np.cos(lat1) * np.sin(lon1)
    z1 = np.sin(lat1)

    x2 = np.cos(lat2) * np.cos(lon2)
    y2 = np.cos(lat2) * np.sin(lon2)
    z2 = np.sin(lat2)

    dot = x1 * x2 + y1 * y2 + z1 * z2
    dot = np.clip(dot, -1.0, 1.0)

    angle = np.arccos(dot)
    chord = 2.0 * radius_m * np.sin(0.5 * angle)

    return float(chord)


def tunnel_distance_sphere(
    lon1_deg: Number,
    lat1_deg: Number,
    lon2_deg: Number,
    lat2_deg: Number,
    elev1_m: Number = 0.0,
    elev2_m: Number = 0.0,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute 3-D straight-line (tunnel) distance between two points on a
    spherical Earth, including elevation.

    This helper converts both geographic points to Cartesian coordinates
    using :func:`wgs84_to_cartesian_sphere` and returns the Euclidean
    distance between the two 3-D locations. Compared to
    :func:`chord_distance`, this function also accounts explicitly for
    the elevation of each point above the reference sphere.

    Parameters
    ----------
    lon1_deg, lat1_deg : float
        Longitude and latitude of the first point in degrees.
    lon2_deg, lat2_deg : float
        Longitude and latitude of the second point in degrees.
    elev1_m, elev2_m : float, optional
        Elevation of the two points in meters above the reference
        sphere. Defaults to 0.0.
    radius_m : float, optional
        Sphere radius in meters. Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    float
        3-D Euclidean distance between the two points in meters.
    """
    x1, y1, z1 = wgs84_to_cartesian_sphere(
        lon1_deg, lat1_deg, elev1_m, radius_m
    )
    x2, y2, z2 = wgs84_to_cartesian_sphere(
        lon2_deg, lat2_deg, elev2_m, radius_m
    )

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    return float(np.sqrt(dx * dx + dy * dy + dz * dz))


def forward_azimuth(
    lon1_deg: Number,
    lat1_deg: Number,
    lon2_deg: Number,
    lat2_deg: Number,
) -> float:
    """
    Compute forward azimuth from point 1 to point 2.

    The azimuth is measured clockwise from geographic north (0°)
    and returned in degrees in the range [0, 360).

    Parameters
    ----------
    lon1_deg, lat1_deg : float
        Longitude and latitude of the start point in degrees.
    lon2_deg, lat2_deg : float
        Longitude and latitude of the end point in degrees.

    Returns
    -------
    float
        Forward azimuth from point 1 to point 2 in degrees.
    """
    lon1 = lon1_deg * DEG2RAD
    lat1 = lat1_deg * DEG2RAD
    lon2 = lon2_deg * DEG2RAD
    lat2 = lat2_deg * DEG2RAD

    dlon = lon2 - lon1

    x = np.sin(dlon) * np.cos(lat2)
    y = (
        np.cos(lat1) * np.sin(lat2)
        - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    )

    az = np.arctan2(x, y) * RAD2DEG
    az_norm = (az + 360.0) % 360.0

    return float(az_norm)


def geocentric_radius(lat_deg: ArrayLike) -> ArrayLike:
    """
    Geocentric radius at a given geodetic latitude (WGS84 ellipsoid).

    This function computes the distance from the Earth's center to the
    ellipsoid surface at the specified geodetic latitude, using the
    WGS84 semi-major and semi-minor axes.

    Parameters
    ----------
    lat_deg : float or array_like
        Geodetic latitude in degrees (positive north).

    Returns
    -------
    float or numpy.ndarray
        Geocentric radius in meters, with the same shape as the input.
    """
    lat_rad = np.asarray(lat_deg) * DEG2RAD
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)

    a2 = EQUATORIAL_EARTH_RADIUS_M ** 2
    b2 = POLAR_EARTH_RADIUS_M ** 2

    num = (a2 * a2 * cos_lat * cos_lat) + (b2 * b2 * sin_lat * sin_lat)
    den = (a2 * cos_lat * cos_lat) + (b2 * sin_lat * sin_lat)

    radius = np.sqrt(num / den)
    return radius


# -----------------------------------------------------------------------------
# Coordinate transforms
# -----------------------------------------------------------------------------


def wgs84_to_cartesian_sphere(
    lon_deg: ArrayLike,
    lat_deg: ArrayLike,
    elevation_m: ArrayLike = 0.0,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert WGS84 coordinates to Cartesian coordinates (spherical Earth).

    Parameters
    ----------
    lon_deg, lat_deg : float or array_like
        Longitude and latitude in degrees.
    elevation_m : float or array_like, optional
        Elevation above the reference sphere in meters. Default is 0.
    radius_m : float, optional
        Sphere radius in meters. Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    x, y, z : numpy.ndarray
        Cartesian coordinates in meters in a right–handed system where
        the origin is at the sphere center, the z–axis points to the
        geographic north pole and the x–axis intersects the equator
        at longitude 0°.
    """
    lon = np.asarray(lon_deg) * DEG2RAD
    lat = np.asarray(lat_deg) * DEG2RAD
    ele = np.asarray(elevation_m)

    rho = radius_m + ele

    cos_lat = np.cos(lat)
    sin_lat = np.sin(lat)
    cos_lon = np.cos(lon)
    sin_lon = np.sin(lon)

    x = rho * cos_lat * cos_lon
    y = rho * cos_lat * sin_lon
    z = rho * sin_lat

    return x, y, z


def cartesian_sphere_to_wgs84(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Cartesian coordinates (spherical Earth) to WGS84.

    Parameters
    ----------
    x, y, z : float or array_like
        Cartesian coordinates in meters in the same system used by
        :func:`wgs84_to_cartesian_sphere`.
    radius_m : float, optional
        Sphere radius in meters used as reference for the elevation.

    Returns
    -------
    lon_deg, lat_deg, elevation_m : numpy.ndarray
        Longitude and latitude in degrees, and elevation above the
        reference sphere in meters.
    """
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)
    z_arr = np.asarray(z)

    rho = np.sqrt(x_arr * x_arr + y_arr * y_arr + z_arr * z_arr)
    lat = np.arcsin(z_arr / rho)
    lon = np.arctan2(y_arr, x_arr)

    lon_deg = lon * RAD2DEG
    lat_deg = lat * RAD2DEG
    elevation_m = rho - radius_m

    return lon_deg, lat_deg, elevation_m


def wgs84_to_sinusoidal(
    lon_deg: ArrayLike,
    lat_deg: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert WGS84 coordinates to a global sinusoidal projection.

    This projection is equal–area and is used here as a simple tool
    to estimate areas on the sphere.

    Parameters
    ----------
    lon_deg, lat_deg : float or array_like
        Longitude and latitude in degrees.
    radius_m : float, optional
        Sphere radius in meters. Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    x, y : numpy.ndarray
        Projected coordinates in meters.
    """
    lon = np.asarray(lon_deg) * DEG2RAD
    lat = np.asarray(lat_deg) * DEG2RAD

    x = radius_m * lon * np.cos(lat)
    y = radius_m * lat

    return x, y


def sinusoidal_to_wgs84(
    x: ArrayLike,
    y: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert coordinates from the sinusoidal projection to WGS84.

    Parameters
    ----------
    x, y : float or array_like
        Projected coordinates in meters.
    radius_m : float, optional
        Sphere radius in meters used in the projection.

    Returns
    -------
    lon_deg, lat_deg : numpy.ndarray
        Longitude and latitude in degrees.
    """
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    lat = y_arr / radius_m
    lon = x_arr / (radius_m * np.cos(lat))

    lon_deg = lon * RAD2DEG
    lat_deg = lat * RAD2DEG

    return lon_deg, lat_deg

