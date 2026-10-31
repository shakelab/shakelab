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
Core geodetic constants and functions.

This module contains low-level geodetic constants and functions that do
not depend on the geometry classes of :mod:`geodeticN`. It is meant to
be small, stable and reusable across the library.

Design choices
--------------
- The WGS84 ellipsoid parameters are defined here once and reused by
  the other modules of the package.
- Functions accept scalars or array-like inputs. Scalar inputs return
  Python floats; array inputs return NumPy arrays.
- Spherical distances are computed from the haversine term, which is
  numerically stable at short distances and clipped against rounding
  near antipodal points.

Conventions
-----------
- Geographic coordinates are given as (longitude, latitude) in decimal
  degrees. Longitudes may use either the ``[-180, 180)`` or the
  ``[0, 360)`` convention; angular computations are periodic and do
  not depend on it.
- Distances and elevations are in meters; elevation is positive upward.
- Angles in internal computations are in radians.
"""

from __future__ import annotations

from numbers import Real
import re
from typing import Any, Tuple, Union

import numpy as np


__all__ = [
    "MEAN_EARTH_RADIUS_M",
    "WGS84_A",
    "WGS84_F",
    "WGS84_B",
    "WGS84_E2",
    "WGS84_EP2",
    "EQUATORIAL_EARTH_RADIUS_M",
    "POLAR_EARTH_RADIUS_M",
    "DEG2RAD",
    "RAD2DEG",
    "deg_to_rad",
    "rad_to_deg",
    "normalize_longitude",
    "relative_longitude",
    "dms2dec",
    "great_circle_distance",
    "chord_distance",
    "tunnel_distance_sphere",
    "forward_azimuth",
    "geocentric_radius",
    "wgs84_to_cartesian_sphere",
    "cartesian_sphere_to_wgs84",
    "wgs84_to_sinusoidal",
    "sinusoidal_to_wgs84",
]


Number = Union[float, int]
ArrayLike = Union[Number, np.ndarray]


# ---------------------------------------------------------------------------
# Geodetic constants
# ---------------------------------------------------------------------------

#: Mean Earth radius in meters (IUGG value, close to the WGS84 mean).
MEAN_EARTH_RADIUS_M: float = 6_371_008.8

#: WGS84 semi-major axis in meters.
WGS84_A: float = 6_378_137.0

#: WGS84 flattening.
WGS84_F: float = 1.0 / 298.257223563

#: WGS84 semi-minor axis in meters.
WGS84_B: float = WGS84_A * (1.0 - WGS84_F)

#: WGS84 first eccentricity squared.
WGS84_E2: float = WGS84_F * (2.0 - WGS84_F)

#: WGS84 second eccentricity squared.
WGS84_EP2: float = WGS84_E2 / (1.0 - WGS84_E2)

#: WGS84 equatorial radius in meters (alias of :data:`WGS84_A`).
EQUATORIAL_EARTH_RADIUS_M: float = WGS84_A

#: WGS84 polar radius in meters (alias of :data:`WGS84_B`).
POLAR_EARTH_RADIUS_M: float = WGS84_B

#: Degrees-to-radians conversion factor.
DEG2RAD: float = np.pi / 180.0

#: Radians-to-degrees conversion factor.
RAD2DEG: float = 180.0 / np.pi


# ---------------------------------------------------------------------------
# Angle helpers
# ---------------------------------------------------------------------------

def deg_to_rad(angle_deg: ArrayLike) -> ArrayLike:
    """
    Convert degrees to radians.

    Parameters
    ----------
    angle_deg
        Angle in degrees, scalar or array-like.

    Returns
    -------
    float or numpy.ndarray
        Angle in radians, with the shape of the input.
    """
    return _to_output(np.asarray(angle_deg, dtype=float) * DEG2RAD)


def rad_to_deg(angle_rad: ArrayLike) -> ArrayLike:
    """
    Convert radians to degrees.

    Parameters
    ----------
    angle_rad
        Angle in radians, scalar or array-like.

    Returns
    -------
    float or numpy.ndarray
        Angle in degrees, with the shape of the input.
    """
    return _to_output(np.asarray(angle_rad, dtype=float) * RAD2DEG)


def normalize_longitude(lon_deg: ArrayLike) -> ArrayLike:
    """
    Normalize longitude to the range ``[-180, 180)``.

    Parameters
    ----------
    lon_deg
        Longitude in degrees, in any convention.

    Returns
    -------
    float or numpy.ndarray
        Normalized longitude in degrees.
    """
    lon = np.asarray(lon_deg, dtype=float)

    return _to_output((lon + 180.0) % 360.0 - 180.0)


def relative_longitude(
    lon_deg: ArrayLike,
    ref_lon_deg: ArrayLike,
) -> ArrayLike:
    """
    Return the longitude difference ``lon - ref`` in ``[-180, 180)``.

    The result is independent of the longitude convention of the inputs
    (``[-180, 180)`` or ``[0, 360)``) and is continuous across the
    antimeridian around ``ref_lon_deg``.

    Parameters
    ----------
    lon_deg
        Longitude in degrees.
    ref_lon_deg
        Reference longitude in degrees.

    Returns
    -------
    float or numpy.ndarray
        Wrapped longitude difference in degrees.
    """
    lon = np.asarray(lon_deg, dtype=float)
    ref = np.asarray(ref_lon_deg, dtype=float)

    return _to_output((lon - ref + 180.0) % 360.0 - 180.0)


def dms2dec(value: Any) -> Union[float, np.ndarray]:
    """
    Convert degrees-minutes-seconds notation to decimal degrees.

    Plain numbers and decimal strings (e.g. ``"12.5"``) are accepted as
    well as DMS strings with arbitrary separators, such as::

        "12°34'56.7\""
        "12 34 56.7N"
        "-12:34:56.7"
        "12°34'S"

    A leading ``-`` applies to the whole angle (``"-12:30"`` is
    -12.5). A trailing hemisphere letter (N, S, E, W) sets the sign and
    overrides any explicit sign. A comma is accepted as decimal
    separator.

    Parameters
    ----------
    value
        Angle as a number, a string, or a sequence of those.

    Returns
    -------
    float or numpy.ndarray
        Angle in decimal degrees; an array for sequence inputs.

    Raises
    ------
    TypeError
        If the input is neither a number, a string nor a sequence.
    ValueError
        If a string cannot be parsed or has more than three fields.
    """
    if isinstance(value, str):
        return _dms_string_to_decimal(value)

    if isinstance(value, (list, tuple, np.ndarray)):
        return np.array([dms2dec(item) for item in value], dtype=float)

    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("DMS value must be a number or a string.")

    return float(value)


# ---------------------------------------------------------------------------
# Distances and azimuth on the sphere
# ---------------------------------------------------------------------------

def great_circle_distance(
    lon1_deg: ArrayLike,
    lat1_deg: ArrayLike,
    lon2_deg: ArrayLike,
    lat2_deg: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> ArrayLike:
    """
    Compute the great-circle distance between points on a sphere.

    The haversine formula is used, with the haversine term clipped to
    ``[0, 1]`` so that nearly antipodal points never yield NaN.

    Parameters
    ----------
    lon1_deg, lat1_deg
        Longitude and latitude of the first point(s) in degrees.
    lon2_deg, lat2_deg
        Longitude and latitude of the second point(s) in degrees.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    float or numpy.ndarray
        Great-circle distance in meters.
    """
    hav = _haversine_term(lon1_deg, lat1_deg, lon2_deg, lat2_deg)
    angle = 2.0 * np.arctan2(np.sqrt(hav), np.sqrt(1.0 - hav))

    return _to_output(radius_m * angle)


def chord_distance(
    lon1_deg: ArrayLike,
    lat1_deg: ArrayLike,
    lon2_deg: ArrayLike,
    lat2_deg: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> ArrayLike:
    """
    Compute the chord (straight-line) distance between points on a
    sphere.

    The chord is computed as ``2 R sqrt(h)``, where ``h`` is the
    haversine term, which keeps full precision at short distances.

    Parameters
    ----------
    lon1_deg, lat1_deg
        Longitude and latitude of the first point(s) in degrees.
    lon2_deg, lat2_deg
        Longitude and latitude of the second point(s) in degrees.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    float or numpy.ndarray
        Chord distance in meters.
    """
    hav = _haversine_term(lon1_deg, lat1_deg, lon2_deg, lat2_deg)

    return _to_output(2.0 * radius_m * np.sqrt(hav))


def tunnel_distance_sphere(
    lon1_deg: ArrayLike,
    lat1_deg: ArrayLike,
    lon2_deg: ArrayLike,
    lat2_deg: ArrayLike,
    elev1_m: ArrayLike = 0.0,
    elev2_m: ArrayLike = 0.0,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> ArrayLike:
    """
    Compute the 3-D straight-line distance between points on a
    spherical Earth, including elevation.

    This is the coordinate-based counterpart of
    :func:`shakelab.libutils.geodeticN.distance.tunnel_distance_sphere`,
    which takes point objects.

    Parameters
    ----------
    lon1_deg, lat1_deg
        Longitude and latitude of the first point(s) in degrees.
    lon2_deg, lat2_deg
        Longitude and latitude of the second point(s) in degrees.
    elev1_m, elev2_m
        Elevations above the reference sphere in meters.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    float or numpy.ndarray
        3-D Euclidean distance in meters.
    """
    x1, y1, z1 = wgs84_to_cartesian_sphere(
        lon1_deg, lat1_deg, elev1_m, radius_m,
    )
    x2, y2, z2 = wgs84_to_cartesian_sphere(
        lon2_deg, lat2_deg, elev2_m, radius_m,
    )

    return _to_output(np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2
                              + (z2 - z1) ** 2))


def forward_azimuth(
    lon1_deg: ArrayLike,
    lat1_deg: ArrayLike,
    lon2_deg: ArrayLike,
    lat2_deg: ArrayLike,
) -> ArrayLike:
    """
    Compute the forward azimuth from point 1 to point 2 on a sphere.

    Parameters
    ----------
    lon1_deg, lat1_deg
        Longitude and latitude of the start point(s) in degrees.
    lon2_deg, lat2_deg
        Longitude and latitude of the end point(s) in degrees.

    Returns
    -------
    float or numpy.ndarray
        Azimuth in degrees, clockwise from north, in ``[0, 360)``.
    """
    lon1 = np.asarray(lon1_deg, dtype=float) * DEG2RAD
    lat1 = np.asarray(lat1_deg, dtype=float) * DEG2RAD
    lon2 = np.asarray(lon2_deg, dtype=float) * DEG2RAD
    lat2 = np.asarray(lat2_deg, dtype=float) * DEG2RAD

    dlon = lon2 - lon1

    x = np.sin(dlon) * np.cos(lat2)
    y = (
        np.cos(lat1) * np.sin(lat2)
        - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    )

    azimuth = np.arctan2(x, y) * RAD2DEG

    return _to_output((azimuth + 360.0) % 360.0)


# ---------------------------------------------------------------------------
# Ellipsoid helpers
# ---------------------------------------------------------------------------

def geocentric_radius(lat_deg: ArrayLike) -> ArrayLike:
    """
    Return the WGS84 geocentric radius at a geodetic latitude.

    Parameters
    ----------
    lat_deg
        Geodetic latitude in degrees (positive north).

    Returns
    -------
    float or numpy.ndarray
        Distance from the Earth's center to the ellipsoid surface, in
        meters.
    """
    lat = np.asarray(lat_deg, dtype=float) * DEG2RAD
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    a2 = WGS84_A ** 2
    b2 = WGS84_B ** 2

    num = (a2 * cos_lat) ** 2 + (b2 * sin_lat) ** 2
    den = a2 * cos_lat ** 2 + b2 * sin_lat ** 2

    return _to_output(np.sqrt(num / den))


# ---------------------------------------------------------------------------
# Coordinate transforms (spherical Earth and sinusoidal projection)
# ---------------------------------------------------------------------------

def wgs84_to_cartesian_sphere(
    lon_deg: ArrayLike,
    lat_deg: ArrayLike,
    elevation_m: ArrayLike = 0.0,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert geographic coordinates to Cartesian (spherical Earth).

    Parameters
    ----------
    lon_deg, lat_deg
        Longitude and latitude in degrees.
    elevation_m
        Elevation above the reference sphere in meters.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    x, y, z : numpy.ndarray
        Cartesian coordinates in meters, in a right-handed system with
        origin at the sphere center, z toward the north pole and x
        through longitude 0 on the equator.
    """
    lon = np.asarray(lon_deg, dtype=float) * DEG2RAD
    lat = np.asarray(lat_deg, dtype=float) * DEG2RAD
    rho = radius_m + np.asarray(elevation_m, dtype=float)

    x = rho * np.cos(lat) * np.cos(lon)
    y = rho * np.cos(lat) * np.sin(lon)
    z = rho * np.sin(lat)

    return x, y, z


def cartesian_sphere_to_wgs84(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Cartesian coordinates (spherical Earth) to geographic.

    Parameters
    ----------
    x, y, z
        Cartesian coordinates in meters, as returned by
        :func:`wgs84_to_cartesian_sphere`.
    radius_m
        Sphere radius in meters, used as reference for the elevation.

    Returns
    -------
    lon_deg, lat_deg, elevation_m : numpy.ndarray
        Longitude in ``[-180, 180]`` and latitude in degrees, and
        elevation above the reference sphere in meters.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    z_arr = np.asarray(z, dtype=float)

    rho = np.sqrt(x_arr ** 2 + y_arr ** 2 + z_arr ** 2)
    lat = np.arctan2(z_arr, np.hypot(x_arr, y_arr))
    lon = np.arctan2(y_arr, x_arr)

    return lon * RAD2DEG, lat * RAD2DEG, rho - radius_m


def wgs84_to_sinusoidal(
    lon_deg: ArrayLike,
    lat_deg: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
    lon0_deg: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Project geographic coordinates with the sinusoidal projection.

    The projection is equal-area. Longitudes are taken relative to the
    central meridian ``lon0_deg`` and wrapped to ``[-180, 180)``, so the
    result does not depend on the longitude convention and is
    continuous across the antimeridian around ``lon0_deg``.

    Parameters
    ----------
    lon_deg, lat_deg
        Longitude and latitude in degrees.
    radius_m
        Sphere radius in meters.
    lon0_deg
        Central meridian in degrees.

    Returns
    -------
    x, y : numpy.ndarray
        Projected coordinates in meters.
    """
    dlon = np.asarray(relative_longitude(lon_deg, lon0_deg)) * DEG2RAD
    lat = np.asarray(lat_deg, dtype=float) * DEG2RAD

    x = radius_m * dlon * np.cos(lat)
    y = radius_m * lat

    return x, y


def sinusoidal_to_wgs84(
    x: ArrayLike,
    y: ArrayLike,
    radius_m: float = MEAN_EARTH_RADIUS_M,
    lon0_deg: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Invert the sinusoidal projection.

    Parameters
    ----------
    x, y
        Projected coordinates in meters.
    radius_m
        Sphere radius in meters used in the projection.
    lon0_deg
        Central meridian in degrees used in the projection.

    Returns
    -------
    lon_deg, lat_deg : numpy.ndarray
        Longitude (``lon0_deg`` plus an offset in ``[-180, 180]``, i.e.
        in the convention of ``lon0_deg``) and latitude in degrees. At
        the poles the longitude is undefined and ``lon0_deg`` is
        returned.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    lat = y_arr / radius_m
    cos_lat = np.cos(lat)

    with np.errstate(divide="ignore", invalid="ignore"):
        dlon = np.where(
            np.abs(cos_lat) > 1e-15,
            x_arr / (radius_m * cos_lat),
            0.0,
        )

    return lon0_deg + dlon * RAD2DEG, lat * RAD2DEG


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _to_output(value: Any) -> ArrayLike:
    """
    Return a float for 0-d results and a float array otherwise.
    """
    array = np.asarray(value, dtype=float)

    if array.ndim == 0:
        return float(array)

    return array


def _haversine_term(
    lon1_deg: ArrayLike,
    lat1_deg: ArrayLike,
    lon2_deg: ArrayLike,
    lat2_deg: ArrayLike,
) -> np.ndarray:
    """
    Return ``sin^2(theta / 2)`` for the central angle ``theta``.

    The value is clipped to ``[0, 1]`` against rounding errors.
    """
    lon1 = np.asarray(lon1_deg, dtype=float) * DEG2RAD
    lat1 = np.asarray(lat1_deg, dtype=float) * DEG2RAD
    lon2 = np.asarray(lon2_deg, dtype=float) * DEG2RAD
    lat2 = np.asarray(lat2_deg, dtype=float) * DEG2RAD

    sin_dlat = np.sin(0.5 * (lat2 - lat1))
    sin_dlon = np.sin(0.5 * (lon2 - lon1))

    hav = sin_dlat ** 2 + np.cos(lat1) * np.cos(lat2) * sin_dlon ** 2

    return np.clip(hav, 0.0, 1.0)


def _dms_string_to_decimal(value: str) -> float:
    """
    Parse a single DMS or decimal string.
    """
    text = value.strip().upper().replace(",", ".")

    if not text:
        raise ValueError("Empty DMS string.")

    hemisphere_sign = None

    if text[-1] in ("N", "S", "E", "W"):
        hemisphere_sign = -1.0 if text[-1] in ("S", "W") else 1.0
        text = text[:-1].strip()

    sign = 1.0

    if text.startswith("-"):
        sign = -1.0
        text = text[1:]
    elif text.startswith("+"):
        text = text[1:]

    parts = [part for part in re.split(r"[^\d.]+", text) if part]

    if not 1 <= len(parts) <= 3:
        raise ValueError(f"Cannot parse DMS string: {value!r}")

    try:
        numbers = [float(part) for part in parts]
    except ValueError:
        raise ValueError(f"Cannot parse DMS string: {value!r}") from None

    degrees = numbers[0]

    if len(numbers) > 1:
        degrees += numbers[1] / 60.0

    if len(numbers) > 2:
        degrees += numbers[2] / 3600.0

    if hemisphere_sign is not None:
        sign = hemisphere_sign

    return sign * degrees
