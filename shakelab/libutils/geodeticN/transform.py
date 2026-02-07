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
Coordinate transformation utilities for ShakeLab.

This module provides basic tools to convert between geodetic WGS84
coordinates (longitude, latitude, elevation) and local Cartesian
frames. Two types of local Cartesian frames are supported:

- **ECEF-shifted**: the point is converted to global ECEF coordinates
  and optionally translated so that a given reference point is the
  origin. Axes remain aligned with the global Earth-Centered
  Earth-Fixed frame.

- **ENU (East–North–Up)**: the point is expressed in a local
  topocentric frame defined around a reference geodetic point.
  Axes are oriented toward East, North, and Up respectively. This
  frame is particularly useful for local seismic, geodetic and
  engineering applications, as it cleanly separates horizontal and
  vertical components.

The functions :func:`wgs_to_metric` and :func:`metric_to_wgs` expose
these options through the optional ``frame`` argument.

Refactoring note
----------------
The local reference point is stored only in :class:`MetricFrame`.
Legacy APIs that passed the reference point separately have been
removed in this module revision.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple, Union


from shakelab.libutils.geodeticN.primitives import WgsPoint


# ============================================================================
# WGS84 ellipsoid parameters
# ============================================================================


WGS84_A: float = 6378137.0              # semi-major axis [m]
WGS84_F: float = 1.0 / 298.257223563    # flattening [-]
WGS84_E2: float = WGS84_F * (2.0 - WGS84_F)  # first ecc. squared

WGS84_B: float = WGS84_A * (1.0 - WGS84_F)   # semi-minor axis [m]
WGS84_EP2: float = (WGS84_A**2 - WGS84_B**2) / (WGS84_B**2)


# ============================================================================
# Basic metric containers
# ============================================================================


GeoLike = Union[WgsPoint, Tuple[float, float, float]]


@dataclass(frozen=True)
class MetricFrame:
    """
    Definition of a local Cartesian reference frame.

    A MetricFrame defines both the origin and the orientation of a
    local Cartesian system used to represent geodetic positions in
    meters.

    The reference point is normalized at initialization time, so
    ``ref_geo`` is always stored internally as a numeric triple
    ``(lon, lat, h)`` with lon/lat in degrees and h in meters.

    Attributes
    ----------
    ref_geo : tuple
        Reference geodetic point defining the origin of the local
        frame, stored as ``(lon, lat, h)`` with lon/lat in degrees and
        h in meters.
    orientation : {"ecef", "enu"}
        Orientation of the Cartesian axes:

        - ``"ecef"``: translated ECEF frame (axes aligned with global
          ECEF; origin at ``ref_geo``).
        - ``"enu"``: local East-North-Up frame defined at ``ref_geo``.
    """

    ref_geo: GeoLike
    orientation: str = "ecef"

    def __post_init__(self) -> None:
        lon, lat, h = _coerce_geo(self.ref_geo)
        object.__setattr__(self, "ref_geo", (lon, lat, h))

        orient = self.orientation.lower().strip()
        if orient not in ("ecef", "enu"):
            raise ValueError(
                "MetricFrame.orientation must be 'ecef' or 'enu'."
            )
        object.__setattr__(self, "orientation", orient)


@dataclass
class MetricPoint:
    """
    A point expressed in a local 3-D Cartesian frame.

    The coordinates (x, y, z) represent the position of a point in
    meters in a Cartesian system. The definition of the local
    frame (origin and axis orientation) can be carried explicitly via
    the optional :class:`MetricFrame` attribute.

    Attributes
    ----------
    x, y, z : float
        Cartesian coordinates in meters.
    frame : MetricFrame or None
        Optional frame definition associated with this point.
        If None, coordinates are interpreted as absolute ECEF.
    """

    x: float
    y: float
    z: float
    frame: Optional[MetricFrame] = None

    def to_wgs(self, frame=None) -> "WgsPoint":
        """
        Convert this MetricPoint to a WgsPoint.

        Parameters
        ----------
        frame : MetricFrame or None, optional
            Frame in which the metric coordinates are expressed.

            If None and this point carries a non-null ``self.frame``,
            that frame is used. Otherwise, the coordinates are
            interpreted as absolute ECEF.

        Returns
        -------
        WgsPoint
            Geodetic coordinates (longitude, latitude, elevation) in
            the WGS84 system.
        """
        # Local import to avoid circular dependencies (transform imports
        # WgsPoint from primitives).
        from shakelab.libutils.geodeticN.transform import metric_to_wgs

        return metric_to_wgs(self, frame=frame)


MetricLike = Union[MetricPoint, Tuple[float, float, float]]


# ============================================================================
# Internal coercion helpers
# ============================================================================


def _coerce_geo(geo: GeoLike) -> Tuple[float, float, float]:
    """
    Normalize a geodetic input to (lon, lat, h) floats.

    Parameters
    ----------
    geo : WgsPoint or tuple
        Either a :class:`WgsPoint` instance or a tuple
        ``(lon, lat, h)`` with lon/lat in degrees and h in meters.

    Returns
    -------
    lon, lat, h : float
        Normalized geodetic coordinates.
    """
    if isinstance(geo, WgsPoint):
        return float(geo.longitude), float(geo.latitude), float(geo.elevation)

    lon, lat, h = geo
    return float(lon), float(lat), float(h)


def _coerce_metric(metric: MetricLike) -> Tuple[float, float, float]:
    """
    Normalize a metric input to a 3-element tuple ``(x, y, z)``.

    Parameters
    ----------
    metric : MetricPoint or tuple
        Either a :class:`MetricPoint` instance or a tuple
        ``(x, y, z)`` in meters.

    Returns
    -------
    x, y, z : float
        Normalized Cartesian coordinates in meters.
    """
    if isinstance(metric, MetricPoint):
        return float(metric.x), float(metric.y), float(metric.z)

    x, y, z = metric
    return float(x), float(y), float(z)


def _coerce_frame(frame: Optional[MetricFrame]) -> Optional[MetricFrame]:
    """
    Validate a MetricFrame argument.

    Parameters
    ----------
    frame : MetricFrame or None
        Frame to validate. If None, coordinates are interpreted as
        absolute ECEF.

    Returns
    -------
    MetricFrame or None
        The validated frame.

    Raises
    ------
    TypeError
        If frame is not a MetricFrame and not None.
    """
    if frame is None or isinstance(frame, MetricFrame):
        return frame
    raise TypeError("frame must be a MetricFrame or None.")


# ============================================================================
# Core WGS84 <-> ECEF transforms (absolute)
# ============================================================================


def _geo_to_ecef_raw(
    lon: float,
    lat: float,
    h: float,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to absolute ECEF (WGS84).

    Parameters
    ----------
    lon, lat : float
        Longitude and latitude in degrees.
    h : float
        Ellipsoidal height above WGS84 in meters.

    Returns
    -------
    x, y, z : float
        Absolute ECEF coordinates in meters.
    """
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)

    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)

    x = (n + h) * cos_lat * math.cos(lon_rad)
    y = (n + h) * cos_lat * math.sin(lon_rad)
    z = (n * (1.0 - WGS84_E2) + h) * sin_lat

    return x, y, z


def geo_to_ecef(
    geo: GeoLike,
    frame: Optional[MetricFrame] = None,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to ECEF, with optional local shift.

    If ``frame`` is None, this function returns absolute ECEF
    coordinates. If ``frame`` is provided, its orientation must be
    ``"ecef"`` and the returned coordinates are expressed in a local
    translated ECEF system whose origin is ``frame.ref_geo``.

    Parameters
    ----------
    geo : WgsPoint or tuple
        Geodetic coordinates of the point to be converted.
    frame : MetricFrame or None, optional
        Local frame definition. Only ``orientation="ecef"`` is valid.

    Returns
    -------
    x, y, z : float
        ECEF coordinates in meters (absolute if frame is None, local
        translated if frame is provided).

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not "ecef".
    """
    frame = _coerce_frame(frame)

    lon, lat, h = _coerce_geo(geo)
    x, y, z = _geo_to_ecef_raw(lon, lat, h)

    if frame is None:
        return x, y, z

    if frame.orientation != "ecef":
        raise ValueError("geo_to_ecef requires orientation='ecef'.")

    lon0, lat0, h0 = _coerce_geo(frame.ref_geo)
    x0, y0, z0 = _geo_to_ecef_raw(lon0, lat0, h0)

    return x - x0, y - y0, z - z0


def ecef_to_geo(
    ecef: MetricLike,
    frame: Optional[MetricFrame] = None,
) -> Tuple[float, float, float]:
    """
    Convert ECEF coordinates to geodetic coordinates.

    If ``frame`` is None, the input is interpreted as absolute ECEF.
    If ``frame`` is provided, its orientation must be ``"ecef"`` and the
    input is interpreted as local translated ECEF coordinates relative
    to ``frame.ref_geo``.

    Parameters
    ----------
    ecef : MetricPoint or tuple
        ECEF coordinates to be converted.
    frame : MetricFrame or None, optional
        Local frame definition. Only ``orientation="ecef"`` is valid.

    Returns
    -------
    lon, lat, h : float
        Geodetic coordinates (longitude, latitude in degrees; height in
        meters).

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not "ecef".
    """
    frame = _coerce_frame(frame)

    x, y, z = _coerce_metric(ecef)

    if frame is not None:
        if frame.orientation != "ecef":
            raise ValueError("ecef_to_geo requires orientation='ecef'.")
        lon0, lat0, h0 = _coerce_geo(frame.ref_geo)
        x0, y0, z0 = _geo_to_ecef_raw(lon0, lat0, h0)
        x += x0
        y += y0
        z += z0

    p = math.hypot(x, y)
    if p == 0.0:
        lat_rad = math.copysign(math.pi / 2.0, z)
        lon_rad = 0.0
        h = abs(z) - WGS84_B
        return math.degrees(lon_rad), math.degrees(lat_rad), float(h)

    theta = math.atan2(z * WGS84_A, p * WGS84_B)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)

    lat_rad = math.atan2(
        z + WGS84_EP2 * WGS84_B * sin_theta**3,
        p - WGS84_E2 * WGS84_A * cos_theta**3,
    )
    lon_rad = math.atan2(y, x)

    sin_lat = math.sin(lat_rad)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    h = p / math.cos(lat_rad) - n

    return math.degrees(lon_rad), math.degrees(lat_rad), float(h)


# ============================================================================
# Core WGS84 <-> ENU transforms
# ============================================================================


def geo_to_enu(
    geo: GeoLike,
    frame: MetricFrame,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to local ENU coordinates.

    Parameters
    ----------
    geo : WgsPoint or tuple
        Geodetic coordinates of the point to be converted.
    frame : MetricFrame
        Local frame definition. Must have ``orientation="enu"``.

    Returns
    -------
    e, n, u : float
        Local East, North and Up coordinates in meters.

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not "enu".
    """
    if frame.orientation != "enu":
        raise ValueError("geo_to_enu requires orientation='enu'.")

    lon, lat, h = _coerce_geo(geo)
    lon0, lat0, h0 = _coerce_geo(frame.ref_geo)

    x, y, z = _geo_to_ecef_raw(lon, lat, h)
    x0, y0, z0 = _geo_to_ecef_raw(lon0, lat0, h0)

    dx = x - x0
    dy = y - y0
    dz = z - z0

    lat0_rad = math.radians(lat0)
    lon0_rad = math.radians(lon0)

    sin_lat0 = math.sin(lat0_rad)
    cos_lat0 = math.cos(lat0_rad)
    sin_lon0 = math.sin(lon0_rad)
    cos_lon0 = math.cos(lon0_rad)

    e = -sin_lon0 * dx + cos_lon0 * dy
    n = (
        -sin_lat0 * cos_lon0 * dx
        -sin_lat0 * sin_lon0 * dy
        + cos_lat0 * dz
    )
    u = (
        cos_lat0 * cos_lon0 * dx
        + cos_lat0 * sin_lon0 * dy
        + sin_lat0 * dz
    )

    return e, n, u


def enu_to_geo(
    enu: MetricLike,
    frame: MetricFrame,
) -> Tuple[float, float, float]:
    """
    Convert local ENU coordinates back to geodetic.

    Parameters
    ----------
    enu : MetricPoint or tuple
        Local ENU coordinates to be converted (E, N, U) in meters.
    frame : MetricFrame
        Local frame definition. Must have ``orientation="enu"``.

    Returns
    -------
    lon, lat, h : float
        Geodetic coordinates (longitude, latitude in degrees; height in
        meters).

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not "enu".
    """
    if frame.orientation != "enu":
        raise ValueError("enu_to_geo requires orientation='enu'.")

    e, n, u = _coerce_metric(enu)
    lon0, lat0, h0 = _coerce_geo(frame.ref_geo)

    x0, y0, z0 = _geo_to_ecef_raw(lon0, lat0, h0)

    lat0_rad = math.radians(lat0)
    lon0_rad = math.radians(lon0)

    sin_lat0 = math.sin(lat0_rad)
    cos_lat0 = math.cos(lat0_rad)
    sin_lon0 = math.sin(lon0_rad)
    cos_lon0 = math.cos(lon0_rad)

    dx = (
        -sin_lon0 * e
        -sin_lat0 * cos_lon0 * n
        + cos_lat0 * cos_lon0 * u
    )
    dy = (
        cos_lon0 * e
        -sin_lat0 * sin_lon0 * n
        + cos_lat0 * sin_lon0 * u
    )
    dz = cos_lat0 * n + sin_lat0 * u

    x = x0 + dx
    y = y0 + dy
    z = z0 + dz

    return ecef_to_geo((x, y, z), frame=None)


# ============================================================================
# High-level helpers using WgsPoint / MetricPoint
# ============================================================================


def wgs_to_metric(
    point: WgsPoint,
    frame: Optional[MetricFrame] = None,
) -> MetricPoint:
    """
    Convert a WgsPoint to a MetricPoint, in a chosen frame.

    Parameters
    ----------
    point : WgsPoint
        Geodetic coordinates to be converted.
    frame : MetricFrame or None, optional
        Local frame definition. If None, absolute ECEF coordinates are
        returned.

    Returns
    -------
    MetricPoint
        Cartesian coordinates in meters, in the requested frame.
    """
    frame = _coerce_frame(frame)

    if frame is None:
        x, y, z = geo_to_ecef(point, frame=None)
        return MetricPoint(float(x), float(y), float(z), frame=None)

    if frame.orientation == "ecef":
        x, y, z = geo_to_ecef(point, frame=frame)
        return MetricPoint(float(x), float(y), float(z), frame=frame)

    if frame.orientation == "enu":
        e, n, u = geo_to_enu(point, frame=frame)
        return MetricPoint(float(e), float(n), float(u), frame=frame)

    raise ValueError(
        f"Unknown frame '{frame.orientation}' (expected 'ecef' or 'enu')."
    )


def metric_to_wgs(
    point: MetricPoint,
    frame: Optional[MetricFrame] = None,
) -> WgsPoint:
    """
    Convert a MetricPoint to a WgsPoint, in a chosen frame.

    If ``frame`` is None and the input point carries a non-null
    ``point.frame``, the point frame is used.

    Parameters
    ----------
    point : MetricPoint
        Cartesian coordinates to be converted.
    frame : MetricFrame or None, optional
        Frame in which the point coordinates are expressed. If None,
        the point's own frame is used (if present); otherwise the
        coordinates are interpreted as absolute ECEF.

    Returns
    -------
    WgsPoint
        Geodetic coordinates (lon, lat, elev) in the WGS84 system.

    Raises
    ------
    ValueError
        If the effective frame is inconsistent.
    """
    if frame is None:
        frame = point.frame
    frame = _coerce_frame(frame)

    if frame is not None and point.frame is not None:
        if point.frame != frame:
            raise ValueError(
                "metric_to_wgs: point.frame is incompatible with the "
                "provided frame."
            )

    if frame is None:
        lon, lat, h = ecef_to_geo(point, frame=None)
        return WgsPoint(
            longitude=float(lon),
            latitude=float(lat),
            elevation=float(h),
        )

    if frame.orientation == "ecef":
        lon, lat, h = ecef_to_geo(point, frame=frame)
    elif frame.orientation == "enu":
        lon, lat, h = enu_to_geo(point, frame=frame)
    else:
        raise ValueError(
            f"Unknown frame '{frame.orientation}' (expected 'ecef' or 'enu')."
        )

    return WgsPoint(
        longitude=float(lon),
        latitude=float(lat),
        elevation=float(h),
    )
