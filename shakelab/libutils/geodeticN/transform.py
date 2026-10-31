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
Coordinate transformations between geodetic and Cartesian frames.

This module converts geodetic WGS84 coordinates (longitude, latitude,
ellipsoidal height) to and from:

- **ECEF** (Earth-Centered, Earth-Fixed) coordinates, absolute or
  translated so that a reference point is the origin;
- **ENU** (East-North-Up) local topocentric coordinates around a
  reference point, convenient for local seismic and engineering work.

The high-level helpers :func:`wgs_to_metric` and :func:`metric_to_wgs`
select ECEF or ENU through an optional :class:`MetricFrame`.

Design choices
--------------
- The local reference point is stored only in :class:`MetricFrame`.
- Local metric coordinates are obtained through :class:`MetricFrame`
  (translated ECEF or ENU); map projections such as UTM are
  intentionally not provided.
- ECEF to geodetic uses Bowring's closed form for the latitude and a
  height formula that remains stable near the poles.

Conventions
-----------
- Geodetic triples are ``(lon, lat, h)``: longitude and latitude in
  degrees, ellipsoidal height in meters (positive upward). Two-element
  tuples ``(lon, lat)`` are accepted with ``h = 0``.
- Longitudes may use the ``[-180, 180)`` or ``[0, 360)`` convention;
  returned longitudes are in ``[-180, 180]``.
- Cartesian coordinates are in meters.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Sequence, Tuple, Union

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.utilities import (
    WGS84_A,
    WGS84_B,
    WGS84_E2,
    WGS84_EP2,
    WGS84_F,
)


__all__ = [
    "WGS84_A",
    "WGS84_F",
    "WGS84_B",
    "WGS84_E2",
    "WGS84_EP2",
    "MetricFrame",
    "MetricPoint",
    "geo_to_ecef",
    "ecef_to_geo",
    "geo_to_enu",
    "enu_to_geo",
    "wgs_to_metric",
    "metric_to_wgs",
]


GeoLike = Union[WgsPoint, Sequence[float]]


# ---------------------------------------------------------------------------
# Metric containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MetricFrame:
    """
    Definition of a local Cartesian reference frame.

    Attributes
    ----------
    ref_geo
        Reference geodetic point defining the origin of the frame: a
        :class:`WgsPoint` or a ``(lon, lat[, h])`` sequence. It is
        normalized at initialization to a ``(lon, lat, h)`` tuple of
        floats.
    orientation
        Orientation of the axes: ``"ecef"`` (translated ECEF, axes
        aligned with global ECEF) or ``"enu"`` (local East-North-Up).

    Raises
    ------
    TypeError
        If ``orientation`` is not a string.
    ValueError
        If ``orientation`` is not ``"ecef"`` or ``"enu"``.
    """

    ref_geo: GeoLike
    orientation: str = "ecef"

    def __post_init__(self) -> None:
        object.__setattr__(self, "ref_geo", _coerce_geo(self.ref_geo))

        if not isinstance(self.orientation, str):
            raise TypeError("MetricFrame.orientation must be a string.")

        orientation = self.orientation.lower().strip()

        if orientation not in ("ecef", "enu"):
            raise ValueError(
                "MetricFrame.orientation must be 'ecef' or 'enu'."
            )

        object.__setattr__(self, "orientation", orientation)


@dataclass
class MetricPoint:
    """
    Point in a 3-D Cartesian frame.

    Attributes
    ----------
    x, y, z
        Cartesian coordinates in meters.
    frame
        Frame in which the coordinates are expressed. ``None`` means
        absolute ECEF.
    """

    x: float
    y: float
    z: float
    frame: Optional[MetricFrame] = None

    def to_wgs(self, frame: Optional[MetricFrame] = None) -> WgsPoint:
        """
        Convert the point to geodetic coordinates.

        Parameters
        ----------
        frame
            Frame of the coordinates. If omitted, :attr:`frame` is used.

        Returns
        -------
        WgsPoint
            Geodetic point.
        """
        return metric_to_wgs(self, frame=frame)


MetricLike = Union[MetricPoint, Sequence[float]]


# ---------------------------------------------------------------------------
# Absolute and translated ECEF
# ---------------------------------------------------------------------------

def geo_to_ecef(
    geo: GeoLike,
    frame: Optional[MetricFrame] = None,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to ECEF, optionally translated.

    Parameters
    ----------
    geo
        Geodetic point: :class:`WgsPoint` or ``(lon, lat[, h])``.
    frame
        If given, it must have ``orientation="ecef"`` and the result is
        translated so that ``frame.ref_geo`` is the origin.

    Returns
    -------
    x, y, z : float
        ECEF coordinates in meters.

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not ``"ecef"``.
    """
    frame = _coerce_frame(frame)
    x, y, z = _geo_to_ecef_raw(*_coerce_geo(geo))

    if frame is None:
        return x, y, z

    if frame.orientation != "ecef":
        raise ValueError("geo_to_ecef requires orientation='ecef'.")

    x0, y0, z0 = _geo_to_ecef_raw(*frame.ref_geo)

    return x - x0, y - y0, z - z0


def ecef_to_geo(
    ecef: MetricLike,
    frame: Optional[MetricFrame] = None,
) -> Tuple[float, float, float]:
    """
    Convert ECEF coordinates to geodetic coordinates.

    Parameters
    ----------
    ecef
        ECEF coordinates: :class:`MetricPoint` or ``(x, y, z)``.
    frame
        If given, it must have ``orientation="ecef"`` and the input is
        interpreted relative to ``frame.ref_geo``.

    Returns
    -------
    lon, lat, h : float
        Longitude and latitude in degrees, ellipsoidal height in meters.

    Raises
    ------
    ValueError
        If ``frame.orientation`` is not ``"ecef"``.
    """
    frame = _coerce_frame(frame)
    x, y, z = _coerce_metric(ecef)

    if frame is not None:
        if frame.orientation != "ecef":
            raise ValueError("ecef_to_geo requires orientation='ecef'.")

        x0, y0, z0 = _geo_to_ecef_raw(*frame.ref_geo)
        x += x0
        y += y0
        z += z0

    p = math.hypot(x, y)

    if p == 0.0:
        lat = math.copysign(math.pi / 2.0, z)
        return 0.0, math.degrees(lat), abs(z) - WGS84_B

    # Bowring's closed-form latitude.
    theta = math.atan2(z * WGS84_A, p * WGS84_B)
    lat = math.atan2(
        z + WGS84_EP2 * WGS84_B * math.sin(theta) ** 3,
        p - WGS84_E2 * WGS84_A * math.cos(theta) ** 3,
    )
    lon = math.atan2(y, x)

    # Height formula stable at all latitudes (no division by cos(lat)).
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    h = (
        p * cos_lat
        + z * sin_lat
        - WGS84_A * math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    )

    return math.degrees(lon), math.degrees(lat), h


# ---------------------------------------------------------------------------
# Local East-North-Up
# ---------------------------------------------------------------------------

def geo_to_enu(
    geo: GeoLike,
    frame: MetricFrame,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to local ENU coordinates.

    Parameters
    ----------
    geo
        Geodetic point: :class:`WgsPoint` or ``(lon, lat[, h])``.
    frame
        Local frame with ``orientation="enu"``.

    Returns
    -------
    e, n, u : float
        East, North and Up coordinates in meters.

    Raises
    ------
    TypeError
        If ``frame`` is not a :class:`MetricFrame`.
    ValueError
        If ``frame.orientation`` is not ``"enu"``.
    """
    frame = _require_enu_frame(frame, "geo_to_enu")

    x, y, z = _geo_to_ecef_raw(*_coerce_geo(geo))
    x0, y0, z0 = _geo_to_ecef_raw(*frame.ref_geo)

    dx = x - x0
    dy = y - y0
    dz = z - z0

    sin_lat0, cos_lat0, sin_lon0, cos_lon0 = _frame_trig(frame)

    e = -sin_lon0 * dx + cos_lon0 * dy
    n = (
        -sin_lat0 * cos_lon0 * dx
        - sin_lat0 * sin_lon0 * dy
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
    Convert local ENU coordinates to geodetic coordinates.

    Parameters
    ----------
    enu
        ENU coordinates: :class:`MetricPoint` or ``(e, n, u)``.
    frame
        Local frame with ``orientation="enu"``.

    Returns
    -------
    lon, lat, h : float
        Longitude and latitude in degrees, ellipsoidal height in meters.

    Raises
    ------
    TypeError
        If ``frame`` is not a :class:`MetricFrame`.
    ValueError
        If ``frame.orientation`` is not ``"enu"``.
    """
    frame = _require_enu_frame(frame, "enu_to_geo")

    e, n, u = _coerce_metric(enu)
    x0, y0, z0 = _geo_to_ecef_raw(*frame.ref_geo)

    sin_lat0, cos_lat0, sin_lon0, cos_lon0 = _frame_trig(frame)

    dx = (
        -sin_lon0 * e
        - sin_lat0 * cos_lon0 * n
        + cos_lat0 * cos_lon0 * u
    )
    dy = (
        cos_lon0 * e
        - sin_lat0 * sin_lon0 * n
        + cos_lat0 * sin_lon0 * u
    )
    dz = cos_lat0 * n + sin_lat0 * u

    return ecef_to_geo((x0 + dx, y0 + dy, z0 + dz), frame=None)


# ---------------------------------------------------------------------------
# High-level helpers
# ---------------------------------------------------------------------------

def wgs_to_metric(
    point: GeoLike,
    frame: Optional[MetricFrame] = None,
) -> MetricPoint:
    """
    Convert a geodetic point to a :class:`MetricPoint`.

    Parameters
    ----------
    point
        Geodetic point: :class:`WgsPoint` or ``(lon, lat[, h])``.
    frame
        Output frame. ``None`` returns absolute ECEF coordinates.

    Returns
    -------
    MetricPoint
        Cartesian coordinates in meters, carrying ``frame``.
    """
    frame = _coerce_frame(frame)

    if frame is None or frame.orientation == "ecef":
        x, y, z = geo_to_ecef(point, frame=frame)
    else:
        x, y, z = geo_to_enu(point, frame=frame)

    return MetricPoint(float(x), float(y), float(z), frame=frame)


def metric_to_wgs(
    point: MetricPoint,
    frame: Optional[MetricFrame] = None,
) -> WgsPoint:
    """
    Convert a :class:`MetricPoint` to a :class:`WgsPoint`.

    Parameters
    ----------
    point
        Cartesian point.
    frame
        Frame of the coordinates. If omitted, ``point.frame`` is used;
        if that is also ``None``, absolute ECEF is assumed.

    Returns
    -------
    WgsPoint
        Geodetic point.

    Raises
    ------
    ValueError
        If ``frame`` and ``point.frame`` are both given and differ.
    """
    if frame is None:
        frame = point.frame

    frame = _coerce_frame(frame)

    if (
        frame is not None
        and point.frame is not None
        and point.frame != frame
    ):
        raise ValueError(
            "metric_to_wgs: point.frame is incompatible with the "
            "provided frame."
        )

    if frame is None or frame.orientation == "ecef":
        lon, lat, h = ecef_to_geo(point, frame=frame)
    else:
        lon, lat, h = enu_to_geo(point, frame=frame)

    return WgsPoint(longitude=lon, latitude=lat, elevation=h)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _coerce_geo(geo: GeoLike) -> Tuple[float, float, float]:
    """
    Normalize a geodetic input to a ``(lon, lat, h)`` tuple of floats.
    """
    if isinstance(geo, WgsPoint):
        return (
            float(geo.longitude),
            float(geo.latitude),
            float(geo.elevation),
        )

    values = tuple(geo)

    if len(values) == 2:
        return float(values[0]), float(values[1]), 0.0

    if len(values) == 3:
        return float(values[0]), float(values[1]), float(values[2])

    raise ValueError("Geodetic tuples must be (lon, lat) or (lon, lat, h).")


def _coerce_metric(metric: MetricLike) -> Tuple[float, float, float]:
    """
    Normalize a metric input to an ``(x, y, z)`` tuple of floats.
    """
    if isinstance(metric, MetricPoint):
        return float(metric.x), float(metric.y), float(metric.z)

    x, y, z = metric

    return float(x), float(y), float(z)


def _coerce_frame(frame: Optional[MetricFrame]) -> Optional[MetricFrame]:
    """
    Check that ``frame`` is a :class:`MetricFrame` or None.
    """
    if frame is None or isinstance(frame, MetricFrame):
        return frame

    raise TypeError("frame must be a MetricFrame or None.")


def _require_enu_frame(frame: MetricFrame, name: str) -> MetricFrame:
    """
    Check that ``frame`` is an ENU :class:`MetricFrame`.
    """
    if not isinstance(frame, MetricFrame):
        raise TypeError(f"{name} requires a MetricFrame.")

    if frame.orientation != "enu":
        raise ValueError(f"{name} requires orientation='enu'.")

    return frame


def _frame_trig(frame: MetricFrame) -> Tuple[float, float, float, float]:
    """
    Return sin/cos of the latitude and longitude of the frame origin.
    """
    lon0, lat0, _ = frame.ref_geo
    lat0 = math.radians(lat0)
    lon0 = math.radians(lon0)

    return math.sin(lat0), math.cos(lat0), math.sin(lon0), math.cos(lon0)


def _geo_to_ecef_raw(
    lon: float,
    lat: float,
    h: float,
) -> Tuple[float, float, float]:
    """
    Convert geodetic coordinates to absolute ECEF (WGS84).
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
