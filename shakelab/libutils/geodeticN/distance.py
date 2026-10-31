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
Seismological distances between geographic points.

This module provides the distances used in seismology between an event
and a site:

- :func:`epicentral_distance`: great-circle (surface) distance;
- :func:`hypocentral_distance`: 3-D straight-line distance including
  elevation, on the WGS84 ellipsoid or on a sphere;
- :func:`tunnel_distance_sphere` and :func:`tunnel_distance_ellipsoid`:
  the two 3-D straight-line distances used by
  :func:`hypocentral_distance`.

Design choices
--------------
- The functions accept any object exposing ``longitude`` and
  ``latitude`` (and optionally ``elevation``), such as
  :class:`WgsPoint`.
- The epicentral distance is computed on a sphere of mean radius, while
  the hypocentral distance uses the WGS84 ellipsoid by default. The two
  models differ by a few tenths of a percent; use
  ``hypocentral_distance(..., approx="sphere")`` for a distance
  consistent with the epicentral one.

Conventions
-----------
- Longitude and latitude in degrees, in any longitude convention.
- Elevation in meters, positive upward; a hypocentre at 10 km depth has
  ``elevation = -10000``.
- Distances in meters.
"""

from __future__ import annotations

from typing import Any

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import geo_to_ecef
from shakelab.libutils.geodeticN import utilities


__all__ = [
    "epicentral_distance",
    "hypocentral_distance",
    "tunnel_distance_sphere",
    "tunnel_distance_ellipsoid",
]


# ---------------------------------------------------------------------------
# Public distances
# ---------------------------------------------------------------------------

def epicentral_distance(
    hypocenter: Any,
    station: Any,
    radius_m: float = utilities.MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute the epicentral (great-circle) distance.

    Elevations are ignored.

    Parameters
    ----------
    hypocenter
        Event location: an object with ``longitude`` and ``latitude``.
    station
        Site location: an object with ``longitude`` and ``latitude``.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    float
        Epicentral distance in meters.

    Raises
    ------
    TypeError
        If an input does not expose ``longitude`` and ``latitude``.
    """
    hyp = _as_wgs_point(hypocenter)
    sta = _as_wgs_point(station)

    return utilities.great_circle_distance(
        hyp.longitude,
        hyp.latitude,
        sta.longitude,
        sta.latitude,
        radius_m=radius_m,
    )


def hypocentral_distance(
    hypocenter: Any,
    station: Any,
    approx: str = "ellipsoid",
    radius_m: float = utilities.MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute the hypocentral (3-D straight-line) distance.

    Parameters
    ----------
    hypocenter
        Event location: an object with ``longitude``, ``latitude`` and
        optionally ``elevation`` (default 0).
    station
        Site location, as ``hypocenter``.
    approx
        Earth model: ``"ellipsoid"`` (WGS84, default) or ``"sphere"``.
    radius_m
        Sphere radius in meters, used when ``approx="sphere"``.

    Returns
    -------
    float
        Hypocentral distance in meters.

    Raises
    ------
    TypeError
        If an input does not expose ``longitude`` and ``latitude``.
    ValueError
        If ``approx`` is not ``"ellipsoid"`` or ``"sphere"``.
    """
    hyp = _as_wgs_point(hypocenter)
    sta = _as_wgs_point(station)
    model = str(approx).lower().strip()

    if model == "ellipsoid":
        return tunnel_distance_ellipsoid(hyp, sta)

    if model == "sphere":
        return tunnel_distance_sphere(hyp, sta, radius_m=radius_m)

    raise ValueError(
        "Invalid 'approx' option. Expected 'ellipsoid' or 'sphere', "
        f"got: {approx!r}"
    )


def tunnel_distance_sphere(
    point1: Any,
    point2: Any,
    radius_m: float = utilities.MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute the 3-D straight-line distance on a spherical Earth.

    This is the point-based counterpart of
    :func:`shakelab.libutils.geodeticN.utilities.tunnel_distance_sphere`.

    Parameters
    ----------
    point1, point2
        Objects with ``longitude``, ``latitude`` and optionally
        ``elevation``.
    radius_m
        Sphere radius in meters.

    Returns
    -------
    float
        Euclidean distance in meters, including elevation.
    """
    p1 = _as_wgs_point(point1)
    p2 = _as_wgs_point(point2)

    return utilities.tunnel_distance_sphere(
        p1.longitude,
        p1.latitude,
        p2.longitude,
        p2.latitude,
        elev1_m=p1.elevation,
        elev2_m=p2.elevation,
        radius_m=radius_m,
    )


def tunnel_distance_ellipsoid(point1: Any, point2: Any) -> float:
    """
    Compute the 3-D straight-line distance on the WGS84 ellipsoid.

    Parameters
    ----------
    point1, point2
        Objects with ``longitude``, ``latitude`` and optionally
        ``elevation``.

    Returns
    -------
    float
        Euclidean distance in meters, including elevation.
    """
    x1, y1, z1 = geo_to_ecef(_as_wgs_point(point1))
    x2, y2, z2 = geo_to_ecef(_as_wgs_point(point2))

    return ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _as_wgs_point(obj: Any) -> WgsPoint:
    """
    Convert an object with ``longitude``/``latitude`` to a WgsPoint.

    ``elevation`` defaults to 0 when missing.
    """
    if isinstance(obj, WgsPoint):
        return obj

    try:
        lon = obj.longitude
        lat = obj.latitude
    except AttributeError:
        raise TypeError(
            f"Object of type {type(obj).__name__} does not expose "
            "'longitude' and 'latitude' attributes."
        ) from None

    return WgsPoint(
        longitude=float(lon),
        latitude=float(lat),
        elevation=float(getattr(obj, "elevation", 0.0)),
    )
