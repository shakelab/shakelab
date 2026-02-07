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
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Distance helpers for geodetic coordinates.

This module provides high-level wrappers to compute seismology-oriented
distances between geographic objects, consistently with the ShakeLab
geodetic conventions.

The main helpers are:

* :func:`epicentral_distance`:
  great-circle (surface) distance between event and station;
* :func:`hypocentral_distance`:
  thin wrapper around the tunnel-distance functions;
* :func:`tunnel_distance_sphere`:
  straight-line 3-D distance on a spherical Earth (with elevation);
* :func:`tunnel_distance_ellipsoid`:
  straight-line 3-D distance on the WGS84 ellipsoid (with elevation).

All geodetic coordinates follow the convention ``(lon, lat, h)``, i.e.
longitude first, then latitude, then ellipsoidal height in meters.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import geo_to_ecef
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
    wgs84_to_cartesian_sphere,
)


# ============================================================================
# Internal helpers
# ============================================================================


def _as_wgs_point(obj: Any) -> WgsPoint:
    """
    Convert a generic object with attributes
    (longitude, latitude, elevation) into a WgsPoint.

    Parameters
    ----------
    obj : Any
        An object carrying geodetic information. Must expose
        attributes 'longitude' and 'latitude'. Elevation defaults to 0.

    Returns
    -------
    WgsPoint
        A WgsPoint instance with the canonical ShakeLab attributes.
    """
    try:
        lon = obj.longitude
        lat = obj.latitude
    except AttributeError:
        raise TypeError(
            f"Object of type {type(obj)} does not expose "
            "'longitude' and 'latitude' attributes."
        )

    elev = getattr(obj, "elevation", 0.0)

    return WgsPoint(
        longitude=float(lon),
        latitude=float(lat),
        elevation=float(elev),
    )


# ============================================================================
# Seismological distance wrapper functions
# ============================================================================


def epicentral_distance(
    hypocenter: Any,
    station: Any,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute epicentral distance between two geographic points.

    This function returns the surface (great-circle) distance between
    two points defined by longitude and latitude. It accepts any
    objects that provide ``.lon`` and ``.lat`` attributes, such as
    :class:`WgsPoint`, without requiring a specific class type.

    Parameters
    ----------
    hypocenter : object
        Object with ``lon`` and ``lat`` attributes representing the
        event location.
    station : object
        Object with ``lon`` and ``lat`` attributes representing the
        station location.
    radius_m : float, optional
        Sphere radius in meters. Default is the mean Earth radius.

    Returns
    -------
    float
        Epicentral distance in meters along the spherical surface.
    """
    return great_circle_distance(
        hypocenter.longitude,
        hypocenter.latitude,
        station.longitude,
        station.latitude,
        radius_m=radius_m,
    )


def hypocentral_distance(
    hypocenter: Any,
    station: Any,
    approx: str = "ellipsoid",
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute hypocentral distance between an event and a station.

    Parameters
    ----------
    hypocenter : object
        Object convertible to :class:`WgsPoint` via ``_as_wgs_point``.
    station : object
        Object convertible to :class:`WgsPoint` via ``_as_wgs_point``.
    approx : str, optional
        Method to use for the 3-D distance. Valid options:
        * "ellipsoid"  → use :func:`tunnel_distance_ellipsoid`
        * "sphere"     → use :func:`tunnel_distance_sphere`
        Default is "ellipsoid".
    radius_m : float, optional
        Sphere radius in meters used when ``approx="sphere"``.
        Default is :data:`MEAN_EARTH_RADIUS_M`.

    Returns
    -------
    float
        Hypocentral distance in meters.

    Raises
    ------
    ValueError
        If ``approx`` is not one of {"ellipsoid", "sphere"}.
    """
    hyp_pt = _as_wgs_point(hypocenter)
    sta_pt = _as_wgs_point(station)

    approx = approx.lower().strip()

    if approx == "ellipsoid":
        return tunnel_distance_ellipsoid(hyp_pt, sta_pt)

    if approx == "sphere":
        return tunnel_distance_sphere(hyp_pt, sta_pt, radius_m=radius_m)

    raise ValueError(
        "Invalid 'approx' option. Expected 'ellipsoid' or 'sphere', "
        f"got: {approx!r}"
    )


# ============================================================================
# Metric distance wrapper functions
# ============================================================================


def tunnel_distance_sphere(
    point1: WgsPoint,
    point2: WgsPoint,
    radius_m: float = MEAN_EARTH_RADIUS_M,
) -> float:
    """
    Compute 3-D straight-line distance between two points using a
    spherical Earth with elevation included.

    Parameters
    ----------
    point1, point2 : WgsPoint
        Geographic points, with ``lon``, ``lat`` and ``elev``.
    radius_m : float, optional
        Sphere radius in meters. Default is the mean Earth radius.

    Returns
    -------
    float
        Euclidean distance in meters (includes elevation).
    """
    x1, y1, z1 = wgs84_to_cartesian_sphere(
        point1.longitude,
        point1.latitude,
        point1.elevation,
        radius_m,
    )
    x2, y2, z2 = wgs84_to_cartesian_sphere(
        point2.longitude,
        point2.latitude,
        point2.elevation,
        radius_m,
    )

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    return float(math.sqrt(dx * dx + dy * dy + dz * dz))


def tunnel_distance_ellipsoid(
    point1: WgsPoint,
    point2: WgsPoint,
) -> float:
    """
    Compute 3-D straight-line distance between two points on the WGS84
    ellipsoid, including elevation.

    Parameters
    ----------
    point1, point2 : WgsPoint
        Points with ``lon``, ``lat`` and ``elev``.

    Returns
    -------
    float
        Euclidean distance in meters.
    """
    x1, y1, z1 = geo_to_ecef(point1, frame=None)
    x2, y2, z2 = geo_to_ecef(point2, frame=None)

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    return float(math.sqrt(dx * dx + dy * dy + dz * dz))

