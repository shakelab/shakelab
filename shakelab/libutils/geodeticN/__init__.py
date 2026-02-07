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
ShakeLab Geodetic Toolkit (geodeticN)
====================================

Public API initialization.

This module re-exports selected classes and functions from the
submodules (utilities, distance, primitives, mesh, sampling,
transform) so they can be imported directly from
:mod:`shakelab.libutils.geodeticN`.

Example
-------
>>> from shakelab.libutils.geodeticN import WgsPoint, epicentral_distance
>>> p1 = WgsPoint(longitude=13.0, latitude=46.0)
>>> p2 = WgsPoint(longitude=14.0, latitude=46.0)
>>> d = epicentral_distance(p1, p2)
"""

from .utilities import (
    MEAN_EARTH_RADIUS_M,
    EQUATORIAL_EARTH_RADIUS_M,
    POLAR_EARTH_RADIUS_M,
    DEG2RAD,
    RAD2DEG,
    deg_to_rad,
    rad_to_deg,
    normalize_longitude,
    great_circle_distance,
    chord_distance,
)

from .distance import (
    epicentral_distance,
    hypocentral_distance,
)

from .primitives import (
    WgsPoint,
    WgsPolygon,
    WgsMesh,
    PolygonKind,
)

from .mesh import (
    mesh_fibonacci_global,
    mesh_fibonacci_spherical_cap,
    mesh_fibonacci_bbox,
    mesh_cartesian,
    mesh_sinusoidal,
)

from .sampling import (
    sample_uniform_in_polygon_xy,
    sample_uniform_in_polygon_wgs84,
)

from .transform import (
    MetricFrame,
    MetricPoint,
    geo_to_ecef,
    ecef_to_geo,
    geo_to_enu,
    enu_to_geo,
    wgs_to_metric,
    metric_to_wgs,
)


__all__ = [
    # utilities / core constants
    "MEAN_EARTH_RADIUS_M",
    "EQUATORIAL_EARTH_RADIUS_M",
    "POLAR_EARTH_RADIUS_M",
    "DEG2RAD",
    "RAD2DEG",
    "deg_to_rad",
    "rad_to_deg",
    "normalize_longitude",
    "great_circle_distance",
    "chord_distance",
    # distance helpers
    "epicentral_distance",
    "hypocentral_distance",
    # geometric primitives
    "WgsPoint",
    "WgsPolygon",
    "WgsMesh",
    "PolygonKind",
    # mesh generators
    "mesh_fibonacci_global",
    "mesh_fibonacci_spherical_cap",
    "mesh_fibonacci_bbox",
    "mesh_cartesian",
    "mesh_sinusoidal",
    # sampling utilities
    "sample_uniform_in_polygon_xy",
    "sample_uniform_in_polygon_wgs84",
    # coordinate transforms
    "MetricFrame",
    "MetricPoint",
    "geo_to_ecef",
    "ecef_to_geo",
    "geo_to_enu",
    "enu_to_geo",
    "wgs_to_metric",
    "metric_to_wgs",
]
