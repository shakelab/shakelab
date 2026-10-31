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
ShakeLab geodetic toolkit.

The package collects geodetic constants and functions
(:mod:`.utilities`), geometric primitives (:mod:`.primitives`),
seismological distances (:mod:`.distance`), coordinate transforms
(:mod:`.transform`), mesh generators (:mod:`.mesh`), random sampling
(:mod:`.sampling`) and GeoJSON/CSV input and output (:mod:`.gio`).

The most used names are re-exported here, e.g.::

    >>> from shakelab.libutils.geodeticN import WgsPoint, epicentral_distance
    >>> p1 = WgsPoint(longitude=13.0, latitude=46.0)
    >>> p2 = WgsPoint(longitude=14.0, latitude=46.0)
    >>> distance = epicentral_distance(p1, p2)

Conventions
-----------
- Geographic coordinates are ``(longitude, latitude)`` in degrees, in
  the ``[-180, 180)`` or ``[0, 360)`` convention.
- Elevations and distances are in meters; elevation is positive upward.
"""

from shakelab.libutils.geodeticN.utilities import (
    DEG2RAD,
    EQUATORIAL_EARTH_RADIUS_M,
    MEAN_EARTH_RADIUS_M,
    POLAR_EARTH_RADIUS_M,
    RAD2DEG,
    chord_distance,
    deg_to_rad,
    dms2dec,
    forward_azimuth,
    great_circle_distance,
    normalize_longitude,
    rad_to_deg,
    relative_longitude,
)
from shakelab.libutils.geodeticN.primitives import (
    PolygonKind,
    WgsMesh,
    WgsPoint,
    WgsPolygon,
)
from shakelab.libutils.geodeticN.distance import (
    epicentral_distance,
    hypocentral_distance,
)
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    MetricPoint,
    ecef_to_geo,
    enu_to_geo,
    geo_to_ecef,
    geo_to_enu,
    metric_to_wgs,
    wgs_to_metric,
)
from shakelab.libutils.geodeticN.mesh import (
    mesh_cartesian,
    mesh_fibonacci_bbox,
    mesh_fibonacci_global,
    mesh_fibonacci_spherical_cap,
    mesh_sinusoidal,
)
from shakelab.libutils.geodeticN.sampling import (
    sample_uniform_in_polygon_wgs84,
    sample_uniform_in_polygon_xy,
)
from shakelab.libutils.geodeticN.gio import (
    read_csv,
    read_geojson,
    write_csv,
    write_geojson,
)


__all__ = [
    # utilities
    "MEAN_EARTH_RADIUS_M",
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
    "forward_azimuth",
    # primitives
    "WgsPoint",
    "WgsPolygon",
    "WgsMesh",
    "PolygonKind",
    # distances
    "epicentral_distance",
    "hypocentral_distance",
    # transforms
    "MetricFrame",
    "MetricPoint",
    "geo_to_ecef",
    "ecef_to_geo",
    "geo_to_enu",
    "enu_to_geo",
    "wgs_to_metric",
    "metric_to_wgs",
    # meshes
    "mesh_fibonacci_global",
    "mesh_fibonacci_spherical_cap",
    "mesh_fibonacci_bbox",
    "mesh_cartesian",
    "mesh_sinusoidal",
    # sampling
    "sample_uniform_in_polygon_xy",
    "sample_uniform_in_polygon_wgs84",
    # input and output
    "read_geojson",
    "write_geojson",
    "read_csv",
    "write_csv",
]
