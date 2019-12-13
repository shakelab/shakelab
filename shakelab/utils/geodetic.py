# ============================================================================
#
# Copyright (C) 2019, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ============================================================================

import numpy as np

# CONSTANTS
MEAN_EARTH_RADIUS = 6371008.8
EQUATORIAL_EARTH_RADIUS = 6378137.0
POLAR_EARTH_RADIUS = 6356752.3


def wgs_to_xyz_sphere(lat, lon, ele):
    """
    Convert WGS84 coordinates to cartesian using
    spherical earth approximation.
    """

    phi = np.radians(90. - lat)
    theta = np.radians(lon)
    rho = ele + MEAN_EARTH_RADIUS

    x = rho * np.sin(phi) * np.cos(theta)
    y = rho * np.sin(phi) * np.sin(theta)
    z = rho * np.cos(phi)

    return round(x, 3), round(y, 3), round(z, 3)

def wgs_to_xyz_ellipsoid(lat, lon, ele):
    """
    Convert WGS84 coordinates to cartesian using
    ellipsoidal earth approximation.
    """

    phi = np.radians(lat)
    theta = np.radians(lon)

    eer2 = EQUATORIAL_EARTH_RADIUS**2
    per2 = POLAR_EARTH_RADIUS**2

    rho = eer2/np.sqrt(eer2 * np.cos(phi) + per2 * np.sin(phi))

    x = (rho + ele) * np.cos(phi) * np.cos(theta)
    y = (rho + ele) * np.cos(phi) * np.sin(theta)
    z = ((per2/eer2) * rho + ele) * np.sin(phi)

    return round(x, 3), round(y, 3), round(z, 3)

def geocentric_radius(lat):
    """
    Compute the distance from the Earth's center
    to a point on the spheroid surface at given
    geodetic latitude.
    """

    phi = np.radians(lat)

    a = EQUATORIAL_EARTH_RADIUS**2
    b = POLAR_EARTH_RADIUS**2
    c = (b/a) * np.tan(phi)**2

    radius = np.sqrt((a + b * c)/(1 + c))

    return round(radius, 3)

def tunnel_distance(lat1, lon1, ele1, lat2, lon2, ele2):
    """
    Compute the linear distance (circle chord) between
    two points on the earth ellipsoid surface.
    """

    x1, y1, z1 = wgs_to_xyz_ellipsoid(lat1, lon1, ele1)
    x2, y2, z2 = wgs_to_xyz_ellipsoid(lat2, lon2, ele2)

    distance = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    return round(distance, 3)

def circle_distance(lat1, lon1, lat2, lon2):
    """
    Compute the great circle distance between two
    points on the spherical earth surface.

    Modified from: Salvador Dali (?)
    http://stackoverflow.com/users/1090562/salvador-dali
    """

    p = np.pi/180

    c1 = np.cos((lat2 - lat1) * p)
    c2 = np.cos(lat1 * p)
    c3 = np.cos(lat2 * p)
    c4 = np.cos((lon2 - lon1) * p)

    a = 0.5 - c1/2 + c2 * c3 * (1 - c4) / 2

    distance = 2 * MEAN_EARTH_RADIUS * np.arcsin(np.sqrt(a))

    return round(distance, 3)

def circle_distance_to_test(lat1, lon1, lat2, lon2):
    """
    Compute the great circle distance between two
    points on the spherical earth surface.
    """
    
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    distance = MEAN_EARTH_RADIUS * c
    
    return round(distance, 3)











