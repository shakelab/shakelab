# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2020, ShakeLab Developers.
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
#
# ****************************************************************************
"""
This modules includes a set of tools for geographical data manipulation.
"""

__all__ = ['read_geometry', 'write_geometry',
           'WgsPoint', 'WgsPolygon']

import numpy as np
import json

import matplotlib as _mat

# GEODETIC CONSTANTS
MEAN_EARTH_RADIUS = 6371008.8
EQUATORIAL_EARTH_RADIUS = 6378137.0
POLAR_EARTH_RADIUS = 6356752.3
DEG_TO_M = 111195.
NDIGITS = 4


def read_geometry(geometry_file):
    """
    """
    collection = []

    with open(geometry_file, 'r') as f:
        data = json.load(f)

        for feature in data['features']:
            ftype = feature['geometry']['type']
            fcoor = feature['geometry']['coordinates']
            fatt = feature['properties']

            if ftype == 'Point':
                item = WgsPoint()
                item.latitude = fcoor[1]
                item.longitue = fcoor[0]
                item.attributes = fatt
                collection.append(item)

            if ftype == "MultiPoint":
                for fpart in fcoor:
                    item = WgsPoint()
                    item.latitude = fcoor[1]
                    item.longitue = fcoor[0]
                    item.attributes = fatt
                    collection.append(item)

            if ftype == "LineString":
                pass

            if ftype == "MultiLineString":
                pass

            if ftype == 'Polygon':
                item = WgsPolygon()
                item.from_list(fcoor[0])
                item.attributes = fatt
                collection.append(item)

            if ftype == 'MultiPolygon':
                for fpart in fcoor:
                    item = WgsPolygon()
                    item.from_list(fpart[0])
                    item.attributes = fatt
                    collection.append(item)

    return collection

def write_geometry(collection, geometry_file, decimals=NDIGITS):
    """
    """

    data = {
        "type": "FeatureCollection",
        "features": []
        }

    for item in collection:
        feature = {
            'geometry': {},
            'properties': {}
            }

        if isinstance(item, WgsPoint):
            lat = np.round(item.latitude, decimals)
            lon = np.round(item.longitude, decimals)
            feature['geometry'] = {
                'type': 'Point',
                'coordinates': [lon, lat],
                'properties': {}
                }

        if isinstance(item, WgsPolygon):
            lat, lon = item.to_array()
            lat = np.round(lat, decimals)
            lon = np.round(lon, decimals)
            feature['geometry'] = {
                'type': 'Polygon',
                'coordinates': [[[x,y] for x,y in zip(lon, lat)]],
                'properties': {}
                }

        data['features'].append(feature)

    with open(geometry_file, 'w') as f:
        json.dump(data, f)


class WgsPoint():
    """
    Class to represent a geographical location of a point in the WGS84
    coordinate system.

    Attributes:
        latitude (float):
            The latitude of the point in decimal degrees.

        longitude (float):
            The longitude of the point in decimal degrees.

        elevation (float):
            The elevation of the point in meters, referred to the sea level.

        attributes (dict):
            Optional container for general meta-information

    Methods:
        circle_distance:

        tunnel_distance:
    """

    def __init__(self, latitude, longitude, elevation=0.,
                       laterror=0., lonerror=0., eleerror=0.):
        """
        Initialise the class by passing the coordinates of the point.
        If no elevation is provided, 0 m.a.s.l. is assumed.

        Args:
            latitude (float):
                The latitude of the point in decimal degrees.

            longitude (float):
                The longitude of the point in decimal degrees.

            elevation (float = 0.):
                The elevation of the point in meters, referred
                to the sea level.
        """
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation

        self.laterror = laterror
        self.lonerror = lonerror
        self.eleerror = eleerror

        self.attributes = {}

    def __call__(self):
        """
        Returns the coordinates of the point.
        """
        return self.latitude, self.longitude, self.elevation

    def __sub__(self, point):
        """
        Returns the tunnel distance in meters.
        """
        if isinstance(point, WgsPoint):
            return self.tunnel_distance(point)

    def circle_distance(self, point):
        """
        Compute the great-circle distance between current location and a
        given geographical point assuming a spherical earth.
        
        Args:
            point (WgsPoint):
                The second geographical location for the calculation.

        Returns:
            The circle distance between points in meters.

        Notes:
            Point elevation is not considered, therefore the two points
            (or their projection) are assumed to be on the earth surface.

        To do:
            calculation for ellipsoid should be implemented.
        """
        return circle_distance(self.latitude, self.longitude,
                               point.latitude, point.longitude)

    def tunnel_distance(self, point, approx='sphere'):
        """
        Compute the tunnel distance between the current location and a
        given geographical location. Calculation can be done assuming
        spherical or ellipsoidal earth.

        Args:
            point (WgsPoint):
                The second geographical location for the calculation.

            approx (str = 'sphere', 'ellipsoid'):
                The earth's form approximation used as reference frame.

        Returns:
            The circle tunnel between points in meters.
        """
        if approx == 'sphere':
            handle = tunnel_distance_sphere
        if approx == 'ellipsoid':
            handle = tunnel_distance_ellipsoid

        return handle(self.latitude, self.longitude, self.elevation,
                      point.latitude, point.longitude, point.elevation)


class WgsPolygon():
    """
    A polygon in geographical coordinates.
    Vertexes are in a list of WgsPoints.
    """

    def __init__(self, points=None):
        self.points = []
        if points is not None:
            self.from_list(points)
        self.attributes = {}

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            point = self.points[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return point

    def __call__(self):

        return self.to_array()

    def add(self, point):

        self.points.append(point)

    def from_array(self, latitude, longitude):

        for lat, lon in zip(latitude, longitude):
            self.add(WgsPoint(lat, lon))

    def from_list(self, points):

        for p in points:
            self.add(WgsPoint(p[0], p[1]))

    def to_array(self):
        """
        Export latitude and longitude as numpy arrays.
        """

        lat = np.array([v.latitude for v in self.points])
        lon = np.array([v.longitude for v in self.points])

        return lat, lon

    def to_list(self):
        """
        Export a list of tuples with geographical coordinates.
        """

        return [(v.latitude, v.longitude) for v in self.points]

    def get_bounds(self):

        lat, lon = self.to_array()

        return (min(lat), max(lat)), (min(lon), max(lon))

    def get_area(self):
        """
        Note: area is in square kilometers
        """

        # Coordinate conversion using equal area projection
        lat, lon = self.to_array()
        x, y = wgs_to_xy_sinproj(lat, lon)

        return 1e-6 * polygon_area_shoelace(x, y)

    def contains(self, point):
        """
        point is a WgsPoint object
        return boolean
        """

        poly_lat, poly_lon = self.to_array()

        return contains(poly_lon, poly_lat,
                        point.longitude, point.latitude)

    def create_mesh(self, delta, meters=False, mesh_type='cartesian'):

        bnd = self.get_bounds()
        if mesh_type == 'spherical':
            grd_lat, grd_lon = spherical_mesh(delta, meters,
                                              latlim=bnd[0],
                                              lonlim=bnd[1])

        if mesh_type == 'cartesian':
            grd_lat, grd_lon = cartesian_mesh(delta, meters,
                                              latlim=bnd[0],
                                              lonlim=bnd[1])

        mesh = WgsMesh()
        poly_lat, poly_lon = self.to_array()
        for lat, lon in zip(grd_lat, grd_lon):
            if contains(poly_lon, poly_lat, lon, lat):
                mesh.add(WgsPoint(lat, lon))

        return mesh

class WgsMesh():
    """
    """

    def __init__(self, ):
        self.points = []
        self.attributes = {}

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            point = self.points[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return point

    def add(self, point):
        self.points.append(point)

    def to_array(self):
        """
        Export latitude and longitude as numpy arrays.
        """

        lat = np.array([v.latitude for v in self.points])
        lon = np.array([v.longitude for v in self.points])

        return lat, lon

    def intersect(self, polygon):

        points = []
        for p in self.points:
            if polygon.contains(p):
                points.append(p)
        self.points = points

    def create_mesh(self, delta, meters=False, polygon=None,
                    latlim=(-90, 90), lonlim=(-180, 180),
                    mesh_type='cartesian'):

        if polygon is not None:
            latlim, lonlim = polygon.get_bounds()

        if mesh_type == 'spherical':
            grd_lat, grd_lon = spherical_mesh(delta, meters=meters,
                                              latlim=latlim,
                                              lonlim=lonlim)

        if mesh_type == 'cartesian':
            grd_lat, grd_lon = cartesian_mesh(delta, meters=meters,
                                              latlim=latlim,
                                              lonlim=lonlim)

        for lat, lon in zip(grd_lat, grd_lon):
            self.add(WgsPoint(lat, lon))

        if polygon is not None:
            self.intersect(polygon)

    def to_csv(self, csv_file):

        with open(csv_file, 'w') as f:
            f.write('latitude, longitude\n')
            for p in self.points:
                f.write('{0},{1}\n'.format(p.latitude, p.longitude))

# ----------------------------------------------------------------------------
# Geometric functions (cartesian x, y)

def contains(polygon_x, polygon_y, x, y):
    """
    Evaluate if a point is inside a polygon.
    Modified from: ????
    """

    n = len(list(zip(polygon_x, polygon_y)))

    x0 = polygon_x[0]
    y0 = polygon_y[0]

    result = False

    for i in range(n + 1):
        x1 = polygon_x[i % n]
        y1 = polygon_y[i % n]

        if min(y0, y1) < y <= max(y0, y1):
            if x <= max(x0, x1):
                if y0 != y1:
                    xints = (y-y0) * (x1-x0) / (y1-y0) + x0
                if x0 == x1 or x <= xints:
                    result = not result

        x0, y0 = x1, y1

    return result

def polygon_area(x, y):
    """
    Calculates the area of an arbitrary polygon given its verticies.
    Modified from Joe Kington
    """

    area = 0.0
    for i in range(-1, len(x)-1):
        area += x[i] * (y[i+1] - y[i-1])

    return abs(area) / 2.0

def polygon_area_shoelace(x, y):
    """
    Using Shoelace formula to compute area.
    """

    a = np.dot(x, np.roll(y, 1))
    b = np.dot(y, np.roll(x, 1))

    return 0.5*np.abs(a - b)

# ----------------------------------------------------------------------------
# Geodetic functions

def wgs_to_xy_sinproj(lat, lon):
    """
    Approximate conversion using sinusoidal projection.
    """

    y = np.radians(lat) * MEAN_EARTH_RADIUS
    x = np.radians(lon) * MEAN_EARTH_RADIUS * np.cos(np.radians(lat))

    return np.round(x, NDIGITS), np.round(y, NDIGITS)

def xy_to_wgs_sinproj(x, y):
    """
    Approximate conversion using sinusoidal projection.
    """

    lat = np.degrees(y / MEAN_EARTH_RADIUS)
    lon = np.degrees(x / (MEAN_EARTH_RADIUS * np.cos(np.radians(lat))))

    return np.round(lat, NDIGITS), np.round(lon, NDIGITS)

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

    return np.round(x, NDIGITS), np.round(y, NDIGITS), np.round(z, NDIGITS)

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

    return np.round(x, NDIGITS), np.round(y, NDIGITS), np.round(z, NDIGITS)

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

    return np.round(radius, NDIGITS)

def tunnel_distance_sphere(lat1, lon1, ele1, lat2, lon2, ele2):
    """
    Compute the linear distance (circle chord) between
    two points on the spherical earth surface.
    """

    x1, y1, z1 = wgs_to_xyz_sphere(lat1, lon1, ele1)
    x2, y2, z2 = wgs_to_xyz_sphere(lat2, lon2, ele2)

    distance = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    return np.round(distance, NDIGITS)

def tunnel_distance_ellipsoid(lat1, lon1, ele1, lat2, lon2, ele2):
    """
    Compute the linear distance (circle chord) between
    two points on the earth ellipsoid surface.
    """

    x1, y1, z1 = wgs_to_xyz_ellipsoid(lat1, lon1, ele1)
    x2, y2, z2 = wgs_to_xyz_ellipsoid(lat2, lon2, ele2)

    distance = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    return np.round(distance, NDIGITS)

def circle_distance_to_test(lat1, lon1, lat2, lon2):
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

    return np.round(distance, NDIGITS)

def circle_distance(lat1, lon1, lat2, lon2):
    """
    Compute the great circle distance between two
    points on the spherical earth surface using the
    Haversine formula
    Modified from: Chris Veness
    """
    
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = (lon2 - lon1)/2
    dlat = (lat2 - lat1)/2
    
    a = np.sin(dlat)**2 + np.sin(dlon)**2 * np.cos(lat1) * np.cos(lat2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    distance = MEAN_EARTH_RADIUS * c
    
    return np.round(distance, NDIGITS)

def wgs_to_azimuth(lat1, lon1, lat2, lon2):
    """
    Compute the azimuth (to north) of a segment by
    its vertex wgs coordinates.
    """

    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlon = lon2 - lon1
    a = np.sin(dlon)*np.cos(lat2)
    b = np.cos(lat1)*np.sin(lat2)
    c = np.sin(lat1)*np.cos(lat2)*np.cos(dlon)

    azimuth = np.degrees(np.arctan2(a, b-c))

    if azimuth < 0.:
        azimuth += 360.

    return azimuth

# ----------------------------------------------------------------------------

def spherical_mesh(delta, meters=False,
                   latlim=(-90, 90), lonlim=(-180, 180)):
    """
    Produce a spherical mesh using golder spiral algorithm.
    Delta is the average distance between nearby points
    (by default in degrees).

    Modified from Chris Drost.
    """

    if meters:
        # Distance is in meters
        num_pts = np.rint(4 * np.pi * MEAN_EARTH_RADIUS**2 / delta**2)
    else:
        # Distance is in degrees
        num_pts = np.rint(4 * np.pi / np.radians(delta)**2)

    indices = np.arange(0, num_pts, dtype=float) + 0.5

    phi = np.arccos(1 - 2 * (indices / num_pts))
    theta = np.pi * (1 + 5**0.5) * indices

    # Conversion to wgs84
    lat = np.degrees(phi) - 90.
    lon = np.degrees(unwrap(theta))

    i = (lat >= latlim[0]) & (lat <= latlim[1])
    j = (lon >= lonlim[0]) & (lon <= lonlim[1])
    lat = lat[i & j]
    lon = lon[i & j]

    return np.round(lat, NDIGITS), np.round(lon, NDIGITS)

def unwrap(angle):
    """
    Unwrap phase angle.
    """

    return angle - (2 * np.pi) * ((angle + np.pi)//(2 * np.pi))

def cartesian_mesh(delta, meters=False,
                   latlim=(-90, 90), lonlim=(-180, 180)):
    """
    Create a spherical mesh from a cartesian grid using
    sinusoidal projection.
    Such approach is inconsistent at the new-day line and should
    be used only locally and within longitude -180 to 180 degrees.
    Delta is in degrees if not speficied otherwise.
    """

    if not meters:
        delta *= DEG_TO_M

    r = MEAN_EARTH_RADIUS

    xlen =  np.pi * r
    ylen =  np.pi * r / 2

    xmax = xlen - (xlen % delta)
    ymax = ylen - (ylen % delta)

    xrng = np.arange(-xmax, xmax, delta)
    yrng = np.arange(-ymax, ymax, delta)

    ylim = np.radians(latlim) * r
    ysel = yrng[(yrng >= ylim[0]) & (yrng <= ylim[1])]

    lat = np.array([])
    lon = np.array([])

    for ys in ysel:
        xlim = np.radians(lonlim) * r * np.cos(ys / r)
        xsel = xrng[(xrng >= xlim[0]) & (xrng <= xlim[1])]

        lat0, lon0 = xy_to_wgs_sinproj(xsel, ys*np.ones(len(xsel)))

        lat = np.append(lat, lat0)
        lon = np.append(lon, lon0)

    return np.round(lat, NDIGITS), np.round(lon, NDIGITS)
