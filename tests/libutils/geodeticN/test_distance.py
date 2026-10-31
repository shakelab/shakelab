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
Unit tests for :mod:`shakelab.libutils.geodeticN.distance`.

The tests check the epicentral and hypocentral distances, the two
tunnel distances and the handling of generic point objects.

Conventions
-----------
- Hypocentres below the surface have negative elevation in meters.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from shakelab.libutils.geodeticN import utilities
from shakelab.libutils.geodeticN.distance import (
    epicentral_distance,
    hypocentral_distance,
    tunnel_distance_ellipsoid,
    tunnel_distance_sphere,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint


HYPOCENTRE = WgsPoint(longitude=12.0, latitude=45.0, elevation=-10000.0)
STATION = WgsPoint(longitude=13.0, latitude=45.0, elevation=0.0)


def test_epicentral_matches_great_circle():
    """The epicentral distance is the great-circle distance."""
    expected = utilities.great_circle_distance(12.0, 45.0, 13.0, 45.0)
    assert epicentral_distance(HYPOCENTRE, STATION) == pytest.approx(
        expected, rel=1e-12
    )


def test_epicentral_ignores_elevation():
    """Elevation does not change the epicentral distance."""
    surface = WgsPoint(longitude=12.0, latitude=45.0)
    assert epicentral_distance(HYPOCENTRE, STATION) == pytest.approx(
        epicentral_distance(surface, STATION), rel=1e-14
    )


def test_hypocentral_exceeds_epicentral():
    """A deep hypocentre is farther than its epicentre."""
    assert hypocentral_distance(HYPOCENTRE, STATION) > epicentral_distance(
        HYPOCENTRE, STATION
    )


@pytest.mark.parametrize("approx", ["ellipsoid", "sphere", " Sphere "])
def test_hypocentral_vertical_offset(approx):
    """A station above the hypocentre is at the depth distance."""
    above = WgsPoint(longitude=12.0, latitude=45.0, elevation=0.0)
    distance = hypocentral_distance(HYPOCENTRE, above, approx=approx)
    assert distance == pytest.approx(10000.0, abs=1e-6)


def test_sphere_and_ellipsoid_are_close():
    """The two Earth models agree within one percent."""
    sphere = tunnel_distance_sphere(HYPOCENTRE, STATION)
    ellipsoid = tunnel_distance_ellipsoid(HYPOCENTRE, STATION)
    assert sphere > 0.0
    assert ellipsoid > 0.0
    assert abs(sphere - ellipsoid) / ellipsoid < 0.01


def test_sphere_model_is_consistent_with_epicentral():
    """At the surface and short range, the spherical chord is the arc."""
    point1 = WgsPoint(longitude=12.0, latitude=45.0)
    point2 = WgsPoint(longitude=12.0, latitude=45.01)
    chord = hypocentral_distance(point1, point2, approx="sphere")
    arc = epicentral_distance(point1, point2)
    assert chord == pytest.approx(arc, rel=1e-7)


def test_point_and_coordinate_sphere_tunnels_agree():
    """The point-based tunnel distance wraps the coordinate-based one."""
    expected = utilities.tunnel_distance_sphere(
        12.0, 45.0, 13.0, 45.0, elev1_m=-10000.0, elev2_m=0.0,
    )
    assert tunnel_distance_sphere(HYPOCENTRE, STATION) == pytest.approx(
        expected, rel=1e-14
    )


def test_hypocentral_rejects_invalid_model():
    """Unknown Earth models raise ValueError."""
    with pytest.raises(ValueError, match="approx"):
        hypocentral_distance(HYPOCENTRE, STATION, approx="invalid-model")


def test_generic_point_objects_are_accepted():
    """Objects with longitude and latitude work; elevation defaults to 0."""
    station = SimpleNamespace(longitude=13.0, latitude=45.0)
    assert hypocentral_distance(HYPOCENTRE, station) == pytest.approx(
        hypocentral_distance(HYPOCENTRE, STATION), rel=1e-14
    )


def test_objects_without_coordinates_are_rejected():
    """Objects without longitude and latitude raise TypeError."""
    with pytest.raises(TypeError, match="longitude"):
        epicentral_distance(HYPOCENTRE, object())


def test_longitude_convention_does_not_matter():
    """Longitudes in [0, 360) give the same distances."""
    west = WgsPoint(longitude=-10.0, latitude=40.0, elevation=-5000.0)
    east = WgsPoint(longitude=350.0, latitude=40.0, elevation=-5000.0)
    site = WgsPoint(longitude=-9.0, latitude=40.5)
    assert hypocentral_distance(east, site) == pytest.approx(
        hypocentral_distance(west, site), rel=1e-12
    )
    assert epicentral_distance(east, site) == pytest.approx(
        epicentral_distance(west, site), rel=1e-12
    )
