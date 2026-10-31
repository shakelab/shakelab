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
Unit tests for :mod:`shakelab.libutils.geodeticN.utilities`.

The tests cover the WGS84 constants, angle and longitude helpers, DMS
parsing, distances and azimuth on the sphere, the geocentric radius
and the spherical and sinusoidal coordinate transforms.

Conventions
-----------
- ``R`` is the mean Earth radius used as default by the distances.
- Longitudes are tested in both the ``[-180, 180)`` and the
  ``[0, 360)`` conventions.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from shakelab.libutils.geodeticN.utilities import (
    EQUATORIAL_EARTH_RADIUS_M,
    MEAN_EARTH_RADIUS_M,
    POLAR_EARTH_RADIUS_M,
    WGS84_A,
    WGS84_B,
    WGS84_E2,
    WGS84_EP2,
    cartesian_sphere_to_wgs84,
    chord_distance,
    deg_to_rad,
    dms2dec,
    forward_azimuth,
    geocentric_radius,
    great_circle_distance,
    normalize_longitude,
    rad_to_deg,
    relative_longitude,
    sinusoidal_to_wgs84,
    tunnel_distance_sphere,
    wgs84_to_cartesian_sphere,
    wgs84_to_sinusoidal,
)


R = MEAN_EARTH_RADIUS_M


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_wgs84_constants_are_consistent():
    """Polar radius and eccentricities follow from axis and flattening."""
    assert POLAR_EARTH_RADIUS_M == pytest.approx(6356752.314245, abs=1e-6)
    assert POLAR_EARTH_RADIUS_M == WGS84_B
    assert EQUATORIAL_EARTH_RADIUS_M == WGS84_A
    difference = WGS84_A ** 2 - WGS84_B ** 2
    assert WGS84_E2 == pytest.approx(difference / WGS84_A ** 2, rel=1e-12)
    assert WGS84_EP2 == pytest.approx(difference / WGS84_B ** 2, rel=1e-12)


# ---------------------------------------------------------------------------
# Angles and longitudes
# ---------------------------------------------------------------------------

def test_angle_conversions():
    """Degrees and radians convert both ways, for scalars and arrays."""
    assert deg_to_rad(180.0) == pytest.approx(math.pi)
    assert rad_to_deg(math.pi) == pytest.approx(180.0)
    assert isinstance(deg_to_rad(90.0), float)
    np.testing.assert_allclose(deg_to_rad([0.0, 90.0]), [0.0, math.pi / 2])


@pytest.mark.parametrize("lon, expected", [
    (190.0, -170.0),
    (-190.0, 170.0),
    (360.0, 0.0),
    (-540.0, -180.0),
    (180.0, -180.0),
    (359.5, -0.5),
])
def test_normalize_longitude(lon, expected):
    """Longitudes are wrapped to [-180, 180)."""
    assert normalize_longitude(lon) == pytest.approx(expected)


def test_normalize_longitude_arrays():
    """Arrays are normalized element-wise."""
    result = normalize_longitude([0.0, 360.0, 190.0])
    np.testing.assert_allclose(result, [0.0, 0.0, -170.0])


@pytest.mark.parametrize("lon, ref, expected", [
    (350.0, -10.0, 0.0),
    (-170.0, 170.0, 20.0),
    (170.0, -170.0, -20.0),
    (10.0, 350.0, 20.0),
])
def test_relative_longitude(lon, ref, expected):
    """Differences are wrapped and independent of the convention."""
    assert relative_longitude(lon, ref) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# DMS parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text, expected", [
    ("12", 12.0),
    (" 45.25 ", 45.25),
    ("12 30", 12.5),
    ("12 30 0", 12.5),
    ("12,5", 12.5),
    ("+12 30", 12.5),
    ("12°30'00\"N", 12.5),
    ("12°30'S", -12.5),
    ("-12°30'N", 12.5),
    ("-12:30:00", -12.5),
    ("-0:30:00", -0.5),
    ("12°34'56.7\"", 12.0 + 34.0 / 60.0 + 56.7 / 3600.0),
    ("12 34 56.7W", -(12.0 + 34.0 / 60.0 + 56.7 / 3600.0)),
])
def test_dms2dec_strings(text, expected):
    """DMS strings are parsed; a leading sign applies to the whole angle."""
    assert dms2dec(text) == pytest.approx(expected, rel=1e-12)


def test_dms2dec_numbers_and_sequences():
    """Numbers pass through; sequences return arrays."""
    assert dms2dec(12) == 12.0
    assert isinstance(dms2dec(12), float)
    result = dms2dec(["12 30", -1.5, "1 0 36S"])
    np.testing.assert_allclose(result, [12.5, -1.5, -1.01])


@pytest.mark.parametrize("text", ["", "N", "abc", "12 30 15 10", "1.2.3"])
def test_dms2dec_rejects_invalid_strings(text):
    """Unparsable strings or too many fields raise ValueError."""
    with pytest.raises(ValueError):
        dms2dec(text)


@pytest.mark.parametrize("value", [None, True, object()])
def test_dms2dec_rejects_invalid_types(value):
    """Values that are neither numbers nor strings raise TypeError."""
    with pytest.raises(TypeError):
        dms2dec(value)


# ---------------------------------------------------------------------------
# Distances and azimuth
# ---------------------------------------------------------------------------

def test_great_circle_quarter_equator():
    """A quarter of the equator has the expected length."""
    distance = great_circle_distance(0.0, 0.0, 90.0, 0.0)
    assert distance == pytest.approx(0.5 * math.pi * R, rel=1e-12)


def test_great_circle_zero_and_symmetry():
    """Distance is zero for equal points and symmetric otherwise."""
    assert great_circle_distance(12.0, 45.0, 12.0, 45.0) == 0.0
    forward = great_circle_distance(12.0, 45.0, 13.0, 46.0)
    backward = great_circle_distance(13.0, 46.0, 12.0, 45.0)
    assert forward == pytest.approx(backward, rel=1e-14)


@pytest.mark.parametrize("lon1, lat1, lon2, lat2", [
    (0.0, 0.0, 180.0, 0.0),
    (10.0, 20.0, -170.0, -20.0),
    (45.3, -12.7, -134.7, 12.7),
    (0.0, 90.0, 0.0, -90.0),
])
def test_great_circle_antipodal_points(lon1, lat1, lon2, lat2):
    """Antipodal points give half a circumference, never NaN."""
    distance = great_circle_distance(lon1, lat1, lon2, lat2)
    assert math.isfinite(distance)
    assert distance == pytest.approx(math.pi * R, rel=1e-9)


def test_great_circle_longitude_convention():
    """Longitudes in [0, 360) give the same distance."""
    assert great_circle_distance(350.0, 10.0, 20.0, 15.0) == pytest.approx(
        great_circle_distance(-10.0, 10.0, 20.0, 15.0), rel=1e-12
    )


def test_great_circle_vectorized():
    """Array inputs return an array of distances."""
    result = great_circle_distance(0.0, 0.0, [0.0, 90.0, 180.0], 0.0)
    assert isinstance(result, np.ndarray)
    np.testing.assert_allclose(result, [0.0, 0.5 * math.pi * R, math.pi * R])


def test_chord_is_shorter_than_arc_and_matches_tunnel():
    """The chord is shorter than the arc and equals the surface tunnel."""
    arc = great_circle_distance(0.0, 0.0, 60.0, 30.0)
    chord = chord_distance(0.0, 0.0, 60.0, 30.0)
    assert chord < arc
    assert chord == pytest.approx(
        tunnel_distance_sphere(0.0, 0.0, 60.0, 30.0), rel=1e-10
    )


def test_chord_antipodal():
    """The chord between antipodal points is the diameter."""
    assert chord_distance(0.0, 0.0, 180.0, 0.0) == pytest.approx(2.0 * R)


def test_chord_short_distance_precision():
    """The chord keeps full precision at millimetre distances."""
    chord = chord_distance(12.0, 45.0, 12.0, 45.0 + 1e-8)
    assert chord == pytest.approx(R * math.radians(1e-8), rel=1e-6)


def test_tunnel_distance_sphere_vertical_offset():
    """Points on the same vertical are separated by the height difference."""
    distance = tunnel_distance_sphere(
        12.0, 45.0, 12.0, 45.0, elev1_m=-10000.0, elev2_m=0.0,
    )
    assert distance == pytest.approx(10000.0, abs=1e-6)


@pytest.mark.parametrize("lon2, lat2, expected", [
    (90.0, 0.0, 90.0),
    (0.0, 10.0, 0.0),
    (0.0, -10.0, 180.0),
    (-90.0, 0.0, 270.0),
    (350.0, 0.0, 270.0),
])
def test_forward_azimuth(lon2, lat2, expected):
    """Azimuths from the origin follow the compass directions."""
    assert forward_azimuth(0.0, 0.0, lon2, lat2) == pytest.approx(expected)


def test_forward_azimuth_vectorized():
    """Array inputs return an array of azimuths."""
    result = forward_azimuth(0.0, 0.0, [90.0, 0.0], [0.0, 10.0])
    np.testing.assert_allclose(result, [90.0, 0.0], atol=1e-12)


# ---------------------------------------------------------------------------
# Ellipsoid helpers
# ---------------------------------------------------------------------------

def test_geocentric_radius():
    """The radius goes from the equatorial to the polar value."""
    assert geocentric_radius(0.0) == pytest.approx(WGS84_A, rel=1e-12)
    assert geocentric_radius(90.0) == pytest.approx(WGS84_B, rel=1e-12)
    assert WGS84_A > geocentric_radius(45.0) > WGS84_B
    assert geocentric_radius([0.0, 90.0]).shape == (2,)


# ---------------------------------------------------------------------------
# Spherical and sinusoidal transforms
# ---------------------------------------------------------------------------

def test_cartesian_sphere_round_trip():
    """Spherical Cartesian coordinates convert back to the input."""
    lon = np.array([0.0, 12.0, -120.0])
    lat = np.array([0.0, 45.0, -60.0])
    elevation = np.array([0.0, 1000.0, -5000.0])
    x, y, z = wgs84_to_cartesian_sphere(lon, lat, elevation)
    lon2, lat2, elevation2 = cartesian_sphere_to_wgs84(x, y, z)
    np.testing.assert_allclose(lon2, lon, atol=1e-10)
    np.testing.assert_allclose(lat2, lat, atol=1e-10)
    np.testing.assert_allclose(elevation2, elevation, atol=1e-6)


def test_sinusoidal_round_trip():
    """The sinusoidal projection is inverted exactly."""
    lon = np.array([-170.0, 0.0, 12.5, 179.0])
    lat = np.array([-60.0, 0.0, 45.0, 10.0])
    x, y = wgs84_to_sinusoidal(lon, lat)
    lon2, lat2 = sinusoidal_to_wgs84(x, y)
    np.testing.assert_allclose(lon2, lon, atol=1e-9)
    np.testing.assert_allclose(lat2, lat, atol=1e-9)


def test_sinusoidal_wraps_longitudes():
    """Longitudes in [0, 360) project like their [-180, 180) equivalent."""
    np.testing.assert_allclose(
        wgs84_to_sinusoidal(350.0, 10.0),
        wgs84_to_sinusoidal(-10.0, 10.0),
    )


def test_sinusoidal_central_meridian():
    """A central meridian keeps the antimeridian region continuous."""
    x, y = wgs84_to_sinusoidal([170.0, -170.0], [0.0, 0.0], lon0_deg=180.0)
    np.testing.assert_allclose(x, R * np.radians([-10.0, 10.0]))
    lon, _ = sinusoidal_to_wgs84(x, y, lon0_deg=180.0)
    np.testing.assert_allclose(lon, [170.0, 190.0])


def test_sinusoidal_pole_is_finite():
    """At the pole the inverse returns the central meridian."""
    lon, lat = sinusoidal_to_wgs84(0.0, R * math.pi / 2.0, lon0_deg=15.0)
    assert float(lon) == 15.0
    assert float(lat) == pytest.approx(90.0)
