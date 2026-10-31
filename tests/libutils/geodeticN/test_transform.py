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
Unit tests for :mod:`shakelab.libutils.geodeticN.transform`.

The tests check the ECEF and ENU transforms (absolute and local), the
high-level metric helpers and the frame validation.

Conventions
-----------
- Round trips are checked to 1e-8 degrees (about 1 mm) and 1 mm in
  height, unless stated otherwise.
"""

from __future__ import annotations

import pytest

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    WGS84_A,
    WGS84_B,
    MetricFrame,
    MetricPoint,
    ecef_to_geo,
    enu_to_geo,
    geo_to_ecef,
    geo_to_enu,
    metric_to_wgs,
    wgs_to_metric,
)


REFERENCE = WgsPoint(longitude=13.0, latitude=46.0, elevation=200.0)


def _assert_same_geo(actual, expected, angle_tol=1e-8, height_tol=1e-3):
    """Compare two (lon, lat, h) triples."""
    assert actual[0] == pytest.approx(expected[0], abs=angle_tol)
    assert actual[1] == pytest.approx(expected[1], abs=angle_tol)
    assert actual[2] == pytest.approx(expected[2], abs=height_tol)


def _geo(point):
    """Return a WgsPoint as a (lon, lat, h) tuple."""
    return point.longitude, point.latitude, point.elevation


# ---------------------------------------------------------------------------
# ECEF
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geo", [
    (0.0, 0.0, 0.0),
    (12.0, 45.0, 1000.0),
    (-123.5, 60.0, -500.0),
    (179.9, -33.0, 25.0),
])
def test_geo_ecef_round_trip(geo):
    """geo_to_ecef and ecef_to_geo invert each other."""
    _assert_same_geo(ecef_to_geo(geo_to_ecef(geo)), geo)


@pytest.mark.parametrize("geo", [
    (10.0, 89.9999, 100.0),
    (-45.0, -89.99999, -500.0),
    (0.0, 90.0, 0.0),
    (30.0, 89.0, 8000.0),
])
def test_geo_ecef_round_trip_near_poles(geo):
    """The height stays accurate close to the poles."""
    lon, lat, h = ecef_to_geo(geo_to_ecef(geo))
    assert lat == pytest.approx(geo[1], abs=1e-8)
    assert h == pytest.approx(geo[2], abs=1e-3)


def test_ecef_reference_values():
    """Equator and pole map to the ellipsoid axes."""
    x, y, z = geo_to_ecef((0.0, 0.0, 0.0))
    assert (x, y, z) == pytest.approx((WGS84_A, 0.0, 0.0), abs=1e-6)
    x, y, z = geo_to_ecef((0.0, 90.0, 0.0))
    assert (x, y, z) == pytest.approx((0.0, 0.0, WGS84_B), abs=1e-6)


def test_ecef_to_geo_on_the_polar_axis():
    """Points on the polar axis are handled explicitly."""
    lon, lat, h = ecef_to_geo((0.0, 0.0, WGS84_B + 10.0))
    assert (lon, lat, h) == pytest.approx((0.0, 90.0, 10.0), abs=1e-9)


def test_two_element_tuples_have_zero_height():
    """(lon, lat) tuples are accepted with h = 0."""
    assert geo_to_ecef((12.0, 45.0)) == pytest.approx(
        geo_to_ecef((12.0, 45.0, 0.0))
    )


def test_ecef_longitude_convention():
    """Longitudes in [0, 360) give the same ECEF coordinates."""
    assert geo_to_ecef((350.0, 10.0, 0.0)) == pytest.approx(
        geo_to_ecef((-10.0, 10.0, 0.0)), abs=1e-6
    )


# ---------------------------------------------------------------------------
# Metric helpers and frames
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("point", [
    WgsPoint(0.0, 0.0, 0.0),
    WgsPoint(13.0, 46.0, 200.0),
    WgsPoint(-75.0, -30.0, 1500.0),
])
def test_wgs_metric_round_trip_absolute(point):
    """Without a frame, metric coordinates are absolute ECEF."""
    metric = wgs_to_metric(point)
    assert metric.frame is None
    _assert_same_geo(_geo(metric_to_wgs(metric)), _geo(point))


def test_local_ecef_frame():
    """A local ECEF frame maps the reference point to the origin."""
    frame = MetricFrame(ref_geo=REFERENCE, orientation="ecef")
    origin = wgs_to_metric(REFERENCE, frame=frame)
    assert (origin.x, origin.y, origin.z) == pytest.approx(
        (0.0, 0.0, 0.0), abs=1e-6
    )
    point = WgsPoint(13.01, 46.02, 250.0)
    back = metric_to_wgs(wgs_to_metric(point, frame=frame))
    _assert_same_geo(_geo(back), _geo(point))


@pytest.mark.parametrize("point", [
    REFERENCE,
    WgsPoint(13.01, 46.02, 250.0),
    WgsPoint(12.98, 45.99, 150.0),
])
def test_geo_enu_round_trip(point):
    """geo_to_enu and enu_to_geo invert each other."""
    frame = MetricFrame(ref_geo=REFERENCE, orientation="enu")
    enu = geo_to_enu(point, frame=frame)
    _assert_same_geo(enu_to_geo(enu, frame=frame), _geo(point))


@pytest.mark.parametrize("point", [
    REFERENCE,
    WgsPoint(13.02, 45.98, 150.0),
    WgsPoint(12.99, 46.03, 300.0),
])
def test_wgs_metric_round_trip_enu(point):
    """High-level helpers round-trip in an ENU frame."""
    frame = MetricFrame(ref_geo=REFERENCE, orientation="enu")
    metric = wgs_to_metric(point, frame=frame)
    assert isinstance(metric, MetricPoint)
    assert metric.frame == frame
    _assert_same_geo(_geo(metric.to_wgs()), _geo(point))


def test_wgs_metric_matches_geo_to_enu():
    """wgs_to_metric with an ENU frame matches geo_to_enu."""
    frame = MetricFrame(ref_geo=REFERENCE, orientation="enu")
    point = WgsPoint(13.015, 46.005, 220.0)
    metric = wgs_to_metric(point, frame=frame)
    assert (metric.x, metric.y, metric.z) == pytest.approx(
        geo_to_enu(point, frame=frame), abs=1e-6
    )


def test_enu_axes_point_east_north_up():
    """A point to the east has positive East and negligible North."""
    frame = MetricFrame(ref_geo=(13.0, 46.0, 0.0), orientation="enu")
    east, north, up = geo_to_enu((13.001, 46.0, 0.0), frame=frame)
    assert east > 0.0
    assert abs(north) < 1e-3 * east
    _, north, _ = geo_to_enu((13.0, 46.001, 0.0), frame=frame)
    assert north > 0.0
    _, _, up = geo_to_enu((13.0, 46.0, 50.0), frame=frame)
    assert up == pytest.approx(50.0, abs=1e-6)


def test_frame_orientation_is_normalized():
    """The orientation is case and whitespace insensitive."""
    frame = MetricFrame(ref_geo=(13.0, 46.0), orientation=" ENU ")
    assert frame.orientation == "enu"
    assert frame.ref_geo == (13.0, 46.0, 0.0)


@pytest.mark.parametrize("orientation, error", [
    ("ned", ValueError),
    (1, TypeError),
])
def test_invalid_frame_orientation(orientation, error):
    """Unknown orientations are rejected."""
    with pytest.raises(error):
        MetricFrame(ref_geo=REFERENCE, orientation=orientation)


def test_frame_orientation_mismatch():
    """ECEF and ENU helpers reject frames of the other kind."""
    enu = MetricFrame(ref_geo=REFERENCE, orientation="enu")
    ecef = MetricFrame(ref_geo=REFERENCE, orientation="ecef")
    with pytest.raises(ValueError, match="ecef"):
        geo_to_ecef(REFERENCE, frame=enu)
    with pytest.raises(ValueError, match="enu"):
        geo_to_enu(REFERENCE, frame=ecef)
    with pytest.raises(TypeError):
        geo_to_enu(REFERENCE, frame=None)


def test_metric_to_wgs_rejects_incompatible_frames():
    """A point frame different from the explicit frame is an error."""
    enu = MetricFrame(ref_geo=REFERENCE, orientation="enu")
    other = MetricFrame(ref_geo=(0.0, 0.0, 0.0), orientation="enu")
    metric = wgs_to_metric(REFERENCE, frame=enu)
    with pytest.raises(ValueError, match="incompatible"):
        metric_to_wgs(metric, frame=other)
