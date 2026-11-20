"""Analytical regression tests for rectangular rupture geometry."""

import numpy as np
import pytest

from shakelab.groundmotion.rupture import RectangularRupture
from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import geo_to_enu


def _rupture(dip=90.0):
    return RectangularRupture(
        WgsPoint(13.0, 46.0, -2000.0),
        strike=0.0, dip=dip, length=10.0, width=6.0,
    )


def test_vertical_vertices_and_area():
    rupture = _rupture()
    assert rupture.vertices_enu == pytest.approx(np.array([
        [0, -5, 0], [0, 5, 0], [0, 5, -6], [0, -5, -6],
    ]), abs=1e-12)
    assert rupture.ztor == pytest.approx(2.0)
    assert rupture.area == pytest.approx(60.0)


def test_eastward_dip_and_vertex_round_trip():
    rupture = _rupture(30.0)
    vertices = rupture.vertices_enu
    assert vertices[3] - vertices[0] == pytest.approx(
        [3 * np.sqrt(3), 0, -3],
    )
    for site, expected in zip(rupture.vertices, vertices):
        actual = np.array(geo_to_enu(site, rupture.frame)) / 1000
        assert actual == pytest.approx(expected, abs=1e-7)


def test_hypocentre_is_inside_and_ztor_is_preserved():
    hypocentre = WgsPoint(13.0, 46.0, -5000.0)
    rupture = RectangularRupture.from_hypocentre(
        hypocentre, strike=120.0, dip=30.0, length=20.0,
        width=12.0, ztor=2.0, along_fraction=0.25,
    )
    assert rupture.down_fraction == pytest.approx(0.5)
    assert rupture.ztor == pytest.approx(2)
    assert rupture.point_enu(0.25, 0.5) == pytest.approx([0, 0, 0])
    assert rupture.distances([hypocentre])["rrup"] == pytest.approx([0])


def test_horizontal_requires_explicit_down_fraction():
    kwargs = dict(strike=0, dip=0, length=10, width=6,
                  ztor=2, along_fraction=0.5)
    reference = WgsPoint(13, 46, -2000)
    with pytest.raises(ValueError, match="down_fraction"):
        RectangularRupture.from_hypocentre(reference, **kwargs)
    rupture = RectangularRupture.from_hypocentre(
        reference, down_fraction=0.25, **kwargs,
    )
    assert rupture.point_enu(0.5, 0.25) == pytest.approx([0, 0, 0])


@pytest.mark.parametrize("depth", [1, 9])
def test_hypocentre_outside_depth_range(depth):
    with pytest.raises(ValueError, match="outside"):
        RectangularRupture.from_hypocentre(
            WgsPoint(13, 46, -1000 * depth),
            strike=0, dip=90, length=10, width=6,
            ztor=2, along_fraction=0.5,
        )


@pytest.mark.parametrize("override", [
    {"strike": 361}, {"dip": -1}, {"width": 0},
    {"length": float("nan")}, {"along_fraction": -0.1},
    {"down_fraction": 1.1}, {"down_fraction": 1.0},
])
def test_invalid_geometry(override):
    kwargs = dict(strike=0, dip=90, length=10, width=6)
    kwargs.update(override)
    with pytest.raises(ValueError):
        RectangularRupture(WgsPoint(13, 46, -2000), **kwargs)


def test_vertical_distances_analytical():
    rupture = _rupture()
    sites = np.array([[0, 0, 2], [3, 0, 2], [-3, 0, 2],
                      [3, 9, 2], [0, 0, -3]])
    result = rupture.distances_enu(sites)
    assert result["rrup"] == pytest.approx(
        [2, np.sqrt(13), np.sqrt(13), np.sqrt(29), 0],
    )
    assert result["rjb"] == pytest.approx([0, 3, 3, 5, 0])
    assert result["rx"] == pytest.approx([0, 3, -3, 3, 0])
    assert result["ry0"] == pytest.approx([0, 0, 0, 4, 0])


def test_dipping_plane_footprint_and_perpendicular_distance():
    rupture = _rupture(30.0)
    result = rupture.distances_enu(np.array([[2, 0, 0]]))
    assert result["rjb"] == pytest.approx([0])
    assert result["rrup"] == pytest.approx([1])
    assert result["rx"] == pytest.approx([2])


def test_strike_rotation_preserves_distances():
    first = _rupture(30)
    second = RectangularRupture(
        first.reference, strike=90, dip=30, length=10, width=6,
    )
    sites = np.array([[3, 9, 2], [-3, 0, 2]])
    rotated = sites[:, [1, 0, 2]].copy()
    rotated[:, 1] *= -1
    for key, distances in first.distances_enu(sites).items():
        assert second.distances_enu(rotated)[key] == pytest.approx(distances)


def test_empty_and_invalid_site_arrays():
    rupture = _rupture()
    assert rupture.distances([])["rrup"].size == 0
    with pytest.raises(ValueError, match="shape"):
        rupture.distances_enu(np.zeros((3, 2)))
    with pytest.raises(ValueError, match="finite"):
        rupture.distances_enu(np.array([[0, 0, np.nan]]))
