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
Unit tests for :mod:`shakelab.libutils.geodeticN.mesh`.

The tests check the node counts and bounds of the Fibonacci, Cartesian
and sinusoidal generators, the alignment of the sinusoidal grid, the
handling of longitude ranges crossing the antimeridian or given in the
``[0, 360)`` convention, and the composition with polygon clipping.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from shakelab.libutils.geodeticN.mesh import (
    mesh_cartesian,
    mesh_fibonacci_bbox,
    mesh_fibonacci_global,
    mesh_fibonacci_spherical_cap,
    mesh_sinusoidal,
)
from shakelab.libutils.geodeticN.primitives import WgsMesh, WgsPolygon
from shakelab.libutils.geodeticN.utilities import (
    MEAN_EARTH_RADIUS_M,
    great_circle_distance,
    normalize_longitude,
)


def _node_set(mesh, wrap=False):
    """Return the nodes of a mesh as a set of rounded (lon, lat)."""
    lon, lat, _ = mesh.to_arrays()

    if wrap:
        lon = normalize_longitude(lon)

    return set(zip(np.round(lon, 6), np.round(lat, 6)))


def _assert_in_box(mesh, lon_min, lon_max, lat_min, lat_max, tol=1e-6):
    """Check that all nodes are inside a box."""
    lon, lat, _ = mesh.to_arrays()
    assert np.all(lon >= lon_min - tol)
    assert np.all(lon <= lon_max + tol)
    assert np.all(lat >= lat_min - tol)
    assert np.all(lat <= lat_max + tol)


# ---------------------------------------------------------------------------
# Fibonacci meshes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("spacing", [0.0, -10.0])
def test_fibonacci_global_rejects_invalid_spacing(spacing):
    """The spacing must be positive."""
    with pytest.raises(ValueError, match="spacing_km"):
        mesh_fibonacci_global(spacing)


def test_fibonacci_global_count_and_bounds():
    """The node count follows the target area per node."""
    mesh = mesh_fibonacci_global(spacing_km=200.0)
    expected = round(4.0 * math.pi * MEAN_EARTH_RADIUS_M ** 2 / 200e3 ** 2)
    assert isinstance(mesh, WgsMesh)
    assert len(mesh) == expected
    _assert_in_box(mesh, -180.0, 180.0, -90.0, 90.0, tol=0.0)


def test_fibonacci_cap_nodes_within_radius():
    """All cap nodes lie within the cap radius."""
    mesh = mesh_fibonacci_spherical_cap(13.0, 46.0, 100.0, 20.0)
    assert len(mesh) > 0
    lon, lat, _ = mesh.to_arrays()
    distance = great_circle_distance(13.0, 46.0, lon, lat)
    assert np.all(distance <= 100e3 * 1.001)


@pytest.mark.parametrize("max_distance, spacing", [(0.0, 20.0), (100.0, 0.0)])
def test_fibonacci_cap_rejects_invalid_parameters(max_distance, spacing):
    """Radius and spacing must be positive."""
    with pytest.raises(ValueError):
        mesh_fibonacci_spherical_cap(13.0, 46.0, max_distance, spacing)


def test_fibonacci_cap_keeps_center_convention():
    """Longitudes stay close to the center, in its convention."""
    mesh = mesh_fibonacci_spherical_cap(350.0, 0.0, 200.0, 50.0)
    lon, _, _ = mesh.to_arrays()
    assert np.all(np.abs(lon - 350.0) < 3.0)


def test_fibonacci_bbox_nodes_inside_box():
    """Nodes of the bbox variant are inside the box."""
    mesh = mesh_fibonacci_bbox(12.0, 14.0, 45.0, 47.0, spacing_km=20.0)
    assert len(mesh) > 0
    _assert_in_box(mesh, 12.0, 14.0, 45.0, 47.0)


def test_fibonacci_bbox_across_antimeridian():
    """A box from 170° to -170° keeps nodes in [170, 190]."""
    mesh = mesh_fibonacci_bbox(170.0, -170.0, -5.0, 5.0, spacing_km=50.0)
    assert len(mesh) > 0
    _assert_in_box(mesh, 170.0, 190.0, -5.0, 5.0)


@pytest.mark.parametrize("kwargs", [
    {"spacing_km": 0.0},
    {"spacing_km": 10.0, "lat_min": 2.0, "lat_max": 1.0},
])
def test_fibonacci_bbox_rejects_invalid_input(kwargs):
    """Spacing and latitude bounds are validated."""
    params = {"lon_min": 0.0, "lon_max": 1.0, "lat_min": 0.0,
              "lat_max": 1.0}
    params.update(kwargs)
    with pytest.raises(ValueError):
        mesh_fibonacci_bbox(**params)


# ---------------------------------------------------------------------------
# Cartesian grid
# ---------------------------------------------------------------------------

def test_cartesian_grid_count_and_bounds():
    """A 2° x 2° box with 1° steps has 9 nodes."""
    mesh = mesh_cartesian(-1.0, 1.0, -1.0, 1.0, 1.0, 1.0)
    assert len(mesh) == 9
    assert mesh.bounding_box() == (-1.0, 1.0, -1.0, 1.0)


def test_cartesian_grid_has_no_drift():
    """Fractional steps reach the upper bound exactly once."""
    mesh = mesh_cartesian(0.0, 1.0, 0.0, 1.0, 0.1, 0.1)
    assert len(mesh) == 121
    lon_max = mesh.bounding_box()[1]
    assert lon_max == pytest.approx(1.0, abs=1e-12)


def test_cartesian_grid_across_antimeridian():
    """Ranges crossing the antimeridian are continuous."""
    mesh = mesh_cartesian(170.0, -170.0, 0.0, 0.0, 5.0, 1.0)
    lon, _, _ = mesh.to_arrays()
    np.testing.assert_allclose(lon, [170.0, 175.0, 180.0, 185.0, 190.0])


@pytest.mark.parametrize("lon_step, lat_step", [(0.0, 1.0), (1.0, -1.0)])
def test_cartesian_grid_rejects_invalid_steps(lon_step, lat_step):
    """Steps must be positive."""
    with pytest.raises(ValueError):
        mesh_cartesian(0.0, 1.0, 0.0, 1.0, lon_step, lat_step)


# ---------------------------------------------------------------------------
# Sinusoidal grid
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("spacing", [0.0, -10.0])
def test_sinusoidal_rejects_invalid_spacing(spacing):
    """The spacing must be positive."""
    with pytest.raises(ValueError):
        mesh_sinusoidal(spacing_km=spacing)


def test_sinusoidal_rejects_invalid_latitudes():
    """Latitude bounds must be ordered and within [-90, 90]."""
    with pytest.raises(ValueError):
        mesh_sinusoidal(100.0, lat_min=10.0, lat_max=0.0)
    with pytest.raises(ValueError):
        mesh_sinusoidal(100.0, lat_max=95.0)


def test_sinusoidal_global_extent():
    """The global grid covers the whole sphere."""
    mesh = mesh_sinusoidal(spacing_km=500.0)
    lon_min, lon_max, lat_min, lat_max = mesh.bounding_box()
    assert lon_min <= -150.0
    assert lon_max >= 150.0
    assert lat_min <= -80.0
    assert lat_max >= 80.0


def test_sinusoidal_global_density():
    """The node count matches the target area per node."""
    mesh = mesh_sinusoidal(spacing_km=300.0)
    expected = 4.0 * math.pi * MEAN_EARTH_RADIUS_M ** 2 / 300e3 ** 2
    assert len(mesh) == pytest.approx(expected, rel=0.02)


def test_sinusoidal_window_is_a_subset_of_the_global_grid():
    """A window contains exactly the global nodes inside it."""
    global_mesh = mesh_sinusoidal(spacing_km=300.0)
    window = mesh_sinusoidal(300.0, 10.0, 20.0, 40.0, 50.0)
    assert len(window) > 0
    _assert_in_box(window, 10.0, 20.0, 40.0, 50.0)

    expected = {
        node for node in _node_set(global_mesh)
        if 10.0 <= node[0] <= 20.0 and 40.0 <= node[1] <= 50.0
    }
    assert _node_set(window) == expected


def test_sinusoidal_window_across_antimeridian():
    """A window from 170° to -170° returns longitudes in [170, 190]."""
    window = mesh_sinusoidal(100.0, 170.0, -170.0, -10.0, 10.0)
    assert len(window) > 0
    _assert_in_box(window, 170.0, 190.0, -10.0, 10.0)

    global_mesh = mesh_sinusoidal(spacing_km=100.0, lat_min=-10.0,
                                  lat_max=10.0)
    lon, lat, _ = global_mesh.to_arrays()
    offset = (lon - 170.0) % 360.0
    expected = int(np.sum(offset <= 20.0 + 1e-9))
    assert len(window) == expected


def test_sinusoidal_window_longitude_convention():
    """A window in [0, 360) has the same nodes as in [-180, 180)."""
    east = mesh_sinusoidal(100.0, 350.0, 370.0, -5.0, 5.0)
    west = mesh_sinusoidal(100.0, -10.0, 10.0, -5.0, 5.0)
    assert len(east) == len(west) > 0
    assert _node_set(east, wrap=True) == _node_set(west, wrap=True)
    _assert_in_box(east, 350.0, 370.0, -5.0, 5.0)


def test_sinusoidal_mesh_clipped_to_polygon():
    """A mesh inside a polygon is obtained by clipping."""
    polygon = WgsPolygon.from_lonlat([
        (10.0, 40.0), (14.0, 40.0), (12.0, 44.0),
    ])
    lon_min, lon_max, lat_min, lat_max = polygon.bounding_box()
    mesh = mesh_sinusoidal(50.0, lon_min, lon_max, lat_min, lat_max)
    clipped = mesh.clip(polygon)
    assert 0 < len(clipped) < len(mesh)
    lon, lat, _ = clipped.to_arrays()
    assert np.all(polygon.contains_points(lon, lat))
