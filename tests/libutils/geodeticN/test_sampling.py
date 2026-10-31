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
Unit tests for :mod:`shakelab.libutils.geodeticN.sampling`.

The tests check argument validation, that samples lie inside planar and
geographic polygons (also concave, multipart and crossing the
antimeridian), that they are statistically uniform, and the
reproducibility granted by the ``rng`` argument.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.libutils.geodeticN.primitives import (
    WgsMesh,
    WgsPolygon,
    _points_in_ring,
)
from shakelab.libutils.geodeticN.sampling import (
    sample_uniform_in_polygon_wgs84,
    sample_uniform_in_polygon_xy,
)


TRIANGLE = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
L_SHAPE = np.array([
    [0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0],
])
SQUARE = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]


# ---------------------------------------------------------------------------
# Planar sampling
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("vertices, n_samples, error", [
    (TRIANGLE, 0, ValueError),
    (TRIANGLE, -10, ValueError),
    (TRIANGLE, 2.5, TypeError),
    (np.array([1.0, 2.0]), 10, ValueError),
    (np.array([[0.0, 0.0], [1.0, 0.0]]), 10, ValueError),
])
def test_xy_rejects_invalid_arguments(vertices, n_samples, error):
    """Counts and vertex arrays are validated."""
    with pytest.raises(error):
        sample_uniform_in_polygon_xy(vertices, n_samples)


def test_xy_triangle_samples():
    """Samples are inside the triangle, with centroid near (1/3, 1/3)."""
    points = sample_uniform_in_polygon_xy(TRIANGLE, 2000, rng=0)
    assert points.shape == (2000, 2)
    assert np.all(_points_in_ring(
        TRIANGLE[:, 0], TRIANGLE[:, 1], points[:, 0], points[:, 1],
    ))
    np.testing.assert_allclose(points.mean(axis=0), [1 / 3, 1 / 3],
                               atol=0.03)


def test_xy_concave_polygon():
    """Samples avoid the notch of an L-shaped polygon."""
    points = sample_uniform_in_polygon_xy(L_SHAPE, 3000, rng=1)
    in_notch = (points[:, 0] > 1.0) & (points[:, 1] > 1.0)
    assert not np.any(in_notch)
    # The L-shape area is 3; each unit square holds about a third.
    in_corner = (points[:, 0] < 1.0) & (points[:, 1] < 1.0)
    assert np.mean(in_corner) == pytest.approx(1 / 3, abs=0.04)


def test_xy_closed_ring_is_accepted():
    """A closing vertex is ignored."""
    closed = np.vstack([TRIANGLE, TRIANGLE[:1]])
    points = sample_uniform_in_polygon_xy(closed, 100, rng=2)
    assert points.shape == (100, 2)


def test_xy_reproducibility():
    """Equal seeds give equal samples; generators are accepted."""
    first = sample_uniform_in_polygon_xy(TRIANGLE, 50, rng=7)
    second = sample_uniform_in_polygon_xy(TRIANGLE, 50, rng=7)
    np.testing.assert_array_equal(first, second)
    third = sample_uniform_in_polygon_xy(
        TRIANGLE, 50, rng=np.random.default_rng(7),
    )
    np.testing.assert_array_equal(first, third)


def test_xy_global_seed_reproducibility():
    """Without rng, numpy.random.seed makes results reproducible."""
    np.random.seed(3)
    first = sample_uniform_in_polygon_xy(TRIANGLE, 20)
    np.random.seed(3)
    second = sample_uniform_in_polygon_xy(TRIANGLE, 20)
    np.testing.assert_array_equal(first, second)


def test_invalid_random_source():
    """Unsupported random sources raise TypeError."""
    with pytest.raises(TypeError, match="rng"):
        sample_uniform_in_polygon_xy(TRIANGLE, 10, rng="seed")


# ---------------------------------------------------------------------------
# Geographic sampling
# ---------------------------------------------------------------------------

def test_wgs84_rejects_invalid_counts():
    """The number of samples must be positive."""
    polygon = WgsPolygon.from_lonlat(SQUARE)
    with pytest.raises(ValueError):
        sample_uniform_in_polygon_wgs84(polygon, 0)


def test_wgs84_empty_polygon():
    """An empty polygon gives an empty mesh."""
    mesh = sample_uniform_in_polygon_wgs84(WgsPolygon(), 10)
    assert isinstance(mesh, WgsMesh)
    assert len(mesh) == 0


def test_wgs84_square_samples():
    """Samples are inside the polygon and centered on it."""
    polygon = WgsPolygon.from_lonlat(SQUARE)
    mesh = sample_uniform_in_polygon_wgs84(polygon, 1000, rng=0)
    lon, lat, _ = mesh.to_arrays()
    assert len(mesh) == 1000
    assert np.all(polygon.contains_points(lon, lat))
    assert abs(lon.mean()) < 0.1
    assert abs(lat.mean()) < 0.1


def test_wgs84_multipolygon_samples():
    """Parts receive samples in proportion to their area."""
    polygon = WgsPolygon.from_parts([
        [(-2.0, -1.0), (-1.0, -1.0), (-1.0, 0.0), (-2.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)],
    ])
    mesh = sample_uniform_in_polygon_wgs84(polygon, 2000, rng=0)
    lon, lat, _ = mesh.to_arrays()
    assert len(mesh) == 2000
    assert np.all(polygon.contains_points(lon, lat))
    assert np.mean(lon < 0.0) == pytest.approx(0.5, abs=0.05)


def test_wgs84_across_antimeridian():
    """Samples in a polygon crossing the antimeridian stay inside it."""
    polygon = WgsPolygon.from_lonlat([
        (170.0, -5.0), (-170.0, -5.0), (-170.0, 5.0), (170.0, 5.0),
    ])
    mesh = sample_uniform_in_polygon_wgs84(polygon, 500, rng=0)
    lon, lat, _ = mesh.to_arrays()
    assert np.all(polygon.contains_points(lon, lat))
    assert np.all((lon >= 170.0) & (lon <= 190.0))
    assert lon.mean() == pytest.approx(180.0, abs=1.0)
