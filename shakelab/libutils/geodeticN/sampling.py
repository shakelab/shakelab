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
Random sampling of points inside polygons.

The module draws points uniformly distributed inside polygons, either in
planar coordinates (:func:`sample_uniform_in_polygon_xy`) or in
geographic coordinates (:func:`sample_uniform_in_polygon_wgs84`).

Design choices
--------------
- Planar sampling triangulates the polygon vertices (Delaunay), draws
  candidate points uniformly in the triangles and keeps those inside the
  polygon (rejection sampling), so concave polygons are supported.
  Candidates are drawn and tested in vectorized batches.
- Geographic sampling works in the equal-area sinusoidal projection
  centered on the mean longitude of each part (as
  :meth:`WgsPolygon.contains_points`), so that uniformity in the plane
  is uniformity in area; parts receive samples in proportion to their
  area.
- The random source is the ``rng`` argument: a
  :class:`numpy.random.Generator`, an integer seed, or ``None`` to use
  the global NumPy random state (so that ``numpy.random.seed`` keeps
  results reproducible).

Conventions
-----------
- Geographic coordinates are ``(lon, lat)`` in degrees; sampled
  longitudes are in the convention of the first vertex of their part.
- Planar coordinates and distances are in meters.
"""

from __future__ import annotations

from numbers import Integral
from typing import Any, List, Optional, Union

import numpy as np
from scipy.spatial import Delaunay

from shakelab.libutils.geodeticN.primitives import (
    WgsMesh,
    WgsPolygon,
    _open_ring,
    _points_in_ring,
    _reference_longitude,
    _signed_area,
)
from shakelab.libutils.geodeticN.utilities import (
    sinusoidal_to_wgs84,
    wgs84_to_sinusoidal,
)


__all__ = [
    "sample_uniform_in_polygon_xy",
    "sample_uniform_in_polygon_wgs84",
]


RandomSource = Optional[Union[np.random.Generator, int]]


# ---------------------------------------------------------------------------
# Planar polygons
# ---------------------------------------------------------------------------

def sample_uniform_in_polygon_xy(
    vertices_xy: np.ndarray,
    n_samples: int,
    rng: RandomSource = None,
) -> np.ndarray:
    """
    Draw points uniformly inside a planar polygon.

    Parameters
    ----------
    vertices_xy
        Polygon vertices, shape ``(n_vertices, 2)``. The polygon must be
        simple (non self-intersecting); a closing vertex is allowed.
    n_samples
        Number of points to draw.
    rng
        Random source (see the module design choices).

    Returns
    -------
    numpy.ndarray
        Sampled points, shape ``(n_samples, 2)``.

    Raises
    ------
    ValueError
        If ``n_samples`` is not positive, the vertices have a wrong
        shape, there are fewer than three distinct vertices, or the
        triangulation is degenerate.
    """
    _require_positive_count(n_samples)
    generator = _random_source(rng)

    vertices = np.asarray(vertices_xy, dtype=float)

    if vertices.ndim != 2 or vertices.shape[1] != 2:
        raise ValueError("vertices_xy must have shape (n, 2).")

    if len(vertices) > 1 and np.array_equal(vertices[0], vertices[-1]):
        vertices = vertices[:-1]

    if vertices.shape[0] < 3:
        raise ValueError("Polygon must have at least three vertices.")

    triangles = vertices[Delaunay(vertices).simplices]
    areas = _triangle_areas(triangles)

    if not np.all(areas > 0.0):
        raise ValueError("Degenerate triangles found in triangulation.")

    cumulative = np.cumsum(areas)
    batch = max(64, 2 * n_samples)
    accepted: List[np.ndarray] = []
    count = 0

    while count < n_samples:
        index = np.searchsorted(
            cumulative,
            generator.random(batch) * cumulative[-1],
        )
        index = np.minimum(index, len(triangles) - 1)
        candidates = _sample_in_triangles(triangles[index], generator)

        inside = _points_in_ring(
            vertices[:, 0],
            vertices[:, 1],
            candidates[:, 0],
            candidates[:, 1],
        )
        accepted.append(candidates[inside])
        count += int(inside.sum())

    return np.vstack(accepted)[:n_samples]


# ---------------------------------------------------------------------------
# Geographic polygons
# ---------------------------------------------------------------------------

def sample_uniform_in_polygon_wgs84(
    polygon: WgsPolygon,
    n_samples: int,
    rng: RandomSource = None,
) -> WgsMesh:
    """
    Draw points uniformly (in area) inside a geographic polygon.

    Parameters
    ----------
    polygon
        Polygon or multipolygon.
    n_samples
        Number of points to draw.
    rng
        Random source (see the module design choices).

    Returns
    -------
    WgsMesh
        Sampled points. An empty polygon gives an empty mesh.

    Raises
    ------
    ValueError
        If ``n_samples`` is not positive.
    """
    _require_positive_count(n_samples)
    generator = _random_source(rng)

    parts = []
    areas = []

    for part in polygon.parts:
        ring = _open_ring(part)

        if len(ring) < 3:
            continue

        lon0 = _reference_longitude(ring)
        x, y = wgs84_to_sinusoidal(
            [p.longitude for p in ring],
            [p.latitude for p in ring],
            lon0_deg=lon0,
        )
        parts.append((np.column_stack((x, y)), lon0))
        areas.append(abs(_signed_area(x, y)))

    total = float(np.sum(areas)) if areas else 0.0

    if total <= 0.0:
        return WgsMesh()

    counts = generator.multinomial(n_samples, np.array(areas) / total)

    lons: List[np.ndarray] = []
    lats: List[np.ndarray] = []

    for (vertices, lon0), count in zip(parts, counts):
        if count <= 0:
            continue

        points = sample_uniform_in_polygon_xy(
            vertices, int(count), rng=generator,
        )
        lon, lat = sinusoidal_to_wgs84(
            points[:, 0], points[:, 1], lon0_deg=lon0,
        )
        lons.append(lon)
        lats.append(lat)

    return WgsMesh.from_arrays(np.concatenate(lons), np.concatenate(lats))


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _random_source(rng: RandomSource) -> Any:
    """
    Return an object with ``random`` and ``multinomial`` methods.

    ``None`` returns the global :mod:`numpy.random` state, an integer is
    used as seed of a new :class:`numpy.random.Generator`.
    """
    if rng is None:
        return np.random

    if isinstance(rng, np.random.Generator):
        return rng

    if isinstance(rng, Integral) and not isinstance(rng, bool):
        return np.random.default_rng(int(rng))

    if rng is np.random:
        return rng

    raise TypeError("rng must be a numpy Generator, an integer or None.")


def _require_positive_count(n_samples: int) -> None:
    """
    Raise ValueError if ``n_samples`` is not a positive integer.
    """
    if isinstance(n_samples, bool) or not isinstance(n_samples, Integral):
        raise TypeError("n_samples must be an integer.")

    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")


def _sample_in_triangles(
    triangles: np.ndarray,
    generator: Any,
) -> np.ndarray:
    """
    Draw one uniform point in each triangle of shape (n, 3, 2).
    """
    n = len(triangles)
    r1 = np.sqrt(generator.random(n))[:, None]
    r2 = generator.random(n)[:, None]

    return (
        (1.0 - r1) * triangles[:, 0]
        + r1 * (1.0 - r2) * triangles[:, 1]
        + r1 * r2 * triangles[:, 2]
    )


def _triangle_areas(triangles: np.ndarray) -> np.ndarray:
    """
    Return the areas of planar triangles of shape (n, 3, 2).
    """
    edge1 = triangles[:, 1] - triangles[:, 0]
    edge2 = triangles[:, 2] - triangles[:, 0]

    # Explicit 2D cross product: np.cross no longer accepts 2D vectors
    # in NumPy 2.
    cross = edge1[:, 0] * edge2[:, 1] - edge1[:, 1] * edge2[:, 0]

    return 0.5 * np.abs(cross)
