# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
Unit tests for :mod:`shakelab.exploration.yarray.response`.

Analytic responses (two stations, regular grids), invariances,
coordinate uncertainty, automatic grid, limits, search domains,
velocity conversion, arrays read from file and errors.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.response import (
    ArrayGeometry,
    ArrayResponse,
    ResponseAnalysis,
    ResponseGrid,
)

DATA = Path(__file__).parent / "data"
SEVEN = np.array([[0, 0], [30, 0], [0, 40], [25, 35], [-20, 10],
                  [10, -30], [-35, -20]], dtype=float)


def _sesame():
    """SESAME array (its SAC headers hold local coordinates)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="STLO")
        return ArrayData(DATA / "sesame_array.json")


def _square(count=4, spacing=10.0):
    """Regular square grid of stations."""
    axis = np.arange(count) * spacing
    x, y = np.meshgrid(axis, axis)
    return np.column_stack((x.ravel(), y.ravel()))


# ----------------------------------------------------------------------------
# Geometry and grid
# ----------------------------------------------------------------------------


def test_geometry_summary():
    geometry = ArrayGeometry(_square(3, 10.0))
    assert geometry.n_stations == 9
    assert geometry.rank == 2 and geometry.is_two_dimensional
    assert geometry.minimum_spacing == pytest.approx(10.0)
    assert geometry.maximum_aperture == pytest.approx(np.hypot(20, 20))
    assert geometry.anisotropy_ratio == pytest.approx(1.0)
    assert np.allclose(geometry.centroid, [10.0, 10.0])
    assert geometry.rayleigh_proxy == pytest.approx(
        2 * np.pi / geometry.maximum_aperture
    )
    assert geometry.nyquist_proxy == pytest.approx(np.pi / 10.0)


def test_geometry_linear_and_errors():
    line = ArrayGeometry([[0, 0], [10, 0], [25, 0]])
    assert line.is_linear and line.anisotropy_ratio == np.inf

    with pytest.raises(ValueError):
        ArrayGeometry([[0, 0]])

    with pytest.raises(ValueError):
        ArrayGeometry([[0, 0], [0, 0], [1, 1]])

    with pytest.raises(ValueError):
        ArrayGeometry([[0, 0], [np.nan, 1]])

    with pytest.raises(ValueError):
        ArrayGeometry("not coordinates")


def test_geometry_from_array():
    array = _sesame()
    geometry = ArrayGeometry(array)
    assert np.allclose(geometry.coordinates, array.coordinates[:, :2])


def test_response_grid():
    grid = ResponseGrid(1.0, 11)
    assert grid.values[grid.center_index] == 0.0
    assert grid.step == pytest.approx(0.2)

    with pytest.raises(ValueError):
        ResponseGrid(1.0, 10)

    with pytest.raises(ValueError):
        ResponseGrid(-1.0, 11)

    with pytest.raises(TypeError):
        ResponseGrid(1.0, 11.0)


def test_automatic_grid():
    geometry = ArrayGeometry(SEVEN)
    automatic = ArrayResponse(SEVEN).grid
    assert automatic.size % 2 == 1
    assert automatic.maximum == pytest.approx(
        1.25 * 2 * np.pi / geometry.minimum_spacing
    )
    assert automatic.step <= geometry.rayleigh_proxy / 12

    extent = ArrayResponse(SEVEN, grid_max=0.3).grid
    assert extent.maximum == 0.3
    assert extent.step == pytest.approx(automatic.step, rel=0.05)
    assert ArrayResponse(SEVEN, grid_size=31).grid.maximum == \
        automatic.maximum

    with pytest.raises(ValueError, match="grid_size"):
        ArrayResponse(SEVEN, grid_max=50.0)

    with pytest.raises(ValueError):
        ArrayResponse(SEVEN, grid_max=-1.0)


# ----------------------------------------------------------------------------
# Response
# ----------------------------------------------------------------------------


def test_two_stations_analytic():
    spacing = 20.0
    response = ArrayResponse([[0, 0], [spacing, 0]],
                              grid_max=0.5, grid_size=101)
    kx, ky = response.grid.mesh()
    assert np.allclose(response.power, np.cos(kx * spacing / 2) ** 2)
    assert response.power[50, 50] == 1.0


def test_translation_invariance():
    first = ArrayResponse(SEVEN, grid_max=0.3, grid_size=61)
    second = ArrayResponse(SEVEN + [1000.0, -500.0],
                            grid_max=0.3, grid_size=61)
    assert np.allclose(first.power, second.power, atol=1e-12)


def test_rotation_by_ninety_degrees():
    first = ArrayResponse(SEVEN, grid_max=0.3, grid_size=61)
    rotated = SEVEN @ np.array([[0.0, -1.0], [1.0, 0.0]])
    second = ArrayResponse(rotated, grid_max=0.3, grid_size=61)
    assert np.allclose(np.rot90(first.power), second.power, atol=1e-12) or \
        np.allclose(np.rot90(first.power, -1), second.power, atol=1e-12)


def test_expected_response():
    grid = {"grid_max": 0.3, "grid_size": 61}
    zero = ArrayResponse(SEVEN, coordinate_covariances=np.zeros((7, 2, 2)),
                         **grid)
    assert np.allclose(zero.expected_power, zero.power, atol=1e-12)

    noisy = ArrayResponse(
        SEVEN, **grid,
        coordinate_covariances=np.repeat(np.eye(3)[None] * 25.0, 7, axis=0),
    )
    center = noisy.grid.center_index
    assert noisy.expected_power[center, center] == 1.0
    remote = np.hypot(*noisy.grid.mesh()) > 0.15
    assert np.all(noisy.expected_power[remote]
                  <= noisy.power[remote] + 1.0 / 7 + 1e-12)
    assert noisy.expected_power[remote].mean() < noisy.power[remote].mean()

    assert noisy.selected_power(expected=True) is noisy.expected_power

    with pytest.raises(ValueError, match="expected"):
        ArrayResponse(SEVEN, **grid).selected_power(expected=True)

    with pytest.raises(ValueError):
        ArrayResponse(SEVEN, coordinate_covariances=np.ones((7, 2, 3)))

    with pytest.raises(ValueError):
        bad = np.repeat(np.array([[[1.0, 0.0], [0.0, -1.0]]]), 7, axis=0)
        ArrayResponse(SEVEN, coordinate_covariances=bad)


def test_response_of_array_data():
    array = _sesame()
    response = ArrayResponse(array)
    assert response.geometry.n_stations == 9
    assert response.expected_power is None


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------


def test_regular_grid_aliases():
    spacing = 10.0
    response = ArrayResponse(_square(4, spacing),
                              grid_max=1.0, grid_size=201)
    analysis = response.analyze()
    assert isinstance(analysis, ResponseAnalysis)
    period = 2 * np.pi / spacing
    assert analysis.ambiguity.peak_sidelobe_level == pytest.approx(1.0,
                                                                   abs=5e-3)
    assert analysis.ambiguity.peak_distance == pytest.approx(period,
                                                             rel=0.01)
    assert analysis.ambiguity.onset_distance < period
    assert analysis.limits.aliasing_radius == pytest.approx(
        analysis.ambiguity.onset_distance / 2
    )
    # The main lobe is symmetric under the 90-degree symmetry.
    radii = analysis.main_lobe.radii
    assert radii[0] == pytest.approx(radii[90], rel=1e-6)
    assert analysis.limits.resolution_radius == \
        analysis.main_lobe.worst_radius


def test_main_lobe_and_level():
    response = ArrayResponse(SEVEN)
    narrow = response.analyze(resolution_level=0.9)
    wide = response.analyze(0.3)
    assert narrow.main_lobe.worst_radius < wide.main_lobe.worst_radius
    lobe = wide.main_lobe
    assert lobe.mask[response.grid.center_index, response.grid.center_index]
    assert lobe.best_radius <= lobe.equivalent_radius <= lobe.worst_radius
    assert lobe.worst_resolution_separation == 2 * lobe.worst_radius
    proxy = response.geometry.rayleigh_proxy
    assert 0.2 * proxy < lobe.worst_radius < 2.0 * proxy


def test_analysis_errors():
    line = ArrayResponse([[0, 0], [10, 0], [25, 0]],
                          grid_max=1.0, grid_size=51)

    with pytest.raises(ValueError, match="linear"):
        line.analyze()

    response = ArrayResponse(SEVEN)

    with pytest.raises(ValueError):
        response.analyze(resolution_level=1.0)

    with pytest.raises(ValueError, match="boundary"):
        ArrayResponse(SEVEN, grid_max=0.01, grid_size=11).analyze()

    with pytest.raises(TypeError):
        response.analyze(azimuth_size=4.5)


def test_search_domain():
    response = ArrayResponse(_square(4, 10.0),
                              grid_max=1.6, grid_size=201)
    analysis = response.analyze()
    safe = response.search_domain(0.9 * analysis.limits.aliasing_radius)
    assert safe.is_safe and safe.maximum_ambiguity == 0.0
    assert safe.difference_radius == pytest.approx(2 * safe.maximum_wavenumber)
    unsafe = response.search_domain(0.4)
    assert not unsafe.is_safe
    assert unsafe.maximum_ambiguity == pytest.approx(1.0, abs=5e-3)
    assert unsafe.ambiguity_peaks

    with pytest.raises(ValueError, match="grid_max"):
        response.search_domain(1.0)


def test_velocity_limits():
    analysis = ArrayResponse(SEVEN).analyze()
    limits = analysis.velocity_limits([2.0, 4.0])
    assert np.allclose(limits.minimum, 2 * np.pi * np.array([2.0, 4.0])
                       / analysis.limits.aliasing_radius)
    assert np.allclose(limits.maximum, 2 * np.pi * np.array([2.0, 4.0])
                       / analysis.limits.resolution_radius)
    assert limits.resolution_level == 0.5

    with pytest.raises(ValueError):
        analysis.velocity_limits([4.0, 2.0])


def test_velocity_limits_without_upper_limit():
    rng = np.random.default_rng(3)
    coordinates = rng.uniform(-50, 50, (10, 2))
    analysis = ArrayResponse(coordinates, grid_max=0.2,
                             grid_size=81).analyze()
    assert analysis.limits.aliasing_radius is None

    with pytest.raises(ValueError, match="upper"):
        analysis.velocity_limits([1.0])


def test_info(capsys):
    response = ArrayResponse(SEVEN)
    analysis = response.analyze()
    response.info()
    response.geometry.info()
    analysis.info()
    response.search_domain(0.03).info()
    text = capsys.readouterr().out
    assert "ArrayResponse" in text and "usable wavenumbers" in text
    assert "SearchDomainAnalysis" in text
