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
Tests for 1D ray geometry and travel times.

Straight-ray distances, Snell's law and image-source solutions provide
independent reference values.
"""

from __future__ import annotations

from math import hypot, inf, sin, sqrt

import numpy as np
import pytest

from shakelab.location import ttsolver as tt


@pytest.fixture
def layered_model():
    """Three layers with increasing velocities above a half-space."""
    return tt.LayeredModel1D(
        [2000, 2000, inf], [4000, 5000, 6000], [2400, 3000, 3600],
    )


def _assert_path(ray, layers, depths):
    """Check path continuity, segment travel times and slowness."""
    segments = ray["segments"]
    assert [segment.layer for segment in segments] == layers
    np.testing.assert_allclose(
        [segment.z_start for segment in segments], depths[:-1],
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        [segment.z_end for segment in segments], depths[1:],
        rtol=0, atol=1e-10,
    )
    for segment in segments:
        assert segment.p**2 + segment.q**2 == pytest.approx(
            1 / segment.v**2, rel=1e-12, abs=0,
        )
    assert sum(segment.t for segment in segments) == pytest.approx(
        ray["t"], rel=1e-12, abs=0,
    )


def test_model_depths_and_interfaces():
    """Layer indices are right-continuous at internal interfaces."""
    model = tt.LayeredModel1D([10, 20], [4000, 5000], [2300, 2800])

    assert model.num_layers == 2
    assert model.total_thickness() == 30
    np.testing.assert_array_equal(model.z_top, [0, 10])
    np.testing.assert_array_equal(model.z_bottom, [10, 30])
    assert [model.layer_at(depth) for depth in (0, 9.999, 10, 29.999)] == [
        0, 0, 1, 1,
    ]
    for depth in (-1, 30):
        with pytest.raises(ValueError):
            model.layer_at(depth)


@pytest.mark.parametrize(("thickness", "vp", "vs"), [
    ([10, 20], [4000], [2300, 2800]),
    ([10, -5], [4000, 5000], [2300, 2800]),
    ([10, 20], [0, 5000], [2300, 2800]),
])
def test_model_rejects_invalid_layers(thickness, vp, vs):
    """Layer arrays must agree in length and contain positive values."""
    with pytest.raises(ValueError):
        tt.LayeredModel1D(thickness, vp, vs)


def test_load_velocity_profile(tmp_path):
    """Read SI columns in their documented order, including density."""
    path = tmp_path / "profile.csv"
    path.write_text(
        "thickness (m),vp (m/s),vs (m/s),density (kg/m3)\n"
        "1500,4000,2300,2400\n2500,5000,3000,2600\n",
        encoding="ascii",
    )
    assert tt.read_1d_model_profile(path) == (
        [1500, 2500], [4000, 5000], [2300, 3000], [2400, 2600],
    )


@pytest.mark.parametrize(("phase", "velocity"), [("P", 4000), ("S", 2300)])
def test_direct_traveltime_in_a_half_space(phase, velocity):
    """T equals the straight-line distance divided by phase velocity."""
    model = tt.LayeredModel1D([inf], [4000], [2300])
    offsets = np.array([0, 500, 1000, 2000])
    expected = np.hypot(offsets, 1000) / velocity
    np.testing.assert_allclose(
        tt.tt_direct_curve(model, offsets, 1000, phase=phase),
        expected, rtol=1e-12, atol=0,
    )
    assert tt.tt_direct_surface(
        model, 500, 1000, phase=phase,
    ) == pytest.approx(expected[1], rel=1e-12)


@pytest.mark.parametrize("offset", [0.0, 3500.0])
def test_direct_path_crosses_layers_in_both_directions(layered_model, offset):
    """Reciprocity reverses the path and preserves its travel time."""
    down = tt.direct_path(layered_model, offset, 2500, 4700)
    up = tt.direct_path(layered_model, offset, 4700, 2500)
    _assert_path(down, [1, 2], [2500, 4000, 4700])
    _assert_path(up, [2, 1], [4700, 4000, 2500])
    assert down["t"] == pytest.approx(up["t"], abs=1e-12)
    assert sin(down["theta_source"]) == pytest.approx(
        down["p"] * 5000, abs=1e-12,
    )
    if offset == 0:
        assert up["t"] == pytest.approx(700 / 6000 + 1500 / 5000)


def test_vertical_path_through_three_layers():
    """A vertical path sums each layer's thickness divided by velocity."""
    model = tt.LayeredModel1D([10, 20, 30], [4000, 5000, 6000],
                             [2300, 2800, 3300])
    ray = tt.direct_path(model, 0, 5, 45)
    _assert_path(ray, [0, 1, 2], [5, 10, 30, 45])
    assert ray["t"] == pytest.approx(5 / 4000 + 20 / 5000 + 15 / 6000)


def test_head_wave_matches_critical_refraction(layered_model):
    """Check the critical slowness, intercept time and head-wave path."""
    offset = 20000
    times, branches = tt.solve_traveltime(
        layered_model, [offset], 2500, 3000, return_branches=True,
    )
    ray = min(branches[0], key=lambda branch: branch["t"])
    expected = offset / 6000 + 2500 * sqrt(5000**-2 - 6000**-2)
    assert ray["kind"] == "head"
    assert ray["p"] == pytest.approx(1 / 6000)
    assert times[0] == pytest.approx(expected, rel=1e-12)
    _assert_path(ray, [1, 2, 1], [2500, 4000, 4000, 3000])


@pytest.mark.parametrize(("phase", "slow", "fast"), [
    ("P", 4000, 6000), ("S", 2300, 3500),
])
def test_first_arrival_is_continuous_at_interface(phase, slow, fast):
    """Check the critical limit and reciprocity on both sides."""
    model = tt.LayeredModel1D([4000, inf], [4000, 6000], [2300, 3500])
    offset = 50000
    limit = offset / fast + 4000 * sqrt(slow**-2 - fast**-2)

    for depth, kind in ((3999.99, "head"), (4000, "head"),
                        (4000.01, "direct")):
        times = []
        for source, receiver in ((depth, 0), (0, depth)):
            time, branches = tt.solve_traveltime(
                model, offset, source, receiver, phase=phase,
                include_diving=False, return_branches=True,
            )
            ray = min(branches[0], key=lambda branch: branch["t"])
            assert ray["kind"] == kind
            assert float(time) == pytest.approx(
                limit, rel=0, abs=1e-11 if depth == 4000 else 1e-5,
            )
            times.append(float(time))
        assert times[0] == pytest.approx(times[1], rel=0, abs=1e-12)


@pytest.mark.parametrize(("phase", "slow", "fast"), [
    ("P", 4000, 6000), ("S", 2300, 3500),
])
@pytest.mark.parametrize("depth", [4000.1, 4000.000001])
def test_nearly_horizontal_direct_ray_matches_snell(phase, slow, fast, depth):
    """Prescribe the fast-layer segment and recover both ray legs."""
    model = tt.LayeredModel1D([4000, inf], [4000, 6000], [2300, 3500])
    fast_offset = 45000
    fast_length = hypot(fast_offset, depth - 4000)
    sin_slow = (slow / fast) * fast_offset / fast_length
    cos_slow = sqrt(1 - sin_slow**2)
    offset = fast_offset + 4000 * sin_slow / cos_slow
    expected = fast_length / fast + 4000 / (slow * cos_slow)

    ray = tt.direct_path(model, offset, depth, 0, phase=phase, tol=1e-10)
    _assert_path(ray, [1, 0], [depth, 4000, 0])
    assert ray["t"] == pytest.approx(expected, rel=0, abs=3e-9)
    assert ray["segments"][0].q == pytest.approx(
        (depth - 4000) / (fast * fast_length), rel=2e-10, abs=0,
    )


def test_same_depth_curve_switches_from_horizontal_to_head(layered_model):
    """The horizontal direct ray wins nearby and the head wave farther out."""
    offsets = np.array([0, 5000, 30000])
    far_head = 30000 / 6000 + 3000 * sqrt(5000**-2 - 6000**-2)
    np.testing.assert_allclose(
        tt.same_depth_curve(layered_model, offsets, 2500),
        [0, 1, far_head], rtol=1e-12, atol=0,
    )
    assert tt.same_depth_best(
        layered_model, 5000, 2500,
    )["kind"] == "direct_horizontal"
    assert tt.same_depth_best(layered_model, 30000, 2500)["kind"] == "head"


def test_traveltime_preserves_offset_shape(layered_model):
    """Return times in the input shape and branch lists in offset order."""
    offsets = np.array([[0, 1000], [2000, 3000]])
    times, branches = tt.solve_traveltime(
        layered_model, offsets, 2500, 3000, arrivals="direct",
        return_branches=True,
    )
    assert times.shape == (2, 2)
    np.testing.assert_allclose(
        times, np.hypot(offsets, 500) / 5000, rtol=1e-12, atol=0,
    )
    assert [row[0]["x_target"] for row in branches] == list(offsets.flat)
    assert all(len(row) == 1 and row[0]["kind"] == "direct"
               for row in branches)


@pytest.mark.parametrize("offset", [0.0, 3000.0])
def test_legacy_reflection_matches_two_straight_legs(layered_model, offset):
    """The legacy branch is the reflection from the deeper interface."""
    branches = tt.compute_branches(layered_model, offset, 2500, 3000)
    ray = next(branch for branch in branches if branch["kind"] == "diving")
    assert ray["t"] == pytest.approx(hypot(offset, 2500) / 5000, abs=1e-6)
    _assert_path(ray, [1, 1], [2500, 4000, 3000])


@pytest.mark.parametrize(("phase", "velocity"), [("P", 5000), ("S", 3000)])
@pytest.mark.parametrize("offset", [0.0, 3000.0])
def test_reflection_at_equal_depth(layered_model, phase, velocity, offset):
    """A flat reflector yields the travel time to the image source."""
    ray = tt.reflected_path(
        layered_model, offset, 2500, 2500,
        interface_index=2, phase=phase, tol=1e-10,
    )
    assert ray["kind"] == "reflected"
    assert ray["reflection_depth"] == 4000
    assert ray["t"] == pytest.approx(hypot(offset, 3000) / velocity, abs=1e-9)
    _assert_path(ray, [1, 1], [2500, 4000, 2500])
    assert ray["segments"][0].x == pytest.approx(offset / 2, abs=1e-6)


@pytest.mark.parametrize(("offset", "source_depth", "index"), [
    (1000, 0, 0), (1000, 4000, 2), (-1000, 0, 2),
])
def test_reflection_rejects_invalid_geometry(
    layered_model, offset, source_depth, index,
):
    """Reject invalid offsets, interfaces and source positions."""
    with pytest.raises(ValueError):
        tt.reflected_path(
            layered_model, offset, source_depth, interface_index=index,
        )


def test_reflection_crosses_layers_with_one_ray_parameter(layered_model):
    """Recover an independently prescribed Snell parameter and path."""
    p = 1 / 10000
    # Source is in the top layer, receiver and reflector are deeper.
    heights = np.array([1000, 2000, 1500])
    speeds = np.array([4000, 5000, 5000])
    cosines = np.sqrt(1 - (p * speeds)**2)
    offset = float(np.sum(heights * p * speeds / cosines))
    expected_time = float(np.sum(heights / (speeds * cosines)))

    ray = tt.reflected_path(
        layered_model, offset, 1000, 2500, interface_index=2, tol=1e-10,
    )
    _assert_path(ray, [0, 1, 1], [1000, 2000, 4000, 2500])
    assert ray["p"] == pytest.approx(p, rel=1e-9)
    assert ray["t"] == pytest.approx(expected_time, rel=1e-9)


def test_first_arrival_curve_with_heads(layered_model):
    """Check the direct-to-head crossover and the critical offset."""
    offsets = [0, 1000, 20000]
    head = 20000 / 6000 + 2500 * sqrt(5000**-2 - 6000**-2)
    np.testing.assert_allclose(
        tt.curve_with_heads(layered_model, offsets, 2500, 3000,
                            include_diving=False),
        [0.1, hypot(1000, 500) / 5000, head], rtol=1e-12, atol=0,
    )
    branches = tt.compute_branches(
        layered_model, 1000, 2500, 3000, include_diving=False,
    )
    assert [branch["kind"] for branch in branches] == ["direct"]
