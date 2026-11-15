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
"""Tests for the Cartesian travel-time engines and their gradients."""

from __future__ import annotations

from math import hypot, inf, sqrt

import numpy as np
import pytest

from shakelab.libutils.geodeticN import MetricFrame, MetricPoint
from shakelab.location import (
    CartHypocenter,
    CartStation,
    HomogeneousEngine,
    LayeredModel1D,
    LayeredTravelTimeEngine,
    ttsolver,
)


@pytest.mark.parametrize("receiver_depth", [0.0, 3000.0, 5000.0])
@pytest.mark.parametrize(("phase", "velocity"), [("P", 6000), ("S", 3500)])
def test_one_layer_traveltime_and_gradients(receiver_depth, phase, velocity):
    """Check the analytic gradient and its finite-difference estimate."""
    frame = MetricFrame((13.0, 46.0, 800.0), "enu")
    engine = LayeredTravelTimeEngine(
        LayeredModel1D([inf], [6000], [3500]), fd_step=1.0,
    )
    source = CartHypocenter(MetricPoint(0, 0, -3800, frame), 0)
    receiver = CartStation(
        "A", MetricPoint(4000, 1000, -receiver_depth - 800, frame),
        receiver_depth,
    )
    delta_z = receiver_depth - 3000
    length = hypot(4000, 1000, delta_z)
    gradient = np.array([-4000, -1000, delta_z]) / (length * velocity)
    assert engine.traveltime(phase, source, receiver) == pytest.approx(
        length / velocity, abs=1e-12,
    )
    np.testing.assert_allclose(
        engine.gradients(phase, source, receiver), gradient, atol=1e-12,
    )
    np.testing.assert_allclose(
        engine.gradients_numerical(phase, source, receiver),
        gradient, atol=1e-10,
    )

@pytest.mark.parametrize(("phase", "velocity"), [("p", 6000), ("S", 3500)])
def test_homogeneous_traveltime_and_gradient(phase, velocity):
    """Differentiate the straight-line travel time with respect to source."""
    engine = HomogeneousEngine(6000, 3500)
    source = CartHypocenter(MetricPoint(0, 0, 0), 0)
    station = CartStation("STA", MetricPoint(3000, 4000, 12000))
    assert engine.traveltime(phase, source, station) == pytest.approx(
        13000 / velocity, rel=1e-12,
    )
    np.testing.assert_allclose(
        engine.gradients(phase, source, station),
        -np.array([3000, 4000, 12000]) / (13000 * velocity),
        rtol=1e-12, atol=0,
    )


def test_homogeneous_limits_and_errors():
    """Return zero at coincidence and reject unknown phases and speeds."""
    engine = HomogeneousEngine(6000, 3500)
    point = MetricPoint(1, 2, 3)
    source = CartHypocenter(point, 0)
    station = CartStation("STA", point)
    assert engine.traveltime("P", source, station) == 0
    assert engine.gradients("P", source, station) == (0, 0, 0)
    with pytest.raises(ValueError, match="Unsupported phase"):
        engine.traveltime("Lg", source, station)
    with pytest.raises(ValueError, match="positive"):
        HomogeneousEngine(6000, 0)


def test_first_arrival_head_wave_omits_legacy_reflections(monkeypatch):
    """Use the head wave and its gradient without tracing a reflection."""
    model = LayeredModel1D([4000, inf], [5000, 6000], [3000, 3600])
    engine = LayeredTravelTimeEngine(model, fd_step=1)
    source = CartHypocenter(MetricPoint(0, 0, -2500), 0)
    station = CartStation("STA", MetricPoint(20000, 0, -3000), 3000)

    def fail_on_reflection(*args, **kwargs):
        pytest.fail("The first-arrival engine requested a reflection.")

    monkeypatch.setattr(
        ttsolver, "_solve_reflected_branch", fail_on_reflection,
    )
    vertical_slowness = sqrt(5000**-2 - 6000**-2)
    assert engine.traveltime("P", source, station) == pytest.approx(
        20000 / 6000 + 2500 * vertical_slowness, rel=1e-12,
    )
    gradient = [-1 / 6000, 0, vertical_slowness]
    np.testing.assert_allclose(
        engine.gradients("P", source, station), gradient, rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        engine.gradients_numerical("P", source, station),
        gradient, rtol=0, atol=1e-10,
    )


def test_layered_gradient_falls_back_to_finite_differences(monkeypatch):
    """Recover the derivative when the analytical branch is unavailable."""
    model = LayeredModel1D([inf], [6000], [3500])
    engine = LayeredTravelTimeEngine(model, fd_step=1)
    source = CartHypocenter(MetricPoint(0, 0, -4000), 0)
    station = CartStation("STA", MetricPoint(3000, 0, 0))
    monkeypatch.setattr(
        engine, "gradients_analytical", lambda *args: (np.nan,) * 3,
    )
    np.testing.assert_allclose(
        engine.gradients("P", source, station),
        [-0.6 / 6000, 0, -0.8 / 6000], rtol=0, atol=1e-10,
    )


@pytest.mark.parametrize("depth", [4000.0, 4000.000001])
def test_layered_gradient_at_grazing_takeoff(depth):
    """Use the one-sided ray slowness on and below the interface."""
    model = LayeredModel1D([4000, inf], [4000, 6000], [2300, 3500])
    engine = LayeredTravelTimeEngine(model, tol=1e-10)
    source = CartHypocenter(MetricPoint(0, 0, -depth), 0)
    fast_offset = 45000
    fast_length = hypot(fast_offset, depth - 4000)
    sin_slow = (4000 / 6000) * fast_offset / fast_length
    offset = fast_offset + 4000 * sin_slow / sqrt(1 - sin_slow**2)
    station = CartStation("STA", MetricPoint(offset, 0, 0))

    dtx, dty, dtz = engine.gradients("P", source, station)
    assert dtx == pytest.approx(-fast_offset / (6000 * fast_length))
    assert dty == 0
    assert dtz == pytest.approx(
        -(depth - 4000) / (6000 * fast_length), rel=2e-10, abs=0,
    )


def test_layered_model_bounds_and_frame():
    """The ENU datum sets model depth and both endpoints share one frame."""
    model = LayeredModel1D([5000], [6000], [3500])
    engine = LayeredTravelTimeEngine(model)
    frame = MetricFrame((13.0, 46.0, 200.0), "enu")
    source = CartHypocenter(MetricPoint(0, 0, -200, frame), 0)
    assert engine.contains_source(source)
    assert not engine.contains_source(
        CartHypocenter(MetricPoint(0, 0, -5200, frame), 0),
    )
    assert not engine.contains_source(
        CartHypocenter(MetricPoint(0, 0, -199, frame), 0),
    )
    station = CartStation("STA", MetricPoint(2000, 0, 0))
    with pytest.raises(ValueError, match="common ENU frame"):
        engine.traveltime("P", source, station)


def test_layered_gradient_method_validation():
    """The engine rejects unsupported gradient and step options."""
    model = LayeredModel1D([inf], [6000], [3500])
    source = CartHypocenter(MetricPoint(0, 0, -4000), 0)
    station = CartStation("STA", MetricPoint(3000, 0, 0))
    engine = LayeredTravelTimeEngine(model, gradient_method="unknown")
    with pytest.raises(ValueError, match="Unsupported gradient_method"):
        engine.gradients("P", source, station)

    engine = LayeredTravelTimeEngine(model, gradient_method="numerical",
                                     fd_step=0)
    with pytest.raises(ValueError, match="fd_step must be positive"):
        engine.gradients("P", source, station)
