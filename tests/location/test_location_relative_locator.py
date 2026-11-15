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
"""Check the relative locator with independently calculated arrivals."""

from __future__ import annotations

from itertools import combinations
from math import dist, inf

import numpy as np
import pytest

from shakelab.libutils.geodeticN import (
    MetricFrame, MetricPoint, WgsPoint, wgs_to_metric,
)
from shakelab.libutils.timeN import Date
from shakelab.location import (
    CartHypocenter, CartStation, HomogeneousEngine, Hypocenter, LayeredModel1D,
    LayeredTravelTimeEngine, PhasePick, Station,
)
from shakelab.location.relative_locator import (
    RelativePhaseObservation, build_relative_design_matrices,
    estimate_origin_time, locate_event_relative,
)


@pytest.mark.parametrize("layered", [False, True])
def test_relative_location_and_origin_time(layered):
    """Recover the hypocenter from differences and t0 from arrivals."""
    ref = WgsPoint(13.05, 46.05, 0.0)
    truth = WgsPoint(13.05, 46.05, -8000.0)
    origin = Date("2024-03-27T21:19:37.123456789Z")
    stations = {
        str(i): Station(str(i), WgsPoint(lon, lat, -100))
        for i, (lon, lat) in enumerate([
            (13.00, 46.00), (13.10, 46.00), (13.00, 46.10),
            (13.10, 46.10), (13.05, 46.05),
        ])
    }
    initial = Hypocenter(WgsPoint(13.06, 46.04, -9000), Date())
    velocities = {"P": 6000.0, "S": 3500.0}
    if layered:
        engine = LayeredTravelTimeEngine(
            LayeredModel1D([inf], [velocities["P"]], [velocities["S"]]),
        )
    else:
        engine = HomogeneousEngine(velocities["P"], velocities["S"])

    frame = MetricFrame(ref, engine.metric_frame)
    source_point = wgs_to_metric(truth, frame)
    cart_stations = {}
    times = {}
    for code, station in stations.items():
        point = wgs_to_metric(station.location, frame)
        depth = -point.z if layered else 100
        cart_stations[code] = CartStation(code, point, depth)
        length = dist(
            (point.x, point.y, point.z),
            (source_point.x, source_point.y, source_point.z),
        )
        for phase, velocity in velocities.items():
            times[code, phase] = origin + length / velocity
    picks = [PhasePick(code, phase, t) for (code, phase), t in times.items()]
    observations = [
        RelativePhaseObservation(
            i, j, phase, times[i, phase] - times[j, phase],
        )
        for i, j in combinations(stations, 2)
        for phase in velocities
    ]
    located, result, used_ref = locate_event_relative(
        observations, stations, engine, initial, ref_geo=ref,
        damping=0.0, tol_space_m=1e-4,
    )
    assert result.converged
    assert result.num_observations == len(observations)
    assert result.rms < 1e-8
    assert used_ref == ref
    point = wgs_to_metric(located.location, frame)
    assert dist((point.x, point.y, point.z),
                (source_point.x, source_point.y, source_point.z)) < 1e-3
    time_result = estimate_origin_time(
        picks, cart_stations, engine, result.hypocenter,
    )
    assert abs(time_result.t0 - origin) < 1e-8
    assert time_result.rms < 1e-8


def test_relative_design_matrix_uses_station_order():
    """Use station-i minus station-j for travel times and derivatives."""
    source = CartHypocenter(MetricPoint(0, 0, 0), 0)
    engine = HomogeneousEngine(5000, 2500)
    stations = {
        "A": CartStation("A", MetricPoint(5000, 0, 0)),
        "B": CartStation("B", MetricPoint(0, 10000, 0)),
    }
    observations = [
        RelativePhaseObservation("A", "B", "P", -0.75, weight=2),
        RelativePhaseObservation("B", "A", "P", 0.75, weight=3),
    ]
    design, residuals, weights = build_relative_design_matrices(
        observations, stations, engine, source,
    )
    np.testing.assert_allclose(
        design, [[-1 / 5000, 1 / 5000, 0], [1 / 5000, -1 / 5000, 0]],
        rtol=0, atol=1e-15,
    )
    np.testing.assert_allclose(residuals, [0.25, -0.25], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(weights, [2, 3])


def test_origin_time_uses_squared_weights():
    """Estimate the weighted mean of the origin inferred by each pick."""
    reference = Date("2024-03-27T21:19:37.123456789Z")
    engine = HomogeneousEngine(5000, 2500)
    source = CartHypocenter(MetricPoint(0, 0, 0), 0)
    stations = {"STA": CartStation("STA", MetricPoint(5000, 0, 0))}
    picks = [
        PhasePick("STA", "P", reference + 1),
        PhasePick("STA", "S", reference + 3, weight=2),
    ]
    result = estimate_origin_time(picks, stations, engine, source)
    assert result.t0 == reference + 0.8
    np.testing.assert_allclose(
        result.residuals, [-0.8, 0.2], rtol=0, atol=1e-12,
    )
    assert result.weighted_rms == pytest.approx(np.hypot(0.8, 0.4) / 2**0.5)


def test_origin_time_rejects_missing_times_and_zero_weights():
    """Require observed times and at least one positive squared weight."""
    reference = Date("2024-03-27T21:19:37Z")
    engine = HomogeneousEngine(5000, 2500)
    source = CartHypocenter(MetricPoint(0, 0, 0), 0)
    stations = {"STA": CartStation("STA", MetricPoint(5000, 0, 0))}
    with pytest.raises(ValueError, match="timed picks"):
        estimate_origin_time(
            [PhasePick("STA", "P", None)], stations, engine, source,
        )
    with pytest.raises(ValueError, match="sum of weights"):
        estimate_origin_time(
            [PhasePick("STA", "P", reference, weight=0)],
            stations, engine, source,
        )
