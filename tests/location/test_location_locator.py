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
"""Tests for absolute earthquake location and its weighted inversion."""

from __future__ import annotations

from math import dist, hypot, inf, sqrt
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.geodeticN import (
    MetricFrame,
    MetricPoint,
    WgsPoint,
    geo_to_ecef,
    metric_to_wgs,
    wgs_to_metric,
)
from shakelab.libutils.timeN import Date
from shakelab.location import (
    CartHypocenter,
    CartStation,
    HomogeneousEngine,
    Hypocenter,
    LayeredModel1D,
    LayeredTravelTimeEngine,
    PhasePick,
    Station,
    build_design_matrices,
    build_normal_equations,
    compute_covariance,
    geiger_locate,
    locate_event,
)

from shakelab.location.locator import prepare_picks_for_location
from shakelab.location.locutils import load_station_list
from shakelab.location.phases import (
    load_picks, normalize_phases, read_pun_origin,
)
from shakelab.location.ttsolver import read_1d_model_profile


@pytest.mark.parametrize("max_iter", [1, 20])
def test_locate_event_homogeneous(max_iter: int) -> None:
    """Recover an analytic source and report residuals at the final point."""
    truth = WgsPoint(13.05, 46.05, -8000.0)
    origin = Date("2024-03-27T21:19:37.123456789Z")
    initial = Hypocenter(WgsPoint(13.07, 46.04, -9000.0), origin + 0.5)
    ref_geo = WgsPoint(13.05, 46.05, 220.0)
    stations = {
        "A": Station("A", WgsPoint(13.00, 46.00, 200.0)),
        "B": Station("B", WgsPoint(13.00, 46.10, 220.0)),
        "C": Station("C", WgsPoint(13.10, 46.00, 210.0)),
        "D": Station("D", WgsPoint(13.10, 46.10, 230.0)),
        "E": Station("E", WgsPoint(13.05, 46.05, 250.0)),
    }
    velocities = {"P": 6000.0, "S": 3500.0}
    truth_xyz = geo_to_ecef(truth)
    picks = []
    for code, station in stations.items():
        length = dist(truth_xyz, geo_to_ecef(station.location))
        for phase, velocity in velocities.items():
            picks.append(PhasePick(code, phase, origin + length / velocity))

    estimated, result, used_ref = locate_event(
        picks,
        stations,
        HomogeneousEngine(vp=velocities["P"], vs=velocities["S"]),
        initial,
        ref_geo=ref_geo,
        damping=0.0,
        max_iter=max_iter,
        tol_space_m=1e-5,
        tol_time_s=1e-8,
    )

    assert used_ref == ref_geo
    assert result.hypocenter.point.frame == MetricFrame(ref_geo, "ecef")
    assert result.time_reference == initial.t0
    assert result.num_picks == len(picks)
    assert result.converged == (max_iter > 1)

    estimated_xyz = geo_to_ecef(estimated.location)
    expected = []
    for pick in sorted(picks, key=lambda pick: pick.time):
        station_xyz = geo_to_ecef(stations[pick.station].location)
        time = dist(estimated_xyz, station_xyz) / velocities[pick.phase]
        expected.append((pick.time - estimated.t0) - time)

    # Date rounds the reported origin to the nearest nanosecond.
    np.testing.assert_allclose(result.residuals, expected, rtol=0, atol=1e-9)
    assert result.rms == pytest.approx(
        float(np.sqrt(np.mean(np.square(expected)))), abs=1e-9,
    )

    if max_iter > 1:
        assert dist(estimated_xyz, truth_xyz) < 1e-4
        assert abs(estimated.t0 - origin) < 1e-8


def test_normal_equations_use_inverse_standard_deviation_weights():
    """Apply inverse-variance weights and add squared damping."""
    design = np.array([[1.0, 2.0], [3.0, 1.0], [1.0, 4.0]])
    residuals = np.array([1.0, -2.0, 3.0])
    weights = np.array([1.0, 4.0, 0.5])
    damping = 0.1
    normal, rhs, data = build_normal_equations(
        design, residuals, weights, damping,
    )
    inverse_variance = np.diag(weights**2)
    np.testing.assert_allclose(data, design.T @ inverse_variance @ design)
    np.testing.assert_allclose(rhs, design.T @ inverse_variance @ residuals)
    np.testing.assert_allclose(normal, data + damping**2 * np.eye(2))


@pytest.mark.parametrize("initial_depth", [2200.0, 3800.0])
def test_recovers_source_from_two_layer_first_arrivals(initial_depth):
    """Recover the known source using near and far P/S picks."""
    ref_geo = WgsPoint(13.05, 46.05, 350.0)
    frame = MetricFrame(ref_geo, "enu")
    model_top = -ref_geo.elevation

    source_east, source_north, source_depth = 1200.0, -800.0, 3000.0
    source = MetricPoint(
        source_east, source_north, model_top - source_depth, frame,
    )
    origin = Date("2024-03-27T21:19:37.123456789Z")
    initial = Hypocenter(
        metric_to_wgs(MetricPoint(500, 600, model_top - initial_depth, frame)),
        origin + 0.4,
    )

    thickness = 6000.0
    velocities = {"P": (4500.0, 6500.0), "S": (2500.0, 3650.0)}
    engine = LayeredTravelTimeEngine(
        LayeredModel1D([thickness, inf], velocities["P"], velocities["S"]),
        tol=1e-10,
    )

    # Endpoints are in the upper layer. Near stations see the direct
    # arrival first, distant stations see the critically refracted wave.
    receivers = {
        "N1": (3000.0, 2000.0, 100.0),
        "N2": (-5000.0, 4000.0, 200.0),
        "N3": (5500.0, -6500.0, 300.0),
        "N4": (-3000.0, -4000.0, 400.0),
        "F1": (38000.0, 13000.0, 250.0),
        "F2": (-29000.0, 24000.0, 200.0),
        "F3": (-18000.0, -40000.0, 100.0),
    }
    stations = {}
    picks = []
    for code, (east, north, depth) in receivers.items():
        location = metric_to_wgs(
            MetricPoint(east, north, model_top - depth, frame),
        )
        stations[code] = Station(code, location)

        offset = hypot(east - source_east, north - source_north)
        interface_legs = 2 * thickness - source_depth - depth
        for phase, (upper, lower) in velocities.items():
            direct_time = hypot(offset, source_depth - depth) / upper

            # Snell's law for critical refraction. A head wave exists
            # when both inclined legs fit inside the horizontal offset.
            critical_sin = upper / lower
            critical_cos = sqrt(1 - critical_sin**2)
            critical_offset = interface_legs * critical_sin / critical_cos
            head_time = inf
            if offset >= critical_offset:
                intercept_time = interface_legs * critical_cos / upper
                head_time = offset / lower + intercept_time

            assert (direct_time < head_time) == code.startswith("N")
            first_arrival = origin + min(direct_time, head_time)
            picks.append(PhasePick(code, phase, first_arrival))

    estimated, result, used_ref = locate_event(
        picks, stations, engine, initial, ref_geo=ref_geo,
        damping=0.0, max_iter=30, tol_space_m=1e-4, tol_time_s=1e-7,
    )

    assert used_ref == ref_geo
    assert result.converged
    assert result.num_picks == len(picks)
    assert result.hypocenter.point.frame == frame
    recovered = wgs_to_metric(estimated.location, frame)
    assert dist(
        (recovered.x, recovered.y, recovered.z),
        (source.x, source.y, source.z),
    ) < 0.01
    assert abs(estimated.t0 - origin) < 1e-6
    assert result.rms < 1e-6


@pytest.mark.parametrize(("require_time", "sort_by_time", "indices"), [
    (True, True, [3, 0]),
    (True, False, [0, 3]),
    (False, True, [2, 3, 0]),
    (False, False, [0, 2, 3]),
])
def test_prepare_picks_filters_and_sorts(require_time, sort_by_time, indices):
    """Filter unavailable stations and apply the requested time ordering."""
    origin = Date("2024-03-27T21:19:37Z")
    stations = {"STA": Station("STA", WgsPoint(13.0, 46.0, 100.0))}
    picks = [
        PhasePick("STA", "S", origin + 2),
        PhasePick("UNKNOWN", "P", origin),
        PhasePick("STA", "P", None),
        PhasePick("STA", "P", origin + 1),
    ]
    original = picks.copy()
    assert prepare_picks_for_location(
        picks, stations, require_time=require_time, sort_by_time=sort_by_time,
    ) == [original[index] for index in indices]
    assert picks == original
    with pytest.raises(ValueError, match="No valid picks"):
        prepare_picks_for_location(picks[1:3], stations)


def test_absolute_design_matrix_uses_relative_time_and_pick_weight():
    """Arrival residuals retain sub-microsecond offsets at modern dates."""
    origin = Date("2024-03-27T21:19:37.123456789Z")
    engine = HomogeneousEngine(5000, 2500)
    source = CartHypocenter(MetricPoint(0, 0, 0), 0.5)
    stations = {"STA": CartStation("STA", MetricPoint(3000, 4000, 0))}
    picks = [
        PhasePick("STA", "P", origin + 1.5000002, weight=2),
        PhasePick("STA", "S", origin + 2.5),
    ]
    design, residuals, weights = build_design_matrices(
        picks, stations, engine, source, time_reference=origin,
    )
    np.testing.assert_allclose(
        design, [[-0.6 / 5000, -0.8 / 5000, 0, 1],
                 [-0.6 / 2500, -0.8 / 2500, 0, 1]],
        rtol=0, atol=1e-15,
    )
    np.testing.assert_allclose(residuals, [2e-7, 0], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(weights, [2, 1])


def test_covariance_uses_weighted_residual_variance():
    """Scale the inverse undamped normal matrix by weighted RSS / dof."""
    design = np.array([[1., 0.], [1., 0.], [0., 1.], [0., 1.]])
    residuals = np.array([1., -1., 2., -2.])
    weights = np.array([2., 2., 1., 1.])
    normal = np.diag([8., 2.])
    covariance, errors = compute_covariance(design, residuals, weights, normal)
    np.testing.assert_allclose(covariance, [[1, 0], [0, 4]], rtol=0, atol=0)
    np.testing.assert_allclose(errors, [1, 2], rtol=0, atol=0)


def test_singular_absolute_location_reports_no_covariance():
    """An undamped rank-deficient system reports non-convergence."""
    origin = Date("2024-03-27T21:19:37Z")
    engine = HomogeneousEngine(5000, 2500)
    initial = CartHypocenter(MetricPoint(0, 0, 0), 0)
    stations = {"STA": CartStation("STA", MetricPoint(5000, 0, 0))}
    picks = [PhasePick("STA", "P", origin + 1.5)]
    result = geiger_locate(
        picks, stations, engine, initial, damping=0, time_reference=origin,
    )
    assert not result.converged
    assert result.cov_matrix is None
    assert result.num_iterations == 1
    np.testing.assert_allclose(result.residuals, [0.5], rtol=0, atol=1e-12)


@pytest.mark.parametrize(("model_name", "rms_s"), [
    ("homogeneous", 0.339),
    ("layered", 0.564),
])
def test_socchieve_relocation_regression(model_name, rms_s):
    """
    Reproduce the workspace comparison for 47 recorded P/S arrivals.

    The RMS is a regression reference at millisecond precision for the
    supplied models and picks. The PUN origin provides the initial
    geometry, with a deliberate spatial and temporal perturbation.
    """
    data = Path(__file__).resolve().parents[2] / "examples/location/data"
    event = data / "Picking" / "20240327211937"
    stations = load_station_list(data / "Stations" / "Station_list_OGS.csv")
    origin = read_pun_origin(event.with_suffix(".pun"))
    initial = Hypocenter(
        WgsPoint(
            origin.longitude - 0.04, origin.latitude + 0.03,
            origin.location.elevation - 2000,
        ),
        origin.time + 0.5,
    )
    reference = WgsPoint(origin.longitude, origin.latitude, 0.0)

    raw = normalize_phases(load_picks(event.with_suffix(".obs")))
    earliest = {}
    for pick in prepare_picks_for_location(raw, stations):
        if pick.phase in ("P", "S"):
            earliest.setdefault((pick.station, pick.phase), pick)
    picks = list(earliest.values())

    if model_name == "layered":
        thickness, vp, vs, _ = read_1d_model_profile(
            data / "Model" / "friuli_velocity_profile_custom.csv",
        )
        engine = LayeredTravelTimeEngine(LayeredModel1D(thickness, vp, vs))
    else:
        engine = HomogeneousEngine(5900, 3300)

    _, result, _ = locate_event(
        picks, stations, engine, initial, ref_geo=reference, damping=0.05,
        max_iter=30, tol_space_m=10, tol_time_s=0.005,
    )
    assert result.converged
    assert result.num_picks == 47
    assert np.isfinite(result.residuals).all()
    assert result.rms == pytest.approx(rms_s, abs=0.001)
