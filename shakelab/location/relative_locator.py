# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
#
# ****************************************************************************
"""
Geiger hypocenter location from relative travel times.

This module implements a damped, weighted Geiger location algorithm
that operates on *relative* arrival times between station pairs for
one or more seismic phases (typically P and S).

The key idea is that for an event with origin time t0:

    t_i^obs = t0 + T_i(x)

where T_i(x) is the travel time from the source x to station i.
The time difference between two stations i and j is:

    Δt_ij^obs = t_i^obs - t_j^obs = T_i(x) - T_j(x)

The origin time t0 cancels out, so relative times constrain only
the source position (x, y, z), not the origin time. A separate
helper is provided to estimate t0 afterwards from absolute picks.

The core inversion operates in local Cartesian coordinates
(`CartStation`, `CartHypocenter`). A high-level wrapper
(:func:`locate_event_relative`) is provided to work directly with
geographic coordinates through the user-facing `Station` and
`Hypocenter` classes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence, Optional, Tuple

import numpy as np

from shakelab.location.locator import (  # type: ignore
    Station,
    Hypocenter,
    CartStation,
    CartHypocenter,
    TravelTimeEngine,
    build_normal_equations,
    compute_covariance,
    _auto_ref_geo,
)
from shakelab.libutils.geodeticN import (
    MetricFrame, MetricPoint, WgsPoint, metric_to_wgs, wgs_to_metric,
)
from shakelab.libutils.timeN import Date
from shakelab.location.phases import PhasePick


# ---------------------------------------------------------------------------
# Data structures for relative-time observations
# ---------------------------------------------------------------------------


@dataclass
class RelativePhaseObservation:
    """Relative arrival time between a pair of stations for one phase.

    Parameters
    ----------
    station_i : str
        Identifier of the first station (i).
    station_j : str
        Identifier of the second station (j).
    phase : str
        Seismic phase identifier (e.g. "P", "S").
    dt_obs : float
        Observed time difference Δt_ij = t_i - t_j [s].
        The sign convention must be consistent across all observations.
    weight : float, optional
        Data weight. A common choice is weight = 1 / sigma, where sigma
        is the expected uncertainty on Δt_ij [s].
    """

    station_i: str
    station_j: str
    phase: str
    dt_obs: float
    weight: float = 1.0


@dataclass
class RelativeLocationResult:
    """Result of a relative-time hypocenter inversion (x, y, z only)."""

    hypocenter: CartHypocenter
    covariance: Optional[np.ndarray]
    std_x: Optional[float]
    std_y: Optional[float]
    std_z: Optional[float]
    rms: float
    weighted_rms: float
    num_iterations: int
    num_observations: int
    converged: bool
    residuals: np.ndarray


@dataclass
class OriginTimeResult:
    """Result of a separate origin-time estimation."""

    t0: Date
    rms: float
    weighted_rms: float
    residuals: np.ndarray


# ---------------------------------------------------------------------------
# Design matrix for relative times (Cartesian)
# ---------------------------------------------------------------------------


def build_relative_design_matrices(
    obs: Sequence[RelativePhaseObservation],
    stations: Dict[str, CartStation],
    engine: TravelTimeEngine,
    source: CartHypocenter,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build design matrix G, residual vector r and weights w for Δt.

    For each observation Δt_ij^obs = t_i - t_j:

        Δt_ij^pred = T_i(x) - T_j(x)

    Linearizing around the current source x_k and keeping only spatial
    derivatives:

        r_ij = Δt_ij^obs - Δt_ij^pred
             ≈ (∇T_i - ∇T_j) · δx

    Parameters
    ----------
    obs : sequence of RelativePhaseObservation
        Relative-time observations between station pairs.
    stations : dict
        Mapping station_id -> CartStation.
    engine : TravelTimeEngine
        Travel-time engine providing traveltimes and gradients.
    source : CartHypocenter
        Current hypocenter estimate (x, y, z). The t0 field is ignored.

    Returns
    -------
    G : ndarray, shape (N, 3)
        Design matrix with columns [dΔT/dx, dΔT/dy, dΔT/dz].
    r : ndarray, shape (N,)
        Residual vector: observed - predicted time differences [s].
    w : ndarray, shape (N,)
        Data weights.
    """
    n = len(obs)
    G = np.zeros((n, 3), dtype=float)
    r = np.zeros(n, dtype=float)
    w = np.zeros(n, dtype=float)

    for k, ob in enumerate(obs):
        sta_i = stations[ob.station_i]
        sta_j = stations[ob.station_j]

        # Predicted traveltimes from current source
        ti = engine.traveltime(ob.phase, source, sta_i)
        tj = engine.traveltime(ob.phase, source, sta_j)
        dt_pred = ti - tj

        # Gradients of traveltime at i and j
        dti_dx, dti_dy, dti_dz = engine.gradients(
            ob.phase,
            source,
            sta_i,
        )
        dtj_dx, dtj_dy, dtj_dz = engine.gradients(
            ob.phase,
            source,
            sta_j,
        )

        # Gradient of ΔT = T_i - T_j
        G[k, 0] = dti_dx - dtj_dx
        G[k, 1] = dti_dy - dtj_dy
        G[k, 2] = dti_dz - dtj_dz

        r[k] = ob.dt_obs - dt_pred
        w[k] = ob.weight

    return G, r, w


# ---------------------------------------------------------------------------
# Geiger-style inversion for relative times (Cartesian)
# ---------------------------------------------------------------------------


def geiger_locate_relative(
    obs: Sequence[RelativePhaseObservation],
    stations: Dict[str, CartStation],
    engine: TravelTimeEngine,
    initial: CartHypocenter,
    damping: float = 0.1,
    max_iter: int = 20,
    tol_space: float = 1.0,
    space_scale: float = 1000.0,
) -> RelativeLocationResult:
    """Iterative Geiger location using only relative times (Δt).

    This solver estimates only the source position (x, y, z). The
    origin time t0 is not constrained and is left as a dummy value
    in the returned CartHypocenter (set to 0.0).

    Parameters
    ----------
    obs : sequence of RelativePhaseObservation
        Relative-time observations between station pairs.
    stations : dict
        Mapping station_id -> CartStation (Cartesian coordinates).
    engine : TravelTimeEngine
        Travel-time engine providing traveltimes and gradients.
    initial : CartHypocenter
        Initial hypocenter estimate (x, y, z). The t0 field is ignored.
    damping : float, optional
        Damping parameter for the least-squares update.
    max_iter : int, optional
        Maximum number of iterations.
    tol_space : float, optional
        Spatial convergence threshold on the update norm [same units
        as coordinates].
    space_scale : float, optional
        Scaling factor applied to spatial parameters when forming the
        design matrix. The default 1000.0 keeps solver magnitudes in
        line with the previous km-based implementation while internal
        coordinates are stored in meters.

    Returns
    -------
    RelativeLocationResult
        Object containing the final hypocenter (x, y, z), covariance
        (if available), residuals and iteration statistics.
    """
    m = np.array([initial.point.x, initial.point.y, initial.point.z])
    converged = False
    iterations = 0

    for iterations in range(1, max_iter + 1):
        source = CartHypocenter(
            MetricPoint(*m, frame=initial.point.frame), t0=initial.t0,
        )
        G, r, w = build_relative_design_matrices(obs, stations, engine, source)
        if not (np.all(np.isfinite(G)) and np.all(np.isfinite(r))):
            break
        A, b, _ = build_normal_equations(G * space_scale, r, w, damping)

        try:
            dm = space_scale * np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(dm)):
            break
        m += dm
        if np.linalg.norm(dm) < tol_space:
            converged = True
            break

    hypocenter = CartHypocenter(
        MetricPoint(*m, frame=initial.point.frame), t0=initial.t0,
    )
    G, residuals, weights = build_relative_design_matrices(
        obs, stations, engine, hypocenter,
    )
    cov = stds = None
    if np.all(np.isfinite(G)) and np.all(np.isfinite(residuals)):
        weighted_G = G * weights[:, None]
        cov, stds = compute_covariance(
            G, residuals, weights, weighted_G.T @ weighted_G,
        )
    num_obs = residuals.size
    rms = float(np.sqrt(np.mean(residuals**2))) if num_obs else 0.0
    wrms = float(np.sqrt(np.mean((weights * residuals)**2))) if num_obs else 0.0

    if stds is not None and stds.size >= 3:
        std_x = float(stds[0])
        std_y = float(stds[1])
        std_z = float(stds[2])
    else:
        std_x = std_y = std_z = None

    result = RelativeLocationResult(
        hypocenter=hypocenter,
        covariance=cov,
        std_x=std_x,
        std_y=std_y,
        std_z=std_z,
        rms=rms,
        weighted_rms=wrms,
        num_iterations=iterations,
        num_observations=int(num_obs),
        converged=converged,
        residuals=residuals,
    )

    return result


# ---------------------------------------------------------------------------
# Geographic wrapper (analogous to locate_event, but without t0)
# ---------------------------------------------------------------------------


def locate_event_relative(
    obs: Sequence[RelativePhaseObservation],
    stations: Dict[str, Station],
    engine: TravelTimeEngine,
    initial: Hypocenter,
    ref_geo: Optional[WgsPoint] = None,
    damping: float = 0.1,
    max_iter: int = 20,
    tol_space_m: float = 1.0,
) -> Tuple[Hypocenter, RelativeLocationResult, WgsPoint]:
    """Locate an event from relative times in geographic coordinates.

    This is the relative-time analogue of :func:`locator.locate_event`,
    but it solves only for the spatial coordinates (lat, lon, depth).
    The origin time in the returned Hypocenter is an unset Date.

    Parameters
    ----------
    obs : sequence of RelativePhaseObservation
        Relative-time observations between station pairs.
    stations : dict
        Mapping station_id -> Station (geographic).
    engine : TravelTimeEngine
        Travel-time engine compatible with the Cartesian locator.
    initial : Hypocenter
        Initial hypocenter guess in geographic coordinates.
    ref_geo : WgsPoint, optional
        Reference geodetic point (latitude [deg], longitude [deg],
        height [m, positive upward]) used to define the local Cartesian
        system. If None, it is inferred from stations.
    damping : float, optional
        Damping factor for the Geiger iterations.
    max_iter : int, optional
        Maximum number of iterations.
    tol_space_m : float, optional
        Spatial convergence tolerance in meters.

    Returns
    -------
    hypo_geo : Hypocenter
        Estimated hypocenter in geographic coordinates. The t0 field
        is unset because relative times do not constrain it.
    result : RelativeLocationResult
        Full result of the Cartesian relative inversion.
    ref_geo_out : WgsPoint
        Reference geodetic point actually used for the conversion.
    """
    if not stations:
        raise ValueError("At least one station is required for location.")

    if ref_geo is None:
        ref_geo = _auto_ref_geo(stations=stations, initial=initial)

    orientation = getattr(engine, "metric_frame", "ecef")
    frame = MetricFrame(ref_geo, orientation)
    cart_stations = {}
    for st_id, station in stations.items():
        point = wgs_to_metric(station.location, frame)
        elevation = station.location.elevation
        if frame.orientation == "enu":
            elevation = point.z + frame.ref_geo[2]
        cart_stations[st_id] = CartStation(st_id, point, max(0.0, -elevation))

    initial_cart = CartHypocenter(wgs_to_metric(initial.location, frame), 0.0)

    # Run the Cartesian relative-time locator
    result = geiger_locate_relative(
        obs=obs,
        stations=cart_stations,
        engine=engine,
        initial=initial_cart,
        damping=damping,
        max_iter=max_iter,
        tol_space=tol_space_m,
    )

    hypo_geo = Hypocenter(metric_to_wgs(result.hypocenter.point), Date())

    return hypo_geo, result, ref_geo


# ---------------------------------------------------------------------------
# Separate origin-time estimation from absolute picks (Cartesian)
# ---------------------------------------------------------------------------


def estimate_origin_time(
    picks: Sequence[PhasePick],
    stations: Dict[str, CartStation],
    engine: TravelTimeEngine,
    source: CartHypocenter,
) -> OriginTimeResult:
    """Estimate origin time t0 from absolute picks and a fixed hypocenter.

    Given an already determined source position (x, y, z) and a set of
    absolute arrival times t_i^obs, this function solves the simple
    weighted least-squares problem:

        t_i^obs = t0 + T_i(x) + ε_i

    where T_i(x) is the traveltime from source to station i.

    With weights w_i (typically 1 / sigma_i), the solution is:

        t0_hat = sum_i w_i^2 (t_i^obs - T_i) / sum_i w_i^2

    Parameters
    ----------
    picks : sequence of PhasePick
        Phase picks with observed arrival times and weights.
    stations : dict
        Mapping station_id -> CartStation.
    engine : TravelTimeEngine
        Travel-time engine providing traveltimes.
    source : CartHypocenter
        Hypocenter with fixed (x, y, z). The current t0 is ignored.

    Returns
    -------
    OriginTimeResult
        Estimated origin time and basic residual statistics.
    """
    if not picks:
        raise ValueError("At least one pick is required to estimate t0.")

    reference = picks[0].time
    if reference is None:
        raise ValueError("Origin-time estimation requires timed picks.")
    y = []
    w = []
    for pick in picks:
        if pick.time is None:
            raise ValueError("Origin-time estimation requires timed picks.")
        station = stations[pick.station]
        tt = engine.traveltime(pick.phase, source, station)
        if not np.isfinite(tt):
            raise ValueError(f"No finite traveltime for station {pick.station}.")
        y.append((pick.time - reference) - tt)
        w.append(pick.weight if pick.weight is not None else 1.0)

    y_arr = np.asarray(y, dtype=float)
    w_arr = np.asarray(w, dtype=float)

    w2 = w_arr * w_arr
    denom = float(np.sum(w2))
    if denom <= 0.0:
        raise ValueError("Non-positive sum of weights when estimating t0.")

    t0_hat = float(np.sum(w2 * y_arr) / denom)

    residuals = y_arr - t0_hat
    rms = float(np.sqrt(np.mean(residuals**2)))
    wr = w_arr * residuals
    wrms = float(np.sqrt(np.mean(wr**2)))

    return OriginTimeResult(
        t0=reference + t0_hat,
        rms=rms,
        weighted_rms=wrms,
        residuals=residuals,
    )
