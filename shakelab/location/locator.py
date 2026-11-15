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
# ShakeLab is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Geiger hypocenter location with modular travel-time engine.

This module implements a damped, weighted Geiger location algorithm,
with optional covariance estimation, in a way that is independent
from the specific velocity model. The user must provide a concrete
TravelTimeEngine implementation.

Design overview
---------------
The core inversion operates in a local Cartesian coordinate system
obtained from WGS84 geodetic coordinates via the projection tools
in :mod:`shakelab.libutils.geodeticN`. The travel-time engine selects
the frame: translated ECEF for homogeneous travel times and local
East-North-Up (ENU) for a horizontally layered model.

In this frame, Cartesian coordinates are represented by the
lightweight :class:`MetricPoint` container, which stores ``x, y, z``
in meters. The internal classes :class:`CartStation` and
:class:`CartHypocenter` wrap station and hypocenter locations in
the selected local frame and are the objects effectively
used by the Geiger algorithm.

User-facing code interacts instead with geographic coordinates
through the :class:`Station` and :class:`Hypocenter` classes, which
are based on :class:`shakelab.libutils.geodeticN.WgsPoint` (lat, lon,
elevation in meters) and on :class:`shakelab.libutils.timeN.Date`
(origin time).

The high-level function :func:`locate_event` performs all necessary
conversions:

    WgsPoint  →  MetricPoint
              →  Geiger inversion
              →  MetricPoint  →  WgsPoint

This isolates the numerical inversion from the details of geographic
coordinate handling and enforces a consistent representation
throughout the location subsystem.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence, Tuple, Optional, List

import numpy as np

from shakelab.libutils.timeN import Date

from shakelab.libutils.geodeticN import (  # type: ignore
    WgsPoint,
    MetricFrame,
    MetricPoint,
    wgs_to_metric,
    metric_to_wgs,
)
from shakelab.location.ttengine import (
    TravelTimeEngine,
    HomogeneousEngine,
    LayeredTravelTimeEngine,
)
from shakelab.location.phases import PhasePick  # type: ignore


# ---------------------------------------------------------------------------
# Geographic data structures (user-facing)
# ---------------------------------------------------------------------------


@dataclass
class Station:
    """
    Seismic station in geographic coordinates.

    Parameters
    ----------
    station_id : str
        Identifier of the station (network code, name, etc.).
    location : WgsPoint
        Station geographic coordinates in WGS84. The embedded
        latitude and longitude are expressed in degrees, while the
        elevation is expressed in meters and is positive upward
        (e.g. above sea level or ellipsoid). The same vertical
        datum must be used consistently for all stations and for
        the reference point used by the locator.

    Notes
    -----
    Older versions of this module stored the latitude, longitude
    and elevation as separate attributes ``lat``, ``lon`` and
    ``elev``. For convenience, read-only properties are still
    provided that forward to the internal :class:`WgsPoint`
    instance.
    """

    station_id: str
    location: WgsPoint

    @property
    def lat(self) -> float:
        """Latitude in degrees (read-only convenience alias)."""
        return self.location.latitude

    @property
    def lon(self) -> float:
        """Longitude in degrees (read-only convenience alias)."""
        return self.location.longitude

    @property
    def elev(self) -> float:
        """Station elevation in meters, positive upward."""
        return self.location.elevation


@dataclass
class Hypocenter:
    """
    Earthquake hypocenter in geographic coordinates.

    Parameters
    ----------
    location : WgsPoint
        Geographic hypocentral position in WGS84. As for stations,
        latitude and longitude are expressed in degrees, and
        ``location.elevation`` is expressed in meters and is
        positive upward. Hypocentral depth is therefore encoded as
        a **negative elevation**:

            depth_m = -location.elevation
    t0 : Date
        Absolute origin time represented as a :class:`Date`.

    Notes
    -----
    Earlier versions exposed latitude, longitude, depth (in meters,
    positive downward) and origin time (in seconds) as separate
    attributes. To ease the transition, read-only properties are
    provided that emulate this API:

    * ``lat`` and ``lon`` forward to ``location.latitude`` and
      ``location.longitude``;
    * ``depth`` returns the hypocentral depth in meters (positive
      downward), computed as ``-location.elevation``;
    * ``depth_km`` returns the depth in kilometers;
    * ``t0_seconds`` returns the origin time in seconds and can be
      used in external code that still expects a float.
    """

    location: WgsPoint
    t0: Date

    @property
    def lat(self) -> float:
        """Hypocenter latitude in degrees (read-only convenience)."""
        return self.location.latitude

    @property
    def lon(self) -> float:
        """Hypocenter longitude in degrees (read-only convenience)."""
        return self.location.longitude

    @property
    def depth(self) -> float:
        """Hypocentral depth in meters (positive downward)."""
        return -self.location.elevation

    @property
    def depth_km(self) -> float:
        """Hypocentral depth in kilometers (positive downward)."""
        return -self.location.elevation / 1000.0

    @property
    def t0_seconds(self) -> float:
        """
        Origin time in seconds, convenience alias for
        ``t0.seconds``.
        """
        return self.t0.seconds


# ---------------------------------------------------------------------------
# Basic data structures (Cartesian, internal)
# ---------------------------------------------------------------------------


@dataclass
class CartStation:
    """
    Seismic station expressed in local Cartesian coordinates.

    This internal representation stores the station position as a
    :class:`MetricPoint` in a local Cartesian frame defined by
    :func:`wgs_to_metric`. In typical ShakeLab workflows this is a
    translated ECEF frame, but it may also be a local ENU frame
    depending on how the conversion function is called.

    Units are meters.

    Attributes
    ----------
    station_id : str
        Station identifier.
    point : MetricPoint
        Cartesian coordinates in meters in the chosen local frame.
    depth : float
        Physical station depth in meters, positive downward. This is
        used by the layered engine. A value of zero places the station
        at the top of the velocity model.
    """

    station_id: str
    point: MetricPoint
    depth: float = 0.0


@dataclass
class CartHypocenter:
    """
    Earthquake hypocenter in Cartesian coordinates.

    Parameters
    ----------
    point : MetricPoint
        Cartesian source location (meters) in the same local frame
        used for the stations (see
        :class:`CartStation` and :func:`wgs_to_metric`).
    t0 : float
        Origin time in seconds relative to the inversion's time
        reference, or seconds from the Date epoch if no reference is set.

    Notes
    -----
    The Geiger inversion operates on a model vector
    ``[x, y, z, t0]`` constructed from this class. Only the
    Cartesian position and the origin time are needed at this
    stage. The physical interpretation in terms of latitude,
    longitude and depth is recovered by projecting the final
    solution back to geodetic coordinates via
    :func:`metric_to_wgs`.
    """

    point: MetricPoint
    t0: float


@dataclass
class LocationResult:
    """
    Result of a Geiger hypocenter location.

    Attributes
    ----------
    hypocenter : CartHypocenter
        Final hypocenter estimate in Cartesian coordinates.
    cov_matrix : ndarray or None
        Covariance matrix of the four model parameters
        (x, y, z, t0), if available.
    std_x, std_y, std_z, std_t0 : float or None
        Standard deviations of the model parameters, derived from
        the covariance matrix when available.
    rms : float
        Root-mean-square residual (unweighted).
    weighted_rms : float
        Root-mean-square residual with data weights applied.
    num_iterations : int
        Number of Geiger iterations performed.
    num_picks : int
        Number of picks effectively used in the inversion.
    converged : bool
        True if both spatial and temporal convergence criteria
        were satisfied; False otherwise.
    residuals : ndarray
        Final residual vector used to compute RMS statistics.
    time_reference : Date or None
        Reference for ``hypocenter.t0``. Add the offset to this date to
        recover the origin time. None denotes the Date epoch.
    """

    hypocenter: CartHypocenter
    cov_matrix: Optional[np.ndarray]
    std_x: Optional[float]
    std_y: Optional[float]
    std_z: Optional[float]
    std_t0: Optional[float]
    rms: float
    weighted_rms: float
    num_iterations: int
    num_picks: int
    converged: bool
    residuals: np.ndarray
    time_reference: Optional[Date] = None


# ---------------------------------------------------------------------------
# Core linear algebra utilities
# ---------------------------------------------------------------------------


def _arrival_seconds(time: Date, reference: Optional[Date]) -> float:
    return time.seconds if reference is None else time - reference


def build_design_matrices(
    picks: Sequence[PhasePick],
    stations: Dict[str, CartStation],
    engine: TravelTimeEngine,
    source: CartHypocenter,
    *,
    time_reference: Optional[Date] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build design matrix G, residual vector r and weights w.

    Parameters
    ----------
    picks : sequence of PhasePick
        Phase picks for which traveltimes are computed. These are
        instances of :class:`shakelab.location.phases.PhasePick`.
        All picks must carry a non-None ``time`` attribute; this
        condition is enforced by
        :func:`prepare_picks_for_location`.
    stations : dict
        Mapping ``station_id -> CartStation``. The station identifier
        must match the ``station`` field of each PhasePick.
    engine : TravelTimeEngine
        Travel-time engine providing traveltimes and gradients.
    source : CartHypocenter
        Current hypocenter estimate.
    time_reference : Date, optional
        Reference for observation times and ``source.t0``.
        Defaults to the Date epoch.

    Returns
    -------
    G : ndarray, shape (N, 4)
        Design matrix with partial derivatives of traveltime with
        respect to x, y, z, t0.
    r : ndarray, shape (N,)
        Data residuals: observed_time - predicted_time.
    w : ndarray, shape (N,)
        Data weights.

    Notes
    -----
    If ``pick.weight`` is None, a default
    weight of 1.0 is used. Picks with an unknown station identifier
    should have been removed beforehand by
    :func:`prepare_picks_for_location`; if any remain, a KeyError
    is raised.
    """
    n = len(picks)
    G = np.zeros((n, 4), dtype=float)
    r = np.zeros(n, dtype=float)
    w = np.zeros(n, dtype=float)

    for i, pick in enumerate(picks):
        if pick.time is None:
            raise ValueError(
                "PhasePick without time encountered in "
                "build_design_matrices. Use "
                "prepare_picks_for_location(require_time=True) "
                "before calling the locator."
            )

        station_id = pick.station
        if station_id not in stations:
            raise KeyError(
                f"Station {station_id!r} in picks is not present in the "
                "stations dictionary."
            )

        station = stations[station_id]

        t_pred = engine.traveltime(pick.phase, source, station)
        dt_dx, dt_dy, dt_dz = engine.gradients(pick.phase, source, station)

        G[i, 0] = dt_dx
        G[i, 1] = dt_dy
        G[i, 2] = dt_dz
        G[i, 3] = 1.0

        t_obs = _arrival_seconds(pick.time, time_reference)
        r[i] = t_obs - (source.t0 + t_pred)

        weight = pick.weight if pick.weight is not None else 1.0
        w[i] = float(weight)

    return G, r, w


def build_normal_equations(
    G: np.ndarray,
    r: np.ndarray,
    w: np.ndarray,
    damping: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build normal equations for the linearized Geiger step.

    Parameters
    ----------
    G : ndarray, shape (N, 4)
        Design matrix.
    r : ndarray, shape (N,)
        Data residuals.
    w : ndarray, shape (N,)
        Data weights.
    damping : float
        Damping parameter (Levenberg-Marquardt style).

    Returns
    -------
    A : ndarray, shape (4, 4)
        Normal matrix including damping:
        A = G^T W^2 G + lambda^2 I.
    b : ndarray, shape (4,)
        Right-hand side: b = G^T W^2 r.
    A_data : ndarray, shape (4, 4)
        Data-only normal matrix: A_data = G^T W^2 G.
        This is used for covariance estimation.
    """
    weighted_G = G * w[:, None]
    a_data = weighted_G.T @ weighted_G
    a = a_data + (damping ** 2) * np.eye(G.shape[1])
    b = weighted_G.T @ (w * r)
    return a, b, a_data


def compute_covariance(
    G: np.ndarray,
    r: np.ndarray,
    w: np.ndarray,
    A_data: np.ndarray,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Estimate covariance matrix and parameter standard deviations.

    Parameters
    ----------
    G : ndarray, shape (N, 4)
        Design matrix.
    r : ndarray, shape (N,)
        Data residuals.
    w : ndarray, shape (N,)
        Data weights.
    A_data : ndarray, shape (4, 4)
        Data-only normal matrix (without damping).

    Returns
    -------
    cov : ndarray or None
        Covariance matrix of the model parameters, or None if the
        matrix cannot be inverted.
    std : ndarray or None
        Standard deviations (square root of the diagonal of cov),
        or None if cov is None.
    """
    dof = max(G.shape[0] - G.shape[1], 1)
    weighted_residuals = w * r
    residual_variance = (weighted_residuals @ weighted_residuals) / dof

    try:
        cov = np.linalg.inv(A_data)
    except np.linalg.LinAlgError:
        return None, None

    cov *= residual_variance
    std = np.sqrt(np.diag(cov))

    return cov, std


# ---------------------------------------------------------------------------
# Helper for preparing picks
# ---------------------------------------------------------------------------


def prepare_picks_for_location(
    picks: Sequence[PhasePick],
    stations: Dict[str, Station],
    require_time: bool = True,
    sort_by_time: bool = True,
) -> List[PhasePick]:
    """
    Prepare and validate phase picks before running the locator.

    This helper centralizes the basic sanity checks on the input
    picks and returns a cleaned list that is safe to pass to the
    Geiger inversion. It performs:

    * optional filtering of picks without a valid time stamp;
    * filtering of picks whose station code is not present in the
      ``stations`` dictionary;
    * optional sorting of picks by observed time.

    Parameters
    ----------
    picks : sequence of PhasePick
        Raw phase picks as produced by the reading/parsing stage.
        These are instances of
        :class:`shakelab.location.phases.PhasePick`.
    stations : dict
        Mapping ``station_id -> Station``. The station identifier
        must match the ``station`` field of each PhasePick.
    require_time : bool, optional
        If ``True`` (default), picks with ``time is None`` are
        discarded. If ``False``, picks without time are kept, but
        note that the current implementation of the locator
        requires a valid time for traveltime residuals, so such
        picks will not be usable.
    sort_by_time : bool, optional
        If ``True`` (default), the returned list is sorted by
        arrival time. Picks without time
        (if retained via ``require_time=False``) are placed at the
        beginning of the list.

    Returns
    -------
    list of PhasePick
        Cleaned and optionally sorted list of picks that are
        consistent with the provided station dictionary.

    Raises
    ------
    ValueError
        If no valid picks remain after filtering.

    Notes
    -----
    The function does not modify the input picks in place. All
    filtering and sorting are applied to a separate list. The
    data weights ``pick.weight`` are not altered here; missing
    weights are handled downstream (by assuming a default of 1.0
    when building the design matrix).
    """
    cleaned: List[PhasePick] = []

    for pick in picks:
        if pick.station not in stations:
            continue
        if require_time and pick.time is None:
            continue
        cleaned.append(pick)

    if not cleaned:
        raise ValueError(
            "No valid picks available for location after filtering. "
            "Check that station codes are consistent and that picks "
            "carry a valid time stamp."
        )

    if sort_by_time:
        timed = [pick for pick in cleaned if pick.time is not None]
        timed.sort(key=lambda pick: pick.time)
        cleaned = [pick for pick in cleaned if pick.time is None] + timed

    return cleaned


# ---------------------------------------------------------------------------
# Geiger location
# ---------------------------------------------------------------------------


def geiger_locate(
    picks: Sequence[PhasePick],
    stations: Dict[str, CartStation],
    engine: TravelTimeEngine,
    initial: CartHypocenter,
    damping: float = 0.1,
    max_iter: int = 20,
    tol_space: float = 1.0,
    tol_time: float = 0.01,
    space_scale: float = 1000.0,
    *,
    time_reference: Optional[Date] = None,
) -> LocationResult:
    """
    Iterative Geiger hypocenter location with covariance estimate.

    Parameters
    ----------
    picks : sequence of PhasePick
        Phase picks with observed arrival times and weights. These
        should normally be the output of
        :func:`prepare_picks_for_location`.
    stations : dict
        Mapping station_id -> CartStation.
    engine : TravelTimeEngine
        Travel-time engine providing traveltimes and gradients.
    initial : CartHypocenter
        Initial hypocenter guess in Cartesian coordinates.
    damping : float, optional
        Damping factor (Levenberg-Marquardt style).
    max_iter : int, optional
        Maximum number of iterations.
    tol_space : float, optional
        Spatial convergence tolerance in the same units as station
        coordinates.
    tol_time : float, optional
        Temporal convergence tolerance in seconds.
    space_scale : float, optional
        Scaling factor applied to the spatial parameters in the
        linear system to improve conditioning. Conceptually, the
        inversion is performed in scaled variables

            x' = x / space_scale,  y' = y / space_scale,
            z' = z / space_scale,  t0' = t0,

        and the updates are mapped back to meters at the end of
        each iteration.
    time_reference : Date, optional
        Reference for observation times and ``initial.t0``.
        Defaults to the Date epoch.

    Returns
    -------
    LocationResult
        Result structure containing the final hypocenter estimate,
        covariance (if available), residuals and iteration statistics.
    """
    # Model vector in physical units (meters, seconds)
    m = np.array(
        [initial.point.x, initial.point.y, initial.point.z, initial.t0],
        dtype=float,
    )

    converged = False
    iter_count = 0

    for _ in range(max_iter):
        iter_count += 1

        source_point = MetricPoint(
            x=float(m[0]),
            y=float(m[1]),
            z=float(m[2]),
            frame=initial.point.frame,
        )
        source = CartHypocenter(point=source_point, t0=float(m[3]))

        # Design matrix and residuals in physical units
        G, r, w = build_design_matrices(
            picks, stations, engine, source, time_reference=time_reference,
        )
        if not (np.all(np.isfinite(G)) and np.all(np.isfinite(r))):
            break

        # Build scaled design matrix for the solver:
        # columns for x,y,z are multiplied so that the corresponding
        # parameters are effectively expressed in kilometers-like
        # units, improving numerical conditioning.
        G_solver = G.copy()
        if space_scale != 1.0:
            G_solver[:, 0:3] *= space_scale

        A, b, _ = build_normal_equations(G_solver, r, w, damping)

        try:
            dm_scaled = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            # Singular or ill-conditioned system: abort iterations
            break

        # Map updates back to physical units:
        # dm_scaled refers to the scaled parameters m' = [x', y', z', t0']
        dm = dm_scaled.copy()
        if space_scale != 1.0:
            dm[0:3] *= space_scale

        if isinstance(engine, LayeredTravelTimeEngine):
            # Backtrack when a full step leaves the model or increases
            # the weighted residual norm.
            cost = np.linalg.norm(w * r)
            for _ in range(12):
                trial = m + dm
                point = MetricPoint(*trial[:3], frame=initial.point.frame)
                candidate = CartHypocenter(point, float(trial[3]))
                if engine.contains_source(candidate):
                    trial_residuals = np.empty_like(r)
                    for i, pick in enumerate(picks):
                        time = engine.traveltime(
                            pick.phase, candidate, stations[pick.station],
                        )
                        observed = _arrival_seconds(pick.time, time_reference)
                        trial_residuals[i] = observed - (trial[3] + time)
                    if np.linalg.norm(w * trial_residuals) <= cost:
                        break
                dm *= 0.5
            else:
                break

        m_new = m + dm

        dx = float(dm[0])
        dy = float(dm[1])
        dz = float(dm[2])
        dt = float(dm[3])

        step_space = float(np.sqrt(dx * dx + dy * dy + dz * dz))
        step_time = abs(dt)

        m = m_new
        if step_space <= tol_space and step_time <= tol_time:
            converged = True
            break

    final_point = MetricPoint(
        x=float(m[0]),
        y=float(m[1]),
        z=float(m[2]),
        frame=initial.point.frame,
    )
    hypo_cart = CartHypocenter(point=final_point, t0=float(m[3]))
    G, residuals, w = build_design_matrices(
        picks, stations, engine, hypo_cart, time_reference=time_reference,
    )

    weighted_G = G * w[:, None]
    cov, std = compute_covariance(G, residuals, w, weighted_G.T @ weighted_G)
    if std is not None:
        std_x, std_y, std_z, std_t0 = [float(value) for value in std]
    else:
        std_x = std_y = std_z = std_t0 = None

    rms = float(np.sqrt(np.mean(residuals ** 2))) if residuals.size else 0.0
    wrms = (
        float(np.sqrt(np.mean((w * residuals) ** 2)))
        if residuals.size else 0.0
    )

    return LocationResult(
        hypocenter=hypo_cart,
        cov_matrix=cov,
        std_x=std_x,
        std_y=std_y,
        std_z=std_z,
        std_t0=std_t0,
        rms=rms,
        weighted_rms=wrms,
        num_iterations=iter_count,
        num_picks=len(picks),
        converged=converged,
        residuals=residuals,
        time_reference=time_reference,
    )


# ---------------------------------------------------------------------------
# Geographic wrapper utilities
# ---------------------------------------------------------------------------


def _auto_ref_geo(
    stations: Dict[str, Station],
    initial: Hypocenter,
) -> WgsPoint:
    """
    Compute a reasonable reference geodetic point for the locator.

    The reference point defines the origin of the local Cartesian
    system used internally by the Geiger algorithm. In this helper,
    it is chosen as the (unweighted) average of station coordinates:

        lat0 = mean(lat_i)
        lon0 = mean(lon_i)
        h0   = mean(h_i)

    where lat_i, lon_i and h_i are the latitude, longitude and
    elevation of each station's :class:`WgsPoint` location. The
    initial hypocenter is currently not included in the average, but
    it is kept as an explicit argument for potential future
    refinements (e.g. weighting the reference toward the expected
    source region, or snapping to the initial epicenter).

    Parameters
    ----------
    stations : dict
        Mapping ``station_id -> Station`` in geographic coordinates.
        Each station contains a :class:`WgsPoint` with latitude,
        longitude and elevation in meters (positive upward). At
        least one station must be provided.
    initial : Hypocenter
        Initial hypocenter guess in geographic coordinates. The
        current implementation does not use it explicitly, but it
        is passed through to preserve API flexibility and to
        document that the reference point conceptually belongs to
        the source–station geometry.

    Returns
    -------
    WgsPoint
        Reference geodetic point in WGS84 (latitude [deg], longitude
        [deg], elevation [m, positive upward]) to be used as origin
        for the local Cartesian system (via :func:`wgs_to_metric`).

    Raises
    ------
    ValueError
        If no stations are provided.

    Notes
    -----
    The precise choice of the reference point does not affect the
    physical meaning of the solution, as long as the transformation
    from geographic to Cartesian coordinates is consistent for all
    stations and for the hypocenter. However, choosing a reference
    that lies close to the cluster of stations helps reduce the
    magnitude of Cartesian coordinates and can improve numerical
    conditioning in some cases.
    """
    if not stations:
        raise ValueError("At least one station is required to infer ref_geo.")

    lats = [st.location.latitude for st in stations.values()]
    lons = [st.location.longitude for st in stations.values()]
    elevs = [st.location.elevation for st in stations.values()]

    lat0 = float(sum(lats) / len(lats))
    lon0 = float(sum(lons) / len(lons))
    h0 = float(sum(elevs) / len(elevs))

    return WgsPoint(
        longitude=lon0,
        latitude=lat0,
        elevation=h0,
    )


def locate_event(
    picks: Sequence[PhasePick],
    stations: Dict[str, Station],
    engine: TravelTimeEngine,
    initial: Hypocenter,
    ref_geo: Optional[WgsPoint] = None,
    damping: float = 0.1,
    max_iter: int = 20,
    tol_space_m: float = 1.0,
    tol_time_s: float = 0.01,
) -> Tuple[Hypocenter, LocationResult, WgsPoint]:
    """
    Locate an earthquake from picks in geographic coordinates.

    This function is a user-friendly wrapper around the Cartesian
    Geiger locator. The user provides stations and an initial
    hypocenter in geographic coordinates; internally, coordinates
    are converted to a local Cartesian system using
    :func:`wgs_to_metric`, the Geiger inversion is performed in
    that system, and the final solution is converted back to
    geographic coordinates with :func:`metric_to_wgs`.

    The function is deliberately kept stateless. All information
    required for the inversion is passed explicitly: picks, station
    geometry, travel-time engine, and initial hypocenter.
    Pick times and origin time are absolute dates.

    Parameters
    ----------
    picks : sequence of PhasePick
        Phase picks with observed arrival times and weights. The
        low-level locator only needs the station identifier,
        observed time and weight; any additional metadata are
        ignored. The input picks are first passed through
        :func:`prepare_picks_for_location`.
    stations : dict
        Mapping ``station_id -> Station`` in geographic coordinates.
        Each station carries a :class:`WgsPoint` with latitude,
        longitude and elevation in meters (positive upward).
    engine : TravelTimeEngine
        Travel-time engine compatible with the Cartesian locator.
    initial : Hypocenter
        Initial hypocenter guess in geographic coordinates. The
        origin time is the reference for the relative times used
        in the inversion.
    ref_geo : WgsPoint, optional
        Reference geodetic point (latitude [deg], longitude [deg],
        elevation [m, positive upward]) used to define the local
        Cartesian system. If ``None``, it is automatically inferred
        from the station coordinates by :func:`_auto_ref_geo`.
    damping : float, optional
        Damping factor (Levenberg-Marquardt style) for the Geiger
        iterations.
    max_iter : int, optional
        Maximum number of iterations.
    tol_space_m : float, optional
        Spatial convergence tolerance in meters.
    tol_time_s : float, optional
        Temporal convergence tolerance in seconds.

    Returns
    -------
    hypo_geo : Hypocenter
        Estimated hypocenter in geographic coordinates. The
        location is returned as a :class:`WgsPoint` (with elevation
        in meters), and the origin time as a :class:`Date`.
    result : LocationResult
        Full result of the Cartesian Geiger inversion, including
        covariance (if available) and residuals. Its Cartesian
        hypocenter time is relative to ``result.time_reference``.
    ref_geo_out : WgsPoint
        Reference geodetic point actually used for the conversion.
        This is useful to reproduce the Cartesian coordinates or
        to interpret the covariance matrix externally.

    Notes
    -----
    Hypocentral depth is not stored explicitly in the returned
    :class:`Hypocenter`. Instead, the vertical coordinate is the
    elevation in meters (positive upward), and depth can be
    obtained as ``depth_m = -hypo_geo.location.elevation``.
    Internally, the projection functions operate on
    :class:`WgsPoint` and :class:`MetricPoint` via
    :func:`wgs_to_metric` and :func:`metric_to_wgs`; the
    :class:`MetricPoint` carries the engine's frame. For ENU, the
    layered model starts at the plane ``up = -ref_geo.elevation``.
    """
    if not stations:
        raise ValueError("At least one station is required for location.")

    # Clean and validate picks before running the locator
    picks_prepared = prepare_picks_for_location(picks, stations)

    # Determine reference geographic point if not provided
    if ref_geo is None:
        ref_geo = _auto_ref_geo(stations, initial)

    orientation = getattr(engine, "metric_frame", "ecef")
    frame = MetricFrame(ref_geo=ref_geo, orientation=orientation)

    # Convert stations from geographic to Cartesian coordinates
    cart_stations: Dict[str, CartStation] = {}
    for st_id, st in stations.items():
        mp = wgs_to_metric(st.location, frame=frame)

        elevation = st.location.elevation
        if frame.orientation == "enu":
            elevation = mp.z + frame.ref_geo[2]
        depth = max(0.0, -elevation)

        cart_stations[st_id] = CartStation(
            station_id=st_id,
            point=mp,
            depth=depth,
        )

    # Convert initial hypocenter to Cartesian coordinates
    mp0 = wgs_to_metric(initial.location, frame=frame)
    initial_cart = CartHypocenter(
        point=mp0,
        t0=0.0,
    )

    # Run the Cartesian Geiger locator
    result = geiger_locate(
        picks=picks_prepared,
        stations=cart_stations,
        engine=engine,
        initial=initial_cart,
        damping=damping,
        max_iter=max_iter,
        tol_space=tol_space_m,
        tol_time=tol_time_s,
        time_reference=initial.t0,
    )

    hypo_geo = Hypocenter(
        location=metric_to_wgs(result.hypocenter.point),
        t0=initial.t0 + result.hypocenter.t0,
    )

    return hypo_geo, result, ref_geo


# ---------------------------------------------------------------------------
# Initial hypocentre model guess
# ---------------------------------------------------------------------------

def locate_from_trilateration(stations, phases, depth=10000):
    """
    """
    pass
    # return hypo_init_guess

def locate_from_grid_search():
    """
    """
    pass
    # return hypo_init_guess
