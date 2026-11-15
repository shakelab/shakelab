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
Travel-time engines for hypocenter location.

This module defines the abstract :class:`TravelTimeEngine` interface and a
few simple concrete implementations that operate in a local Cartesian
coordinate system:

* :class:`HomogeneousEngine` implements straight-ray propagation in a
  homogeneous medium with constant P- and S-wave velocities, mainly
  intended for testing.
* :class:`LayeredTravelTimeEngine` wraps the 1D layered traveltime
  solver in :mod:`shakelab.location.traveltime` to provide a simple
  engine for horizontally layered models.

These engines are used by the Geiger locator in
:mod:`shakelab.location.locator`, but can in principle be re-used by
other modules requiring traveltimes and their spatial gradients.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, Optional, Tuple, Union

import math
import numpy as np

from shakelab.location.ttsolver import (  # type: ignore
    LayeredModel1D,
    solve_traveltime,
    TOLERANCE,
    MAX_ITER,
)


class TravelTimeEngine:
    """
    Abstract interface for traveltime and gradient computation.

    Concrete implementations must provide traveltimes and spatial
    gradients for arbitrary Cartesian source and receiver
    coordinates. The Cartesian coordinates must be expressed in
    the same local frame for both source and stations.
    """

    metric_frame = "ecef"

    def traveltime(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> float:
        """
        Return traveltime [s] for the specified phase.

        Parameters
        ----------
        phase : str
            Phase identifier (e.g. ``"P"`` or ``"S"``).
        source : CartHypocenter
            Source location in the local Cartesian system.
        station : CartStation
            Station location in the same Cartesian system.

        Returns
        -------
        float
            Traveltime in seconds.
        """
        raise NotImplementedError

    def gradients(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Tuple[float, float, float]:
        """
        Return spatial traveltime gradients (dT/dx, dT/dy, dT/dz).

        The derivatives are taken with respect to the Cartesian
        coordinates of the source.

        Parameters
        ----------
        phase : str
            Phase identifier (e.g. ``"P"`` or ``"S"``).
        source : CartHypocenter
            Source location in the local Cartesian system.
        station : CartStation
            Station location in the same Cartesian system.

        Returns
        -------
        (float, float, float)
            Traveltime derivatives (dT/dx, dT/dy, dT/dz), in seconds
            per meter.
        """
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Simple homogeneous travel-time engine
# ---------------------------------------------------------------------------


class HomogeneousEngine(TravelTimeEngine):
    """
    Very simple homogeneous medium travel-time engine.

    This implementation is intended only for testing the inversion
    machinery. It assumes constant P- and S-wave velocities and
    straight rays in a homogeneous medium.
    """

    def __init__(self, vp: float, vs: float):
        """
        Parameters
        ----------
        vp : float
            P-wave velocity (m/s).
        vs : float
            S-wave velocity (m/s).
        """
        if vp <= 0.0 or vs <= 0.0:
            raise ValueError("vp and vs must be positive.")

        self.vp = float(vp)
        self.vs = float(vs)

    def _phase_velocity(self, phase: str) -> float:
        """
        Return phase velocity for the requested phase.

        Parameters
        ----------
        phase : str
            Phase identifier (``"P"`` or ``"S"``).

        Returns
        -------
        float
            Phase velocity in m/s.

        Raises
        ------
        ValueError
            If an unsupported phase is requested.
        """
        ph = phase.upper()
        if ph == "P":
            return self.vp
        if ph == "S":
            return self.vs

        raise ValueError(f"Unsupported phase {phase!r} for HomogeneousEngine.")

    def traveltime(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> float:
        """
        Traveltime in a homogeneous medium.

        For a homogeneous medium with straight rays, the traveltime is
        simply:

            T = R / v,

        where R is the Euclidean distance between source and station and
        v is the phase velocity.

        Parameters
        ----------
        phase : str
            Phase identifier (``"P"`` or ``"S"``).
        source : CartHypocenter
            Source location in the local Cartesian system.
        station : CartStation
            Station location in the same Cartesian system.

        Returns
        -------
        float
            Traveltime in seconds.
        """
        v = self._phase_velocity(phase)

        dx = station.point.x - source.point.x
        dy = station.point.y - source.point.y
        dz = station.point.z - source.point.z

        dist = float(np.sqrt(dx * dx + dy * dy + dz * dz))
        return dist / v

    def gradients(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Tuple[float, float, float]:
        """
        Spatial gradients of traveltime in a homogeneous medium.

        For T = R / v, with

            R = ||s - x|| = sqrt(dx^2 + dy^2 + dz^2),

        where ``dx = sx - x``, ``dy = sy - y`` and ``dz = sz - z``,
        the gradients become

            dT/dx = -(dx) / (v * R),
            dT/dy = -(dy) / (v * R),
            dT/dz = -(dz) / (v * R).

        Parameters
        ----------
        phase : str
            Phase identifier (``"P"`` or ``"S"``).
        source : CartHypocenter
            Source location in the local Cartesian system.
        station : CartStation
            Station location in the same Cartesian system.

        Returns
        -------
        (float, float, float)
            Traveltime derivatives (dT/dx, dT/dy, dT/dz), in seconds
            per meter.
        """
        v = self._phase_velocity(phase)

        dx = station.point.x - source.point.x
        dy = station.point.y - source.point.y
        dz = station.point.z - source.point.z

        dist = float(np.sqrt(dx * dx + dy * dy + dz * dz))
        if dist == 0.0:
            # Traveltime is zero and gradients are undefined; return 0
            # to avoid division by zero. This situation should not
            # occur in practical location problems.
            return 0.0, 0.0, 0.0

        coeff = 1.0 / (v * dist)

        dtx = -coeff * dx
        dty = -coeff * dy
        dtz = -coeff * dz

        return dtx, dty, dtz


# ---------------------------------------------------------------------------
# Layered traveltime engine
# ---------------------------------------------------------------------------


@dataclass
class LayeredTravelTimeEngine(TravelTimeEngine):
    """
    Travel-time engine for horizontally layered media.

    This engine wraps the 1D layered solver provided by
    :mod:`shakelab.location.ttsolver` and exposes traveltimes and
    spatial gradients in the local Cartesian system used by the
    locator. First arrivals use direct branches and head waves by default.

    Parameters
    ----------
    model : LayeredModel1D
        1D velocity model with layer thickness and Vp/Vs profiles.
        Depth is measured in meters, positive downward from the free
        surface.
    arrivals : str, optional
        Arrival mode passed to :func:`solve_traveltime`. Accepted
        values are:

        * ``"direct"`` – only direct branches are considered;
        * ``"all"`` – direct branches and head waves are considered,
          with legacy interface reflections available by explicit opt-in.

        The default is ``"all"``.
    include_diving : bool, optional
        Include legacy interface reflections when ``arrivals="all"``.
        The default is ``False``.
    tol : float, optional
        Tolerance for the internal root-finding iterations of the
        solver, in meters. The default is :data:`TOLERANCE`.
    max_iter : int, optional
        Maximum number of iterations per ray. The default is
        :data:`MAX_ITER`.
    max_heads : int or None, optional
        Placeholder for future use. At present this value is ignored
        and kept only for API stability.
    require_convergence : bool, optional
        If ``True``, only branches that satisfy the convergence and
        residual criteria of the solver are kept. The default is
        ``True`` and should normally not be changed.
    same_depth_switch : float, optional
        Depth threshold (in meters) that was used in the original
        implementation to activate special same-depth logic. It is
        currently kept only for API compatibility and is not passed to
        :func:`solve_traveltime`.
    fd_step : float, optional
        Step size (in meters) used for finite-difference evaluation
        of traveltime gradients. The default is ``50.0``.
    gradient_method : {"analytical", "numerical"}, optional
        Method used by :meth:`gradients`. If ``"analytical"``,
        gradients are derived from the ray parameter and local
        velocity at the source (eikonal relation). If the analytical
        computation fails, the engine silently falls back to the
        numerical scheme. If ``"numerical"``, derivatives are always
        evaluated by finite differences.

    Notes
    -----
    * The engine works in the local Cartesian system used by the
      locator (see :class:`CartStation` and :class:`CartHypocenter`).
      Traveltimes are computed from the *horizontal* separation
      between source and receiver and their depths.
    * Geographic locations use a common ENU frame, with model depth
      zero at the plane ``up = -frame.ref_geo[2]``. Legacy Cartesian
      inputs without a frame use ``z = 0`` as the top of the model.
      A station can specify its model depth explicitly.
    """

    model: LayeredModel1D
    arrivals: str = "all"
    include_diving: bool = False
    tol: float = TOLERANCE
    max_iter: int = MAX_ITER
    max_heads: Optional[int] = None
    require_convergence: bool = True
    same_depth_switch: float = 1.0
    fd_step: float = 50.0
    gradient_method: str = "analytical"
    # Preferred metric frame for this engine (used by higher-level
    # code such as the locator, if it chooses to honour it).
    metric_frame: str = "enu"

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _height_above_model_top(
        obj: Union[CartHypocenter, CartStation],
    ) -> float:
        frame = obj.point.frame
        reference_height = frame.ref_geo[2] if frame is not None else 0.0
        return float(obj.point.z) + reference_height

    def _get_depth(self, obj: Union[CartHypocenter, CartStation]) -> float:
        """Return model depth in meters, positive downward."""
        depth = getattr(obj, "depth", None)
        if depth is None:
            depth = -self._height_above_model_top(obj)
        return max(0.0, float(depth))

    def contains_source(self, source: CartHypocenter) -> bool:
        """Return whether the source is inside the layered velocity model."""
        depth = -self._height_above_model_top(source)
        return 0.0 <= depth < self.model.z_bottom[-1]

    @staticmethod
    def _horizontal_offset(
        source: CartHypocenter,
        station: CartStation,
    ) -> float:
        frame = source.point.frame
        if frame != station.point.frame or (
            frame is not None and frame.orientation != "enu"
        ):
            raise ValueError(
                "Layered travel times require a common ENU frame."
            )
        return math.hypot(
            station.point.x - source.point.x,
            station.point.y - source.point.y,
        )

    def _single_traveltime(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> float:
        """
        Compute traveltime for a single phase and station.

        This is a thin wrapper around :func:`solve_traveltime`
        operating on one source–receiver pair.

        Parameters
        ----------
        phase : str
            Phase identifier (e.g. ``"P"`` or ``"S"``).
        source : CartHypocenter
            Source location in the local Cartesian system.
        station : CartStation
            Station location in the same Cartesian system.

        Returns
        -------
        float
            Traveltime in seconds. Returns ``NaN`` if the solver
            fails or if the requested phase cannot be computed for
            the given geometry/model.
        """
        phase_u = phase.upper()
        if phase_u not in ("P", "S"):
            return math.nan

        offset = self._horizontal_offset(source, station)

        z_s = self._get_depth(source)
        z_r = self._get_depth(station)

        try:
            tt_arr, _ = solve_traveltime(
                model=self.model,
                offsets=[offset],
                source_depth=z_s,
                receiver_depth=z_r,
                phase=phase_u,
                arrivals=self.arrivals,
                return_branches=False,
                tol=self.tol,
                max_iter=self.max_iter,
                require_convergence=self.require_convergence,
                include_diving=self.include_diving,
            )
        except Exception:
            return math.nan

        arr = np.asarray(tt_arr, dtype=float)
        if arr.size == 0:
            return math.nan

        t_val = float(arr.flat[0])
        if not math.isfinite(t_val):
            return math.nan

        return t_val

    def _single_best_branch(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Optional[Dict]:
        """
        Compute the fastest branch for a single source–receiver pair.

        This helper wraps :func:`solve_traveltime` with
        ``return_branches=True`` and returns the branch dictionary
        corresponding to the minimum traveltime, or ``None`` if no
        valid branch can be found.
        """
        phase_u = phase.upper()
        if phase_u not in ("P", "S"):
            return None

        offset = self._horizontal_offset(source, station)

        z_s = self._get_depth(source)
        z_r = self._get_depth(station)

        try:
            tt_arr, branches_all = solve_traveltime(
                model=self.model,
                offsets=[offset],
                source_depth=z_s,
                receiver_depth=z_r,
                phase=phase_u,
                arrivals=self.arrivals,
                return_branches=True,
                tol=self.tol,
                max_iter=self.max_iter,
                require_convergence=self.require_convergence,
                include_diving=self.include_diving,
            )
        except Exception:
            return None

        if not branches_all or not branches_all[0]:
            return None

        branches = branches_all[0]
        try:
            best = min(
                branches,
                key=lambda b: b.get("t", math.inf),
            )
        except Exception:
            return None

        t_best = float(best.get("t", math.nan))
        if not math.isfinite(t_best):
            return None

        return best

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def traveltime(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> float:
        """
        Compute the theoretical traveltime between source and station.

        Parameters
        ----------
        phase : str
            Seismic phase identifier (case-insensitive). Only ``"P"``
            and ``"S"`` are currently supported.
        source : CartHypocenter
            Hypocenter in local Cartesian coordinates.
        station : CartStation
            Station in local Cartesian coordinates.

        Returns
        -------
        float
            Traveltime in seconds, or ``NaN`` if the solver fails
            for the given geometry.
        """
        return self._single_traveltime(phase, source, station)

    # ------------------------------------------------------------------
    # Gradients
    # ------------------------------------------------------------------

    def gradients_numerical(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Tuple[float, float, float]:
        """
        Numerical finite-difference approximation of traveltime
        gradients.

        Derivatives are evaluated by perturbing the source position
        in the local Cartesian system and calling
        :meth:`_single_traveltime`. The station is kept fixed.

        Parameters
        ----------
        phase : str
            Seismic phase identifier (``"P"`` or ``"S"``).
        source : CartHypocenter
            Hypocenter in local Cartesian coordinates.
        station : CartStation
            Station in local Cartesian coordinates.

        Returns
        -------
        dtx, dty, dtz : float
            Partial derivatives of traveltime with respect to the
            source coordinates (x, y, z), in s/m. If the base
            traveltime or any of the perturbed traveltimes cannot
            be computed, the corresponding component is set to
            ``NaN``.
        """
        # Base traveltime
        t0 = self._single_traveltime(phase, source, station)
        if not math.isfinite(t0):
            return math.nan, math.nan, math.nan

        h = float(self.fd_step)
        if h <= 0.0:
            raise ValueError("fd_step must be positive.")

        def _perturb(dx: float, dy: float, dz: float) -> float:
            point = replace(
                source.point,
                x=source.point.x + dx,
                y=source.point.y + dy,
                z=source.point.z + dz,
            )
            return self._single_traveltime(
                phase, replace(source, point=point), station,
            )

        # dT/dx
        t_plus = _perturb(h, 0.0, 0.0)
        t_minus = _perturb(-h, 0.0, 0.0)
        if math.isfinite(t_plus) and math.isfinite(t_minus):
            dtx = (t_plus - t_minus) / (2.0 * h)
        else:
            dtx = math.nan

        # dT/dy
        t_plus = _perturb(0.0, h, 0.0)
        t_minus = _perturb(0.0, -h, 0.0)
        if math.isfinite(t_plus) and math.isfinite(t_minus):
            dty = (t_plus - t_minus) / (2.0 * h)
        else:
            dty = math.nan

        # dT/dz
        t_plus = _perturb(0.0, 0.0, h)
        t_minus = _perturb(0.0, 0.0, -h)
        if math.isfinite(t_plus) and math.isfinite(t_minus):
            dtz = (t_plus - t_minus) / (2.0 * h)
        else:
            dtz = math.nan

        return dtx, dty, dtz

    def gradients_analytical(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Tuple[float, float, float]:
        """
        Analytical traveltime gradients based on the ray parameter.

        Project the slowness components of the first ray segment onto
        the local Cartesian axes. At an interface, use the one-sided
        derivative from the layer occupied by that segment.

        If the branch or ray parameter cannot be determined, the method
        returns ``NaN`` for all components.

        Parameters
        ----------
        phase : str
            Seismic phase identifier (``"P"`` or ``"S"``).
        source : CartHypocenter
            Hypocenter in local Cartesian coordinates.
        station : CartStation
            Station in local Cartesian coordinates.

        Returns
        -------
        dtx, dty, dtz : float
            Partial derivatives of traveltime with respect to the
            source coordinates (x, y, z), in s/m. If the analytical
            computation fails, all components are set to ``NaN``.
        """
        phase_u = phase.upper()
        if phase_u not in ("P", "S"):
            return math.nan, math.nan, math.nan

        branch = self._single_best_branch(phase_u, source, station)
        if branch is None:
            return math.nan, math.nan, math.nan

        p = float(branch.get("p", math.nan))
        if not math.isfinite(p):
            return math.nan, math.nan, math.nan

        z_s = self._get_depth(source)
        z_r = self._get_depth(station)
        kind = branch["kind"]
        if kind in ("head", "diving"):
            vertical_sign = 1.0
        elif z_s == z_r:
            vertical_sign = 0.0
        else:
            vertical_sign = 1.0 if z_r > z_s else -1.0

        q = branch["segments"][0].q

        dx = source.point.x - station.point.x
        dy = source.point.y - station.point.y
        offset = math.hypot(dx, dy)
        horizontal_scale = p / offset if offset else 0.0
        if self._height_above_model_top(source) > 0.0:
            vertical_sign = 0.0
        return dx * horizontal_scale, dy * horizontal_scale, vertical_sign * q

    def gradients(
        self,
        phase: str,
        source: CartHypocenter,
        station: CartStation,
    ) -> Tuple[float, float, float]:
        """
        Dispatch method for traveltime gradients.

        This method selects either the analytical or numerical
        implementation based on :attr:`gradient_method`. If the
        analytical computation fails (e.g. invalid branch or local
        velocity), it automatically falls back to the numerical
        finite-difference scheme.

        Parameters
        ----------
        phase : str
            Seismic phase identifier (``"P"`` or ``"S"``).
        source : CartHypocenter
            Hypocenter in local Cartesian coordinates.
        station : CartStation
            Station in local Cartesian coordinates.

        Returns
        -------
        dtx, dty, dtz : float
            Partial derivatives of traveltime with respect to the
            source coordinates (x, y, z), in s/m. If all attempts
            fail, the components are set to ``NaN``.
        """
        method = (self.gradient_method or "").strip().lower()

        if method.startswith("a"):
            dtx, dty, dtz = self.gradients_analytical(
                phase,
                source,
                station,
            )
            if (
                math.isfinite(dtx)
                and math.isfinite(dty)
                and math.isfinite(dtz)
            ):
                return dtx, dty, dtz
            # Fallback to numerical if analytical failed
            return self.gradients_numerical(phase, source, station)

        if method.startswith("n"):
            return self.gradients_numerical(phase, source, station)

        raise ValueError(
            "Unsupported gradient_method "
            f"'{self.gradient_method}'. Use 'analytical' or 'numerical'."
        )
