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
Layered Traveltime Solver (Extended Version)

This is an improved and partially refactored version of the original
tsolver.py module. It keeps the public API fully compatible with the
original implementation while adding:

  • full ray-path decomposition (segments) for all branch types
  • unified RaySegment structure
  • single-interface reflected paths in constant-velocity layers
  • better geometrical consistency for head-wave paths
  • improved numerical stability near critical angles
  • reorganized internal helpers (R3 restructuring level)

This file is intended to support analytical gradients for the upcoming
LayeredTravelTimeEngine in ShakeLab.

TO DO:
  • remove compatibility with legacy code
  • add switch to turn off segments
  • refactor for speed improvement
"""


from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple

import numpy as np

# Publicly exported names
__all__ = [
    "LayeredModel1D",
    "solve_traveltime",
    "compute_branches",
    "solve_with_heads",
    "curve_with_heads",
    "same_depth_branches",
    "same_depth_best",
    "same_depth_curve",
    "tt_direct",
    "direct_path",
    "reflected_path",
    "tt_direct_curve",
    "tt_direct_surface",
    "read_1d_model_profile",
]

# Solver defaults
TOLERANCE = 1e-6
MAX_ITER = 60


# ---------------------------------------------------------------------------
# Ray-segment container
# ---------------------------------------------------------------------------

@dataclass
class RaySegment:
    """
    Representation of a single segment of a ray path.

    For direct and diving rays, a segment corresponds to a straight path
    portion inside a single layer, defined by:
      • d_vert   vertical distance       (m)
      • v        velocity in layer        (m/s)
      • p        ray parameter            (s/m)
      • q        vertical slowness sqrt(1/v^2 - p^2)
      • x        horizontal projection    (m)
      • t        segment traveltime       (s)

    For head-wave paths:
      • type == "horizontal" is used for the along-interface segment
        with:
          d_vert = 0
          q = 0
          p = 1 / v_head
          x = horizontal length (L)
          t = L / v_head
    """

    type: str                 # "direct", "diving", "down", "up", "horizontal"
    layer: Optional[int]      # layer index, None for horizontal segment
    d_vert: float             # vertical distance (m)
    v: float                  # velocity (m/s)
    p: float                  # ray parameter (s/m)
    q: float                  # vertical slowness sqrt(...), 0 for horizontal
    x: float                  # horizontal projection (m)
    t: float                  # segment travel time (s)
    z_start: Optional[float] = None   # optional depth info for future use
    z_end: Optional[float] = None     # optional depth info for future use

    def as_dict(self) -> Dict:
        """Return as plain dictionary (for export/debug)."""
        return {
            "type": self.type,
            "layer": self.layer,
            "d_vert": self.d_vert,
            "v": self.v,
            "p": self.p,
            "q": self.q,
            "x": self.x,
            "t": self.t,
            "z_start": self.z_start,
            "z_end": self.z_end,
        }

# ---------------------------------------------------------------------------
# 1D Layered Earth Model
# ---------------------------------------------------------------------------

class LayeredModel1D:
    """
    Representation of a horizontally layered 1D model.

    Each layer is defined by:
        • thickness      (m)
        • vp             (P-wave velocity, m/s)
        • vs             (S-wave velocity, m/s)

    Depth increases downward. The last layer may have infinite thickness.

    This class is intentionally kept lightweight, because the traveltime
    solver relies on fast repeated geometric queries such as:
        • locating the layer containing a given depth
        • retrieving depth to top/bottom of layers
        • obtaining velocity arrays

    Attributes
    ----------
    thickness : ndarray
        Layer thicknesses (m).
    vp : ndarray
        P-wave velocities (m/s).
    vs : ndarray
        S-wave velocities (m/s).
    z_top : ndarray
        Depth to the top of each layer.
    z_bottom : ndarray
        Depth to the bottom of each layer.

    Notes
    -----
    • The model does not include density, Qp/Qs or gradients inside layers.
      Velocities are constant in each layer.
    • Velocities must be strictly positive.
    """

    def __init__(self, thickness, vp, vs):
        self.thickness = np.asarray(thickness, dtype=float)
        self.vp = np.asarray(vp, dtype=float)
        self.vs = np.asarray(vs, dtype=float)

        if not (self.thickness.size ==
                self.vp.size == self.vs.size):
            raise ValueError(
                "thickness, vp and vs must have the same length."
            )

        if np.any(self.thickness <= 0.0):
            raise ValueError("Layer thicknesses must be positive.")

        if np.any(self.vp <= 0.0) or np.any(self.vs <= 0.0):
            raise ValueError("Velocities must be positive.")

        # Precompute depth to top and bottom of layers
        self.z_top = np.zeros_like(self.thickness)
        if self.thickness.size > 1:
            self.z_top[1:] = np.cumsum(self.thickness[:-1])

        self.z_bottom = self.z_top + self.thickness

    # ------------------------------------------------------------------
    # Model geometry helpers
    # ------------------------------------------------------------------

    @property
    def num_layers(self) -> int:
        return self.thickness.size

    def layer_at(self, depth: float) -> int:
        """
        Return index of layer containing the given depth.

        Depth is positive downward. If last layer has infinite thickness,
        depths below model remain in that layer.
        """
        depth = float(depth)
        if depth < 0.0:
            raise ValueError("Depth must be non-negative.")

        for k in range(self.num_layers):
            if depth < self.z_bottom[k]:
                return k

        # Depth below all finite layers:
        if math.isinf(self.thickness[-1]):
            return self.num_layers - 1

        raise ValueError("Depth below bottom of finite model.")

    def total_thickness(self) -> float:
        """Return total thickness (may be infinite)."""
        return float(np.sum(self.thickness))


# ---------------------------------------------------------------------------
# Read model from CSV (unchanged API)
# ---------------------------------------------------------------------------

def read_1d_model_profile(filename: str) -> Tuple[
    List[float], List[float], List[float], List[float]
]:
    """
    Read layered model from CSV.

    Expected columns:
        "thickness (m)", "vp (m/s)", "vs (m/s)", "density (kg/m3)"

    Density is returned but not used by this solver.

    Use :func:`shakelab.location.profiles.read_model` to load a model
    directly from a JSON or CSV file.
    """
    from .profiles import _read_csv_columns

    with open(filename, "r", newline="") as f:
        return _read_csv_columns(f, with_density=True)


# ---------------------------------------------------------------------------
# Helper functions: geometry and numerical primitives
# ---------------------------------------------------------------------------

def _segment_lengths(model: LayeredModel1D,
                     z1: float,
                     z2: float) -> np.ndarray:
    """
    Return vertical path length within each layer between depths z1 and z2.

    This returns an array d_vert[k], with the amount of vertical distance
    travelled inside layer k when connecting z1 → z2 with a straight segment.

    No turning is assumed here (pure monotonic vertical path).
    """
    z1 = float(z1)
    z2 = float(z2)
    zmin = min(z1, z2)
    zmax = max(z1, z2)

    d_vert = np.zeros(model.num_layers, dtype=float)
    for k in range(model.num_layers):
        top = model.z_top[k]
        bottom = model.z_bottom[k]
        seg_top = max(zmin, top)
        seg_bottom = min(zmax, bottom)
        if seg_bottom > seg_top:
            d_vert[k] = seg_bottom - seg_top
    return d_vert


# ---------------------------------------------------------------------------
# Numerical stability helper
# ---------------------------------------------------------------------------

def _safe_q(p: float, v: float) -> float:
    """
    Return the vertical direction cosine with a numerical floor.

    If arguments would make the square root invalid or numerically unstable,
    use a small relative floor. This matches the behavior of the original
    solver but is encapsulated in a cleaner helper.
    """
    pv = p * v
    arg = 1.0 - pv * pv

    # Relative stability floor to avoid q → 0 causing blow-ups in dx, dt.
    eps_rel = 1e-12 + 1e-10 * abs(pv * pv)
    if arg <= -eps_rel:
        return math.nan
    arg = max(arg, eps_rel)

    return math.sqrt(arg)


# ---------------------------------------------------------------------------
# Core function: evaluate x(p), t(p) and build per-layer segments
# ---------------------------------------------------------------------------

def _xt(
    p: float,
    d_vert: np.ndarray,
    v_arr: np.ndarray,
    z_start: float,
    direction: str,
    *,
    descending: bool,
    cosines: Optional[np.ndarray] = None,
) -> Tuple[float, float, List[RaySegment]]:
    """
    Evaluate horizontal offset and traveltime for a ray with parameter p.

    Parameters
    ----------
    p : float
        Ray parameter (s/m).
    d_vert : ndarray
        Vertical distance in each crossed layer (m).
    v_arr : ndarray
        Velocities in crossed layers (m/s).
    z_start : float
        Starting depth of this segment (m).
    direction : str
        Segment type marker.
    descending : bool
        Whether depth increases along the ray.
    cosines : ndarray, optional
        Precomputed vertical direction cosines indexed by model layer.

    Returns
    -------
    x_total : float
        Total horizontal projection (m).
    t_total : float
        Total travel time (s).
    segments : list of RaySegment
        One segment per crossed layer.
    """
    x_total = 0.0
    t_total = 0.0
    segments: List[RaySegment] = []

    z_curr = float(z_start)

    layers = np.flatnonzero(d_vert > 0.0)
    if not descending:
        layers = layers[::-1]
    for k in layers:
        dz = d_vert[k]
        v = float(v_arr[k])
        q = _safe_q(p, v) if cosines is None else float(cosines[k])

        if not math.isfinite(q):
            return math.inf, math.inf, []

        # Horizontal projection for this layer
        x_k = dz * p * v / q
        t_k = dz / (v * q)

        x_total += x_k
        t_total += t_k

        z_next = z_curr + (dz if descending else -dz)

        segments.append(
            RaySegment(
                type=direction,
                layer=int(k),
                d_vert=dz,
                v=v,
                p=p,
                q=q / v,
                x=x_k,
                t=t_k,
                z_start=z_curr,
                z_end=z_next,
            )
        )

        z_curr = z_next

    return x_total, t_total, segments


# ---------------------------------------------------------------------------
# Helper: detect direction of travel (used for z_start / z_end)
# ---------------------------------------------------------------------------

def z2_is_deeper(z_start: float, z_curr: float) -> bool:
    """
    Return True if ray is travelling downward, False if upward.

    In a layered model this is simply determined by comparing depths.
    Only used to set approximate z_start/z_end for segments.
    """
    return z_curr >= z_start

# ---------------------------------------------------------------------------
# Direct branch solver (rewritten with full ray-path segments)
# ---------------------------------------------------------------------------

def _straight_direct_branch(
    layer: int, velocity: float, offset: float, z_s: float, z_r: float,
) -> Dict:
    """Build the exact direct ray for endpoints in a single layer."""
    vertical = abs(z_r - z_s)
    length = math.hypot(offset, vertical)
    p = offset / (length * velocity) if length else 0.0
    q = vertical / (length * velocity) if length else 0.0
    time = length / velocity
    kind = "direct" if vertical else "direct_horizontal"
    return {
        "kind": kind,
        "t": time,
        "p": p,
        "theta_source": math.atan2(offset, vertical),
        "x": offset,
        "x_target": offset,
        "dx": 0.0,
        "n_iter": 0,
        "converged": True,
        "segments": [
            RaySegment(
                type=kind, layer=layer, d_vert=vertical, v=velocity,
                p=p, q=q, x=offset, t=time, z_start=z_s, z_end=z_r,
            ),
        ],
    }


def _solve_direct_basic(model: LayeredModel1D,
                        offset: float,
                        z_s: float,
                        z_r: float,
                        phase: str,
                        tol: float = TOLERANCE,
                        max_iter: int = MAX_ITER) -> Optional[Dict]:
    """
    Solve the direct refracted branch between (z_s) and (z_r).

    Bisect the ray angle tangent in the fastest crossed layer to match
    horizontal offset, including nearly horizontal rays.
    Returns ray-path segments for full geometrical reconstruction.
    """

    offset = float(offset)
    if offset < 0.0:
        raise ValueError("Offset must be non-negative.")

    # Select velocity array
    phase_u = phase.upper()
    if phase_u == "P":
        v_arr_full = model.vp
    elif phase_u == "S":
        v_arr_full = model.vs
    else:
        raise ValueError("phase must be 'P' or 'S'.")

    # Build vertical crossing pattern
    d_vert_full = _segment_lengths(model, z_s, z_r)
    crossed = d_vert_full > 0.0
    layers = np.flatnonzero(crossed)
    if z_s == z_r or (
        layers.size == 1 and math.isclose(
            d_vert_full[layers[0]], abs(z_r - z_s), abs_tol=1e-9,
        )
    ):
        layer = model.layer_at(z_s) if z_s == z_r else int(layers[0])
        return _straight_direct_branch(
            layer, float(v_arr_full[layer]), offset, z_s, z_r,
        )
    if not layers.size:
        return None

    d_vert = d_vert_full[crossed]
    v_arr = v_arr_full[crossed]

    if offset == 0.0:
        t_vert = float(np.sum(d_vert / v_arr))
        _, _, segs = _xt(
            0.0, d_vert_full, v_arr_full, z_s, "direct",
            descending=z_r > z_s,
        )
        return {
            "kind": "direct",
            "t": t_vert,
            "p": 0.0,
            "theta_source": 0.0,
            "x": 0.0,
            "x_target": 0.0,
            "dx": 0.0,
            "n_iter": 0,
            "converged": True,
            "segments": segs,
        }

    v_max = float(np.max(v_arr))
    ratios = v_arr / v_max
    contrast = np.sqrt((1.0 - ratios) * (1.0 + ratios))
    # The other layers have smaller or equal ray angle tangents.
    tan_lo = offset / float(np.sum(d_vert))
    tan_hi = offset / float(np.sum(d_vert[v_arr == v_max]))
    tol_x = tol * max(1.0, abs(offset))

    for it in range(max_iter):
        tangent = math.sqrt(tan_lo) * math.sqrt(tan_hi)
        cos_fast = 1.0 / math.hypot(1.0, tangent)
        sin_fast = tangent * cos_fast
        # Retain cos_fast even when sin_fast rounds to one.
        cosines = np.hypot(cos_fast, sin_fast * contrast)
        x_mid = float(np.sum(d_vert * (sin_fast * ratios / cosines)))
        if abs(x_mid - offset) <= tol_x:
            break
        if x_mid < offset:
            tan_lo = tangent
        else:
            tan_hi = tangent
    else:
        return None

    p = sin_fast / v_max
    cosines_full = np.zeros_like(d_vert_full)
    cosines_full[crossed] = cosines
    x, time, segments = _xt(
        p, d_vert_full, v_arr_full, z_s, "direct",
        descending=z_r > z_s, cosines=cosines_full,
    )

    return {
        "kind": "direct",
        "t": float(time),
        "p": p,
        "theta_source": math.atan2(p, segments[0].q),
        "x": float(x),
        "x_target": float(offset),
        "dx": float(x - offset),
        "n_iter": it + 1,
        "converged": True,
        "segments": segments,
    }


# ---------------------------------------------------------------------------
# Public direct-path wrapper (API-compatible)
# ---------------------------------------------------------------------------

def direct_path(model: LayeredModel1D,
                offset: float,
                source_depth: float,
                receiver_depth: float = 0.0,
                phase: str = "P",
                tol: float = TOLERANCE,
                max_iter: int = MAX_ITER) -> Dict:
    """
    Solve direct-branch path. Returns full branch dictionary with segments.
    """
    branch = _solve_direct_basic(
        model=model,
        offset=offset,
        z_s=source_depth,
        z_r=receiver_depth,
        phase=phase,
        tol=tol,
        max_iter=max_iter,
    )
    if branch is None:
        raise ValueError("Direct branch unreachable or non-convergent.")
    return branch


# ---------------------------------------------------------------------------
# Traveltime-only wrappers (API unchanged)
# ---------------------------------------------------------------------------

def tt_direct(model: LayeredModel1D,
              offset: float,
              source_depth: float,
              receiver_depth: float = 0.0,
              phase: str = "P",
              tol: float = TOLERANCE,
              max_iter: int = MAX_ITER) -> float:
    """Return traveltime for direct branch only."""
    return direct_path(
        model=model,
        offset=offset,
        source_depth=source_depth,
        receiver_depth=receiver_depth,
        phase=phase,
        tol=tol,
        max_iter=max_iter,
    )["t"]


def tt_direct_surface(model: LayeredModel1D,
                      offset: float,
                      source_depth: float,
                      phase: str = "P",
                      tol: float = TOLERANCE,
                      max_iter: int = MAX_ITER) -> float:
    """Direct traveltime with receiver at the surface."""
    return tt_direct(
        model=model,
        offset=offset,
        source_depth=source_depth,
        receiver_depth=0.0,
        phase=phase,
        tol=tol,
        max_iter=max_iter,
    )


def tt_direct_curve(model: LayeredModel1D,
                    offsets,
                    source_depth: float,
                    receiver_depth: float = 0.0,
                    phase: str = "P",
                    tol: float = TOLERANCE,
                    max_iter: int = MAX_ITER) -> np.ndarray:
    """Vectorized direct-branch traveltime curve."""
    offs = np.asarray(offsets, dtype=float)
    t = np.empty_like(offs, dtype=float)
    it = np.nditer(offs, flags=["multi_index"])
    for x in it:
        idx = it.multi_index
        t[idx] = tt_direct(
            model=model,
            offset=float(x),
            source_depth=source_depth,
            receiver_depth=receiver_depth,
            phase=phase,
            tol=tol,
            max_iter=max_iter,
        )
    return t

# ---------------------------------------------------------------------------
# Single-interface reflected branch solver
# ---------------------------------------------------------------------------

def _solve_reflected_branch(model: LayeredModel1D,
                         offset: float,
                         z_s: float,
                         z_r: float,
                         phase: str,
                         iface_index: int,
                         tol: float = TOLERANCE,
                         max_iter: int = MAX_ITER) -> Optional[Dict]:
    """Solve a same-phase reflection from a deeper interface."""

    # Geometry check
    if iface_index <= 0 or iface_index >= model.num_layers:
        return None

    z_iface = float(model.z_top[iface_index])
    if not (z_iface > z_s and z_iface > z_r):
        return None

    # Select velocity array for chosen phase
    phase_u = phase.upper()
    if phase_u == "P":
        v_arr_full = model.vp
    else:
        v_arr_full = model.vs

    # Build descending and ascending vertical segments
    d_down_full = _segment_lengths(model, z_s, z_iface)
    d_up_full = _segment_lengths(model, z_r, z_iface)

    crossed_down = d_down_full > 0.0
    crossed_up = d_up_full > 0.0

    if not np.any(crossed_down) or not np.any(crossed_up):
        return None

    d_down = d_down_full[crossed_down]
    v_down = v_arr_full[crossed_down]
    d_up = d_up_full[crossed_up]
    v_up = v_arr_full[crossed_up]

    # Special case: zero offset → vertical composite path
    if offset == 0.0:
        t_vert = float(np.sum(d_down / v_down) + np.sum(d_up / v_up))
        _, _, seg_down = _xt(
            0.0, d_down_full, v_arr_full, z_s, "down", descending=True,
        )
        _, _, seg_up = _xt(
            0.0, d_up_full, v_arr_full, z_iface, "up", descending=False,
        )

        return {
            "kind": "reflected",
            "iface": int(iface_index),
            "t": t_vert,
            "p": 0.0,
            "theta_source": 0.0,
            "x": 0.0,
            "x_target": 0.0,
            "dx": 0.0,
            "n_iter": 0,
            "converged": True,
            "reflection_depth": z_iface,
            "segments": seg_down + seg_up,
        }

    # General case: bisection on p
    p_lo = 0.0
    p_hi = 1.0 / float(max(np.max(v_down), np.max(v_up)))
    tol_x = tol * max(1.0, abs(offset))

    converged = False
    best_segments = []
    p_mid = None
    x_mid = None
    t_mid = None

    for it in range(max_iter):
        p_mid = 0.5 * (p_lo + p_hi)

        # Evaluate descending
        x_d, t_d, seg_down = _xt(
            p=p_mid,
            d_vert=d_down_full,
            v_arr=v_arr_full,
            z_start=z_s,
            direction="down",
            descending=True,
        )
        if not math.isfinite(x_d):
            break

        # Evaluate ascending
        x_u, t_u, seg_up = _xt(
            p=p_mid,
            d_vert=d_up_full,
            v_arr=v_arr_full,
            z_start=z_iface,
            direction="up",
            descending=False,
        )
        if not math.isfinite(x_u):
            break

        x_mid = x_d + x_u
        t_mid = t_d + t_u

        if abs(x_mid - offset) <= tol_x:
            converged = True
            best_segments = seg_down + seg_up
            break

        if x_mid < offset:
            p_lo = p_mid
        else:
            p_hi = p_mid

        best_segments = seg_down + seg_up

    if not converged:
        return None

    # Take-off angle
    k_s = model.layer_at(z_s)
    v_s = float(v_arr_full[k_s])
    pv_s = max(-1.0, min(1.0, p_mid * v_s))
    theta_s = math.asin(pv_s)

    return {
        "kind": "reflected",
        "iface": int(iface_index),
        "t": float(t_mid),
        "p": float(p_mid),
        "theta_source": float(theta_s),
        "x": float(x_mid),
        "x_target": float(offset),
        "dx": float(x_mid - offset),
        "n_iter": it + 1,
        "converged": converged,
        "reflection_depth": z_iface,
        "segments": best_segments,
    }


def reflected_path(
    model: LayeredModel1D,
    offset: float,
    source_depth: float,
    receiver_depth: float = 0.0,
    *,
    interface_index: int,
    phase: str = "P",
    tol: float = TOLERANCE,
    max_iter: int = MAX_ITER,
) -> dict:
    """
    Compute a same-phase reflection from a specified interface.

    Parameters
    ----------
    model : LayeredModel1D
        Model with constant P and S velocities in each layer.
    offset : float
        Non-negative horizontal distance in meters.
    source_depth : float
        Source depth in meters, positive downward.
    receiver_depth : float, optional
        Receiver depth in meters, positive downward.
    interface_index : int
        Index of the layer whose top is the reflector, with
        ``1`` denoting the interface below the first layer.
        Both endpoints must be strictly above the interface.
    phase : str, optional
        Wave phase, either ``"P"`` (P-to-P) or ``"S"`` (S-to-S).
    tol : float, optional
        Offset residual tolerance scaled by max(1 m, offset).
    max_iter : int, optional
        Maximum number of bisection iterations.

    Returns
    -------
    branch : dict
        Branch of kind ``"reflected"`` with travel time ``t`` in seconds,
        ray parameter ``p`` in s/m, ``iface``, ``reflection_depth`` in
        meters and a list of ray ``segments``.

    Raises
    ------
    TypeError
        If the interface index is not an integer.
    ValueError
        If the geometry or options are invalid, or the solve fails.

    Notes
    -----
    Source and receiver can share a depth. The result describes
    ray geometry and travel time. Reflection amplitudes and mode
    conversions are outside the scope of this function.
    """
    if isinstance(interface_index, (bool, np.bool_)):
        raise TypeError("interface_index must be an integer.")
    interface_index = operator.index(interface_index)
    if not 1 <= interface_index < model.num_layers:
        raise ValueError("interface_index must identify an internal boundary.")

    geometry = (offset, source_depth, receiver_depth)
    if not all(math.isfinite(value) and value >= 0 for value in geometry):
        raise ValueError(
            "Offset and endpoint depths must be finite and non-negative."
        )
    if not max(source_depth, receiver_depth) < model.z_top[interface_index]:
        raise ValueError("Source and receiver must be above the reflector.")
    if phase.upper() not in ("P", "S"):
        raise ValueError("phase must be 'P' or 'S'.")
    if not math.isfinite(tol) or tol <= 0.0 or max_iter <= 0:
        raise ValueError("tol and max_iter must be positive.")

    branch = _solve_reflected_branch(
        model, offset, source_depth, receiver_depth, phase,
        interface_index, tol=tol, max_iter=max_iter,
    )
    if branch is None:
        raise ValueError("Reflected branch did not converge.")
    return branch


# ---------------------------------------------------------------------------
# Head-wave solver (rewritten with explicit segments)
# ---------------------------------------------------------------------------

def _solve_head_at(model: LayeredModel1D,
                   offset: float,
                   z_s: float,
                   z_r: float,
                   phase: str,
                   iface_index: int) -> Optional[Dict]:
    """
    Solve the head-wave branch generated at interface iface_index.

    Returns full ray-path segments:
        • down-slant (source → interface)
        • horizontal (along interface at critical slowness)
        • up-slant (interface → receiver)

    Valid only if:
        – v_lower > v_upper  (velocity increase)
        – both endpoints are at or above the interface, with at least
          one strictly above it
        – offset > x_down + x_up   (positive horizontal residual)
    """

    # Sanity
    if offset <= 0.0:
        return None

    if iface_index <= 0 or iface_index >= model.num_layers:
        return None

    z_iface = float(model.z_top[iface_index])

    if max(z_s, z_r) > z_iface or min(z_s, z_r) == z_iface:
        return None

    phase_u = phase.upper()
    v_arr_full = model.vp if phase_u == "P" else model.vs

    v_upper = float(v_arr_full[iface_index - 1])
    v_lower = float(v_arr_full[iface_index])

    # Must have velocity increase
    if v_lower <= v_upper:
        return None

    # Critical slowness
    p_c = 1.0 / v_lower

    d_down_full = _segment_lengths(model, z_s, z_iface)
    d_up_full = _segment_lengths(model, z_r, z_iface)
    crossed = (d_down_full > 0.0) | (d_up_full > 0.0)
    crossed_velocities = v_arr_full[crossed]
    if np.any(crossed_velocities >= v_lower):
        return None
    ratios = crossed_velocities / v_lower
    cosines = np.zeros_like(d_down_full)
    cosines[crossed] = np.sqrt((1.0 - ratios) * (1.0 + ratios))

    x_down, t_down, seg_down = _xt(
        p=p_c,
        d_vert=d_down_full,
        v_arr=v_arr_full,
        z_start=z_s,
        direction="down",
        descending=True,
        cosines=cosines,
    )

    x_up, t_up, seg_up = _xt(
        p=p_c,
        d_vert=d_up_full,
        v_arr=v_arr_full,
        z_start=z_iface,
        direction="up",
        descending=False,
        cosines=cosines,
    )

    # Horizontal residual
    L = float(offset) - (x_down + x_up)
    if L <= 0.0:
        return None

    # Time of horizontal interface segment
    t_h = L / v_lower
    x_h = L

    seg_h = RaySegment(
        type="horizontal",
        layer=iface_index,
        d_vert=0.0,
        v=v_lower,
        p=p_c,
        q=0.0,
        x=x_h,
        t=t_h,
        z_start=z_iface,
        z_end=z_iface,
    )

    segments = seg_down + [seg_h] + seg_up

    return {
        "kind": "head",
        "iface": int(iface_index),
        "t": float(t_down + t_h + t_up),
        "p": float(p_c),
        "theta_source": math.atan2(p_c, segments[0].q),
        "x": float(offset),
        "x_target": float(offset),
        "dx": 0.0,
        "n_iter": 0,
        "converged": True,
        "segments": segments,
    }


# ---------------------------------------------------------------------------
# Enumerate all valid head-wave branches for a given offset
# ---------------------------------------------------------------------------

def _solve_all_heads(model: LayeredModel1D,
                     offset: float,
                     z_s: float,
                     z_r: float,
                     phase: str) -> List[Dict]:
    """
    Build all head-wave branches generated at velocity-increasing interfaces.
    """
    phase_u = phase.upper()
    v_arr = model.vp if phase_u == "P" else model.vs

    branches = []
    for iface in range(1, model.num_layers):
        if v_arr[iface] > v_arr[iface - 1]:
            br = _solve_head_at(
                model=model,
                offset=offset,
                z_s=z_s,
                z_r=z_r,
                phase=phase,
                iface_index=iface,
            )
            if br is not None:
                branches.append(br)

    return branches

# ---------------------------------------------------------------------------
# Combined branch computation (direct + diving + head)
# ---------------------------------------------------------------------------

def compute_branches(model: LayeredModel1D,
                     offset: float,
                     source_depth: float,
                     receiver_depth: float = 0.0,
                     phase: str = "P",
                     include_direct: bool = True,
                     include_diving: bool = True,
                     include_heads: bool = True,
                     tol: float = TOLERANCE,
                     max_iter: int = MAX_ITER,
                     require_convergence: bool = True) -> List[Dict]:
    """
    Build candidate branches for a single offset.

    The legacy ``include_diving`` option adds interface reflections
    under the branch label ``"diving"``. Use :func:`reflected_path` to
    request a specific reflection, including at equal endpoint depths.

    Returns a list of branch dictionaries. Each branch has:
        kind: "direct" | "diving" | "head"
        t: travel time
        p: ray parameter
        theta_source
        x, x_target, dx
        n_iter
        converged
        segments: List[RaySegment]
        extra fields:
            iface        (for diving and head)
            turning_depth (for diving)
    """

    branches = []

    # -------------------------------------------------------------
    # 1) Direct branch
    # -------------------------------------------------------------
    if include_direct:
        direct = _solve_direct_basic(
            model=model,
            offset=offset,
            z_s=source_depth,
            z_r=receiver_depth,
            phase=phase,
            tol=tol,
            max_iter=max_iter,
        )
        if direct is not None:
            branches.append(direct)

    # -------------------------------------------------------------
    # 2) Reflections through the legacy include_diving option
    # -------------------------------------------------------------
    if include_diving:
        # Keep the legacy enumeration rule for this option.
        if source_depth != receiver_depth:
            for iface in range(1, model.num_layers):
                z_iface = model.z_top[iface]
                # only deeper than both source/receiver
                if z_iface > max(source_depth, receiver_depth):
                    br = _solve_reflected_branch(
                        model=model,
                        offset=offset,
                        z_s=source_depth,
                        z_r=receiver_depth,
                        phase=phase,
                        iface_index=iface,
                        tol=tol,
                        max_iter=max_iter,
                    )
                    if br is not None:
                        br["kind"] = "diving"
                        br["turning_depth"] = br.pop("reflection_depth")
                        if offset == 0.0:
                            for segment in br["segments"]:
                                segment.type = "diving"
                        branches.append(br)

    # -------------------------------------------------------------
    # 3) Head waves
    # -------------------------------------------------------------
    if include_heads and offset > 0.0:
        heads = _solve_all_heads(
            model=model,
            offset=offset,
            z_s=source_depth,
            z_r=receiver_depth,
            phase=phase,
        )
        branches.extend(heads)

    # -------------------------------------------------------------
    # Filtering for convergence / residual quality
    # -------------------------------------------------------------
    if require_convergence:
        tol_x = tol * max(1.0, abs(float(offset)))

        filtered = []
        for b in branches:
            if not b.get("converged", True):
                continue
            if abs(b.get("dx", 0.0)) > tol_x:
                continue
            filtered.append(b)

        branches = filtered

    return branches

# ---------------------------------------------------------------------------
# Fastest-arrival solver for one offset
# ---------------------------------------------------------------------------

def solve_with_heads(model: LayeredModel1D,
                     offset: float,
                     source_depth: float,
                     receiver_depth: float = 0.0,
                     phase: str = "P",
                     tol: float = TOLERANCE,
                     max_iter: int = MAX_ITER,
                     require_convergence: bool = True,
                     include_diving: bool = True) -> Dict:
    """
    Select the fastest branch among direct, diving and head waves.

    Uses compute_branches() internally, which already provides R3
    segments and stability filtering.
    """

    branches = compute_branches(
        model=model,
        offset=offset,
        source_depth=source_depth,
        receiver_depth=receiver_depth,
        phase=phase,
        include_direct=True,
        include_diving=include_diving,
        include_heads=True,
        tol=tol,
        max_iter=max_iter,
        require_convergence=require_convergence,
    )

    if not branches:
        raise ValueError("No valid branch for this geometry.")

    return min(branches, key=lambda b: b["t"])


# ---------------------------------------------------------------------------
# Curve with heads (vectorized fastest-arrival for many offsets)
# ---------------------------------------------------------------------------

def curve_with_heads(model: LayeredModel1D,
                     offsets,
                     source_depth: float,
                     receiver_depth: float = 0.0,
                     phase: str = "P",
                     tol: float = TOLERANCE,
                     max_iter: int = MAX_ITER,
                     require_convergence: bool = True,
                     include_diving: bool = True) -> np.ndarray:
    """
    Vectorized version of solve_with_heads(), returns T(x) for all offsets.
    """
    offs = np.asarray(offsets, dtype=float)
    t = np.empty_like(offs, dtype=float)

    it = np.nditer(offs, flags=["multi_index"])
    for x in it:
        idx = it.multi_index
        best = solve_with_heads(
            model=model,
            offset=float(x),
            source_depth=source_depth,
            receiver_depth=receiver_depth,
            phase=phase,
            tol=tol,
            max_iter=max_iter,
            require_convergence=require_convergence,
            include_diving=include_diving,
        )
        t[idx] = best["t"]

    return t


# ---------------------------------------------------------------------------
# Same-depth geometry helpers
# ---------------------------------------------------------------------------

def same_depth_branches(model: LayeredModel1D,
                        offset: float,
                        depth: float,
                        phase: str = "P") -> List[Dict]:
    """
    Build candidate branches for source and receiver at same depth.

    Includes:
      - direct horizontal branch
      - head waves (with full R3 segment structure)
    """

    return compute_branches(
        model=model, offset=offset, source_depth=depth, receiver_depth=depth,
        phase=phase, include_diving=False,
    )


# ---------------------------------------------------------------------------
# Fastest branch for same-depth geometry
# ---------------------------------------------------------------------------

def same_depth_best(model: LayeredModel1D,
                    offset: float,
                    depth: float,
                    phase: str = "P") -> Dict:
    """Select the fastest same-depth branch."""
    branches = same_depth_branches(
        model=model,
        offset=offset,
        depth=depth,
        phase=phase,
    )
    if not branches:
        raise ValueError("No branch for same-depth configuration.")
    return min(branches, key=lambda b: b["t"])


# ---------------------------------------------------------------------------
# Same-depth curve (vectorized)
# ---------------------------------------------------------------------------

def same_depth_curve(model: LayeredModel1D,
                     offsets,
                     depth: float,
                     phase: str = "P") -> np.ndarray:
    """Vectorized version of same_depth_best() for T(x) curve."""
    offs = np.asarray(offsets, dtype=float)
    t = np.empty_like(offs, dtype=float)

    it = np.nditer(offs, flags=["multi_index"])
    for x in it:
        idx = it.multi_index
        best = same_depth_best(
            model=model,
            offset=float(x),
            depth=depth,
            phase=phase,
        )
        t[idx] = best["t"]

    return t

# ---------------------------------------------------------------------------
# High-level traveltime solver (public API)
# ---------------------------------------------------------------------------

def solve_traveltime(model: LayeredModel1D,
                     offsets,
                     source_depth: float,
                     receiver_depth: float = 0.0,
                     phase: str = "P",
                     arrivals: str = "all",
                     return_branches: bool = False,
                     tol: float = TOLERANCE,
                     max_iter: int = MAX_ITER,
                     require_convergence: bool = True,
                     include_diving: bool = True) -> Tuple[np.ndarray, List]:
    """
    Compute traveltime(s) for given offset(s) and a single source depth.

    Parameters
    ----------
    model : LayeredModel1D
        Layered velocity model.
    offsets : float or array-like
        Horizontal distance(s) between source and receiver(s), meters.
    source_depth : float
        Source depth (m, positive downward).
    receiver_depth : float
        Receiver depth (m, normally 0 for surface).
    phase : str
        'P' or 'S'.
    arrivals : {'direct', 'all'}
        • 'direct' → only direct arrival (may fail if not reachable)
        • 'all'    → fastest arrival among direct / diving / head
    return_branches : bool
        If True, return also the list of branches for each offset.
        If False, return only the traveltimes.
    tol : float
        Horizontal residual tolerance for convergence.
    max_iter : int
        Maximum iterations for direct and reflected ray searches.
    require_convergence : bool
        Filter branches with residual dx above tolerance.
    include_diving : bool
        Whether to include interface reflections with the legacy label
        ``"diving"``.

    Returns
    -------
    tt : ndarray
        Traveltimes for each offset.
    branches_all : list
        If return_branches=True, a list of lists of branches.
        Otherwise, an empty list.
    """

    offs = np.asarray(offsets, dtype=float)
    shape = offs.shape
    offs_flat = offs.ravel()

    tt_flat = np.zeros_like(offs_flat, dtype=float)
    branches_all = [] if return_branches else None

    for i, x in enumerate(offs_flat):

        if arrivals.lower() == "direct":
            # Direct arrival only
            b = _solve_direct_basic(
                model=model,
                offset=float(x),
                z_s=source_depth,
                z_r=receiver_depth,
                phase=phase,
                tol=tol,
                max_iter=max_iter,
            )
            if b is None or not b.get("converged", True):
                raise ValueError(
                    f"Direct branch did not converge for offset {x} m."
                )
            tt_flat[i] = b["t"]
            if return_branches:
                branches_all.append([b])

        else:
            # Full set: direct + diving + head
            branches = compute_branches(
                model=model,
                offset=float(x),
                source_depth=source_depth,
                receiver_depth=receiver_depth,
                phase=phase,
                include_direct=True,
                include_diving=include_diving,
                include_heads=True,
                tol=tol,
                max_iter=max_iter,
                require_convergence=require_convergence,
            )

            if not branches:
                raise ValueError(
                    f"No valid branch for offset {x} m."
                )

            best = min(branches, key=lambda b: b["t"])
            tt_flat[i] = best["t"]

            if return_branches:
                branches_all.append(branches)

    tt = tt_flat.reshape(shape)
    return tt, branches_all if return_branches else []

# ---------------------------------------------------------------------------
# Footer, version tag, and retro-compatibility aliases
# ---------------------------------------------------------------------------

__version__ = "R3-extended-2025"

# Retro-compatibility:
# (the original solver exposed some names that now map to rewritten functions)

# Direct solvers
direct = tt_direct
direct_surface = tt_direct_surface
direct_curve = tt_direct_curve

# Full solver alias
traveltime = solve_traveltime

# Same-depth aliases
horizontal = same_depth_best
horizontal_curve = same_depth_curve


# Explicitly define exported API
__all__ = [
    # Main model class
    "LayeredModel1D",

    # Ray-segment data structure
    "RaySegment",

    # High-level solver
    "solve_traveltime",
    "traveltime",

    # Branch solvers
    "compute_branches",
    "solve_with_heads",
    "direct_path",
    "reflected_path",
    "direct",
    "tt_direct",
    "direct_surface",
    "direct_curve",

    # Same-depth solvers
    "same_depth_branches",
    "same_depth_best",
    "horizontal",
    "same_depth_curve",
    "horizontal_curve",

    # Curve solvers
    "curve_with_heads",

    # I/O
    "read_1d_model_profile",

    # Meta
    "__version__",
]
