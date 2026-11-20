"""Planar rectangular rupture geometry in a local WGS84 ENU frame.

Public lengths, depths and distances are in kilometres; WgsPoint elevation
and Cartesian coordinates returned by geodeticN are in metres. The finite
surface is a plane in one local tangent frame, not a curved Earth surface.
Distances are exact for that Cartesian representation. For large ruptures
or distant sites, the tangent-plane approximation needs reassessment.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, isclose, isfinite, radians, sin, sqrt
from numbers import Real

import numpy as np

from shakelab.groundmotion.scaling import ScalingEstimate

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    enu_to_geo,
    geo_to_enu,
)


def _number(value: float, name: str) -> float:
    """Return a finite real scalar."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number.")
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite.")
    return value


@dataclass(frozen=True)
class RectangularRupture:
    """Rectangular rupture anchored to an explicitly located WGS84 point.

    Args:
        reference: Point on the rupture, defining its local ENU frame.
        strike: Clockwise azimuth from North in degrees.
        dip: Inclination from horizontal in degrees, dipping to the
            right of strike.
        length: Along-strike length in km.
        width: Down-dip width in km.
        along_fraction: Reference position from the first upper vertex
            toward the second, as a fraction of length.
        down_fraction: Reference position from the upper edge toward
            the lower edge, as a fraction of width.

    By default the reference is the midpoint of the upper edge.
    Fractions are geometry coordinates, not probabilities.
    """

    reference: WgsPoint
    strike: float
    dip: float
    length: float
    width: float
    along_fraction: float = 0.5
    down_fraction: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.reference, WgsPoint):
            raise TypeError("reference must be a WgsPoint.")
        for name in ("strike", "dip", "length", "width",
                     "along_fraction", "down_fraction"):
            object.__setattr__(
                self, name, _number(getattr(self, name), name),
            )
        if not 0.0 <= self.strike <= 360.0:
            raise ValueError("strike must be in [0, 360] degrees.")
        if not 0.0 <= self.dip <= 90.0:
            raise ValueError("dip must be in [0, 90] degrees.")
        if self.length <= 0.0 or self.width <= 0.0:
            raise ValueError("length and width must be strictly positive.")
        if not 0.0 <= self.along_fraction <= 1.0:
            raise ValueError("along_fraction must be in [0, 1].")
        if not 0.0 <= self.down_fraction <= 1.0:
            raise ValueError("down_fraction must be in [0, 1].")
        if self.ztor < -1e-9:
            raise ValueError("Rupture upper edge is above the local surface.")

    @classmethod
    def from_hypocentre(
        cls,
        hypocentre: WgsPoint,
        *,
        strike: float,
        dip: float,
        length: float,
        width: float,
        ztor: float,
        along_fraction: float,
        down_fraction: float | None = None,
    ) -> RectangularRupture:
        """Construct a plane containing the hypocentre.

        ztor is the depth of the upper edge in the hypocentre's local
        tangent frame. For positive dip, the down-dip position follows
        from hypocentral depth and ztor. For horizontal faults it must
        be specified explicitly because depth cannot constrain it.
        """
        if not isinstance(hypocentre, WgsPoint):
            raise TypeError("hypocentre must be a WgsPoint.")
        dip = _number(dip, "dip")
        width = _number(width, "width")
        ztor = _number(ztor, "ztor")
        if not 0.0 <= dip <= 90.0 or width <= 0.0 or ztor < 0.0:
            raise ValueError("Invalid dip, width or ztor.")
        depth = -hypocentre.elevation / 1000.0
        relief = width * sin(radians(dip))
        if dip == 0.0:
            if abs(depth - ztor) > 1e-9:
                raise ValueError("Horizontal rupture must match hypocentre depth.")
            if down_fraction is None:
                raise ValueError("Horizontal rupture requires down_fraction.")
        else:
            fraction = (depth - ztor) / relief
            if fraction < -1e-9 or fraction > 1.0 + 1e-9:
                raise ValueError("Hypocentre is outside the rupture depth range.")
            fraction = min(1.0, max(0.0, fraction))
            if down_fraction is not None:
                supplied = _number(down_fraction, "down_fraction")
                if abs(supplied - fraction) > 1e-9:
                    raise ValueError("down_fraction is inconsistent with ztor.")
            down_fraction = fraction
        return cls(
            hypocentre, strike, dip, length, width,
            along_fraction, down_fraction,
        )

    @property
    def frame(self) -> MetricFrame:
        """Local East-North-Up frame centred on reference."""
        return MetricFrame(self.reference, orientation="enu")

    @property
    def strike_vector(self) -> np.ndarray:
        """Unit vector along strike."""
        angle = radians(self.strike)
        return np.array([sin(angle), cos(angle), 0.0])

    @property
    def cross_strike_vector(self) -> np.ndarray:
        """Horizontal unit vector to the right of strike."""
        angle = radians(self.strike)
        return np.array([cos(angle), -sin(angle), 0.0])

    @property
    def dip_vector(self) -> np.ndarray:
        """Unit vector down the rupture plane."""
        angle = radians(self.dip)
        return (cos(angle) * self.cross_strike_vector
                + np.array([0.0, 0.0, -sin(angle)]))

    @property
    def ztor(self) -> float:
        """Upper-edge depth in the local frame, in km."""
        return (-self.reference.elevation / 1000.0
                - self.down_fraction * self.width * sin(radians(self.dip)))

    @property
    def area(self) -> float:
        """Rupture area in square kilometres."""
        return self.length * self.width

    def point_enu(self, along: float, down: float) -> np.ndarray:
        """Return a point on the rectangle as ENU coordinates in km.

        along and down are fractions measured from the first upper
        vertex. Both must lie in [0, 1].
        """
        along = _number(along, "along")
        down = _number(down, "down")
        if not 0.0 <= along <= 1.0 or not 0.0 <= down <= 1.0:
            raise ValueError("Rupture coordinates must lie in [0, 1].")
        return ((along - self.along_fraction) * self.length
                * self.strike_vector
                + (down - self.down_fraction) * self.width
                * self.dip_vector)

    def point(self, along: float, down: float) -> WgsPoint:
        """Return a WGS84 point on the rectangle."""
        xyz = self.point_enu(along, down) * 1000.0
        return WgsPoint(*enu_to_geo(xyz.tolist(), self.frame))

    @property
    def vertices_enu(self) -> np.ndarray:
        """Four vertices in km: top-start, top-end, bottom-end, bottom-start."""
        return np.array([
            self.point_enu(0.0, 0.0),
            self.point_enu(1.0, 0.0),
            self.point_enu(1.0, 1.0),
            self.point_enu(0.0, 1.0),
        ])

    @property
    def vertices(self) -> tuple[WgsPoint, ...]:
        """Four WGS84 vertices in the same order as vertices_enu."""
        return tuple(self.point(a, d) for a, d in (
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),
        ))

    def distances_enu(self, sites: np.ndarray) -> dict[str, np.ndarray]:
        """Calculate finite-rupture distances for local ENU sites in km.

        Args:
            sites: Array of shape (n, 3), in this rupture's local frame.

        Returns:
            Mapping of arrays in km: rrup (minimum 3-D distance to the
            rectangle), rjb (minimum horizontal distance to its
            projected footprint), rx (signed cross-strike distance from
            the upper edge's infinite line, positive on its right),
            ry0 (horizontal along-strike distance beyond either end).
        """
        sites = np.asarray(sites, dtype=float)
        if sites.ndim != 2 or sites.shape[1] != 3:
            raise ValueError("sites must have shape (n, 3).")
        if not np.isfinite(sites).all():
            raise ValueError("sites must contain finite coordinates.")

        delta = sites - self.point_enu(0.0, 0.0)
        strike = self.strike_vector
        dip = self.dip_vector
        along = delta @ strike
        down = delta @ dip
        cross = delta @ self.cross_strike_vector

        nearest = (np.clip(along, 0.0, self.length)[:, None] * strike
                   + np.clip(down, 0.0, self.width)[:, None] * dip)
        rrup = np.linalg.norm(delta - nearest, axis=1)

        ry0 = np.maximum(np.maximum(-along, along - self.length), 0.0)
        horizontal_width = self.width * cos(radians(self.dip))
        excess = np.maximum(
            np.maximum(-cross, cross - horizontal_width), 0.0,
        )
        return {"rrup": rrup, "rjb": np.hypot(ry0, excess),
                "rx": cross, "ry0": ry0}

    def distances(self, sites: list[WgsPoint]) -> dict[str, np.ndarray]:
        """Calculate distances in km for WGS84 sites, using their elevations."""
        if not all(isinstance(site, WgsPoint) for site in sites):
            raise TypeError("sites must contain WgsPoint objects.")
        frame = self.frame
        coordinates = np.array(
            [geo_to_enu(site, frame) for site in sites], dtype=float,
        ).reshape(-1, 3) / 1000.0
        return self.distances_enu(coordinates)


@dataclass(frozen=True)
class RuptureBuildResult:
    """Generated geometry and provenance of resolved source parameters."""

    rupture: RectangularRupture
    source_parameters: dict[str, float]
    parameter_origins: dict[str, str]
    scaling_model: str | None = None
    scaling_style: str | None = None
    scaling_estimate: ScalingEstimate | None = None


def build_rectangular_rupture(
    hypocentre: WgsPoint,
    *,
    strike: float,
    dip: float,
    magnitude: float | None = None,
    rake: float | None = None,
    area: float | None = None,
    length: float | None = None,
    width: float | None = None,
    aspect_ratio: float | None = None,
    ztor: float | None = None,
    along_fraction: float | None = None,
    down_fraction: float | None = None,
    upper_seismogenic_depth: float = 0.0,
    lower_seismogenic_depth: float | None = None,
    scaling_model: str | None = None,
    scaling_style: str = "auto",
    scaling_epsilon: float = 0.0,
    scaling_strict: bool = False,
) -> RuptureBuildResult:
    """Resolve a rectangular rupture from incomplete source information.

    Depths and sizes are in km, area in km**2, angles in degrees.
    Magnitude must be Mw when scaling is used. Scaling is opt-in and
    requires a model name. At least one shape constraint is required
    (length, width or aspect_ratio), unless the depth and nucleation
    constraints determine width. Redundant explicit inputs must agree.

    Aspect ratio is a target when neither dimension is fixed: width
    can be adjusted to fit depth bounds while preserving area. When
    either dimension is fixed, an explicit aspect ratio is a constraint.
    The default nucleation is central; only an unspecified down-dip
    fraction can be shifted to accommodate the depth bounds.

    Parameter origins distinguish explicit, scaling, derived, default
    and constrained values. No focal angles or scaling model are
    inferred. An omitted upper depth is zero; the lower depth is
    unbounded unless provided.
    """
    from shakelab.groundmotion.scaling import (
        create_scaling_relation,
        fault_style_from_rake,
    )

    if not isinstance(hypocentre, WgsPoint):
        raise TypeError("hypocentre must be a WgsPoint.")
    strike = _number(strike, "strike")
    dip = _number(dip, "dip")
    if not 0 <= strike <= 360 or not 0 <= dip <= 90:
        raise ValueError("Invalid strike or dip.")
    if rake is not None:
        rake = _number(rake, "rake")
        fault_style_from_rake(rake)
    if magnitude is not None:
        magnitude = _number(magnitude, "magnitude")
    if not isinstance(scaling_strict, bool):
        raise TypeError("scaling_strict must be boolean.")
    epsilon = _number(scaling_epsilon, "scaling_epsilon")
    if not isinstance(scaling_style, str):
        raise TypeError("scaling_style must be a string.")

    values = {}
    origins = {}
    for name, value in (
        ("area", area), ("length", length), ("width", width),
        ("aspect_ratio", aspect_ratio), ("ztor", ztor),
        ("along_fraction", along_fraction),
        ("down_fraction", down_fraction),
    ):
        if value is not None:
            values[name] = _number(value, name)
            origins[name] = "explicit"
    area = values.get("area")
    length = values.get("length")
    width = values.get("width")
    aspect_ratio = values.get("aspect_ratio")
    ztor = values.get("ztor")
    along_fraction = values.get("along_fraction")
    down_fraction = values.get("down_fraction")
    for name in ("area", "length", "width", "aspect_ratio"):
        if name in values and values[name] <= 0:
            raise ValueError(f"{name} must be strictly positive.")
    for name in ("along_fraction", "down_fraction"):
        if name in values and not 0 <= values[name] <= 1:
            raise ValueError(f"{name} must be in [0, 1].")

    upper = _number(upper_seismogenic_depth, "upper_seismogenic_depth")
    lower = (float("inf") if lower_seismogenic_depth is None
             else _number(lower_seismogenic_depth,
                          "lower_seismogenic_depth"))
    depth = -hypocentre.elevation / 1000
    if upper < 0 or lower <= upper:
        raise ValueError("Invalid seismogenic depth bounds.")
    if not upper <= depth <= lower:
        raise ValueError("Hypocentre is outside the seismogenic layer.")
    if ztor is not None and not upper <= ztor <= depth:
        raise ValueError("ztor is outside the permitted depth range.")

    # Width interval implied by depth and nucleation constraints.
    slope = sin(radians(dip))
    min_width, max_width = 0.0, float("inf")
    fixed_width = None
    if slope == 0:
        if ztor is not None and not isclose(
            ztor, depth, rel_tol=0, abs_tol=1e-9,
        ):
            raise ValueError("Horizontal rupture requires ztor = depth.")
    elif ztor is not None:
        min_width = (depth - ztor) / slope
        max_width = (lower - ztor) / slope
        if down_fraction is not None:
            if down_fraction == 0:
                if depth != ztor:
                    raise ValueError("down_fraction conflicts with ztor.")
            else:
                fixed_width = min_width / down_fraction
                if fixed_width <= 0:
                    raise ValueError("ztor and down_fraction imply zero width.")
    elif down_fraction is not None:
        if down_fraction > 0:
            max_width = (depth - upper) / (slope * down_fraction)
        if down_fraction < 1:
            max_width = min(
                max_width, (lower - depth) /
                (slope * (1 - down_fraction)),
            )
    else:
        max_width = (lower - upper) / slope

    def agrees(a, b):
        return isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)

    if fixed_width is not None:
        if width is not None and not agrees(width, fixed_width):
            raise ValueError("width conflicts with ztor and down_fraction.")
        if width is None:
            width = fixed_width
            origins["width"] = "derived"

    # An explicit aspect ratio with a known dimension fixes the other.
    if aspect_ratio is not None:
        if length is not None and width is not None:
            if not agrees(length / width, aspect_ratio):
                raise ValueError("aspect_ratio conflicts with dimensions.")
        elif length is not None:
            width = length / aspect_ratio
            origins["width"] = "derived"
        elif width is not None:
            length = aspect_ratio * width
            origins["length"] = "derived"

    estimate = None
    used_style = None
    if length is not None and width is not None:
        if area is not None and not agrees(area, length * width):
            raise ValueError("area conflicts with length and width.")
        if area is None:
            area = length * width
            origins["area"] = "derived"
    else:
        if area is None:
            if scaling_model is None or magnitude is None:
                raise ValueError("Supply area or magnitude and scaling_model.")
            used_style = scaling_style
            if used_style == "auto":
                used_style = (
                    "all" if scaling_model == "HanksBakun2002"
                    else fault_style_from_rake(rake)
                )
            model = create_scaling_relation(
                scaling_model, style=used_style, strict=scaling_strict,
            )
            estimate = model.magnitude_to_area(magnitude)
            area = estimate.at_epsilon(epsilon)
            origins["area"] = "scaling"
        if length is not None:
            width = area / length
            origins["width"] = "derived"
        elif width is not None:
            length = area / width
            origins["length"] = "derived"
        elif aspect_ratio is not None:
            target = sqrt(area / aspect_ratio)
            width = min(max(target, min_width), max_width)
            if width <= 0 or min_width > max_width:
                raise ValueError("No width satisfies the depth constraints.")
            length = area / width
            origins["width"] = (
                "derived" if agrees(width, target) else "constrained"
            )
            origins["length"] = origins["width"]
        else:
            raise ValueError("Supply length, width or aspect_ratio.")

    if not all(isfinite(v) and v > 0 for v in (area, length, width)):
        raise ValueError(
            "Resolved dimensions and area must be finite and positive."
        )
    if (width <= 0 or width < min_width - 1e-9
            or width > max_width + 1e-9):
        raise ValueError("Rupture width conflicts with depth constraints.")
    relief = width * slope
    if along_fraction is None:
        along_fraction = 0.5
        origins["along_fraction"] = "default"
    if ztor is None:
        if down_fraction is not None:
            ztor = depth - down_fraction * relief
            origins["ztor"] = "derived"
        else:
            target = depth - relief / 2
            ztor = max(upper, min(target, lower - relief))
            origins["ztor"] = (
                "derived" if agrees(ztor, target) else "constrained"
            )
            if slope == 0:
                down_fraction = 0.5
                origins["down_fraction"] = "default"
    if slope == 0 and down_fraction is None:
        down_fraction = 0.5
        origins["down_fraction"] = "default"

    # Remove harmless round-off at the surface.
    if -1e-9 <= ztor < 0:
        ztor = 0.0
    rupture = RectangularRupture.from_hypocentre(
        hypocentre, strike=strike, dip=dip, length=length, width=width,
        ztor=ztor, along_fraction=along_fraction,
        down_fraction=down_fraction,
    )
    if rupture.ztor < upper - 1e-9 or (
        rupture.ztor + relief > lower + 1e-9
    ):
        raise ValueError("Rupture extends outside the seismogenic layer.")

    resolved = dict(
        strike=strike, dip=dip, area=rupture.area,
        length=rupture.length, width=rupture.width, ztor=rupture.ztor,
        aspect_ratio=rupture.length / rupture.width,
        along_fraction=rupture.along_fraction,
        down_fraction=rupture.down_fraction,
    )
    origins.update(strike="explicit", dip="explicit")
    origins.setdefault("down_fraction", "derived")
    if aspect_ratio is None:
        origins["aspect_ratio"] = "derived"
    elif not agrees(resolved["aspect_ratio"], aspect_ratio):
        origins["aspect_ratio"] = "constrained"
    if rake is not None:
        resolved["rake"] = rake
        origins["rake"] = "explicit"
    return RuptureBuildResult(
        rupture, resolved, origins,
        scaling_model if estimate is not None else None,
        used_style, estimate,
    )
