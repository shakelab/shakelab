"""Planar rectangular rupture geometry in a local WGS84 ENU frame.

Public lengths, depths and distances are in kilometres; WgsPoint elevation
and Cartesian coordinates returned by geodeticN are in metres. The finite
surface is a plane in one local tangent frame, not a curved Earth surface.
Distances are exact for that Cartesian representation. For large ruptures
or distant sites, the tangent-plane approximation needs reassessment.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, radians, sin
from numbers import Real

import numpy as np

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
