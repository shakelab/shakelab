"""Empirical magnitude scaling and consistent rectangular dimensions.

Wells and Coppersmith (1994), BSSA 84(4), 974-1002, Table 2A.
doi:10.1785/BSSA0840040974

Magnitude is Mw. Length and width are in km, area in km**2.
Forward regressions predict log10(size); reverse regressions predict Mw.
The regressions are independent fits, not algebraic inverses.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log10, sqrt
from numbers import Real
import warnings


def _scalar(value: float, name: str) -> float:
    """Validate a finite scalar."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number.")
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite.")
    return value


def _positive(value: float, name: str) -> float:
    value = _scalar(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be strictly positive.")
    return value


def fault_style_from_rake(rake: float | None) -> str:
    """Classify rake with +/-45 and +/-135 degree boundaries.

    Boundary values are assigned to strike-slip. None selects all styles.
    This is a classification convention, not a regression in WC1994.
    """
    if rake is None:
        return "all"
    rake = _scalar(rake, "rake")
    if not -180 <= rake <= 180:
        raise ValueError("rake must be in [-180, 180] degrees.")
    if 45 < rake < 135:
        return "reverse"
    if -135 < rake < -45:
        return "normal"
    return "strike-slip"


@dataclass(frozen=True)
class ScalingEstimate:
    """Prediction and residual standard deviation in its regression space."""

    median: float
    sigma: float
    sigma_space: str
    extrapolated: bool

    def at_epsilon(self, epsilon: float = 0.0) -> float:
        """Evaluate a signed number of residual standard deviations."""
        epsilon = _scalar(epsilon, "epsilon")
        if self.sigma_space == "log10":
            return self.median * 10.0 ** (epsilon * self.sigma)
        return self.median + epsilon * self.sigma


@dataclass(frozen=True)
class RuptureDimensions:
    """Consistent rectangular dimensions in km, with area in km**2."""

    length: float
    width: float

    def __post_init__(self) -> None:
        for name in ("length", "width"):
            object.__setattr__(
                self, name, _positive(getattr(self, name), name),
            )

    @property
    def area(self) -> float:
        return self.length * self.width

    @property
    def aspect_ratio(self) -> float:
        return self.length / self.width

    @classmethod
    def from_area(
        cls,
        area: float,
        *,
        aspect_ratio: float | None = None,
        length: float | None = None,
        width: float | None = None,
        max_width: float | None = None,
    ) -> RuptureDimensions:
        """Preserve area using exactly one explicit shape constraint.

        When an aspect ratio is used, a maximum width caps the derived
        width and extends length to preserve area. Explicit dimensions
        are never changed to enforce the cap.
        """
        area = _positive(area, "area")
        if sum(x is not None for x in
               (aspect_ratio, length, width)) != 1:
            raise ValueError(
                "Specify exactly one of aspect_ratio, length or width."
            )
        if aspect_ratio is not None:
            ratio = _positive(aspect_ratio, "aspect_ratio")
            width = sqrt(area / ratio)
        elif length is not None:
            width = area / _positive(length, "length")
        else:
            width = _positive(width, "width")
        if max_width is not None:
            limit = _positive(max_width, "max_width")
            if width > limit:
                if aspect_ratio is None:
                    raise ValueError(
                        "Specified dimension conflicts with max_width."
                    )
                width = limit
        return cls(area / width, width)


class WellsCoppersmith1994:
    """Independent regressions from WC1994 Table 2A.

    Select all, strike-slip, reverse or normal. Predictions outside the
    published calibration range warn by default, or raise with strict=True.
    Residual sigma is in base-10 logarithmic units in forward regressions
    and magnitude units in reverse regressions.
    """

    # Forward rows: intercept, slope, residual sigma, Mw min, Mw max.
    # Reverse rows: intercept, slope, residual sigma, size min, size max;
    # width uses Mw bounds because Table 2A repeats the length ranges.
    _FORWARD = {
        "length": {
            "all": (-2.44, .59, .16, 4.8, 8.1),
            "strike-slip": (-2.57, .62, .15, 4.8, 8.1),
            "reverse": (-2.42, .58, .16, 4.8, 7.6),
            "normal": (-1.88, .50, .17, 5.2, 7.3),
        },
        "width": {
            "all": (-1.01, .32, .15, 4.8, 8.1),
            "strike-slip": (-.76, .27, .14, 4.8, 8.1),
            "reverse": (-1.61, .41, .15, 4.8, 7.6),
            "normal": (-1.14, .35, .12, 5.2, 7.3),
        },
        "area": {
            "all": (-3.49, .91, .24, 4.8, 7.9),
            "strike-slip": (-3.42, .90, .22, 4.8, 7.9),
            "reverse": (-3.99, .98, .26, 4.8, 7.6),
            "normal": (-2.87, .82, .22, 5.2, 7.3),
        },
    }
    _REVERSE = {
        "length": {
            "all": (4.38, 1.49, .26, 1.1, 350),
            "strike-slip": (4.33, 1.49, .24, 1.5, 350),
            "reverse": (4.49, 1.49, .26, 1.1, 80),
            "normal": (4.34, 1.54, .31, 3.8, 63),
        },
        "width": {
            "all": (4.06, 2.25, .41, 4.8, 8.1),
            "strike-slip": (3.80, 2.59, .45, 4.8, 8.1),
            "reverse": (4.37, 1.95, .32, 4.8, 7.6),
            "normal": (4.04, 2.11, .31, 5.2, 7.3),
        },
        "area": {
            "all": (4.07, .98, .24, 2.2, 5184),
            "strike-slip": (3.98, 1.02, .23, 3, 5184),
            "reverse": (4.33, .90, .25, 2.2, 2400),
            "normal": (3.93, 1.02, .25, 19, 900),
        },
    }

    def __init__(self, style: str = "all", *, strict: bool = False):
        if style not in self._FORWARD["area"]:
            raise ValueError(f"Unknown fault style: {style!r}.")
        if not isinstance(strict, bool):
            raise TypeError("strict must be boolean.")
        self.style = style
        self.strict = strict

    def _range_check(self, value, lower, upper, name):
        outside = not lower <= value <= upper
        if outside:
            message = (
                f"WC1994 {self.style} {name}: outside calibration "
                f"range [{lower}, {upper}]."
            )
            if self.strict:
                raise ValueError(message)
            warnings.warn(message, UserWarning, stacklevel=3)
        return outside

    def estimate(self, magnitude: float, quantity: str) -> ScalingEstimate:
        """Predict median length, width or area from Mw."""
        magnitude = _scalar(magnitude, "magnitude")
        if quantity not in self._FORWARD:
            raise ValueError(f"Unknown scaling quantity: {quantity!r}.")
        a, b, sigma, lower, upper = self._FORWARD[quantity][self.style]
        outside = self._range_check(magnitude, lower, upper, quantity)
        return ScalingEstimate(
            10.0 ** (a + b * magnitude), sigma, "log10", outside,
        )

    def magnitude_to_area(self, magnitude: float) -> ScalingEstimate:
        return self.estimate(magnitude, "area")

    def magnitude_to_length(self, magnitude: float) -> ScalingEstimate:
        return self.estimate(magnitude, "length")

    def magnitude_to_width(self, magnitude: float) -> ScalingEstimate:
        return self.estimate(magnitude, "width")

    def estimate_magnitude(
        self, size: float, quantity: str,
    ) -> ScalingEstimate:
        """Use the independent regression of Mw on a rupture parameter."""
        size = _positive(size, quantity)
        if quantity not in self._REVERSE:
            raise ValueError(f"Unknown scaling quantity: {quantity!r}.")
        a, b, sigma, lower, upper = self._REVERSE[quantity][self.style]
        magnitude = a + b * log10(size)
        # WC1994's printed RW range repeats RLD's dimension range.
        # Do not present it as an observed width calibration interval.
        diagnostic = magnitude if quantity == "width" else size
        outside = self._range_check(diagnostic, lower, upper, quantity)
        return ScalingEstimate(magnitude, sigma, "Mw", outside)

    def area_to_magnitude(self, area: float) -> ScalingEstimate:
        return self.estimate_magnitude(area, "area")

    def length_to_magnitude(self, length: float) -> ScalingEstimate:
        return self.estimate_magnitude(length, "length")

    def width_to_magnitude(self, width: float) -> ScalingEstimate:
        return self.estimate_magnitude(width, "width")
