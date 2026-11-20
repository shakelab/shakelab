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
    sigma: float | None
    sigma_space: str
    extrapolated: bool | None

    def at_epsilon(self, epsilon: float = 0.0) -> float:
        """Evaluate a signed number of residual standard deviations."""
        epsilon = _scalar(epsilon, "epsilon")
        if epsilon == 0.0:
            return self.median
        if self.sigma is None:
            raise ValueError("Residual dispersion is not available.")
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


class ThingbaijamEtAl2017(WellsCoppersmith1994):
    """Table 1 of Thingbaijam, Mai and Goda (2017).

    DOI: 10.1785/0120170017. Forward predictions use the published
    log10 regressions. Reverse predictions invert the orthogonal
    regression line and transform the residual by the slope; they
    are not separately fitted conditional magnitude regressions.
    No pooled all-style regression is available.
    """

    _FORWARD = {
        "length": {
            "reverse": (-2.693, .614, .083, 5.59, 7.69),
            "normal": (-1.722, .485, .128, 5.86, 8.39),
            "strike-slip": (-2.943, .681, .151, 5.38, 8.70),
            "interface": (-2.412, .583, .107, 6.68, 9.19),
        },
        "width": {
            "reverse": (-1.669, .435, .087, 5.59, 7.69),
            "normal": (-.829, .323, .128, 5.86, 8.39),
            "strike-slip": (-.543, .261, .105, 5.38, 8.70),
            "interface": (-.880, .366, .099, 6.68, 9.19),
        },
        "area": {
            "reverse": (-4.362, 1.049, .121, 5.59, 7.69),
            "normal": (-2.551, .808, .181, 5.86, 8.39),
            "strike-slip": (-3.486, .942, .184, 5.38, 8.70),
            "interface": (-3.292, .949, .150, 6.68, 9.19),
        },
    }
    _SIZE_RANGES = {
        "length": {
            "reverse": (4.9, 108), "normal": (9, 262.5),
            "strike-slip": (6, 580), "interface": (29.2, 1420),
        },
        "width": {
            "reverse": (4.8, 45), "normal": (6, 112.5),
            "strike-slip": (6.5, 50), "interface": (29.2, 260),
        },
        "area": {
            "reverse": (23.5, 4860), "normal": (54, 29531.3),
            "strike-slip": (39, 29000), "interface": (852.6, 318080),
        },
    }

    def __init__(self, style: str, *, strict: bool = False):
        super().__init__(style, strict=strict)

    def _range_check(self, value, lower, upper, name):
        outside = not lower <= value <= upper
        if outside:
            message = (
                f"TMG2017 {self.style} {name}: outside calibration "
                f"range [{lower}, {upper}]."
            )
            if self.strict:
                raise ValueError(message)
            warnings.warn(message, UserWarning, stacklevel=3)
        return outside

    def estimate_magnitude(
        self, size: float, quantity: str,
    ) -> ScalingEstimate:
        """Invert the fitted line and transform log10 residual units."""
        size = _positive(size, quantity)
        if quantity not in self._FORWARD:
            raise ValueError(f"Unknown scaling quantity: {quantity!r}.")
        a, b, sigma, _, _ = self._FORWARD[quantity][self.style]
        lower, upper = self._SIZE_RANGES[quantity][self.style]
        outside = self._range_check(size, lower, upper, quantity)
        return ScalingEstimate(
            (log10(size) - a) / b, sigma / abs(b), "Mw", outside,
        )


class _AreaScaling:
    """Common interface for deterministic magnitude-area laws.

    A None sigma or extrapolated flag indicates unavailable information,
    not zero dispersion or confirmed validity.
    """

    def estimate(self, magnitude: float, quantity: str) -> ScalingEstimate:
        if quantity != "area":
            raise ValueError(f"Unsupported scaling quantity: {quantity!r}.")
        return self.magnitude_to_area(magnitude)

    def estimate_magnitude(
        self, size: float, quantity: str,
    ) -> ScalingEstimate:
        if quantity != "area":
            raise ValueError(f"Unsupported scaling quantity: {quantity!r}.")
        return self.area_to_magnitude(size)


class _Leonard2014(_AreaScaling):
    """Magnitude-area subset of Leonard (2014), BSSA 104, 2953-2965.

    DOI: 10.1785/0120140087. Median constants follow the OpenQuake
    magnitude-area convention. No residual dispersion or quantitative
    calibration limits are assigned here.
    """

    _CONSTANTS: dict[str, float] = {}

    def __init__(self, style: str = "all", *, strict: bool = False):
        if style not in self._CONSTANTS:
            raise ValueError(f"Unknown fault style: {style!r}.")
        if not isinstance(strict, bool):
            raise TypeError("strict must be boolean.")
        if strict:
            raise ValueError(
                "Strict calibration checks are unavailable for this model."
            )
        self.style = style
        self.strict = strict

    def magnitude_to_area(self, magnitude: float) -> ScalingEstimate:
        magnitude = _scalar(magnitude, "magnitude")
        area = 10.0 ** (magnitude - self._CONSTANTS[self.style])
        return ScalingEstimate(area, None, "log10", None)

    def area_to_magnitude(self, area: float) -> ScalingEstimate:
        area = _positive(area, "area")
        magnitude = log10(area) + self._CONSTANTS[self.style]
        return ScalingEstimate(magnitude, None, "Mw", None)


class Leonard2014_Interplate(_Leonard2014):
    """Tectonically active crustal faults (not a subduction-interface law)."""

    _CONSTANTS = {
        "all": 3.995, "strike-slip": 3.99,
        "normal": 4.00, "reverse": 4.00,
    }


class Leonard2014_StableContinental(_Leonard2014):
    """Stable continental regions."""

    _CONSTANTS = {
        "all": 4.185, "strike-slip": 4.18,
        "normal": 4.19, "reverse": 4.19,
    }


class HanksBakun2002(_AreaScaling):
    """Bilinear continental-earthquake magnitude-area relation.

    BSSA 92, 1841-1846, doi:10.1785/0120010148.
    Published rounded intercepts give a tiny discontinuity at 537 km**2.
    For inversion use the lower-branch magnitude at 537 km**2 as the
    crossover, preserving the published coefficients on both branches.
    No residual dispersion or quantitative calibration bounds are
    assigned. Intercept standard errors are not prediction sigmas.
    """

    crossover_area = 537.0
    crossover_magnitude = 3.98 + log10(crossover_area)

    def __init__(self, style: str = "all", *, strict: bool = False):
        if style != "all":
            raise ValueError("HanksBakun2002 supports only style='all'.")
        if not isinstance(strict, bool):
            raise TypeError("strict must be boolean.")
        if strict:
            raise ValueError(
                "Strict calibration checks are unavailable for this model."
            )
        self.style = style
        self.strict = strict

    def magnitude_to_area(self, magnitude: float) -> ScalingEstimate:
        magnitude = _scalar(magnitude, "magnitude")
        if magnitude <= self.crossover_magnitude:
            area = 10.0 ** (magnitude - 3.98)
        else:
            area = 10.0 ** (0.75 * (magnitude - 3.07))
        return ScalingEstimate(area, None, "log10", None)

    def area_to_magnitude(self, area: float) -> ScalingEstimate:
        area = _positive(area, "area")
        if area <= self.crossover_area:
            magnitude = 3.98 + log10(area)
        else:
            magnitude = 3.07 + (4.0 / 3.0) * log10(area)
        return ScalingEstimate(magnitude, None, "Mw", None)


_SCALING_RELATIONS = {
    "Leonard2014_Interplate": Leonard2014_Interplate,
    "Leonard2014_StableContinental": Leonard2014_StableContinental,
    "HanksBakun2002": HanksBakun2002,
    "WellsCoppersmith1994": WellsCoppersmith1994,
    "ThingbaijamEtAl2017": ThingbaijamEtAl2017,
}


def available_scaling_relations() -> tuple[str, ...]:
    """Return identifiers accepted by create_scaling_relation."""
    return tuple(sorted(_SCALING_RELATIONS))


def create_scaling_relation(name: str, **kwargs):
    """Instantiate a named scaling relation without external dependencies."""
    if not isinstance(name, str) or name not in _SCALING_RELATIONS:
        raise ValueError(
            f"Unknown scaling relation {name!r}; available: "
            + ", ".join(available_scaling_relations())
        )
    return _SCALING_RELATIONS[name](**kwargs)
