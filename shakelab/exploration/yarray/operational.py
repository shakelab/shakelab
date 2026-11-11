# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Operational assessment of an array (survey planning).

Combines the geometric response (:mod:`.response`) and the Monte Carlo
performance of an estimator (:mod:`.performance`) into a
frequency-velocity-azimuth envelope: for every frequency ``f``, phase
velocity ``c`` and propagation azimuth, the experiment of
:meth:`ArrayResponse.performance
<shakelab.exploration.yarray.response.ArrayResponse.performance>` is
run with the wavenumber ``k = 2 pi f / c`` and the configuration of the
user (:meth:`ArrayResponse.assess
<shakelab.exploration.yarray.response.ArrayResponse.assess>`). The
envelope holds the probabilities of LOCAL, ALIAS, OTHER and BOUNDARY
outcomes, the local wavevector RMSE, the exact relative velocity RMSE,
the quantization error of the search grid and a usability mask.

There is no opaque quality score: a point is usable when it meets
explicit criteria (the ``min_*`` and ``max_*`` parameters of
``assess``), and otherwise it
gets the first limiting factor found (:class:`LimitingFactor`). The
recommendations (:class:`Recommendation`) explain the limits and
suggest processing changes (search range, grid sampling, SNR,
realizations); the geometry itself is not optimized.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

import numpy as np

from shakelab.exploration.yarray.config import PerformanceConfig
from shakelab.exploration.yarray.performance import (
    PerformanceResult,
    PlaneWaveSource,
    _check,
    _context,
    _evaluate,
    _search,
)
from shakelab.exploration.yarray.response import ArrayResponse, _readonly

__all__ = [
    "LimitingFactor",
    "OperationalAssessment",
    "Recommendation",
    "RecommendationSeverity",
    "RecommendationType",
]

# Criteria of a usable point and their defaults (see ArrayResponse.assess).
_CRITERIA = ("min_local", "max_alias", "max_other", "max_boundary",
             "max_velocity_rmse", "max_grid_error")


def _axis(value: Any, name: str) -> np.ndarray:
    """Return an increasing axis of positive values."""
    axis = np.asarray(value, dtype=float)

    if (
        axis.ndim != 1 or axis.size == 0 or not np.all(np.isfinite(axis))
        or np.any(axis <= 0.0) or np.any(np.diff(axis) <= 0.0)
    ):
        raise ValueError(f"{name} must be increasing positive values")

    return axis


class LimitingFactor(Enum):
    """Main limitation of an envelope point."""

    NONE = "none"
    RESOLUTION = "resolution"
    ALIASING = "aliasing"
    LOW_SNR = "low_snr"
    SEARCH_BOUNDARY = "search_boundary"
    GRID_QUANTIZATION = "grid_quantization"
    OTHER_FAILURE = "other_failure"


class RecommendationType(Enum):
    """Category of a recommendation."""

    ARRAY_ANISOTROPY = "array_anisotropy"
    SEARCH_RANGE = "search_range"
    RADIAL_SAMPLING = "radial_sampling"
    REALIZATIONS = "realizations"
    SNR = "snr"


class RecommendationSeverity(Enum):
    """Severity of a recommendation."""

    INFO = "info"
    WARNING = "warning"
    CRITICAL = "critical"


def _azimuths(value: Any) -> np.ndarray:
    """Return distinct azimuths in ``[0, 360)`` (12 uniform if None)."""
    if value is None:
        return np.linspace(0.0, 360.0, 12, endpoint=False)

    azimuths = np.mod(np.asarray(value, dtype=float), 360.0)

    if (
        azimuths.ndim != 1 or azimuths.size == 0
        or not np.all(np.isfinite(azimuths))
        or np.unique(azimuths).size != azimuths.size
    ):
        raise ValueError("azimuths must be distinct finite values")

    return azimuths


def _criteria(values: dict) -> dict:
    """Return validated criteria."""
    criteria = {}

    for name in _CRITERIA:
        value = float(values[name])

        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")

        if name in ("min_local", "max_alias", "max_other",
                    "max_boundary") and value > 1.0:
            raise ValueError(f"{name} cannot exceed one")

        criteria[name] = value

    return criteria


@dataclass(frozen=True)
class Recommendation:
    """
    Operational recommendation.

    Attributes
    ----------
    category : RecommendationType
        Category.
    severity : RecommendationSeverity
        Severity.
    message : str
        Explanation.
    value : float, optional
        Quantitative target, if any.
    unit : str, optional
        Unit of ``value``.
    """

    category: RecommendationType
    severity: RecommendationSeverity
    message: str
    value: Optional[float] = None
    unit: Optional[str] = None


@dataclass(frozen=True, eq=False)
class OperationalAssessment:
    """
    Frequency-velocity-azimuth envelope of an array.

    The envelope arrays have shape ``(n_frequencies, n_velocities,
    n_azimuths)``; the ``conservative_*`` properties reduce them over
    the azimuths (worst case).

    Attributes
    ----------
    response : ArrayResponse
        Geometric response.
    config : PerformanceConfig
        Configuration of the experiments (a copy).
    frequencies, velocities, azimuths : numpy.ndarray
        Axes of the envelope: frequencies in Hz, phase velocities and
        propagation azimuths in degrees clockwise from north.
    criteria : dict
        Criteria of a usable point (``min_local``, ``max_alias``,
        ``max_other``, ``max_boundary``, ``max_velocity_rmse``,
        ``max_grid_error``; see :meth:`ArrayResponse.assess
        <shakelab.exploration.yarray.response.ArrayResponse.assess>`).
    conservative_azimuth : bool
        Rule of :attr:`conservative_usable` (all the azimuths, or at
        least half).
    resolution_radius : float
        Resolution radius of the response at the classification level.
    local_probability, alias_probability : numpy.ndarray
        Outcome probabilities.
    other_probability, boundary_probability : numpy.ndarray
        Outcome probabilities.
    local_wavevector_rmse : numpy.ndarray
        LOCAL wavevector RMSE (NaN without LOCAL estimates).
    relative_velocity_rmse : numpy.ndarray
        LOCAL velocity RMSE relative to the true velocity.
    grid_quantization_error : numpy.ndarray
        Distance to the nearest grid wavevector (NaN outside the grid).
    usable : numpy.ndarray
        True where all the criteria are met.
    limiting_factor : numpy.ndarray
        :class:`LimitingFactor` of each point.
    recommendations : tuple of Recommendation
        Recommendations.
    """

    response: ArrayResponse
    config: PerformanceConfig
    frequencies: np.ndarray
    velocities: np.ndarray
    azimuths: np.ndarray
    criteria: dict
    conservative_azimuth: bool
    resolution_radius: float
    local_probability: np.ndarray
    alias_probability: np.ndarray
    other_probability: np.ndarray
    boundary_probability: np.ndarray
    local_wavevector_rmse: np.ndarray
    relative_velocity_rmse: np.ndarray
    grid_quantization_error: np.ndarray
    usable: np.ndarray
    limiting_factor: np.ndarray
    recommendations: tuple = ()

    @property
    def conservative_usable(self) -> np.ndarray:
        """Usability of the frequency-velocity points."""
        if self.conservative_azimuth:
            return np.all(self.usable, axis=2)

        return np.mean(self.usable, axis=2) >= 0.5

    @property
    def conservative_local_probability(self) -> np.ndarray:
        """Smallest LOCAL probability over the azimuths."""
        return np.min(self.local_probability, axis=2)

    @property
    def conservative_alias_probability(self) -> np.ndarray:
        """Largest ALIAS probability over the azimuths."""
        return np.max(self.alias_probability, axis=2)

    @property
    def conservative_relative_velocity_rmse(self) -> np.ndarray:
        """Largest relative velocity RMSE over the azimuths."""
        values = self.relative_velocity_rmse
        result = np.full(values.shape[:2], np.nan)
        valid = np.any(np.isfinite(values), axis=2)
        result[valid] = np.nanmax(values[valid], axis=1)

        return result

    def usable_velocity_range(self, frequency_index: int
                              ) -> Optional[tuple[float, float]]:
        """
        Return the extent of the usable velocities at a frequency.

        The interval spans all the usable samples; velocities between
        separate usable islands are not necessarily usable.

        Parameters
        ----------
        frequency_index : int
            Index of the frequency.

        Returns
        -------
        tuple of float or None
            Smallest and largest usable velocity, None if none.
        """
        indices = np.flatnonzero(self.conservative_usable[frequency_index])

        if indices.size == 0:
            return None

        return (float(self.velocities[indices[0]]),
                float(self.velocities[indices[-1]]))

    def info(self) -> None:
        """Print a summary of the assessment."""
        geometry = self.response.geometry
        config = self.config
        usable = self.conservative_usable
        print(
            f"OperationalAssessment: {geometry.n_stations} stations "
            f"(spacing {geometry.minimum_spacing:.6g} to "
            f"{geometry.maximum_aperture:.6g}, anisotropy "
            f"{geometry.anisotropy_ratio:.4g}), "
            f"{config.method.value}, SNR {config.snr_db:g} dB, "
            f"{config.realizations} realizations"
        )
        print(
            f"  envelope {self.frequencies.size} frequencies x "
            f"{self.velocities.size} velocities x {self.azimuths.size} "
            f"azimuths, {np.count_nonzero(usable)} of {usable.size} "
            "points usable"
        )

        for index, frequency in enumerate(self.frequencies):
            interval = self.usable_velocity_range(index)
            text = (
                "none" if interval is None
                else f"{interval[0]:.6g} to {interval[1]:.6g}"
            )
            print(f"  {frequency:g} Hz: {text}")

        for item in self.recommendations:
            print(f"  [{item.severity.value}] {item.message}")


def _limiting_factor(result: PerformanceResult, frequency: float,
                     criteria: dict, radius: float) -> LimitingFactor:
    """Return the first limiting factor of an experiment."""
    outcome = result.outcome

    if outcome.boundary_probability > criteria["max_boundary"]:
        return LimitingFactor.SEARCH_BOUNDARY

    fraction = (
        result.grid_quantization_error / radius if radius > 0.0
        else float("inf")
    )

    if fraction > criteria["max_grid_error"]:
        return LimitingFactor.GRID_QUANTIZATION

    if outcome.alias_probability > criteria["max_alias"]:
        return LimitingFactor.ALIASING

    if outcome.other_probability > criteria["max_other"]:
        return LimitingFactor.OTHER_FAILURE

    if outcome.local_probability < criteria["min_local"]:
        return LimitingFactor.LOW_SNR

    velocity = result.velocity_statistics(frequency)

    if velocity is None:
        return LimitingFactor.OTHER_FAILURE

    if velocity.relative_rmse > criteria["max_velocity_rmse"]:
        return LimitingFactor.RESOLUTION

    return LimitingFactor.NONE


def _recommendations(assessment: OperationalAssessment) -> tuple:
    """Return the recommendations of an assessment."""
    items = []
    criteria = assessment.criteria
    anisotropy = assessment.response.geometry.anisotropy_ratio

    def add(category, severity, message, value=None):
        items.append(Recommendation(category, severity, message, value))

    if anisotropy > 2.0:
        add(RecommendationType.ARRAY_ANISOTROPY,
            RecommendationSeverity.WARNING,
            "The geometry is strongly anisotropic: the resolution varies "
            "with the propagation azimuth.", float(anisotropy))

    if not np.any(assessment.conservative_usable):
        add(RecommendationType.SEARCH_RANGE,
            RecommendationSeverity.CRITICAL,
            "No frequency-velocity point meets all the criteria: review "
            "search range, SNR, averaging or geometry.")

    if np.any(assessment.conservative_alias_probability
              > criteria["max_alias"]):
        add(RecommendationType.SEARCH_RANGE,
            RecommendationSeverity.WARNING,
            "Part of the domain is limited by geometric aliasing: restrict "
            "the high wavenumbers searched where physically acceptable.")

    fraction = assessment.grid_quantization_error / max(
        assessment.resolution_radius, np.finfo(float).eps
    )

    if np.any(fraction[np.isfinite(fraction)]
              > criteria["max_grid_error"]):
        add(RecommendationType.RADIAL_SAMPLING,
            RecommendationSeverity.WARNING,
            "The quantization of the search grid is not negligible with "
            "respect to the resolution of the array: sample radius and/or "
            "azimuth more finely.")

    if np.max(assessment.other_probability) > \
            criteria["max_other"]:
        add(RecommendationType.SNR, RecommendationSeverity.INFO,
            "Noise-driven OTHER outcomes occur in part of the domain: test "
            "higher SNR or more averaging with "
            "ArrayResponse.performance_sweep().")

    if np.min(assessment.local_probability) < \
            criteria["min_local"]:
        add(RecommendationType.REALIZATIONS, RecommendationSeverity.INFO,
            "The LOCAL probability is below the criterion in part of the "
            "domain: test more Fourier realizations before changing the "
            "geometry.")

    if not items:
        add(RecommendationType.SEARCH_RANGE, RecommendationSeverity.INFO,
            "The domain meets the criteria: no processing change is "
            "indicated.")

    return tuple(items)


def _assess(response: ArrayResponse, config: PerformanceConfig,
            frequencies: Any, velocities: Any, azimuths: Any,
            criteria: dict, conservative_azimuth: bool
            ) -> OperationalAssessment:
    """Return the envelope (see :meth:`ArrayResponse.assess`)."""
    if not isinstance(response, ArrayResponse):
        raise TypeError("response must be an ArrayResponse")

    if not isinstance(config, PerformanceConfig):
        raise TypeError("config must be a PerformanceConfig")

    if not isinstance(conservative_azimuth, bool):
        raise TypeError("conservative_azimuth must be a bool")

    config.check_complete()
    config = config.copy()
    criteria = _criteria(criteria)
    frequencies = _readonly(_axis(frequencies, "frequencies"))
    velocities = _readonly(_axis(velocities, "velocities"))
    azimuths = _readonly(_azimuths(azimuths))
    grid = config.grid
    shape = (frequencies.size, velocities.size, azimuths.size)
    values = {
        name: np.full(shape, np.nan) for name in (
            "local_probability", "alias_probability", "other_probability",
            "boundary_probability", "local_wavevector_rmse",
            "relative_velocity_rmse", "grid_quantization_error",
        )
    }
    usable = np.zeros(shape, dtype=bool)
    factors = np.empty(shape, dtype=object)
    context = _context(response, config)
    search = _search(response, grid)
    radius = context.analysis.limits.resolution_radius
    base = config.seed or 0
    point = 0

    for index in np.ndindex(shape):
        i, j, m = index
        frequency = float(frequencies[i])
        wavenumber = 2.0 * np.pi * frequency / velocities[j]
        seed = base + point
        point += 1

        if not grid.minimum <= wavenumber <= grid.maximum:
            for name in ("local_probability", "alias_probability",
                         "other_probability"):
                values[name][index] = 0.0

            values["boundary_probability"][index] = 1.0
            factors[index] = LimitingFactor.SEARCH_BOUNDARY
            continue

        source = PlaneWaveSource(float(wavenumber), float(azimuths[m]),
                                 frequency=frequency)
        experiment = config.copy(seed=seed)
        _check(response, source, experiment)
        result = _evaluate(response, source, experiment, context, search)
        outcome = result.outcome

        values["local_probability"][index] = outcome.local_probability
        values["alias_probability"][index] = outcome.alias_probability
        values["other_probability"][index] = outcome.other_probability
        values["boundary_probability"][index] = \
            outcome.boundary_probability
        values["grid_quantization_error"][index] = \
            result.grid_quantization_error

        if result.local_statistics is not None:
            values["local_wavevector_rmse"][index] = \
                result.local_statistics.rmse

        velocity = result.velocity_statistics(frequency)

        if velocity is not None:
            values["relative_velocity_rmse"][index] = velocity.relative_rmse

        factors[index] = _limiting_factor(result, frequency, criteria,
                                          radius)
        usable[index] = factors[index] is LimitingFactor.NONE

    assessment = OperationalAssessment(
        response=response,
        config=config,
        frequencies=frequencies,
        velocities=velocities,
        azimuths=azimuths,
        criteria=criteria,
        conservative_azimuth=conservative_azimuth,
        resolution_radius=float(radius),
        usable=_readonly(usable, bool),
        limiting_factor=_readonly(factors, object),
        **{name: _readonly(value) for name, value in values.items()},
    )
    object.__setattr__(assessment, "recommendations",
                       _recommendations(assessment))

    return assessment
