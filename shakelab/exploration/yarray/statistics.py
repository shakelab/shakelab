# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
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
# You should have received a copy of the GNU Affero General Public
# License with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Statistical post-processing of YArray spectral maxima.

This module operates on the compact maxima contained in a
:class:`BeamformingResult`. It intentionally does not depend on
Matplotlib or other graphical backends.

The first implementation provides conditional empirical densities of
spectral maxima as a function of frequency. Statistics can be evaluated
in radial angular wavenumber or logarithmic slowness coordinates.

The density at each analysis frequency is estimated as a weighted
Gaussian kernel mixture. By default, kernel bandwidth is tied to the
effective radial resolution of the beamforming search grid rather than
to the global spread of maxima at that frequency. This prevents
multimodal samples from artificially broadening every individual mode.

A robust adaptive bandwidth and an explicit fixed bandwidth remain
available as alternatives. If future YArray results provide positional
standard errors, these can be folded into each kernel by convolution
in quadrature.

The resulting quantity is an empirical conditional density, not a
Bayesian posterior probability for the true dispersion curve.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from .beamforming import (
    BeamformingResult,
    GridQuantity,
    SpectralMaximum,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


__all__ = [
    "BandwidthMode",
    "DensityCoordinate",
    "DensityPeaks",
    "DensitySummary",
    "DensityWeighting",
    "MaximumDensity",
    "density_peaks",
    "density_summary",
    "maximum_density",
]


class BandwidthMode(Enum):
    """Strategy used to determine Gaussian density bandwidth."""

    RESOLUTION = "resolution"
    ADAPTIVE = "adaptive"
    FIXED = "fixed"


class DensityCoordinate(Enum):
    """Coordinate used for maximum-density estimation."""

    WAVENUMBER = "wavenumber"
    LOG_SLOWNESS = "log_slowness"
    LOG_VELOCITY = "log_velocity"


class DensityWeighting(Enum):
    """Weight assigned to individual spectral maxima."""

    UNIFORM = "uniform"
    POWER = "power"
    RELATIVE_POWER = "relative_power"


@dataclass(frozen=True)
class DensityPeaks:
    """
    Local maxima of the conditional density at each frequency.

    Parameters
    ----------
    frequencies
        Frequencies represented by the density estimate.
    values
        Density-coordinate values of selected local maxima, with shape
        ``(n_frequencies, n_peaks)``. Missing peaks are represented by
        ``NaN``.
    relative_density
        Peak heights normalized by the global density maximum at the
        same frequency. Values therefore lie in ``[0, 1]``.
    density
        Absolute conditional-density values at the selected maxima.
    prominence
        Peak prominence in absolute conditional-density units.
    count
        Number of retained local maxima at each frequency.

    Notes
    -----
    Peaks are detected independently at each frequency and sorted by
    decreasing density amplitude. Their column index is therefore an
    amplitude order, not a tracked modal identity across frequency.
    """

    frequencies: FloatArray
    values: FloatArray
    relative_density: FloatArray
    density: FloatArray
    prominence: FloatArray
    count: IntArray


@dataclass(frozen=True)
class DensitySummary:
    """
    Frequency-dependent summary of a maximum-density estimate.

    Parameters
    ----------
    frequencies
        Frequencies represented by the summary.
    mode
        Coordinate of the global density maximum at each frequency.
    median
        Coordinate dividing the conditional density into equal halves.
    lower
        Lower conditional quantile.
    upper
        Upper conditional quantile.
    lower_probability
        Cumulative probability associated with ``lower``.
    upper_probability
        Cumulative probability associated with ``upper``.

    Notes
    -----
    The median and quantile interval summarize the complete conditional
    density. In a strongly multimodal distribution they may fall
    between distinct modes. The ``mode`` is therefore generally more
    appropriate for visual ridge diagnostics.
    """

    frequencies: FloatArray
    mode: FloatArray
    median: FloatArray
    lower: FloatArray
    upper: FloatArray
    lower_probability: float
    upper_probability: float


@dataclass(frozen=True)
class MaximumDensity:
    """
    Conditional empirical density of spectral maxima.

    Parameters
    ----------
    frequencies
        Sorted frequencies represented in the density.
    values
        Common coordinate grid used for every frequency.
    density
        Conditional density with shape
        ``(n_frequencies, n_values)``. Each non-empty row integrates
        to one over ``values``.
    counts
        Number of selected maxima contributing at each frequency.
    bandwidths
        Gaussian KDE bandwidth used at each frequency before optional
        positional uncertainty is folded into each kernel.
    bandwidth_mode
        Strategy used to estimate the bandwidths.
    bandwidth_factor
        Multiplicative factor applied to resolution-based bandwidths.
    component
        Processed component represented by the density.
    coordinate
        Statistical coordinate: angular wavenumber, logarithmic
        slowness, or logarithmic velocity.
    weighting
        Maximum-weighting strategy.
    ranks
        Optional one-based maximum ranks retained in the estimate.
    """

    frequencies: FloatArray
    values: FloatArray
    density: FloatArray
    counts: IntArray
    bandwidths: FloatArray
    bandwidth_mode: BandwidthMode
    bandwidth_factor: float
    component: str
    coordinate: DensityCoordinate
    weighting: DensityWeighting
    ranks: tuple[int, ...] | None

    def __post_init__(self) -> None:
        if self.frequencies.ndim != 1:
            raise ValueError(
                "frequencies must be one-dimensional."
            )

        if self.values.ndim != 1:
            raise ValueError(
                "values must be one-dimensional."
            )

        expected_shape = (
            self.frequencies.size,
            self.values.size,
        )

        if self.density.shape != expected_shape:
            raise ValueError(
                "density shape is inconsistent with its axes."
            )

        if self.counts.shape != self.frequencies.shape:
            raise ValueError(
                "counts shape is inconsistent with frequencies."
            )

        if self.bandwidths.shape != self.frequencies.shape:
            raise ValueError(
                "bandwidths shape is inconsistent with frequencies."
            )

    @property
    def velocity(self) -> FloatArray:
        """
        Return the velocity grid implied by log-slowness coordinates.

        Raises
        ------
        ValueError
            If the density coordinate is not logarithmic slowness.
        """
        if self.coordinate is not DensityCoordinate.LOG_SLOWNESS:
            raise ValueError(
                "Velocity is defined here only for log-slowness "
                "density coordinates."
            )

        return np.exp(
            -self.values
        )

    @property
    def slowness(self) -> FloatArray:
        """
        Return the slowness grid implied by log-slowness coordinates.

        Raises
        ------
        ValueError
            If the density coordinate is not logarithmic slowness.
        """
        if self.coordinate is not DensityCoordinate.LOG_SLOWNESS:
            raise ValueError(
                "Slowness is defined here only for log-slowness "
                "density coordinates."
            )

        return np.exp(
            self.values
        )

    def relative_density(
        self,
        *,
        floor: float = 0.0,
    ) -> FloatArray:
        """
        Return density normalized independently at each frequency.

        Each non-empty row is divided by its own maximum so that the
        strongest value equals one. This representation is useful for
        visualization because weak but coherent secondary modes are not
        suppressed by large absolute density differences among
        frequencies.

        Parameters
        ----------
        floor
            Optional lower clipping level in the interval ``[0, 1)``.
            A positive floor is convenient for logarithmic display.
        """
        if (
            not np.isfinite(floor)
            or not 0.0 <= floor < 1.0
        ):
            raise ValueError(
                "floor must be finite and in the range [0, 1)."
            )

        maximum = np.max(
            self.density,
            axis=1,
        )

        relative = np.zeros_like(
            self.density
        )

        valid = maximum > 0.0

        relative[valid] = (
            self.density[valid]
            / maximum[valid, None]
        )

        if floor > 0.0:
            relative[valid] = np.maximum(
                relative[valid],
                floor,
            )

        return relative

    def summary(
        self,
        *,
        interval: tuple[float, float] = (
            0.16,
            0.84,
        ),
    ) -> DensitySummary:
        """Return mode, median, and a conditional quantile interval."""
        return density_summary(
            self,
            interval=interval,
        )

    def peaks(
        self,
        *,
        number: int = 4,
        minimum_relative_density: float = 0.0,
        minimum_prominence: float = 0.0,
        minimum_separation: int = 1,
    ) -> DensityPeaks:
        """Return independently detected density maxima by frequency."""
        return density_peaks(
            self,
            number=number,
            minimum_relative_density=minimum_relative_density,
            minimum_prominence=minimum_prominence,
            minimum_separation=minimum_separation,
        )


def maximum_density(
    result: BeamformingResult,
    *,
    component: str | None = None,
    ranks: Iterable[int] | None = None,
    coordinate: DensityCoordinate | str = (
        DensityCoordinate.WAVENUMBER
    ),
    weighting: DensityWeighting | str = (
        DensityWeighting.UNIFORM
    ),
    size: int = 256,
    bandwidth_mode: BandwidthMode | str = BandwidthMode.RESOLUTION,
    bandwidth: float | None = None,
    bandwidth_factor: float = 1.0,
    use_uncertainty: bool = True,
    normalize: bool = True,
) -> MaximumDensity:
    """
    Estimate an empirical density from extracted spectral maxima.

    Parameters
    ----------
    result
        Beamforming result containing extracted maxima.
    component
        Component to analyse. If omitted, it is inferred only when the
        result contains maxima for exactly one component.
    ranks
        Optional one-based maximum ranks to retain. Rank acts only as a
        selection criterion; it is never used as a statistical weight.
    coordinate
        ``"wavenumber"``, ``"log_slowness"``, or
        ``"log_velocity"``.
    weighting
        ``"uniform"``, ``"power"``, or ``"relative_power"``.
        Relative power is normalized by the strongest extracted
        maximum in the same frequency, block, and component.
    size
        Number of samples in the common density-coordinate grid.
    bandwidth_mode
        ``"resolution"`` ties the kernel width to the local radial
        sampling of the beamforming search grid and is the default.
        ``"adaptive"`` uses a robust weighted Silverman estimate.
        ``"fixed"`` uses the explicit value supplied by ``bandwidth``.
    bandwidth
        Explicit Gaussian kernel bandwidth required by ``"fixed"``.
    bandwidth_factor
        Positive multiplier applied to resolution-based bandwidths.
        Values near one correspond to approximately one local radial
        grid step.
    use_uncertainty
        Fold available positional standard errors into individual
        kernels. Current YArray maxima store power uncertainty but not
        coordinate uncertainty, so this option presently leaves the
        kernel widths unchanged.
    normalize
        Normalize every non-empty frequency row to unit integral.

    Returns
    -------
    MaximumDensity
        Conditional empirical density and associated metadata.

    Notes
    -----
    Power uncertainty is not positional uncertainty and is therefore
    not used to broaden the density kernels.

    Future maxima exposing ``wavenumber_standard_error`` or
    ``slowness_standard_error`` are supported automatically.
    """
    if not isinstance(
        result,
        BeamformingResult,
    ):
        raise TypeError(
            "result must be a BeamformingResult instance."
        )

    coordinate = DensityCoordinate(
        coordinate
    )

    weighting = DensityWeighting(
        weighting
    )

    if size < 16:
        raise ValueError(
            "Density grid size must be at least 16."
        )

    bandwidth_mode = BandwidthMode(
        bandwidth_mode
    )

    if (
        not np.isfinite(bandwidth_factor)
        or bandwidth_factor <= 0.0
    ):
        raise ValueError(
            "bandwidth_factor must be finite and positive."
        )

    if bandwidth_mode is BandwidthMode.FIXED:
        if bandwidth is None:
            raise ValueError(
                "Fixed bandwidth mode requires bandwidth."
            )

        if (
            not np.isfinite(bandwidth)
            or bandwidth <= 0.0
        ):
            raise ValueError(
                "bandwidth must be finite and positive."
            )

    elif bandwidth is not None:
        raise ValueError(
            "bandwidth can be supplied only with fixed bandwidth mode."
        )

    selected_ranks = _validated_ranks(
        ranks
    )

    maxima = tuple(
        result.maxima
    )

    if not maxima:
        raise ValueError(
            "The beamforming result contains no spectral maxima."
        )

    selected_component = _resolve_component(
        maxima,
        component,
    )

    selected = tuple(
        maximum
        for maximum in maxima
        if maximum.component == selected_component
        and (
            selected_ranks is None
            or maximum.rank in selected_ranks
        )
    )

    if not selected:
        raise ValueError(
            "No spectral maxima match the requested selection."
        )

    frequencies = np.unique(
        np.asarray(
            [
                maximum.frequency
                for maximum in selected
            ],
            dtype=float,
        )
    )

    groups = [
        tuple(
            maximum
            for maximum in selected
            if maximum.frequency == frequency
        )
        for frequency in frequencies
    ]

    group_values = [
        np.asarray(
            [
                _coordinate_value(
                    maximum,
                    coordinate,
                )
                for maximum in group
            ],
            dtype=float,
        )
        for group in groups
    ]

    group_weights = [
        _maximum_weights(
            group,
            maxima,
            weighting,
        )
        for group in groups
    ]

    bandwidths = np.asarray(
        [
            _density_bandwidth(
                result,
                float(frequency),
                values,
                weights,
                coordinate,
                bandwidth_mode,
                bandwidth,
                bandwidth_factor,
            )
            for frequency, values, weights in zip(
                frequencies,
                group_values,
                group_weights,
            )
        ],
        dtype=float,
    )

    values = _density_grid(
        group_values,
        bandwidths,
        size,
    )

    density = np.zeros(
        (
            frequencies.size,
            values.size,
        ),
        dtype=float,
    )

    counts = np.asarray(
        [
            len(group)
            for group in groups
        ],
        dtype=np.int64,
    )

    for index, (
        group,
        samples,
        weights,
        empirical_bandwidth,
    ) in enumerate(
        zip(
            groups,
            group_values,
            group_weights,
            bandwidths,
        )
    ):
        uncertainty = np.asarray(
            [
                _coordinate_standard_error(
                    maximum,
                    coordinate,
                )
                if use_uncertainty
                else 0.0
                for maximum in group
            ],
            dtype=float,
        )

        sigma = np.sqrt(
            empirical_bandwidth**2
            + uncertainty**2
        )

        density[index] = _gaussian_mixture(
            values,
            samples,
            weights,
            sigma,
        )

        if normalize:
            integral = np.trapz(
                density[index],
                values,
            )

            if integral > 0.0:
                density[index] /= integral

    return MaximumDensity(
        frequencies=frequencies,
        values=values,
        density=density,
        counts=counts,
        bandwidths=bandwidths,
        bandwidth_mode=bandwidth_mode,
        bandwidth_factor=float(
            bandwidth_factor
        ),
        component=selected_component,
        coordinate=coordinate,
        weighting=weighting,
        ranks=selected_ranks,
    )


def density_peaks(
    density: MaximumDensity,
    *,
    number: int = 4,
    minimum_relative_density: float = 0.0,
    minimum_prominence: float = 0.0,
    minimum_separation: int = 1,
) -> DensityPeaks:
    """
    Extract the strongest local density maxima at each frequency.

    Parameters
    ----------
    density
        Maximum-density estimate.
    number
        Maximum number of local maxima retained per frequency.
    minimum_relative_density
        Minimum peak height relative to the strongest density maximum
        at the same frequency.
    minimum_prominence
        Minimum prominence relative to the strongest density maximum
        at the same frequency.
    minimum_separation
        Minimum separation between retained peaks in density-grid
        samples.

    Returns
    -------
    DensityPeaks
        Independently detected maxima sorted by decreasing density
        amplitude at every frequency.

    Notes
    -----
    No continuity constraint is imposed between adjacent frequencies.
    This function therefore identifies candidate modes but does not
    perform modal tracking or curve extraction.
    """
    if not isinstance(
        density,
        MaximumDensity,
    ):
        raise TypeError(
            "density must be a MaximumDensity instance."
        )

    if isinstance(number, bool) or number < 1:
        raise ValueError(
            "number must be a positive integer."
        )

    if not 0.0 <= minimum_relative_density <= 1.0:
        raise ValueError(
            "minimum_relative_density must be in the range [0, 1]."
        )

    if not 0.0 <= minimum_prominence <= 1.0:
        raise ValueError(
            "minimum_prominence must be in the range [0, 1]."
        )

    if (
        isinstance(
            minimum_separation,
            bool,
        )
        or minimum_separation < 1
    ):
        raise ValueError(
            "minimum_separation must be a positive integer."
        )

    shape = (
        density.frequencies.size,
        number,
    )

    peak_values = np.full(
        shape,
        np.nan,
        dtype=float,
    )

    peak_density = np.full(
        shape,
        np.nan,
        dtype=float,
    )

    relative_density = np.full(
        shape,
        np.nan,
        dtype=float,
    )

    prominence = np.full(
        shape,
        np.nan,
        dtype=float,
    )

    count = np.zeros(
        density.frequencies.size,
        dtype=np.int64,
    )

    for frequency_index, row in enumerate(
        density.density
    ):
        row_maximum = float(
            np.max(
                row
            )
        )

        if not np.isfinite(row_maximum) or row_maximum <= 0.0:
            continue

        indices = _local_maximum_indices(
            row
        )

        if indices.size == 0:
            indices = np.asarray(
                [
                    int(
                        np.argmax(
                            row
                        )
                    )
                ],
                dtype=np.int64,
            )

        relative = (
            row[indices]
            / row_maximum
        )

        prominences = np.asarray(
            [
                _peak_prominence(
                    row,
                    int(index),
                )
                / row_maximum
                for index in indices
            ],
            dtype=float,
        )

        valid = (
            relative
            >= minimum_relative_density
        ) & (
            prominences
            >= minimum_prominence
        )

        indices = indices[
            valid
        ]

        relative = relative[
            valid
        ]

        prominences = prominences[
            valid
        ]

        if indices.size == 0:
            continue

        order = np.argsort(
            row[indices]
        )[::-1]

        indices = indices[
            order
        ]

        relative = relative[
            order
        ]

        prominences = prominences[
            order
        ]

        selected: list[int] = []
        selected_relative: list[float] = []
        selected_prominence: list[float] = []

        for index, rel, prom in zip(
            indices,
            relative,
            prominences,
        ):
            if any(
                abs(
                    int(index)
                    - previous
                )
                < minimum_separation
                for previous in selected
            ):
                continue

            selected.append(
                int(index)
            )

            selected_relative.append(
                float(rel)
            )

            selected_prominence.append(
                float(prom)
            )

            if len(selected) >= number:
                break

        retained = len(
            selected
        )

        if retained == 0:
            continue

        selected_indices = np.asarray(
            selected,
            dtype=np.int64,
        )

        peak_values[
            frequency_index,
            :retained,
        ] = density.values[
            selected_indices
        ]

        peak_density[
            frequency_index,
            :retained,
        ] = row[
            selected_indices
        ]

        relative_density[
            frequency_index,
            :retained,
        ] = np.asarray(
            selected_relative,
            dtype=float,
        )

        prominence[
            frequency_index,
            :retained,
        ] = np.asarray(
            selected_prominence,
            dtype=float,
        )

        count[
            frequency_index
        ] = retained

    return DensityPeaks(
        frequencies=density.frequencies.copy(),
        values=peak_values,
        relative_density=relative_density,
        density=peak_density,
        prominence=prominence,
        count=count,
    )


def density_summary(
    density: MaximumDensity,
    *,
    interval: tuple[float, float] = (
        0.16,
        0.84,
    ),
) -> DensitySummary:
    """
    Summarize each conditional density by mode and quantiles.

    Parameters
    ----------
    density
        Maximum-density estimate.
    interval
        Lower and upper cumulative probabilities. The default
        ``(0.16, 0.84)`` is analogous to a central 68-percent interval
        for a Gaussian distribution.

    Returns
    -------
    DensitySummary
        Frequency-dependent mode, median, and conditional interval.
    """
    if not isinstance(
        density,
        MaximumDensity,
    ):
        raise TypeError(
            "density must be a MaximumDensity instance."
        )

    if len(interval) != 2:
        raise ValueError(
            "interval must contain exactly two probabilities."
        )

    lower_probability = float(
        interval[0]
    )

    upper_probability = float(
        interval[1]
    )

    if not (
        0.0
        < lower_probability
        < 0.5
        < upper_probability
        < 1.0
    ):
        raise ValueError(
            "interval must satisfy 0 < lower < 0.5 < upper < 1."
        )

    mode = np.full(
        density.frequencies.size,
        np.nan,
        dtype=float,
    )

    median = np.full_like(
        mode,
        np.nan,
    )

    lower = np.full_like(
        mode,
        np.nan,
    )

    upper = np.full_like(
        mode,
        np.nan,
    )

    for index, row in enumerate(
        density.density
    ):
        integral = float(
            np.trapz(
                row,
                density.values,
            )
        )

        if not np.isfinite(integral) or integral <= 0.0:
            continue

        mode[index] = density.values[
            int(
                np.argmax(
                    row
                )
            )
        ]

        cumulative = _density_cumulative(
            density.values,
            row,
        )

        lower[index] = np.interp(
            lower_probability,
            cumulative,
            density.values,
        )

        median[index] = np.interp(
            0.5,
            cumulative,
            density.values,
        )

        upper[index] = np.interp(
            upper_probability,
            cumulative,
            density.values,
        )

    return DensitySummary(
        frequencies=density.frequencies.copy(),
        mode=mode,
        median=median,
        lower=lower,
        upper=upper,
        lower_probability=lower_probability,
        upper_probability=upper_probability,
    )


def _density_cumulative(
    values: FloatArray,
    density: FloatArray,
) -> FloatArray:
    """Return a normalized CDF from density samples on a regular grid."""
    increments = (
        0.5
        * (
            density[1:]
            + density[:-1]
        )
        * np.diff(
            values
        )
    )

    cumulative = np.concatenate(
        (
            np.array(
                [0.0],
                dtype=float,
            ),
            np.cumsum(
                increments
            ),
        )
    )

    total = cumulative[-1]

    if total <= 0.0:
        return cumulative

    return cumulative / total


def _validated_ranks(
    ranks: Iterable[int] | None,
) -> tuple[int, ...] | None:
    """Return a validated tuple of unique one-based ranks."""
    if ranks is None:
        return None

    values = tuple(
        int(rank)
        for rank in ranks
    )

    if not values:
        raise ValueError(
            "ranks cannot be empty."
        )

    if any(
        rank < 1
        for rank in values
    ):
        raise ValueError(
            "ranks must contain positive one-based integers."
        )

    return tuple(
        dict.fromkeys(
            values
        )
    )


def _resolve_component(
    maxima: tuple[SpectralMaximum, ...],
    component: str | None,
) -> str:
    """Resolve or validate the component represented by the density."""
    components = tuple(
        sorted(
            {
                maximum.component
                for maximum in maxima
            }
        )
    )

    if component is None:
        if len(components) != 1:
            names = ", ".join(
                components
            )

            raise ValueError(
                "component must be specified because the result "
                f"contains multiple components: {names}."
            )

        return components[0]

    if component not in components:
        names = ", ".join(
            components
        )

        raise ValueError(
            f"Unknown component {component!r}. "
            f"Available components: {names}."
        )

    return component


def _coordinate_value(
    maximum: SpectralMaximum,
    coordinate: DensityCoordinate,
) -> float:
    """Return one maximum in the selected statistical coordinate."""
    if coordinate is DensityCoordinate.WAVENUMBER:
        value = maximum.wavenumber

    elif coordinate is DensityCoordinate.LOG_SLOWNESS:
        slowness = maximum.slowness

        if (
            not np.isfinite(slowness)
            or slowness <= 0.0
        ):
            raise ValueError(
                "Log-slowness density requires finite positive "
                "slowness values."
            )

        value = np.log(
            slowness
        )

    else:
        velocity = maximum.velocity

        if (
            not np.isfinite(velocity)
            or velocity <= 0.0
        ):
            raise ValueError(
                "Log-velocity density requires finite positive "
                "velocity values."
            )

        value = np.log(
            velocity
        )

    if not np.isfinite(value):
        raise ValueError(
            "Maximum coordinate contains a non-finite value."
        )

    return float(
        value
    )


def _maximum_weights(
    group: tuple[SpectralMaximum, ...],
    all_maxima: tuple[SpectralMaximum, ...],
    weighting: DensityWeighting,
) -> FloatArray:
    """Return normalized statistical weights for one frequency group."""
    if weighting is DensityWeighting.UNIFORM:
        weights = np.ones(
            len(group),
            dtype=float,
        )

    elif weighting is DensityWeighting.POWER:
        weights = np.asarray(
            [
                max(
                    0.0,
                    float(maximum.power),
                )
                for maximum in group
            ],
            dtype=float,
        )

    else:
        weights = np.asarray(
            [
                _relative_power(
                    maximum,
                    all_maxima,
                )
                for maximum in group
            ],
            dtype=float,
        )

    total = float(
        np.sum(weights)
    )

    if not np.isfinite(total) or total <= 0.0:
        return np.full(
            len(group),
            1.0 / len(group),
            dtype=float,
        )

    return weights / total


def _relative_power(
    maximum: SpectralMaximum,
    all_maxima: tuple[SpectralMaximum, ...],
) -> float:
    """
    Return power relative to the strongest maximum in one spectrum.

    The denominator is evaluated before rank filtering so the weight
    retains the meaning of power relative to the strongest extracted
    maximum of that frequency, block, and component.
    """
    reference = max(
        (
            item.power
            for item in all_maxima
            if item.frequency == maximum.frequency
            and item.block == maximum.block
            and item.component == maximum.component
        ),
        default=0.0,
    )

    if not np.isfinite(reference) or reference <= 0.0:
        return 0.0

    return max(
        0.0,
        float(maximum.power),
    ) / reference


def _density_bandwidth(
    result: BeamformingResult,
    frequency: float,
    values: FloatArray,
    weights: FloatArray,
    coordinate: DensityCoordinate,
    mode: BandwidthMode,
    bandwidth: float | None,
    factor: float,
) -> float:
    """Return the bandwidth selected by the requested strategy."""
    if mode is BandwidthMode.FIXED:
        return float(
            bandwidth
        )

    if mode is BandwidthMode.ADAPTIVE:
        return _automatic_bandwidth(
            values,
            weights,
        )

    return (
        factor
        * _resolution_bandwidth(
            result,
            frequency,
            values,
            coordinate,
        )
    )


def _resolution_bandwidth(
    result: BeamformingResult,
    frequency: float,
    samples: FloatArray,
    coordinate: DensityCoordinate,
) -> float:
    """
    Estimate density bandwidth from local radial search-grid spacing.

    The configured radial grid is transformed to the statistical
    coordinate at the current frequency. For each observed maximum, the
    nearest grid sample is found and its local one-step spacing is
    estimated from adjacent grid values. The median local spacing of
    the observed maxima defines the density bandwidth.

    This choice reflects numerical resolving power of the search grid
    and does not grow merely because several distinct modal populations
    coexist at one frequency.
    """
    grid = _coordinate_grid(
        result,
        frequency,
        coordinate,
    )

    grid = np.unique(
        np.sort(
            grid
        )
    )

    if grid.size < 2:
        raise ValueError(
            "At least two distinct radial grid samples are required "
            "for resolution-based density estimation."
        )

    local_spacing = np.empty(
        grid.size,
        dtype=float,
    )

    local_spacing[0] = (
        grid[1]
        - grid[0]
    )

    local_spacing[-1] = (
        grid[-1]
        - grid[-2]
    )

    if grid.size > 2:
        local_spacing[1:-1] = 0.5 * (
            grid[2:]
            - grid[:-2]
        )

    indices = _nearest_sorted_indices(
        grid,
        samples,
    )

    spacing = local_spacing[
        indices
    ]

    spacing = spacing[
        np.isfinite(spacing)
        & (
            spacing > 0.0
        )
    ]

    if spacing.size == 0:
        raise ValueError(
            "Cannot estimate a positive search-grid resolution."
        )

    return float(
        np.median(
            spacing
        )
    )


def _coordinate_grid(
    result: BeamformingResult,
    frequency: float,
    coordinate: DensityCoordinate,
) -> FloatArray:
    """Return the configured radial grid in a statistical coordinate."""
    grid = result.config.grid

    if coordinate is DensityCoordinate.WAVENUMBER:
        values = grid.wavenumbers(
            frequency
        )

    elif coordinate is DensityCoordinate.LOG_SLOWNESS:
        slowness = grid.slowness(
            frequency
        )

        valid = (
            np.isfinite(slowness)
            & (
                slowness > 0.0
            )
        )

        values = np.log(
            slowness[
                valid
            ]
        )

    else:
        velocity = grid.velocity(
            frequency
        )

        valid = (
            np.isfinite(velocity)
            & (
                velocity > 0.0
            )
        )

        values = np.log(
            velocity[
                valid
            ]
        )

    values = np.asarray(
        values,
        dtype=float,
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if values.size < 2:
        raise ValueError(
            "The configured radial grid cannot be transformed into "
            "the requested density coordinate."
        )

    return values


def _nearest_sorted_indices(
    grid: FloatArray,
    samples: FloatArray,
) -> IntArray:
    """Return nearest indices on an increasing one-dimensional grid."""
    positions = np.searchsorted(
        grid,
        samples,
        side="left",
    )

    positions = np.clip(
        positions,
        1,
        grid.size - 1,
    )

    left = positions - 1
    right = positions

    choose_right = (
        samples - grid[left]
        > grid[right] - samples
    )

    return np.where(
        choose_right,
        right,
        left,
    ).astype(
        np.int64
    )


def _automatic_bandwidth(
    values: FloatArray,
    weights: FloatArray,
) -> float:
    """
    Return a robust weighted Silverman-style bandwidth estimate.

    A conservative scale fallback is used for one-point or degenerate
    frequency groups.
    """
    values = np.asarray(
        values,
        dtype=float,
    )

    weights = np.asarray(
        weights,
        dtype=float,
    )

    count = values.size

    if count == 0:
        raise ValueError(
            "Cannot estimate bandwidth from an empty sample."
        )

    center = float(
        np.sum(
            weights * values
        )
    )

    variance = float(
        np.sum(
            weights
            * (
                values - center
            ) ** 2
        )
    )

    standard_deviation = np.sqrt(
        max(
            variance,
            0.0,
        )
    )

    q25, q75 = _weighted_quantiles(
        values,
        weights,
        (
            0.25,
            0.75,
        ),
    )

    robust_scale = (
        q75 - q25
    ) / 1.349

    candidates = [
        value
        for value in (
            standard_deviation,
            robust_scale,
        )
        if np.isfinite(value)
        and value > 0.0
    ]

    if candidates:
        scale = min(
            candidates
        )

        effective_count = 1.0 / float(
            np.sum(
                weights**2
            )
        )

        bandwidth = (
            0.9
            * scale
            * effective_count ** (-0.2)
        )

        if np.isfinite(bandwidth) and bandwidth > 0.0:
            return float(
                bandwidth
            )

    magnitude = max(
        1.0,
        abs(center),
    )

    return float(
        0.02 * magnitude
    )


def _weighted_quantiles(
    values: FloatArray,
    weights: FloatArray,
    quantiles: tuple[float, ...],
) -> FloatArray:
    """Return weighted quantiles of one-dimensional samples."""
    order = np.argsort(
        values
    )

    sorted_values = values[
        order
    ]

    sorted_weights = weights[
        order
    ]

    cumulative = np.cumsum(
        sorted_weights
    )

    cumulative /= cumulative[-1]

    return np.interp(
        np.asarray(
            quantiles,
            dtype=float,
        ),
        cumulative,
        sorted_values,
    )


def _density_grid(
    groups: list[FloatArray],
    bandwidths: FloatArray,
    size: int,
) -> FloatArray:
    """Return a common grid covering all kernels."""
    lower = min(
        float(
            np.min(values)
            - 4.0 * bandwidth
        )
        for values, bandwidth in zip(
            groups,
            bandwidths,
        )
    )

    upper = max(
        float(
            np.max(values)
            + 4.0 * bandwidth
        )
        for values, bandwidth in zip(
            groups,
            bandwidths,
        )
    )

    if not np.isfinite(lower) or not np.isfinite(upper):
        raise ValueError(
            "Cannot construct a finite density grid."
        )

    if upper <= lower:
        scale = max(
            1.0,
            abs(lower),
        )

        lower -= 0.01 * scale
        upper += 0.01 * scale

    return np.linspace(
        lower,
        upper,
        size,
        dtype=float,
    )


def _coordinate_standard_error(
    maximum: SpectralMaximum,
    coordinate: DensityCoordinate,
) -> float:
    """
    Return available positional uncertainty in density coordinates.

    Current :class:`SpectralMaximum` objects do not yet expose
    wavenumber/slowness standard errors, so this function normally
    returns zero. The attribute-based implementation is intentionally
    forward compatible with future uncertainty propagation.
    """
    wavenumber_error = getattr(
        maximum,
        "wavenumber_standard_error",
        None,
    )

    if coordinate is DensityCoordinate.WAVENUMBER:
        return _validated_optional_error(
            wavenumber_error
        )

    if coordinate is DensityCoordinate.LOG_VELOCITY:
        velocity_error = getattr(
            maximum,
            "velocity_standard_error",
            None,
        )

        if velocity_error is not None:
            error = _validated_optional_error(
                velocity_error
            )

            velocity = maximum.velocity

            if velocity <= 0.0:
                return 0.0

            return error / velocity

        if wavenumber_error is not None:
            error = _validated_optional_error(
                wavenumber_error
            )

            if maximum.wavenumber <= 0.0:
                return 0.0

            return error / maximum.wavenumber

        return 0.0

    slowness_error = getattr(
        maximum,
        "slowness_standard_error",
        None,
    )

    if slowness_error is not None:
        error = _validated_optional_error(
            slowness_error
        )

        slowness = maximum.slowness

        if slowness <= 0.0:
            return 0.0

        return error / slowness

    if wavenumber_error is not None:
        error = _validated_optional_error(
            wavenumber_error
        )

        if maximum.wavenumber <= 0.0:
            return 0.0

        return error / maximum.wavenumber

    return 0.0


def _validated_optional_error(
    value: object,
) -> float:
    """Return a valid non-negative optional standard error."""
    if value is None:
        return 0.0

    error = float(
        value
    )

    if not np.isfinite(error) or error < 0.0:
        return 0.0

    return error


def _local_maximum_indices(
    values: FloatArray,
) -> IntArray:
    """Return local maxima, including valid one-sided edge maxima."""
    values = np.asarray(
        values,
        dtype=float,
    )

    if values.size == 1:
        return np.asarray(
            [0],
            dtype=np.int64,
        )

    interior = np.flatnonzero(
        (
            values[1:-1]
            > values[:-2]
        )
        & (
            values[1:-1]
            >= values[2:]
        )
    ) + 1

    edges: list[int] = []

    if values[0] > values[1]:
        edges.append(
            0
        )

    if values[-1] > values[-2]:
        edges.append(
            values.size - 1
        )

    if not edges:
        return interior.astype(
            np.int64
        )

    return np.concatenate(
        (
            np.asarray(
                edges,
                dtype=np.int64,
            ),
            interior.astype(
                np.int64
            ),
        )
    )


def _peak_prominence(
    values: FloatArray,
    index: int,
) -> float:
    """
    Estimate one-dimensional peak prominence without SciPy dependency.

    The reference level is the higher of the minimum values encountered
    to the left and right of the peak. At an edge, the available side
    determines the reference.
    """
    peak = float(
        values[index]
    )

    if values.size == 1:
        return peak

    if index == 0:
        reference = float(
            np.min(
                values[1:]
            )
        )

    elif index == values.size - 1:
        reference = float(
            np.min(
                values[:-1]
            )
        )

    else:
        left = float(
            np.min(
                values[:index]
            )
        )

        right = float(
            np.min(
                values[index + 1:]
            )
        )

        reference = max(
            left,
            right,
        )

    return max(
        0.0,
        peak - reference,
    )


def _gaussian_mixture(
    grid: FloatArray,
    samples: FloatArray,
    weights: FloatArray,
    sigma: FloatArray,
) -> FloatArray:
    """Evaluate a weighted Gaussian mixture on a common grid."""
    sigma = np.maximum(
        np.asarray(
            sigma,
            dtype=float,
        ),
        np.finfo(float).eps,
    )

    offset = (
        grid[:, None]
        - samples[None, :]
    ) / sigma[None, :]

    kernels = (
        np.exp(
            -0.5 * offset**2
        )
        / (
            np.sqrt(
                2.0 * np.pi
            )
            * sigma[None, :]
        )
    )

    return kernels @ weights
