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
Statistics of spectral maxima.

The maxima of many covariance blocks (long recordings, or the results
of the real-time monitor joined) scatter around the dispersion curves
of the waves crossing the array. :func:`maximum_density` estimates, at
each frequency, the density of the maxima along the radial coordinate
(wavenumber, log-slowness or log-velocity) as a weighted mixture of
Gaussian kernels; :meth:`MaximumDensity.summary` gives mode, median
and a quantile interval, :meth:`MaximumDensity.peaks` the local maxima
of the density (candidate modes).

Design choices
--------------
By default (``bandwidth_mode="resolution"``) the kernel width is the
local spacing of the radial search grid in the chosen coordinate: it
reflects the resolving power of the analysis and does not widen when
several modes coexist at one frequency, as a bandwidth from the spread
of the samples would. ``"adaptive"`` (weighted Silverman rule) and
``"fixed"`` are alternatives.

The density is an empirical conditional density of the maxima at each
frequency, not a probability of the true dispersion curve. Peaks are
found independently at each frequency (no modal tracking).

References
----------
Silverman, B. W. (1986). Density Estimation for Statistics and Data
Analysis. Chapman and Hall, Section 3.4.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

import numpy as np
from scipy.signal import find_peaks

from shakelab.exploration.yarray.beamforming import (
    BeamformingResult,
    SpectralMaximum,
)

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
    """Width of the density kernels."""

    RESOLUTION = "resolution"
    ADAPTIVE = "adaptive"
    FIXED = "fixed"


class DensityCoordinate(Enum):
    """Radial coordinate of the density."""

    WAVENUMBER = "wavenumber"
    LOG_SLOWNESS = "log_slowness"
    LOG_VELOCITY = "log_velocity"


class DensityWeighting(Enum):
    """Weight of each maximum."""

    UNIFORM = "uniform"
    POWER = "power"
    RELATIVE_POWER = "relative_power"


def _enum(value: Any, kind: type, name: str) -> Enum:
    """Return an enumeration member from a member or its value."""
    if isinstance(value, kind):
        return value

    for member in kind:
        if member.value == value:
            return member

    names = ", ".join(repr(member.value) for member in kind)
    raise ValueError(f"{name} must be one of {names}")


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class DensitySummary:
    """
    Mode, median and quantile interval of the density at each frequency.

    Attributes
    ----------
    frequencies : numpy.ndarray
        Frequencies in Hz.
    mode, median, lower, upper : numpy.ndarray
        Coordinate of the density maximum, of the median and of the
        lower and upper quantiles (NaN where there are no maxima).
    probabilities : tuple of float
        Cumulative probabilities of ``lower`` and ``upper``.

    Notes
    -----
    Median and interval describe the whole density: with several modes
    they may fall between them, where the mode does not.
    """

    frequencies: np.ndarray
    mode: np.ndarray
    median: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    probabilities: tuple


@dataclass(frozen=True, eq=False)
class DensityPeaks:
    """
    Local maxima of the density at each frequency.

    Attributes
    ----------
    frequencies : numpy.ndarray
        Frequencies in Hz.
    values : numpy.ndarray
        Coordinates of the peaks ``(n_frequencies, number)``, by
        decreasing height; NaN where there are fewer peaks.
    relative_density : numpy.ndarray
        Heights divided by the largest density at the frequency.
    prominence : numpy.ndarray
        Prominences divided by the largest density at the frequency.
    count : numpy.ndarray
        Number of peaks at each frequency.

    Notes
    -----
    The column of a peak is its rank by height at that frequency, not a
    mode followed across frequencies.
    """

    frequencies: np.ndarray
    values: np.ndarray
    relative_density: np.ndarray
    prominence: np.ndarray
    count: np.ndarray


@dataclass(frozen=True, eq=False)
class MaximumDensity:
    """
    Density of the spectral maxima at each frequency.

    Attributes
    ----------
    frequencies : numpy.ndarray
        Frequencies with maxima, in Hz.
    values : numpy.ndarray
        Coordinate grid, common to all frequencies.
    density : numpy.ndarray
        Density ``(n_frequencies, n_values)``; each row integrates to
        one over ``values``.
    counts : numpy.ndarray
        Number of maxima at each frequency.
    bandwidths : numpy.ndarray
        Kernel width at each frequency, in the coordinate.
    component : str
        Component of the maxima.
    coordinate : DensityCoordinate
        Coordinate of ``values``: wavenumber (rad per distance unit),
        natural logarithm of slowness or of velocity.
    weighting : DensityWeighting
        Weight of the maxima.
    ranks : tuple of int or None
        Ranks of the maxima used (None for all).
    """

    frequencies: np.ndarray
    values: np.ndarray
    density: np.ndarray
    counts: np.ndarray
    bandwidths: np.ndarray
    component: str
    coordinate: DensityCoordinate
    weighting: DensityWeighting
    ranks: Optional[tuple]

    @property
    def velocity(self) -> np.ndarray:
        """
        Phase velocity of the grid (logarithmic coordinates only).

        Raises
        ------
        ValueError
            For the wavenumber coordinate, whose velocity depends on
            the frequency.
        """
        if self.coordinate is DensityCoordinate.LOG_VELOCITY:
            return np.exp(self.values)

        if self.coordinate is DensityCoordinate.LOG_SLOWNESS:
            return np.exp(-self.values)

        raise ValueError("the velocity of a wavenumber grid depends on f")

    def relative_density(self, floor: float = 0.0) -> np.ndarray:
        """
        Return the density divided by its maximum at each frequency.

        Parameters
        ----------
        floor : float, optional
            Smallest value, in ``[0, 1)`` (e.g. for a logarithmic
            colour scale). Default 0.

        Returns
        -------
        numpy.ndarray
            Values in ``[floor, 1]``, shape of :attr:`density`.
        """
        if not 0.0 <= floor < 1.0:
            raise ValueError("floor must be in [0, 1)")

        peak = np.max(self.density, axis=1, keepdims=True)
        relative = np.divide(self.density, peak, out=np.zeros_like(
            self.density), where=peak > 0.0)
        return np.maximum(relative, floor)

    def summary(self, interval: tuple = (0.16, 0.84)) -> DensitySummary:
        """
        Return mode, median and a quantile interval.

        Parameters
        ----------
        interval : tuple of float, optional
            Lower and upper cumulative probabilities. Default (0.16,
            0.84), the central 68% of a Gaussian.

        Returns
        -------
        DensitySummary
            Summary at each frequency.
        """
        return density_summary(self, interval)

    def peaks(
        self,
        number: int = 4,
        minimum_relative_density: float = 0.0,
        minimum_prominence: float = 0.0,
        minimum_separation: int = 1,
    ) -> DensityPeaks:
        """
        Return the local maxima of the density.

        Parameters
        ----------
        number : int, optional
            Peaks kept at each frequency. Default 4.
        minimum_relative_density : float, optional
            Smallest height relative to the largest. Default 0.
        minimum_prominence : float, optional
            Smallest prominence relative to the largest height.
            Default 0.
        minimum_separation : int, optional
            Smallest distance between peaks, in grid samples. Default 1.

        Returns
        -------
        DensityPeaks
            Peaks at each frequency.
        """
        return density_peaks(
            self, number, minimum_relative_density, minimum_prominence,
            minimum_separation,
        )

    def info(self) -> None:
        """Print a summary of the density."""
        print(
            f"MaximumDensity: {self.component}, {self.coordinate.value}, "
            f"{self.weighting.value}, {len(self.frequencies)} "
            f"frequencies, {int(np.sum(self.counts))} maxima, grid of "
            f"{len(self.values)} values"
        )


# ----------------------------------------------------------------------------
# Density
# ----------------------------------------------------------------------------


def maximum_density(
    result: BeamformingResult,
    component: Optional[str] = None,
    ranks: Optional[Iterable] = None,
    coordinate: Any = DensityCoordinate.WAVENUMBER,
    weighting: Any = DensityWeighting.UNIFORM,
    size: int = 256,
    bandwidth_mode: Any = BandwidthMode.RESOLUTION,
    bandwidth: Optional[float] = None,
    bandwidth_factor: float = 1.0,
) -> MaximumDensity:
    """
    Estimate the density of the spectral maxima at each frequency.

    Parameters
    ----------
    result : BeamformingResult
        Result with maxima (e.g. joined from many windows).
    component : str, optional
        Component; required if several were processed.
    ranks : iterable of int, optional
        Ranks of the maxima used (a selection, not a weight). Default:
        all.
    coordinate : DensityCoordinate or str, optional
        ``"wavenumber"`` (default), ``"log_slowness"`` or
        ``"log_velocity"`` (natural logarithms).
    weighting : DensityWeighting or str, optional
        ``"uniform"`` (default), ``"power"``, or ``"relative_power"``
        (power divided by that of the strongest maximum of the same
        block, frequency and component).
    size : int, optional
        Number of values of the common grid (at least 16). Default 256.
    bandwidth_mode : BandwidthMode or str, optional
        ``"resolution"`` (default): local spacing of the search grid in
        the coordinate; ``"adaptive"``: weighted Silverman rule;
        ``"fixed"``: ``bandwidth``.
    bandwidth : float, optional
        Kernel width with ``"fixed"``, in the coordinate.
    bandwidth_factor : float, optional
        Factor applied to the resolution bandwidth. Default 1.

    Returns
    -------
    MaximumDensity
        Density and its parameters.

    Raises
    ------
    TypeError
        If ``result`` is not a BeamformingResult.
    ValueError
        If no maximum is selected, the component is ambiguous, or an
        option is invalid.
    """
    if not isinstance(result, BeamformingResult):
        raise TypeError("result must be a BeamformingResult")

    coordinate = _enum(coordinate, DensityCoordinate, "coordinate")
    weighting = _enum(weighting, DensityWeighting, "weighting")
    bandwidth_mode = _enum(bandwidth_mode, BandwidthMode, "bandwidth_mode")

    if int(size) < 16:
        raise ValueError("size must be at least 16")

    if not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0.0:
        raise ValueError("bandwidth_factor must be positive")

    if (bandwidth_mode is BandwidthMode.FIXED) != (bandwidth is not None):
        raise ValueError("bandwidth is required by, and only by, 'fixed'")

    if bandwidth is not None and not float(bandwidth) > 0.0:
        raise ValueError("bandwidth must be positive")

    ranks = None if ranks is None else tuple(sorted({int(r) for r in ranks}))

    if ranks is not None and (not ranks or ranks[0] < 1):
        raise ValueError("ranks must be positive integers")

    component = _component(result, component)
    selected = [
        item for item in result.maxima
        if item.component == component
        and (ranks is None or item.rank in ranks)
    ]

    if not selected:
        raise ValueError("no maximum matches the selection")

    frequencies = np.unique([item.frequency for item in selected])
    reference = _reference_powers(result.maxima)
    samples, weights, widths = [], [], []

    for frequency in frequencies:
        group = [item for item in selected if item.frequency == frequency]
        values = np.array([_coordinate(item, coordinate) for item in group])
        weight = _weights(group, weighting, reference)
        samples.append(values)
        weights.append(weight)

        if bandwidth_mode is BandwidthMode.FIXED:
            widths.append(float(bandwidth))
        elif bandwidth_mode is BandwidthMode.ADAPTIVE:
            widths.append(_silverman(values, weight))
        else:
            widths.append(bandwidth_factor * _grid_spacing(
                result, float(frequency), values, coordinate
            ))

    widths = np.array(widths)
    lower = min(np.min(v) - 4.0 * w for v, w in zip(samples, widths))
    upper = max(np.max(v) + 4.0 * w for v, w in zip(samples, widths))
    grid = np.linspace(lower, upper, int(size))
    density = np.array([
        _mixture(grid, values, weight, width)
        for values, weight, width in zip(samples, weights, widths)
    ])
    density /= np.trapezoid(density, grid, axis=1)[:, None]

    return MaximumDensity(
        frequencies=frequencies,
        values=grid,
        density=density,
        counts=np.array([len(values) for values in samples]),
        bandwidths=widths,
        component=component,
        coordinate=coordinate,
        weighting=weighting,
        ranks=ranks,
    )


def density_summary(
    density: MaximumDensity,
    interval: tuple = (0.16, 0.84),
) -> DensitySummary:
    """
    Return mode, median and a quantile interval of a density.

    Parameters
    ----------
    density : MaximumDensity
        Density of the maxima.
    interval : tuple of float, optional
        Lower and upper cumulative probabilities, with
        ``0 < lower < 0.5 < upper < 1``. Default (0.16, 0.84).

    Returns
    -------
    DensitySummary
        Summary at each frequency.

    Raises
    ------
    ValueError
        If the interval is invalid.
    """
    if not isinstance(density, MaximumDensity):
        raise TypeError("density must be a MaximumDensity")

    lower_p, upper_p = (float(value) for value in interval)

    if not 0.0 < lower_p < 0.5 < upper_p < 1.0:
        raise ValueError("interval must satisfy 0 < lower < 0.5 < upper < 1")

    grid = density.values
    steps = 0.5 * (density.density[:, 1:] + density.density[:, :-1]) * (
        np.diff(grid)
    )
    cumulative = np.concatenate(
        [np.zeros((len(steps), 1)), np.cumsum(steps, axis=1)], axis=1
    )
    cumulative /= cumulative[:, -1:]

    def quantile(probability: float) -> np.ndarray:
        return np.array([np.interp(probability, row, grid)
                         for row in cumulative])

    return DensitySummary(
        frequencies=density.frequencies.copy(),
        mode=grid[np.argmax(density.density, axis=1)],
        median=quantile(0.5),
        lower=quantile(lower_p),
        upper=quantile(upper_p),
        probabilities=(lower_p, upper_p),
    )


def density_peaks(
    density: MaximumDensity,
    number: int = 4,
    minimum_relative_density: float = 0.0,
    minimum_prominence: float = 0.0,
    minimum_separation: int = 1,
) -> DensityPeaks:
    """
    Return the local maxima of a density at each frequency.

    Peaks and their prominences are those of
    ``scipy.signal.find_peaks``, on each row padded with zeros so that
    maxima at the ends of the grid are found too.

    Parameters
    ----------
    density : MaximumDensity
        Density of the maxima.
    number : int, optional
        Peaks kept at each frequency, by decreasing height. Default 4.
    minimum_relative_density : float, optional
        Smallest height relative to the largest, in ``[0, 1]``.
    minimum_prominence : float, optional
        Smallest prominence relative to the largest height, in
        ``[0, 1]``.
    minimum_separation : int, optional
        Smallest distance between peaks in grid samples (the higher
        peak is kept). Default 1.

    Returns
    -------
    DensityPeaks
        Peaks at each frequency.

    Raises
    ------
    ValueError
        If an option is out of range.
    """
    if not isinstance(density, MaximumDensity):
        raise TypeError("density must be a MaximumDensity")

    if int(number) < 1 or int(minimum_separation) < 1:
        raise ValueError("number and minimum_separation must be >= 1")

    for value in (minimum_relative_density, minimum_prominence):
        if not 0.0 <= value <= 1.0:
            raise ValueError("relative thresholds must be in [0, 1]")

    shape = (len(density.frequencies), int(number))
    values = np.full(shape, np.nan)
    relative = np.full(shape, np.nan)
    prominence = np.full(shape, np.nan)
    count = np.zeros(shape[0], dtype=int)

    for row_index, row in enumerate(density.density):
        top = np.max(row)
        padded = np.concatenate([[0.0], row, [0.0]])
        index, properties = find_peaks(
            padded,
            height=minimum_relative_density * top,
            prominence=max(minimum_prominence * top, 0.0),
            distance=int(minimum_separation),
        )
        order = np.argsort(padded[index])[::-1][:int(number)]
        found = order.size
        values[row_index, :found] = density.values[index[order] - 1]
        relative[row_index, :found] = padded[index[order]] / top
        prominence[row_index, :found] = (
            properties["prominences"][order] / top
        )
        count[row_index] = found

    return DensityPeaks(
        frequencies=density.frequencies.copy(),
        values=values,
        relative_density=relative,
        prominence=prominence,
        count=count,
    )


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _component(result: BeamformingResult, component: Optional[str]) -> str:
    """Component of the density: given, or the only one with maxima."""
    available = sorted({item.component for item in result.maxima})

    if not available:
        raise ValueError("the result has no maxima")

    if component is None:
        if len(available) > 1:
            raise ValueError("choose the component: " + ", ".join(available))

        return available[0]

    if component not in available:
        raise ValueError(
            f"no maxima of component {component!r} ({', '.join(available)})"
        )

    return component


def _coordinate(maximum: SpectralMaximum, coordinate: DensityCoordinate
                ) -> float:
    """Value of a maximum in the density coordinate."""
    if coordinate is DensityCoordinate.WAVENUMBER:
        return maximum.wavenumber

    if maximum.wavenumber <= 0.0:
        raise ValueError("logarithmic coordinates need positive wavenumbers")

    if coordinate is DensityCoordinate.LOG_SLOWNESS:
        return float(np.log(maximum.slowness))

    return float(np.log(maximum.velocity))


def _reference_powers(maxima: Iterable) -> dict:
    """Largest power of each block, frequency and component."""
    reference: dict = {}

    for item in maxima:
        key = (item.component, item.frequency, item.block)
        reference[key] = max(reference.get(key, 0.0), item.power)

    return reference


def _weights(group: list, weighting: DensityWeighting, reference: dict
             ) -> np.ndarray:
    """Normalized weights of the maxima of one frequency."""
    if weighting is DensityWeighting.UNIFORM:
        weights = np.ones(len(group))
    elif weighting is DensityWeighting.POWER:
        weights = np.array([max(item.power, 0.0) for item in group])
    else:
        weights = np.array([
            max(item.power, 0.0)
            / reference[(item.component, item.frequency, item.block)]
            if reference[(item.component, item.frequency, item.block)] > 0
            else 0.0
            for item in group
        ])

    total = np.nansum(weights)

    if not np.isfinite(total) or total <= 0.0:
        return np.full(len(group), 1.0 / len(group))

    return np.nan_to_num(weights) / total


def _grid_spacing(
    result: BeamformingResult,
    frequency: float,
    samples: np.ndarray,
    coordinate: DensityCoordinate,
) -> float:
    """Median local spacing of the search grid at the samples."""
    k = result.config.grid.wavenumbers(frequency)

    if coordinate is DensityCoordinate.WAVENUMBER:
        grid = k
    else:
        k = k[k > 0.0]
        grid = np.log(k / (2.0 * np.pi * frequency))

        if coordinate is DensityCoordinate.LOG_VELOCITY:
            grid = -grid

    grid = np.unique(grid)
    spacing = np.gradient(grid)
    index = np.argmin(np.abs(grid[None, :] - samples[:, None]), axis=1)
    return float(np.median(spacing[index]))


def _silverman(values: np.ndarray, weights: np.ndarray) -> float:
    """Weighted Silverman bandwidth, robust to outliers."""
    centre = np.sum(weights * values)
    deviation = np.sqrt(np.sum(weights * (values - centre) ** 2))
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    q25, q75 = np.interp([0.25, 0.75], cumulative / cumulative[-1],
                         values[order])
    scales = [s for s in (deviation, (q75 - q25) / 1.349) if s > 0.0]

    if not scales:
        return 0.02 * max(1.0, abs(centre))

    effective = 1.0 / np.sum(weights ** 2)
    return float(0.9 * min(scales) * effective ** -0.2)


def _mixture(grid: np.ndarray, samples: np.ndarray, weights: np.ndarray,
             width: float) -> np.ndarray:
    """Weighted mixture of Gaussian kernels on a grid."""
    offset = (grid[:, None] - samples[None, :]) / width
    kernels = np.exp(-0.5 * offset ** 2) / (np.sqrt(2.0 * np.pi) * width)
    return kernels @ weights
