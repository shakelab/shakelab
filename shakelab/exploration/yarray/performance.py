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
Statistical performance of the estimators (Monte Carlo).

Controlled experiments with one scalar plane wave in spatially white,
circular complex Gaussian noise. For each trial, independent Fourier
realizations::

    x_m = s_m a(k0) + n_m,    s_m ~ CN(0, P_s),    n_m ~ CN(0, P_n I),

are averaged into a sample cross-spectral matrix
(:func:`shakelab.signals.crossspectral.cross_spectral_matrix`), the
estimator (:func:`.estimators.beam_power`, Bartlett, Capon or MUSIC)
is evaluated on the search grid and its maximum ``k_hat`` is
classified from the geometric response (:mod:`.response`) through the
error ``dk = k_hat - k0``:

``BOUNDARY``
    the maximum lies on the radial boundary of the search grid;
``LOCAL``
    ``dk`` lies in the main lobe (centred on the true wavevector);
``ALIAS``
    ``dk`` lies in a remote region of the response above the
    ambiguity level (an alias family);
``OTHER``
    none of the above (noise-driven remote maxima).

``LOCAL`` means geometrically local, not accurate to a tolerance; a
remote noise maximum is an alias only if it belongs to a geometric
ambiguity family. The SNR is per station, ``P_s / P_n``; the
realizations are independent, so their effective number equals their
number.

:func:`evaluate_performance` runs one experiment,
:func:`performance_sweep` a grid of SNR and realizations. The
experiments are reproducible through ``MonteCarloConfig.seed``. The
search grid must be of wavenumbers (``quantity="wavenumber"``): the
experiment does not depend on the frequency, which enters only the
conversion to velocity (:meth:`PerformanceResult.velocity_statistics`).
Correlated windows, coloured noise, several sources, sampled
coordinate errors and RTBF are not modelled.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from enum import IntEnum
from statistics import NormalDist
from typing import Any, Optional

import numpy as np
from scipy.ndimage import label

from shakelab.exploration.yarray.config import GridQuantity, SearchGrid
from shakelab.exploration.yarray.estimators import (
    BeamformingMethod,
    beam_power,
    steering_vectors,
)
from shakelab.exploration.yarray.response import (
    ArrayResponse,
    ResponseAnalysis,
    analyze_response,
)
from shakelab.signals.crossspectral import cross_spectral_matrix

__all__ = [
    "AliasFamily",
    "ClassificationConfig",
    "EstimatorConfig",
    "LocalStatistics",
    "MonteCarloConfig",
    "OutcomeClass",
    "OutcomeSummary",
    "PerformanceConfig",
    "PerformanceResult",
    "PerformanceSamples",
    "PerformanceSweepResult",
    "PlaneWaveSource",
    "VelocityStatistics",
    "evaluate_performance",
    "performance_sweep",
    "theoretical_cross_spectrum",
]

_NO_ALIAS = -1


def _readonly(value: Any, dtype: Any = float) -> np.ndarray:
    """Return a read-only copy of an array."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _integer(value: Any, name: str, minimum: int) -> int:
    """Return an integer not smaller than ``minimum``."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")

    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")

    return int(value)


def _finite(value: Any, name: str) -> float:
    """Return a finite real number."""
    value = float(value)

    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")

    return value


def _method(value: Any) -> BeamformingMethod:
    """Return a scalar estimator from a member or its value."""
    if not isinstance(value, BeamformingMethod):
        for member in BeamformingMethod:
            if member.value == value:
                value = member
                break
        else:
            raise ValueError(f"unknown method {value!r}; the experiments "
                             "support bartlett, capon and music")

    return value


class _Replace:
    """Mixin returning modified copies of a frozen dataclass."""

    def replace(self, **changes: Any) -> Any:
        """
        Return a copy with some fields changed.

        Parameters
        ----------
        **changes
            New field values.

        Returns
        -------
        object
            Validated copy.
        """
        return dataclasses.replace(self, **changes)


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------


class OutcomeClass(IntEnum):
    """Mutually exclusive outcomes of a trial."""

    LOCAL = 0
    ALIAS = 1
    OTHER = 2
    BOUNDARY = 3


@dataclass(frozen=True)
class PlaneWaveSource:
    """
    Scalar stochastic plane wave.

    Parameters
    ----------
    wavenumber : float
        Angular wavenumber (non-negative).
    azimuth : float
        Propagation azimuth in degrees clockwise from north.
    signal_power : float, optional
        Expected Fourier power ``E[|s|**2]`` at each station. Default 1.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    wavenumber: float
    azimuth: float
    signal_power: float = 1.0

    def __post_init__(self) -> None:
        wavenumber = _finite(self.wavenumber, "wavenumber")
        power = _finite(self.signal_power, "signal_power")

        if wavenumber < 0.0:
            raise ValueError("wavenumber cannot be negative")

        if power <= 0.0:
            raise ValueError("signal_power must be positive")

        object.__setattr__(self, "wavenumber", wavenumber)
        object.__setattr__(self, "signal_power", power)
        object.__setattr__(
            self, "azimuth", float(np.mod(_finite(self.azimuth, "azimuth"),
                                          360.0))
        )

    @classmethod
    def from_velocity(
        cls,
        velocity: float,
        frequency: float,
        azimuth: float,
        signal_power: float = 1.0,
    ) -> PlaneWaveSource:
        """
        Return the plane wave of a phase velocity at a frequency.

        Parameters
        ----------
        velocity : float
            Phase velocity.
        frequency : float
            Frequency in Hz.
        azimuth : float
            Propagation azimuth in degrees clockwise from north.
        signal_power : float, optional
            Fourier power at each station. Default is 1.

        Returns
        -------
        PlaneWaveSource
            Source with wavenumber ``2 pi f / c``.
        """
        velocity = _finite(velocity, "velocity")
        frequency = _finite(frequency, "frequency")

        if velocity <= 0.0 or frequency <= 0.0:
            raise ValueError("velocity and frequency must be positive")

        return cls(2.0 * np.pi * frequency / velocity, azimuth, signal_power)

    @property
    def kx(self) -> float:
        """East wavenumber."""
        return float(self.wavenumber * np.sin(np.radians(self.azimuth)))

    @property
    def ky(self) -> float:
        """North wavenumber."""
        return float(self.wavenumber * np.cos(np.radians(self.azimuth)))

    @property
    def wavevector(self) -> np.ndarray:
        """Wavevector ``[kx, ky]``."""
        return np.array([self.kx, self.ky])


@dataclass(frozen=True)
class EstimatorConfig(_Replace):
    """
    Estimator of the experiments.

    Parameters
    ----------
    method : BeamformingMethod or str, optional
        ``"bartlett"`` (default), ``"capon"`` or ``"music"``.
    diagonal_loading : float, optional
        Relative diagonal loading of Capon. Default is 0.001.
    music_sources : int, optional
        Signal eigenvectors excluded by MUSIC. Default is 1.
    strict : bool, optional
        Require more realizations than stations for Capon and MUSIC.
        Default is True.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    method: Any = BeamformingMethod.BARTLETT
    diagonal_loading: float = 0.001
    music_sources: int = 1
    strict: bool = True

    def __post_init__(self) -> None:
        loading = _finite(self.diagonal_loading, "diagonal_loading")

        if loading < 0.0:
            raise ValueError("diagonal_loading cannot be negative")

        if not isinstance(self.strict, bool):
            raise TypeError("strict must be a bool")

        object.__setattr__(self, "method", _method(self.method))
        object.__setattr__(self, "diagonal_loading", loading)
        object.__setattr__(self, "music_sources",
                           _integer(self.music_sources, "music_sources", 1))


@dataclass(frozen=True)
class MonteCarloConfig(_Replace):
    """
    Monte Carlo sampling.

    Parameters
    ----------
    trials : int, optional
        Number of trials (at least two). Default is 1000.
    realizations : int, optional
        Independent Fourier realizations of each sample matrix.
        Default is 40.
    seed : int, optional
        Seed of the random generator; None for a non-reproducible run.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    trials: int = 1000
    realizations: int = 40
    seed: Optional[int] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "trials",
                           _integer(self.trials, "trials", 2))
        object.__setattr__(self, "realizations",
                           _integer(self.realizations, "realizations", 1))

        if self.seed is not None:
            object.__setattr__(self, "seed", _integer(self.seed, "seed", 0))


@dataclass(frozen=True)
class ClassificationConfig(_Replace):
    """
    Geometric classification of the trials.

    Parameters
    ----------
    resolution_level : float, optional
        Response level of the LOCAL region. Default is 0.5.
    ambiguity_level : float, optional
        Response level of the ALIAS regions. Default is 0.5.
    boundary_margin : int, optional
        Radial samples next to each radial boundary of the search grid
        classified as BOUNDARY (the azimuth is periodic). Default 0.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    resolution_level: float = 0.5
    ambiguity_level: float = 0.5
    boundary_margin: int = 0

    def __post_init__(self) -> None:
        for name in ("resolution_level", "ambiguity_level"):
            value = _finite(getattr(self, name), name)

            if not 0.0 < value < 1.0:
                raise ValueError(f"{name} must be between 0 and 1 "
                                 "(excluded)")

            object.__setattr__(self, name, value)

        object.__setattr__(
            self, "boundary_margin",
            _integer(self.boundary_margin, "boundary_margin", 0),
        )


@dataclass(frozen=True)
class PerformanceConfig(_Replace):
    """
    Configuration of a performance experiment.

    Parameters
    ----------
    snr_db : float
        SNR per station in decibels: the noise power is
        ``signal_power / 10**(snr_db / 10)``.
    grid : SearchGrid
        Search grid of wavenumbers (``quantity="wavenumber"``).
    estimator : EstimatorConfig, optional
        Estimator. Default: Bartlett.
    monte_carlo : MonteCarloConfig, optional
        Trials, realizations and seed.
    classification : ClassificationConfig, optional
        Levels of the classification.

    Raises
    ------
    ValueError
        If the grid is not of wavenumbers or the boundary margin
        leaves no interior samples.
    """

    snr_db: float
    grid: SearchGrid
    estimator: EstimatorConfig = EstimatorConfig()
    monte_carlo: MonteCarloConfig = MonteCarloConfig()
    classification: ClassificationConfig = ClassificationConfig()

    def __post_init__(self) -> None:
        object.__setattr__(self, "snr_db", _finite(self.snr_db, "snr_db"))

        if not isinstance(self.grid, SearchGrid):
            raise TypeError("grid must be a SearchGrid")

        if self.grid.quantity is not GridQuantity.WAVENUMBER:
            raise ValueError("the performance experiments need a grid of "
                             "wavenumbers")

        for name, kind in (
            ("estimator", EstimatorConfig),
            ("monte_carlo", MonteCarloConfig),
            ("classification", ClassificationConfig),
        ):
            if not isinstance(getattr(self, name), kind):
                raise TypeError(f"{name} must be a {kind.__name__}")

        if 2 * self.classification.boundary_margin >= self.grid.size:
            raise ValueError("boundary_margin leaves no interior radial "
                             "samples")

    @property
    def snr(self) -> float:
        """Linear SNR per station."""
        return float(10.0 ** (self.snr_db / 10.0))


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class AliasFamily:
    """
    Remote connected region of the response (geometric alias family).

    Attributes
    ----------
    index : int
        Index from 0, by increasing onset distance.
    onset_distance : float
        Smallest wavevector displacement in the region.
    peak_kx, peak_ky : float
        Displacement of the strongest grid point of the region.
    peak_power : float
        Largest response of the region.
    """

    index: int
    onset_distance: float
    peak_kx: float
    peak_ky: float
    peak_power: float

    @property
    def peak_distance(self) -> float:
        """Length of the peak displacement."""
        return float(np.hypot(self.peak_kx, self.peak_ky))

    @property
    def peak_azimuth(self) -> float:
        """Azimuth of the peak displacement, clockwise from north."""
        return float(np.mod(np.degrees(np.arctan2(self.peak_kx,
                                                   self.peak_ky)), 360.0))


@dataclass(frozen=True, eq=False)
class PerformanceSamples:
    """
    Estimates and classes of the trials.

    Attributes
    ----------
    kx, ky : numpy.ndarray
        Estimated wavevectors.
    power : numpy.ndarray
        Power of the maxima.
    outcome : numpy.ndarray
        :class:`OutcomeClass` values.
    alias_region : numpy.ndarray
        Index of the alias family for ALIAS trials, -1 otherwise.
    """

    kx: np.ndarray
    ky: np.ndarray
    power: np.ndarray
    outcome: np.ndarray
    alias_region: np.ndarray

    def __post_init__(self) -> None:
        if self.kx.ndim != 1:
            raise ValueError("the samples must be one-dimensional")

        for name in ("ky", "power", "outcome", "alias_region"):
            if getattr(self, name).shape != self.kx.shape:
                raise ValueError(f"{name} does not match kx")

    @property
    def size(self) -> int:
        """Number of trials."""
        return int(self.kx.size)

    @property
    def wavenumber(self) -> np.ndarray:
        """Estimated wavenumbers."""
        return np.hypot(self.kx, self.ky)

    @property
    def azimuth(self) -> np.ndarray:
        """Estimated propagation azimuths."""
        return np.mod(np.degrees(np.arctan2(self.kx, self.ky)), 360.0)

    def mask(self, outcome: Any) -> np.ndarray:
        """
        Return the trials of an outcome.

        Parameters
        ----------
        outcome : OutcomeClass or str
            Outcome (e.g. ``"local"``).

        Returns
        -------
        numpy.ndarray
            Boolean mask.
        """
        return self.outcome == int(_outcome(outcome))


def _outcome(value: Any) -> OutcomeClass:
    """Return an outcome from a member or its (case-insensitive) name."""
    if isinstance(value, OutcomeClass):
        return value

    if isinstance(value, str) and value.upper() in OutcomeClass.__members__:
        return OutcomeClass[value.upper()]

    raise ValueError(f"unknown outcome {value!r}")


@dataclass(frozen=True)
class OutcomeSummary:
    """
    Counts and probabilities of the outcomes.

    Attributes
    ----------
    trials : int
        Number of trials.
    local_count, alias_count, other_count, boundary_count : int
        Trials of each outcome.
    """

    trials: int
    local_count: int
    alias_count: int
    other_count: int
    boundary_count: int

    def __post_init__(self) -> None:
        total = (self.local_count + self.alias_count + self.other_count
                 + self.boundary_count)

        if total != self.trials:
            raise ValueError("the counts do not sum to the trials")

    def count(self, outcome: Any) -> int:
        """Return the trials of an outcome (OutcomeClass or name)."""
        outcome = _outcome(outcome)
        return int(getattr(self, f"{outcome.name.lower()}_count"))

    def probability(self, outcome: Any) -> float:
        """Return the probability of an outcome."""
        return self.count(outcome) / self.trials

    def standard_error(self, outcome: Any) -> float:
        """Return the binomial standard error of a probability."""
        probability = self.probability(outcome)
        return float(np.sqrt(probability * (1.0 - probability)
                             / self.trials))

    def wilson_interval(self, outcome: Any,
                        confidence: float = 0.95) -> tuple[float, float]:
        """
        Return the Wilson score interval of a probability.

        Unlike the symmetric normal interval, it behaves well near
        zero and one.

        Parameters
        ----------
        outcome : OutcomeClass or str
            Outcome.
        confidence : float, optional
            Confidence level. Default is 0.95.

        Returns
        -------
        tuple of float
            Lower and upper bounds.
        """
        confidence = _finite(confidence, "confidence")

        if not 0.0 < confidence < 1.0:
            raise ValueError("confidence must be between 0 and 1")

        count = self.count(outcome)
        n = float(self.trials)
        p = count / n
        z = NormalDist().inv_cdf(0.5 + 0.5 * confidence)
        denominator = 1.0 + z ** 2 / n
        center = (p + z ** 2 / (2.0 * n)) / denominator
        half = z / denominator * np.sqrt(p * (1.0 - p) / n
                                         + z ** 2 / (4.0 * n ** 2))

        # The bounds are exactly 0 and 1 at the extremes (rounding).
        low = 0.0 if count == 0 else max(0.0, center - half)
        high = 1.0 if count == self.trials else min(1.0, center + half)

        return float(low), float(high)

    @property
    def local_probability(self) -> float:
        """Probability of LOCAL estimates."""
        return self.probability(OutcomeClass.LOCAL)

    @property
    def alias_probability(self) -> float:
        """Probability of ALIAS estimates."""
        return self.probability(OutcomeClass.ALIAS)

    @property
    def other_probability(self) -> float:
        """Probability of OTHER estimates."""
        return self.probability(OutcomeClass.OTHER)

    @property
    def boundary_probability(self) -> float:
        """Probability of BOUNDARY estimates."""
        return self.probability(OutcomeClass.BOUNDARY)


@dataclass(frozen=True, eq=False)
class LocalStatistics:
    """
    Wavevector errors of the LOCAL estimates.

    Attributes
    ----------
    count : int
        Number of LOCAL estimates.
    bias : numpy.ndarray
        Mean Cartesian error ``[dkx, dky]``.
    covariance : numpy.ndarray
        Covariance ``(2, 2)`` of the errors.
    rmse : float
        Root mean square error length.
    radial_bias, radial_standard_deviation : float
        Error along the true wavevector.
    transverse_bias, transverse_standard_deviation : float
        Error across the true wavevector.
    """

    count: int
    bias: np.ndarray
    covariance: np.ndarray
    rmse: float
    radial_bias: float
    radial_standard_deviation: float
    transverse_bias: float
    transverse_standard_deviation: float

    @property
    def standard_deviation_kx(self) -> float:
        """Standard deviation of the east wavenumber."""
        return float(np.sqrt(self.covariance[0, 0]))

    @property
    def standard_deviation_ky(self) -> float:
        """Standard deviation of the north wavenumber."""
        return float(np.sqrt(self.covariance[1, 1]))


@dataclass(frozen=True)
class VelocityStatistics:
    """
    Phase velocities of the LOCAL estimates at a frequency.

    Each estimate is converted exactly, ``c = 2 pi f / k`` (no
    linearization).

    Attributes
    ----------
    frequency : float
        Frequency in Hz.
    count : int
        Number of estimates.
    true_velocity : float
        Velocity of the source.
    bias, standard_deviation, rmse : float
        Error statistics.
    median, lower_quantile, upper_quantile : float
        Median and central quantile interval.
    minimum, maximum : float
        Extremes.
    """

    frequency: float
    count: int
    true_velocity: float
    bias: float
    standard_deviation: float
    rmse: float
    median: float
    lower_quantile: float
    upper_quantile: float
    minimum: float
    maximum: float

    @property
    def relative_rmse(self) -> float:
        """RMSE relative to the true velocity."""
        return self.rmse / self.true_velocity


@dataclass(frozen=True, eq=False)
class PerformanceResult:
    """
    Result of a performance experiment.

    Attributes
    ----------
    source : PlaneWaveSource
        True plane wave.
    config : PerformanceConfig
        Configuration.
    analysis : ResponseAnalysis
        Response analysis used to classify the trials.
    alias_families : tuple of AliasFamily
        Alias families of the response.
    samples : PerformanceSamples
        Estimates of the trials.
    outcome : OutcomeSummary
        Outcome counts.
    local_statistics : LocalStatistics or None
        Errors of the LOCAL estimates (None if there are none).
    global_bias, global_covariance : numpy.ndarray
        Error statistics of all the estimates.
    global_rmse : float
        RMSE of all the estimates.
    alias_counts : numpy.ndarray
        Trials of each alias family.
    """

    source: PlaneWaveSource
    config: PerformanceConfig
    analysis: ResponseAnalysis
    alias_families: tuple
    samples: PerformanceSamples
    outcome: OutcomeSummary
    local_statistics: Optional[LocalStatistics]
    global_bias: np.ndarray
    global_covariance: np.ndarray
    global_rmse: float
    alias_counts: np.ndarray

    @property
    def response(self) -> ArrayResponse:
        """Geometric response."""
        return self.analysis.response

    @property
    def noise_power(self) -> float:
        """Noise power per station."""
        return float(self.source.signal_power / self.config.snr)

    @property
    def nearest_search_wavevector(self) -> np.ndarray:
        """Grid wavevector nearest to the source."""
        kx, ky = _search_wavevectors(self.config.grid)
        distance = (kx - self.source.kx) ** 2 + (ky - self.source.ky) ** 2
        row, column = np.unravel_index(int(np.argmin(distance)),
                                       distance.shape)

        return np.array([kx[row, column], ky[row, column]])

    @property
    def grid_quantization_error(self) -> float:
        """
        Distance between the source and the nearest grid wavevector.

        A floor of the error set by the grid, not by the estimator.
        """
        return float(np.linalg.norm(self.nearest_search_wavevector
                                    - self.source.wavevector))

    @property
    def reached_alias_families(self) -> tuple:
        """Alias families selected by at least one trial."""
        return tuple(
            family for family, count in zip(self.alias_families,
                                             self.alias_counts)
            if count > 0
        )

    @property
    def alias_probabilities(self) -> np.ndarray:
        """Probability of each alias family."""
        return self.alias_counts / self.config.monte_carlo.trials

    @property
    def local_wavenumbers(self) -> np.ndarray:
        """Wavenumbers of the LOCAL estimates."""
        return self.samples.wavenumber[self.samples.mask(OutcomeClass.LOCAL)]

    @property
    def effective_realizations(self) -> float:
        """Effective realizations (independent: their number)."""
        return float(self.config.monte_carlo.realizations)

    def velocity_statistics(
        self,
        frequency: float,
        central_probability: float = 0.95,
    ) -> Optional[VelocityStatistics]:
        """
        Return the phase velocities of the LOCAL estimates.

        Parameters
        ----------
        frequency : float
            Frequency in Hz of ``c = 2 pi f / k``.
        central_probability : float, optional
            Probability of the quantile interval. Default is 0.95.

        Returns
        -------
        VelocityStatistics or None
            Statistics, or None without LOCAL estimates.

        Raises
        ------
        ValueError
            If a value is out of range or the source wavenumber is
            zero.
        """
        frequency = _finite(frequency, "frequency")
        central_probability = _finite(central_probability,
                                      "central_probability")

        if frequency <= 0.0:
            raise ValueError("frequency must be positive")

        if not 0.0 < central_probability < 1.0:
            raise ValueError("central_probability must be between 0 and 1")

        if self.source.wavenumber <= 0.0:
            raise ValueError("the velocity needs a positive source "
                             "wavenumber")

        k = self.local_wavenumbers
        k = k[np.isfinite(k) & (k > 0.0)]

        if k.size == 0:
            return None

        velocity = 2.0 * np.pi * frequency / k
        truth = 2.0 * np.pi * frequency / self.source.wavenumber
        error = velocity - truth
        alpha = 0.5 * (1.0 - central_probability)
        lower, upper = np.quantile(velocity, [alpha, 1.0 - alpha])

        return VelocityStatistics(
            frequency=frequency,
            count=int(velocity.size),
            true_velocity=float(truth),
            bias=float(np.mean(error)),
            standard_deviation=(
                float(np.std(velocity, ddof=1)) if velocity.size > 1
                else float("nan")
            ),
            rmse=float(np.sqrt(np.mean(error ** 2))),
            median=float(np.median(velocity)),
            lower_quantile=float(lower),
            upper_quantile=float(upper),
            minimum=float(np.min(velocity)),
            maximum=float(np.max(velocity)),
        )

    def info(self) -> None:
        """Print a summary of the experiment."""
        config = self.config
        outcome = self.outcome
        low, high = outcome.wilson_interval(OutcomeClass.LOCAL)
        print(
            f"PerformanceResult: {config.estimator.method.value}, "
            f"SNR {config.snr_db:g} dB, {config.monte_carlo.trials} "
            f"trials of {config.monte_carlo.realizations} realizations, "
            f"source k={self.source.wavenumber:.6g} azimuth "
            f"{self.source.azimuth:g}"
        )
        print(
            f"  local {outcome.local_probability:.4f} "
            f"[{low:.4f}, {high:.4f}], alias "
            f"{outcome.alias_probability:.4f}, other "
            f"{outcome.other_probability:.4f}, boundary "
            f"{outcome.boundary_probability:.4f}"
        )

        if self.local_statistics is None:
            print("  no local estimates")
        else:
            local = self.local_statistics
            print(
                f"  local errors: rmse {local.rmse:.6g}, radial "
                f"{local.radial_bias:.3g} +/- "
                f"{local.radial_standard_deviation:.3g}, transverse "
                f"{local.transverse_bias:.3g} +/- "
                f"{local.transverse_standard_deviation:.3g}"
            )

        print(
            f"  grid quantization {self.grid_quantization_error:.6g}, "
            f"global rmse {self.global_rmse:.6g}, "
            f"{len(self.alias_families)} alias families "
            f"({len(self.reached_alias_families)} reached)"
        )


@dataclass(frozen=True, eq=False)
class PerformanceSweepResult:
    """
    Experiments on a grid of SNR and realizations.

    Attributes
    ----------
    snr_db : numpy.ndarray
        SNR values in decibels (rows).
    realizations : numpy.ndarray
        Realization counts (columns).
    results : tuple of PerformanceResult
        Results in row-major order.
    """

    snr_db: np.ndarray
    realizations: np.ndarray
    results: tuple

    def __post_init__(self) -> None:
        if len(self.results) != self.snr_db.size * self.realizations.size:
            raise ValueError("the results do not match the sweep axes")

    @property
    def shape(self) -> tuple[int, int]:
        """``(n_snr, n_realizations)``."""
        return int(self.snr_db.size), int(self.realizations.size)

    def result(self, snr_index: int, realizations_index: int
               ) -> PerformanceResult:
        """Return the result of a sweep point."""
        return self.results[snr_index * self.realizations.size
                            + realizations_index]

    def _matrix(self, getter: Any) -> np.ndarray:
        """Return a quantity of every point as a matrix."""
        return np.array([getter(result) for result in self.results],
                        dtype=float).reshape(self.shape)

    @property
    def local_probability(self) -> np.ndarray:
        """LOCAL probabilities."""
        return self._matrix(lambda r: r.outcome.local_probability)

    @property
    def alias_probability(self) -> np.ndarray:
        """ALIAS probabilities."""
        return self._matrix(lambda r: r.outcome.alias_probability)

    @property
    def other_probability(self) -> np.ndarray:
        """OTHER probabilities."""
        return self._matrix(lambda r: r.outcome.other_probability)

    @property
    def boundary_probability(self) -> np.ndarray:
        """BOUNDARY probabilities."""
        return self._matrix(lambda r: r.outcome.boundary_probability)

    def _local(self, name: str) -> np.ndarray:
        """Return a LOCAL statistic of every point (NaN if none)."""
        return self._matrix(
            lambda r: float("nan") if r.local_statistics is None
            else getattr(r.local_statistics, name)
        )

    @property
    def local_rmse(self) -> np.ndarray:
        """LOCAL wavevector RMSE."""
        return self._local("rmse")

    @property
    def radial_standard_deviation(self) -> np.ndarray:
        """LOCAL radial standard deviations."""
        return self._local("radial_standard_deviation")

    @property
    def transverse_standard_deviation(self) -> np.ndarray:
        """LOCAL transverse standard deviations."""
        return self._local("transverse_standard_deviation")

    def info(self) -> None:
        """Print a summary of the sweep."""
        local = self.local_probability
        first = self.results[0]
        print(
            f"PerformanceSweepResult: {first.config.estimator.method.value}"
            f", {self.shape[0]} SNR ({self.snr_db.min():g} to "
            f"{self.snr_db.max():g} dB) x {self.shape[1]} realizations "
            f"({self.realizations.min()} to {self.realizations.max()}), "
            f"local probability {local.min():.3f} to {local.max():.3f}"
        )


# ----------------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class _Context:
    """Quantities shared by the experiments on one response."""

    analysis: ResponseAnalysis
    lobe: np.ndarray
    labels: np.ndarray
    families: tuple
    levels: tuple


def _search_wavevectors(grid: SearchGrid) -> tuple[np.ndarray, np.ndarray]:
    """Return the wavevectors of a wavenumber grid."""
    return grid.wavevectors(1.0)


def _remote_regions(response: ArrayResponse,
                    level: float) -> tuple[np.ndarray, tuple]:
    """Return the labels (-1 outside) and the alias families."""
    power = response.power
    components, count = label(power >= level,
                              structure=np.ones((3, 3), dtype=int))
    center = response.grid.center_index
    central = int(components[center, center])
    kx, ky = response.grid.mesh()
    radial = np.hypot(kx, ky)
    candidates = []

    for component in range(1, count + 1):
        if component == central:
            continue

        mask = components == component
        flat = np.flatnonzero(mask)
        winner = int(flat[int(np.argmax(power.ravel()[flat]))])
        row, column = np.unravel_index(winner, power.shape)
        candidates.append((float(np.min(radial[mask])), component,
                           float(kx[row, column]), float(ky[row, column]),
                           float(power[row, column])))

    candidates.sort(key=lambda item: item[0])
    labels = np.full(components.shape, _NO_ALIAS, dtype=np.int64)
    families = []

    for index, (onset, component, px, py, peak) in enumerate(candidates):
        labels[components == component] = index
        families.append(AliasFamily(index, onset, px, py, peak))

    return labels, tuple(families)


def _context(response: ArrayResponse,
             classification: ClassificationConfig) -> _Context:
    """Return the classification quantities of a response."""
    analysis = analyze_response(
        response,
        resolution_level=classification.resolution_level,
        ambiguity_level=classification.ambiguity_level,
    )
    labels, families = _remote_regions(response,
                                       classification.ambiguity_level)

    return _Context(
        analysis=analysis,
        lobe=np.asarray(analysis.main_lobe.mask, dtype=bool),
        labels=labels,
        families=families,
        levels=(classification.resolution_level,
                classification.ambiguity_level),
    )


def _check(response: ArrayResponse, source: PlaneWaveSource,
           config: PerformanceConfig) -> None:
    """Check the inputs of an experiment."""
    if not isinstance(response, ArrayResponse):
        raise TypeError("response must be an ArrayResponse")

    if not isinstance(source, PlaneWaveSource):
        raise TypeError("source must be a PlaneWaveSource")

    if not isinstance(config, PerformanceConfig):
        raise TypeError("config must be a PerformanceConfig")

    if not response.geometry.is_two_dimensional:
        raise ValueError("the performance experiments need a "
                         "two-dimensional array")

    grid = config.grid
    tolerance = 16.0 * np.finfo(float).eps * max(1.0, abs(grid.maximum))

    if not (grid.minimum - tolerance <= source.wavenumber
            <= grid.maximum + tolerance):
        raise ValueError("the source wavenumber is outside the search grid")

    if (
        config.estimator.method is BeamformingMethod.MUSIC
        and config.estimator.music_sources >= response.geometry.n_stations
    ):
        raise ValueError("music_sources must be smaller than the number of "
                         "stations")


def _complex_normal(rng: np.random.Generator, shape: tuple,
                    power: float) -> np.ndarray:
    """Return circular complex Gaussian samples of a given power."""
    scale = np.sqrt(power / 2.0)
    return scale * (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape))


def _source_steering(response: ArrayResponse,
                     source: PlaneWaveSource) -> np.ndarray:
    """Return the steering vector of the source."""
    return np.asarray(steering_vectors(
        response.geometry.centered_coordinates, np.asarray(source.kx),
        np.asarray(source.ky),
    ), dtype=complex).reshape(-1)


def theoretical_cross_spectrum(
    response: ArrayResponse,
    source: PlaneWaveSource,
    config: PerformanceConfig,
) -> np.ndarray:
    """
    Return the expected cross-spectral matrix of an experiment.

    ``R = P_s a(k0) a(k0)^H + P_n I`` (same convention as
    :func:`shakelab.signals.crossspectral.cross_spectral_matrix`).

    Parameters
    ----------
    response : ArrayResponse
        Geometric response (its centered coordinates are used).
    source : PlaneWaveSource
        Plane wave.
    config : PerformanceConfig
        Configuration (SNR).

    Returns
    -------
    numpy.ndarray
        Matrix ``(n, n)``.
    """
    _check(response, source, config)
    steering = _source_steering(response, source)
    matrix = (
        source.signal_power * np.outer(steering, steering.conj())
        + source.signal_power / config.snr * np.eye(steering.size)
    )

    return _readonly(matrix, complex)


def _error_statistics(kx: np.ndarray, ky: np.ndarray,
                      source: PlaneWaveSource) -> tuple:
    """Return Cartesian bias, covariance and RMSE of the errors."""
    error = np.column_stack((kx - source.kx, ky - source.ky))
    covariance = (
        np.cov(error, rowvar=False, ddof=1) if error.shape[0] > 1
        else np.full((2, 2), np.nan)
    )

    return (np.mean(error, axis=0), np.asarray(covariance),
            float(np.sqrt(np.mean(np.sum(error ** 2, axis=1)))))


def _local_statistics(samples: PerformanceSamples,
                      source: PlaneWaveSource) -> Optional[LocalStatistics]:
    """Return the error statistics of the LOCAL estimates."""
    mask = samples.mask(OutcomeClass.LOCAL)
    count = int(np.count_nonzero(mask))

    if count == 0:
        return None

    kx, ky = samples.kx[mask], samples.ky[mask]
    bias, covariance, rmse = _error_statistics(kx, ky, source)
    nan = float("nan")
    radial_bias = transverse_bias = radial_std = transverse_std = nan

    if source.wavenumber > 0.0:
        error = np.column_stack((kx - source.kx, ky - source.ky))
        radial = source.wavevector / source.wavenumber
        along = error @ radial
        across = error @ np.array([radial[1], -radial[0]])
        radial_bias = float(np.mean(along))
        transverse_bias = float(np.mean(across))

        if count > 1:
            radial_std = float(np.std(along, ddof=1))
            transverse_std = float(np.std(across, ddof=1))

    return LocalStatistics(
        count=count,
        bias=_readonly(bias),
        covariance=_readonly(covariance),
        rmse=rmse,
        radial_bias=radial_bias,
        radial_standard_deviation=radial_std,
        transverse_bias=transverse_bias,
        transverse_standard_deviation=transverse_std,
    )


def _evaluate(response: ArrayResponse, source: PlaneWaveSource,
              config: PerformanceConfig, context: _Context,
              search: tuple) -> PerformanceResult:
    """Run one experiment with precomputed quantities."""
    search_kx, search_ky, search_steering = search
    estimator = config.estimator
    monte_carlo = config.monte_carlo
    grid = config.grid
    margin = config.classification.boundary_margin
    steering = _source_steering(response, source)
    noise = source.signal_power / config.snr
    rng = np.random.default_rng(monte_carlo.seed)
    trials = monte_carlo.trials
    count = monte_carlo.realizations
    estimates = np.empty((trials, 3))
    outcomes = np.empty(trials, dtype=np.int64)
    regions = np.full(trials, _NO_ALIAS, dtype=np.int64)
    rgrid = response.grid

    for trial in range(trials):
        amplitude = _complex_normal(rng, (count, 1), source.signal_power)
        values = amplitude * steering + _complex_normal(
            rng, (count, steering.size), noise
        )
        power = beam_power(
            cross_spectral_matrix(values), search_steering, count,
            estimator.method, strict=estimator.strict,
            diagonal_loading=estimator.diagonal_loading,
            music_sources=estimator.music_sources,
        )
        safe = np.where(np.isfinite(power), power, -np.inf)
        row, column = np.unravel_index(int(np.argmax(safe)), power.shape)
        kx = float(search_kx[row, column])
        ky = float(search_ky[row, column])
        estimates[trial] = kx, ky, power[row, column]

        if column <= margin or column >= grid.size - 1 - margin:
            outcomes[trial] = OutcomeClass.BOUNDARY
            continue

        dx, dy = kx - source.kx, ky - source.ky

        if abs(dx) > rgrid.maximum or abs(dy) > rgrid.maximum:
            outcomes[trial] = OutcomeClass.OTHER
            continue

        c = int(np.clip(np.rint((dx + rgrid.maximum) / rgrid.step), 0,
                        rgrid.size - 1))
        r = int(np.clip(np.rint((dy + rgrid.maximum) / rgrid.step), 0,
                        rgrid.size - 1))

        if context.lobe[r, c]:
            outcomes[trial] = OutcomeClass.LOCAL
        elif context.labels[r, c] >= 0:
            outcomes[trial] = OutcomeClass.ALIAS
            regions[trial] = context.labels[r, c]
        else:
            outcomes[trial] = OutcomeClass.OTHER

    samples = PerformanceSamples(
        kx=_readonly(estimates[:, 0]),
        ky=_readonly(estimates[:, 1]),
        power=_readonly(estimates[:, 2]),
        outcome=_readonly(outcomes, np.int64),
        alias_region=_readonly(regions, np.int64),
    )
    alias_counts = np.bincount(regions[regions >= 0],
                               minlength=len(context.families))
    bias, covariance, rmse = _error_statistics(samples.kx, samples.ky,
                                               source)
    summary = OutcomeSummary(
        trials=trials,
        **{
            f"{member.name.lower()}_count":
                int(np.count_nonzero(outcomes == member))
            for member in OutcomeClass
        },
    )

    return PerformanceResult(
        source=source,
        config=config,
        analysis=context.analysis,
        alias_families=context.families,
        samples=samples,
        outcome=summary,
        local_statistics=_local_statistics(samples, source),
        global_bias=_readonly(bias),
        global_covariance=_readonly(covariance),
        global_rmse=rmse,
        alias_counts=_readonly(alias_counts, np.int64),
    )


def _search(response: ArrayResponse, grid: SearchGrid) -> tuple:
    """Return the search wavevectors and their steering vectors."""
    kx, ky = _search_wavevectors(grid)
    steering = steering_vectors(response.geometry.centered_coordinates,
                                kx, ky)

    return kx, ky, steering


def evaluate_performance(
    response: ArrayResponse,
    source: PlaneWaveSource,
    config: PerformanceConfig,
) -> PerformanceResult:
    """
    Evaluate an estimator by Monte Carlo simulation.

    Parameters
    ----------
    response : ArrayResponse
        Geometric response of a two-dimensional array; its nominal
        coordinates generate the wave and steer the estimator.
    source : PlaneWaveSource
        True plane wave (inside the radial range of the search grid).
    config : PerformanceConfig
        SNR, estimator, search grid, sampling and classification.

    Returns
    -------
    PerformanceResult
        Estimates, outcome probabilities and error statistics.

    Raises
    ------
    ValueError
        If the array is linear, the source is outside the grid or the
        estimator needs more realizations (Capon and MUSIC with
        ``strict``).

    Notes
    -----
    The coordinate covariances of ``response`` are not sampled.
    """
    _check(response, source, config)

    return _evaluate(response, source, config,
                     _context(response, config.classification),
                     _search(response, config.grid))


def performance_sweep(
    response: ArrayResponse,
    source: PlaneWaveSource,
    config: PerformanceConfig,
    snr_db: Any,
    realizations: Any,
    seed_stride: int = 1000003,
) -> PerformanceSweepResult:
    """
    Evaluate an estimator on a grid of SNR and realizations.

    ``config`` is a template: SNR, realizations and seed change at
    each point, everything else is kept.

    Parameters
    ----------
    response : ArrayResponse
        Geometric response.
    source : PlaneWaveSource
        True plane wave.
    config : PerformanceConfig
        Template configuration.
    snr_db : array-like
        SNR values in decibels.
    realizations : array-like
        Realization counts (positive integers).
    seed_stride : int, optional
        Seed increment between points (seed of point ``i`` is
        ``seed + i * seed_stride``, with seed 0 if None). Default is
        1000003.

    Returns
    -------
    PerformanceSweepResult
        Results of every point.
    """
    snr = np.asarray(snr_db, dtype=float)
    raw = np.asarray(realizations)

    if snr.ndim != 1 or snr.size == 0 or not np.all(np.isfinite(snr)):
        raise ValueError("snr_db must be a sequence of finite values")

    if (
        raw.ndim != 1 or raw.size == 0
        or not np.all(np.equal(raw, np.floor(raw))) or np.any(raw < 1)
    ):
        raise ValueError("realizations must be a sequence of positive "
                         "integers")

    counts = raw.astype(np.int64)
    seed_stride = _integer(seed_stride, "seed_stride", 1)
    base = config.monte_carlo.seed or 0
    _check(response, source, config)
    context = _context(response, config.classification)
    search = _search(response, config.grid)
    results = []

    for value in snr:
        for count in counts:
            point = config.replace(
                snr_db=float(value),
                monte_carlo=config.monte_carlo.replace(
                    realizations=int(count),
                    seed=base + len(results) * seed_stride,
                ),
            )
            _check(response, source, point)
            results.append(_evaluate(response, source, point, context,
                                     search))

    return PerformanceSweepResult(snr_db=_readonly(snr),
                                  realizations=_readonly(counts, np.int64),
                                  results=tuple(results))
