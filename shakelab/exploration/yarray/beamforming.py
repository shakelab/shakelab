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
Frequency-wavenumber beamforming for seismic arrays.

The functions in this module operate on immutable ``ArrayData``
snapshots. File access, live buffering, waveform synchronization,
and sensor-orientation correction are handled outside the numerical
processing layer.

The reference coordinate system is:

    x = East
    y = North
    z = Up

Propagation azimuths are measured clockwise from North.

Reusable windowing, Fourier-transform, and cross-spectral-matrix
operations are implemented in :mod:`.spectral`. This module contains
the beamforming-specific orchestration layer: component preparation,
search-grid adaptation, Rayleigh joint RTBF, uncertainty propagation,
and spectral-maximum extraction. Shared scalar steering and estimator
kernels are implemented in :mod:`.estimators`.

Elementary windows may use either fixed duration or a fixed number of
cycles. Fourier realizations can be obtained by direct transforms at
arbitrary frequencies or, for fixed windows, by a shared RFFT with
nearest-bin or complex-interpolated frequency selection. Consecutive
elementary windows are grouped into covariance blocks by the spectral
layer, allowing both statistically stable matrix estimates and
time-dependent array-processing results.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import maximum_filter

from .estimators import (
    BeamformingMethod,
    bartlett_power_statistics,
    beam_power,
    capon_power_statistics,
    hermitian,
    regularized_matrix,
    steering_vectors,
)
from .data import EAST, NORTH, VERTICAL
from .spectral import (
    FFTSelection,
    FrequencyConfig,
    FrequencyWindows,
    SpectralComponents,
    SpectralConfig,
    SpectralMatrixBlocks,
    TransformMethod,
    WindowConfig,
    WindowType,
    cross_spectrum,
    effective_realizations as spectral_effective_realizations,
    matrix_blocks as build_matrix_blocks,
    spectral_slices,
)

if TYPE_CHECKING:
    from .data import ArrayData


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


class ComponentMode(Enum):
    """Component combinations supported by the beamformer."""

    VERTICAL = "vertical"
    RADIAL = "radial"
    TRANSVERSE = "transverse"
    THREE_COMPONENT_SEPARATE = "three_component_separate"
    RAYLEIGH_JOINT = "rayleigh_joint"


class GridQuantity(Enum):
    """Physical quantity sampled along the radial search axis."""

    VELOCITY = "velocity"
    SLOWNESS = "slowness"
    WAVENUMBER = "wavenumber"


class GridScale(Enum):
    """Sampling scale used along the radial search axis."""

    LINEAR = "linear"
    LOG = "log"


class OutputMode(Enum):
    """Amount of information retained after processing."""

    MAXIMA = "maxima"
    SPECTRA = "spectra"
    BOTH = "both"


class StatisticalUncertainty(Enum):
    """Statistical uncertainty treatment for stacked windows."""

    NONE = "none"
    EMPIRICAL = "empirical"


class GeometricUncertainty(Enum):
    """Treatment of station-coordinate uncertainty."""

    NONE = "none"
    EXPECTED_COHERENCE = "expected_coherence"


@dataclass(frozen=True)
class SearchGrid:
    """
    Polar search grid for frequency-wavenumber analysis.

    Parameters
    ----------
    minimum
        Minimum velocity, slowness, or wavenumber.
    maximum
        Maximum velocity, slowness, or wavenumber.
    size
        Number of samples along the radial axis.
    quantity
        Physical quantity represented by the radial values.
    scale
        Linear or logarithmic sampling.
    azimuth_size
        Number of uniformly spaced azimuths between 0 and 360
        degrees.
    azimuths
        Optional explicit azimuth array in degrees clockwise from
        North. If provided, ``azimuth_size`` is ignored.

    Notes
    -----
    Velocity is expressed in distance units per second, slowness in
    seconds per distance unit, and wavenumber in radians per distance
    unit. Station coordinates must use the corresponding distance
    unit.
    """

    minimum: float
    maximum: float
    size: int
    quantity: GridQuantity = GridQuantity.VELOCITY
    scale: GridScale = GridScale.LOG
    azimuth_size: int = 180
    azimuths: ArrayLike | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.quantity, GridQuantity):
            raise TypeError(
                "quantity must be a GridQuantity value."
            )

        if not isinstance(self.scale, GridScale):
            raise TypeError("scale must be a GridScale value.")

        if not np.isfinite(self.minimum):
            raise ValueError("Grid minimum must be finite.")

        if not np.isfinite(self.maximum):
            raise ValueError("Grid maximum must be finite.")

        if self.minimum >= self.maximum:
            raise ValueError(
                "Grid maximum must be greater than its minimum."
            )

        if self.size < 2:
            raise ValueError("Grid size must be at least two.")

        if self.quantity in {
            GridQuantity.VELOCITY,
            GridQuantity.SLOWNESS,
        }:
            if self.minimum <= 0.0:
                raise ValueError(
                    "Velocity and slowness values must be positive."
                )

        if self.quantity is GridQuantity.WAVENUMBER:
            if self.minimum < 0.0:
                raise ValueError(
                    "Wavenumber values cannot be negative."
                )

            if self.scale is GridScale.LOG:
                raise ValueError(
                    "Logarithmic wavenumber sampling is not supported."
                )

        if self.scale is GridScale.LOG and self.minimum <= 0.0:
            raise ValueError(
                "Logarithmic grid values must be positive."
            )

        if self.azimuths is None:
            if self.azimuth_size < 1:
                raise ValueError(
                    "azimuth_size must be at least one."
                )
        else:
            azimuths = np.asarray(self.azimuths, dtype=float)

            if azimuths.ndim != 1 or azimuths.size == 0:
                raise ValueError(
                    "azimuths must be a non-empty one-dimensional "
                    "array."
                )

            if not np.all(np.isfinite(azimuths)):
                raise ValueError(
                    "azimuths must contain only finite values."
                )

            azimuths = np.mod(azimuths, 360.0)

            if np.unique(azimuths).size != azimuths.size:
                raise ValueError(
                    "azimuths must contain unique directions."
                )

            object.__setattr__(self, "azimuths", azimuths.copy())

    @property
    def values(self) -> FloatArray:
        """Return the radial search values."""
        if self.scale is GridScale.LINEAR:
            return np.linspace(
                self.minimum,
                self.maximum,
                self.size,
                dtype=float,
            )

        return np.geomspace(
            self.minimum,
            self.maximum,
            self.size,
            dtype=float,
        )

    @property
    def azimuth_values(self) -> FloatArray:
        """Return propagation azimuths in degrees."""
        if self.azimuths is not None:
            return np.asarray(self.azimuths, dtype=float).copy()

        return np.linspace(
            0.0,
            360.0,
            self.azimuth_size,
            endpoint=False,
            dtype=float,
        )

    def wavenumbers(self, frequency: float) -> FloatArray:
        """
        Convert radial search values to angular wavenumbers.

        Parameters
        ----------
        frequency
            Analysis frequency in Hz.
        """
        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError("Frequency must be positive.")

        values = self.values

        if self.quantity is GridQuantity.VELOCITY:
            return 2.0 * np.pi * frequency / values

        if self.quantity is GridQuantity.SLOWNESS:
            return 2.0 * np.pi * frequency * values

        return values.copy()

    def velocity(self, frequency: float) -> FloatArray:
        """Return phase velocities corresponding to the grid."""
        if self.quantity is GridQuantity.VELOCITY:
            return self.values

        wavenumbers = self.wavenumbers(frequency)

        velocity = np.full(wavenumbers.shape, np.inf, dtype=float)
        nonzero = wavenumbers > 0.0

        velocity[nonzero] = (
            2.0 * np.pi * frequency / wavenumbers[nonzero]
        )

        return velocity

    def slowness(self, frequency: float) -> FloatArray:
        """Return phase slownesses corresponding to the grid."""
        if self.quantity is GridQuantity.SLOWNESS:
            return self.values

        velocity = self.velocity(frequency)

        slowness = np.zeros(velocity.shape, dtype=float)
        finite = np.isfinite(velocity) & (velocity > 0.0)

        slowness[finite] = 1.0 / velocity[finite]

        return slowness

    def wavevectors(
        self,
        frequency: float,
    ) -> tuple[FloatArray, FloatArray]:
        """
        Return kx and ky arrays over the polar grid.

        Returned arrays have shape ``(n_azimuths, n_values)``.
        """
        wavenumbers = self.wavenumbers(frequency)
        azimuths = np.deg2rad(self.azimuth_values)

        kx = np.sin(azimuths)[:, None] * wavenumbers[None, :]
        ky = np.cos(azimuths)[:, None] * wavenumbers[None, :]

        return kx, ky


@dataclass(frozen=True)
class EllipticityGrid:
    """
    Rayleigh ellipticity search grid sampled uniformly in angle.

    The signed ellipticity ratio ``e`` is parameterized through the
    ellipticity angle

        chi = arctan(e)

    so that

        e = tan(chi).

    Sampling uniformly in ``chi`` avoids imposing an arbitrary linear
    bound directly on ``e`` and gives a naturally bounded search
    variable. Angles are expressed in degrees.

    The interval must remain strictly inside (-90, 90) degrees because
    the tangent diverges at the endpoints.
    """

    minimum_angle: float = -85.0
    maximum_angle: float = 85.0
    size: int = 61

    def __post_init__(self) -> None:
        if not np.isfinite(self.minimum_angle):
            raise ValueError(
                "Ellipticity minimum angle must be finite."
            )

        if not np.isfinite(self.maximum_angle):
            raise ValueError(
                "Ellipticity maximum angle must be finite."
            )

        if not -90.0 < self.minimum_angle < 90.0:
            raise ValueError(
                "Ellipticity minimum angle must lie inside "
                "(-90, 90) degrees."
            )

        if not -90.0 < self.maximum_angle < 90.0:
            raise ValueError(
                "Ellipticity maximum angle must lie inside "
                "(-90, 90) degrees."
            )

        if self.minimum_angle >= self.maximum_angle:
            raise ValueError(
                "Ellipticity maximum angle must be greater than "
                "minimum angle."
            )

        if self.size < 3:
            raise ValueError(
                "Ellipticity grid size must be at least three."
            )

    @property
    def angles(self) -> FloatArray:
        """Return uniformly sampled signed ellipticity angles."""
        return np.linspace(
            self.minimum_angle,
            self.maximum_angle,
            self.size,
            dtype=float,
        )

    @property
    def values(self) -> FloatArray:
        """Return signed ellipticity ratios corresponding to angles."""
        return np.tan(
            np.deg2rad(
                self.angles
            )
        )


@dataclass(frozen=True)
class BeamformingConfig:
    """
    Complete array-processing configuration.

    Parameters
    ----------
    frequencies
        Analysis-frequency definition. A one-dimensional array keeps
        the historical direct-transform behavior. A
        :class:`spectral.FrequencyConfig` enables range-based
        definitions and native-bin RFFT sampling.
    grid
        Polar velocity, slowness, or wavenumber grid.
    window
        Elementary-window and covariance-block configuration defined
        by :class:`spectral.WindowConfig`.
    method
        Beam-power estimator.
    components
        Component combination to process.
    diagonal_loading
        Relative diagonal-loading coefficient. Loading is scaled by
        the mean diagonal power of the cross-spectral matrix.
    music_sources
        Number of signal eigenvectors excluded from the MUSIC noise
        subspace.
    output
        Whether to retain maxima, spectra, or both.
    maxima
        Maximum number of local maxima retained for each component,
        frequency, and covariance block.
    peak_neighborhood
        Neighborhood size used for local-maximum detection, expressed
        as ``(azimuth_samples, radial_samples)``.
    minimum_relative_power
        Minimum peak power relative to the global spectrum maximum.
    strict
        If True, require more Fourier realizations than matrix
        dimensions for high-resolution methods.
    transform
        Direct arbitrary-frequency transform or fixed-window RFFT.
    fft_selection
        Frequency selection used with ``RFFT``. ``NEAREST`` uses real
        FFT bins, while ``INTERP`` uses local complex cubic
        interpolation at the requested target frequencies.
    fft_padding
        Optional RFFT zero-padding factor. ``None`` uses the spectral
        layer default: one for nearest-bin selection and eight for
        complex interpolation.

    Notes
    -----
    ``RFFT`` requires fixed-duration windows because all processed
    frequencies must share the same elementary-window layout. Direct
    processing supports both fixed-duration and cycle-based windows.
    """

    frequencies: ArrayLike | FrequencyConfig
    grid: SearchGrid
    window: WindowConfig = WindowConfig()
    method: BeamformingMethod = BeamformingMethod.CAPON
    components: ComponentMode = ComponentMode.VERTICAL
    diagonal_loading: float = 1e-3
    music_sources: int = 1
    ellipticity_grid: EllipticityGrid = EllipticityGrid()
    output: OutputMode = OutputMode.MAXIMA
    maxima: int = 1
    peak_neighborhood: tuple[int, int] = (3, 3)
    minimum_relative_power: float = 0.0
    statistical_uncertainty: StatisticalUncertainty = (
        StatisticalUncertainty.NONE
    )
    geometric_uncertainty: GeometricUncertainty = (
        GeometricUncertainty.NONE
    )
    strict: bool = True
    transform: TransformMethod = TransformMethod.DIRECT
    fft_selection: FFTSelection = FFTSelection.NEAREST
    fft_padding: int | None = None

    def __post_init__(self) -> None:
        if isinstance(
            self.frequencies,
            FrequencyConfig,
        ):
            frequency_config = self.frequencies

        else:
            frequencies = np.asarray(
                self.frequencies,
                dtype=float,
            )

            if frequencies.ndim != 1 or frequencies.size == 0:
                raise ValueError(
                    "frequencies must be a non-empty one-dimensional "
                    "array or a FrequencyConfig instance."
                )

            if not np.all(
                np.isfinite(frequencies)
            ):
                raise ValueError(
                    "frequencies must contain only finite values."
                )

            if np.any(
                frequencies <= 0.0
            ):
                raise ValueError(
                    "frequencies must be strictly positive."
                )

            if np.any(
                np.diff(frequencies) <= 0.0
            ):
                raise ValueError(
                    "frequencies must be strictly increasing."
                )

            frequencies = frequencies.copy()

            object.__setattr__(
                self,
                "frequencies",
                frequencies,
            )

            frequency_config = FrequencyConfig.from_values(
                frequencies
            )

        if not isinstance(
            self.grid,
            SearchGrid,
        ):
            raise TypeError(
                "grid must be a SearchGrid instance."
            )

        if not isinstance(
            self.window,
            WindowConfig,
        ):
            raise TypeError(
                "window must be a WindowConfig instance."
            )

        if not isinstance(
            self.transform,
            TransformMethod,
        ):
            raise TypeError(
                "transform must be a TransformMethod value."
            )

        if not isinstance(
            self.fft_selection,
            FFTSelection,
        ):
            raise TypeError(
                "fft_selection must be an FFTSelection value."
            )

        if not isinstance(
            self.method,
            BeamformingMethod,
        ):
            raise TypeError(
                "method must be a BeamformingMethod value."
            )

        if not isinstance(
            self.components,
            ComponentMode,
        ):
            raise TypeError(
                "components must be a ComponentMode value."
            )

        if not isinstance(
            self.ellipticity_grid,
            EllipticityGrid,
        ):
            raise TypeError(
                "ellipticity_grid must be an EllipticityGrid instance."
            )

        if not isinstance(
            self.output,
            OutputMode,
        ):
            raise TypeError(
                "output must be an OutputMode value."
            )

        if not isinstance(
            self.statistical_uncertainty,
            StatisticalUncertainty,
        ):
            raise TypeError(
                "statistical_uncertainty must be a "
                "StatisticalUncertainty value."
            )

        if not isinstance(
            self.geometric_uncertainty,
            GeometricUncertainty,
        ):
            raise TypeError(
                "geometric_uncertainty must be a "
                "GeometricUncertainty value."
            )

        SpectralConfig(
            frequencies=frequency_config,
            window=self.window,
            transform=self.transform,
            fft_selection=self.fft_selection,
            fft_padding=self.fft_padding,
        )

        if (
            self.statistical_uncertainty
            is StatisticalUncertainty.EMPIRICAL
            and self.method
            not in {
                BeamformingMethod.BARTLETT,
                BeamformingMethod.CAPON,
            }
        ):
            raise ValueError(
                "Empirical beam-power uncertainty is currently "
                "implemented for Bartlett and Capon beamforming."
            )

        if (
            self.components is ComponentMode.RAYLEIGH_JOINT
            and self.method is not BeamformingMethod.CAPON
        ):
            raise ValueError(
                "Rayleigh joint processing currently implements "
                "high-resolution RTBF and requires Capon."
            )

        if (
            self.components is ComponentMode.RAYLEIGH_JOINT
            and self.statistical_uncertainty
            is not StatisticalUncertainty.NONE
        ):
            raise ValueError(
                "Statistical uncertainty is not yet implemented "
                "for Rayleigh joint processing."
            )

        if (
            self.components is ComponentMode.RAYLEIGH_JOINT
            and self.geometric_uncertainty
            is not GeometricUncertainty.NONE
        ):
            raise ValueError(
                "Geometric uncertainty is not yet implemented "
                "for Rayleigh joint processing."
            )

        if (
            self.geometric_uncertainty
            is GeometricUncertainty.EXPECTED_COHERENCE
            and self.method is not BeamformingMethod.BARTLETT
        ):
            raise ValueError(
                "Expected-coherence coordinate uncertainty is "
                "currently implemented only for Bartlett "
                "beamforming."
            )

        if not np.isfinite(
            self.diagonal_loading
        ):
            raise ValueError(
                "diagonal_loading must be finite."
            )

        if self.diagonal_loading < 0.0:
            raise ValueError(
                "diagonal_loading cannot be negative."
            )

        if self.music_sources < 1:
            raise ValueError(
                "music_sources must be at least one."
            )

        if self.maxima < 1:
            raise ValueError(
                "maxima must be at least one."
            )

        if len(
            self.peak_neighborhood
        ) != 2:
            raise ValueError(
                "peak_neighborhood must contain two integers."
            )

        if any(
            value < 1
            for value in self.peak_neighborhood
        ):
            raise ValueError(
                "Peak-neighborhood values must be positive."
            )

        if not (
            0.0
            <= self.minimum_relative_power
            <= 1.0
        ):
            raise ValueError(
                "minimum_relative_power must be in the range [0, 1]."
            )

    @property
    def spectral_config(self) -> SpectralConfig:
        """
        Return the spectral configuration used by the beamformer.

        Existing explicit frequency arrays are converted lazily to a
        :class:`FrequencyConfig`, preserving the historical public
        constructor while allowing the spectral layer to operate on a
        single coherent configuration object.
        """
        if isinstance(
            self.frequencies,
            FrequencyConfig,
        ):
            frequencies = self.frequencies
        else:
            frequencies = FrequencyConfig.from_values(
                self.frequencies
            )

        return SpectralConfig(
            frequencies=frequencies,
            window=self.window,
            transform=self.transform,
            fft_selection=self.fft_selection,
            fft_padding=self.fft_padding,
        )


@dataclass(frozen=True)
class SpectralMaximum:
    """
    Local maximum extracted from one beam-power spectrum.

    The canonical propagation coordinates stored in the result are the
    radial angular wavenumber and propagation azimuth. Velocity,
    slowness, and Cartesian wavevector components are derived
    properties and are therefore not duplicated in the stored state.

    ``rank`` is one-based and orders local maxima by decreasing power
    within the same block, frequency, and component.
    """

    block: int
    time: object
    frequency: float
    component: str
    rank: int
    wavenumber: float
    azimuth: float
    power: float
    power_standard_error: float | None = None
    ellipticity: float | None = None
    ellipticity_angle: float | None = None

    @property
    def velocity(self) -> float:
        """Return phase velocity."""
        if self.wavenumber == 0.0:
            return float("inf")

        return 2.0 * np.pi * self.frequency / self.wavenumber

    @property
    def slowness(self) -> float:
        """Return phase slowness."""
        return self.wavenumber / (2.0 * np.pi * self.frequency)

    @property
    def kx(self) -> float:
        """Return East Cartesian angular wavenumber."""
        return float(
            self.wavenumber
            * np.sin(np.deg2rad(self.azimuth))
        )

    @property
    def ky(self) -> float:
        """Return North Cartesian angular wavenumber."""
        return float(
            self.wavenumber
            * np.cos(np.deg2rad(self.azimuth))
        )


@dataclass(frozen=True)
class FrequencySpectrum:
    """
    Beam-power spectra for one frequency and component.

    Parameters
    ----------
    frequency
        Analysis frequency.
    component
        Name of the processed component.
    times
        Representative times of covariance blocks.
    azimuths
        Propagation azimuths in degrees.
    grid_values
        Radial search values in the configured physical quantity.
    quantity
        Meaning of ``grid_values``.
    power
        Beam-power array with shape
        ``(n_blocks, n_azimuths, n_values)``.
    """

    frequency: float
    component: str
    times: tuple[object, ...]
    azimuths: FloatArray
    grid_values: FloatArray
    quantity: GridQuantity
    n_realizations: int
    effective_realizations: float
    power: FloatArray
    power_standard_error: FloatArray | None = None
    ellipticity: FloatArray | None = None
    ellipticity_angle: FloatArray | None = None

    def wavenumbers(self) -> FloatArray:
        """Return angular wavenumbers for this spectrum."""
        if self.quantity is GridQuantity.WAVENUMBER:
            return self.grid_values.copy()

        if self.quantity is GridQuantity.VELOCITY:
            return (
                2.0
                * np.pi
                * self.frequency
                / self.grid_values
            )

        return (
            2.0
            * np.pi
            * self.frequency
            * self.grid_values
        )

    def wavevectors(
        self,
    ) -> tuple[FloatArray, FloatArray]:
        """Return kx and ky arrays for plotting or export."""
        wavenumbers = self.wavenumbers()
        azimuths = np.deg2rad(self.azimuths)

        kx = np.sin(azimuths)[:, None] * wavenumbers[None, :]
        ky = np.cos(azimuths)[:, None] * wavenumbers[None, :]

        return kx, ky


@dataclass(frozen=True)
class BeamformingResult:
    """Complete result of one array-processing operation."""

    maxima: tuple[SpectralMaximum, ...]
    spectra: tuple[FrequencySpectrum, ...] | None
    config: BeamformingConfig


def _processing_directions(
    mode: ComponentMode,
) -> tuple[tuple[FloatArray, ...], tuple[str, ...]]:
    """Return canonical directions and labels required by a mode."""
    if mode is ComponentMode.VERTICAL:
        return (VERTICAL,), ("vertical",)

    if mode in {
        ComponentMode.RADIAL,
        ComponentMode.TRANSVERSE,
    }:
        return (EAST, NORTH), ("east", "north")

    if mode in {
        ComponentMode.THREE_COMPONENT_SEPARATE,
        ComponentMode.RAYLEIGH_JOINT,
    }:
        return (EAST, NORTH, VERTICAL), (
            "east",
            "north",
            "vertical",
        )

    raise ValueError(
        f"Unsupported component mode: {mode!r}."
    )


def _prepare_waveforms(
    array: ArrayData,
    mode: ComponentMode,
) -> dict[str, FloatArray]:
    """
    Project the input snapshot once onto the required directions.

    Returns
    -------
    dict
        Canonical component matrices with shape
        ``(n_stations, n_samples)``.
    """
    directions, labels = _processing_directions(mode)
    matrix = array.matrix(directions)

    return {
        label: matrix[:, index, :]
        for index, label in enumerate(labels)
    }


def _geometric_attenuation(
    coordinate_covariances: FloatArray,
    frequency: float,
    grid: SearchGrid,
) -> FloatArray:
    """Return Gaussian phase-coherence attenuation factors.

    Station-coordinate errors are assumed independent between
    stations and Gaussian with the covariance matrices supplied by
    ``ArrayData``. Only the horizontal 2 x 2 covariance block enters
    the present 2-D frequency-wavenumber steering model.

    The returned array has shape
    ``(n_azimuths, n_values, n_stations)`` and contains

    ``exp(-0.5 * k.T @ Sigma_i @ k)``.
    """
    covariance = np.asarray(
        coordinate_covariances,
        dtype=float,
    )

    if covariance.ndim != 3 or covariance.shape[1:] != (3, 3):
        raise ValueError(
            "Coordinate covariances must have shape "
            "(n_stations, 3, 3)."
        )

    if not np.all(np.isfinite(covariance)):
        raise ValueError(
            "Coordinate covariances must contain only finite values."
        )

    kx, ky = grid.wavevectors(frequency)
    sigma_xx = covariance[:, 0, 0]
    sigma_xy = covariance[:, 0, 1]
    sigma_yy = covariance[:, 1, 1]

    phase_variance = (
        kx[:, :, None] ** 2 * sigma_xx[None, None, :]
        + 2.0
        * kx[:, :, None]
        * ky[:, :, None]
        * sigma_xy[None, None, :]
        + ky[:, :, None] ** 2 * sigma_yy[None, None, :]
    )

    # Roundoff in positive-semidefinite covariance matrices can
    # produce tiny negative values.
    phase_variance = np.maximum(phase_variance, 0.0)

    return np.exp(-0.5 * phase_variance)


def _search_steering_vectors(
    coordinates: FloatArray,
    frequency: float,
    grid: SearchGrid,
) -> ComplexArray:
    """
    Return shared estimator steering vectors on a YArray SearchGrid.

    This adapter keeps grid parametrization in the beamforming layer
    while the low-level steering kernel remains independent of
    frequency, velocity, slowness, and SearchGrid.
    """
    kx, ky = grid.wavevectors(
        frequency
    )

    return steering_vectors(
        coordinates,
        kx,
        ky,
    )


def _estimator_power(
    matrix: ComplexArray,
    steering: ComplexArray,
    n_realizations: int,
    config: BeamformingConfig,
    attenuation: FloatArray | None = None,
) -> FloatArray:
    """Evaluate the configured shared scalar estimator."""
    return beam_power(
        matrix,
        steering,
        n_realizations,
        config.method,
        strict=config.strict,
        diagonal_loading=config.diagonal_loading,
        music_sources=config.music_sources,
        attenuation=attenuation,
    )


def _horizontal_matrix(
    blocks: SpectralMatrixBlocks,
    azimuth: float,
    *,
    transverse: bool,
) -> ComplexArray:
    """
    Construct radial or transverse cross-spectral matrix.

    Parameters
    ----------
    blocks
        Horizontal cross-spectral matrix blocks.
    azimuth
        Propagation azimuth in degrees clockwise from North.
    transverse
        If True, construct the transverse matrix. Otherwise construct
        the radial matrix.
    """
    if blocks.ee is None:
        raise ValueError(
            "East cross-spectral matrix is not available."
        )

    if blocks.nn is None:
        raise ValueError(
            "North cross-spectral matrix is not available."
        )

    if blocks.en is None:
        raise ValueError(
            "East-North cross-spectral matrix is not available."
        )

    angle = np.deg2rad(azimuth)

    radial_east = np.sin(angle)
    radial_north = np.cos(angle)

    if transverse:
        weight_east = np.cos(angle)
        weight_north = -np.sin(angle)
    else:
        weight_east = radial_east
        weight_north = radial_north

    ne = blocks.en.conj().T

    matrix = (
        weight_east**2 * blocks.ee
        + weight_north**2 * blocks.nn
        + weight_east
        * weight_north
        * (blocks.en + ne)
    )

    return hermitian(matrix)


def _horizontal_power(
    blocks: SpectralMatrixBlocks,
    steering: ComplexArray,
    azimuths: FloatArray,
    n_realizations: int,
    config: BeamformingConfig,
    attenuation: FloatArray | None = None,
    *,
    transverse: bool,
) -> FloatArray:
    """Evaluate radial or transverse power for all grid points."""
    n_azimuths = steering.shape[0]
    n_values = steering.shape[1]

    power = np.empty(
        (n_azimuths, n_values),
        dtype=float,
    )

    for index, azimuth in enumerate(azimuths):
        matrix = _horizontal_matrix(
            blocks,
            float(azimuth),
            transverse=transverse,
        )

        direction_steering = steering[
            index:index + 1
        ]

        direction_attenuation = (
            None
            if attenuation is None
            else attenuation[index:index + 1]
        )

        direction_power = _estimator_power(
            matrix,
            direction_steering,
            n_realizations,
            config,
            attenuation=direction_attenuation,
        )

        power[index] = direction_power[0]

    return power


def _component_power_standard_error(
    component: str,
    spectral: SpectralComponents,
    block_start: int,
    windows_per_block: int,
    steering: ComplexArray,
    azimuths: FloatArray,
    effective_realizations: float,
    config: BeamformingConfig,
    attenuation: FloatArray | None = None,
) -> FloatArray:
    """Estimate empirical beam-power standard error.

    Bartlett uncertainty is computed directly from the dispersion of
    elementary beam powers. Capon uncertainty uses a first-order delta
    method applied to the stacked spectral-matrix estimator.
    """
    block = slice(
        block_start,
        block_start + windows_per_block,
    )

    def estimate(
        coefficients: ComplexArray,
        direction_steering: ComplexArray,
        direction_attenuation: FloatArray | None,
    ) -> FloatArray:
        if config.method is BeamformingMethod.BARTLETT:
            _, standard_error = bartlett_power_statistics(
                coefficients,
                direction_steering,
                effective_realizations,
                attenuation=direction_attenuation,
            )
            return standard_error

        if config.method is BeamformingMethod.CAPON:
            if direction_attenuation is not None:
                raise ValueError(
                    "Geometric uncertainty is not implemented for "
                    "Capon beamforming."
                )

            _, standard_error = capon_power_statistics(
                coefficients,
                direction_steering,
                effective_realizations,
                config.diagonal_loading,
            )
            return standard_error

        raise ValueError(
            "Empirical beam-power uncertainty is not implemented "
            f"for {config.method.value} beamforming."
        )

    if component == "vertical":
        if spectral.vertical is None:
            raise ValueError(
                "Vertical Fourier coefficients are not available."
            )

        return estimate(
            spectral.vertical[block],
            steering,
            attenuation,
        )

    if spectral.east is None or spectral.north is None:
        raise ValueError(
            "Horizontal Fourier coefficients are not available."
        )

    east = spectral.east[block]
    north = spectral.north[block]

    standard_error = np.empty(
        steering.shape[:-1],
        dtype=float,
    )

    for index, azimuth in enumerate(azimuths):
        angle = np.deg2rad(float(azimuth))

        if component == "radial":
            weight_east = np.sin(angle)
            weight_north = np.cos(angle)

        elif component == "transverse":
            weight_east = np.cos(angle)
            weight_north = -np.sin(angle)

        else:
            raise ValueError(
                f"Unknown output component: {component!r}."
            )

        coefficients = (
            weight_east * east
            + weight_north * north
        )

        direction_attenuation = (
            None
            if attenuation is None
            else attenuation[index:index + 1]
        )

        direction_error = estimate(
            coefficients,
            steering[index:index + 1],
            direction_attenuation,
        )

        standard_error[index] = direction_error[0]

    return standard_error

def _component_names(
    mode: ComponentMode,
) -> tuple[str, ...]:
    """Return output component names for a processing mode."""
    if mode is ComponentMode.VERTICAL:
        return ("vertical",)

    if mode is ComponentMode.RADIAL:
        return ("radial",)

    if mode is ComponentMode.TRANSVERSE:
        return ("transverse",)

    if mode is ComponentMode.THREE_COMPONENT_SEPARATE:
        return (
            "vertical",
            "radial",
            "transverse",
        )

    if mode is ComponentMode.RAYLEIGH_JOINT:
        return ("rayleigh_joint",)

    raise ValueError(
        f"Unsupported component mode: {mode!r}."
    )


def _rayleigh_joint_coefficients(
    spectral: SpectralComponents,
    block_start: int,
    windows_per_block: int,
    azimuth: float,
) -> ComplexArray:
    """Build radial-vertical Fourier realizations for one azimuth."""
    if (
        spectral.east is None
        or spectral.north is None
        or spectral.vertical is None
    ):
        raise ValueError(
            "Rayleigh joint processing requires East, North, and "
            "vertical Fourier coefficients."
        )

    block = slice(
        block_start,
        block_start + windows_per_block,
    )

    angle = np.deg2rad(float(azimuth))

    radial = (
        np.sin(angle) * spectral.east[block]
        + np.cos(angle) * spectral.north[block]
    )

    vertical = spectral.vertical[block]

    return np.concatenate(
        (
            radial,
            vertical,
        ),
        axis=1,
    )


def _rayleigh_joint_power(
    spectral: SpectralComponents,
    block_start: int,
    windows_per_block: int,
    steering: ComplexArray,
    azimuths: FloatArray,
    config: BeamformingConfig,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Evaluate high-resolution Rayleigh three-component beam power.

    For each propagation azimuth, the horizontal Fourier coefficients
    are projected onto the radial direction and stacked with the
    vertical coefficients, forming a 2N-channel realization matrix.

    The Capon steering vector follows the high-resolution RTBF
    radial-shift convention of Wathelet et al. (2018):

        E_Rh = [-j e q, q]^T

    where ``e`` is signed Rayleigh ellipticity and ``q`` is the
    N-station plane-wave steering vector.

    The robust signed-ellipticity objective is

        P_Rs = P_Rh * P_Rz = e**2 * P_Rh**2.

    The method returns the maximum objective over the configured
    ellipticity grid and the signed ellipticity at that maximum.
    """
    if config.method is not BeamformingMethod.CAPON:
        raise ValueError(
            "Rayleigh joint processing requires Capon."
        )

    ellipticity_angles = config.ellipticity_grid.angles
    ellipticity = config.ellipticity_grid.values

    # e = 0 makes the alternative vertical-shift power singular and
    # carries no signed-ellipticity information.
    nonzero = ~np.isclose(
        ellipticity,
        0.0,
        rtol=0.0,
        atol=np.finfo(float).eps,
    )

    ellipticity_angles = ellipticity_angles[
        nonzero
    ]
    ellipticity = ellipticity[
        nonzero
    ]

    n_azimuths, n_values, n_stations = steering.shape

    if windows_per_block <= 2 * n_stations and config.strict:
        raise ValueError(
            "Strict Rayleigh joint Capon processing requires more "
            "Fourier realizations than the 2N matrix dimension: "
            f"got {windows_per_block} for dimension {2 * n_stations}."
        )

    power = np.empty(
        (n_azimuths, n_values),
        dtype=float,
    )
    best_ellipticity = np.empty(
        (n_azimuths, n_values),
        dtype=float,
    )
    best_ellipticity_angle = np.empty(
        (n_azimuths, n_values),
        dtype=float,
    )

    for azimuth_index, azimuth in enumerate(azimuths):
        coefficients = _rayleigh_joint_coefficients(
            spectral,
            block_start,
            windows_per_block,
            float(azimuth),
        )

        matrix = cross_spectrum(
            coefficients
        )
        regularized = regularized_matrix(
            matrix,
            config.diagonal_loading,
        )

        q = steering[azimuth_index]

        radial = (
            -1j
            * ellipticity[:, None, None]
            * q[None, :, :]
        )
        vertical = np.broadcast_to(
            q[None, :, :],
            radial.shape,
        )

        joint = np.concatenate(
            (
                radial,
                vertical,
            ),
            axis=2,
        )

        flat = joint.reshape(
            -1,
            2 * n_stations,
        )

        solved = np.linalg.solve(
            regularized,
            flat.T,
        )

        denominator = np.real(
            np.einsum(
                "gi,ig->g",
                flat.conj(),
                solved,
                optimize=True,
            )
        )

        p_rh = np.full(
            denominator.shape,
            np.nan,
            dtype=float,
        )

        valid = denominator > 0.0
        p_rh[valid] = 1.0 / denominator[valid]

        p_rh = p_rh.reshape(
            ellipticity.size,
            n_values,
        )

        objective = (
            ellipticity[:, None] ** 2
            * p_rh**2
        )

        safe = np.where(
            np.isfinite(objective),
            objective,
            -np.inf,
        )

        winner = np.argmax(
            safe,
            axis=0,
        )

        columns = np.arange(
            n_values,
        )

        power[azimuth_index] = objective[
            winner,
            columns,
        ]
        best_ellipticity[azimuth_index] = ellipticity[
            winner
        ]
        best_ellipticity_angle[azimuth_index] = (
            ellipticity_angles[
                winner
            ]
        )

    return (
        power,
        best_ellipticity,
        best_ellipticity_angle,
    )


def _component_power(
    component: str,
    blocks: SpectralMatrixBlocks,
    steering: ComplexArray,
    azimuths: FloatArray,
    n_realizations: int,
    config: BeamformingConfig,
    attenuation: FloatArray | None = None,
    *,
    spectral: SpectralComponents | None = None,
    block_start: int | None = None,
) -> FloatArray | tuple[
    FloatArray,
    FloatArray,
    FloatArray,
]:
    """Compute one component spectrum."""
    if component == "rayleigh_joint":
        if spectral is None or block_start is None:
            raise ValueError(
                "Rayleigh joint processing requires spectral "
                "coefficients and block_start."
            )

        return _rayleigh_joint_power(
            spectral,
            block_start,
            n_realizations,
            steering,
            azimuths,
            config,
        )

    if component == "vertical":
        if blocks.zz is None:
            raise ValueError(
                "Vertical cross-spectral matrix is not available."
            )

        return _estimator_power(
            blocks.zz,
            steering,
            n_realizations,
            config,
            attenuation=attenuation,
        )

    if component == "radial":
        return _horizontal_power(
            blocks,
            steering,
            azimuths,
            n_realizations,
            config,
            attenuation=attenuation,
            transverse=False,
        )

    if component == "transverse":
        return _horizontal_power(
            blocks,
            steering,
            azimuths,
            n_realizations,
            config,
            attenuation=attenuation,
            transverse=True,
        )

    raise ValueError(
        f"Unknown output component: {component!r}."
    )


def _block_time(
    starttime: object,
    delta: float,
    windows: FrequencyWindows,
    block_start: int,
    windows_per_block: int,
) -> object:
    """Return the temporal midpoint of a covariance block."""
    first_window = block_start
    last_window = (
        block_start
        + windows_per_block
        - 1
    )

    first_sample = int(
        windows.starts[first_window]
    )

    last_sample = int(
        windows.starts[last_window]
        + windows.size
        - 1
    )

    midpoint = 0.5 * (
        first_sample + last_sample
    )

    return starttime + midpoint * delta


def _local_maxima(
    spectrum: FloatArray,
    count: int,
    neighborhood: tuple[int, int],
    minimum_relative_power: float,
) -> tuple[tuple[int, int], ...]:
    """Return indices of strongest local maxima."""
    spectrum = np.asarray(spectrum, dtype=float)

    if spectrum.ndim != 2:
        raise ValueError(
            "Peak picking requires a two-dimensional spectrum."
        )

    finite = np.isfinite(spectrum)

    if not np.any(finite):
        return ()

    work = np.where(
        finite,
        spectrum,
        -np.inf,
    )

    filtered = maximum_filter(
        work,
        size=neighborhood,
        mode=("wrap", "nearest"),
    )

    candidates = finite & (work == filtered)

    maximum = np.max(work[finite])

    if minimum_relative_power > 0.0:
        candidates &= (
            work
            >= minimum_relative_power * maximum
        )

    indices = np.argwhere(candidates)

    if indices.size == 0:
        return ()

    values = work[
        indices[:, 0],
        indices[:, 1],
    ]

    order = np.argsort(values)[::-1]
    indices = indices[order[:count]]

    return tuple(
        (int(row), int(column))
        for row, column in indices
    )


def _make_maximum(
    block: int,
    time: object,
    frequency: float,
    component: str,
    rank: int,
    row: int,
    column: int,
    power: FloatArray,
    grid: SearchGrid,
    power_standard_error: FloatArray | None = None,
    ellipticity: FloatArray | None = None,
    ellipticity_angle: FloatArray | None = None,
) -> SpectralMaximum:
    """Construct one spectral maximum from grid indices."""
    azimuth = float(grid.azimuth_values[row])
    wavenumber = float(grid.wavenumbers(frequency)[column])

    return SpectralMaximum(
        block=block,
        time=time,
        frequency=float(frequency),
        component=component,
        rank=rank,
        wavenumber=wavenumber,
        azimuth=azimuth,
        power=float(power[row, column]),
        power_standard_error=(
            None
            if power_standard_error is None
            else float(power_standard_error[row, column])
        ),
        ellipticity=(
            None
            if ellipticity is None
            else float(ellipticity[row, column])
        ),
        ellipticity_angle=(
            None
            if ellipticity_angle is None
            else float(ellipticity_angle[row, column])
        ),
    )


def beamform(
    array: ArrayData,
    config: BeamformingConfig,
) -> BeamformingResult:
    """
    Perform frequency-wavenumber beamforming.

    Parameters
    ----------
    array
        Immutable and synchronized array-data snapshot.
    config
        Beamforming configuration.

    Returns
    -------
    BeamformingResult
        Extracted local maxima and, optionally, complete spectra.

    Notes
    -----
    A separate elementary-window layout is generated for every
    frequency when frequency-dependent windowing is selected.
    Consequently, the number and temporal position of covariance
    blocks may vary with frequency.

    This implementation supports vertical, radial, transverse,
    separate three-component processing, and joint radial-vertical
    high-resolution RTBF processing.
    """
    component_names = _component_names(
        config.components
    )

    waveforms = _prepare_waveforms(
        array,
        config.components,
    )
    coordinates = array.coordinates

    geometric_enabled = (
        config.geometric_uncertainty
        is GeometricUncertainty.EXPECTED_COHERENCE
    )
    coordinate_covariances = (
        array.coordinate_covariances
        if geometric_enabled
        else None
    )

    if geometric_enabled and coordinate_covariances is None:
        raise ValueError(
            "Expected-coherence coordinate uncertainty requires "
            "coordinate covariance matrices in ArrayData."
        )

    statistical_enabled = (
        config.statistical_uncertainty
        is StatisticalUncertainty.EMPIRICAL
    )

    maxima_output: list[SpectralMaximum] = []
    spectra: list[FrequencySpectrum] = []

    slices = spectral_slices(
        waveforms,
        array.npts,
        array.delta,
        config.spectral_config,
    )

    for spectral_slice in slices:
        frequency = spectral_slice.frequency
        windows = spectral_slice.windows
        spectral = spectral_slice.components

        steering = _search_steering_vectors(
            coordinates,
            frequency,
            config.grid,
        )

        attenuation = (
            _geometric_attenuation(
                coordinate_covariances,
                frequency,
                config.grid,
            )
            if coordinate_covariances is not None
            else None
        )

        azimuths = config.grid.azimuth_values
        grid_values = config.grid.values

        component_buffers = {
            component: np.empty(
                (
                    windows.block_starts.size,
                    azimuths.size,
                    grid_values.size,
                ),
                dtype=float,
            )
            for component in component_names
        }

        ellipticity_buffers = (
            {
                component: np.empty(
                    (
                        windows.block_starts.size,
                        azimuths.size,
                        grid_values.size,
                    ),
                    dtype=float,
                )
                for component in component_names
            }
            if config.components is ComponentMode.RAYLEIGH_JOINT
            else None
        )

        ellipticity_angle_buffers = (
            {
                component: np.empty(
                    (
                        windows.block_starts.size,
                        azimuths.size,
                        grid_values.size,
                    ),
                    dtype=float,
                )
                for component in component_names
            }
            if config.components is ComponentMode.RAYLEIGH_JOINT
            else None
        )

        standard_error_buffers = (
            {
                component: np.empty(
                    (
                        windows.block_starts.size,
                        azimuths.size,
                        grid_values.size,
                    ),
                    dtype=float,
                )
                for component in component_names
            }
            if statistical_enabled
            else None
        )

        effective_realizations = spectral_effective_realizations(
            windows,
            config.window.windows_per_block,
        )

        block_times: list[object] = []

        for block_index, block_start in enumerate(
            windows.block_starts
        ):
            block_start = int(block_start)

            time = _block_time(
                array.starttime,
                array.delta,
                windows,
                block_start,
                config.window.windows_per_block,
            )

            block_times.append(time)

            matrix_blocks = build_matrix_blocks(
                spectral,
                block_start,
                config.window.windows_per_block,
            )

            for component in component_names:
                power_result = _component_power(
                    component,
                    matrix_blocks,
                    steering,
                    azimuths,
                    config.window.windows_per_block,
                    config,
                    attenuation=attenuation,
                    spectral=spectral,
                    block_start=block_start,
                )

                if isinstance(power_result, tuple):
                    (
                        power,
                        ellipticity,
                        ellipticity_angle,
                    ) = power_result
                else:
                    power = power_result
                    ellipticity = None
                    ellipticity_angle = None

                component_buffers[component][
                    block_index
                ] = power

                if (
                    ellipticity_buffers is not None
                    and ellipticity is not None
                ):
                    ellipticity_buffers[
                        component
                    ][block_index] = ellipticity

                if (
                    ellipticity_angle_buffers is not None
                    and ellipticity_angle is not None
                ):
                    ellipticity_angle_buffers[
                        component
                    ][block_index] = ellipticity_angle

                power_standard_error = None

                if statistical_enabled:
                    power_standard_error = (
                        _component_power_standard_error(
                            component,
                            spectral,
                            block_start,
                            config.window.windows_per_block,
                            steering,
                            azimuths,
                            effective_realizations,
                            config,
                            attenuation=attenuation,
                        )
                    )
                    standard_error_buffers[component][
                        block_index
                    ] = power_standard_error

                if config.output in {
                    OutputMode.MAXIMA,
                    OutputMode.BOTH,
                }:
                    maxima = _local_maxima(
                        power,
                        config.maxima,
                        config.peak_neighborhood,
                        config.minimum_relative_power,
                    )

                    for rank, (row, column) in enumerate(
                        maxima,
                        start=1,
                    ):
                        maxima_output.append(
                            _make_maximum(
                                block=block_index,
                                time=time,
                                frequency=frequency,
                                component=component,
                                rank=rank,
                                row=row,
                                column=column,
                                power=power,
                                grid=config.grid,
                                power_standard_error=(
                                    power_standard_error
                                ),
                                ellipticity=ellipticity,
                                ellipticity_angle=(
                                    ellipticity_angle
                                ),
                            )
                        )

        if config.output in {
            OutputMode.SPECTRA,
            OutputMode.BOTH,
        }:
            for component in component_names:
                spectra.append(
                    FrequencySpectrum(
                        frequency=frequency,
                        component=component,
                        times=tuple(block_times),
                        azimuths=azimuths.copy(),
                        grid_values=grid_values.copy(),
                        quantity=config.grid.quantity,
                        n_realizations=(
                            config.window.windows_per_block
                        ),
                        effective_realizations=(
                            effective_realizations
                        ),
                        power=component_buffers[
                            component
                        ],
                        power_standard_error=(
                            None
                            if standard_error_buffers is None
                            else standard_error_buffers[component]
                        ),
                        ellipticity=(
                            None
                            if ellipticity_buffers is None
                            else ellipticity_buffers[component]
                        ),
                        ellipticity_angle=(
                            None
                            if ellipticity_angle_buffers is None
                            else ellipticity_angle_buffers[
                                component
                            ]
                        ),
                    )
                )

    spectra_output: tuple[FrequencySpectrum, ...] | None

    if config.output in {
        OutputMode.SPECTRA,
        OutputMode.BOTH,
    }:
        spectra_output = tuple(spectra)
    else:
        spectra_output = None

    return BeamformingResult(
        maxima=tuple(maxima_output),
        spectra=spectra_output,
        config=config,
    )
