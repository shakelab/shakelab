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
Reusable spectral primitives for YArray processing.

This module contains the time-frequency operations shared by
frequency-wavenumber beamforming, future FDD processing, and other
array methods requiring Fourier realizations or cross-spectral
matrices.

Two Fourier-transform strategies are supported:

``DIRECT``
    Evaluate the discrete-time Fourier transform directly at arbitrary
    requested frequencies. This is the natural strategy for
    frequency-dependent windows defined by a number of signal cycles.

``RFFT``
    Evaluate one real FFT for each fixed-duration elementary window
    and reuse it for all requested frequencies. Requested frequencies
    may be mapped to the nearest native FFT bins or obtained by local
    complex cubic interpolation between native bins.

The RFFT implementation never requires all native FFT bins to be
retained for all windows. Each elementary-window transform is reduced
immediately to the selected analysis frequencies.

Multitapering, frequency smoothing, and FDD-specific logic are not
implemented yet. The public abstractions are designed so these
features can be added without changing the beamforming layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal.windows import tukey


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]
IntArray = NDArray[np.int64]


__all__ = [
    "FFTSelection",
    "FrequencyConfig",
    "FrequencyScale",
    "FrequencyWindows",
    "SpectralComponents",
    "SpectralConfig",
    "SpectralMatrixBlocks",
    "SpectralMatrixEstimate",
    "SpectralSlice",
    "TransformMethod",
    "WindowConfig",
    "WindowType",
    "block_start_indices",
    "cross_spectrum",
    "effective_realizations",
    "estimate_spectral_matrix",
    "fourier_coefficients",
    "frequency_windows",
    "matrix_blocks",
    "normalized_tukey",
    "rfft_coefficients",
    "spectral_components",
    "spectral_slices",
    "window_sample_count",
    "window_start_indices",
]


class WindowType(Enum):
    """Window-length parametrization."""

    FIXED = "fixed"
    CYCLES = "cycles"


class TransformMethod(Enum):
    """Fourier-transform strategy."""

    DIRECT = "direct"
    RFFT = "rfft"


class FFTSelection(Enum):
    """Selection of requested frequencies from a native RFFT grid."""

    NEAREST = "nearest"
    INTERP = "interp"


class FrequencyScale(Enum):
    """Sampling scale used to generate target frequencies."""

    LINEAR = "linear"
    LOG = "log"


@dataclass(frozen=True)
class WindowConfig:
    """
    Configuration of elementary windows and covariance blocks.

    Parameters
    ----------
    type
        Fixed-duration or frequency-dependent windowing.
    length
        Window duration in seconds when ``type`` is ``FIXED``.
        Number of signal cycles when ``type`` is ``CYCLES``.
    taper
        Tukey-window shape parameter between zero and one.
    overlap
        Fractional overlap between consecutive elementary windows.
    windows_per_block
        Number of elementary windows stacked to estimate one
        cross-spectral matrix.
    block_overlap
        Fractional overlap between consecutive covariance blocks,
        expressed in terms of elementary windows.
    """

    type: WindowType = WindowType.CYCLES
    length: float = 20.0
    taper: float = 0.1
    overlap: float = 0.0
    windows_per_block: int = 20
    block_overlap: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.type, WindowType):
            raise TypeError(
                "type must be a WindowType value."
            )

        if not np.isfinite(self.length) or self.length <= 0.0:
            raise ValueError(
                "Window length must be positive."
            )

        if not 0.0 <= self.taper <= 1.0:
            raise ValueError(
                "Window taper must be in the range [0, 1]."
            )

        if not 0.0 <= self.overlap < 1.0:
            raise ValueError(
                "Window overlap must be in the range [0, 1)."
            )

        if self.windows_per_block < 1:
            raise ValueError(
                "windows_per_block must be at least one."
            )

        if not 0.0 <= self.block_overlap < 1.0:
            raise ValueError(
                "Block overlap must be in the range [0, 1)."
            )


@dataclass(frozen=True)
class FrequencyConfig:
    """
    Analysis-frequency definition.

    Frequencies may be supplied explicitly through ``values`` or
    generated from ``minimum``, ``maximum``, ``size``, and ``scale``.

    For ``RFFT`` processing, a range with ``size=None`` selects every
    native FFT bin between ``minimum`` and ``maximum``. When ``size``
    is specified, it defines target frequencies which are subsequently
    handled according to :class:`FFTSelection`.

    Parameters
    ----------
    values
        Explicit positive and strictly increasing target frequencies.
    minimum
        Minimum frequency for generated targets or native-bin ranges.
    maximum
        Maximum frequency for generated targets or native-bin ranges.
    size
        Number of generated target frequencies. ``None`` is permitted
        only for ``RFFT`` and means use every native bin in the range.
    scale
        Linear or logarithmic target-frequency spacing.
    """

    values: ArrayLike | None = None
    minimum: float | None = None
    maximum: float | None = None
    size: int | None = None
    scale: FrequencyScale = FrequencyScale.LOG

    def __post_init__(self) -> None:
        if not isinstance(self.scale, FrequencyScale):
            raise TypeError(
                "scale must be a FrequencyScale value."
            )

        if self.values is not None:
            if (
                self.minimum is not None
                or self.maximum is not None
                or self.size is not None
            ):
                raise ValueError(
                    "Explicit frequency values cannot be combined "
                    "with minimum, maximum, or size."
                )

            values = np.asarray(
                self.values,
                dtype=float,
            )

            _validate_frequency_array(
                values
            )

            object.__setattr__(
                self,
                "values",
                values.copy(),
            )

            return

        if self.minimum is None or self.maximum is None:
            raise ValueError(
                "Frequency ranges require minimum and maximum."
            )

        if (
            not np.isfinite(self.minimum)
            or self.minimum <= 0.0
        ):
            raise ValueError(
                "Frequency minimum must be finite and positive."
            )

        if (
            not np.isfinite(self.maximum)
            or self.maximum <= self.minimum
        ):
            raise ValueError(
                "Frequency maximum must be finite and greater than "
                "the minimum."
            )

        if self.size is not None and self.size < 2:
            raise ValueError(
                "Frequency size must be at least two."
            )

        if (
            self.scale is FrequencyScale.LOG
            and self.minimum <= 0.0
        ):
            raise ValueError(
                "Logarithmic frequency sampling requires a positive "
                "minimum."
            )

    @classmethod
    def from_values(
        cls,
        values: ArrayLike,
    ) -> "FrequencyConfig":
        """Construct an explicit frequency configuration."""
        return cls(
            values=values
        )

    @property
    def is_explicit(self) -> bool:
        """Return whether explicit target frequencies are stored."""
        return self.values is not None

    def targets(self) -> FloatArray:
        """
        Return explicit or generated target frequencies.

        Raises
        ------
        ValueError
            If this is an unsampled range with ``size=None``.
        """
        if self.values is not None:
            return np.asarray(
                self.values,
                dtype=float,
            ).copy()

        if self.size is None:
            raise ValueError(
                "This frequency range has no target size. Native FFT "
                "bins must be resolved from the window length."
            )

        if self.scale is FrequencyScale.LINEAR:
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


@dataclass(frozen=True)
class SpectralConfig:
    """
    Spectral-estimation configuration.

    Parameters
    ----------
    frequencies
        Analysis-frequency definition.
    window
        Elementary-window and covariance-block configuration.
    transform
        Direct arbitrary-frequency transform or fixed-window RFFT.
    fft_selection
        Frequency selection used with ``RFFT``. ``NEAREST`` maps each
        target to a real FFT bin. ``INTERP`` returns the exact target
        frequencies using local complex cubic interpolation.
    fft_padding
        Optional zero-padding factor used by ``RFFT``. ``None`` selects
        an adaptive default: one for ``NEAREST`` and eight for
        ``INTERP``. Zero padding increases frequency-grid density but
        does not increase physical spectral resolution.
    """

    frequencies: FrequencyConfig
    window: WindowConfig = WindowConfig()
    transform: TransformMethod = TransformMethod.DIRECT
    fft_selection: FFTSelection = FFTSelection.NEAREST
    fft_padding: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(
            self.frequencies,
            FrequencyConfig,
        ):
            raise TypeError(
                "frequencies must be a FrequencyConfig instance."
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

        if self.fft_padding is not None:
            if (
                isinstance(
                    self.fft_padding,
                    bool,
                )
                or not isinstance(
                    self.fft_padding,
                    int,
                )
                or self.fft_padding < 1
            ):
                raise ValueError(
                    "fft_padding must be a positive integer or None."
                )

        if (
            self.transform is TransformMethod.RFFT
            and self.window.type is not WindowType.FIXED
        ):
            raise ValueError(
                "RFFT processing requires fixed-duration windows."
            )

        if (
            self.transform is TransformMethod.DIRECT
            and not self.frequencies.is_explicit
            and self.frequencies.size is None
        ):
            raise ValueError(
                "Direct processing requires explicit frequencies or "
                "a frequency range with a finite size."
            )


    @property
    def resolved_fft_padding(self) -> int:
        """
        Return the effective RFFT zero-padding factor.

        Nearest-bin processing defaults to the native FFT grid.
        Interpolated processing defaults to eight-fold zero padding,
        which strongly reduces interpolation error while preserving
        the same physical window bandwidth.
        """
        if self.fft_padding is not None:
            return self.fft_padding

        if (
            self.frequencies.values is None
            and self.frequencies.size is None
        ):
            return 1

        if self.fft_selection is FFTSelection.INTERP:
            return 8

        return 1


@dataclass(frozen=True)
class FrequencyWindows:
    """
    Elementary-window layout for one analysis frequency.

    Attributes
    ----------
    size
        Number of samples in each elementary window.
    starts
        Start sample of every complete elementary window.
    block_starts
        Start elementary-window index of every complete covariance
        block.
    taper
        Mean-square-normalized taper applied to each elementary
        window.
    """

    size: int
    starts: IntArray
    block_starts: IntArray
    taper: FloatArray


@dataclass(frozen=True)
class SpectralComponents:
    """
    Fourier coefficients of canonical components.

    Available matrices have shape ``(n_windows, n_stations)``.
    """

    vertical: ComplexArray | None = None
    east: ComplexArray | None = None
    north: ComplexArray | None = None

    @property
    def n_windows(self) -> int:
        """Return the number of Fourier realizations."""
        sizes = {
            component.shape[0]
            for component in (
                self.vertical,
                self.east,
                self.north,
            )
            if component is not None
        }

        if not sizes:
            return 0

        if len(sizes) != 1:
            raise ValueError(
                "Spectral components contain inconsistent numbers "
                "of Fourier realizations."
            )

        return sizes.pop()


@dataclass(frozen=True)
class SpectralSlice:
    """
    Fourier realizations and window layout for one analysis frequency.

    ``frequency`` is always the frequency actually represented by the
    coefficients. With nearest-bin RFFT selection it therefore stores
    the selected native FFT frequency rather than the ideal target.
    """

    frequency: float
    windows: FrequencyWindows
    components: SpectralComponents


@dataclass(frozen=True)
class SpectralMatrixEstimate:
    """
    Cross-spectral matrix estimate and sampling uncertainty.

    ``variance`` is the element-wise sample variance of the complex
    elementary cross products. It combines real and imaginary
    dispersion through the squared complex modulus. The uncertainty
    of the stacked mean is obtained by dividing this variance by the
    effective number of independent realizations.
    """

    mean: ComplexArray
    variance: FloatArray | None
    n_realizations: int
    effective_realizations: float

    @property
    def standard_error(self) -> FloatArray | None:
        """Return element-wise standard error of the matrix mean."""
        if self.variance is None:
            return None

        return np.sqrt(
            self.variance / self.effective_realizations
        )


@dataclass(frozen=True)
class SpectralMatrixBlocks:
    """
    Canonical cross-spectral matrices for one covariance block.

    The matrices use the fixed East-North-Vertical component naming
    adopted by YArray. Missing components are represented by ``None``.
    """

    zz: ComplexArray | None = None
    ee: ComplexArray | None = None
    nn: ComplexArray | None = None
    en: ComplexArray | None = None
    ez: ComplexArray | None = None
    nz: ComplexArray | None = None


@dataclass(frozen=True)
class _RFFTPlan:
    """Resolved RFFT frequency-selection plan."""

    frequencies: FloatArray
    indices: IntArray | None = None
    positions: FloatArray | None = None


def _validate_frequency_array(
    values: FloatArray,
) -> None:
    """Validate an explicit frequency array."""
    if values.ndim != 1 or values.size == 0:
        raise ValueError(
            "Frequency values must be a non-empty one-dimensional "
            "array."
        )

    if not np.all(
        np.isfinite(values)
    ):
        raise ValueError(
            "Frequency values must contain only finite numbers."
        )

    if np.any(
        values <= 0.0
    ):
        raise ValueError(
            "Frequency values must be strictly positive."
        )

    if np.any(
        np.diff(values) <= 0.0
    ):
        raise ValueError(
            "Frequency values must be strictly increasing."
        )


def window_sample_count(
    frequency: float,
    delta: float,
    config: WindowConfig,
) -> int:
    """Return the number of samples in one elementary window."""
    if not isinstance(
        config,
        WindowConfig,
    ):
        raise TypeError(
            "config must be a WindowConfig instance."
        )

    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(
            "frequency must be finite and positive."
        )

    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError(
            "delta must be finite and positive."
        )

    if config.type is WindowType.FIXED:
        duration = config.length
    else:
        duration = config.length / frequency

    sample_count = int(
        np.floor(
            duration / delta
        )
    )

    if sample_count < 2:
        raise ValueError(
            f"Window at {frequency:g} Hz contains fewer than "
            "two samples."
        )

    return sample_count


def window_start_indices(
    npts: int,
    window_size: int,
    overlap: float,
) -> IntArray:
    """Return start indices of complete elementary windows."""
    if npts < 1:
        raise ValueError(
            "npts must be at least one."
        )

    if not 0.0 <= overlap < 1.0:
        raise ValueError(
            "Window overlap must be in the range [0, 1)."
        )

    if window_size < 2:
        raise ValueError(
            "Window size must be at least two samples."
        )

    if window_size > npts:
        return np.empty(
            0,
            dtype=np.int64,
        )

    step = max(
        1,
        int(
            np.floor(
                window_size
                * (1.0 - overlap)
            )
        ),
    )

    return np.arange(
        0,
        npts - window_size + 1,
        step,
        dtype=np.int64,
    )


def block_start_indices(
    n_windows: int,
    windows_per_block: int,
    overlap: float,
) -> IntArray:
    """Return start indices of complete covariance blocks."""
    if n_windows < 0:
        raise ValueError(
            "n_windows cannot be negative."
        )

    if windows_per_block < 1:
        raise ValueError(
            "windows_per_block must be at least one."
        )

    if not 0.0 <= overlap < 1.0:
        raise ValueError(
            "Block overlap must be in the range [0, 1)."
        )

    if windows_per_block > n_windows:
        return np.empty(
            0,
            dtype=np.int64,
        )

    step = max(
        1,
        int(
            np.floor(
                windows_per_block
                * (1.0 - overlap)
            )
        ),
    )

    return np.arange(
        0,
        n_windows - windows_per_block + 1,
        step,
        dtype=np.int64,
    )


def normalized_tukey(
    size: int,
    alpha: float,
) -> FloatArray:
    """Return a Tukey taper normalized to unit mean-square energy."""
    if size < 2:
        raise ValueError(
            "Taper size must be at least two samples."
        )

    if (
        not np.isfinite(alpha)
        or not 0.0 <= alpha <= 1.0
    ):
        raise ValueError(
            "Taper alpha must be finite and in the range [0, 1]."
        )

    window = tukey(
        size,
        alpha=alpha,
        sym=False,
    ).astype(
        float
    )

    energy = np.mean(
        window**2
    )

    if not np.isfinite(energy) or energy <= 0.0:
        raise ValueError(
            "The selected taper has non-positive energy."
        )

    return window / np.sqrt(
        energy
    )


def frequency_windows(
    npts: int,
    delta: float,
    frequency: float,
    config: WindowConfig,
) -> FrequencyWindows:
    """Generate elementary-window and covariance-block indices."""
    window_size = window_sample_count(
        frequency,
        delta,
        config,
    )

    starts = window_start_indices(
        npts,
        window_size,
        config.overlap,
    )

    if starts.size == 0:
        raise ValueError(
            f"No complete elementary window is available at "
            f"{frequency:g} Hz."
        )

    block_starts = block_start_indices(
        starts.size,
        config.windows_per_block,
        config.block_overlap,
    )

    if block_starts.size == 0:
        raise ValueError(
            f"Only {starts.size} elementary windows are available "
            f"at {frequency:g} Hz, but "
            f"{config.windows_per_block} are required per block."
        )

    return FrequencyWindows(
        size=window_size,
        starts=starts,
        block_starts=block_starts,
        taper=normalized_tukey(
            window_size,
            config.taper,
        ),
    )


def effective_realizations(
    windows: FrequencyWindows,
    n_realizations: int,
) -> float:
    """
    Estimate the effective number of independent windows.

    The correction accounts for correlation introduced solely by
    overlap of the deterministic taper. For window lag ``q``, the
    normalized taper overlap is squared because beam-power estimates
    depend quadratically on the Fourier coefficients.
    """
    if not isinstance(
        windows,
        FrequencyWindows,
    ):
        raise TypeError(
            "windows must be a FrequencyWindows instance."
        )

    if n_realizations < 1:
        raise ValueError(
            "n_realizations must be at least one."
        )

    if n_realizations > windows.starts.size:
        raise ValueError(
            "n_realizations cannot exceed the available windows."
        )

    if n_realizations == 1:
        return 1.0

    if windows.starts.size > 1:
        step = int(
            windows.starts[1]
            - windows.starts[0]
        )
    else:
        step = windows.size

    energy = float(
        np.dot(
            windows.taper,
            windows.taper,
        )
    )

    denominator = 1.0

    for lag in range(
        1,
        n_realizations,
    ):
        shift = lag * step

        if shift >= windows.size:
            break

        correlation = float(
            np.dot(
                windows.taper[:-shift],
                windows.taper[shift:],
            )
            / energy
        )

        denominator += (
            2.0
            * (
                1.0
                - lag / n_realizations
            )
            * correlation**2
        )

    return n_realizations / denominator


def fourier_coefficients(
    data: FloatArray,
    frequency: float,
    delta: float,
    windows: FrequencyWindows,
) -> ComplexArray:
    """
    Evaluate one arbitrary Fourier frequency for all windows.

    Parameters
    ----------
    data
        Matrix with shape ``(n_stations, n_samples)``.
    frequency
        Analysis frequency in Hz.
    delta
        Sampling interval in seconds.
    windows
        Elementary-window layout.

    Returns
    -------
    numpy.ndarray
        Complex coefficients with shape
        ``(n_windows, n_stations)``.
    """
    data = _validated_waveform_matrix(
        data
    )

    if not isinstance(
        windows,
        FrequencyWindows,
    ):
        raise TypeError(
            "windows must be a FrequencyWindows instance."
        )

    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(
            "frequency must be finite and positive."
        )

    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError(
            "delta must be finite and positive."
        )

    sample_time = (
        np.arange(
            windows.size,
            dtype=float,
        )
        * delta
    )

    kernel = (
        windows.taper
        * np.exp(
            -2j
            * np.pi
            * frequency
            * sample_time
        )
    )

    coefficients = np.empty(
        (
            windows.starts.size,
            data.shape[0],
        ),
        dtype=np.complex128,
    )

    scale = 2.0 / windows.size

    for index, start in enumerate(
        windows.starts
    ):
        stop = (
            int(start)
            + windows.size
        )

        segment = data[
            :,
            int(start):stop,
        ]

        coefficients[index] = (
            scale
            * (
                segment
                @ kernel
            )
        )

    return coefficients


def rfft_coefficients(
    data: FloatArray,
    frequencies: FloatArray,
    delta: float,
    windows: FrequencyWindows,
    *,
    selection: FFTSelection,
    padding: int = 1,
) -> ComplexArray:
    """
    Evaluate selected frequencies from fixed-window real FFTs.

    Parameters
    ----------
    data
        Matrix with shape ``(n_stations, n_samples)``.
    frequencies
        Frequencies returned by the resolved RFFT plan.
    delta
        Sampling interval in seconds.
    windows
        Common fixed-duration elementary-window layout.
    selection
        Nearest-bin extraction or complex cubic interpolation.
    padding
        Positive integer zero-padding factor applied before the RFFT.

    Returns
    -------
    numpy.ndarray
        Complex coefficients with shape
        ``(n_windows, n_frequencies, n_stations)``.
    """
    data = _validated_waveform_matrix(
        data
    )

    frequencies = np.asarray(
        frequencies,
        dtype=float,
    )

    _validate_frequency_array(
        frequencies
    )

    if not isinstance(
        windows,
        FrequencyWindows,
    ):
        raise TypeError(
            "windows must be a FrequencyWindows instance."
        )

    if not isinstance(
        selection,
        FFTSelection,
    ):
        raise TypeError(
            "selection must be an FFTSelection value."
        )

    if (
        isinstance(
            padding,
            bool,
        )
        or not isinstance(
            padding,
            int,
        )
        or padding < 1
    ):
        raise ValueError(
            "padding must be a positive integer."
        )

    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError(
            "delta must be finite and positive."
        )

    nfft = (
        windows.size
        * padding
    )

    native = np.fft.rfftfreq(
        nfft,
        d=delta,
    )

    nyquist = native[-1]

    if frequencies[-1] > nyquist:
        raise ValueError(
            "Requested RFFT frequency exceeds Nyquist."
        )

    coefficients = np.empty(
        (
            windows.starts.size,
            frequencies.size,
            data.shape[0],
        ),
        dtype=np.complex128,
    )

    scale = 2.0 / windows.size

    if selection is FFTSelection.NEAREST:
        indices = _nearest_bin_indices(
            native,
            frequencies,
        )

        if not np.allclose(
            native[indices],
            frequencies,
            rtol=0.0,
            atol=8.0
            * np.finfo(float).eps
            * max(
                1.0,
                float(nyquist),
            ),
        ):
            raise ValueError(
                "Nearest-bin RFFT coefficients must be requested at "
                "resolved native FFT frequencies."
            )

    else:
        indices = None
        positions = frequencies / (
            native[1] - native[0]
        )

    for window_index, start in enumerate(
        windows.starts
    ):
        stop = (
            int(start)
            + windows.size
        )

        segment = data[
            :,
            int(start):stop,
        ]

        tapered = (
            segment
            * windows.taper[None, :]
        )

        spectrum = (
            scale
            * np.fft.rfft(
                tapered,
                n=nfft,
                axis=1,
            )
        )

        if selection is FFTSelection.NEAREST:
            selected = spectrum[
                :,
                indices,
            ]
        else:
            selected = _complex_cubic_uniform(
                spectrum,
                positions,
            )

        coefficients[window_index] = (
            selected.T
        )

    return coefficients


def spectral_components(
    waveforms: dict[str, FloatArray],
    frequency: float,
    delta: float,
    windows: FrequencyWindows,
) -> SpectralComponents:
    """Compute direct Fourier coefficients for one frequency."""
    waveforms = _validated_waveforms(
        waveforms
    )

    return SpectralComponents(
        vertical=_optional_direct_component(
            waveforms.get(
                "vertical"
            ),
            frequency,
            delta,
            windows,
        ),
        east=_optional_direct_component(
            waveforms.get(
                "east"
            ),
            frequency,
            delta,
            windows,
        ),
        north=_optional_direct_component(
            waveforms.get(
                "north"
            ),
            frequency,
            delta,
            windows,
        ),
    )


def spectral_slices(
    waveforms: dict[str, FloatArray],
    npts: int,
    delta: float,
    config: SpectralConfig,
) -> Iterator[SpectralSlice]:
    """
    Yield Fourier realizations one analysis frequency at a time.

    The beamforming layer uses this generator without needing to know
    whether coefficients were produced by direct transforms or by one
    shared fixed-window RFFT calculation.
    """
    waveforms = _validated_waveforms(
        waveforms
    )

    if not isinstance(
        config,
        SpectralConfig,
    ):
        raise TypeError(
            "config must be a SpectralConfig instance."
        )

    if config.transform is TransformMethod.DIRECT:
        yield from _direct_spectral_slices(
            waveforms,
            npts,
            delta,
            config,
        )
        return

    yield from _rfft_spectral_slices(
        waveforms,
        npts,
        delta,
        config,
    )


def estimate_spectral_matrix(
    coefficients_a: ComplexArray,
    coefficients_b: ComplexArray | None = None,
    *,
    effective_realizations: float | None = None,
    estimate_variance: bool = False,
) -> SpectralMatrixEstimate:
    """
    Estimate a cross-spectral matrix and optional uncertainty.

    The element-wise complex sample variance is evaluated without
    materializing the elementary outer-product matrices. For
    ``z_m = a_m b_m^H``,

    ``Var(z) = L / (L - 1) * (E[|z|^2] - |E[z]|^2)``.
    """
    coefficients_a = np.asarray(
        coefficients_a,
        dtype=np.complex128,
    )

    if coefficients_b is None:
        coefficients_b = coefficients_a
    else:
        coefficients_b = np.asarray(
            coefficients_b,
            dtype=np.complex128,
        )

    if coefficients_a.ndim != 2:
        raise ValueError(
            "Fourier coefficients must be two-dimensional."
        )

    if coefficients_a.shape != coefficients_b.shape:
        raise ValueError(
            "Coefficient matrices must have equal shape."
        )

    n_realizations = coefficients_a.shape[0]

    if n_realizations == 0:
        raise ValueError(
            "At least one Fourier realization is required."
        )

    if effective_realizations is None:
        effective_realizations = float(
            n_realizations
        )

    if (
        not np.isfinite(
            effective_realizations
        )
        or effective_realizations <= 0.0
        or effective_realizations > n_realizations
    ):
        raise ValueError(
            "effective_realizations must be finite, positive, and "
            "not greater than n_realizations."
        )

    mean = (
        coefficients_a.T
        @ coefficients_b.conj()
        / n_realizations
    )

    variance = None

    if estimate_variance:
        if n_realizations < 2:
            raise ValueError(
                "At least two realizations are required to estimate "
                "cross-spectral variance."
            )

        power_a = np.abs(
            coefficients_a
        ) ** 2

        power_b = np.abs(
            coefficients_b
        ) ** 2

        mean_squared_magnitude = (
            power_a.T
            @ power_b
            / n_realizations
        )

        variance = (
            n_realizations
            / (
                n_realizations
                - 1
            )
            * (
                mean_squared_magnitude
                - np.abs(mean) ** 2
            )
        )

        variance = np.maximum(
            variance,
            0.0,
        )

    return SpectralMatrixEstimate(
        mean=mean,
        variance=variance,
        n_realizations=n_realizations,
        effective_realizations=float(
            effective_realizations
        ),
    )


def cross_spectrum(
    coefficients_a: ComplexArray,
    coefficients_b: ComplexArray | None = None,
) -> ComplexArray:
    """Return the stacked cross-spectral matrix mean."""
    return estimate_spectral_matrix(
        coefficients_a,
        coefficients_b,
    ).mean


def matrix_blocks(
    spectral: SpectralComponents,
    block_start: int,
    windows_per_block: int,
) -> SpectralMatrixBlocks:
    """Construct canonical cross-spectral matrices for one block."""
    if not isinstance(
        spectral,
        SpectralComponents,
    ):
        raise TypeError(
            "spectral must be a SpectralComponents instance."
        )

    if block_start < 0:
        raise ValueError(
            "block_start cannot be negative."
        )

    if windows_per_block < 1:
        raise ValueError(
            "windows_per_block must be at least one."
        )

    if (
        block_start
        + windows_per_block
        > spectral.n_windows
    ):
        raise ValueError(
            "Requested covariance block exceeds the available "
            "Fourier realizations."
        )

    block = slice(
        block_start,
        block_start + windows_per_block,
    )

    vertical = (
        None
        if spectral.vertical is None
        else spectral.vertical[block]
    )

    east = (
        None
        if spectral.east is None
        else spectral.east[block]
    )

    north = (
        None
        if spectral.north is None
        else spectral.north[block]
    )

    return SpectralMatrixBlocks(
        zz=(
            None
            if vertical is None
            else cross_spectrum(
                vertical
            )
        ),
        ee=(
            None
            if east is None
            else cross_spectrum(
                east
            )
        ),
        nn=(
            None
            if north is None
            else cross_spectrum(
                north
            )
        ),
        en=(
            None
            if east is None
            or north is None
            else cross_spectrum(
                east,
                north,
            )
        ),
        ez=(
            None
            if east is None
            or vertical is None
            else cross_spectrum(
                east,
                vertical,
            )
        ),
        nz=(
            None
            if north is None
            or vertical is None
            else cross_spectrum(
                north,
                vertical,
            )
        ),
    )


def _validated_waveform_matrix(
    data: FloatArray,
) -> FloatArray:
    """Return a validated two-dimensional waveform matrix."""
    data = np.asarray(
        data,
        dtype=float,
    )

    if data.ndim != 2:
        raise ValueError(
            "Input waveform data must be two-dimensional."
        )

    if not np.all(
        np.isfinite(data)
    ):
        raise ValueError(
            "Input waveform data must contain only finite values."
        )

    return data


def _validated_waveforms(
    waveforms: dict[str, FloatArray],
) -> dict[str, FloatArray]:
    """Validate canonical waveform-component matrices."""
    if not isinstance(
        waveforms,
        dict,
    ):
        raise TypeError(
            "waveforms must be a dictionary."
        )

    unknown = set(
        waveforms
    ) - {
        "vertical",
        "east",
        "north",
    }

    if unknown:
        names = ", ".join(
            sorted(
                unknown
            )
        )

        raise ValueError(
            f"Unsupported waveform components: {names}."
        )

    if not waveforms:
        raise ValueError(
            "At least one waveform component is required."
        )

    validated = {
        name: _validated_waveform_matrix(
            data
        )
        for name, data in waveforms.items()
    }

    shapes = {
        data.shape
        for data in validated.values()
    }

    if len(shapes) != 1:
        raise ValueError(
            "All waveform components must have equal shape."
        )

    return validated


def _optional_direct_component(
    data: FloatArray | None,
    frequency: float,
    delta: float,
    windows: FrequencyWindows,
) -> ComplexArray | None:
    """Return direct coefficients for one optional component."""
    if data is None:
        return None

    return fourier_coefficients(
        data,
        frequency,
        delta,
        windows,
    )


def _direct_spectral_slices(
    waveforms: dict[str, FloatArray],
    npts: int,
    delta: float,
    config: SpectralConfig,
) -> Iterator[SpectralSlice]:
    """Yield direct-transform spectral slices."""
    frequencies = config.frequencies.targets()

    for frequency in frequencies:
        frequency = float(
            frequency
        )

        windows = frequency_windows(
            npts,
            delta,
            frequency,
            config.window,
        )

        components = spectral_components(
            waveforms,
            frequency,
            delta,
            windows,
        )

        yield SpectralSlice(
            frequency=frequency,
            windows=windows,
            components=components,
        )


def _rfft_spectral_slices(
    waveforms: dict[str, FloatArray],
    npts: int,
    delta: float,
    config: SpectralConfig,
) -> Iterator[SpectralSlice]:
    """Yield fixed-window spectral slices from shared real FFTs."""
    reference_frequency = _rfft_reference_frequency(
        config.frequencies
    )

    windows = frequency_windows(
        npts,
        delta,
        reference_frequency,
        config.window,
    )

    padding = config.resolved_fft_padding

    plan = _resolve_rfft_plan(
        config.frequencies,
        windows.size,
        delta,
        config.fft_selection,
        padding,
    )

    arrays: dict[str, ComplexArray] = {}

    for name, data in waveforms.items():
        arrays[name] = rfft_coefficients(
            data,
            plan.frequencies,
            delta,
            windows,
            selection=config.fft_selection,
            padding=padding,
        )

    for index, frequency in enumerate(
        plan.frequencies
    ):
        yield SpectralSlice(
            frequency=float(
                frequency
            ),
            windows=windows,
            components=SpectralComponents(
                vertical=_slice_rfft_component(
                    arrays.get(
                        "vertical"
                    ),
                    index,
                ),
                east=_slice_rfft_component(
                    arrays.get(
                        "east"
                    ),
                    index,
                ),
                north=_slice_rfft_component(
                    arrays.get(
                        "north"
                    ),
                    index,
                ),
            ),
        )


def _slice_rfft_component(
    coefficients: ComplexArray | None,
    index: int,
) -> ComplexArray | None:
    """Return one frequency from optional multi-frequency coefficients."""
    if coefficients is None:
        return None

    return coefficients[
        :,
        index,
        :,
    ]


def _rfft_reference_frequency(
    config: FrequencyConfig,
) -> float:
    """Return a valid frequency used only to build fixed windows."""
    if config.values is not None:
        return float(
            np.asarray(
                config.values,
                dtype=float,
            )[0]
        )

    return float(
        config.minimum
    )


def _resolve_rfft_plan(
    config: FrequencyConfig,
    window_size: int,
    delta: float,
    selection: FFTSelection,
    padding: int,
) -> _RFFTPlan:
    """Resolve native-bin or interpolated RFFT analysis frequencies."""
    native = np.fft.rfftfreq(
        window_size * padding,
        d=delta,
    )

    if native.size < 2:
        raise ValueError(
            "RFFT requires at least two frequency bins."
        )

    nyquist = float(
        native[-1]
    )

    if config.values is None and config.size is None:
        mask = (
            native >= config.minimum
        ) & (
            native <= config.maximum
        )

        frequencies = native[
            mask
        ]

        frequencies = frequencies[
            frequencies > 0.0
        ]

        if frequencies.size == 0:
            raise ValueError(
                "No native RFFT bins fall inside the requested "
                "frequency range."
            )

        return _RFFTPlan(
            frequencies=frequencies.copy(),
        )

    targets = config.targets()

    if targets[-1] > nyquist:
        raise ValueError(
            f"Maximum requested frequency {targets[-1]:g} Hz exceeds "
            f"the Nyquist frequency {nyquist:g} Hz."
        )

    if selection is FFTSelection.INTERP:
        return _RFFTPlan(
            frequencies=targets.copy(),
            positions=(
                targets
                / (
                    native[1]
                    - native[0]
                )
            ),
        )

    indices = _nearest_bin_indices(
        native,
        targets,
    )

    indices = _unique_preserve_order(
        indices
    )

    frequencies = native[
        indices
    ]

    frequencies = frequencies[
        frequencies > 0.0
    ]

    if frequencies.size == 0:
        raise ValueError(
            "Nearest-bin selection resolved only to the zero-frequency "
            "RFFT bin."
        )

    return _RFFTPlan(
        frequencies=frequencies.copy(),
        indices=indices,
    )


def _nearest_bin_indices(
    native: FloatArray,
    targets: FloatArray,
) -> IntArray:
    """Return indices of nearest values on a sorted native grid."""
    positions = np.searchsorted(
        native,
        targets,
        side="left",
    )

    positions = np.clip(
        positions,
        1,
        native.size - 1,
    )

    left = positions - 1
    right = positions

    choose_right = (
        targets - native[left]
        > native[right] - targets
    )

    return np.where(
        choose_right,
        right,
        left,
    ).astype(
        np.int64
    )


def _unique_preserve_order(
    values: IntArray,
) -> IntArray:
    """Return unique integer values while preserving input order."""
    _, indices = np.unique(
        values,
        return_index=True,
    )

    return values[
        np.sort(
            indices
        )
    ]


def _complex_cubic_uniform(
    spectrum: ComplexArray,
    positions: FloatArray,
) -> ComplexArray:
    """
    Interpolate a complex RFFT on its uniform native-bin coordinate.

    A local four-point cubic Lagrange interpolant is used for interior
    targets. Linear complex interpolation is used within one bin of
    either boundary. Exact native-bin targets are returned without
    interpolation. In normal YArray use this interpolation operates on
    an adaptively zero-padded RFFT grid to limit complex interpolation
    error.

    Parameters
    ----------
    spectrum
        Complex RFFT values with shape ``(n_stations, n_bins)``.
    positions
        Fractional native-bin positions of requested frequencies.

    Returns
    -------
    numpy.ndarray
        Complex values with shape
        ``(n_stations, n_frequencies)``.
    """
    spectrum = np.asarray(
        spectrum,
        dtype=np.complex128,
    )

    positions = np.asarray(
        positions,
        dtype=float,
    )

    if spectrum.ndim != 2:
        raise ValueError(
            "RFFT spectrum must be two-dimensional."
        )

    if positions.ndim != 1:
        raise ValueError(
            "Interpolation positions must be one-dimensional."
        )

    maximum = spectrum.shape[1] - 1

    if np.any(
        positions < 0.0
    ) or np.any(
        positions > maximum
    ):
        raise ValueError(
            "Interpolation positions fall outside the RFFT grid."
        )

    output = np.empty(
        (
            spectrum.shape[0],
            positions.size,
        ),
        dtype=np.complex128,
    )

    rounded = np.rint(
        positions
    ).astype(
        np.int64
    )

    exact = np.isclose(
        positions,
        rounded,
        rtol=0.0,
        atol=32.0
        * np.finfo(float).eps
        * max(
            1.0,
            float(maximum),
        ),
    )

    if np.any(
        exact
    ):
        output[
            :,
            exact,
        ] = spectrum[
            :,
            rounded[exact],
        ]

    pending = np.flatnonzero(
        ~exact
    )

    if pending.size == 0:
        return output

    local_positions = positions[
        pending
    ]

    base = np.floor(
        local_positions
    ).astype(
        np.int64
    )

    interior = (
        base >= 1
    ) & (
        base + 2 <= maximum
    )

    if np.any(
        interior
    ):
        selected = pending[
            interior
        ]

        x = (
            positions[selected]
            - base[interior]
        )

        i0 = base[interior] - 1
        i1 = base[interior]
        i2 = base[interior] + 1
        i3 = base[interior] + 2

        w0 = (
            -x
            * (x - 1.0)
            * (x - 2.0)
            / 6.0
        )

        w1 = (
            (x + 1.0)
            * (x - 1.0)
            * (x - 2.0)
            / 2.0
        )

        w2 = (
            -(x + 1.0)
            * x
            * (x - 2.0)
            / 2.0
        )

        w3 = (
            (x + 1.0)
            * x
            * (x - 1.0)
            / 6.0
        )

        output[
            :,
            selected,
        ] = (
            spectrum[
                :,
                i0,
            ]
            * w0[None, :]
            + spectrum[
                :,
                i1,
            ]
            * w1[None, :]
            + spectrum[
                :,
                i2,
            ]
            * w2[None, :]
            + spectrum[
                :,
                i3,
            ]
            * w3[None, :]
        )

    boundary = ~interior

    if np.any(
        boundary
    ):
        selected = pending[
            boundary
        ]

        left = np.floor(
            positions[selected]
        ).astype(
            np.int64
        )

        left = np.clip(
            left,
            0,
            maximum - 1,
        )

        right = left + 1

        fraction = (
            positions[selected]
            - left
        )

        output[
            :,
            selected,
        ] = (
            spectrum[
                :,
                left,
            ]
            * (
                1.0
                - fraction[None, :]
            )
            + spectrum[
                :,
                right,
            ]
            * fraction[None, :]
        )

    return output
