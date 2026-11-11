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
Configuration of frequency-wavenumber processing.

:class:`BeamformingConfig` collects every parameter of the analysis:
analysis frequencies (:class:`FrequencyConfig`), windows and covariance
blocks (:class:`WindowConfig`), search grid (:class:`SearchGrid`),
ellipticity grid of the Rayleigh joint processing
(:class:`EllipticityGrid`), estimator, components, output and
uncertainties. Its parameters can be given singly, with the names of
the options of the command line (``fmin``, ``grid_max``,
``window_length``...), or grouped in the objects, which are immutable
(:meth:`replace` returns a modified copy); the configuration itself
changes with :meth:`BeamformingConfig.update`. A configuration is read
from and written to JSON: ``BeamformingConfig("capon.json")``,
:meth:`BeamformingConfig.read`, :meth:`BeamformingConfig.write`.

Options with a fixed set of values are enumerations; their string
values are accepted everywhere (``method="capon"`` as
``method=BeamformingMethod.CAPON``).

JSON format
-----------
Version 1, one object with the sections ``frequencies`` (``values``,
or ``minimum``, ``maximum``, ``size``, ``scale``), ``spectral``
(``transform``, ``fft_selection``, ``fft_padding``), ``method``
(``name``, ``diagonal_loading``, ``music_sources``), ``components``,
``grid`` (``quantity``, ``scale``, ``minimum``, ``maximum``, ``size``,
``azimuth_size`` or ``azimuths``), ``window`` (``type``, ``length``,
``taper``, ``overlap``, ``windows_per_block``, ``block_overlap``),
``ellipticity`` (``minimum_angle``, ``maximum_angle``, ``size``),
``output`` (``spectra``, ``maxima``, ``neighborhood``,
``minimum_relative_power``; the ``mode`` of the August files, maxima,
spectra or both, is read as ``spectra``), ``uncertainty`` (``statistical``,
``geometric``) and ``strict``. Missing optional fields take the
defaults of the classes; unknown sections are ignored, so that other
applications (e.g. the graphics of the command-line interface) can
keep their options in the same file.
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Optional

import numpy as np

from shakelab.exploration.yarray.estimators import BeamformingMethod
from shakelab.signals.windowing import window_size

__all__ = [
    "BeamformingConfig",
    "BeamformingMethod",
    "ComponentMode",
    "EllipticityGrid",
    "FFTSelection",
    "FrequencyConfig",
    "FrequencyScale",
    "GeometricUncertainty",
    "GridQuantity",
    "GridScale",
    "SearchGrid",
    "StatisticalUncertainty",
    "TransformMethod",
    "WindowConfig",
    "WindowType",
]

# Version of the JSON format.
_VERSION = 1


# ----------------------------------------------------------------------------
# Enumerations
# ----------------------------------------------------------------------------


class FrequencyScale(Enum):
    """Spacing of generated analysis frequencies."""

    LINEAR = "linear"
    LOG = "log"


class WindowType(Enum):
    """Length of the windows: fixed duration or number of cycles."""

    FIXED = "fixed"
    CYCLES = "cycles"


class TransformMethod(Enum):
    """Fourier transform: direct at any frequency, or one RFFT."""

    DIRECT = "direct"
    RFFT = "rfft"


class FFTSelection(Enum):
    """Frequencies from an RFFT: nearest bins or interpolated."""

    NEAREST = "nearest"
    INTERP = "interp"


class GridQuantity(Enum):
    """Quantity sampled along the radial axis of the search grid."""

    VELOCITY = "velocity"
    SLOWNESS = "slowness"
    WAVENUMBER = "wavenumber"


class GridScale(Enum):
    """Spacing of the radial axis of the search grid."""

    LINEAR = "linear"
    LOG = "log"


class ComponentMode(Enum):
    """Components processed."""

    VERTICAL = "vertical"
    RADIAL = "radial"
    TRANSVERSE = "transverse"
    THREE_COMPONENT_SEPARATE = "three_component_separate"
    RAYLEIGH_JOINT = "rayleigh_joint"


class StatisticalUncertainty(Enum):
    """Statistical uncertainty of the beam power."""

    NONE = "none"
    EMPIRICAL = "empirical"


class GeometricUncertainty(Enum):
    """Effect of the uncertainty of the station coordinates."""

    NONE = "none"
    EXPECTED_COHERENCE = "expected_coherence"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _enum(value: Any, kind: type, name: str) -> Enum:
    """Return an enumeration member from a member or its value."""
    if isinstance(value, kind):
        return value

    if isinstance(value, str):
        for member in kind:
            if member.value == value.strip().lower():
                return member

    names = ", ".join(repr(member.value) for member in kind)
    raise ValueError(f"{name} must be one of {names}, got {value!r}")


def _real(value: Any, name: str) -> float:
    """Return a finite real number."""
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a real number")

    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")

    return float(value)


def _integer(value: Any, name: str, minimum: int) -> int:
    """Return an integer not smaller than ``minimum``."""
    if (
        not isinstance(value, Integral)
        or isinstance(value, (bool, np.bool_))
    ):
        raise TypeError(f"{name} must be an integer")

    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")

    return int(value)


def _fraction(value: Any, name: str) -> float:
    """Return a fraction in ``[0, 1)``."""
    value = _real(value, name)

    if not 0.0 <= value < 1.0:
        raise ValueError(f"{name} must be in [0, 1)")

    return value


def _readonly(value: Any) -> np.ndarray:
    """Return a read-only float copy."""
    array = np.array(value, dtype=float)
    array.setflags(write=False)
    return array


def _set(instance: Any, **values: Any) -> None:
    """Set attributes of a frozen dataclass."""
    for name, value in values.items():
        object.__setattr__(instance, name, value)


class _Replace:
    """Mixin providing :meth:`replace` to frozen dataclasses."""

    def replace(self, **changes: Any) -> Any:
        """
        Return a copy with some fields changed.

        Parameters
        ----------
        **changes
            New values of the fields, by name.

        Returns
        -------
        object
            New validated instance of the same class.
        """
        return dataclasses.replace(self, **changes)


# ----------------------------------------------------------------------------
# Frequencies, windows and grids
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class FrequencyConfig(_Replace):
    """
    Analysis frequencies.

    Either explicit ``values``, or a range from ``minimum`` to
    ``maximum`` with ``size`` frequencies spaced as ``scale``. With the
    RFFT transform a range without ``size`` means every FFT bin in the
    range.

    Parameters
    ----------
    values : array-like, optional
        Explicit frequencies in Hz, positive and increasing.
    minimum, maximum : float, optional
        Range of the frequencies, in Hz.
    size : int, optional
        Number of frequencies of the range (at least two).
    scale : FrequencyScale or str, optional
        ``"log"`` (default) or ``"linear"`` spacing.

    Raises
    ------
    ValueError
        If values and range are mixed or incomplete, or invalid.
    """

    values: Any = None
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    size: Optional[int] = None
    scale: Any = FrequencyScale.LOG

    def __post_init__(self) -> None:
        _set(self, scale=_enum(self.scale, FrequencyScale, "scale"))

        if self.values is not None:
            if (self.minimum, self.maximum, self.size) != (None,) * 3:
                raise ValueError(
                    "give either frequency values or a range, not both"
                )

            values = np.atleast_1d(np.asarray(self.values, dtype=float))

            if (
                values.ndim != 1
                or not np.all(np.isfinite(values))
                or np.any(values <= 0.0)
                or np.any(np.diff(values) <= 0.0)
            ):
                raise ValueError(
                    "frequency values must be positive and increasing"
                )

            _set(self, values=_readonly(values))
            return

        if self.minimum is None or self.maximum is None:
            raise ValueError("a frequency range needs minimum and maximum")

        minimum = _real(self.minimum, "minimum")
        maximum = _real(self.maximum, "maximum")

        if not 0.0 < minimum < maximum:
            raise ValueError("the frequency range must be 0 < min < max")

        size = None if self.size is None else _integer(self.size, "size", 2)
        _set(self, minimum=minimum, maximum=maximum, size=size)

    @property
    def is_explicit(self) -> bool:
        """Whether the frequencies are given explicitly."""
        return self.values is not None

    @property
    def is_sampled(self) -> bool:
        """Whether the frequencies are defined without the FFT grid."""
        return self.values is not None or self.size is not None

    def targets(self) -> np.ndarray:
        """
        Return the requested frequencies.

        Returns
        -------
        numpy.ndarray
            Frequencies in Hz.

        Raises
        ------
        ValueError
            For a range without ``size`` (every FFT bin), whose
            frequencies depend on the window length.
        """
        if self.values is not None:
            return np.array(self.values)

        if self.size is None:
            raise ValueError(
                "a range without size selects every FFT bin: the "
                "frequencies depend on the window"
            )

        if self.scale is FrequencyScale.LINEAR:
            return np.linspace(self.minimum, self.maximum, self.size)

        return np.geomspace(self.minimum, self.maximum, self.size)

    def to_dict(self) -> dict:
        """
        Return the JSON representation.

        Returns
        -------
        dict
            Section ``frequencies`` of the JSON format.
        """
        if self.values is not None:
            return {"values": [float(value) for value in self.values]}

        return {
            "minimum": self.minimum,
            "maximum": self.maximum,
            "size": self.size,
            "scale": self.scale.value,
        }


@dataclass(frozen=True, eq=False)
class WindowConfig(_Replace):
    """
    Windows and covariance blocks.

    Each cross-spectral matrix averages ``windows_per_block``
    consecutive windows; consecutive blocks give the evolution in time.

    Parameters
    ----------
    type : WindowType or str, optional
        ``"cycles"`` (default): ``length`` is a number of cycles of the
        analysis frequency (windows shorter at higher frequency).
        ``"fixed"``: ``length`` is a duration in seconds.
    length : float, optional
        Number of cycles or duration. Default is 20.
    taper : float, optional
        Fraction of each window tapered by the Tukey window, in
        ``[0, 1]``. Default is 0.1.
    overlap : float, optional
        Fractional overlap between windows, in ``[0, 1)``. Default 0.5.
    windows_per_block : int, optional
        Windows averaged in each block. Default is 40.
    block_overlap : float, optional
        Fractional overlap between blocks, in windows. Default 0.5.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    type: Any = WindowType.CYCLES
    length: float = 20.0
    taper: float = 0.1
    overlap: float = 0.5
    windows_per_block: int = 40
    block_overlap: float = 0.5

    def __post_init__(self) -> None:
        length = _real(self.length, "length")

        if length <= 0.0:
            raise ValueError("the window length must be positive")

        taper = _real(self.taper, "taper")

        if not 0.0 <= taper <= 1.0:
            raise ValueError("taper must be in [0, 1]")

        _set(
            self,
            type=_enum(self.type, WindowType, "type"),
            length=length,
            taper=taper,
            overlap=_fraction(self.overlap, "overlap"),
            windows_per_block=_integer(
                self.windows_per_block, "windows_per_block", 1
            ),
            block_overlap=_fraction(self.block_overlap, "block_overlap"),
        )

    def size(self, delta: float, frequency: Optional[float] = None) -> int:
        """
        Return the number of samples of a window.

        Parameters
        ----------
        delta : float
            Sampling interval in seconds.
        frequency : float, optional
            Analysis frequency, required for windows in cycles.

        Returns
        -------
        int
            Number of samples.
        """
        if self.type is WindowType.FIXED:
            return window_size(delta, duration=self.length)

        if frequency is None:
            raise ValueError("windows in cycles need the frequency")

        return window_size(delta, cycles=self.length, frequency=frequency)

    def to_dict(self) -> dict:
        """
        Return the JSON representation.

        Returns
        -------
        dict
            Section ``window`` of the JSON format.
        """
        return {
            "type": self.type.value,
            "length": self.length,
            "taper": self.taper,
            "overlap": self.overlap,
            "windows_per_block": self.windows_per_block,
            "block_overlap": self.block_overlap,
        }


@dataclass(frozen=True, eq=False)
class SearchGrid(_Replace):
    """
    Polar grid of the wavevectors searched.

    Parameters
    ----------
    minimum, maximum : float
        Range of the radial quantity: velocity (m/s), slowness (s/m) or
        angular wavenumber (rad/m), with distances in the units of the
        station coordinates.
    size : int
        Number of radial values (at least two).
    quantity : GridQuantity or str, optional
        ``"velocity"`` (default), ``"slowness"`` or ``"wavenumber"``.
    scale : GridScale or str, optional
        ``"log"`` or ``"linear"`` spacing. Default: log for velocities
        and slownesses, linear for wavenumbers (the only spacing
        allowed for them).
    azimuth_size : int, optional
        Number of propagation azimuths, uniform in ``[0, 360)``.
        Default is 180.
    azimuths : array-like, optional
        Explicit propagation azimuths in degrees clockwise from north;
        replace ``azimuth_size``.

    Raises
    ------
    ValueError
        If a value is out of range.
    """

    minimum: float
    maximum: float
    size: int
    quantity: Any = GridQuantity.VELOCITY
    scale: Any = None
    azimuth_size: int = 180
    azimuths: Any = None

    def __post_init__(self) -> None:
        quantity = _enum(self.quantity, GridQuantity, "quantity")
        scale = self.scale

        if scale is None:
            scale = (
                GridScale.LINEAR if quantity is GridQuantity.WAVENUMBER
                else GridScale.LOG
            )

        scale = _enum(scale, GridScale, "scale")
        minimum = _real(self.minimum, "minimum")
        maximum = _real(self.maximum, "maximum")

        if minimum >= maximum:
            raise ValueError("the grid maximum must exceed the minimum")

        if quantity is GridQuantity.WAVENUMBER:
            if minimum < 0.0:
                raise ValueError("wavenumbers cannot be negative")

            if scale is GridScale.LOG:
                raise ValueError("wavenumbers are sampled linearly")

        elif minimum <= 0.0:
            raise ValueError("velocities and slownesses must be positive")

        azimuths = self.azimuths

        if azimuths is not None:
            azimuths = np.mod(np.asarray(azimuths, dtype=float), 360.0)

            if (
                azimuths.ndim != 1
                or azimuths.size == 0
                or not np.all(np.isfinite(azimuths))
                or np.unique(azimuths).size != azimuths.size
            ):
                raise ValueError("azimuths must be distinct finite values")

            azimuths = _readonly(azimuths)

        _set(
            self,
            minimum=minimum,
            maximum=maximum,
            size=_integer(self.size, "size", 2),
            quantity=quantity,
            scale=scale,
            azimuth_size=(
                azimuths.size if azimuths is not None
                else _integer(self.azimuth_size, "azimuth_size", 1)
            ),
            azimuths=azimuths,
        )

    @property
    def values(self) -> np.ndarray:
        """Radial values of the grid."""
        if self.scale is GridScale.LINEAR:
            return np.linspace(self.minimum, self.maximum, self.size)

        return np.geomspace(self.minimum, self.maximum, self.size)

    @property
    def azimuth_values(self) -> np.ndarray:
        """Propagation azimuths in degrees clockwise from north."""
        if self.azimuths is not None:
            return np.array(self.azimuths)

        return np.linspace(0.0, 360.0, self.azimuth_size, endpoint=False)

    def wavenumbers(self, frequency: float) -> np.ndarray:
        """
        Return the angular wavenumbers of the radial values.

        Parameters
        ----------
        frequency : float
            Frequency in Hz.

        Returns
        -------
        numpy.ndarray
            Angular wavenumbers ``2 pi f / c`` (rad per distance unit).
        """
        frequency = _real(frequency, "frequency")

        if frequency <= 0.0:
            raise ValueError("frequency must be positive")

        if self.quantity is GridQuantity.VELOCITY:
            return 2.0 * np.pi * frequency / self.values

        if self.quantity is GridQuantity.SLOWNESS:
            return 2.0 * np.pi * frequency * self.values

        return self.values

    def wavevectors(self, frequency: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the east and north wavenumbers of the grid.

        Parameters
        ----------
        frequency : float
            Frequency in Hz.

        Returns
        -------
        kx, ky : numpy.ndarray
            Arrays ``(n_azimuths, n_values)``: ``k sin(azimuth)`` and
            ``k cos(azimuth)``.
        """
        k = self.wavenumbers(frequency)
        azimuth = np.radians(self.azimuth_values)

        return np.sin(azimuth)[:, None] * k, np.cos(azimuth)[:, None] * k

    def to_dict(self) -> dict:
        """
        Return the JSON representation.

        Returns
        -------
        dict
            Section ``grid`` of the JSON format.
        """
        data = {
            "quantity": self.quantity.value,
            "scale": self.scale.value,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "size": self.size,
        }

        if self.azimuths is None:
            data["azimuth_size"] = self.azimuth_size
        else:
            data["azimuths"] = [float(value) for value in self.azimuths]

        return data


@dataclass(frozen=True, eq=False)
class EllipticityGrid(_Replace):
    """
    Signed Rayleigh ellipticities searched by the joint processing.

    The ellipticity ``e = tan(chi)`` is sampled uniformly in the angle
    ``chi``, strictly inside (-90, 90) degrees.

    Parameters
    ----------
    minimum_angle, maximum_angle : float, optional
        Range of ``chi`` in degrees. Default is -85 to 85.
    size : int, optional
        Number of angles (at least three). Default is 61.

    Raises
    ------
    ValueError
        If the range is not inside (-90, 90) or is reversed.
    """

    minimum_angle: float = -85.0
    maximum_angle: float = 85.0
    size: int = 61

    def __post_init__(self) -> None:
        lower = _real(self.minimum_angle, "minimum_angle")
        upper = _real(self.maximum_angle, "maximum_angle")

        if not -90.0 < lower < upper < 90.0:
            raise ValueError(
                "ellipticity angles must satisfy -90 < min < max < 90"
            )

        _set(
            self,
            minimum_angle=lower,
            maximum_angle=upper,
            size=_integer(self.size, "size", 3),
        )

    @property
    def angles(self) -> np.ndarray:
        """Ellipticity angles in degrees."""
        return np.linspace(self.minimum_angle, self.maximum_angle, self.size)

    @property
    def values(self) -> np.ndarray:
        """Signed ellipticities ``tan(angle)``."""
        return np.tan(np.radians(self.angles))

    def to_dict(self) -> dict:
        """
        Return the JSON representation.

        Returns
        -------
        dict
            Section ``ellipticity`` of the JSON format.
        """
        return {
            "minimum_angle": self.minimum_angle,
            "maximum_angle": self.maximum_angle,
            "size": self.size,
        }


# ----------------------------------------------------------------------------
# Complete configuration
# ----------------------------------------------------------------------------


# Objects of the configuration, their classes and JSON sections.
_OBJECTS = {
    "frequencies": (FrequencyConfig, "frequencies"),
    "grid": (SearchGrid, "grid"),
    "window": (WindowConfig, "window"),
    "ellipticity_grid": (EllipticityGrid, "ellipticity"),
}

# Single parameters of the objects: name -> (object, field). The names
# are those of the command-line options (``--grid-max`` -> grid_max).
_FIELDS = {
    "fmin": ("frequencies", "minimum"),
    "fmax": ("frequencies", "maximum"),
    "nfreq": ("frequencies", "size"),
    "frequency_scale": ("frequencies", "scale"),
    "grid_quantity": ("grid", "quantity"),
    "grid_scale": ("grid", "scale"),
    "grid_min": ("grid", "minimum"),
    "grid_max": ("grid", "maximum"),
    "grid_size": ("grid", "size"),
    "azimuth_size": ("grid", "azimuth_size"),
    "azimuths": ("grid", "azimuths"),
    "window_type": ("window", "type"),
    "window_length": ("window", "length"),
    "window_taper": ("window", "taper"),
    "window_overlap": ("window", "overlap"),
    "windows_per_block": ("window", "windows_per_block"),
    "block_overlap": ("window", "block_overlap"),
    "ellipticity_angle_min": ("ellipticity_grid", "minimum_angle"),
    "ellipticity_angle_max": ("ellipticity_grid", "maximum_angle"),
    "ellipticity_angle_size": ("ellipticity_grid", "size"),
}

# Other parameters: name -> (default, JSON section, JSON key); a None
# section is the top level.
_OPTIONS = {
    "method": (BeamformingMethod.CAPON, "method", "name"),
    "components": (ComponentMode.VERTICAL, None, "components"),
    "diagonal_loading": (0.001, "method", "diagonal_loading"),
    "music_sources": (1, "method", "music_sources"),
    "spectra": (False, "output", "spectra"),
    "maxima": (1, "output", "maxima"),
    "peak_neighborhood": ((3, 3), "output", "neighborhood"),
    "minimum_relative_power": (0.0, "output", "minimum_relative_power"),
    "statistical_uncertainty": (
        StatisticalUncertainty.NONE, "uncertainty", "statistical"
    ),
    "geometric_uncertainty": (
        GeometricUncertainty.NONE, "uncertainty", "geometric"
    ),
    "strict": (True, None, "strict"),
    "transform": (TransformMethod.DIRECT, "spectral", "transform"),
    "fft_selection": (FFTSelection.NEAREST, "spectral", "fft_selection"),
    "fft_padding": (None, "spectral", "fft_padding"),
}

_NAMES = frozenset(_OBJECTS) | frozenset(_FIELDS) | frozenset(_OPTIONS)


class BeamformingConfig:
    """
    Complete configuration of frequency-wavenumber processing.

    Every parameter can be given singly, with the names of the options
    of the command line (``--grid-max`` is ``grid_max``), or grouped in
    the objects :class:`FrequencyConfig`, :class:`SearchGrid`,
    :class:`WindowConfig` and :class:`EllipticityGrid` (also as
    dictionaries with the fields of the JSON sections). Single
    parameters apply on top of the objects and of the file.

    Parameters
    ----------
    source : str, path-like or mapping, optional
        JSON file, or its content as a dictionary (see the module
        notes). Default: the default values.
    **parameters
        Parameters (see Notes), replacing those of ``source``.

    Raises
    ------
    TypeError
        If a parameter is unknown or has the wrong type.
    ValueError
        If a value is invalid or the options are incompatible.

    Notes
    -----
    Parameters and defaults:

    - frequencies: ``frequencies`` (FrequencyConfig, dictionary or
      explicit values in Hz), or ``fmin``, ``fmax``, ``nfreq``,
      ``frequency_scale`` (``"log"``). Required.
    - search grid: ``grid`` (SearchGrid or dictionary), or
      ``grid_min``, ``grid_max``, ``grid_size``, ``grid_quantity``
      (``"velocity"``), ``grid_scale``, ``azimuth_size`` (180) or
      ``azimuths``. Required.
    - windows: ``window`` (WindowConfig or dictionary), or
      ``window_type`` (``"cycles"``), ``window_length`` (20),
      ``window_taper`` (0.1), ``window_overlap`` (0.5),
      ``windows_per_block`` (40), ``block_overlap`` (0.5).
    - ``method``: ``"capon"`` (default), ``"bartlett"`` or ``"music"``;
      ``diagonal_loading`` of Capon, relative to the mean power of the
      matrix (0.001); ``music_sources`` excluded from the MUSIC noise
      subspace (1).
    - ``components``: ``"vertical"`` (default), ``"radial"``,
      ``"transverse"``, ``"three_component_separate"`` or
      ``"rayleigh_joint"`` (Capon only); ellipticities of the joint
      processing: ``ellipticity_grid``, or ``ellipticity_angle_min``
      (-85), ``ellipticity_angle_max`` (85), ``ellipticity_angle_size``
      (61).
    - output: the maxima are always kept; ``spectra`` keeps also the
      complete beam-power spectra of every block (False: they can
      take much memory); ``maxima`` kept for each block, frequency and
      component (1), ``peak_neighborhood`` (azimuths, radial values) of
      a local maximum ((3, 3)), ``minimum_relative_power`` of a maximum
      relative to the largest of the spectrum (0).
    - ``statistical_uncertainty``: ``"none"`` or ``"empirical"``
      (Bartlett and Capon); ``geometric_uncertainty``: ``"none"`` or
      ``"expected_coherence"`` (Bartlett).
    - ``strict``: Capon and MUSIC require more windows per block than
      the dimension of the matrix (True).
    - ``transform``: ``"direct"`` or ``"rfft"`` (fixed windows only),
      ``fft_selection`` (``"nearest"`` or ``"interp"``) and
      ``fft_padding`` (1 for nearest, 8 for interp).

    Every parameter is also an attribute, and assigning it is the same
    as :meth:`update`.

    Examples
    --------
    >>> config = BeamformingConfig(fmin=3.0, fmax=10.0, nfreq=12,
    ...                            grid_min=100.0, grid_max=1500.0,
    ...                            grid_size=60, method="capon")
    >>> config = BeamformingConfig("capon.json", window_length=4.0)
    >>> config.update(method="music")
    >>> config.grid_max = 2000.0
    """

    def __init__(self, source: Any = None, **parameters: Any) -> None:
        object.__setattr__(self, "_state", _validated(_defaults()))

        if source is not None:
            self.read(source)

        if parameters:
            self.update(**parameters)

    def __getattr__(self, name: str) -> Any:
        state = self.__dict__.get("_state")

        if state is not None:
            if name in state:
                return state[name]

            if name in _FIELDS:
                owner, field_name = _FIELDS[name]
                value = state[owner]
                return None if value is None else getattr(value, field_name)

        raise AttributeError(
            f"'BeamformingConfig' object has no attribute {name!r}"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        if name not in _NAMES:
            raise AttributeError(f"unknown configuration parameter {name!r}")

        self.update(**{name: value})

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, BeamformingConfig):
            return NotImplemented

        return self.to_dict() == other.to_dict()

    __hash__ = None

    def __repr__(self) -> str:
        state = self._state
        return (
            f"BeamformingConfig({state['method'].value}, "
            f"{state['components'].value}, frequencies="
            f"{'set' if state['frequencies'] else 'missing'}, grid="
            f"{'set' if state['grid'] else 'missing'})"
        )

    # -- changes -----------------------------------------------------------

    def update(self, **parameters: Any) -> None:
        """
        Change parameters.

        The new values are validated together: if any is invalid, the
        configuration is not changed.

        Parameters
        ----------
        **parameters
            Parameters by name (see the class notes).

        Raises
        ------
        TypeError
            If a parameter is unknown or has the wrong type.
        ValueError
            If a value is invalid or the options are incompatible.
        """
        unknown = sorted(set(parameters) - _NAMES)

        if unknown:
            raise TypeError(
                "unknown configuration parameters: " + ", ".join(unknown)
                + "; accepted: " + ", ".join(sorted(_NAMES))
            )

        state = dict(self._state)

        for name, (kind, _) in _OBJECTS.items():
            if name in parameters:
                state[name] = _as_object(name, kind, parameters[name])

        changes: dict[str, dict] = {}

        for name, (owner, field_name) in _FIELDS.items():
            if name in parameters:
                changes.setdefault(owner, {})[field_name] = parameters[name]

        for owner, values in changes.items():
            state[owner] = _merged(owner, state[owner], values)

        for name in _OPTIONS:
            if name in parameters:
                state[name] = parameters[name]

        object.__setattr__(self, "_state", _validated(state))

    def copy(self) -> BeamformingConfig:
        """
        Return an independent copy.

        Returns
        -------
        BeamformingConfig
            Copy (the objects of the configuration are immutable and
            shared).
        """
        new = BeamformingConfig()
        object.__setattr__(new, "_state", dict(self._state))
        return new

    @property
    def complete(self) -> bool:
        """Whether frequencies and grid are given."""
        return (
            self._state["frequencies"] is not None
            and self._state["grid"] is not None
        )

    def check_complete(self) -> None:
        """
        Raise if the configuration cannot be used for processing.

        Raises
        ------
        ValueError
            If frequencies or grid are missing.
        """
        missing = [name for name in ("frequencies", "grid")
                   if self._state[name] is None]

        if missing:
            raise ValueError(
                "the configuration needs " + " and ".join(missing)
                + " (e.g. fmin, fmax, nfreq and grid_min, grid_max, "
                "grid_size)"
            )

    # -- JSON --------------------------------------------------------------

    def read(self, source: Any) -> None:
        """
        Replace the configuration with that of a JSON file.

        Missing fields take the default values; unknown sections are
        ignored (e.g. ``graphics`` of the command line).

        Parameters
        ----------
        source : str, path-like or mapping
            JSON file, or its content as a dictionary.

        Raises
        ------
        TypeError
            If a section has the wrong type.
        ValueError
            If the file is not valid JSON, the version is not supported
            or a value is invalid.
        """
        if isinstance(source, Mapping):
            data = source
        else:
            path = Path(source).expanduser()

            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path}: invalid JSON at line {error.lineno}, column "
                    f"{error.colno}: {error.msg}"
                ) from error

        parameters = _from_json(data)
        object.__setattr__(self, "_state", _validated(_defaults()))
        self.update(**parameters)

    def write(self, target: Any, owrite: bool = False) -> None:
        """
        Write the configuration to a JSON file.

        Parameters
        ----------
        target : str or path-like
            JSON file.
        owrite : bool, optional
            Overwrite an existing file. Default is False.

        Raises
        ------
        FileExistsError
            If the file exists and ``owrite`` is False.
        """
        path = Path(target).expanduser()

        if path.exists() and not owrite:
            raise FileExistsError(f"{path} exists (use owrite=True)")

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8"
        )

    def to_dict(self) -> dict:
        """
        Return the JSON representation.

        Returns
        -------
        dict
            Content of a JSON configuration (version 1); frequencies
            and grid only if given.
        """
        state = self._state
        data: dict = {"version": _VERSION}

        for name, (_, section) in _OBJECTS.items():
            if state[name] is not None:
                data[section] = state[name].to_dict()

        for name, (_, section, key) in _OPTIONS.items():
            value = state[name]

            if isinstance(value, Enum):
                value = value.value
            elif isinstance(value, tuple):
                value = list(value)

            if section is None:
                data[key] = value
            else:
                data.setdefault(section, {})[key] = value

        return data

    def info(self) -> None:
        """Print the configuration."""
        print(self.summary())

    def summary(self) -> str:
        """
        Return the configuration as text.

        Returns
        -------
        str
            The JSON representation, indented.
        """
        return "BeamformingConfig\n" + json.dumps(self.to_dict(), indent=2)


# ----------------------------------------------------------------------------
# Configuration helpers
# ----------------------------------------------------------------------------


class _OutputMode(Enum):
    """Output modes of the JSON files of August 2026 (``output.mode``)."""

    MAXIMA = "maxima"
    SPECTRA = "spectra"
    BOTH = "both"


def _defaults() -> dict:
    """State of a configuration with the default values."""
    state = {name: None for name in ("frequencies", "grid")}
    state["window"] = WindowConfig()
    state["ellipticity_grid"] = EllipticityGrid()
    state.update({name: value[0] for name, value in _OPTIONS.items()})
    return state


def _as_object(name: str, kind: type, value: Any) -> Any:
    """Object of the configuration from an object, mapping or values."""
    if value is None and name in ("frequencies", "grid"):
        return None

    if isinstance(value, kind):
        return value

    if isinstance(value, Mapping):
        return kind(**value)

    if name == "frequencies":
        return FrequencyConfig(values=value)

    raise TypeError(f"{name} must be a {kind.__name__} or a dictionary")


def _merged(owner: str, current: Any, values: dict) -> Any:
    """Object with some fields replaced."""
    kind = _OBJECTS[owner][0]
    base = {} if current is None else {
        item.name: getattr(current, item.name)
        for item in dataclasses.fields(current)
    }

    if owner == "frequencies" and set(values) & {
        "minimum", "maximum", "size", "scale"
    }:
        base.pop("values", None)

    if owner == "grid":
        if "azimuth_size" in values:
            base.pop("azimuths", None)

        if "quantity" in values and "scale" not in values:
            base.pop("scale", None)

    return kind(**{**base, **values})


def _validated(state: dict) -> dict:
    """Return a validated state (enumerations, numbers, compatibility)."""
    state = dict(state)

    for name in ("strict", "spectra"):
        if not isinstance(state[name], (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")

    state["spectra"] = bool(state["spectra"])

    neighborhood = tuple(state["peak_neighborhood"])

    if len(neighborhood) != 2:
        raise ValueError("peak_neighborhood must contain two integers")

    loading = _real(state["diagonal_loading"], "diagonal_loading")

    if loading < 0.0:
        raise ValueError("diagonal_loading cannot be negative")

    relative = _real(state["minimum_relative_power"],
                     "minimum_relative_power")

    if not 0.0 <= relative <= 1.0:
        raise ValueError("minimum_relative_power must be in [0, 1]")

    padding = state["fft_padding"]

    if padding is not None:
        padding = _integer(padding, "fft_padding", 1)

    state.update(
        method=_enum(state["method"], BeamformingMethod, "method"),
        components=_enum(state["components"], ComponentMode, "components"),
        diagonal_loading=loading,
        music_sources=_integer(state["music_sources"], "music_sources", 1),
        maxima=_integer(state["maxima"], "maxima", 1),
        peak_neighborhood=tuple(
            _integer(value, "peak_neighborhood", 1) for value in neighborhood
        ),
        minimum_relative_power=relative,
        statistical_uncertainty=_enum(
            state["statistical_uncertainty"], StatisticalUncertainty,
            "statistical_uncertainty",
        ),
        geometric_uncertainty=_enum(
            state["geometric_uncertainty"], GeometricUncertainty,
            "geometric_uncertainty",
        ),
        transform=_enum(state["transform"], TransformMethod, "transform"),
        fft_selection=_enum(
            state["fft_selection"], FFTSelection, "fft_selection"
        ),
        fft_padding=padding,
    )
    _check_compatibility(state)
    return state


def _check_compatibility(state: dict) -> None:
    """Raise for combinations of options that are not supported."""
    method = state["method"]
    joint = state["components"] is ComponentMode.RAYLEIGH_JOINT
    frequencies = state["frequencies"]

    if state["transform"] is TransformMethod.RFFT:
        if state["window"].type is not WindowType.FIXED:
            raise ValueError("the RFFT needs windows of fixed duration")
    elif frequencies is not None and not frequencies.is_sampled:
        raise ValueError(
            "the direct transform needs explicit frequencies or a range "
            "with size"
        )

    if joint and method is not BeamformingMethod.CAPON:
        raise ValueError("the Rayleigh joint processing needs Capon")

    if state["statistical_uncertainty"] is StatisticalUncertainty.EMPIRICAL:
        if method is BeamformingMethod.MUSIC or joint:
            raise ValueError(
                "the empirical uncertainty is available for Bartlett and "
                "Capon on single components"
            )

    if state["geometric_uncertainty"] is not GeometricUncertainty.NONE:
        if method is not BeamformingMethod.BARTLETT:
            raise ValueError(
                "the coordinate uncertainty is available for Bartlett"
            )


def _from_json(data: Any) -> dict:
    """Parameters of a JSON configuration."""
    if not isinstance(data, Mapping):
        raise TypeError("the configuration must be a JSON object")

    version = data.get("version", _VERSION)

    if version != _VERSION:
        raise ValueError(f"unsupported configuration version {version}")

    data = dict(data)

    if isinstance(data.get("method"), str):
        data["method"] = {"name": data["method"]}

    output = data.get("output")

    if isinstance(output, Mapping) and "mode" in output:
        output = dict(output)
        mode = _enum(output.pop("mode"), _OutputMode, "output.mode")
        output.setdefault("spectra", mode is not _OutputMode.MAXIMA)
        data["output"] = output

    parameters = {}

    for name, (_, section) in _OBJECTS.items():
        if section in data:
            parameters[name] = _section(data, section)

    for name, (_, section, key) in _OPTIONS.items():
        values = data if section is None else _section(data, section)

        if key in values:
            parameters[name] = values[key]

    return parameters


def _section(data: Mapping, name: str) -> dict:
    """Return an optional JSON section as a dictionary."""
    value = data.get(name, {})

    if not isinstance(value, Mapping):
        raise TypeError(f"the section {name!r} must be an object")

    return dict(value)
