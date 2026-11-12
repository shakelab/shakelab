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
Frequency-wavenumber beamforming of seismic arrays.

:class:`Beamformer` applies a :class:`~.config.BeamformingConfig` to
an :class:`~.data.ArrayData`, on the window common to its stations or
on a given one, and returns a :class:`BeamformingResult`: the local
maxima of the beam power (:class:`SpectralMaximum`) and, optionally,
the complete spectra (:class:`FrequencySpectrum`), for every analysis
frequency, covariance block and processed component. Results of
successive windows are joined with :meth:`BeamformingResult.append`;
:meth:`BeamformingResult.dispersion` and
:meth:`BeamformingResult.density` summarize the maxima.

Processing
----------
For each analysis frequency the motion of the stations along the
required directions (vertical, east, north) is divided into windows
(:mod:`shakelab.signals.windowing`), transformed
(:mod:`shakelab.signals.crossspectral`) and averaged by blocks of
windows into cross-spectral matrices; the beam power of each block is
evaluated on the polar search grid by the estimators of
:mod:`.estimators`:

- ``vertical``: matrix of the vertical motion;
- ``radial`` and ``transverse``: for each azimuth, matrix of the
  horizontal motion along the propagation direction or normal to it
  (azimuth + 90 degrees);
- ``three_component_separate``: the three above, separately;
- ``rayleigh_joint``: high-resolution three-component beamforming of
  Rayleigh waves (RTBF, Wathelet et al., 2018): radial and vertical
  motion in one matrix of dimension 2n, steering vectors
  ``[-1j * e * a, a]`` with the signed ellipticity ``e``, objective
  ``e**2 * P**2`` maximized over the ellipticity grid (Capon only).

Signed ellipticity
------------------
``e`` is the ratio of the radial to the vertical amplitude, with the
radial component positive in the direction of propagation and the
vertical one positive upwards. With the Fourier kernel
``exp(-2j pi f t)``, ``e > 0`` is radial motion ``e sin(wt)`` with
vertical ``cos(wt)``: **positive for prograde, negative for
retrograde** particle motion, as in Wathelet et al. (2018). Checked on
the fundamental mode of a homogeneous half-space (retrograde at the
surface, ``e = -0.68`` for a Poisson solid). Data whose vertical
component is positive downwards, or whose horizontal components are
reversed, give the opposite sign.

Units
-----
The matrices are cross power spectral densities (one-sided), so that
Bartlett and Capon give the density of a plane wave, in the units of
the data squared per Hz: physical units for data corrected for the
instrument, counts squared per Hz otherwise
(:attr:`BeamformingResult.power_units`). MUSIC gives a pseudo-power
without units; the RTBF objective has the square of the units of the
power.

References
----------
Wathelet, M., Guillier, B., Roux, P., Cornou, C., and Ohrnberger, M.
(2018). Rayleigh wave three-component beamforming: signed ellipticity
assessment from high-resolution frequency-wavenumber processing of
ambient vibration arrays. Geophysical Journal International, 215(1),
507-523.
"""

from __future__ import annotations

import json
import math
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
from scipy.ndimage import maximum_filter

from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    ComponentMode,
    FFTSelection,
    GeometricUncertainty,
    GridQuantity,
    StatisticalUncertainty,
    TransformMethod,
)
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.estimators import (
    BeamformingMethod,
    bartlett_power_statistics,
    beam_power,
    capon_power_statistics,
    check_realizations,
    regularized_matrix,
    steering_vectors,
)
from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import warn_external
from shakelab.signals import fourier
from shakelab.signals.crossspectral import coefficients, cross_spectral_matrix
from shakelab.signals.windowing import (
    effective_realizations,
    taper,
    window_starts,
)

__all__ = [
    "BeamformingResult",
    "BeamformingWarning",
    "Beamformer",
    "DispersionCurve",
    "FrequencySpectrum",
    "SpectralMaximum",
]

# Directions (east, north, up) of the motion used by each mode.
_DIRECTIONS = {
    "east": (1.0, 0.0, 0.0),
    "north": (0.0, 1.0, 0.0),
    "vertical": (0.0, 0.0, 1.0),
}

_MODE_INPUTS = {
    ComponentMode.VERTICAL: ("vertical",),
    ComponentMode.RADIAL: ("east", "north"),
    ComponentMode.TRANSVERSE: ("east", "north"),
    ComponentMode.THREE_COMPONENT_SEPARATE: ("east", "north", "vertical"),
    ComponentMode.RAYLEIGH_JOINT: ("east", "north", "vertical"),
}

_MODE_OUTPUTS = {
    ComponentMode.VERTICAL: ("vertical",),
    ComponentMode.RADIAL: ("radial",),
    ComponentMode.TRANSVERSE: ("transverse",),
    ComponentMode.THREE_COMPONENT_SEPARATE: (
        "vertical", "radial", "transverse",
    ),
    ComponentMode.RAYLEIGH_JOINT: ("rayleigh_joint",),
}


class BeamformingWarning(UserWarning):
    """Frequencies without complete covariance blocks."""


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class SpectralMaximum:
    """
    Local maximum of a beam-power spectrum.

    Attributes
    ----------
    block : int
        Index of the covariance block at this frequency.
    time : Date
        Centre of the block.
    frequency : float
        Frequency in Hz.
    component : str
        Processed component (``"vertical"``, ``"radial"``,
        ``"transverse"``, ``"rayleigh_joint"``).
    rank : int
        Rank by decreasing power in the block (1 is the largest).
    wavenumber : float
        Angular wavenumber (rad per distance unit).
    azimuth : float
        Propagation azimuth, degrees clockwise from north.
    power : float
        Beam power (see :attr:`BeamformingResult.power_units`).
    power_standard_error : float, optional
        Standard error of the power.
    ellipticity : float, optional
        Signed Rayleigh ellipticity, radial over vertical amplitude
        (joint processing; positive prograde, negative retrograde).
    ellipticity_angle : float, optional
        Ellipticity angle ``atan(ellipticity)`` in degrees.
    """

    block: int
    time: Date
    frequency: float
    component: str
    rank: int
    wavenumber: float
    azimuth: float
    power: float
    power_standard_error: Optional[float] = None
    ellipticity: Optional[float] = None
    ellipticity_angle: Optional[float] = None

    @property
    def velocity(self) -> float:
        """Phase velocity ``2 pi f / k``."""
        if self.wavenumber == 0.0:
            return math.inf

        return 2.0 * math.pi * self.frequency / self.wavenumber

    @property
    def slowness(self) -> float:
        """Phase slowness ``k / (2 pi f)``."""
        return self.wavenumber / (2.0 * math.pi * self.frequency)

    @property
    def kx(self) -> float:
        """East wavenumber."""
        return self.wavenumber * math.sin(math.radians(self.azimuth))

    @property
    def ky(self) -> float:
        """North wavenumber."""
        return self.wavenumber * math.cos(math.radians(self.azimuth))


@dataclass(frozen=True, eq=False)
class FrequencySpectrum:
    """
    Beam-power spectra of one frequency and component.

    Attributes
    ----------
    frequency : float
        Frequency in Hz.
    component : str
        Processed component.
    times : tuple of Date
        Centres of the covariance blocks.
    azimuths : numpy.ndarray
        Propagation azimuths of the grid, in degrees.
    grid_values : numpy.ndarray
        Radial values of the grid (see ``quantity``).
    quantity : GridQuantity
        Quantity of ``grid_values``.
    n_realizations : int
        Windows averaged in each block.
    effective_realizations : float
        Effective number of independent windows of each block.
    power : numpy.ndarray
        Beam power ``(n_blocks, n_azimuths, n_values)``.
    power_standard_error : numpy.ndarray, optional
        Standard error of the power, same shape.
    ellipticity, ellipticity_angle : numpy.ndarray, optional
        Signed ellipticity and its angle at each grid point (joint
        processing).
    """

    frequency: float
    component: str
    times: tuple
    azimuths: np.ndarray
    grid_values: np.ndarray
    quantity: GridQuantity
    n_realizations: int
    effective_realizations: float
    power: np.ndarray
    power_standard_error: Optional[np.ndarray] = None
    ellipticity: Optional[np.ndarray] = None
    ellipticity_angle: Optional[np.ndarray] = None

    @property
    def nblocks(self) -> int:
        """Number of covariance blocks."""
        return self.power.shape[0]

    def wavenumbers(self) -> np.ndarray:
        """
        Return the angular wavenumbers of the radial values.

        Returns
        -------
        numpy.ndarray
            Wavenumbers (rad per distance unit).
        """
        if self.quantity is GridQuantity.WAVENUMBER:
            return np.array(self.grid_values)

        if self.quantity is GridQuantity.VELOCITY:
            return 2.0 * np.pi * self.frequency / self.grid_values

        return 2.0 * np.pi * self.frequency * self.grid_values

    def wavevectors(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the east and north wavenumbers of the grid.

        Returns
        -------
        kx, ky : numpy.ndarray
            Arrays ``(n_azimuths, n_values)``.
        """
        k = self.wavenumbers()
        azimuth = np.radians(self.azimuths)
        return np.sin(azimuth)[:, None] * k, np.cos(azimuth)[:, None] * k


@dataclass(frozen=True, eq=False)
class DispersionCurve:
    """
    Statistics of the spectral maxima of one rank at each frequency.

    Attributes
    ----------
    component : str
        Processed component.
    rank : int
        Rank of the maxima.
    frequency : numpy.ndarray
        Frequencies in Hz with at least one maximum.
    velocity : numpy.ndarray
        Median phase velocity of the blocks, in m/s (distance units of
        the coordinates per second).
    velocity_low, velocity_high : numpy.ndarray
        16th and 84th percentiles of the velocity.
    azimuth : numpy.ndarray
        Circular mean of the propagation azimuths, in degrees.
    count : numpy.ndarray
        Number of maxima (blocks) at each frequency.
    ellipticity : numpy.ndarray, optional
        Median signed ellipticity (Rayleigh joint processing;
        positive prograde, negative retrograde).
    """

    component: str
    rank: int
    frequency: np.ndarray
    velocity: np.ndarray
    velocity_low: np.ndarray
    velocity_high: np.ndarray
    azimuth: np.ndarray
    count: np.ndarray
    ellipticity: Optional[np.ndarray] = None

    def info(self) -> None:
        """Print the curve as a table."""
        print(self.summary())

    def summary(self) -> str:
        """
        Return the curve as a text table.

        Returns
        -------
        str
            One line for each frequency.
        """
        lines = [
            f"DispersionCurve: {self.component}, rank {self.rank}",
            "  f (Hz)    v (m/s)   v16 (m/s)  v84 (m/s)  azimuth  blocks"
            + ("  ellipticity" if self.ellipticity is not None else ""),
        ]

        for index, frequency in enumerate(self.frequency):
            line = (
                f"  {frequency:7.3f}  {self.velocity[index]:9.1f}  "
                f"{self.velocity_low[index]:9.1f}  "
                f"{self.velocity_high[index]:9.1f}  "
                f"{self.azimuth[index]:7.1f}  {self.count[index]:6d}"
            )

            if self.ellipticity is not None:
                line += f"  {self.ellipticity[index]:11.3f}"

            lines.append(line)

        return "\n".join(lines)


class BeamformingResult:
    """
    Result of frequency-wavenumber processing.

    A result holds the maxima (and optionally the spectra) of the
    covariance blocks of one or more analyses with the same
    configuration: :meth:`append` adds those of a later analysis,
    :meth:`write` with ``append=True`` adds them to a file, which is
    read back joined.

    Parameters
    ----------
    source : str or path-like, optional
        Result file to read (see :meth:`read`). Default: an empty
        result.

    Attributes
    ----------
    maxima : list of SpectralMaximum
        Local maxima.
    spectra : list of FrequencySpectrum or None
        Complete spectra (None unless the configuration has
        ``spectra=True``).
    config : BeamformingConfig or None
        Configuration used (None for an empty result).
    block_times : dict
        Centres of the covariance blocks, by frequency.
    units : str or None
        Units of the data (None if unknown).

    Examples
    --------
    >>> result = beamformer.process()
    >>> result.write("result.pkl")
    >>> result = BeamformingResult("result.pkl")
    >>> result.append(beamformer.process(start, end))
    """

    def __init__(self, source: Any = None) -> None:
        self.maxima: list = []
        self.spectra: Optional[list] = None
        self.config: Optional[BeamformingConfig] = None
        self.block_times: dict = {}
        self.units: Optional[str] = None
        self._block_stations: dict = {}
        self._station_sets: list = []

        if source is not None:
            self.read(source)

    def __repr__(self) -> str:
        if self.config is None:
            return "BeamformingResult(empty)"

        return (
            f"BeamformingResult({self.config.method.value}, "
            f"{self.config.components.value}, {self.nblocks} blocks, "
            f"{len(self.maxima)} maxima)"
        )

    # -- properties -----------------------------------------------------

    @property
    def frequencies(self) -> np.ndarray:
        """Analysis frequencies, in Hz."""
        return np.array(sorted(self.block_times))

    @property
    def components(self) -> tuple[str, ...]:
        """Processed components."""
        if self.config is None:
            return ()

        return _MODE_OUTPUTS[self.config.components]

    @property
    def nblocks(self) -> int:
        """Number of covariance blocks, summed over the frequencies."""
        return sum(len(times) for times in self.block_times.values())

    @property
    def sids(self) -> tuple[str, ...]:
        """Stations processed in at least one block."""
        sids: dict = {}

        for names, _ in self._station_sets:
            sids.update(dict.fromkeys(names))

        return tuple(sids)

    @property
    def coordinates(self) -> np.ndarray:
        """Coordinates ``(n_stations, 3)`` of :attr:`sids`."""
        positions: dict = {}

        for names, coordinates in self._station_sets:
            for name, row in zip(names, coordinates):
                positions.setdefault(name, row)

        if not positions:
            return np.empty((0, 3))

        return np.array([positions[name] for name in self.sids])

    @property
    def power_units(self) -> Optional[str]:
        """Units of the power, None for the MUSIC pseudo-power."""
        if self.config is None or (
            self.config.method is BeamformingMethod.MUSIC
        ):
            return None

        units = f"({self.units or 'unit'})**2/Hz"

        if self.config.components is ComponentMode.RAYLEIGH_JOINT:
            return f"({units})**2"

        return units

    # -- selection ------------------------------------------------------

    def maxima_at(
        self,
        rank: Optional[int] = None,
        component: Optional[str] = None,
        frequency: Optional[float] = None,
    ) -> list[SpectralMaximum]:
        """
        Return the maxima of a rank, component or frequency.

        Parameters
        ----------
        rank : int, optional
            Rank (1 is the largest maximum of each block).
        component : str, optional
            Component.
        frequency : float, optional
            Frequency in Hz (the nearest analysis frequency).

        Returns
        -------
        list of SpectralMaximum
            Selected maxima, in the order of the result.
        """
        if frequency is not None and self.block_times:
            frequency = self._nearest_frequency(frequency)

        return [
            item for item in self.maxima
            if (rank is None or item.rank == rank)
            and (component is None or item.component == component)
            and (frequency is None or item.frequency == frequency)
        ]

    def spectrum(
        self,
        frequency: float,
        component: Optional[str] = None,
    ) -> FrequencySpectrum:
        """
        Return the spectra of the nearest analysis frequency.

        Parameters
        ----------
        frequency : float
            Frequency in Hz.
        component : str, optional
            Component; default is the first processed.

        Returns
        -------
        FrequencySpectrum
            Spectra of all blocks.

        Raises
        ------
        ValueError
            If the result has no spectra.
        KeyError
            If the component was not processed.
        """
        if not self.spectra:
            raise ValueError(
                "the result has no spectra (configuration with "
                "spectra=False)"
            )

        frequency = self._nearest_frequency(frequency)
        component = component or self.components[0]

        for item in self.spectra:
            if item.frequency == frequency and item.component == component:
                return item

        raise KeyError(component)

    def stations_at(self, frequency: float, block: int) -> tuple[str, ...]:
        """
        Return the stations processed in a covariance block.

        Parameters
        ----------
        frequency : float
            Frequency in Hz (the nearest analysis frequency).
        block : int
            Index of the block at that frequency (as
            :attr:`SpectralMaximum.block`).

        Returns
        -------
        tuple of str
            Station identifiers.
        """
        frequency = self._nearest_frequency(frequency)
        return self._station_sets[self._block_stations[frequency][block]][0]

    # -- analysis -------------------------------------------------------

    def dispersion(
        self,
        rank: int = 1,
        component: Optional[str] = None,
    ) -> DispersionCurve:
        """
        Return the statistics of the maxima of a rank at each frequency.

        Parameters
        ----------
        rank : int, optional
            Rank of the maxima. Default 1 (the largest of each block).
        component : str, optional
            Component; required if several were processed.

        Returns
        -------
        DispersionCurve
            Median velocity, percentiles, mean azimuth and number of
            blocks at each frequency with maxima.

        Raises
        ------
        ValueError
            If the component is ambiguous or not processed, or there
            are no maxima of that rank.
        """
        component = self._component(component)
        groups: dict = {}

        for item in self.maxima_at(rank=rank, component=component):
            groups.setdefault(item.frequency, []).append(item)

        if not groups:
            raise ValueError(f"no maximum of rank {rank} for {component}")

        frequencies = sorted(groups)
        columns: dict = {name: [] for name in (
            "velocity", "velocity_low", "velocity_high", "azimuth",
            "count", "ellipticity",
        )}

        for frequency in frequencies:
            items = groups[frequency]
            velocity = np.array([item.velocity for item in items])
            low, median, high = np.percentile(velocity, (16.0, 50.0, 84.0))
            angles = np.radians([item.azimuth for item in items])
            columns["velocity"].append(median)
            columns["velocity_low"].append(low)
            columns["velocity_high"].append(high)
            columns["azimuth"].append(np.degrees(np.arctan2(
                np.mean(np.sin(angles)), np.mean(np.cos(angles))
            )) % 360.0)
            columns["count"].append(len(items))
            values = [item.ellipticity for item in items
                      if item.ellipticity is not None]
            columns["ellipticity"].append(
                float(np.median(values)) if values else np.nan
            )

        ellipticity = np.array(columns.pop("ellipticity"))

        return DispersionCurve(
            component=component,
            rank=int(rank),
            frequency=np.array(frequencies),
            ellipticity=(
                None if np.all(np.isnan(ellipticity)) else ellipticity
            ),
            **{
                name: np.array(values, dtype=int if name == "count"
                               else float)
                for name, values in columns.items()
            },
        )

    def density(self, **options: Any) -> Any:
        """
        Return the density of the maxima at each frequency.

        Parameters
        ----------
        **options
            Options of
            :func:`~shakelab.exploration.yarray.statistics.maximum_density`
            (``component``, ``ranks``, ``coordinate``, ``weighting``,
            ``bandwidth_mode``...).

        Returns
        -------
        MaximumDensity
            Density, with its summary and peaks.
        """
        from shakelab.exploration.yarray.statistics import maximum_density

        return maximum_density(self, **options)

    # -- joining and files ----------------------------------------------

    def append(self, other: Any) -> None:
        """
        Add the blocks of other results.

        Blocks are renumbered; maxima and spectra keep the order of the
        results. The stations of each block are kept
        (:meth:`stations_at`).

        Parameters
        ----------
        other : BeamformingResult or iterable of BeamformingResult
            Results with the same configuration, frequencies and units.

        Raises
        ------
        TypeError
            If an element is not a BeamformingResult.
        ValueError
            If the results are not compatible (nothing is added).
        """
        items = (
            [other] if isinstance(other, BeamformingResult) else list(other)
        )

        for item in items:
            if not isinstance(item, BeamformingResult):
                raise TypeError("only BeamformingResult can be appended")

        items = [item for item in items if item.config is not None]
        reference = self if self.config is not None else (
            items[0] if items else None
        )

        for item in items:
            _check_compatible(reference, item)

        for item in items:
            if item is self:
                item = item.copy()

            self._join(item)

    def copy(self) -> BeamformingResult:
        """
        Return a copy of the result.

        Returns
        -------
        BeamformingResult
            New result; maxima and spectra are immutable and shared.
        """
        new = BeamformingResult()
        new._join(self)
        return new

    def read(self, source: Any) -> None:
        """
        Add the results of a file.

        Parameters
        ----------
        source : str or path-like
            File written by :meth:`write` (trusted origin only: pickle).

        Raises
        ------
        ValueError
            If the file is not a result file or its pieces are not
            compatible with each other or with this result.
        """
        from shakelab.exploration.yarray.io import read_result_states

        self.append([
            BeamformingResult._from_state(state)
            for state in read_result_states(source)
        ])

    def write(
        self,
        target: Any,
        append: bool = False,
        owrite: bool = False,
    ) -> None:
        """
        Write the result to a file (pickle).

        Parameters
        ----------
        target : str or path-like
            File; folders are created.
        append : bool, optional
            Add the result to the file (created if missing), to
            accumulate successive analyses; the file must hold results
            with the same configuration. Default False.
        owrite : bool, optional
            Overwrite an existing file (ignored with ``append``).
            Default False.

        Raises
        ------
        ValueError
            If the result is empty, or ``append`` and the file holds
            results that are not compatible.
        FileExistsError
            If the file exists and neither ``append`` nor ``owrite`` is
            set.
        """
        from shakelab.exploration.yarray.io import (
            read_result_states,
            write_result_state,
        )

        if self.config is None:
            raise ValueError("an empty result cannot be written")

        path = Path(target).expanduser()

        if append and path.exists() and path.stat().st_size > 0:
            first = BeamformingResult._from_state(
                read_result_states(path, first=True)[0]
            )
            _check_compatible(first, self, context=f"{path}: ")

        write_result_state(self._state(), path, append=append, owrite=owrite)

    def info(self) -> None:
        """Print a summary of the result."""
        print(self.summary())

    def summary(self) -> str:
        """
        Return a summary of the result as text.

        Returns
        -------
        str
            Configuration, frequencies, blocks and strongest maxima.
        """
        config = self.config

        if config is None:
            return "BeamformingResult: empty"

        frequencies = self.frequencies
        lines = [
            f"BeamformingResult: {config.method.value}, "
            f"{config.components.value}, {len(self.sids)} stations",
            f"  frequencies : {len(frequencies)}"
            + (
                f" ({frequencies[0]:.4g} - {frequencies[-1]:.4g} Hz)"
                if len(frequencies) else ""
            ),
            f"  blocks      : {self.nblocks} (all frequencies)",
            f"  maxima      : {len(self.maxima)}",
            f"  spectra     : {'yes' if self.spectra else 'no'}",
            f"  power units : {self.power_units}",
        ]

        for component in self.components:
            best = self.maxima_at(rank=1, component=component)

            if best:
                velocity = np.median([item.velocity for item in best])
                lines.append(
                    f"  {component}: median rank-1 velocity "
                    f"{velocity:.4g} m/s"
                )

        return "\n".join(lines)

    # -- helpers --------------------------------------------------------

    @classmethod
    def _from_parts(
        cls,
        maxima: list,
        spectra: Optional[list],
        config: BeamformingConfig,
        block_times: dict,
        sids: tuple,
        coordinates: np.ndarray,
        units: Optional[str],
    ) -> BeamformingResult:
        """Result of one analysis (one set of stations)."""
        result = cls()
        result.maxima = list(maxima)
        result.spectra = None if spectra is None else list(spectra)
        result.config = config.copy()
        result.block_times = {
            frequency: tuple(times) for frequency, times in block_times.items()
        }
        result.units = units
        coordinates = np.array(coordinates, dtype=float)
        coordinates.setflags(write=False)
        result._station_sets = [(tuple(sids), coordinates)]
        result._block_stations = {
            frequency: (0,) * len(times)
            for frequency, times in block_times.items()
        }
        return result

    @classmethod
    def _from_state(cls, state: dict) -> BeamformingResult:
        """Result from the state of a file."""
        result = cls()
        result.maxima = list(state["maxima"])
        result.spectra = (
            None if state["spectra"] is None else list(state["spectra"])
        )
        result.config = BeamformingConfig(state["config"])
        result.block_times = dict(state["block_times"])
        result.units = state["units"]
        result._station_sets = [
            (tuple(sids), np.array(coordinates, dtype=float))
            for sids, coordinates in state["station_sets"]
        ]
        result._block_stations = dict(state["block_stations"])
        return result

    def _state(self) -> dict:
        """State written to a file."""
        return {
            "config": self.config.to_dict(),
            "maxima": list(self.maxima),
            "spectra": None if self.spectra is None else list(self.spectra),
            "block_times": dict(self.block_times),
            "block_stations": dict(self._block_stations),
            "station_sets": [
                (sids, np.array(coordinates))
                for sids, coordinates in self._station_sets
            ],
            "units": self.units,
        }

    def _join(self, other: BeamformingResult) -> None:
        """Add the blocks of a compatible result."""
        if other.config is None:
            return

        if self.config is None:
            self.config = other.config.copy()
            self.units = other.units
            self.spectra = None if other.spectra is None else []
            self.block_times = dict.fromkeys(other.block_times, ())
            self._block_stations = dict.fromkeys(other.block_times, ())

        mapping = []

        for sids, coordinates in other._station_sets:
            for index, (names, values) in enumerate(self._station_sets):
                if names == sids and np.array_equal(values, coordinates):
                    mapping.append(index)
                    break
            else:
                mapping.append(len(self._station_sets))
                self._station_sets.append((sids, coordinates))

        offsets = {}

        for frequency, times in other.block_times.items():
            offsets[frequency] = len(self.block_times[frequency])
            self.block_times[frequency] = self.block_times[frequency] + times
            self._block_stations[frequency] = (
                self._block_stations[frequency] + tuple(
                    mapping[index]
                    for index in other._block_stations[frequency]
                )
            )

        self.maxima.extend(
            _replace_block(maximum, offsets[maximum.frequency])
            for maximum in other.maxima
        )

        if other.spectra is not None:
            joined = []

            for spectrum in other.spectra:
                mine = [
                    item for item in self.spectra
                    if item.frequency == spectrum.frequency
                    and item.component == spectrum.component
                ]
                joined.append(
                    _join_spectra(mine + [spectrum]) if mine else spectrum
                )

            self.spectra = joined

    def _component(self, component: Optional[str]) -> str:
        """Return the requested component, or the only one processed."""
        components = self.components

        if component is None:
            if len(components) != 1:
                raise ValueError(
                    "choose the component: " + ", ".join(components)
                )

            return components[0]

        if component not in components:
            raise ValueError(f"component {component!r} was not processed")

        return component

    def _nearest_frequency(self, frequency: float) -> float:
        """Return the analysis frequency nearest to a value."""
        frequencies = self.frequencies

        if not len(frequencies):
            raise ValueError("the result has no frequencies")

        return float(frequencies[np.argmin(np.abs(frequencies - frequency))])


def _check_compatible(
    reference: Optional[BeamformingResult],
    other: BeamformingResult,
    context: str = "",
) -> None:
    """Raise if two results cannot be joined."""
    if reference is None or reference is other:
        return

    if reference.config.to_dict() != other.config.to_dict():
        raise ValueError(f"{context}results with different configurations")

    if set(reference.block_times) != set(other.block_times):
        raise ValueError(f"{context}results with different frequencies")

    if reference.units != other.units:
        raise ValueError(
            f"{context}results in different units ({reference.units}, "
            f"{other.units})"
        )


def _replace_block(maximum: SpectralMaximum, offset: int) -> SpectralMaximum:
    """Return a maximum with the block index shifted."""
    if offset == 0:
        return maximum

    return SpectralMaximum(**{
        **maximum.__dict__, "block": maximum.block + offset
    })


def _join_spectra(items: Sequence) -> FrequencySpectrum:
    """Join the blocks of spectra of the same frequency and component."""
    def stack(name: str) -> Optional[np.ndarray]:
        values = [getattr(item, name) for item in items]
        return None if values[0] is None else np.concatenate(values)

    first = items[0]

    return FrequencySpectrum(
        frequency=first.frequency,
        component=first.component,
        times=tuple(time for item in items for time in item.times),
        azimuths=first.azimuths,
        grid_values=first.grid_values,
        quantity=first.quantity,
        n_realizations=first.n_realizations,
        effective_realizations=first.effective_realizations,
        power=stack("power"),
        power_standard_error=stack("power_standard_error"),
        ellipticity=stack("ellipticity"),
        ellipticity_angle=stack("ellipticity_angle"),
    )


# ----------------------------------------------------------------------------
# Beamformer
# ----------------------------------------------------------------------------


class Beamformer:
    """
    Frequency-wavenumber processing of an array.

    Parameters
    ----------
    array : ArrayData, optional
        Array to process (kept by reference: stations added later are
        processed too).
    config : BeamformingConfig, str or path-like, optional
        Configuration, or a JSON file to read it from. The beamformer
        keeps a copy (:attr:`config`).
    cache : bool, optional
        Keep the steering vectors of the search grid between calls,
        while station coordinates and grid do not change (repeated
        processing, e.g. in real time). Default is False: the memory is
        about ``16 * n_stations * n_grid * n_frequencies`` bytes.

    Examples
    --------
    >>> beamformer = Beamformer(array, "capon.json")
    >>> result = beamformer.process()
    >>> result = beamformer.process("2017-03-01T10:00", "2017-03-01T10:05")
    >>> beamformer.config.update(method="music")
    """

    def __init__(
        self,
        array: Optional[ArrayData] = None,
        config: Any = None,
        cache: bool = False,
    ) -> None:
        self._array: Optional[ArrayData] = None
        self._config: Optional[BeamformingConfig] = None
        self._cache: dict = {}
        self._cache_key: Any = None
        self.cache = bool(cache)

        if array is not None:
            self.array = array

        if config is not None:
            self.config = config

    def __repr__(self) -> str:
        config = self._config
        method = "no config" if config is None else (
            f"{config.method.value}, {config.components.value}"
        )
        stations = 0 if self._array is None else self._array.nsta
        return (
            f"Beamformer({method}, {stations} stations, cache={self.cache})"
        )

    # -- properties -----------------------------------------------------

    @property
    def array(self) -> Optional[ArrayData]:
        """Array processed."""
        return self._array

    @array.setter
    def array(self, array: ArrayData) -> None:
        if not isinstance(array, ArrayData):
            raise TypeError("array must be an ArrayData")

        self._array = array

    @property
    def config(self) -> Optional[BeamformingConfig]:
        """Configuration (the copy of the beamformer, changeable)."""
        return self._config

    @config.setter
    def config(self, config: Any) -> None:
        if isinstance(config, BeamformingConfig):
            config = config.copy()
        elif isinstance(config, (str, os.PathLike)):
            config = BeamformingConfig(config)
        else:
            raise TypeError(
                "config must be a BeamformingConfig or a JSON file"
            )

        self._config = config

    # -- processing -----------------------------------------------------

    def process(
        self,
        starttime: Any = None,
        endtime: Any = None,
        gaps: str = "drop",
    ) -> BeamformingResult:
        """
        Process the array on a time window.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Window, bounds included. Default: the window common to all
            the stations.
        gaps : {"drop", "raise"}, optional
            Stations that do not cover the window: excluded with an
            :class:`~.data.ArrayDataWarning` (default) or an error.

        Returns
        -------
        BeamformingResult
            Maxima and, if configured, spectra.

        Raises
        ------
        ValueError
            If array or configuration are missing or incomplete, the
            window is invalid, the components are not observable at
            every station, the coordinate covariances are missing for
            the geometric uncertainty, or a block has too few windows
            for the estimator.

        Warns
        -----
        BeamformingWarning
            For frequencies whose windows do not fill a covariance
            block (no result at these frequencies).
        """
        if self._array is None:
            raise ValueError("the beamformer has no array")

        if self._config is None:
            raise ValueError("the beamformer has no configuration")

        self._config.check_complete()
        array = self._array.window(starttime, endtime, gaps=gaps)
        steering = None

        if self.cache:
            key = (
                array.coordinates.tobytes(),
                json.dumps(self._config.grid.to_dict()),
            )

            if key != self._cache_key:
                self._cache_key = key
                self._cache = {}

            steering = self._cache

        return _beamform(array, self._config, steering_cache=steering)


def _beamform(
    array: ArrayData,
    config: BeamformingConfig,
    steering_cache: Optional[dict] = None,
) -> BeamformingResult:
    """
    Apply frequency-wavenumber processing to a synchronized array.

    ``steering_cache`` keeps the steering vectors by frequency; it must
    be used with one geometry and grid only (see :class:`Beamformer`).
    """
    inputs = _MODE_INPUTS[config.components]
    motion = np.moveaxis(
        array.matrix([_DIRECTIONS[name] for name in inputs]), 1, 0
    )
    covariances = None

    if config.geometric_uncertainty is not GeometricUncertainty.NONE:
        covariances = array.coordinate_covariances

        if covariances is None:
            raise ValueError(
                "the coordinate uncertainty needs coordinate covariances"
            )

    if config.method is not BeamformingMethod.BARTLETT:
        dimension = array.nsta * (
            2 if config.components is ComponentMode.RAYLEIGH_JOINT else 1
        )
        check_realizations(
            config.window.windows_per_block, dimension, config.method,
            config.strict,
        )

    context = _Context(array, config, inputs, covariances, steering_cache)
    skipped = []

    for frequency, values, starts, weights in _spectral_slices(
        motion, array.delta, config
    ):
        if not context.add(frequency, values, starts, weights):
            skipped.append(frequency)

    if skipped:
        warn_external(
            f"no complete covariance block at {len(skipped)} frequencies "
            f"({min(skipped):.4g} - {max(skipped):.4g} Hz): the data are "
            "too short for the windows",
            BeamformingWarning,
        )

    return context.result()


# ----------------------------------------------------------------------------
# Spectral slices
# ----------------------------------------------------------------------------


def _spectral_slices(motion: np.ndarray, delta: float, config: Any):
    """
    Yield the Fourier coefficients of each analysis frequency.

    Yields ``(frequency, values, starts, weights)`` with ``values`` of
    shape ``(n_windows, n_inputs, n_stations)`` (None if no window
    fits), the window starts and the taper.
    """
    window = config.window
    spec = ("tukey", window.taper)
    nsamp = motion.shape[-1]
    frequencies = config.frequencies

    if config.transform is TransformMethod.DIRECT:
        for frequency in frequencies.targets():
            size = window.size(delta, frequency)
            starts = window_starts(nsamp, size, overlap=window.overlap)

            if starts.size == 0:
                yield float(frequency), None, starts, None
                continue

            _, values = coefficients(
                motion, delta, frequency, size, overlap=window.overlap,
                window=spec,
            )
            yield float(frequency), values[0], starts, taper(size, spec)

        return

    size = window.size(delta)
    padding = config.fft_padding

    if frequencies.is_sampled:
        targets = frequencies.targets()
    else:
        padding = padding or 1
        grid = fourier.frequency_axis(delta, size * padding)
        targets = grid[
            (grid >= frequencies.minimum) & (grid <= frequencies.maximum)
            & (grid > 0.0) & (grid < 0.5 / delta)
        ]

        if targets.size == 0:
            raise ValueError("no FFT bin falls in the frequency range")

    starts = window_starts(nsamp, size, overlap=window.overlap)

    if starts.size == 0:
        for frequency in targets:
            yield float(frequency), None, starts, None

        return

    selection = (
        "interp" if config.fft_selection is FFTSelection.INTERP
        else "nearest"
    )
    out, values = coefficients(
        motion, delta, targets, size, overlap=window.overlap, window=spec,
        transform="rfft", selection=selection, padding=padding,
    )
    weights = taper(size, spec)

    for index, frequency in enumerate(out):
        yield float(frequency), values[index], starts, weights


# ----------------------------------------------------------------------------
# Processing of one frequency
# ----------------------------------------------------------------------------


class _Context:
    """Accumulator of the results of one call of :func:`_beamform`."""

    def __init__(
        self,
        array: ArrayData,
        config: BeamformingConfig,
        inputs: tuple,
        covariances: Optional[np.ndarray],
        cache: Optional[dict],
    ) -> None:
        self.array = array
        self.config = config
        self.inputs = inputs
        self.covariances = covariances
        self.cache = cache
        self.maxima: list = []
        self.spectra: list = []
        self.block_times: dict = {}

    def steering(self, frequency: float) -> np.ndarray:
        """Steering vectors of the grid at a frequency."""
        if self.cache is not None and frequency in self.cache:
            return self.cache[frequency]

        kx, ky = self.config.grid.wavevectors(frequency)
        value = steering_vectors(self.array.coordinates, kx, ky)

        if self.cache is not None:
            self.cache[frequency] = value

        return value

    def add(
        self,
        frequency: float,
        values: Optional[np.ndarray],
        starts: np.ndarray,
        weights: Optional[np.ndarray],
    ) -> bool:
        """Process one frequency; return False if it has no block."""
        config = self.config
        window = config.window
        count = window.windows_per_block
        blocks = window_starts(
            starts.size, count, overlap=window.block_overlap
        )
        outputs = _MODE_OUTPUTS[config.components]
        grid = config.grid
        shape = (blocks.size, grid.azimuth_size, grid.size)

        if blocks.size == 0:
            self.block_times[frequency] = ()
            self._store_spectra(frequency, outputs, shape, {}, (), 1.0)
            return False

        step = int(starts[1] - starts[0]) if starts.size > 1 else len(weights)
        effective = effective_realizations(weights, step, count)
        steering = self.steering(frequency)
        attenuation = (
            None if self.covariances is None
            else _attenuation(self.covariances, frequency, grid)
        )
        coefficient = dict(zip(self.inputs, np.moveaxis(values, 1, 0)))
        buffers = {name: {} for name in outputs}
        times = []
        delta = self.array.delta
        size = len(weights)

        for index, first in enumerate(blocks):
            block = slice(int(first), int(first) + count)
            last = starts[first + count - 1] + size - 1
            centre = 0.5 * (starts[first] + last)
            time = self.array.starttime + float(centre) * delta
            times.append(time)

            for name in outputs:
                fields = self._component(
                    name, coefficient, block, steering, attenuation, effective
                )

                for key, value in fields.items():
                    buffers[name].setdefault(key, np.empty(shape))[index] = (
                        value
                    )

                self._store_maxima(
                    index, time, frequency, name, fields, steering, grid
                )

        self.block_times[frequency] = tuple(times)
        self._store_spectra(
            frequency, outputs, shape, buffers, tuple(times), effective
        )
        return True

    def _component(
        self,
        name: str,
        values: dict,
        block: slice,
        steering: np.ndarray,
        attenuation: Optional[np.ndarray],
        effective: float,
    ) -> dict:
        """Power (and other fields) of one component in one block."""
        config = self.config
        count = config.window.windows_per_block
        statistics = (
            config.statistical_uncertainty
            is StatisticalUncertainty.EMPIRICAL
        )

        if name == "rayleigh_joint":
            return _rayleigh_joint(values, block, steering, config)

        if name == "vertical":
            series = values["vertical"][block]
            power = beam_power(
                cross_spectral_matrix(series), steering, count,
                config.method, config.strict, config.diagonal_loading,
                config.music_sources, attenuation,
            )
            fields = {"power": power}

            if statistics:
                fields["power_standard_error"] = _standard_error(
                    series, steering, attenuation, effective, config
                )

            return fields

        east = values["east"][block]
        north = values["north"][block]
        ee = cross_spectral_matrix(east)
        nn = cross_spectral_matrix(north)
        en = cross_spectral_matrix(east, north)
        azimuths = np.radians(config.grid.azimuth_values)
        power = np.empty(steering.shape[:-1])
        error = np.empty(steering.shape[:-1]) if statistics else None

        for row, angle in enumerate(azimuths):
            if name == "radial":
                we, wn = math.sin(angle), math.cos(angle)
            else:
                we, wn = math.cos(angle), -math.sin(angle)

            matrix = we * we * ee + wn * wn * nn + we * wn * (en + en.conj().T)
            local = None if attenuation is None else attenuation[row:row + 1]
            power[row] = beam_power(
                matrix, steering[row:row + 1], count, config.method,
                config.strict, config.diagonal_loading,
                config.music_sources, local,
            )[0]

            if statistics:
                error[row] = _standard_error(
                    we * east + wn * north, steering[row:row + 1], local,
                    effective, config,
                )[0]

        fields = {"power": power}

        if statistics:
            fields["power_standard_error"] = error

        return fields

    def _store_maxima(
        self,
        block: int,
        time: Date,
        frequency: float,
        component: str,
        fields: dict,
        steering: np.ndarray,
        grid: Any,
    ) -> None:
        """Append the local maxima of one spectrum."""
        config = self.config
        power = fields["power"]
        wavenumbers = grid.wavenumbers(frequency)
        azimuths = grid.azimuth_values

        for rank, (row, column) in enumerate(
            _local_maxima(
                power, config.maxima, config.peak_neighborhood,
                config.minimum_relative_power,
            ),
            start=1,
        ):
            extra = {
                key: float(fields[key][row, column])
                for key in ("power_standard_error", "ellipticity",
                            "ellipticity_angle")
                if key in fields
            }
            self.maxima.append(SpectralMaximum(
                block=block,
                time=time,
                frequency=frequency,
                component=component,
                rank=rank,
                wavenumber=float(wavenumbers[column]),
                azimuth=float(azimuths[row]),
                power=float(power[row, column]),
                **extra,
            ))

    def _store_spectra(
        self,
        frequency: float,
        outputs: tuple,
        shape: tuple,
        buffers: dict,
        times: tuple,
        effective: float,
    ) -> None:
        """Append the spectra of one frequency."""
        if not self.config.spectra:
            return

        grid = self.config.grid

        for name in outputs:
            fields = buffers.get(name, {})
            self.spectra.append(FrequencySpectrum(
                frequency=frequency,
                component=name,
                times=times,
                azimuths=grid.azimuth_values,
                grid_values=grid.values,
                quantity=grid.quantity,
                n_realizations=self.config.window.windows_per_block,
                effective_realizations=effective,
                power=fields.get("power", np.empty(shape)),
                power_standard_error=fields.get("power_standard_error"),
                ellipticity=fields.get("ellipticity"),
                ellipticity_angle=fields.get("ellipticity_angle"),
            ))

    def result(self) -> BeamformingResult:
        """Return the accumulated result."""
        return BeamformingResult._from_parts(
            maxima=self.maxima,
            spectra=(
                self.spectra if self.config.spectra else None
            ),
            config=self.config,
            block_times=self.block_times,
            sids=self.array.sids,
            coordinates=self.array.coordinates,
            units=self.array.units,
        )


def _standard_error(
    series: np.ndarray,
    steering: np.ndarray,
    attenuation: Optional[np.ndarray],
    effective: float,
    config: BeamformingConfig,
) -> np.ndarray:
    """Empirical standard error of Bartlett or Capon power."""
    if config.method is BeamformingMethod.BARTLETT:
        return bartlett_power_statistics(
            series, steering, effective, attenuation
        )[1]

    return capon_power_statistics(
        series, steering, effective, config.diagonal_loading
    )[1]


def _rayleigh_joint(
    values: dict,
    block: slice,
    steering: np.ndarray,
    config: BeamformingConfig,
) -> dict:
    """
    Rayleigh three-component beamforming (RTBF) of one block.

    For each azimuth the radial and vertical coefficients form a matrix
    of dimension 2n; the Capon power ``P`` of the steering vectors
    ``[-1j * e * a, a]`` is evaluated for every signed ellipticity
    ``e`` and the objective ``e**2 * P**2`` is maximized over ``e``
    (Wathelet et al., 2018).
    """
    grid = config.ellipticity_grid
    angles = grid.angles
    ellipticity = grid.values
    keep = np.abs(ellipticity) > np.finfo(float).eps
    angles = angles[keep]
    ellipticity = ellipticity[keep]
    nazimuth, nvalue, nsta = steering.shape
    power = np.empty((nazimuth, nvalue))
    best = np.empty((nazimuth, nvalue))
    best_angle = np.empty((nazimuth, nvalue))
    east = values["east"][block]
    north = values["north"][block]
    vertical = values["vertical"][block]
    columns = np.arange(nvalue)

    for row, angle in enumerate(np.radians(config.grid.azimuth_values)):
        radial = math.sin(angle) * east + math.cos(angle) * north
        matrix = regularized_matrix(
            cross_spectral_matrix(np.concatenate((radial, vertical), axis=1)),
            config.diagonal_loading,
        )
        a = steering[row]
        joint = np.concatenate((
            -1j * ellipticity[:, None, None] * a[None],
            np.broadcast_to(a[None], (ellipticity.size,) + a.shape),
        ), axis=2).reshape(-1, 2 * nsta)
        solved = np.linalg.solve(matrix, joint.T)
        denominator = np.real(np.einsum("gi,ig->g", joint.conj(), solved))
        capon = np.full(denominator.shape, np.nan)
        positive = denominator > 0.0
        capon[positive] = 1.0 / denominator[positive]
        objective = ellipticity[:, None] ** 2 * capon.reshape(-1, nvalue) ** 2
        winner = np.argmax(np.where(np.isfinite(objective), objective,
                                    -np.inf), axis=0)
        power[row] = objective[winner, columns]
        best[row] = ellipticity[winner]
        best_angle[row] = angles[winner]

    return {
        "power": power,
        "ellipticity": best,
        "ellipticity_angle": best_angle,
    }


def _attenuation(
    covariances: np.ndarray,
    frequency: float,
    grid: Any,
) -> np.ndarray:
    """
    Expected coherence of each station at each wavevector.

    For Gaussian coordinate errors of covariance ``S`` (horizontal
    part) the expected phase factor is ``exp(-k^T S k / 2)``; shape
    ``(n_azimuths, n_values, n_stations)``.
    """
    kx, ky = grid.wavevectors(frequency)
    sxx = covariances[:, 0, 0]
    sxy = covariances[:, 0, 1]
    syy = covariances[:, 1, 1]
    variance = (
        kx[..., None] ** 2 * sxx + 2.0 * kx[..., None] * ky[..., None] * sxy
        + ky[..., None] ** 2 * syy
    )
    return np.exp(-0.5 * np.maximum(variance, 0.0))


def _local_maxima(
    power: np.ndarray,
    count: int,
    neighborhood: tuple,
    minimum_relative_power: float,
) -> list[tuple[int, int]]:
    """
    Indices of the strongest local maxima of a polar spectrum.

    The azimuth axis is periodic; NaN values are never maxima.
    """
    finite = np.isfinite(power)

    if not np.any(finite):
        return []

    work = np.where(finite, power, -np.inf)
    peaks = finite & (
        work == maximum_filter(work, size=neighborhood,
                               mode=("wrap", "nearest"))
    )

    if minimum_relative_power > 0.0:
        peaks &= work >= minimum_relative_power * np.max(work[finite])

    index = np.argwhere(peaks)
    order = np.argsort(work[index[:, 0], index[:, 1]], kind="stable")[::-1]

    return [(int(row), int(col)) for row, col in index[order[:count]]]
