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
Frequency domain decomposition (FDD) of array data.

At each frequency the cross-spectral matrix of the channels of the
array (each component of each station) is estimated on blocks of
consecutive windows (:class:`shakelab.signals.crossspectral.
CrossSpectra`) and decomposed into singular values and vectors. Close
to a resonance (a mode of a structure, a two-dimensional resonance of a
valley) one source of motion dominates: the first singular value peaks
and the first singular vector is proportional to the shape of the mode.
The decomposition of every block gives the statistics of the singular
values (geometric mean and logarithmic standard deviation) and of the
shapes.

- :class:`FDD`: processor of an :class:`~.data.ArrayData` with an
  :class:`~.config.FDDConfig`;
- :class:`FDDResult`: singular values and vectors of every block, with
  :meth:`FDDResult.singular_values`, :meth:`FDDResult.peaks` and
  :meth:`FDDResult.mode`;
- :class:`Mode`: shape of a mode at a frequency;
- :func:`mac`: modal assurance criterion between two shapes.

Phase and sign
--------------
A singular vector is defined up to a complex unit factor. Each vector
is rotated to be as real as possible (angle ``-arg(u^T u) / 2``, which
maximizes the norm of its real part) and the vectors of successive
blocks are given the same sign (positive scalar product with the sum
of the previous ones); the shape of a mode is the mean of the real
parts, with its largest component positive.

References
----------
Brincker, R., Zhang, L., and Andersen, P. (2001). Modal identification
of output-only systems using frequency domain decomposition. Smart
Materials and Structures, 10(3), 441-445.

Allemang, R. J., and Brown, D. L. (1982). A correlation coefficient for
modal vector analysis. Proceedings of the 1st International Modal
Analysis Conference, 110-116.

Ermert, L., Poggi, V., Burjanek, J., and Faeh, D. (2014). Fundamental
and higher two-dimensional resonance modes of an Alpine valley.
Geophysical Journal International, 198(2), 795-811.
"""

from __future__ import annotations

import math
import os
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np

from shakelab.exploration.yarray.config import (
    FDDComponents,
    FDDConfig,
    FDDNormalization,
    WindowType,
)
from shakelab.exploration.yarray.data import ArrayData
from shakelab.libutils.utils import warn_external
from shakelab.signals import fourier
from shakelab.signals.crossspectral import CrossSpectra

__all__ = [
    "FDD",
    "FDDResult",
    "FDDWarning",
    "Mode",
    "mac",
]

# Kind of the result states in the result files (see :mod:`.io`).
_KIND = "fdd"


# ----------------------------------------------------------------------------
# Warnings and helpers
# ----------------------------------------------------------------------------


class FDDWarning(UserWarning):
    """Frequencies without complete blocks."""


def mac(first: Any, second: Any) -> float:
    """
    Return the modal assurance criterion of two shapes.

    ``|a^H b|**2 / ((a^H a) (b^H b))``: 1 for proportional shapes
    (also with a complex factor), 0 for orthogonal ones.

    Parameters
    ----------
    first, second : array-like
        Shapes (real or complex) with the same number of channels.

    Returns
    -------
    float
        Value in ``[0, 1]``.

    Raises
    ------
    ValueError
        If the shapes are not one-dimensional, differ in length or one
        is zero.

    Examples
    --------
    >>> mac([1.0, 2.0], [2.0, 4.0])
    1.0
    """
    a = np.asarray(first, dtype=complex)
    b = np.asarray(second, dtype=complex)

    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError("the shapes must be vectors of the same length")

    norm = np.vdot(a, a).real * np.vdot(b, b).real

    if norm == 0.0:
        raise ValueError("the shapes must not be zero")

    return float(min(1.0, abs(np.vdot(a, b)) ** 2 / norm))


def _align(vectors: np.ndarray) -> np.ndarray:
    """
    Return vectors ``(n_blocks, n)`` with normalized phase and sign.

    Each vector is rotated to be as real as possible; the sign of each
    makes its real part agree with the sum of the previous ones, and
    the largest component of the sum is positive.
    """
    if vectors.shape[0] == 0:
        return vectors.copy()

    angle = -0.5 * np.angle(np.sum(vectors * vectors, axis=1))
    out = vectors * np.exp(1j * angle)[:, None]
    reference = out[0].real.copy()

    for index in range(1, out.shape[0]):
        if np.dot(reference, out[index].real) < 0.0:
            out[index] = -out[index]

        reference += out[index].real

    if reference[np.argmax(np.abs(reference))] < 0.0:
        out = -out

    return out


def _log(values: np.ndarray) -> np.ndarray:
    """Natural logarithm of non-negative values (zero to the tiniest)."""
    return np.log(np.maximum(values, np.finfo(float).tiny))


# ----------------------------------------------------------------------------
# Mode
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class Mode:
    """
    Shape of a mode at a frequency.

    Parameters
    ----------
    frequency : float
        Frequency in Hz.
    shape : numpy.ndarray
        Real shape ``(n_channels,)``, unit norm, largest component
        positive (mean of the aligned singular vectors of the blocks).
    shape_std : numpy.ndarray
        Standard deviation of each component among the blocks, on the
        scale of ``shape`` (NaN with one block).
    value : float
        Geometric mean of the singular value among the blocks.
    value_std : float
        Standard deviation of its natural logarithm (NaN with one
        block).
    mac : float
        Mean modal assurance criterion between ``shape`` and the
        singular vectors of the blocks (1 for a stable shape).
    blocks : int
        Number of blocks.
    index : int
        Index of the singular value (0 for the first).
    stations : tuple of str
        Station of each channel.
    components : tuple of str
        Component of each channel.
    coordinates : numpy.ndarray
        Position ``(n_channels, 3)`` of each channel in metres (East,
        North, Up).
    damping : float or None
        Damping ratio (None: not estimated).
    method : str
        Method of the estimate (``"fdd"``).
    """

    frequency: float
    shape: np.ndarray
    shape_std: np.ndarray
    value: float
    value_std: float
    mac: float
    blocks: int
    index: int
    stations: tuple
    components: tuple
    coordinates: np.ndarray
    damping: Optional[float] = None
    method: str = "fdd"

    @property
    def channels(self) -> tuple[str, ...]:
        """Labels of the channels (``"station component"``)."""
        return tuple(
            f"{sid} {name}" for sid, name in zip(self.stations,
                                                 self.components)
        )

    def summary(self) -> str:
        """
        Return the mode as text.

        Returns
        -------
        str
            Frequency, singular value, stability and shape.
        """
        lines = [
            f"Mode at {self.frequency:.4g} Hz (singular value "
            f"{self.index + 1}, {self.method}): value {self.value:.4g}, "
            f"MAC {self.mac:.3f}, {self.blocks} blocks"
            + (
                "" if self.damping is None
                else f", damping {100.0 * self.damping:.2f}%"
            ),
            "  channel                          shape      std",
        ]

        for label, value, spread in zip(self.channels, self.shape,
                                        self.shape_std):
            lines.append(f"  {label:<30s} {value:8.4f} {spread:8.4f}")

        return "\n".join(lines)

    def info(self) -> None:
        """Print the summary (see :meth:`summary`)."""
        print(self.summary())


# ----------------------------------------------------------------------------
# Result
# ----------------------------------------------------------------------------


class FDDResult:
    """
    Result of a frequency domain decomposition.

    A result holds the singular values and the first singular vectors
    of every block of one or more analyses with the same configuration
    and channels: :meth:`append` adds those of a later analysis,
    :meth:`write` with ``append=True`` adds them to a file, which is
    read back joined.

    Parameters
    ----------
    source : str or path-like, optional
        Result file to read (see :meth:`read`). Default: an empty
        result.

    Attributes
    ----------
    config : FDDConfig or None
        Configuration used (None for an empty result).
    block_times : dict
        Centres of the blocks (Date), by frequency.
    values : dict
        Singular values ``(n_blocks, n_channels)`` in decreasing order,
        by frequency.
    vectors : dict
        First singular vectors ``(n_blocks, n_channels, n_vectors)``
        with normalized phase and sign, by frequency.
    stations : tuple of str
        Station of each channel.
    components : tuple of str
        Component of each channel.
    coordinates : numpy.ndarray
        Position ``(n_channels, 3)`` of each channel in metres.
    units : str or None
        Units of the data (None if unknown).

    Examples
    --------
    >>> result = FDD(array, config).process()
    >>> result.peaks(fmin=0.2, fmax=0.8)
    >>> result.mode(0.32).info()
    >>> result.write("fdd.pkl")
    """

    def __init__(self, source: Any = None) -> None:
        self.config: Optional[FDDConfig] = None
        self.block_times: dict = {}
        self.values: dict = {}
        self.vectors: dict = {}
        self.stations: tuple = ()
        self.components: tuple = ()
        self.coordinates: np.ndarray = np.empty((0, 3))
        self.units: Optional[str] = None

        if source is not None:
            self.read(source)

    def __repr__(self) -> str:
        if self.config is None:
            return "FDDResult(empty)"

        return (
            f"FDDResult({self.config.components.value}, "
            f"{self.nchannels} channels, {len(self.block_times)} "
            f"frequencies, {self.nblocks} blocks)"
        )

    # -- properties -----------------------------------------------------

    @property
    def frequencies(self) -> np.ndarray:
        """Analysis frequencies, in Hz."""
        return np.array(sorted(self.block_times))

    @property
    def nblocks(self) -> int:
        """Number of blocks, summed over the frequencies."""
        return sum(len(times) for times in self.block_times.values())

    @property
    def nchannels(self) -> int:
        """Number of channels."""
        return len(self.stations)

    @property
    def channels(self) -> tuple[str, ...]:
        """Labels of the channels (``"station component"``)."""
        return tuple(
            f"{sid} {name}" for sid, name in zip(self.stations,
                                                 self.components)
        )

    @property
    def power_units(self) -> Optional[str]:
        """Units of the singular values."""
        if self.config is None:
            return None

        return f"({self.units or 'unit'})**2/Hz"

    # -- statistics -----------------------------------------------------

    def singular_values(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return the statistics of the singular values.

        Returns
        -------
        frequency : numpy.ndarray
            Frequencies in Hz, ``(n_freq,)``.
        mean : numpy.ndarray
            Geometric mean among the blocks, ``(n_freq, n_channels)``
            (NaN at frequencies without blocks).
        log_std : numpy.ndarray
            Standard deviation of the natural logarithm among the
            blocks (NaN with fewer than two blocks): the band
            ``mean * exp(+/- log_std)``.
        """
        frequency = self.frequencies
        shape = (frequency.size, self.nchannels)
        mean = np.full(shape, np.nan)
        spread = np.full(shape, np.nan)

        for row, value in enumerate(frequency):
            logs = _log(self.values[value])

            if logs.shape[0]:
                mean[row] = np.exp(logs.mean(axis=0))

            if logs.shape[0] > 1:
                spread[row] = logs.std(axis=0, ddof=1)

        return frequency, mean, spread

    def peaks(
        self,
        fmin: Optional[float] = None,
        fmax: Optional[float] = None,
        index: int = 0,
        prominence: float = 2.0,
        count: Optional[int] = None,
    ) -> np.ndarray:
        """
        Return the frequencies of the peaks of a singular value.

        The peaks are the local maxima of the geometric mean of the
        singular value over the analysis frequencies whose prominence
        is at least ``prominence``: the peak is that many times larger
        than the higher of the two minima that separate it from larger
        peaks (or from the ends of the range). The prominence discards
        the small fluctuations of the estimate.

        Parameters
        ----------
        fmin, fmax : float, optional
            Frequency range in Hz. Default: all the frequencies.
        index : int, optional
            Index of the singular value (0, default, for the first).
        prominence : float, optional
            Smallest ratio between a peak and its separating minimum
            (at least one). Default is 2.
        count : int, optional
            Keep only the ``count`` most prominent peaks. Default: all.

        Returns
        -------
        numpy.ndarray
            Frequencies in Hz, increasing.

        Raises
        ------
        ValueError
            If the result is empty, ``index`` is out of range,
            ``prominence`` is smaller than one or ``count`` smaller
            than one.
        """
        from scipy.signal import find_peaks

        if self.config is None:
            raise ValueError("the result is empty")

        if not 0 <= index < self.nchannels:
            raise ValueError(
                f"index must be in [0, {self.nchannels - 1}]"
            )

        if not prominence >= 1.0:
            raise ValueError("prominence must be at least 1")

        if count is not None and count < 1:
            raise ValueError("count must be at least 1")

        frequency, mean, _ = self.singular_values()
        inside = np.isfinite(mean[:, index])

        if fmin is not None:
            inside &= frequency >= fmin

        if fmax is not None:
            inside &= frequency <= fmax

        frequency = frequency[inside]
        rows, properties = find_peaks(
            np.log(mean[inside, index]), prominence=math.log(prominence)
        )

        if count is not None:
            order = np.argsort(properties["prominences"])[::-1][:count]
            rows = rows[order]

        return np.sort(frequency[rows])

    def mode(self, frequency: float, index: int = 0) -> Mode:
        """
        Return the shape at the analysis frequency nearest to a value.

        Parameters
        ----------
        frequency : float
            Frequency in Hz.
        index : int, optional
            Index of the singular vector (0, default, for the first;
            fewer than the ``vectors`` of the configuration).

        Returns
        -------
        Mode
            Shape, its spread and stability, and the singular value.

        Raises
        ------
        ValueError
            If the result is empty, the vector was not kept or the
            frequency has no blocks.
        """
        if self.config is None:
            raise ValueError("the result is empty")

        frequencies = self.frequencies
        nearest = float(
            frequencies[np.argmin(np.abs(frequencies - frequency))]
        )
        vectors = self.vectors[nearest]

        if not 0 <= index < vectors.shape[2]:
            raise ValueError(
                f"vector {index} was not kept ({vectors.shape[2]} kept)"
            )

        if vectors.shape[0] == 0:
            raise ValueError(f"no blocks at {nearest:g} Hz")

        aligned = _align(vectors[:, :, index]).real
        shape = aligned.mean(axis=0)
        norm = float(np.linalg.norm(shape))
        blocks = aligned.shape[0]
        spread = (
            aligned.std(axis=0, ddof=1) if blocks > 1
            else np.full(shape.size, np.nan)
        )
        logs = _log(self.values[nearest][:, index])

        return Mode(
            frequency=nearest,
            shape=_readonly(shape / norm),
            shape_std=_readonly(spread / norm),
            value=float(np.exp(logs.mean())),
            value_std=(float(logs.std(ddof=1)) if blocks > 1
                       else math.nan),
            mac=float(np.mean([
                mac(shape, vectors[row, :, index]) for row in range(blocks)
            ])),
            blocks=blocks,
            index=index,
            stations=self.stations,
            components=self.components,
            coordinates=_readonly(self.coordinates),
        )

    # -- joining --------------------------------------------------------

    def append(self, other: Any) -> None:
        """
        Add the blocks of other results.

        Parameters
        ----------
        other : FDDResult or sequence of FDDResult
            Results with the same configuration and channels.

        Raises
        ------
        TypeError
            If an item is not an FDDResult.
        ValueError
            If a result is not compatible.
        """
        items = [other] if isinstance(other, FDDResult) else list(other)

        for item in items:
            if not isinstance(item, FDDResult):
                raise TypeError("results must be FDDResult")

        for item in items:
            if item.config is None:
                continue

            if self.config is None:
                self._adopt(item)
                continue

            _check_compatible(self, item)

            for frequency in self.block_times:
                self.block_times[frequency] = (
                    self.block_times[frequency] + item.block_times[frequency]
                )
                self.values[frequency] = np.concatenate(
                    [self.values[frequency], item.values[frequency]]
                )
                self.vectors[frequency] = np.concatenate(
                    [self.vectors[frequency], item.vectors[frequency]]
                )

    def copy(self) -> FDDResult:
        """
        Return an independent copy.

        Returns
        -------
        FDDResult
            Copy of the result.
        """
        result = FDDResult()
        result.append(self)
        return result

    # -- files ----------------------------------------------------------

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
            If the file does not hold FDD results or they are not
            compatible.
        """
        from shakelab.exploration.yarray.io import read_result_states

        self.append([
            FDDResult._from_state(state)
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
            with the same configuration and channels. Default False.
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
            first = FDDResult._from_state(
                read_result_states(path, first=True)[0]
            )
            _check_compatible(first, self, context=f"{path}: ")

        write_result_state(self._state(), path, append=append, owrite=owrite)

    # -- summary --------------------------------------------------------

    def summary(self) -> str:
        """
        Return a summary of the result as text.

        Returns
        -------
        str
            Configuration, channels, blocks and strongest peaks.
        """
        if self.config is None:
            return "FDDResult: empty"

        frequency = self.frequencies
        lines = [
            repr(self),
            f"  frequencies: {frequency.min():g} - {frequency.max():g} Hz",
            f"  normalization: {self.config.normalization.value}, "
            f"power units: {self.power_units}",
            textwrap.fill(
                "channels: " + ", ".join(self.channels), width=77,
                initial_indent="  ", subsequent_indent="    ",
            ),
        ]
        peaks = self.peaks(count=5)

        if peaks.size:
            lines.append(
                "  peaks of the first singular value: "
                + ", ".join(f"{value:.4g}" for value in peaks) + " Hz"
            )

        return "\n".join(lines)

    def info(self) -> None:
        """Print a summary of the result."""
        print(self.summary())

    # -- helpers --------------------------------------------------------

    @classmethod
    def _from_parts(
        cls,
        config: FDDConfig,
        block_times: dict,
        values: dict,
        vectors: dict,
        stations: Sequence[str],
        components: Sequence[str],
        coordinates: Any,
        units: Optional[str],
    ) -> FDDResult:
        """Result of one analysis."""
        result = cls()
        result.config = config.copy()
        result.block_times = {
            frequency: tuple(times) for frequency, times in block_times.items()
        }
        result.values = dict(values)
        result.vectors = dict(vectors)
        result.stations = tuple(stations)
        result.components = tuple(components)
        result.coordinates = _readonly(coordinates)
        result.units = units
        return result

    @classmethod
    def _from_state(cls, state: dict) -> FDDResult:
        """Result from the state of a file."""
        if state.get("kind") != _KIND:
            raise ValueError("not an FDD result")

        return cls._from_parts(
            FDDConfig(state["config"]), state["block_times"],
            state["values"], state["vectors"], state["stations"],
            state["components"], state["coordinates"], state["units"],
        )

    def _state(self) -> dict:
        """State written to a file."""
        return {
            "kind": _KIND,
            "config": self.config.to_dict(),
            "block_times": dict(self.block_times),
            "values": dict(self.values),
            "vectors": dict(self.vectors),
            "stations": self.stations,
            "components": self.components,
            "coordinates": np.array(self.coordinates),
            "units": self.units,
        }

    def _adopt(self, other: FDDResult) -> None:
        """Take the content of a result (as a copy)."""
        self.config = other.config.copy()
        self.block_times = dict(other.block_times)
        self.values = {
            key: value.copy() for key, value in other.values.items()
        }
        self.vectors = {
            key: value.copy() for key, value in other.vectors.items()
        }
        self.stations = other.stations
        self.components = other.components
        self.coordinates = other.coordinates
        self.units = other.units


def _check_compatible(
    first: FDDResult,
    other: FDDResult,
    context: str = "",
) -> None:
    """Raise if two results cannot be joined."""
    if first.config != other.config:
        raise ValueError(context + "the results have different configurations")

    if (first.stations, first.components) != (other.stations,
                                               other.components):
        raise ValueError(context + "the results have different channels")

    if not np.allclose(first.coordinates, other.coordinates):
        raise ValueError(
            context + "the results have different station coordinates"
        )

    if first.units != other.units:
        raise ValueError(context + "the results have different units")


def _readonly(value: Any) -> np.ndarray:
    """Return a read-only float copy."""
    array = np.array(value, dtype=float)
    array.setflags(write=False)
    return array


# ----------------------------------------------------------------------------
# Processor
# ----------------------------------------------------------------------------

# Directions (east, north, up) of the fixed components.
_DIRECTIONS = {
    "east": (1.0, 0.0, 0.0),
    "north": (0.0, 1.0, 0.0),
    "vertical": (0.0, 0.0, 1.0),
}

_MODE_COMPONENTS = {
    FDDComponents.VERTICAL: ("vertical",),
    FDDComponents.NORTH: ("north",),
    FDDComponents.EAST: ("east",),
    FDDComponents.RADIAL: ("radial",),
    FDDComponents.TRANSVERSE: ("transverse",),
    FDDComponents.HORIZONTAL: ("east", "north"),
    FDDComponents.THREE_COMPONENT: ("east", "north", "vertical"),
}


class FDD:
    """
    Frequency domain decomposition of a seismic array.

    Parameters
    ----------
    array : ArrayData, optional
        Stations to process (kept by reference: later changes of the
        array are seen by the processor).
    config : FDDConfig, str or path-like, optional
        Configuration (copied), or a JSON file.

    Examples
    --------
    >>> fdd = FDD(array, "fdd.json")
    >>> result = fdd.process()
    >>> fdd.config.update(normalization="coherence")
    """

    def __init__(
        self,
        array: Optional[ArrayData] = None,
        config: Any = None,
    ) -> None:
        self._array: Optional[ArrayData] = None
        self._config: Optional[FDDConfig] = None

        if array is not None:
            self.array = array

        if config is not None:
            self.config = config

    def __repr__(self) -> str:
        config = self._config
        method = "no config" if config is None else (
            f"{config.components.value}, {config.normalization.value}"
        )
        stations = 0 if self._array is None else self._array.nsta
        return f"FDD({method}, {stations} stations)"

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
    def config(self) -> Optional[FDDConfig]:
        """Configuration (the copy of the processor, changeable)."""
        return self._config

    @config.setter
    def config(self, config: Any) -> None:
        if isinstance(config, FDDConfig):
            config = config.copy()
        elif isinstance(config, (str, os.PathLike)):
            config = FDDConfig(config)
        else:
            raise TypeError("config must be an FDDConfig or a JSON file")

        self._config = config

    # -- processing -----------------------------------------------------

    def process(
        self,
        starttime: Any = None,
        endtime: Any = None,
        gaps: str = "drop",
    ) -> FDDResult:
        """
        Decompose the cross-spectral matrices on a time window.

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
        FDDResult
            Singular values and vectors of every block.

        Raises
        ------
        ValueError
            If array or configuration are missing or incomplete, the
            window is invalid, the components are not observable at
            every station, or a channel has no power with
            ``normalization="coherence"``.

        Warns
        -----
        FDDWarning
            For frequencies whose windows do not fill a block (no
            result at these frequencies).
        """
        if self._array is None:
            raise ValueError("the processor has no array")

        if self._config is None:
            raise ValueError("the processor has no configuration")

        config = self._config
        config.check_complete()
        array = self._array.window(starttime, endtime, gaps=gaps)
        names = _MODE_COMPONENTS[config.components]
        directions = [_direction(name, config.azimuth) for name in names]
        motion = array.matrix(directions)
        nsta, ncomp, nsamp = motion.shape
        stations = tuple(sid for sid in array.sids for _ in names)
        components = names * nsta
        coordinates = np.repeat(array.coordinates, ncomp, axis=0)
        delta = array.delta
        frequency = _frequencies(config, delta)
        window = config.window
        length = (
            {"duration": window.length}
            if window.type is WindowType.FIXED
            else {"cycles": window.length}
        )
        spectra = CrossSpectra.from_data(
            motion.reshape(nsta * ncomp, nsamp), delta, frequency,
            overlap=window.overlap, window=("tukey", window.taper),
            windows_per_block=window.windows_per_block,
            block_overlap=window.block_overlap,
            subfrequencies=config.subfrequencies,
            bandwidth=config.bandwidth,
            channels=[f"{a} {b}" for a, b in zip(stations, components)],
            **length,
        )

        if config.normalization is FDDNormalization.COHERENCE:
            spectra = spectra.coherence()

        kept = min(config.vectors, nsta * ncomp)
        block_times, values, vectors, skipped = {}, {}, {}, []

        for index, value in enumerate(spectra.frequency):
            value = float(value)
            singular, singular_vectors = spectra.decompose(index)
            values[value] = singular
            vectors[value] = np.stack(
                [_align(singular_vectors[:, :, k]) for k in range(kept)],
                axis=2,
            ) if singular.shape[0] else np.zeros(
                (0, nsta * ncomp, kept), dtype=complex
            )
            block_times[value] = tuple(
                array.starttime + float(time)
                for time in spectra.times(index)
            )

            if not singular.shape[0]:
                skipped.append(value)

        if skipped:
            warn_external(
                f"no complete block at {len(skipped)} frequencies "
                f"({min(skipped):.4g} - {max(skipped):.4g} Hz): the data "
                "are too short for the windows",
                FDDWarning,
            )

        return FDDResult._from_parts(
            config, block_times, values, vectors, stations, components,
            coordinates, array.units,
        )


def _direction(name: str, azimuth: Optional[float]) -> tuple:
    """Direction (east, north, up) of a component."""
    if name in _DIRECTIONS:
        return _DIRECTIONS[name]

    angle = math.radians(azimuth)

    if name == "radial":
        return (math.sin(angle), math.cos(angle), 0.0)

    return (math.cos(angle), -math.sin(angle), 0.0)


def _frequencies(config: FDDConfig, delta: float) -> np.ndarray:
    """Analysis frequencies (FFT frequencies of fixed windows if unset)."""
    frequencies = config.frequencies

    if frequencies.is_sampled:
        return frequencies.targets()

    if config.window.type is WindowType.CYCLES:
        raise ValueError(
            "windows in cycles need the number of frequencies (nfreq) or "
            "explicit frequencies"
        )

    grid = fourier.frequency_axis(delta, config.window.size(delta))
    grid = grid[
        (grid >= frequencies.minimum) & (grid <= frequencies.maximum)
        & (grid > 0.0) & (grid < 0.5 / delta)
    ]

    if grid.size == 0:
        raise ValueError("no FFT frequency falls in the frequency range")

    return grid
