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
Array data.

The array-processing methods of :mod:`shakelab.exploration.yarray`
work on an :class:`ArrayData`: the waveforms of a set of stations, with
the position of every station in a local metric frame and the
orientation of every recorded component. Each station
(:class:`StationData`) has one, two or three components on a common
time axis.

An ``ArrayData`` is a container filled in one go or progressively:

.. code-block:: python

    array = ArrayData("array.json")              # array file
    array = ArrayData(collection)                # waveforms with metadata
    array = ArrayData()
    array.read("array.json")
    array.add(records, coordinates={...})
    array.add(station)
    array.remove("XX.S01.")

Waveforms come from any ShakeLab source (files read with
:func:`shakelab.signals.io.reader`, the buffer of a SeedLink client, a
FDSN request): a :class:`~shakelab.signals.base.StreamCollection`, a
:class:`~shakelab.signals.base.Stream`, records or ready
:class:`StationData`.

Conventions
-----------
Coordinates are East, North, Up in metres (``x = East``, ``y = North``,
``z = Up``). When they are computed from geographic positions, the
frame is the local ENU frame of :attr:`ArrayData.frame`; explicit
coordinates are in a local frame chosen by the user. Orientations are
unit vectors of the positive direction of each component in the same
frame. Azimuths are in degrees clockwise from north.

Design choices
--------------
The stations of an array need not share their time axis: the time
window common to all of them is taken when the array is processed
(:meth:`ArrayData.window`, used by the beamformer), so that stations
can be added one at a time from different sources. Their sampling
interval and units must be the same; the samples of the processed
window must be simultaneous (no resampling is performed).

The origin of the ENU frame is fixed: given by ``reference`` (or by
the array file), otherwise the position of the first station added
with a geographic position; adding stations does not move it.
:meth:`ArrayData.set_reference` changes it and recomputes the
coordinates.

``StationData`` is immutable; an ``ArrayData`` changes only by adding
and removing stations. The records are stored by reference: the
builders create new records, so that later changes of the input data do
not affect the array.
"""

from __future__ import annotations

import dataclasses
import math
import os
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from numbers import Real
from typing import Any, Optional, Union

import numpy as np

from shakelab.libutils.geodeticN import (
    MetricFrame,
    WgsPoint,
    geo_to_enu,
    relative_longitude,
)
from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import warn_external
from shakelab.signals.base import (
    SNAP_TOLERANCE,
    Record,
    Stream,
    StreamCollection,
)
from shakelab.signals.codes import match_sid

__all__ = [
    "ArrayData",
    "ArrayDataWarning",
    "StationData",
]

# Unit vectors of the canonical components, as (east, north, up).
_COMPONENT_VECTORS = {
    "E": (1.0, 0.0, 0.0),
    "N": (0.0, 1.0, 0.0),
    "Z": (0.0, 0.0, 1.0),
}

_EAST = np.array(_COMPONENT_VECTORS["E"])
_NORTH = np.array(_COMPONENT_VECTORS["N"])
_VERTICAL = np.array(_COMPONENT_VECTORS["Z"])

# Order of the components within a station.
_COMPONENT_ORDER = "ZNE"

# Smallest singular value of a usable orientation matrix.
_ORIENTATION_TOLERANCE = 1e-8

# Relative tolerance on equal sampling intervals.
_SAMPLING_RTOL = 1e-9

# Largest distance, in metres, between the positions of the components
# of one station.
_LOCATION_TOLERANCE = 0.1

# Accepted values of the ``gaps`` arguments.
_GAP_POLICIES = ("drop", "raise")

CoordinateCovariance = Union[None, Real, Sequence, np.ndarray, Mapping]


class ArrayDataWarning(UserWarning):
    """Stations excluded while building an :class:`ArrayData`."""


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _readonly(value: Any) -> np.ndarray:
    """Return a read-only float copy of an array."""
    array = np.array(value, dtype=float)
    array.setflags(write=False)
    return array


def _as_coordinates(value: Any, context: str) -> np.ndarray:
    """Return a 3-D coordinate vector (2-D values get zero height)."""
    coordinates = np.asarray(value, dtype=float)

    if coordinates.shape == (2,):
        coordinates = np.append(coordinates, 0.0)

    if coordinates.shape != (3,):
        raise ValueError(
            f"{context}: coordinates must contain two or three values"
        )

    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{context}: coordinates must be finite")

    return _readonly(coordinates)


def _as_orientation(value: Any, ncomp: int, context: str) -> np.ndarray:
    """Return a normalized orientation matrix ``(ncomp, 3)``."""
    orientation = np.asarray(value, dtype=float)

    if orientation.shape != (ncomp, 3):
        raise ValueError(
            f"{context}: orientation must have shape ({ncomp}, 3), got "
            f"{orientation.shape}"
        )

    if not np.all(np.isfinite(orientation)):
        raise ValueError(f"{context}: orientation must be finite")

    norms = np.linalg.norm(orientation, axis=1)

    if np.any(norms == 0.0):
        raise ValueError(
            f"{context}: every orientation must be a non-zero vector"
        )

    orientation = orientation / norms[:, None]
    singular = np.linalg.svd(orientation, compute_uv=False)

    if singular[-1] <= _ORIENTATION_TOLERANCE:
        raise ValueError(
            f"{context}: the component orientations are not linearly "
            "independent"
        )

    return _readonly(orientation)


def _as_covariance(value: Any, context: str) -> Optional[np.ndarray]:
    """
    Return a 3 x 3 coordinate covariance matrix, or None.

    A scalar is the variance of each horizontal coordinate (isotropic,
    no vertical uncertainty); a 2 x 2 matrix is horizontal.
    """
    if value is None:
        return None

    covariance = np.asarray(value, dtype=float)

    if covariance.ndim == 0:
        covariance = np.diag([float(covariance)] * 2)

    if covariance.shape == (2, 2):
        expanded = np.zeros((3, 3))
        expanded[:2, :2] = covariance
        covariance = expanded

    if covariance.shape != (3, 3):
        raise ValueError(
            f"{context}: a coordinate covariance must be a scalar or a "
            "(2, 2) or (3, 3) matrix"
        )

    if not np.all(np.isfinite(covariance)):
        raise ValueError(f"{context}: the covariance must be finite")

    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-12):
        raise ValueError(f"{context}: the covariance must be symmetric")

    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues = np.linalg.eigvalsh(covariance)
    tolerance = (
        3 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(eigenvalues))))
    )

    if np.any(eigenvalues < -tolerance):
        raise ValueError(
            f"{context}: the covariance must be positive semidefinite"
        )

    return _readonly(covariance)


def _as_direction(value: Any) -> np.ndarray:
    """Return a unit vector (east, north, up)."""
    direction = np.asarray(value, dtype=float)

    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("a direction must contain three finite values")

    norm = np.linalg.norm(direction)

    if norm == 0.0:
        raise ValueError("a direction must be a non-zero vector")

    return direction / norm


def _check_azimuth(azimuth: Any) -> float:
    """Return a finite azimuth in degrees."""
    if (
        not isinstance(azimuth, Real)
        or isinstance(azimuth, (bool, np.bool_))
    ):
        raise TypeError("azimuth must be a real number")

    if not math.isfinite(azimuth):
        raise ValueError("azimuth must be finite")

    return float(azimuth)


def _check_records(sid: str, records: tuple) -> None:
    """Check that the records of a station share their time axis."""
    reference = records[0]

    if reference.delta is None or reference.time is None:
        raise ValueError(
            f"station {sid}: records need a sampling interval and a "
            "start time"
        )

    for record in records:
        if record.nsamp == 0:
            raise ValueError(f"station {sid}: empty record")

        if record.delta is None or record.time is None:
            raise ValueError(
                f"station {sid}: records need a sampling interval and "
                "a start time"
            )

        if not np.all(np.isfinite(record.data)):
            raise ValueError(f"station {sid}: non-finite samples")

        if not math.isclose(
            record.delta, reference.delta, rel_tol=_SAMPLING_RTOL
        ):
            raise ValueError(
                f"station {sid}: records with different sampling "
                "intervals"
            )

        if (
            record.nsamp != reference.nsamp
            or record.starttime != reference.starttime
        ):
            raise ValueError(
                f"station {sid}: records on different time axes "
                f"({record.starttime}, {record.nsamp} samples and "
                f"{reference.starttime}, {reference.nsamp} samples)"
            )


def _common_units(records: Iterable[Record], context: str) -> Optional[str]:
    """Return the common units of records (None if all are unknown)."""
    units = {record.head.units for record in records} - {None}

    if len(units) > 1:
        names = ", ".join(sorted(str(unit) for unit in units))
        raise ValueError(f"{context}: records in different units ({names})")

    return units.pop() if units else None


def _record_orientation(record: Record, sid: str) -> np.ndarray:
    """Orientation of a record, from its metadata or its channel code."""
    if record.orientation is not None:
        return record.orientation.vector

    code = (record.sid or "").split(".")[-1][-1:]

    if code in _COMPONENT_VECTORS:
        return np.array(_COMPONENT_VECTORS[code])

    raise ValueError(
        f"station {sid}: the orientation of channel {record.sid} is "
        "unknown (no metadata and not a Z, N or E component)"
    )


def _component_key(record: Record) -> tuple:
    """Sort key of the components of a station: Z, N, E, then others."""
    code = (record.sid or "").split(".")[-1][-1:]
    index = _COMPONENT_ORDER.find(code) if code else -1
    return (index if index >= 0 else len(_COMPONENT_ORDER), code)


# ----------------------------------------------------------------------------
# Station data
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class StationData:
    """
    Waveforms, position and orientation of one array station.

    Parameters
    ----------
    sid : str
        Station identifier (``NET.STA.LOC`` when built from ShakeLab
        waveforms).
    records : sequence of Record
        One to three records on the same time axis (same start time,
        sampling interval and number of samples).
    coordinates : array-like
        Position ``(x, y[, z])`` in metres in the frame of the array
        (East, North, Up); a missing ``z`` is zero.
    orientation : array-like, optional
        Matrix ``(n_components, 3)`` whose rows are the directions of
        the positive motion of the records, as (east, north, up). They
        are normalized; they must be linearly independent. Default: the
        orientation of each record (:attr:`Record.orientation`), or its
        channel code (Z, N, E).
    coordinate_covariance : float or array-like, optional
        Uncertainty of the coordinates, in square metres: a scalar
        (variance of each horizontal coordinate), or a (2, 2)
        horizontal or (3, 3) covariance matrix.
    location : WgsPoint, optional
        Geographic position from which the coordinates derive (None for
        coordinates given in a local frame).

    Raises
    ------
    TypeError
        If a record is not a Record or ``location`` is not a WgsPoint.
    ValueError
        If the records are not one to three, are empty, contain
        non-finite samples, are not on the same time axis or have
        different units, or if coordinates, orientation or covariance
        are invalid.

    Notes
    -----
    The records are not copied.
    """

    sid: str
    records: Sequence[Record]
    coordinates: Any
    orientation: Any = None
    coordinate_covariance: Any = None
    location: Optional[WgsPoint] = None

    def __post_init__(self) -> None:
        if not isinstance(self.sid, str) or not self.sid.strip():
            raise ValueError("the station identifier must be a string")

        sid = self.sid.strip()
        records = tuple(self.records)

        if not 1 <= len(records) <= 3:
            raise ValueError(
                f"station {sid}: one to three records are required"
            )

        if not all(isinstance(record, Record) for record in records):
            raise TypeError(f"station {sid}: records must be Record")

        if self.location is not None and not isinstance(
            self.location, WgsPoint
        ):
            raise TypeError(f"station {sid}: location must be a WgsPoint")

        _check_records(sid, records)
        _common_units(records, f"station {sid}")

        orientation = self.orientation

        if orientation is None:
            orientation = [_record_orientation(rec, sid) for rec in records]

        context = f"station {sid}"
        object.__setattr__(self, "sid", sid)
        object.__setattr__(self, "records", records)
        object.__setattr__(
            self, "coordinates", _as_coordinates(self.coordinates, context)
        )
        object.__setattr__(
            self,
            "orientation",
            _as_orientation(orientation, len(records), context),
        )
        object.__setattr__(
            self,
            "coordinate_covariance",
            _as_covariance(self.coordinate_covariance, context),
        )

    def __len__(self) -> int:
        """Return the number of components."""
        return len(self.records)

    def __iter__(self) -> Iterator[Record]:
        """Iterate over the records."""
        return iter(self.records)

    # -- properties -----------------------------------------------------

    @property
    def ncomp(self) -> int:
        """Number of components."""
        return len(self.records)

    @property
    def nsamp(self) -> int:
        """Number of samples of each component."""
        return self.records[0].nsamp

    @property
    def delta(self) -> float:
        """Sampling interval in seconds."""
        return float(self.records[0].delta)

    @property
    def rate(self) -> float:
        """Sampling rate in Hz."""
        return 1.0 / self.delta

    @property
    def starttime(self) -> Date:
        """Time of the first sample."""
        return self.records[0].starttime

    @property
    def endtime(self) -> Date:
        """Time of the last sample."""
        return self.records[0].endtime

    @property
    def units(self) -> Optional[str]:
        """Common units of the records, or None if unknown."""
        return _common_units(self.records, f"station {self.sid}")

    @property
    def data(self) -> np.ndarray:
        """Samples as a matrix ``(n_components, n_samples)``."""
        return np.stack([np.asarray(rec.data, dtype=float)
                         for rec in self.records])

    # -- projections ----------------------------------------------------

    def projection_weights(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> np.ndarray:
        """
        Return the weights that combine the components into a direction.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up); it is normalized.
        tolerance : float, optional
            Largest accepted residual of the reconstructed unit vector.

        Returns
        -------
        numpy.ndarray
            One weight per component.

        Raises
        ------
        ValueError
            If the direction cannot be reconstructed from the
            components (e.g. a horizontal direction from a vertical
            sensor).
        """
        target = _as_direction(direction)
        weights = np.linalg.lstsq(self.orientation.T, target, rcond=None)[0]
        residual = np.linalg.norm(self.orientation.T @ weights - target)

        if residual > tolerance:
            raise ValueError(
                f"station {self.sid}: direction {target.tolist()} cannot "
                "be reconstructed from the recorded components"
            )

        return weights

    def supports(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> bool:
        """
        Return whether a direction can be reconstructed.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up).
        tolerance : float, optional
            Largest accepted residual (see :meth:`projection_weights`).

        Returns
        -------
        bool
            True if the direction is observable.
        """
        try:
            self.projection_weights(direction, tolerance)
        except ValueError:
            return False

        return True

    def project(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> np.ndarray:
        """
        Return the motion along a direction.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up).
        tolerance : float, optional
            Largest accepted residual (see :meth:`projection_weights`).

        Returns
        -------
        numpy.ndarray
            Samples ``(n_samples,)``.
        """
        return self.projection_weights(direction, tolerance) @ self.data

    def matrix(
        self,
        directions: Iterable,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> np.ndarray:
        """
        Return the motion along several directions.

        Parameters
        ----------
        directions : iterable of array-like
            Target directions (east, north, up).
        tolerance : float, optional
            Largest accepted residual (see :meth:`projection_weights`).

        Returns
        -------
        numpy.ndarray
            Samples ``(n_directions, n_samples)``.

        Raises
        ------
        ValueError
            If no direction is given or one cannot be reconstructed.
        """
        directions = tuple(directions)

        if not directions:
            raise ValueError("at least one direction is required")

        weights = np.stack([
            self.projection_weights(direction, tolerance)
            for direction in directions
        ])

        return weights @ self.data

    def vertical(self) -> np.ndarray:
        """
        Return the vertical motion (positive upward).

        Returns
        -------
        numpy.ndarray
            Samples ``(n_samples,)``.
        """
        return self.project(_VERTICAL)

    def horizontal(self) -> np.ndarray:
        """
        Return the east and north motion.

        Returns
        -------
        numpy.ndarray
            Samples ``(2, n_samples)``.
        """
        return self.matrix((_EAST, _NORTH))

    def canonical(self) -> np.ndarray:
        """
        Return the east, north and vertical motion.

        Returns
        -------
        numpy.ndarray
            Samples ``(3, n_samples)``.
        """
        return self.matrix((_EAST, _NORTH, _VERTICAL))

    def radial(self, azimuth: float) -> np.ndarray:
        """
        Return the horizontal motion along a propagation azimuth.

        Parameters
        ----------
        azimuth : float
            Propagation azimuth, degrees clockwise from north.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_samples,)``.
        """
        angle = math.radians(_check_azimuth(azimuth))
        return self.project((math.sin(angle), math.cos(angle), 0.0))

    def transverse(self, azimuth: float) -> np.ndarray:
        """
        Return the horizontal motion normal to a propagation azimuth.

        The transverse direction is the radial one rotated by 90
        degrees clockwise (seen from above): azimuth ``azimuth + 90``.

        Parameters
        ----------
        azimuth : float
            Propagation azimuth, degrees clockwise from north.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_samples,)``.
        """
        angle = math.radians(_check_azimuth(azimuth))
        return self.project((math.cos(angle), -math.sin(angle), 0.0))

    # -- conversion -----------------------------------------------------

    def cut(
        self,
        starttime: Any,
        endtime: Any,
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> Optional[StationData]:
        """
        Return the station on a time window.

        Parameters
        ----------
        starttime, endtime : Date or str
            Window, bounds included.
        snap_tolerance : float, optional
            Times within this fraction of a sampling interval from a
            sample coincide with it.

        Returns
        -------
        StationData or None
            Station with new records on the window, None if the records
            do not cover it.
        """
        start = Date(starttime)
        end = Date(endtime)
        tolerance = snap_tolerance * self.delta

        if (
            self.starttime - start > tolerance
            or end - self.endtime > tolerance
        ):
            return None

        return dataclasses.replace(self, records=[
            record.extract(start, end, snap_tolerance=snap_tolerance)
            for record in self.records
        ])


# ----------------------------------------------------------------------------
# Array data
# ----------------------------------------------------------------------------


class ArrayData:
    """
    Stations of a seismic array.

    Parameters
    ----------
    source : str, path-like or waveforms, optional
        An array file (JSON, see :meth:`read`), or the waveforms of the
        stations to add (see :meth:`add`). Default: an empty array.
    reference : MetricFrame, WgsPoint or sequence of float, optional
        Origin of the ENU frame: a frame, or a point ``(lon, lat[,
        elevation])``. Default: the reference of the array file, or the
        position of the first station added with a geographic position.
        For an array in local coordinates it is the geographic position
        of their origin (optional).
    **options
        Options of :meth:`add` (only with waveforms).

    Raises
    ------
    TypeError
        If ``options`` are given with an array file.

    Examples
    --------
    >>> array = ArrayData("array.json")
    >>> array = ArrayData(collection, select="XX.*..HH?")
    >>> array = ArrayData()
    >>> array.add(records, coordinates={"XX.A01.": (0.0, 0.0)})
    """

    def __init__(
        self,
        source: Any = None,
        reference: Any = None,
        **options: Any,
    ) -> None:
        self._stations: dict[str, StationData] = {}
        self._frame: Optional[MetricFrame] = None
        self._excluded: list = []

        if reference is not None:
            self._frame = _frame(reference, [])

        if source is None:
            if options:
                raise TypeError("options need waveforms to add")

        elif _is_path(source):
            if options:
                raise TypeError(
                    "an array file takes no options: "
                    + ", ".join(sorted(options))
                )

            self.read(source)

        else:
            self.add(source, **options)

    def __len__(self) -> int:
        """Return the number of stations."""
        return len(self._stations)

    def __iter__(self) -> Iterator[StationData]:
        """Iterate over the stations."""
        return iter(tuple(self._stations.values()))

    def __contains__(self, sid: Any) -> bool:
        """Return whether a station identifier is in the array."""
        return sid in self._stations

    def __getitem__(self, item: Any) -> Any:
        """Return a station by position or identifier."""
        if isinstance(item, str):
            return self._stations[item]

        return self.stations[item]

    def __repr__(self) -> str:
        if not self._stations:
            return "ArrayData(no stations)"

        return (
            f"ArrayData({len(self)} stations at {self.rate:g} Hz, common "
            f"window {self.starttime} - {self.endtime})"
        )

    # -- filling --------------------------------------------------------

    def read(self, source: Any) -> None:
        """
        Add the stations of an array file.

        The file (JSON) lists the stations with their local coordinates
        and the waveform file of each component, for data without (or
        with unusable) metadata, and optionally the geographic position
        of the origin of the coordinates; see :mod:`.io` for the format.

        Parameters
        ----------
        source : str or path-like
            Array file.

        Raises
        ------
        FileNotFoundError
            If the array file or a waveform file does not exist.
        ValueError
            If the file is invalid, its reference differs from that of
            the array, or a station cannot be added (see :meth:`add`).
        """
        from shakelab.exploration.yarray.io import read_array_file

        stations, reference = read_array_file(source)

        if reference is not None:
            frame = _frame(reference, [])

            if self._frame is not None and self._frame != frame:
                raise ValueError(
                    f"{source}: the reference of the file differs from "
                    "that of the array"
                )

        self._insert(stations, [])

        if reference is not None:
            self._frame = frame

    def add(
        self,
        data: Any,
        select: Any = None,
        starttime: Any = None,
        endtime: Any = None,
        coordinates: Optional[Mapping] = None,
        orientations: Optional[Mapping] = None,
        coordinate_covariance: CoordinateCovariance = None,
        gaps: str = "drop",
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> None:
        """
        Add stations.

        Waveforms are grouped by station (``NET.STA.LOC``, up to three
        channels of one instrument); records of the same channel are
        consecutive segments of it. Each station is cut to the time
        window covered by all its channels, or to the given window.

        Parameters
        ----------
        data : StationData, Record, Stream, StreamCollection or iterable
            Stations, or waveforms with source identifiers
            ``NET.STA.LOC.CHA``. A StreamCollection is not modified.
        select : str or list of str, optional
            Channel pattern(s) ``NET.STA.LOC.CHA`` (see
            :mod:`shakelab.signals.codes`). Default: every channel.
        starttime, endtime : Date or str, optional
            Time window, bounds included. Default: the interval covered
            by all the channels of each station.
        coordinates : mapping, optional
            Station positions ``{sid: (x, y[, z])}`` in metres, in the
            frame of the array, for every added station (``sid`` is
            ``NET.STA.LOC``). Default: the positions of the records
            (:attr:`Record.location`), in the ENU frame of the array.
        orientations : mapping, optional
            Orientation matrices ``{sid: matrix}`` overriding those of
            the metadata for some stations, rows in the order Z, N, E
            and then by channel code.
        coordinate_covariance : float, array-like or mapping, optional
            Coordinate uncertainty, in square metres: one value for all
            stations (see :class:`StationData`) or a mapping by station.
        gaps : {"drop", "raise"}, optional
            Stations whose channels do not cover the window with
            continuous data: excluded with an :class:`ArrayDataWarning`
            (``"drop"``, default; listed in :attr:`excluded`) or an
            error (``"raise"``).
        snap_tolerance : float, optional
            Times within this fraction of a sampling interval from a
            sample coincide with it.

        Raises
        ------
        TypeError
            If ``data`` has the wrong type, or options are given with
            StationData.
        ValueError
            If no channel is selected, a station has channels of
            different instruments or more than three, a station has a
            gap with ``gaps="raise"``, positions or orientations are
            missing, an identifier is already in the array, or the
            sampling interval or the units differ from those of the
            array.
        """
        if gaps not in _GAP_POLICIES:
            raise ValueError("gaps must be 'drop' or 'raise'")

        stations, groups = _split_input(data, select)

        if stations is not None:
            given = {
                "select": select, "starttime": starttime,
                "endtime": endtime, "coordinates": coordinates,
                "orientations": orientations,
                "coordinate_covariance": coordinate_covariance,
            }
            used = sorted(name for name, value in given.items()
                          if value is not None)

            if used:
                raise TypeError(
                    "options apply to waveforms, not to StationData: "
                    + ", ".join(used)
                )

            self._insert(stations, [])
            return

        if coordinates is not None and not isinstance(coordinates, Mapping):
            raise TypeError("coordinates must be a mapping {sid: (x, y[, z])}")

        if orientations is not None and not isinstance(orientations, Mapping):
            raise TypeError("orientations must be a mapping {sid: matrix}")

        channels, excluded = _station_channels(
            groups, starttime, endtime, gaps, snap_tolerance
        )

        if not channels:
            _warn_excluded(excluded)
            raise ValueError("no station covers the time window")

        if coordinates is not None:
            missing = sorted(set(channels) - set(coordinates))

            if missing:
                raise ValueError(
                    "coordinates are missing for stations "
                    + ", ".join(missing)
                )

            positions = {sid: (coordinates[sid], None) for sid in channels}

        else:
            positions = self._geographic_positions(channels)

        stations = [
            StationData(
                sid=sid,
                records=records,
                coordinates=positions[sid][0],
                orientation=_station_value(orientations, sid),
                coordinate_covariance=_station_value(
                    coordinate_covariance, sid
                ),
                location=positions[sid][1],
            )
            for sid, records in channels.items()
        ]
        self._insert(stations, excluded)

    def remove(self, sid: Any) -> None:
        """
        Remove stations.

        Parameters
        ----------
        sid : str or list of str
            Station identifier(s).

        Raises
        ------
        KeyError
            If a station is not in the array (nothing is removed).
        """
        sids = [sid] if isinstance(sid, str) else list(sid)
        missing = [item for item in sids if item not in self._stations]

        if missing:
            raise KeyError(", ".join(str(item) for item in missing))

        for item in sids:
            del self._stations[item]

    def set_reference(self, reference: Any = None) -> None:
        """
        Change the origin of the ENU frame.

        The coordinates of the stations with a geographic position are
        recomputed; for an array in local coordinates the reference is
        only the geographic position of their origin.

        Parameters
        ----------
        reference : MetricFrame, WgsPoint or sequence of float, optional
            New origin. Default: the mean position of the stations
            (centre of the array).

        Raises
        ------
        ValueError
            If only some stations have a geographic position, or the
            default is requested for an array without geographic
            positions.
        """
        stations = list(self._stations.values())
        located = [st for st in stations if st.location is not None]

        if located and len(located) != len(stations):
            missing = sorted(st.sid for st in stations if st.location is None)
            raise ValueError(
                "stations " + ", ".join(missing) + " have no geographic "
                "position: the coordinates cannot be recomputed"
            )

        if not located:
            if reference is None:
                raise ValueError(
                    "the stations have no geographic position: give the "
                    "reference"
                )

            self._frame = _frame(reference, [])
            return

        frame = _frame(reference, [st.location for st in located])
        self._frame = frame
        self._stations = {
            st.sid: dataclasses.replace(
                st, coordinates=geo_to_enu(st.location, frame)
            )
            for st in stations
        }

    def copy(self) -> ArrayData:
        """
        Return a copy of the array.

        Returns
        -------
        ArrayData
            New array with the same stations (immutable, shared).
        """
        return self._derived(list(self._stations.values()), [])

    def window(
        self,
        starttime: Any = None,
        endtime: Any = None,
        gaps: str = "drop",
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> ArrayData:
        """
        Return the array on a common time window.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Window, bounds included. Default: the window common to all
            the stations (:attr:`starttime`, :attr:`endtime`).
        gaps : {"drop", "raise"}, optional
            Stations that do not cover the window: excluded with an
            :class:`ArrayDataWarning` (default) or an error.
        snap_tolerance : float, optional
            Times within this fraction of a sampling interval from a
            sample coincide with it.

        Returns
        -------
        ArrayData
            Synchronized array with new records.

        Raises
        ------
        ValueError
            If the array is empty, the window is empty, a station does
            not cover it with ``gaps="raise"``, no station is left or
            the samples of the stations are not simultaneous.
        """
        if gaps not in _GAP_POLICIES:
            raise ValueError("gaps must be 'drop' or 'raise'")

        if not self._stations:
            raise ValueError("the array has no station")

        start = self.starttime if starttime is None else Date(starttime)
        end = self.endtime if endtime is None else Date(endtime)

        if end < start:
            raise ValueError(
                "the stations have no common time window"
                if starttime is None and endtime is None
                else f"empty time window ({start} - {end})"
            )

        stations = []
        excluded = []

        for station in self._stations.values():
            cut = station.cut(start, end, snap_tolerance)

            if cut is not None:
                stations.append(cut)
                continue

            reason = (
                f"data from {station.starttime} to {station.endtime} do "
                f"not cover {start} - {end}"
            )

            if gaps == "raise":
                raise ValueError(f"station {station.sid}: {reason}")

            excluded.append((station.sid, reason))

        _warn_excluded(excluded)

        if not stations:
            raise ValueError("no station covers the time window")

        reference = stations[0]

        for station in stations[1:]:
            if (
                station.nsamp != reference.nsamp
                or station.starttime != reference.starttime
            ):
                raise ValueError(
                    "stations on different time axes: "
                    f"{reference.sid} from {reference.starttime} "
                    f"({reference.nsamp} samples), {station.sid} from "
                    f"{station.starttime} ({station.nsamp} samples); "
                    "the samples must be simultaneous (no resampling "
                    "is performed)"
                )

        return self._derived(stations, excluded)

    # -- properties -----------------------------------------------------

    @property
    def stations(self) -> tuple[StationData, ...]:
        """Stations, in the order they were added."""
        return tuple(self._stations.values())

    @property
    def frame(self) -> Optional[MetricFrame]:
        """
        Local ENU frame of the coordinates.

        None for an array in local coordinates without a geographic
        reference (see :meth:`set_reference`).
        """
        return self._frame

    @property
    def excluded(self) -> tuple[tuple, ...]:
        """``(sid, reason)`` of the stations left out by the builders."""
        return tuple(self._excluded)

    @property
    def nsta(self) -> int:
        """Number of stations."""
        return len(self._stations)

    @property
    def sids(self) -> tuple[str, ...]:
        """Station identifiers, in the order they were added."""
        return tuple(self._stations)

    @property
    def coordinates(self) -> np.ndarray:
        """Station coordinates ``(n_stations, 3)``, in metres."""
        if not self._stations:
            return np.empty((0, 3))

        return np.stack([st.coordinates for st in self._stations.values()])

    @property
    def coordinate_covariances(self) -> Optional[np.ndarray]:
        """
        Coordinate covariances ``(n_stations, 3, 3)``.

        None if no station has one; zeros for the stations without one
        otherwise.
        """
        values = [st.coordinate_covariance for st in self._stations.values()]

        if all(value is None for value in values):
            return None

        return np.stack([
            np.zeros((3, 3)) if value is None else value for value in values
        ])

    @property
    def synchronized(self) -> bool:
        """Whether all the stations share their time axis."""
        stations = self.stations

        return bool(stations) and all(
            st.starttime == stations[0].starttime
            and st.nsamp == stations[0].nsamp
            for st in stations[1:]
        )

    @property
    def delta(self) -> Optional[float]:
        """Sampling interval in seconds (None for an empty array)."""
        return self.stations[0].delta if self._stations else None

    @property
    def rate(self) -> Optional[float]:
        """Sampling rate in Hz (None for an empty array)."""
        return self.stations[0].rate if self._stations else None

    @property
    def starttime(self) -> Optional[Date]:
        """Start of the window common to all the stations."""
        if not self._stations:
            return None

        return max(st.starttime for st in self._stations.values())

    @property
    def endtime(self) -> Optional[Date]:
        """End of the window common to all the stations."""
        if not self._stations:
            return None

        return min(st.endtime for st in self._stations.values())

    @property
    def nsamp(self) -> int:
        """Number of samples of the common window (0 if empty)."""
        if self.synchronized:
            return self.stations[0].nsamp

        if not self._stations or self.endtime < self.starttime:
            return 0

        return int(round((self.endtime - self.starttime) / self.delta)) + 1

    @property
    def duration(self) -> float:
        """Length of the common window, in seconds."""
        return max(self.nsamp - 1, 0) * (self.delta or 0.0)

    @property
    def units(self) -> Optional[str]:
        """Common units of the records, or None if unknown."""
        return _common_units(
            (rec for st in self._stations.values() for rec in st.records),
            "array",
        )

    # -- projections ----------------------------------------------------

    def supports(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> bool:
        """
        Return whether every station can reconstruct a direction.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up).
        tolerance : float, optional
            Largest accepted residual (see
            :meth:`StationData.projection_weights`).

        Returns
        -------
        bool
            True if all stations observe the direction.
        """
        return not self.unsupported_stations(direction, tolerance)

    def unsupported_stations(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> tuple[str, ...]:
        """
        Return the stations that cannot reconstruct a direction.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up).
        tolerance : float, optional
            Largest accepted residual.

        Returns
        -------
        tuple of str
            Station identifiers.
        """
        return tuple(
            st.sid for st in self._stations.values()
            if not st.supports(direction, tolerance)
        )

    def project(
        self,
        direction: Any,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> np.ndarray:
        """
        Return the motion of every station along a direction.

        Parameters
        ----------
        direction : array-like
            Target direction (east, north, up).
        tolerance : float, optional
            Largest accepted residual.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_stations, n_samples)``.

        Raises
        ------
        ValueError
            If the array is not synchronized (see :meth:`window`) or
            some stations cannot reconstruct the direction.
        """
        self._require_support((direction,), tolerance)
        return np.stack([
            st.project(direction, tolerance) for st in self._stations.values()
        ])

    def matrix(
        self,
        directions: Iterable,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> np.ndarray:
        """
        Return the motion of every station along several directions.

        Parameters
        ----------
        directions : iterable of array-like
            Target directions (east, north, up).
        tolerance : float, optional
            Largest accepted residual.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_stations, n_directions, n_samples)``.

        Raises
        ------
        ValueError
            If no direction is given, the array is not synchronized or
            some stations cannot reconstruct a direction.
        """
        directions = tuple(directions)

        if not directions:
            raise ValueError("at least one direction is required")

        self._require_support(directions, tolerance)
        return np.stack([
            st.matrix(directions, tolerance)
            for st in self._stations.values()
        ])

    def vertical(self) -> np.ndarray:
        """
        Return the vertical motion of every station.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_stations, n_samples)``.
        """
        return self.project(_VERTICAL)

    def horizontal(self) -> np.ndarray:
        """
        Return the east and north motion of every station.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_stations, 2, n_samples)``.
        """
        return self.matrix((_EAST, _NORTH))

    def canonical(self) -> np.ndarray:
        """
        Return the east, north and vertical motion of every station.

        Returns
        -------
        numpy.ndarray
            Samples ``(n_stations, 3, n_samples)``.
        """
        return self.matrix((_EAST, _NORTH, _VERTICAL))

    def info(self) -> None:
        """Print a summary of the array."""
        print(self.summary())

    def summary(self) -> str:
        """
        Return a summary of the array as text.

        Returns
        -------
        str
            Common window, frame, stations and excluded stations.
        """
        frame = (
            "local, no geographic reference" if self._frame is None
            else "ENU from lon {:.6f}, lat {:.6f}, h {:.1f} m".format(
                *self._frame.ref_geo
            )
        )

        if not self._stations:
            return f"ArrayData: no stations\n  frame     : {frame}"

        synchronized = self.synchronized
        lines = [
            f"ArrayData: {self.nsta} stations"
            + ("" if synchronized else " (different time axes)"),
            f"  window    : {self.starttime} - {self.endtime} "
            f"({self.nsamp} samples at {self.rate:g} Hz)",
            f"  units     : {self.units}",
            f"  frame     : {frame}",
        ]

        for station in self._stations.values():
            x, y, z = station.coordinates
            line = (
                f"  {station.sid:<16} {station.ncomp} comp  "
                f"x {x:10.2f}  y {y:10.2f}  z {z:8.2f} m"
            )

            if not synchronized:
                line += f"  {station.starttime} - {station.endtime}"

            lines.append(line)

        for sid, reason in self._excluded:
            lines.append(f"  excluded {sid}: {reason}")

        return "\n".join(lines)

    # -- helpers --------------------------------------------------------

    def _derived(self, stations: list, excluded: list) -> ArrayData:
        """New array with the frame of this one and other stations."""
        array = ArrayData()
        array._frame = self._frame
        array._stations = {st.sid: st for st in stations}
        array._excluded = list(self._excluded) + list(excluded)
        return array

    def _insert(self, stations: list, excluded: list) -> None:
        """Add checked stations and record the excluded ones."""
        for station in stations:
            if not isinstance(station, StationData):
                raise TypeError("stations must be StationData")

        sids = [station.sid for station in stations]
        repeated = sorted({
            sid for sid in sids
            if sid in self._stations or sids.count(sid) > 1
        })

        if repeated:
            raise ValueError(
                "stations already in the array: " + ", ".join(repeated)
            )

        every = list(self._stations.values()) + list(stations)

        if every:
            first = every[0]

            for station in every[1:]:
                if not math.isclose(
                    station.delta, first.delta, rel_tol=_SAMPLING_RTOL
                ):
                    raise ValueError(
                        "stations with different sampling intervals "
                        f"({first.sid}, {station.sid}): no resampling is "
                        "performed"
                    )

            _common_units(
                (rec for station in every for rec in station.records),
                "array",
            )

        for station in stations:
            self._stations[station.sid] = station

        self._excluded.extend(excluded)
        _warn_excluded(excluded)

    def _geographic_positions(self, channels: dict) -> dict:
        """Coordinates and locations of new stations from metadata."""
        locations = {
            sid: _station_location(sid, records)
            for sid, records in channels.items()
        }
        missing = sorted(sid for sid, point in locations.items()
                         if point is None)

        if missing:
            raise ValueError(
                "the location of stations " + ", ".join(missing)
                + " is unknown: give metadata or coordinates"
            )

        if self._frame is None:
            local = sorted(st.sid for st in self._stations.values()
                           if st.location is None)

            if local:
                raise ValueError(
                    "the array is in local coordinates without a "
                    "geographic reference: give the coordinates of the "
                    "new stations, or a reference (set_reference)"
                )

            self._frame = _frame(next(iter(locations.values())), [])

        return {
            sid: (geo_to_enu(point, self._frame), point)
            for sid, point in locations.items()
        }

    def _require_support(self, directions: tuple, tolerance: float) -> None:
        """Raise if the array cannot give the motion along directions."""
        if not self._stations:
            raise ValueError("the array has no station")

        if not self.synchronized:
            raise ValueError(
                "the stations are on different time axes: take a common "
                "window with window()"
            )

        missing = sorted({
            sid for direction in directions
            for sid in self.unsupported_stations(direction, tolerance)
        })

        if missing:
            raise ValueError(
                "the requested direction is not observable at stations "
                + ", ".join(missing)
            )


# ----------------------------------------------------------------------------
# Builders
# ----------------------------------------------------------------------------


def _selected(sid: Optional[str], select: Any) -> bool:
    """Return whether a channel matches the selection."""
    if select is None:
        return True

    patterns = [select] if isinstance(select, str) else list(select)
    return any(match_sid(sid, pattern) for pattern in patterns)


def _station_sid(sid: Optional[str]) -> str:
    """Return ``NET.STA.LOC`` of a channel identifier."""
    parts = (sid or "").split(".")

    if len(parts) != 4:
        raise ValueError(
            f"channel identifier {sid!r} is not NET.STA.LOC.CHA"
        )

    return ".".join(parts[:3])


def _cut_channel(
    segments: list,
    start: Date,
    end: Date,
    snap_tolerance: float,
) -> Optional[Record]:
    """Return the channel cut to the window, or None if incomplete."""
    for record in segments:
        if record.nsamp == 0 or record.delta is None or record.time is None:
            continue

        tolerance = snap_tolerance * record.delta

        if (
            record.starttime - start <= tolerance
            and end - record.endtime <= tolerance
        ):
            return record.extract(start, end, snap_tolerance=snap_tolerance)

    return None


def _frame(reference: Any, locations: list) -> MetricFrame:
    """Return the ENU frame: given, or centred on the stations."""
    if isinstance(reference, MetricFrame):
        if reference.orientation != "enu":
            raise ValueError("the reference frame must be an ENU frame")

        return reference

    if reference is not None:
        if isinstance(reference, WgsPoint):
            reference = _geographic(reference)

        return MetricFrame(reference, "enu")

    lon0 = locations[0].longitude
    longitude = lon0 + float(np.mean([
        relative_longitude(point.longitude, lon0) for point in locations
    ]))
    latitude = float(np.mean([point.latitude for point in locations]))
    elevation = float(np.mean([point.elevation for point in locations]))

    return MetricFrame((longitude, latitude, elevation), "enu")


def _geographic(point: WgsPoint) -> tuple:
    """Return ``(longitude, latitude, elevation)`` of a point."""
    return (point.longitude, point.latitude, point.elevation or 0.0)


def _station_location(sid: str, records: list) -> Optional[WgsPoint]:
    """Common location of the records of a station, or None."""
    points = [rec.location for rec in records if rec.location is not None]

    if not points:
        return None

    frame = MetricFrame(_geographic(points[0]), "enu")

    for point in points[1:]:
        if np.linalg.norm(geo_to_enu(point, frame)) > _LOCATION_TOLERANCE:
            raise ValueError(
                f"station {sid}: the components have different locations"
            )

    return points[0]


def _station_value(value: Any, sid: str) -> Any:
    """Value of a per-station argument (mapping or common value)."""
    if isinstance(value, Mapping):
        return value.get(sid)

    return value


def _is_path(value: Any) -> bool:
    """Return whether a value is a file path."""
    return isinstance(value, (str, os.PathLike))


def _split_input(data: Any, select: Any) -> tuple:
    """
    Classify the data given to :meth:`ArrayData.add`.

    Returns ``(stations, None)`` for StationData, ``(None, groups)``
    for waveforms, where ``groups`` maps ``NET.STA.LOC`` to the list of
    the segments of each selected channel.
    """
    if isinstance(data, StationData):
        return [data], None

    if isinstance(data, StreamCollection):
        segments = [
            stream.records() for stream in data
            if stream.records() and _selected(stream.sid, select)
        ]
        return None, _groups(segments)

    if isinstance(data, Stream):
        items = list(data.records())
    elif isinstance(data, Record):
        items = [data]
    elif isinstance(data, Iterable) and not isinstance(data, (str, bytes)):
        items = list(data)
    else:
        raise TypeError(
            "data must be StationData, Record, Stream, StreamCollection "
            "or an iterable of them"
        )

    if items and all(isinstance(item, StationData) for item in items):
        return items, None

    channels: dict[str, list] = {}

    for record in items:
        if not isinstance(record, Record):
            raise TypeError(
                "data must be StationData or waveforms, not mixed"
            )

        if not record.sid:
            raise ValueError("records need a source identifier")

        if _selected(record.sid, select):
            channels.setdefault(record.sid, []).append(record)

    return None, _groups(list(channels.values()))


def _groups(segments: list) -> dict:
    """Group the segment lists of the channels by station."""
    if not segments:
        raise ValueError("no channel is selected")

    groups: dict[str, list] = {}

    for channel in segments:
        groups.setdefault(_station_sid(channel[0].sid), []).append(channel)

    return groups


def _station_channels(
    groups: dict,
    starttime: Any,
    endtime: Any,
    gaps: str,
    snap_tolerance: float,
) -> tuple[dict, list]:
    """Cut the channels of each station to its window."""
    channels: dict[str, list] = {}
    excluded = []

    for sid in sorted(groups):
        segments = [
            sorted(channel, key=lambda rec: rec.starttime)
            for channel in groups[sid]
        ]
        codes = {(channel[0].sid or "").split(".")[-1][:-1]
                 for channel in segments}

        if len(codes) > 1:
            raise ValueError(
                f"station {sid} has channels of different instruments "
                f"({', '.join(sorted(codes))}): choose them with select"
            )

        if len(segments) > 3:
            raise ValueError(
                f"station {sid} has more than three channels: choose "
                "them with select"
            )

        start = (
            max(channel[0].starttime for channel in segments)
            if starttime is None else Date(starttime)
        )
        end = (
            min(channel[-1].endtime for channel in segments)
            if endtime is None else Date(endtime)
        )
        records = []
        reason = None

        if end < start:
            reason = (
                "the channels have no common time window"
                if starttime is None and endtime is None
                else f"empty time window ({start} - {end})"
            )
        else:
            for channel in segments:
                record = _cut_channel(channel, start, end, snap_tolerance)

                if record is None:
                    reason = (
                        f"channel {channel[0].sid} does not cover "
                        f"{start} - {end} without gaps"
                    )
                    break

                records.append(record)

        if reason is not None:
            if gaps == "raise":
                raise ValueError(f"station {sid}: {reason}")

            excluded.append((sid, reason))
            continue

        channels[sid] = sorted(records, key=_component_key)

    return channels, excluded


def _warn_excluded(excluded: list) -> None:
    """Warn about the stations left out."""
    if excluded:
        warn_external(
            "stations excluded from the array: "
            + "; ".join(f"{sid} ({reason})" for sid, reason in excluded),
            ArrayDataWarning,
        )
