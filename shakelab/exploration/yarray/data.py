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
Synchronized array data.

The array-processing methods of :mod:`shakelab.exploration.yarray`
work on an :class:`ArrayData`: a snapshot of the waveforms of a set of
stations on a common time axis, with the position of every station in
a local metric frame and the orientation of every recorded component.
Each station (:class:`StationData`) has one, two or three components.

An ``ArrayData`` is built from ShakeLab waveforms, whatever their
origin (files read with :func:`shakelab.signals.io.reader`, the buffer
of a SeedLink client, a FDSN request):

- :meth:`ArrayData.from_collection` from a
  :class:`~shakelab.signals.base.StreamCollection`;
- :meth:`ArrayData.from_records` from a sequence of
  :class:`~shakelab.signals.base.Record`.

Both cut the channels on a common time window, group them by station
and take positions and orientations from the channel metadata, unless
they are given explicitly.

Conventions
-----------
Coordinates are East, North, Up in metres (``x = East``, ``y = North``,
``z = Up``). When they are computed from geographic positions, the
frame is the local ENU frame of :attr:`ArrayData.frame` (by default
centred on the mean position of the stations); explicit coordinates
are in a local frame chosen by the user. Orientations are unit vectors
of the positive direction of each component in the same frame.
Azimuths are in degrees clockwise from north.

Design choices
--------------
``StationData`` and ``ArrayData`` are immutable: stations, records,
coordinates and orientations cannot be replaced after construction.
The records are stored by reference: the builders create new records
(cut to the common window), so that later changes of the input data do
not affect the snapshot.
"""

from __future__ import annotations

import math
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
from shakelab.signals.base import SNAP_TOLERANCE, Record, StreamCollection
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

# Accepted values of the ``gaps`` argument of the builders.
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

    Raises
    ------
    TypeError
        If a record is not a Record.
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


# ----------------------------------------------------------------------------
# Array data
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class ArrayData:
    """
    Stations of an array on a common time axis.

    Parameters
    ----------
    stations : sequence of StationData
        Stations, in processing order; identifiers must be unique and
        all records must share start time, sampling interval, number of
        samples and units.
    frame : MetricFrame, optional
        Local ENU frame of the coordinates, when they derive from
        geographic positions. None for coordinates in a local frame
        without a geographic reference.
    excluded : sequence of tuple, optional
        ``(sid, reason)`` of the stations left out by the builders.

    Raises
    ------
    TypeError
        If a station is not a StationData or ``frame`` is not a
        MetricFrame.
    ValueError
        If there are no stations, identifiers are repeated or the
        stations are not on the same time axis.
    """

    stations: Sequence[StationData]
    frame: Optional[MetricFrame] = None
    excluded: Sequence[tuple] = ()

    def __post_init__(self) -> None:
        stations = tuple(self.stations)

        if not stations:
            raise ValueError("an array needs at least one station")

        if not all(isinstance(item, StationData) for item in stations):
            raise TypeError("stations must be StationData")

        if self.frame is not None and not isinstance(self.frame, MetricFrame):
            raise TypeError("frame must be a MetricFrame")

        sids = [station.sid for station in stations]

        if len(sids) != len(set(sids)):
            raise ValueError("station identifiers must be unique")

        reference = stations[0]

        for station in stations[1:]:
            if not math.isclose(
                station.delta, reference.delta, rel_tol=_SAMPLING_RTOL
            ):
                raise ValueError(
                    "stations with different sampling intervals "
                    f"({reference.sid}, {station.sid})"
                )

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

        _common_units(
            (rec for station in stations for rec in station.records),
            "array",
        )

        object.__setattr__(self, "stations", stations)
        object.__setattr__(
            self, "excluded", tuple(tuple(item) for item in self.excluded)
        )

    def __len__(self) -> int:
        """Return the number of stations."""
        return len(self.stations)

    def __iter__(self) -> Iterator[StationData]:
        """Iterate over the stations."""
        return iter(self.stations)

    def __getitem__(self, item: Any) -> Any:
        """Return a station by position or identifier."""
        if isinstance(item, str):
            for station in self.stations:
                if station.sid == item:
                    return station

            raise KeyError(item)

        return self.stations[item]

    def __repr__(self) -> str:
        return (
            f"ArrayData({len(self)} stations, {self.nsamp} samples at "
            f"{self.rate:g} Hz from {self.starttime})"
        )

    # -- builders -------------------------------------------------------

    @classmethod
    def from_collection(
        cls,
        collection: StreamCollection,
        select: Any = None,
        starttime: Any = None,
        endtime: Any = None,
        coordinates: Optional[Mapping] = None,
        reference: Any = None,
        orientations: Optional[Mapping] = None,
        coordinate_covariance: CoordinateCovariance = None,
        gaps: str = "drop",
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> ArrayData:
        """
        Build an array from a stream collection.

        Parameters
        ----------
        collection : StreamCollection
            Waveforms of the array, e.g. from
            :func:`shakelab.signals.io.reader` or the buffer of a
            SeedLink client. It is not modified.
        select : str or list of str, optional
            Channel pattern(s) ``NET.STA.LOC.CHA`` (see
            :mod:`shakelab.signals.codes`). Default: every channel.
        starttime, endtime : Date or str, optional
            Common time window, bounds included. Default: the interval
            covered by all the selected channels.
        coordinates : mapping, optional
            Station positions ``{sid: (x, y[, z])}`` in metres, in a
            local frame, for every station (``sid`` is ``NET.STA.LOC``).
            Default: the positions of the records
            (:attr:`Record.location`), in the ENU frame of
            ``reference``.
        reference : MetricFrame, WgsPoint or sequence of float, optional
            Origin of the ENU frame: a frame, or a point ``(lon, lat[,
            elevation])``. Default: the mean position of the stations.
            With ``coordinates`` it is the geographic origin of the
            given frame (optional).
        orientations : mapping, optional
            Orientation matrices ``{sid: matrix}`` overriding those of
            the metadata for some stations, rows in the order Z, N, E
            and then by channel code.
        coordinate_covariance : float, array-like or mapping, optional
            Coordinate uncertainty, in square metres: one value for all
            stations (see :class:`StationData`) or a mapping by station.
        gaps : {"drop", "raise"}, optional
            Stations whose channels do not cover the whole window with
            continuous data: excluded with an
            :class:`ArrayDataWarning` (``"drop"``, default; listed in
            :attr:`excluded`) or an error (``"raise"``).
        snap_tolerance : float, optional
            Times within this fraction of a sampling interval from a
            sample coincide with it.

        Returns
        -------
        ArrayData
            Array with new records cut to the window.

        Raises
        ------
        TypeError
            If ``collection`` is not a StreamCollection.
        ValueError
            If no channel is selected, the channels of a station have
            different band or instrument codes (use ``select``), the
            window is empty, a station has a gap with
            ``gaps="raise"``, no station is left, positions or
            orientations are missing, or the stations are not sampled
            simultaneously.
        """
        if not isinstance(collection, StreamCollection):
            raise TypeError("collection must be a StreamCollection")

        streams = [
            stream for stream in collection
            if stream.records() and _selected(stream.sid, select)
        ]

        if not streams:
            raise ValueError("no channel is selected")

        groups: dict[str, list] = {}

        for stream in streams:
            groups.setdefault(_station_sid(stream.sid), []).append(
                stream.records()
            )

        return _build(
            groups,
            starttime=starttime,
            endtime=endtime,
            coordinates=coordinates,
            reference=reference,
            orientations=orientations,
            covariance=coordinate_covariance,
            gaps=gaps,
            snap_tolerance=snap_tolerance,
        )

    @classmethod
    def from_records(
        cls,
        records: Iterable[Record],
        starttime: Any = None,
        endtime: Any = None,
        coordinates: Optional[Mapping] = None,
        reference: Any = None,
        orientations: Optional[Mapping] = None,
        coordinate_covariance: CoordinateCovariance = None,
        gaps: str = "drop",
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> ArrayData:
        """
        Build an array from records.

        Records of the same channel (same ``sid``) are taken as
        consecutive segments of it. The arguments are those of
        :meth:`from_collection`.

        Parameters
        ----------
        records : iterable of Record
            Records with a source identifier ``NET.STA.LOC.CHA``.
        starttime, endtime : Date or str, optional
            Common time window, bounds included.
        coordinates : mapping, optional
            Station positions ``{sid: (x, y[, z])}`` in metres.
        reference : MetricFrame, WgsPoint or sequence of float, optional
            Origin of the ENU frame.
        orientations : mapping, optional
            Orientation matrices by station.
        coordinate_covariance : float, array-like or mapping, optional
            Coordinate uncertainty, in square metres.
        gaps : {"drop", "raise"}, optional
            Treatment of incomplete stations.
        snap_tolerance : float, optional
            Time tolerance, as a fraction of a sampling interval.

        Returns
        -------
        ArrayData
            Array with new records cut to the window.

        Raises
        ------
        TypeError
            If an element is not a Record.
        ValueError
            As :meth:`from_collection`; also for records without a
            source identifier.
        """
        channels: dict[str, list] = {}

        for record in records:
            if not isinstance(record, Record):
                raise TypeError("records must be Record")

            if not record.sid:
                raise ValueError("records need a source identifier")

            channels.setdefault(record.sid, []).append(record)

        if not channels:
            raise ValueError("no record is given")

        groups: dict[str, list] = {}

        for sid, segments in channels.items():
            groups.setdefault(_station_sid(sid), []).append(segments)

        return _build(
            groups,
            starttime=starttime,
            endtime=endtime,
            coordinates=coordinates,
            reference=reference,
            orientations=orientations,
            covariance=coordinate_covariance,
            gaps=gaps,
            snap_tolerance=snap_tolerance,
        )

    # -- properties -----------------------------------------------------

    @property
    def nsta(self) -> int:
        """Number of stations."""
        return len(self.stations)

    @property
    def sids(self) -> tuple[str, ...]:
        """Station identifiers, in processing order."""
        return tuple(station.sid for station in self.stations)

    @property
    def coordinates(self) -> np.ndarray:
        """Station coordinates ``(n_stations, 3)``, in metres."""
        return np.stack([station.coordinates for station in self.stations])

    @property
    def coordinate_covariances(self) -> Optional[np.ndarray]:
        """
        Coordinate covariances ``(n_stations, 3, 3)``.

        None if no station has one; zeros for the stations without one
        otherwise.
        """
        values = [station.coordinate_covariance for station in self.stations]

        if all(value is None for value in values):
            return None

        return np.stack([
            np.zeros((3, 3)) if value is None else value for value in values
        ])

    @property
    def nsamp(self) -> int:
        """Number of samples of each record."""
        return self.stations[0].nsamp

    @property
    def delta(self) -> float:
        """Sampling interval in seconds."""
        return self.stations[0].delta

    @property
    def rate(self) -> float:
        """Sampling rate in Hz."""
        return self.stations[0].rate

    @property
    def starttime(self) -> Date:
        """Time of the first sample."""
        return self.stations[0].starttime

    @property
    def endtime(self) -> Date:
        """Time of the last sample."""
        return self.stations[0].endtime

    @property
    def duration(self) -> float:
        """Time between first and last sample, in seconds."""
        return (self.nsamp - 1) * self.delta

    @property
    def units(self) -> Optional[str]:
        """Common units of the records, or None if unknown."""
        return _common_units(
            (rec for station in self.stations for rec in station.records),
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
            station.sid for station in self.stations
            if not station.supports(direction, tolerance)
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
            If some stations cannot reconstruct the direction.
        """
        self._require_support((direction,), tolerance)
        return np.stack([
            station.project(direction, tolerance)
            for station in self.stations
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
            If no direction is given or some stations cannot
            reconstruct one.
        """
        directions = tuple(directions)

        if not directions:
            raise ValueError("at least one direction is required")

        self._require_support(directions, tolerance)
        return np.stack([
            station.matrix(directions, tolerance)
            for station in self.stations
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
            Time axis, frame, stations and excluded stations.
        """
        lines = [
            f"ArrayData: {self.nsta} stations",
            f"  time      : {self.starttime} - {self.endtime} "
            f"({self.nsamp} samples at {self.rate:g} Hz)",
            f"  units     : {self.units}",
            "  frame     : "
            + (
                "local, no geographic reference" if self.frame is None
                else "ENU from lon {:.6f}, lat {:.6f}, h {:.1f} m".format(
                    *self.frame.ref_geo
                )
            ),
        ]

        for station in self.stations:
            x, y, z = station.coordinates
            lines.append(
                f"  {station.sid:<16} {station.ncomp} comp  "
                f"x {x:10.2f}  y {y:10.2f}  z {z:8.2f} m"
            )

        for sid, reason in self.excluded:
            lines.append(f"  excluded {sid}: {reason}")

        return "\n".join(lines)

    def _require_support(self, directions: tuple, tolerance: float) -> None:
        """Raise if some station cannot reconstruct a direction."""
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


def _time_window(groups: dict, starttime: Any, endtime: Any) -> tuple:
    """Return the requested window, or the one common to all channels."""
    starts = [
        min(record.starttime for record in segments)
        for channels in groups.values() for segments in channels
    ]
    ends = [
        max(record.endtime for record in segments)
        for channels in groups.values() for segments in channels
    ]

    start = max(starts) if starttime is None else Date(starttime)
    end = min(ends) if endtime is None else Date(endtime)

    if end < start:
        raise ValueError(
            "the selected channels have no common time window"
            if starttime is None and endtime is None
            else f"empty time window ({start} - {end})"
        )

    return start, end


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


def _build(
    groups: dict,
    starttime: Any,
    endtime: Any,
    coordinates: Optional[Mapping],
    reference: Any,
    orientations: Optional[Mapping],
    covariance: Any,
    gaps: str,
    snap_tolerance: float,
) -> ArrayData:
    """Build an ArrayData from channels grouped by station."""
    if gaps not in _GAP_POLICIES:
        raise ValueError("gaps must be 'drop' or 'raise'")

    if coordinates is not None and not isinstance(coordinates, Mapping):
        raise TypeError("coordinates must be a mapping {sid: (x, y[, z])}")

    if orientations is not None and not isinstance(orientations, Mapping):
        raise TypeError("orientations must be a mapping {sid: matrix}")

    for sid, channels in groups.items():
        codes = {
            (segments[0].sid or "").split(".")[-1][:-1]
            for segments in channels
        }

        if len(codes) > 1:
            raise ValueError(
                f"station {sid} has channels of different instruments "
                f"({', '.join(sorted(codes))}): choose them with select"
            )

        if len(channels) > 3:
            raise ValueError(
                f"station {sid} has more than three channels: choose "
                "them with select"
            )

    start, end = _time_window(groups, starttime, endtime)
    stations: dict[str, list] = {}
    excluded = []

    for sid in sorted(groups):
        records = []

        for segments in groups[sid]:
            segments = sorted(segments, key=lambda rec: rec.starttime)
            record = _cut_channel(segments, start, end, snap_tolerance)

            if record is None:
                reason = (
                    f"channel {segments[0].sid} does not cover "
                    f"{start} - {end} without gaps"
                )

                if gaps == "raise":
                    raise ValueError(f"station {sid}: {reason}")

                excluded.append((sid, reason))
                break

            records.append(record)

        else:
            stations[sid] = sorted(records, key=_component_key)

    if excluded:
        warn_external(
            "stations excluded from the array: "
            + "; ".join(f"{sid} ({reason})" for sid, reason in excluded),
            ArrayDataWarning,
        )

    if not stations:
        raise ValueError("no station covers the time window")

    frame = None

    if coordinates is not None:
        missing = sorted(set(stations) - set(coordinates))

        if missing:
            raise ValueError(
                "coordinates are missing for stations " + ", ".join(missing)
            )

        positions = {sid: coordinates[sid] for sid in stations}

        if reference is not None:
            frame = _frame(reference, [])

    else:
        locations = {
            sid: _station_location(sid, records)
            for sid, records in stations.items()
        }
        missing = sorted(sid for sid, point in locations.items()
                         if point is None)

        if missing:
            raise ValueError(
                "the location of stations " + ", ".join(missing)
                + " is unknown: give metadata or coordinates"
            )

        frame = _frame(reference, list(locations.values()))
        positions = {
            sid: geo_to_enu(point, frame) for sid, point in locations.items()
        }

    return ArrayData(
        stations=[
            StationData(
                sid=sid,
                records=records,
                coordinates=positions[sid],
                orientation=_station_value(orientations, sid),
                coordinate_covariance=_station_value(covariance, sid),
            )
            for sid, records in stations.items()
        ],
        frame=frame,
        excluded=excluded,
    )
