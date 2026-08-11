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
Data structures for seismic-array processing.

This module defines the public data model used by array-processing
algorithms. It does not perform file I/O and does not contain any
beamforming-specific implementation.

The global Cartesian reference system is:

    x = East
    y = North
    z = Up

Each station may contain one, two, or three recorded components.
Component orientations are represented by unit vectors expressed in
the global reference system.

``StationData`` and ``ArrayData`` are structurally immutable snapshots:
stations, records, coordinates, and orientations cannot be replaced
after construction. The contained ``Record`` instances are not copied,
however, and their waveform arrays remain mutable. Callers must
therefore avoid modifying records while an array-processing operation
is running.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from shakelab.signals.base import Record


FloatArray = NDArray[np.float64]


EAST = np.array(
    [1.0, 0.0, 0.0],
    dtype=float,
)
NORTH = np.array(
    [0.0, 1.0, 0.0],
    dtype=float,
)
VERTICAL = np.array(
    [0.0, 0.0, 1.0],
    dtype=float,
)

_DEFAULT_ORIENTATION = np.eye(
    3,
    dtype=float,
)

_ORIENTATION_TOLERANCE = 1e-8
_SAMPLING_RTOL = 1e-9


def _readonly_copy(value: ArrayLike) -> FloatArray:
    """Return a non-writeable floating-point array copy."""
    array = np.array(
        value,
        dtype=float,
        copy=True,
    )

    array.setflags(write=False)

    return array


def _as_coordinates(value: ArrayLike) -> FloatArray:
    """
    Convert coordinates to a three-dimensional vector.

    Two-dimensional coordinates are interpreted as East and North,
    with zero elevation.

    Parameters
    ----------
    value
        Sequence containing two or three coordinate values.

    Returns
    -------
    numpy.ndarray
        Read-only coordinate vector with shape ``(3,)``.
    """
    coordinates = np.asarray(
        value,
        dtype=float,
    )

    if coordinates.ndim != 1:
        raise ValueError(
            "Coordinates must be a one-dimensional array."
        )

    if coordinates.size == 2:
        coordinates = np.concatenate(
            (
                coordinates,
                np.array([0.0]),
            )
        )

    elif coordinates.size != 3:
        raise ValueError(
            "Coordinates must contain either two or three values."
        )

    if not np.all(np.isfinite(coordinates)):
        raise ValueError(
            "Coordinates must contain only finite values."
        )

    return _readonly_copy(coordinates)


def _as_orientation(
    value: ArrayLike | None,
    n_components: int,
) -> FloatArray:
    """
    Convert and validate a component-orientation matrix.

    If no matrix is provided, the first ``n_components`` rows of the
    identity matrix are used. This corresponds to the convention:

        component 1 = East
        component 2 = North
        component 3 = Vertical

    Parameters
    ----------
    value
        Orientation matrix or ``None``.
    n_components
        Number of recorded components.

    Returns
    -------
    numpy.ndarray
        Read-only normalized orientation matrix with shape
        ``(n_components, 3)``.
    """
    if not 1 <= n_components <= 3:
        raise ValueError(
            "A station must contain between one and three components."
        )

    if value is None:
        orientation = _DEFAULT_ORIENTATION[:n_components]
        return _readonly_copy(orientation)

    orientation = np.asarray(
        value,
        dtype=float,
    )

    expected_shape = (
        n_components,
        3,
    )

    if orientation.shape != expected_shape:
        raise ValueError(
            "Orientation must have shape "
            f"{expected_shape}, got {orientation.shape}."
        )

    if not np.all(np.isfinite(orientation)):
        raise ValueError(
            "Orientation must contain only finite values."
        )

    norms = np.linalg.norm(
        orientation,
        axis=1,
    )

    if np.any(norms == 0.0):
        raise ValueError(
            "Each orientation row must define a non-zero direction."
        )

    orientation = orientation / norms[:, None]

    singular_values = np.linalg.svd(
        orientation,
        compute_uv=False,
    )

    if singular_values[-1] <= _ORIENTATION_TOLERANCE:
        raise ValueError(
            "Component orientations must be linearly independent "
            "and numerically well conditioned."
        )

    return _readonly_copy(orientation)


def _as_coordinate_covariance(
    value: ArrayLike | None,
) -> FloatArray | None:
    """
    Convert and validate a coordinate covariance matrix.

    A 2 x 2 covariance matrix is interpreted in the horizontal
    East-North plane and expanded to 3 x 3 with zero vertical
    uncertainty.

    Parameters
    ----------
    value
        Coordinate covariance matrix or ``None``.

    Returns
    -------
    numpy.ndarray or None
        Read-only covariance matrix with shape ``(3, 3)``.
    """
    if value is None:
        return None

    covariance = np.asarray(
        value,
        dtype=float,
    )

    if covariance.shape == (2, 2):
        expanded = np.zeros(
            (3, 3),
            dtype=float,
        )

        expanded[:2, :2] = covariance
        covariance = expanded

    elif covariance.shape != (3, 3):
        raise ValueError(
            "Coordinate covariance must have shape (2, 2) "
            "or (3, 3)."
        )

    if not np.all(np.isfinite(covariance)):
        raise ValueError(
            "Coordinate covariance must contain only finite values."
        )

    if not np.allclose(
        covariance,
        covariance.T,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(
            "Coordinate covariance must be symmetric."
        )

    covariance = 0.5 * (
        covariance + covariance.T
    )

    eigenvalues = np.linalg.eigvalsh(
        covariance
    )

    tolerance = (
        np.finfo(float).eps
        * max(
            1.0,
            float(np.max(np.abs(eigenvalues))),
        )
        * covariance.shape[0]
    )

    if np.any(eigenvalues < -tolerance):
        raise ValueError(
            "Coordinate covariance must be positive semidefinite."
        )

    return _readonly_copy(covariance)


def _as_direction(value: ArrayLike) -> FloatArray:
    """
    Convert a direction to a normalized three-dimensional vector.

    Parameters
    ----------
    value
        Direction expressed in East-North-Up coordinates.

    Returns
    -------
    numpy.ndarray
        Normalized direction vector.
    """
    direction = np.asarray(
        value,
        dtype=float,
    )

    if direction.shape != (3,):
        raise ValueError(
            "Direction must contain exactly three values."
        )

    if not np.all(np.isfinite(direction)):
        raise ValueError(
            "Direction must contain only finite values."
        )

    norm = np.linalg.norm(direction)

    if norm == 0.0:
        raise ValueError(
            "Direction must be non-zero."
        )

    return direction / norm


def _record_unit(record: Record):
    """Return the physical unit associated with a record."""
    header = getattr(
        record,
        "head",
        None,
    )

    if header is None:
        return None

    return getattr(
        header,
        "units",
        None,
    )


def _validate_unit_compatibility(
    records: Sequence[Record],
    context: str,
) -> None:
    """
    Validate known physical units.

    Missing units are allowed. If two or more records provide explicit
    units, all known values must be equal.
    """
    known_units = {
        _record_unit(record)
        for record in records
        if _record_unit(record) is not None
    }

    if len(known_units) > 1:
        values = ", ".join(
            sorted(str(unit) for unit in known_units)
        )

        raise ValueError(
            f"{context} contains incompatible physical units: "
            f"{values}."
        )


@dataclass(frozen=True)
class StationData:
    """
    Waveform data and geometry for one array station.

    Parameters
    ----------
    sid
        Unique station identifier.
    records
        Ordered sequence containing one, two, or three waveform
        records.
    coordinates
        Station coordinates in the global East-North-Up reference
        system. Two-dimensional coordinates are accepted and imply
        zero elevation.
    orientation
        Matrix with shape ``(n_components, 3)``. Each row contains
        the unit vector of the corresponding recorded component in
        the East-North-Up reference system.

        If omitted, the default component order is East, North,
        Vertical.
    coordinate_covariance
        Optional 2 x 2 or 3 x 3 covariance matrix describing
        coordinate uncertainty.

    Notes
    -----
    Record order and orientation-row order must match.

    The ``Record`` objects are stored by reference and are not copied.
    They must not be modified while the station snapshot is being
    processed.
    """

    sid: str
    records: Sequence[Record]
    coordinates: ArrayLike
    orientation: ArrayLike | None = None
    coordinate_covariance: ArrayLike | None = None

    def __post_init__(self) -> None:
        sid = str(self.sid).strip()

        if not sid:
            raise ValueError(
                "Station identifier cannot be empty."
            )

        records = tuple(self.records)

        if not 1 <= len(records) <= 3:
            raise ValueError(
                "A station must contain between one and three records."
            )

        if not all(
            isinstance(record, Record)
            for record in records
        ):
            raise TypeError(
                "Station records must be Record instances."
            )

        coordinates = _as_coordinates(
            self.coordinates
        )

        orientation = _as_orientation(
            self.orientation,
            len(records),
        )

        covariance = _as_coordinate_covariance(
            self.coordinate_covariance
        )

        object.__setattr__(
            self,
            "sid",
            sid,
        )
        object.__setattr__(
            self,
            "records",
            records,
        )
        object.__setattr__(
            self,
            "coordinates",
            coordinates,
        )
        object.__setattr__(
            self,
            "orientation",
            orientation,
        )
        object.__setattr__(
            self,
            "coordinate_covariance",
            covariance,
        )

        self._validate_records()

    def _validate_records(self) -> None:
        """Check internal consistency among component records."""
        reference = self.records[0]

        if reference.nsamp == 0:
            raise ValueError(
                f"Station {self.sid!r} contains an empty record."
            )

        if reference.delta is None:
            raise ValueError(
                f"Station {self.sid!r} has no sampling interval."
            )

        if reference.starttime is None:
            raise ValueError(
                f"Station {self.sid!r} has no start time."
            )

        reference_delta = float(reference.delta)
        reference_npts = int(reference.nsamp)
        reference_starttime = reference.starttime

        if (
            not np.isfinite(reference_delta)
            or reference_delta <= 0.0
        ):
            raise ValueError(
                f"Station {self.sid!r} has an invalid "
                "sampling interval."
            )

        for index, record in enumerate(self.records):
            data = np.asarray(record.data)

            if data.ndim != 1:
                raise ValueError(
                    f"Record {index} of station {self.sid!r} "
                    "is not one-dimensional."
                )

            if data.size == 0:
                raise ValueError(
                    f"Record {index} of station {self.sid!r} "
                    "is empty."
                )

            if not np.all(np.isfinite(data)):
                raise ValueError(
                    f"Record {index} of station {self.sid!r} "
                    "contains non-finite samples."
                )

            if record.delta is None:
                raise ValueError(
                    f"Record {index} of station {self.sid!r} "
                    "has no sampling interval."
                )

            if record.starttime is None:
                raise ValueError(
                    f"Record {index} of station {self.sid!r} "
                    "has no start time."
                )

            if record.nsamp != reference_npts:
                raise ValueError(
                    f"Records of station {self.sid!r} do not have "
                    "the same number of samples."
                )

            if not np.isclose(
                float(record.delta),
                reference_delta,
                rtol=_SAMPLING_RTOL,
                atol=0.0,
            ):
                raise ValueError(
                    f"Records of station {self.sid!r} do not have "
                    "the same sampling interval."
                )

            if record.starttime != reference_starttime:
                raise ValueError(
                    f"Records of station {self.sid!r} do not have "
                    "the same start time."
                )

        _validate_unit_compatibility(
            self.records,
            f"Station {self.sid!r}",
        )

    def __len__(self) -> int:
        """Return the number of recorded components."""
        return len(self.records)

    def __iter__(self) -> Iterator[Record]:
        """Iterate over component records."""
        return iter(self.records)

    @property
    def n_components(self) -> int:
        """Number of recorded components."""
        return len(self.records)

    @property
    def npts(self) -> int:
        """Number of samples in each component."""
        return self.records[0].nsamp

    @property
    def delta(self) -> float:
        """Sampling interval in seconds."""
        return float(self.records[0].delta)

    @property
    def rate(self) -> float:
        """Sampling rate in samples per second."""
        return float(self.records[0].rate)

    @property
    def starttime(self):
        """Start time shared by all station records."""
        return self.records[0].starttime

    @property
    def units(self):
        """
        Return the common known physical unit.

        ``None`` is returned when all component units are undefined.
        """
        for record in self.records:
            unit = _record_unit(record)

            if unit is not None:
                return unit

        return None

    @property
    def data(self) -> FloatArray:
        """
        Return recorded components as a matrix.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_components, n_samples)``. Rows
            follow the order of ``records`` and ``orientation``.
        """
        return np.stack(
            [
                np.asarray(
                    record.data,
                    dtype=float,
                )
                for record in self.records
            ],
            axis=0,
        )

    @property
    def orientation_rank(self) -> int:
        """Rank of the component-orientation matrix."""
        return int(
            np.linalg.matrix_rank(
                self.orientation
            )
        )

    def projection_weights(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> FloatArray:
        """
        Return weights required to reconstruct a target direction.

        Parameters
        ----------
        direction
            Target direction in East-North-Up coordinates.
        tolerance
            Maximum accepted reconstruction residual.

        Returns
        -------
        numpy.ndarray
            One weight for each recorded component.

        Raises
        ------
        ValueError
            If the requested direction is not observable from the
            available components.
        """
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError(
                "Projection tolerance must be finite and non-negative."
            )

        target = _as_direction(direction)

        weights, _, _, _ = np.linalg.lstsq(
            self.orientation.T,
            target,
            rcond=None,
        )

        reconstructed = (
            self.orientation.T @ weights
        )

        residual = np.linalg.norm(
            reconstructed - target
        )

        if residual > tolerance:
            raise ValueError(
                f"Direction {target.tolist()} cannot be reconstructed "
                f"from the available components of station "
                f"{self.sid!r}."
            )

        return weights

    def supports(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> bool:
        """Return whether a target direction can be reconstructed."""
        try:
            self.projection_weights(
                direction,
                tolerance=tolerance,
            )

        except ValueError:
            return False

        return True

    def project(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> FloatArray:
        """
        Project the recorded motion onto a target direction.

        Returns
        -------
        numpy.ndarray
            Projected waveform with shape ``(n_samples,)``.
        """
        weights = self.projection_weights(
            direction,
            tolerance=tolerance,
        )

        return weights @ self.data

    def matrix(
        self,
        directions: Iterable[ArrayLike],
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> FloatArray:
        """
        Project the station onto multiple target directions.

        Parameters
        ----------
        directions
            Iterable of target directions.
        tolerance
            Maximum accepted reconstruction residual.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_directions, n_samples)``.
        """
        direction_list = tuple(directions)

        if not direction_list:
            raise ValueError(
                "At least one projection direction is required."
            )

        weights = np.stack(
            [
                self.projection_weights(
                    direction,
                    tolerance=tolerance,
                )
                for direction in direction_list
            ],
            axis=0,
        )

        return weights @ self.data

    def east(self) -> FloatArray:
        """Return the East component."""
        return self.project(EAST)

    def north(self) -> FloatArray:
        """Return the North component."""
        return self.project(NORTH)

    def vertical(self) -> FloatArray:
        """Return the vertical component."""
        return self.project(VERTICAL)

    def horizontal(self) -> FloatArray:
        """
        Return East and North components.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(2, n_samples)``.
        """
        return self.matrix(
            (
                EAST,
                NORTH,
            )
        )

    def canonical(self) -> FloatArray:
        """
        Return East, North, and vertical components.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(3, n_samples)``.
        """
        return self.matrix(
            (
                EAST,
                NORTH,
                VERTICAL,
            )
        )

    def radial(self, azimuth: float) -> FloatArray:
        """
        Return the radial component for a propagation azimuth.

        Azimuth is expressed in degrees clockwise from North.
        """
        azimuth = float(azimuth)

        if not np.isfinite(azimuth):
            raise ValueError(
                "Azimuth must be finite."
            )

        angle = np.deg2rad(azimuth)

        direction = np.array(
            [
                np.sin(angle),
                np.cos(angle),
                0.0,
            ]
        )

        return self.project(direction)

    def transverse(self, azimuth: float) -> FloatArray:
        """
        Return the transverse component for a propagation azimuth.

        Azimuth is expressed in degrees clockwise from North.
        """
        azimuth = float(azimuth)

        if not np.isfinite(azimuth):
            raise ValueError(
                "Azimuth must be finite."
            )

        angle = np.deg2rad(azimuth)

        direction = np.array(
            [
                np.cos(angle),
                -np.sin(angle),
                0.0,
            ]
        )

        return self.project(direction)


@dataclass(frozen=True)
class ArrayData:
    """
    Immutable collection of stations prepared for array processing.

    Parameters
    ----------
    stations
        Ordered sequence of :class:`StationData` instances.

    Notes
    -----
    Station order is preserved and defines the row order of all
    matrices returned by this class.

    ``ArrayData`` validates the common sampling and time axis during
    construction. Those checks are not repeated by individual
    properties.
    """

    stations: Sequence[StationData] = field(
        default_factory=tuple
    )

    def __post_init__(self) -> None:
        stations = tuple(self.stations)

        if not stations:
            raise ValueError(
                "ArrayData must contain at least one station."
            )

        if not all(
            isinstance(station, StationData)
            for station in stations
        ):
            raise TypeError(
                "All ArrayData elements must be StationData instances."
            )

        object.__setattr__(
            self,
            "stations",
            stations,
        )

        self._validate_station_ids()
        self._validate_common_time_axis()
        self._validate_units()

    def __len__(self) -> int:
        """Return the number of stations."""
        return len(self.stations)

    def __iter__(self) -> Iterator[StationData]:
        """Iterate over stations in processing order."""
        return iter(self.stations)

    def __getitem__(
        self,
        index: int | slice,
    ):
        """Return one or more stations by position."""
        return self.stations[index]

    def _validate_station_ids(self) -> None:
        """Check that station identifiers are unique."""
        station_ids = self.station_ids

        if len(station_ids) != len(set(station_ids)):
            raise ValueError(
                "ArrayData contains duplicate station identifiers."
            )

    def _validate_common_time_axis(self) -> None:
        """Check sampling and temporal consistency among stations."""
        reference = self.stations[0]

        for station in self.stations[1:]:
            if station.npts != reference.npts:
                raise ValueError(
                    "Array stations do not have the same number "
                    "of samples."
                )

            if not np.isclose(
                station.delta,
                reference.delta,
                rtol=_SAMPLING_RTOL,
                atol=0.0,
            ):
                raise ValueError(
                    "Array stations do not have the same sampling "
                    "interval."
                )

            if station.starttime != reference.starttime:
                raise ValueError(
                    "Array stations do not have the same start time."
                )

    def _validate_units(self) -> None:
        """Check known physical units across all station records."""
        records = tuple(
            record
            for station in self.stations
            for record in station.records
        )

        _validate_unit_compatibility(
            records,
            "ArrayData",
        )

    @property
    def n_stations(self) -> int:
        """Number of stations in the array."""
        return len(self.stations)

    @property
    def station_ids(self) -> tuple[str, ...]:
        """Ordered station identifiers."""
        return tuple(
            station.sid
            for station in self.stations
        )

    @property
    def coordinates(self) -> FloatArray:
        """
        Return station coordinates.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_stations, 3)``.
        """
        return np.stack(
            [
                station.coordinates
                for station in self.stations
            ],
            axis=0,
        )

    @property
    def coordinate_covariances(self) -> FloatArray | None:
        """
        Return coordinate covariance matrices.

        Returns
        -------
        numpy.ndarray or None
            Array with shape ``(n_stations, 3, 3)``. If no station
            provides coordinate uncertainty, ``None`` is returned.
            Missing individual covariance matrices are represented by
            zeros when at least one station provides them.
        """
        covariances = tuple(
            station.coordinate_covariance
            for station in self.stations
        )

        if all(
            covariance is None
            for covariance in covariances
        ):
            return None

        return np.stack(
            [
                (
                    np.zeros(
                        (3, 3),
                        dtype=float,
                    )
                    if covariance is None
                    else covariance
                )
                for covariance in covariances
            ],
            axis=0,
        )

    @property
    def npts(self) -> int:
        """Number of samples in each component record."""
        return self.stations[0].npts

    @property
    def delta(self) -> float:
        """Sampling interval shared by all stations."""
        return self.stations[0].delta

    @property
    def rate(self) -> float:
        """Sampling rate shared by all stations."""
        return self.stations[0].rate

    @property
    def starttime(self):
        """Start time shared by all stations."""
        return self.stations[0].starttime

    @property
    def units(self):
        """Return the common known physical unit."""
        for station in self.stations:
            if station.units is not None:
                return station.units

        return None

    def supports(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> bool:
        """Return whether all stations support a target direction."""
        return all(
            station.supports(
                direction,
                tolerance=tolerance,
            )
            for station in self.stations
        )

    def unsupported_stations(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> tuple[str, ...]:
        """
        Return identifiers of stations that cannot reconstruct a
        target direction.
        """
        return tuple(
            station.sid
            for station in self.stations
            if not station.supports(
                direction,
                tolerance=tolerance,
            )
        )

    def project(
        self,
        direction: ArrayLike,
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> FloatArray:
        """
        Project all stations onto one target direction.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_stations, n_samples)``.
        """
        unsupported = self.unsupported_stations(
            direction,
            tolerance=tolerance,
        )

        if unsupported:
            station_list = ", ".join(unsupported)

            raise ValueError(
                "The requested direction is not supported by "
                f"stations: {station_list}."
            )

        return np.stack(
            [
                station.project(
                    direction,
                    tolerance=tolerance,
                )
                for station in self.stations
            ],
            axis=0,
        )

    def matrix(
        self,
        directions: Iterable[ArrayLike],
        *,
        tolerance: float = _ORIENTATION_TOLERANCE,
    ) -> FloatArray:
        """
        Project all stations onto multiple target directions.

        Parameters
        ----------
        directions
            Iterable of target directions in East-North-Up
            coordinates.
        tolerance
            Maximum accepted reconstruction residual.

        Returns
        -------
        numpy.ndarray
            Array with shape
            ``(n_stations, n_directions, n_samples)``.
        """
        direction_list = tuple(directions)

        if not direction_list:
            raise ValueError(
                "At least one projection direction is required."
            )

        return np.stack(
            [
                station.matrix(
                    direction_list,
                    tolerance=tolerance,
                )
                for station in self.stations
            ],
            axis=0,
        )

    def vertical(self) -> FloatArray:
        """
        Return vertical data for all stations.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_stations, n_samples)``.
        """
        return self.project(VERTICAL)

    def horizontal(self) -> FloatArray:
        """
        Return East and North data for all stations.

        Returns
        -------
        numpy.ndarray
            Array with shape
            ``(n_stations, 2, n_samples)``.
        """
        return self.matrix(
            (
                EAST,
                NORTH,
            )
        )

    def three_component(self) -> FloatArray:
        """
        Return East, North, and vertical data for all stations.

        Returns
        -------
        numpy.ndarray
            Array with shape
            ``(n_stations, 3, n_samples)``.
        """
        return self.matrix(
            (
                EAST,
                NORTH,
                VERTICAL,
            )
        )