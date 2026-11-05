# ****************************************************************************
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
# ****************************************************************************
"""
Channel and station metadata.

This module describes the recording channels: where the sensor is
(coordinates and burial depth), how it is oriented, which equipment
records, at which sampling rate, and with which instrument response.
The description is valid in a time interval (an epoch); a collection
holds the epochs of many channels and selects the one applicable to a
source identifier and a time.

Classes
-------
:class:`Orientation`
    Azimuth and dip of a component.
:class:`Equipment`
    Description of a sensor or a datalogger.
:class:`StationEpoch`
    Station (``NET.STA``) valid in a time interval.
:class:`ChannelEpoch`
    Channel (``NET.STA.LOC.CHA``) valid in a time interval, with its
    station and instrument response.
:class:`MetadataCollection`
    Channel and station epochs, read from StationXML, RESP or JSON.

Design choices
--------------
- The model is flat: a collection holds channel epochs, each with a
  reference to its station epoch, rather than a tree of networks,
  stations and channels. Stations are also available on their own
  (:meth:`MetadataCollection.stations`), for maps, arrays and location.
- The instrument response of an epoch is a
  :class:`shakelab.signals.response.StageSet`: the response module is
  the model of the instrument, this module adds identifiers and epochs.
- Formats provide what they contain: StationXML the complete
  description, RESP only the response (no coordinates, orientation or
  equipment).
- The metadata are read as written; inconsistencies are reported by
  :meth:`MetadataCollection.validate`, never corrected.
- Collections store copies (``append``), as the other ShakeLab
  collections.

Conventions
-----------
- Coordinates follow ``docs/conventions.md``: ``location`` is a
  :class:`~shakelab.libutils.geodeticN.WgsPoint` (longitude, latitude,
  elevation in m, positive upward). The elevation is taken as written
  in the metadata (StationXML refers it to the local ground surface,
  usually mean sea level; using it as an ellipsoidal height is an
  approximation).
- ``depth`` is the burial depth of the sensor below ``location``, in
  km, positive downward (StationXML ``Depth`` is in m). The elevation
  of the sensor is ``location.elevation - 1000 * depth``.
- Azimuth in degrees clockwise from north, in ``[0, 360)``; dip in
  degrees, positive downward from the horizontal, in ``[-90, 90]``
  (SEED and StationXML: a vertical component positive upward has dip
  -90).
- Epochs include their start time and exclude their end time; a
  missing bound is open.
"""

from __future__ import annotations

import json as _json
import math
import os as _os
import warnings
from copy import copy as _copy
from copy import deepcopy
from dataclasses import dataclass, fields
from numbers import Real
from typing import Any, Iterator, Optional

import numpy as np

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.libutils.timeN import Date
from shakelab.signals import response as rspm
from shakelab.signals.codes import match_code, parse_sid_pattern

__all__ = [
    "MetadataWarning",
    "Orientation",
    "Equipment",
    "StationEpoch",
    "ChannelEpoch",
    "MetadataCollection",
    "select_response",
    "select_sensitivity",
]


class MetadataWarning(UserWarning):
    """Warning about incomplete or inconsistent channel metadata."""


# ----------------------------------------------------------------------------
# Value objects
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Orientation:
    """
    Orientation of a component.

    Parameters
    ----------
    azimuth
        Degrees clockwise from north; reduced to ``[0, 360)``.
    dip
        Degrees below the horizontal, in ``[-90, 90]``: 0 for a
        horizontal component, -90 for a vertical component positive
        upward.

    Raises
    ------
    TypeError
        If a value is not a real number.
    ValueError
        If a value is not finite or the dip is out of range.
    """

    azimuth: float
    dip: float

    def __post_init__(self) -> None:
        azimuth = _finite(self.azimuth, "azimuth") % 360.0
        dip = _finite(self.dip, "dip")
        if not -90.0 <= dip <= 90.0:
            raise ValueError("the dip must be in [-90, 90] degrees")
        object.__setattr__(self, "azimuth", azimuth)
        object.__setattr__(self, "dip", dip)

    @property
    def vector(self) -> np.ndarray:
        """
        Unit vector of the positive direction, as (east, north, up).
        """
        azimuth = math.radians(self.azimuth)
        dip = math.radians(self.dip)
        return np.array([
            math.cos(dip) * math.sin(azimuth),
            math.cos(dip) * math.cos(azimuth),
            -math.sin(dip),
        ])

    def to_dict(self) -> dict:
        """Return the attributes as a dictionary."""
        return {"azimuth": self.azimuth, "dip": self.dip}


@dataclass(frozen=True)
class Equipment:
    """
    Description of a sensor or a datalogger.

    Attributes follow the StationXML ``Equipment`` type; all are
    optional strings.
    """

    type: Optional[str] = None
    description: Optional[str] = None
    manufacturer: Optional[str] = None
    model: Optional[str] = None
    serial_number: Optional[str] = None

    def __str__(self) -> str:
        name = " ".join(
            item for item in (self.manufacturer, self.model) if item
        )
        return name or self.description or self.type or ""

    def to_dict(self) -> dict:
        """Return the attributes as a dictionary."""
        return {item.name: getattr(self, item.name) for item in fields(self)}


# ----------------------------------------------------------------------------
# Epochs
# ----------------------------------------------------------------------------


class _Epoch:
    """Validity interval shared by station and channel epochs."""

    _starttime: Optional[Date]
    _endtime: Optional[Date]

    @property
    def starttime(self) -> Optional[Date]:
        """Start of validity (included); None if open."""
        return self._starttime

    @starttime.setter
    def starttime(self, value: Any) -> None:
        self._starttime = _date(value)

    @property
    def endtime(self) -> Optional[Date]:
        """End of validity (excluded); None if open."""
        return self._endtime

    @endtime.setter
    def endtime(self, value: Any) -> None:
        self._endtime = _date(value)

    def match(self, time: Any) -> bool:
        """
        Return whether a time falls within the epoch.

        Parameters
        ----------
        time
            Any value accepted by :class:`Date`.
        """
        time = time if isinstance(time, Date) else Date(time)
        return (
            (self._starttime is None or time >= self._starttime)
            and (self._endtime is None or time < self._endtime)
        )

    def overlaps(self, other: _Epoch) -> bool:
        """Return whether two epochs share an interval of time."""
        starts_before_end = (
            self._starttime is None or other._endtime is None
            or self._starttime < other._endtime
        )
        ends_after_start = (
            self._endtime is None or other._starttime is None
            or self._endtime > other._starttime
        )
        return starts_before_end and ends_after_start

    def _interval(self) -> str:
        start = str(self._starttime) if self._starttime else "open"
        end = str(self._endtime) if self._endtime else "open"
        return f"{start} -> {end}"


class StationEpoch(_Epoch):
    """
    Station valid in a time interval.

    Parameters
    ----------
    code
        Station code ``NET.STA``.
    starttime, endtime
        Validity interval (any value accepted by :class:`Date`).
    location
        :class:`WgsPoint` or ``(longitude, latitude[, elevation])``.
    site
        Name of the site.
    description
        Free description.
    """

    def __init__(
        self,
        code: str,
        starttime: Any = None,
        endtime: Any = None,
        location: Any = None,
        site: Optional[str] = None,
        description: Optional[str] = None,
    ) -> None:
        self.code = _code(code, 2, "station code")
        self.starttime = starttime
        self.endtime = endtime
        self.location = location
        self.site = site
        self.description = description

    @property
    def location(self) -> Optional[WgsPoint]:
        """Coordinates of the station."""
        return self._location

    @location.setter
    def location(self, value: Any) -> None:
        self._location = _point(value)

    @property
    def network(self) -> str:
        """Network code."""
        return self.code.split(".")[0]

    @property
    def station(self) -> str:
        """Station code without the network."""
        return self.code.split(".")[1]

    def __repr__(self) -> str:
        return f"StationEpoch({self.code!r}, {self._interval()})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StationEpoch):
            return NotImplemented
        return self.to_dict() == other.to_dict()

    __hash__ = None

    def to_dict(self) -> dict:
        """Return the attributes as a dictionary of plain values."""
        return {
            "code": self.code,
            "starttime": _date_str(self.starttime),
            "endtime": _date_str(self.endtime),
            "location": _point_list(self.location),
            "site": self.site,
            "description": self.description,
        }

    @classmethod
    def from_dict(cls, data: dict) -> StationEpoch:
        """Build a station epoch from :meth:`to_dict` output."""
        return cls(**data)

    def copy(self) -> StationEpoch:
        """Return a copy."""
        return deepcopy(self)


class ChannelEpoch(_Epoch):
    """
    Channel valid in a time interval.

    Parameters
    ----------
    sid
        Source identifier ``NET.STA.LOC.CHA``; None only for the empty
        epoch of a record without identifier (see :attr:`empty`).
    starttime, endtime
        Validity interval (any value accepted by :class:`Date`).
    location
        Coordinates of the sensor: :class:`WgsPoint` or
        ``(longitude, latitude[, elevation])``.
    depth
        Burial depth in km, positive downward.
    orientation
        :class:`Orientation` or ``(azimuth, dip)`` in degrees.
    sample_rate
        Declared sampling rate in Hz.
    sensor, datalogger
        :class:`Equipment` descriptions.
    response
        Instrument response (:class:`StageSet`).
    station
        Station epoch the channel belongs to.

    Attributes
    ----------
    extra : dict
        Arbitrary additional metadata of the channel.

    Notes
    -----
    A record refers to the epoch of its channel through
    ``record.head.meta``. The records of the same epoch share one
    object; the record methods never modify it but replace it with an
    updated copy (:meth:`replace`).
    """

    def __init__(
        self,
        sid: Optional[str] = None,
        starttime: Any = None,
        endtime: Any = None,
        location: Any = None,
        depth: Optional[Real] = None,
        orientation: Any = None,
        sample_rate: Optional[Real] = None,
        sensor: Optional[Equipment] = None,
        datalogger: Optional[Equipment] = None,
        response: Optional[rspm.StageSet] = None,
        station: Optional[StationEpoch] = None,
    ) -> None:
        self.sid = sid
        self.starttime = starttime
        self.endtime = endtime
        self.location = location
        self.depth = depth
        self.orientation = orientation
        self.sample_rate = sample_rate
        self.sensor = sensor
        self.datalogger = datalogger
        self.response = response
        self.station = station
        self.extra: dict = {}

    # -- validated attributes -------------------------------------------

    @property
    def sid(self) -> Optional[str]:
        """Source identifier ``NET.STA.LOC.CHA``."""
        return self._sid

    @sid.setter
    def sid(self, value: Optional[str]) -> None:
        self._sid = (
            None if value is None else _code(value, 4, "source identifier")
        )

    @property
    def location(self) -> Optional[WgsPoint]:
        """Coordinates of the sensor (elevation of the ground)."""
        return self._location

    @location.setter
    def location(self, value: Any) -> None:
        self._location = _point(value)

    @property
    def depth(self) -> Optional[float]:
        """Burial depth in km, positive downward."""
        return self._depth

    @depth.setter
    def depth(self, value: Optional[Real]) -> None:
        self._depth = None if value is None else _finite(value, "depth")

    @property
    def orientation(self) -> Optional[Orientation]:
        """Orientation of the component."""
        return self._orientation

    @orientation.setter
    def orientation(self, value: Any) -> None:
        if value is None or isinstance(value, Orientation):
            self._orientation = value
            return
        try:
            azimuth, dip = value
        except (TypeError, ValueError):
            raise TypeError(
                "orientation must be an Orientation or (azimuth, dip)"
            ) from None
        self._orientation = Orientation(azimuth, dip)

    @property
    def sample_rate(self) -> Optional[float]:
        """Declared sampling rate in Hz."""
        return self._sample_rate

    @sample_rate.setter
    def sample_rate(self, value: Optional[Real]) -> None:
        if value is None:
            self._sample_rate = None
            return
        value = _finite(value, "sample_rate")
        if value < 0:
            raise ValueError("sample_rate must not be negative")
        self._sample_rate = value

    @property
    def response(self) -> Optional[rspm.StageSet]:
        """Instrument response."""
        return self._response

    @response.setter
    def response(self, value: Optional[rspm.StageSet]) -> None:
        if value is not None and not isinstance(value, rspm.StageSet):
            raise TypeError("response must be a StageSet")
        self._response = value

    # -- codes ----------------------------------------------------------

    @property
    def network(self) -> Optional[str]:
        """Network code."""
        return self._part(0)

    @property
    def station_code(self) -> Optional[str]:
        """Station code ``NET.STA``."""
        if self.sid is None:
            return None
        return ".".join(self.sid.split(".")[:2])

    @property
    def location_code(self) -> Optional[str]:
        """Location code (may be empty)."""
        return self._part(2)

    @property
    def channel(self) -> Optional[str]:
        """Channel code."""
        return self._part(3)

    def _part(self, index: int) -> Optional[str]:
        return None if self.sid is None else self.sid.split(".")[index]

    # -- other ----------------------------------------------------------

    @property
    def empty(self) -> bool:
        """
        Whether the epoch describes nothing beyond the identifier.

        An empty epoch (no validity interval, coordinates, orientation,
        sampling rate, equipment, response, station or extra
        information) is the default metadata of a record.
        """
        return not self.extra and all(
            getattr(self, name) is None for name in _CONTENT_FIELDS
        )

    def replace(self, **changes: Any) -> ChannelEpoch:
        """
        Return a copy with some attributes changed.

        Parameters
        ----------
        **changes
            New values of constructor attributes (e.g.
            ``location=(13.0, 46.0)``), validated as in the
            constructor.

        Returns
        -------
        ChannelEpoch
            Shallow copy: the unchanged response and station are shared
            with the original, ``extra`` is copied.

        Raises
        ------
        TypeError
            If an attribute is unknown.
        """
        unknown = set(changes) - _CHANNEL_FIELDS
        if unknown:
            raise TypeError(
                f"unknown channel attributes: {', '.join(sorted(unknown))}"
            )
        epoch = _copy(self)
        epoch.extra = dict(self.extra)
        for key, value in changes.items():
            setattr(epoch, key, value)
        return epoch

    def __repr__(self) -> str:
        return f"ChannelEpoch({self.sid!r}, {self._interval()})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ChannelEpoch):
            return NotImplemented
        return (
            self.to_dict() == other.to_dict()
            and self.station == other.station
        )

    __hash__ = None

    def conflicts(self, other: ChannelEpoch) -> list:
        """
        Attributes defined in both epochs with different values.

        Attributes missing in either epoch are not conflicts, so the
        partial metadata of a data file (e.g. coordinates read from SAC)
        agree with a complete StationXML epoch of the same channel.
        Coordinates are compared within 1e-5 degrees (about 1 m, the
        precision of single-precision formats) and 0.1 m of elevation
        (an elevation of 0, the default of coordinates given without
        it, is not compared),
        depths within 1 mm, orientations within 1e-3 degrees and
        sampling rates within a relative 1e-4.

        Parameters
        ----------
        other : ChannelEpoch
            Epoch to compare.

        Returns
        -------
        list of str
            Names of the conflicting attributes, in definition order.

        Raises
        ------
        TypeError
            If ``other`` is not a ChannelEpoch.
        """
        if not isinstance(other, ChannelEpoch):
            raise TypeError("other must be a ChannelEpoch")

        def close(a: Any, b: Any, tol: float) -> bool:
            return abs(a - b) <= tol

        def same_location(a: WgsPoint, b: WgsPoint) -> bool:
            # Elevation 0 also stands for "not given" (WgsPoint default).
            elevations = (a.elevation or 0.0, b.elevation or 0.0)
            return (
                close(a.longitude, b.longitude, 1e-5)
                and close(a.latitude, b.latitude, 1e-5)
                and (0.0 in elevations
                     or close(a.elevation, b.elevation, 0.1))
            )

        checks = {
            "sid": lambda a, b: a == b,
            "starttime": lambda a, b: a == b,
            "endtime": lambda a, b: a == b,
            "location": same_location,
            "depth": lambda a, b: close(a, b, 1e-6),
            "orientation": lambda a, b: (
                close((a.azimuth - b.azimuth + 180.0) % 360.0, 180.0, 1e-3)
                and close(a.dip, b.dip, 1e-3)
            ),
            "sample_rate": lambda a, b: close(a, b, 1e-4 * abs(b)),
            "sensor": lambda a, b: a == b,
            "datalogger": lambda a, b: a == b,
            "response": lambda a, b: a.to_dict() == b.to_dict(),
            "station": lambda a, b: a == b,
        }
        names = []
        for name, same in checks.items():
            a, b = getattr(self, name), getattr(other, name)
            if a is not None and b is not None and not same(a, b):
                names.append(name)
        if any(
            self.extra[key] != other.extra[key]
            for key in self.extra.keys() & other.extra.keys()
        ):
            names.append("extra")
        return names

    def info(self) -> None:
        """Print a summary of the channel epoch."""
        print(f"{self.sid}  {self._interval()}")
        if self.location is not None:
            loc = self.location
            print(
                f"  location: lon {loc.longitude} lat {loc.latitude} "
                f"elev {loc.elevation} m, depth {self.depth} km"
            )
        if self.orientation is not None:
            print(
                f"  orientation: azimuth {self.orientation.azimuth}, "
                f"dip {self.orientation.dip}"
            )
        if self.sample_rate is not None:
            print(f"  sample rate: {self.sample_rate} Hz")
        for label, item in (
            ("sensor", self.sensor), ("datalogger", self.datalogger)
        ):
            if item is not None:
                print(f"  {label}: {item}")
        if self.response is not None:
            sens = self.response.sensitivity
            text = (
                f"{sens.value} {sens.output_units} per {sens.input_units}"
                f" at {sens.frequency} Hz"
                if sens is not None else "no sensitivity"
            )
            print(f"  response: {len(self.response.stage)} stages, {text}")

    def to_dict(self) -> dict:
        """
        Return the attributes as a dictionary of plain values.

        The station is not included (see
        :meth:`MetadataCollection.write`).
        """
        return {
            "sid": self.sid,
            "starttime": _date_str(self.starttime),
            "endtime": _date_str(self.endtime),
            "location": _point_list(self.location),
            "depth": self.depth,
            "orientation": (
                self.orientation.to_dict() if self.orientation else None
            ),
            "sample_rate": self.sample_rate,
            "sensor": self.sensor.to_dict() if self.sensor else None,
            "datalogger": (
                self.datalogger.to_dict() if self.datalogger else None
            ),
            "response": (
                self.response.to_dict() if self.response else None
            ),
            "extra": deepcopy(self.extra),
        }

    @classmethod
    def from_dict(
        cls, data: dict, station: Optional[StationEpoch] = None
    ) -> ChannelEpoch:
        """Build a channel epoch from :meth:`to_dict` output."""
        data = dict(data)
        extra = data.pop("extra", None) or {}
        orientation = data.pop("orientation", None)
        sensor = data.pop("sensor", None)
        datalogger = data.pop("datalogger", None)
        resp = data.pop("response", None)
        epoch = cls(
            **data,
            orientation=(
                Orientation(**orientation) if orientation else None
            ),
            sensor=Equipment(**sensor) if sensor else None,
            datalogger=Equipment(**datalogger) if datalogger else None,
            response=rspm.StageSet.from_dict(resp) if resp else None,
            station=station,
        )
        epoch.extra = extra
        return epoch

    def copy(self) -> ChannelEpoch:
        """Return a deep copy (including the station epoch)."""
        return deepcopy(self)


_CONTENT_FIELDS = (
    "_starttime", "_endtime", "_location", "_depth", "_orientation",
    "_sample_rate", "sensor", "datalogger", "_response", "station",
)
_CHANNEL_FIELDS = {
    "sid", "starttime", "endtime", "location", "depth", "orientation",
    "sample_rate", "sensor", "datalogger", "response", "station", "extra",
}


# ----------------------------------------------------------------------------
# Collection
# ----------------------------------------------------------------------------


class MetadataCollection:
    """
    Channel and station epochs indexed by identifier.

    Parameters
    ----------
    source
        Optional file to read (see :meth:`read`).
    ftype
        File type; inferred from the file name if None.

    Notes
    -----
    Iteration and ``len`` refer to channel epochs. ``collection[sid]``
    returns the epochs of a channel; :meth:`get` the epoch valid at a
    time; :meth:`select` a sub-collection.
    """

    def __init__(self, source: Any = None, ftype: Optional[str] = None):
        self._channels: dict[str, list[ChannelEpoch]] = {}
        self._stations: list[StationEpoch] = []
        if source is not None:
            self.read(source, ftype=ftype)

    def __len__(self) -> int:
        return sum(len(epochs) for epochs in self._channels.values())

    def __iter__(self) -> Iterator[ChannelEpoch]:
        for epochs in self._channels.values():
            yield from epochs

    def __contains__(self, sid: Any) -> bool:
        return sid in self._channels

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, str):
            if key not in self._channels:
                raise KeyError(f"source identifier not found: {key}")
            return list(self._channels[key])
        return list(self)[key]

    def __repr__(self) -> str:
        return (
            f"MetadataCollection(stations={len(self._stations)}, "
            f"channels={len(self._channels)}, epochs={len(self)})"
        )

    @property
    def sid(self) -> list:
        """Source identifiers of the channels."""
        return list(self._channels)

    # -- building -------------------------------------------------------

    def append(self, item: Any) -> None:
        """
        Append a copy of a channel or station epoch.

        The station of a channel epoch is shared with an equal station
        epoch already in the collection, or added. A channel epoch equal
        to one already present is ignored; an epoch that overlaps a
        different epoch of the same channel is added with a
        :class:`MetadataWarning`.

        Raises
        ------
        TypeError
            If the item is not a ChannelEpoch or a StationEpoch.
        ValueError
            If a channel epoch has no source identifier.
        """
        if isinstance(item, StationEpoch):
            self._add_station(item)
            return
        if not isinstance(item, ChannelEpoch):
            raise TypeError("only channel and station epochs can be added")
        if item.sid is None:
            raise ValueError("a channel epoch needs a source identifier")

        epochs = self._channels.setdefault(item.sid, [])
        for epoch in epochs:
            if epoch == item:
                return
            if epoch.overlaps(item):
                warnings.warn(
                    f"{item.sid}: epoch {item._interval()} overlaps "
                    f"{epoch._interval()}",
                    MetadataWarning,
                    stacklevel=2,
                )

        station = item.station
        new = deepcopy(item, {id(station): None})
        new.station = (
            self._add_station(station) if station is not None else None
        )
        epochs.append(new)
        epochs.sort(key=_start_key)

    def extend(self, other: Any) -> None:
        """
        Append the content of another collection, or of an iterable of
        channel and station epochs (see :meth:`append`).
        """
        if isinstance(other, MetadataCollection):
            for station in other._stations:
                self._add_station(station)
        for item in other:
            self.append(item)

    def _add_station(self, station: StationEpoch) -> StationEpoch:
        """Return the stored station equal to ``station``, adding it."""
        if not isinstance(station, StationEpoch):
            raise TypeError("station must be a StationEpoch")
        for stored in self._stations:
            if stored == station:
                return stored
        stored = station.copy()
        self._stations.append(stored)
        return stored

    # -- selection ------------------------------------------------------

    def get(self, sid: str, time: Any = None) -> ChannelEpoch:
        """
        Return the epoch of a channel valid at a time.

        Parameters
        ----------
        sid
            Source identifier.
        time
            Any value accepted by :class:`Date`. It may be omitted if
            the channel has only one epoch.

        Returns
        -------
        ChannelEpoch
            The stored epoch (not a copy). If several overlapping
            epochs are valid, the first is returned with a
            :class:`MetadataWarning`.

        Raises
        ------
        KeyError
            If the channel is not found or no epoch is valid at the
            time.
        ValueError
            If the time is omitted and the channel has several epochs.
        """
        epochs = self[sid]
        if time is None:
            if len(epochs) > 1:
                raise ValueError(
                    f"{sid} has {len(epochs)} epochs: a time is required"
                )
            return epochs[0]

        valid = [epoch for epoch in epochs if epoch.match(time)]
        if not valid:
            raise KeyError(f"no epoch of {sid} is valid at {time}")
        if len(valid) > 1:
            warnings.warn(
                f"{sid}: {len(valid)} epochs valid at {time}, the first "
                "is used",
                MetadataWarning,
                stacklevel=2,
            )
        return valid[0]

    def select(self, sid: str = "*", time: Any = None) -> MetadataCollection:
        """
        Return a sub-collection.

        Parameters
        ----------
        sid
            Source identifier pattern, matched code by code (see
            :mod:`shakelab.signals.codes`), e.g. ``"IV.*..HN?"``,
            ``"IV.{MODE,CAE}"`` or ``"IV"`` (every channel of IV).
        time
            If given, keep only the epochs valid at that time.

        Returns
        -------
        MetadataCollection
            Copies of the matching epochs, with their stations.
        """
        patterns = parse_sid_pattern(sid)
        result = MetadataCollection()
        for epoch in self:
            if not any(pattern.matches(epoch.sid) for pattern in patterns):
                continue
            if time is not None and not epoch.match(time):
                continue
            result.append(epoch)
        return result

    def stations(self, code: str = "*", time: Any = None) -> list:
        """
        Return the station epochs.

        Parameters
        ----------
        code
            Station code pattern ``NET.STA``, matched code by code (see
            :mod:`shakelab.signals.codes`), e.g. ``"IV.*"`` or ``"IV"``.
        time
            If given, keep only the epochs valid at that time.

        Returns
        -------
        list of StationEpoch
            The stored station epochs (not copies).
        """
        return [
            station for station in self._stations
            if match_code(station.code, code)
            and (time is None or station.match(time))
        ]

    # -- checks and output ----------------------------------------------

    def validate(self, rtol: float = 1e-3) -> list:
        """
        Check the consistency of the metadata.

        The checks are: overlapping epochs of the same channel; the
        response of each epoch (see :meth:`StageSet.validate`); the
        declared sampling rate against that of the response.

        Returns
        -------
        list of str
            Problems found, each prefixed by channel and epoch start.
        """
        issues = []
        for sid, epochs in self._channels.items():
            for first, second in zip(epochs, epochs[1:]):
                if first.overlaps(second):
                    issues.append(
                        f"{sid}: epochs {first._interval()} and "
                        f"{second._interval()} overlap"
                    )
            for epoch in epochs:
                prefix = f"{sid} [{epoch.starttime}]"
                resp = epoch.response
                if resp is None:
                    continue
                issues.extend(
                    f"{prefix}: {issue}"
                    for issue in resp.validate(rtol=rtol)
                )
                if (
                    epoch.sample_rate and resp.sample_rate
                    and not math.isclose(
                        epoch.sample_rate, resp.sample_rate, rel_tol=rtol
                    )
                ):
                    issues.append(
                        f"{prefix}: sampling rate {epoch.sample_rate} Hz "
                        f"differs from that of the response "
                        f"{resp.sample_rate} Hz"
                    )
        return issues

    def info(self) -> None:
        """Print a summary of the stations and channels."""
        print(
            f"{len(self._stations)} stations, {len(self._channels)} "
            f"channels, {len(self)} epochs"
        )
        for station in self._stations:
            loc = station.location
            where = (
                f"lon {loc.longitude} lat {loc.latitude} "
                f"elev {loc.elevation} m" if loc else "no coordinates"
            )
            site = f" ({station.site})" if station.site else ""
            print(f"{station.code}{site}  {station._interval()}  {where}")
        for epoch in self:
            orientation = (
                f"az {epoch.orientation.azimuth:g} "
                f"dip {epoch.orientation.dip:g}"
                if epoch.orientation else "orientation unknown"
            )
            rate = (
                f"{epoch.sample_rate:g} Hz" if epoch.sample_rate else "? Hz"
            )
            response = "response" if epoch.response else "no response"
            print(
                f"  {epoch.sid}  {epoch._interval()}  {rate}  "
                f"{orientation}  {response}"
            )

    def copy(self) -> MetadataCollection:
        """Return a deep copy."""
        return deepcopy(self)

    def read(
        self,
        source: Any,
        ftype: Optional[str] = None,
        append: bool = False,
    ) -> None:
        """
        Read metadata from a StationXML, RESP or JSON file.

        Parameters
        ----------
        source
            File path; for StationXML and RESP also text, bytes or an
            open file (see
            :func:`shakelab.signals.stationxml.read_stationxml` and
            :func:`shakelab.signals.respfile.read_resp`).
        ftype
            ``"xml"`` (also ``"sxml"``), ``"resp"`` or ``"json"``. If
            None, it is inferred from the file name: ``"resp"`` for
            names starting with ``RESP.`` (as written by rdseed) or
            ending with ``.resp``, otherwise the extension.
        append
            If True, add to the current content; otherwise replace it.

        Raises
        ------
        ValueError
            If the file type is unknown.
        """
        if ftype is None:
            ftype = _file_type(source)
        ftype = ftype.lower()

        if ftype in ("xml", "sxml"):
            # Local imports: the readers build metadata objects, so
            # they sit above this module.
            from shakelab.signals.stationxml import read_stationxml

            new = read_stationxml(source)

        elif ftype == "resp":
            from shakelab.signals.respfile import read_resp

            new = read_resp(source)

        elif ftype == "json":
            with open(source, "r") as fobj:
                new = _collection_from_json(_json.load(fobj))

        else:
            raise ValueError(f"Unknown file type: {ftype}")

        if not append:
            self._channels = {}
            self._stations = []
        self.extend(new)

    def write(self, file: Any, ftype: str = "json") -> None:
        """
        Write the metadata to a JSON file.

        Stations are written once; channel epochs refer to them by
        position.

        Raises
        ------
        NotImplementedError
            For file types other than JSON.
        """
        if ftype != "json":
            raise NotImplementedError("Only the json format is supported")

        index = {id(station): i for i, station in enumerate(self._stations)}
        channels = []
        for epoch in self:
            entry = epoch.to_dict()
            entry["station"] = (
                index[id(epoch.station)] if epoch.station else None
            )
            channels.append(entry)

        content = {
            "stations": [station.to_dict() for station in self._stations],
            "channels": channels,
        }
        with open(file, "w") as fobj:
            _json.dump(content, fobj, indent=2)


def _collection_from_json(content: dict) -> MetadataCollection:
    """Build a collection from the content of a JSON file."""
    collection = MetadataCollection()
    stations = [
        StationEpoch.from_dict(item) for item in content.get("stations", [])
    ]
    for station in stations:
        collection.append(station)
    for entry in content.get("channels", []):
        entry = dict(entry)
        position = entry.pop("station", None)
        station = stations[position] if position is not None else None
        collection.append(ChannelEpoch.from_dict(entry, station=station))
    return collection


def _file_type(path: Any) -> str:
    """Return the metadata file type inferred from a file name."""
    name = _os.path.basename(str(path))
    if name.upper().startswith("RESP."):
        return "resp"
    _, ext = _os.path.splitext(name)
    return ext.lstrip(".").lower()


# ----------------------------------------------------------------------------
# Response selection
# ----------------------------------------------------------------------------


def select_response(
    resp: Any,
    sid: Optional[str] = None,
    time: Any = None,
) -> Any:
    """
    Return the response applicable to a channel at a time.

    Parameters
    ----------
    resp
        Stage, stage set, channel epoch or metadata collection.
    sid
        Source identifier of the data.
    time
        Start time of the data.

    Returns
    -------
    StageResponse or StageSet
        Object providing ``transfer_function``. Stages and stage sets
        are returned unchanged.

    Raises
    ------
    TypeError
        If the response type is unsupported.
    ValueError
        If the identifier needed for the selection is missing, if a
        channel epoch belongs to another channel or another time, or if
        it has no response.
    KeyError
        If the channel or the epoch are not available.
    """
    if isinstance(resp, (rspm.StageResponse, rspm.StageSet)):
        return resp

    if isinstance(resp, MetadataCollection):
        if sid is None:
            raise ValueError(
                "source identifier is required to select a response"
            )
        return select_response(resp.get(sid, time), sid=sid, time=time)

    if isinstance(resp, ChannelEpoch):
        if sid is not None and resp.sid != sid:
            raise ValueError(
                f"the metadata of {resp.sid} cannot be applied to {sid}"
            )
        if time is not None and not resp.match(time):
            raise ValueError(
                f"the epoch {resp._interval()} of {resp.sid} is not "
                f"valid at {time}"
            )
        if resp.response is None:
            raise ValueError(
                f"{resp.sid or 'record'}: no instrument response"
            )
        return resp.response

    raise TypeError(f"unsupported response object: {type(resp).__name__}")


def select_sensitivity(
    resp: Any,
    sid: Optional[str] = None,
    time: Any = None,
) -> rspm.Sensitivity:
    """
    Return the overall sensitivity applicable to a channel at a time.

    Parameters
    ----------
    resp
        Gain stage, stage set, channel epoch or metadata collection.
    sid, time
        Source identifier and start time of the data (see
        :func:`select_response`).

    Raises
    ------
    ValueError
        If the stage set declares no sensitivity.
    TypeError
        If the selected response is not a stage set or a gain stage.
    """
    selected = select_response(resp, sid=sid, time=time)

    if isinstance(selected, rspm.StageSet):
        sens = selected.sensitivity
        if sens is None or not sens.value:
            raise ValueError("the response declares no overall sensitivity")
        return sens

    if isinstance(selected, rspm.StageGain):
        return rspm.Sensitivity(selected.sensitivity, selected.frequency)

    raise TypeError(
        f"no sensitivity available from {type(selected).__name__}"
    )


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _finite(value: Any, name: str) -> float:
    """Return a finite float, or raise."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _date(value: Any) -> Optional[Date]:
    """Return a copy of a date, None for an open bound."""
    if value is None or value == "":
        return None
    if isinstance(value, Date):
        return value.copy()
    return Date(value)


def _date_str(value: Optional[Date]) -> Optional[str]:
    return None if value is None else str(value)


def _point(value: Any) -> Optional[WgsPoint]:
    """Return a WgsPoint from a point or a coordinate sequence."""
    if value is None or isinstance(value, WgsPoint):
        return value
    try:
        coordinates = tuple(value)
    except TypeError:
        raise TypeError(
            "location must be a WgsPoint or (longitude, latitude"
            "[, elevation])"
        ) from None
    if len(coordinates) not in (2, 3):
        raise ValueError(
            "location must contain longitude, latitude and optionally "
            "elevation"
        )
    return WgsPoint(*coordinates)


def _point_list(value: Optional[WgsPoint]) -> Optional[list]:
    if value is None:
        return None
    return [value.longitude, value.latitude, value.elevation]


def _code(value: Any, parts: int, name: str) -> str:
    """Check a dotted code with the given number of parts."""
    if not isinstance(value, str) or value.count(".") != parts - 1:
        raise ValueError(
            f"{name} must have {parts} dot-separated parts: {value!r}"
        )
    return value


def _start_key(epoch: _Epoch) -> tuple:
    """Sort key of epochs by start time (open start first)."""
    start = epoch.starttime
    return (start is not None, start)
