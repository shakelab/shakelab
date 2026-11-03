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
Waveform metadata.

This module provides :class:`Header`, the metadata of a waveform
record or spectrum: sampling, start time, identifiers, units and the
description of the recording channel. It depends only on the time
primitives and on :mod:`shakelab.signals.metadata`, so that both
time-domain and frequency-domain objects can use it.

Design choices
--------------
- The source identifier and the description of the channel
  (coordinates, orientation, equipment, instrument response) are a
  :class:`~shakelab.signals.metadata.ChannelEpoch`, ``head.meta``,
  always present (empty by default).
- ``head.meta`` is shared, not copied: the records of the same epoch
  refer to one object, and :meth:`Header.copy` keeps the reference.
  Changes made through the header and record attributes (``sid``,
  ``location``, ...) replace the epoch with an updated copy, so they
  never affect other records. Changing the shared object directly
  (``head.meta.location = ...``) affects all the records referring to
  it.
- ``head.extra`` holds free information about one record (e.g. a
  quality flag); free information about the channel belongs to
  ``head.meta.extra``.

Conventions
-----------
- Times are :class:`~shakelab.libutils.timeN.Date` instances (UTC);
  sampling interval ``delta`` in seconds and sampling rate ``rate`` in
  Hz, with ``rate = 1 / delta``.
- :meth:`Header.info` prints a summary for interactive inspection.
"""

from __future__ import annotations

from copy import deepcopy
from numbers import Real
from typing import Any, Optional

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.signals.metadata import ChannelEpoch

__all__ = ["Header"]


# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------


class Header:
    """
    Metadata associated with a waveform record or spectrum.

    A Header stores the sampling interval, sampling rate, start time,
    event identifier and physical units of a :class:`Record`, and the
    description of its channel (:attr:`meta`), which includes the
    source identifier.

    The sampling interval and sampling rate are kept synchronized:
    ``rate = 1 / delta``.

    Parameters
    ----------
    delta
        Sampling interval in seconds. It must be finite and strictly
        positive.
    time
        Start time of the associated record, in any format accepted by
        :class:`~shakelab.libutils.timeN.Date` (e.g. an ISO 8601
        string). The header stores its own copy.
    sid
        Source identifier ``NET.STA.LOC.CHA`` (FDSN source code).
    eid
        Event identifier.
    units
        Physical units of the waveform samples.
    parent
        Object owning the header, normally a :class:`Record`.
    meta
        Channel epoch (:class:`~shakelab.signals.metadata.ChannelEpoch`)
        of the record, stored by reference. Default: an empty epoch
        with identifier ``sid``.

    Attributes
    ----------
    extra : dict
        Free additional information about this record only.

    Raises
    ------
    TypeError
        If ``meta`` is not a ChannelEpoch.
    ValueError
        If ``sid`` differs from the identifier of ``meta``.

    Notes
    -----
    The parent reference is intentionally excluded from :meth:`copy`
    unless explicitly provided. This prevents copying a Header from
    recursively duplicating its owning Record.
    """

    def __init__(
        self,
        delta: Optional[Real] = None,
        time: Any = None,
        sid: Optional[str] = None,
        eid: Optional[str] = None,
        units: Optional[str] = None,
        parent: Any = None,
        meta: Optional[ChannelEpoch] = None,
    ) -> None:
        self._delta: Optional[float] = None
        self._rate: Optional[float] = None
        self._time: Optional[Date] = None
        self._parent = parent

        if meta is None:
            meta = ChannelEpoch(sid)
        elif not isinstance(meta, ChannelEpoch):
            raise TypeError("meta must be a ChannelEpoch")
        elif sid is not None and sid != meta.sid:
            raise ValueError(
                f"source identifier {sid!r} differs from the metadata "
                f"identifier {meta.sid!r}"
            )
        self._meta = meta

        self.eid = eid
        self.units = units
        self.extra: dict = {}

        if delta is not None:
            self.delta = delta

        if time is not None:
            self.time = time

    def __repr__(self) -> str:
        return (
            f"Header(sid={self.sid!r}, time={self.time}, "
            f"delta={self.delta!r})"
        )

    def __str__(self) -> str:
        meta = self.meta
        lines = [
            f"sampling interval: {self.delta}",
            f"sampling rate: {self.rate}",
            f"starting time: {self.time}",
            f"source identifier: {self.sid}",
            f"event identifier: {self.eid}",
            f"units: {self.units}",
            f"location: {meta.location}",
            f"orientation: {meta.orientation}",
            "response: "
            + ("declared" if meta.response is not None else "none"),
        ]

        return "\n".join(lines)

    @property
    def meta(self) -> ChannelEpoch:
        """
        Description of the channel (shared channel epoch).

        Assigning None replaces it with an empty epoch keeping the
        source identifier.
        """
        return self._meta

    @meta.setter
    def meta(self, value: Optional[ChannelEpoch]) -> None:
        if value is None:
            value = ChannelEpoch(self.sid)
        elif not isinstance(value, ChannelEpoch):
            raise TypeError("meta must be a ChannelEpoch")
        self._meta = value

    @property
    def sid(self) -> Optional[str]:
        """
        Source identifier, stored in :attr:`meta`.

        Assigning a new identifier replaces :attr:`meta` with an
        updated copy (see :meth:`ChannelEpoch.replace`).
        """
        return self._meta.sid

    @sid.setter
    def sid(self, value: Optional[str]) -> None:
        if value != self._meta.sid:
            self._meta = self._meta.replace(sid=value)

    @staticmethod
    def _validate_sampling_value(value: Any, name: str) -> float:
        """
        Validate a sampling interval or sampling rate.

        Parameters
        ----------
        value
            Value to validate.
        name
            Parameter name used in error messages.

        Returns
        -------
        float
            Validated value converted to float.

        Raises
        ------
        TypeError
            If the value is not a real number.
        ValueError
            If the value is not finite or is not strictly positive.
        """
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"{name} must be a real number")

        value = float(value)

        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")

        if value <= 0.0:
            raise ValueError(f"{name} must be strictly positive")

        return value

    @property
    def delta(self) -> Optional[float]:
        """Sampling interval in seconds."""
        return self._delta

    @delta.setter
    def delta(self, value: Optional[Real]) -> None:
        if value is None:
            self._delta = None
            self._rate = None
            return

        value = self._validate_sampling_value(value, "delta")

        self._delta = value
        self._rate = 1.0 / value

    @property
    def rate(self) -> Optional[float]:
        """Sampling rate in samples per second."""
        return self._rate

    @rate.setter
    def rate(self, value: Optional[Real]) -> None:
        if value is None:
            self._rate = None
            self._delta = None
            return

        value = self._validate_sampling_value(value, "rate")

        self._rate = value
        self._delta = 1.0 / value

    @property
    def time(self) -> Optional[Date]:
        """
        Start time of the associated record.

        Any value accepted by :class:`~shakelab.libutils.timeN.Date` can
        be assigned; the header always stores its own copy, so later
        changes to the assigned object do not affect the record.
        """
        return self._time

    @time.setter
    def time(self, value: Any) -> None:
        if value is None:
            self._time = None
            return

        if isinstance(value, Date):
            self._time = value.copy()
            return

        self._time = Date(value)

    @property
    def nsamp(self) -> Optional[int]:
        """
        Number of samples in the parent object, when available.

        Returns
        -------
        int or None
            Number of samples exposed by the parent object, or ``None``
            when no compatible parent is assigned.
        """
        if self._parent is None:
            return None

        return getattr(self._parent, "nsamp", None)

    def info(self) -> None:
        """Print a human-readable description of the header."""
        print(self)

    def copy(self, parent: Any = None) -> Header:
        """
        Return a copy of the header.

        The channel epoch (:attr:`meta`) is shared with the original
        (see the module notes); the other attributes are copied.

        Parameters
        ----------
        parent
            Parent to associate with the copied header. By default, the
            copied header has no parent.

        Returns
        -------
        Header
            Header copy.

        Notes
        -----
        The parent reference is not copied automatically because it
        normally points back to the owning Record. Copying that reference
        with ``deepcopy`` could duplicate the complete object graph or
        leave the copied Header associated with the wrong Record.
        """
        header = Header(
            delta=self.delta,
            time=self.time,
            eid=self.eid,
            units=self.units,
            parent=parent,
            meta=self._meta,
        )

        header.extra = deepcopy(self.extra)

        return header
