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
Time-limited circular waveform streams (experimental).

This module provides :class:`CircularStream`, a specialized
:class:`~shakelab.signals.base.Stream` that retains only the most
recent waveform data within a configurable time window, for real-time
acquisition (e.g. from a SeedLink server).

Design choices
--------------
- Records are stored as independent blocks in a deque. Contiguous
  records are not concatenated automatically, avoiding repeated
  allocation and copying of increasingly large NumPy arrays. A
  continuous segment is obtained on demand with
  ``extract(..., merge=True)``.
- Appended records are copied, so the caller's records are never
  modified when the buffer is trimmed.
- Signal processing (filtering, tapering, resampling, response
  correction, ...) is not applied to the buffer: applied block by block
  it would introduce edge effects at every block boundary. Selections
  (slices and :meth:`extract`) return standard streams, which can be
  processed.

Status
------
Experimental: the API, the internal representation and the performance
may change. Importing the module emits a :class:`FutureWarning`.
"""

from __future__ import annotations

import warnings
from collections import deque
from copy import deepcopy
from math import ceil, isclose, isfinite
from numbers import Real
from typing import Any, Iterable, Optional

from shakelab.signals.base import SNAP_TOLERANCE, Record, Stream

__all__ = ["CircularStream"]

warnings.warn(
    "shakelab.signals.circular is experimental; its API may change",
    FutureWarning,
    stacklevel=2,
)


# ----------------------------------------------------------------------------
# CircularStream
# ----------------------------------------------------------------------------


class CircularStream(Stream):
    """
    Time-limited stream implemented as a deque of waveform blocks.

    A CircularStream retains the samples that fall within
    ``max_duration`` seconds of the latest appended sample.

    Records must be appended in chronological order. Contiguous records
    remain separate blocks, preventing repeated concatenation of large
    waveform arrays.

    Parameters
    ----------
    sid
        Source identifier associated with the stream.
    max_duration
        Maximum temporal span retained by the buffer, in seconds.
        Default is 3600 seconds.

    Notes
    -----
    The retained interval includes both temporal boundaries. Therefore,
    a continuous stream sampled every ``delta`` seconds may contain
    ``floor(max_duration / delta) + 1`` samples between its first and
    final retained sample.

    Gaps between records are preserved. Overlapping records are rejected
    unless ``enforce=True`` is passed to :meth:`append`. In that case,
    existing samples at or after the new record start time are discarded,
    giving precedence to the new record.
    """

    def __init__(
        self,
        sid: Optional[str] = None,
        max_duration: Real = 3600.0,
    ) -> None:
        self.max_duration = self._validate_max_duration(max_duration)

        super().__init__(sid=sid)

    def __repr__(self) -> str:
        return (
            f"CircularStream(sid={self.sid!r}, "
            f"max_duration={self.max_duration!r}, "
            f"records={len(self.record)})"
        )

    def _new_record_container(
        self,
        records: Optional[Iterable[Record]] = None,
    ) -> deque:
        """Create a deque-based record container."""
        if records is None:
            records = ()

        return deque(records)

    @staticmethod
    def _validate_max_duration(value: Any) -> float:
        """
        Validate the maximum buffer duration.

        Raises
        ------
        TypeError
            If the value is not a real number.
        ValueError
            If the value is not finite or strictly positive.
        """
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("max_duration must be a real number")

        value = float(value)

        if not isfinite(value):
            raise ValueError("max_duration must be finite")

        if value <= 0.0:
            raise ValueError("max_duration must be strictly positive")

        return value

    @staticmethod
    def _validate_tolerance(value: Any) -> float:
        """Validate a non-negative snapping tolerance."""
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("snap_tolerance must be a real number")

        value = float(value)

        if not isfinite(value) or value < 0.0:
            raise ValueError(
                "snap_tolerance must be finite and non-negative"
            )

        return value

    @staticmethod
    def _time_difference(time_1: Any, time_0: Any) -> float:
        """Return the difference between two absolute times in seconds."""
        return float(time_1 - time_0)

    def copy(self) -> CircularStream:
        """
        Return an independent copy of the circular stream.

        Returns
        -------
        CircularStream
            Circular stream with independent records and the same
            duration limit.
        """
        stream = self.__class__(
            sid=deepcopy(self.sid),
            max_duration=self.max_duration,
        )

        stream._replace_records(record.copy() for record in self.record)

        return stream

    def append(
        self,
        record: Record,
        enforce: bool = False,
        fillvalue: Any = 0.0,
        merge_tolerance: Optional[Real] = None,
        snap_tolerance: Real = SNAP_TOLERANCE,
    ) -> None:
        """
        Append a waveform block to the circular buffer.

        Parameters
        ----------
        record
            Record to append. The buffer stores a copy.
        enforce
            If ``False``, overlapping records are rejected. If ``True``,
            existing samples at or after the new record start time are
            discarded before the new record is appended. Default is
            ``False``.
        fillvalue
            Accepted for API compatibility with :meth:`Stream.append`.
            Gaps are preserved and are never filled by CircularStream.
        merge_tolerance
            Accepted for API compatibility with :meth:`Stream.append`:
            blocks are never joined in the buffer.
        snap_tolerance
            Times closer than this to a sample, as a fraction of one
            sampling interval, coincide with it when the buffer is
            trimmed. Default is
            :data:`~shakelab.signals.base.SNAP_TOLERANCE`.

        Raises
        ------
        TypeError
            If the input record or arguments have invalid types.
        ValueError
            If source identifiers are incompatible, temporal metadata is
            missing, records are appended out of chronological order, or
            an overlap is found without ``enforce=True``.
        """
        if not isinstance(record, Record):
            raise TypeError("record must be a Record instance")

        if not isinstance(enforce, bool):
            raise TypeError("enforce must be a boolean")

        snap_tolerance = self._validate_tolerance(snap_tolerance)

        if self.sid is None:
            self.sid = record.sid

        if (
            self.sid is not None
            and record.sid is not None
            and self.sid != record.sid
        ):
            raise ValueError(
                "record source identifier does not match the stream"
            )

        if record.nsamp == 0:
            return

        record._require_time_axis()

        record = record.copy()

        if record.sid is None:
            record.sid = self.sid

        if self.record:
            previous = self.record[-1]

            start_offset = self._time_difference(
                record.starttime,
                previous.starttime,
            )

            if start_offset < 0.0:
                raise ValueError(
                    "records must be appended in chronological order"
                )

            if record.starttime <= previous.endtime:
                if not enforce:
                    raise ValueError(
                        "record overlaps the current circular buffer; "
                        "use enforce=True to replace existing samples"
                    )

                self._discard_overlapping_tail(
                    record.starttime,
                    snap_tolerance=snap_tolerance,
                )

        self.record.append(record)
        self._enforce_duration_limit(snap_tolerance=snap_tolerance)

    def _discard_overlapping_tail(
        self,
        starttime: Any,
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> None:
        """
        Remove buffered samples at or after a given start time.

        Complete overlapping records are removed from the right side of
        the deque. If only the final record is partially overlapped, it
        is truncated to retain samples strictly preceding ``starttime``.
        """
        while self.record and self.record[-1].starttime >= starttime:
            self.record.pop()

        if not self.record:
            return

        last = self.record[-1]

        if last.endtime < starttime:
            return

        offset = self._time_difference(starttime, last.starttime)
        position = offset / last.delta
        nearest = round(position)

        if isclose(position, nearest, rel_tol=0.0, abs_tol=snap_tolerance):
            keep_samples = int(nearest)
        else:
            keep_samples = int(ceil(position))

        if keep_samples <= 0:
            self.record.pop()
            return

        if keep_samples >= last.nsamp:
            return

        last.data = last.data[:keep_samples]

    def _enforce_duration_limit(
        self,
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> None:
        """
        Discard samples older than the configured time window.

        Complete expired records are removed with ``popleft()``. At most
        one record, the oldest retained block, is sliced.
        """
        if not self.record:
            return

        latest_time = self.record[-1].endtime
        cutoff = latest_time - self.max_duration

        while self.record and self.record[0].endtime < cutoff:
            self.record.popleft()

        if not self.record:
            return

        first = self.record[0]

        if first.starttime >= cutoff:
            return

        offset = self._time_difference(cutoff, first.starttime)
        position = offset / first.delta
        nearest = round(position)

        if isclose(position, nearest, rel_tol=0.0, abs_tol=snap_tolerance):
            first_index = int(nearest)
        else:
            first_index = int(ceil(position))

        if first_index <= 0:
            return

        if first_index >= first.nsamp:
            self.record.popleft()
            return

        first.data = first.data[first_index:]
        first.time = first.time + first_index * first.delta

    def sort(self, reverse: bool = False) -> None:
        """
        Validate or restore chronological record ordering.

        Parameters
        ----------
        reverse
            Reverse ordering is not supported for a circular stream.

        Raises
        ------
        TypeError
            If ``reverse`` is not a boolean.
        ValueError
            If reverse ordering is requested.

        Notes
        -----
        Circular streams must remain ordered from earliest to latest so
        that append and duration enforcement remain efficient.
        """
        if not isinstance(reverse, bool):
            raise TypeError("reverse must be a boolean")

        if reverse:
            raise ValueError("a CircularStream cannot be reverse-sorted")

        super().sort(reverse=False)

    # ------------------------------------------------------------------------
    # Processing is not applied block by block
    # ------------------------------------------------------------------------

    def _processing_not_supported(self, *args: Any, **kwargs: Any) -> None:
        """
        Refuse block-by-block processing of the buffer.

        Raises
        ------
        NotImplementedError
            Always: processing each block separately would introduce
            edge effects at every block boundary.
        """
        raise NotImplementedError(
            "processing is not applied to a CircularStream block by "
            "block; extract a continuous segment first, e.g. "
            "stream.extract(starttime, endtime, merge=True), and "
            "process the resulting Record"
        )

    remove_mean = _processing_not_supported
    filter = _processing_not_supported
    resample = _processing_not_supported
    taper = _processing_not_supported
    shift = _processing_not_supported
    integrate = _processing_not_supported
    differentiate = _processing_not_supported
    convolve_response = _processing_not_supported
    deconvolve_response = _processing_not_supported
