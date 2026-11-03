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
Core waveform classes and processing utilities.

This module provides the fundamental data structures used throughout
ShakeLab for uniformly sampled time-series analysis.

The main classes are:

- :class:`Header`, storing waveform metadata (defined in
  :mod:`shakelab.signals.header` and re-exported here).
- :class:`Record`, representing a continuous waveform segment together
  with its metadata.
- :class:`Stream`, representing a collection of temporally ordered
  records belonging to the same source.
- :class:`StreamCollection`, grouping multiple streams. With
  ``max_duration`` it holds time-limited circular streams
  (:class:`shakelab.signals.circular.CircularStream`, experimental).

The module also provides basic signal-processing utilities including
filtering, resampling, tapering, Fourier-domain operations, temporal
windowing, convolution, correlation, and response correction.

Waveform objects are designed to be lightweight, largely independent of
external libraries, and suitable for scientific processing pipelines.

Conventions
-----------
- Times are :class:`~shakelab.libutils.timeN.Date` instances (UTC);
  sampling interval ``delta`` in seconds and sampling rate ``rate`` in
  Hz, with ``rate = 1 / delta``.
- Locations are :class:`~shakelab.libutils.geodeticN.WgsPoint`
  instances, given as ``(longitude, latitude[, elevation])``.
- Time tolerances are fractions of the sampling interval:
  :data:`MERGE_TOLERANCE` when joining records,
  :data:`SNAP_TOLERANCE` when matching times to samples.
- The ``info()`` methods print a summary for interactive inspection;
  they are the only functions of the module that print.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from numbers import Integral, Number, Real
from typing import TYPE_CHECKING, Any, Iterator, Optional

import numpy as np
from scipy import integrate, signal
from scipy.signal import windows

from shakelab.signals import fourier
from shakelab.signals import metadata
from shakelab.signals import response
from shakelab.signals import io
from shakelab.signals.header import Header
from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    # Annotations only: spectrum sits above this module.
    from shakelab.signals.spectrum import Spectrum

__all__ = [
    "MERGE_TOLERANCE",
    "SNAP_TOLERANCE",
    "Header",
    "Record",
    "Stream",
    "StreamCollection",
    "tukey_window",
]

MERGE_TOLERANCE = 0.5
"""
Default tolerance for joining contiguous records, as a fraction of the
sampling interval: two records are contiguous when the start of the
second one is within half a sample of the expected time (as in
libmseed). It absorbs the limited time resolution of data formats
(e.g. 0.1 ms in MiniSEED 2).
"""

SNAP_TOLERANCE = 1e-6
"""
Default tolerance for snapping a time to a sample of the sampling grid,
as a fraction of the sampling interval: a time this close to a sample
coincides with it. It only absorbs floating-point rounding.
"""


# ----------------------------------------------------------------------------
# Record
# ----------------------------------------------------------------------------


class Record:
    """
    Continuous waveform segment.

    A Record contains a one-dimensional array of samples and the
    associated metadata stored in a :class:`Header`.

    Parameters
    ----------
    data : array-like, optional
        Waveform samples. Scalar values are converted to one-sample
        arrays. Data must be one-dimensional.
    delta : float, optional
        Sampling interval in seconds.
    time : Date, str or other Date input, optional
        Start time of the record, in any format accepted by
        :class:`~shakelab.libutils.timeN.Date`.
    location : WgsPoint or sequence of float, optional
        Geographic location ``(longitude, latitude[, elevation])``.
    sid : str, optional
        Source identifier.
    eid : str, optional
        Event identifier.
    units : str, optional
        Physical units of the waveform samples.

    Notes
    -----
    Numerical operators return independent Record objects and never
    modify their operands. Signal-processing methods generally operate
    in place unless explicitly documented otherwise.
    """

    def __init__(
        self,
        data: Any = None,
        delta: Optional[Real] = None,
        time: Any = None,
        location: Any = None,
        sid: Optional[str] = None,
        eid: Optional[str] = None,
        units: Optional[str] = None,
    ) -> None:
        self.head = Header(
            delta=delta,
            time=time,
            location=location,
            sid=sid,
            eid=eid,
            units=units,
            parent=self,
        )

        self._data = np.array([], dtype=float)

        if data is not None:
            self.data = data

    def __len__(self) -> int:
        """
        Return the number of samples.
        """
        return self.nsamp

    def __getitem__(self, item: Any) -> Any:
        """
        Return one or more waveform samples.
        """
        return self.data[item]

    def __str__(self) -> str:
        return (
            "starttime {0}, duration {1}s, sampling rate {2} Hz, "
            "{3} samples\n"
        ).format(
            self.starttime,
            self.duration,
            self.rate,
            self.nsamp,
        )

    @staticmethod
    def _is_numeric_scalar(value):
        """
        Return whether a value is a supported numerical scalar.
        """
        return (
            isinstance(value, Number)
            and not isinstance(value, (bool, np.bool_))
        )

    @property
    def data(self) -> np.ndarray:
        """
        One-dimensional waveform data.
        """
        return self._data

    @data.setter
    def data(self, value: Any) -> None:
        array = np.asarray(value)

        if array.ndim == 0:
            array = array.reshape(1)

        if array.ndim != 1:
            raise ValueError(
                "record data must be one-dimensional"
            )

        self._data = array.copy()

    @property
    def sid(self) -> Optional[str]:
        """
        Source identifier.
        """
        return self.head.sid

    @sid.setter
    def sid(self, value: Optional[str]) -> None:
        self.head.sid = value

    @property
    def eid(self) -> Optional[str]:
        """
        Event identifier.
        """
        return self.head.eid

    @eid.setter
    def eid(self, value: Optional[str]) -> None:
        self.head.eid = value

    @property
    def nsamp(self) -> int:
        """
        Number of waveform samples.
        """
        return self.data.size

    @property
    def delta(self) -> Optional[float]:
        """
        Sampling interval in seconds.
        """
        return self.head.delta

    @delta.setter
    def delta(self, value: Optional[Real]) -> None:
        self.head.delta = value

    @property
    def rate(self) -> Optional[float]:
        """
        Sampling rate in samples per second.
        """
        return self.head.rate

    @rate.setter
    def rate(self, value: Optional[Real]) -> None:
        self.head.rate = value

    @property
    def time(self) -> Optional[Date]:
        """
        Start time of the record.
        """
        return self.head.time

    @time.setter
    def time(self, value: Any) -> None:
        self.head.time = value

    @property
    def peak(self) -> Optional[float]:
        """
        Maximum absolute sample amplitude.

        Returns
        -------
        number or None
            Maximum absolute amplitude, or ``None`` for an empty record.
        """
        if not self.nsamp:
            return None

        return np.max(np.abs(self.data))

    def copy(self) -> Record:
        """
        Return an independent copy of the record.

        Returns
        -------
        Record
            Copy containing independent waveform data and metadata.
        """
        record = Record()

        record.head = self.head.copy(parent=record)
        record.data = self.data.copy()

        return record

    def info(self) -> None:
        """
        Print a human-readable description of the record.
        """
        print(self)

    def __add__(self, value: Any) -> Record:
        """
        Return a new record after addition or concatenation.

        A numerical scalar is added to every waveform sample. Adding
        another Record attempts to append it temporally to a copy of the
        current record.

        Parameters
        ----------
        value : Number or Record
            Scalar value or record to add.

        Returns
        -------
        Record
            Independent result record.

        Raises
        ------
        TypeError
            If the operand is unsupported.
        ValueError
            If another Record cannot be appended.
        """
        record = self.copy()

        if self._is_numeric_scalar(value):
            record.data = record.data + value
            return record

        if isinstance(value, Record):
            appended = record.append(
                value.copy(),
            )

            if not appended:
                raise ValueError(
                    "records cannot be concatenated"
                )

            return record

        raise TypeError(
            "unsupported operand type for Record addition: {0}".format(
                type(value).__name__
            )
        )

    def __radd__(self, value: Any) -> Record:
        """
        Return the result of scalar addition.
        """
        return self.__add__(value)

    def __sub__(self, value: Any) -> Record:
        """
        Return a new record after subtracting a numerical scalar.

        Parameters
        ----------
        value : Number
            Scalar value subtracted from every waveform sample.

        Returns
        -------
        Record
            Independent result record.

        Raises
        ------
        TypeError
            If the operand is not a numerical scalar.
        """
        if not self._is_numeric_scalar(value):
            raise TypeError(
                "unsupported operand type for Record subtraction: "
                "{0}".format(type(value).__name__)
            )

        record = self.copy()
        record.data = record.data - value

        return record

    def __rsub__(self, value: Any) -> Record:
        """
        Return a scalar minus every waveform sample.
        """
        if not self._is_numeric_scalar(value):
            raise TypeError(
                "unsupported operand type for Record subtraction: "
                "{0}".format(type(value).__name__)
            )

        record = self.copy()
        record.data = value - record.data

        return record

    def __mul__(self, value: Any) -> Record:
        """
        Return a scaled or response-convolved record.

        Parameters
        ----------
        value : Number or response
            Numerical scale factor, or instrumental response: a stage,
            a stage set, a channel epoch or a metadata collection
            (see :meth:`convolve_response`).

        Returns
        -------
        Record
            Independent result record.

        Raises
        ------
        TypeError
            If the operand is unsupported.
        """
        record = self.copy()

        if self._is_numeric_scalar(value):
            record.data = record.data * value
            return record

        if isinstance(value, _response_types()):
            record.convolve_response(value)
            return record

        raise TypeError(
            "unsupported operand type for Record multiplication: "
            "{0}".format(type(value).__name__)
        )

    def __rmul__(self, value: Any) -> Record:
        """
        Return the result of scalar multiplication.
        """
        return self.__mul__(value)

    def __truediv__(self, value: Any) -> Record:
        """
        Return a scaled or response-deconvolved record.

        Parameters
        ----------
        value : Number or response
            Numerical divisor, or instrumental response: a stage, a
            stage set, a channel epoch or a metadata collection
            (see :meth:`deconvolve_response`).

        Returns
        -------
        Record
            Independent result record.

        Raises
        ------
        TypeError
            If the operand is unsupported.
        ZeroDivisionError
            If the numerical divisor is zero.
        """
        record = self.copy()

        if self._is_numeric_scalar(value):
            if value == 0:
                raise ZeroDivisionError(
                    "cannot divide a Record by zero"
                )

            record.data = record.data / value
            return record

        if isinstance(value, _response_types()):
            record.deconvolve_response(value)
            return record

        raise TypeError(
            "unsupported operand type for Record division: {0}".format(
                type(value).__name__
            )
        )

    @property
    def starttime(self) -> Optional[Date]:
        """
        Time of the first sample.

        Returns
        -------
        Date or None
            Record start time.
        """
        return self.head.time

    @property
    def endtime(self) -> Optional[Date]:
        """
        Time of the final sample.

        Returns
        -------
        Date or None
            Time corresponding to the final waveform sample. ``None`` is
            returned when the record has no start time or contains no
            samples.
        """
        if self.time is None or self.nsamp == 0:
            return None

        return self.time + self.duration

    @property
    def duration(self) -> float:
        """
        Time span between the first and final sample, in seconds.

        Returns
        -------
        float
            Record duration. An empty or one-sample record has zero
            duration.

        Notes
        -----
        Duration is defined as:

        ``(nsamp - 1) * delta``

        It therefore represents the time difference between the first and
        final samples, not the total width of ``nsamp`` sampling intervals.
        """
        if self.nsamp <= 1:
            return 0.0

        if self.delta is None:
            raise ValueError(
                "sampling interval is required to calculate duration"
            )

        return (self.nsamp - 1) * self.delta

    def _require_time_axis(self):
        """
        Validate the metadata required for temporal operations.

        Raises
        ------
        ValueError
            If the sampling interval or start time is undefined.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required for temporal operations"
            )

        if self.time is None:
            raise ValueError(
                "start time is required for temporal operations"
            )

    @staticmethod
    def _time_difference(time_1, time_0):
        """
        Return the difference between two times in seconds.

        Parameters
        ----------
        time_1, time_0 : Date
            Times to subtract.

        Returns
        -------
        float
            ``time_1 - time_0`` in seconds.
        """
        return float(time_1 - time_0)

    def _relative_time(self, value, name):
        """
        Convert a temporal value to seconds from record start.

        Parameters
        ----------
        value : Date, str or Real
            Absolute time or relative time in seconds.
        name : str
            Parameter name used in error messages.

        Returns
        -------
        float
            Time in seconds relative to the first record sample.

        Raises
        ------
        TypeError
            If the temporal representation is unsupported.
        """
        if isinstance(value, Date):
            self._require_time_axis()
            return self._time_difference(value, self.time)

        if isinstance(value, str):
            self._require_time_axis()
            return self._time_difference(Date(value), self.time)

        if (
            isinstance(value, Real)
            and not isinstance(value, (bool, np.bool_))
        ):
            value = float(value)

            if not np.isfinite(value):
                raise ValueError(
                    "{0} must be finite".format(name)
                )

            return value

        raise TypeError(
            "{0} must be a Date, string, or relative time in "
            "seconds".format(name)
        )

    def _relative_time_window(
        self,
        starttime=None,
        endtime=None,
    ):
        """
        Convert a requested time window to record-relative seconds.

        Returns
        -------
        tuple of float
            Relative start and end times.

        Raises
        ------
        ValueError
            If the requested interval is invalid.
        """
        relative_start = 0.0
        relative_end = self.duration

        if starttime is not None:
            relative_start = self._relative_time(
                starttime,
                "starttime",
            )

        if endtime is not None:
            relative_end = self._relative_time(
                endtime,
                "endtime",
            )

        if relative_end < relative_start:
            raise ValueError(
                "endtime precedes starttime"
            )

        return relative_start, relative_end

    def overlaps(self, starttime: Any = None, endtime: Any = None) -> bool:
        """
        Check whether a time window intersects the record.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the requested window.
        endtime : Date, str or float, optional
            End of the requested window.

        Returns
        -------
        bool
            ``True`` if the window contains at least one record sample.
        """
        if self.nsamp == 0:
            return False

        relative_start, relative_end = self._relative_time_window(
            starttime,
            endtime,
        )

        return not (
            relative_end < 0.0
            or relative_start > self.duration
        )

    def time_axis(
        self,
        reference: Any = "relative",
        shift: float = 0.0,
    ) -> np.ndarray:
        """
        Generate the record time axis, in seconds.

        Parameters
        ----------
        reference : {"relative", "absolute"} or Date, optional
            Origin of the time axis:

            - ``"relative"`` (or ``"r"``): the record start time, so the
              axis starts from zero (default);
            - ``"absolute"`` (or ``"a"``): the internal epoch of
              :class:`~shakelab.libutils.timeN.Date`
              (0001-01-01T00:00:00Z);
            - a :class:`~shakelab.libutils.timeN.Date`, or any value
              accepted by it (e.g. an ISO 8601 string): seconds relative
              to that time, useful to align several records.
        shift : float, optional
            Additional shift in seconds.

        Returns
        -------
        numpy.ndarray
            Time coordinate of each waveform sample.

        Raises
        ------
        ValueError
            If the sampling interval is missing, or if an absolute or
            date reference is requested without a record start time.

        Notes
        -----
        Absolute times are about 6e10 seconds today, so in double
        precision they resolve only about 10 microseconds. A date
        reference keeps the full precision: the offset between the
        record start and the reference is computed exactly before
        conversion to float.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required to generate a time axis"
            )

        if (
            not isinstance(shift, Real)
            or isinstance(shift, (bool, np.bool_))
        ):
            raise TypeError("shift must be a real number")

        shift = float(shift)

        if not np.isfinite(shift):
            raise ValueError("shift must be finite")

        axis = np.arange(self.nsamp, dtype=float) * self.delta

        if isinstance(reference, str) and reference in ("relative", "r"):
            return axis + shift

        if self.time is None:
            raise ValueError(
                "start time is required for an absolute or date-"
                "referenced time axis"
            )

        if isinstance(reference, str) and reference in ("absolute", "a"):
            return axis + self.time.to_seconds() + shift

        try:
            origin = Date(reference)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "reference must be 'relative', 'absolute' or a date"
            ) from error

        offset = float(self.time - origin)

        return axis + offset + shift

    @property
    def taxis(self) -> np.ndarray:
        """
        Relative time axis in seconds.
        """
        return self.time_axis()

    def append(
        self,
        record: Record,
        enforce: bool = False,
        fillvalue: float = 0.0,
        merge_tolerance: Optional[float] = None,
    ) -> bool:
        """
        Append another record to the current record.

        Parameters
        ----------
        record : Record
            Record to append.
        enforce : bool, optional
            If ``False``, only directly contiguous records are appended.
            If ``True``, integer-sample gaps are filled and integer-sample
            overlaps are resolved by retaining the current record up to
            the first sample of the new record.
        fillvalue : scalar, optional
            Value used to fill gaps when ``enforce=True``.
        merge_tolerance : float, optional
            Maximum distance from an integer sample offset, expressed as
            a fraction of one sampling interval. Default is
            :data:`MERGE_TOLERANCE` (half a sample).

        Returns
        -------
        bool
            ``True`` when the record was appended, otherwise ``False``.

        Raises
        ------
        TypeError
            If ``record`` is not a Record.
        ValueError
            If metadata or sampling parameters are incompatible.

        Notes
        -----
        The method modifies the current record in place.

        For overlapping records with ``enforce=True``, samples from the
        new record take precedence from its start time onward.
        """
        if not isinstance(record, Record):
            raise TypeError(
                "record must be a Record instance"
            )

        if not isinstance(enforce, bool):
            raise TypeError("enforce must be a boolean")

        if merge_tolerance is None:
            merge_tolerance = MERGE_TOLERANCE

        if (
            not isinstance(merge_tolerance, Real)
            or isinstance(merge_tolerance, (bool, np.bool_))
        ):
            raise TypeError("merge_tolerance must be a real number")

        merge_tolerance = float(merge_tolerance)

        if not np.isfinite(merge_tolerance) or merge_tolerance < 0.0:
            raise ValueError(
                "merge_tolerance must be finite and non-negative"
            )

        if self.nsamp == 0:
            self.head = record.head.copy(parent=self)
            self.data = record.data.copy()
            return True

        if record.nsamp == 0:
            return True

        self._require_time_axis()
        record._require_time_axis()

        if not np.isclose(
            self.delta,
            record.delta,
            rtol=1e-9,
            atol=0.0,
        ):
            raise ValueError(
                "records have incompatible sampling intervals"
            )

        if (
            self.sid is not None
            and record.sid is not None
            and self.sid != record.sid
        ):
            raise ValueError(
                "records have different source identifiers"
            )

        expected_time = self.endtime + self.delta
        offset_time = self._time_difference(
            record.starttime,
            expected_time,
        )
        offset_samples = offset_time / self.delta
        sample_offset = int(round(offset_samples))

        if not np.isclose(
            offset_samples,
            sample_offset,
            rtol=0.0,
            atol=merge_tolerance,
        ):
            if enforce:
                raise ValueError(
                    "record offset is not aligned with the sampling grid"
                )

            return False

        if sample_offset == 0:
            self.data = np.concatenate(
                (self.data, record.data)
            )
            return True

        if not enforce:
            return False

        if sample_offset > 0:
            gap = np.full(
                sample_offset,
                fillvalue,
                dtype=np.result_type(
                    self.data.dtype,
                    record.data.dtype,
                    type(fillvalue),
                ),
            )

            self.data = np.concatenate(
                (self.data, gap, record.data)
            )
            return True

        overlap = -sample_offset

        if overlap >= self.nsamp:
            raise ValueError(
                "new record starts before the current record"
            )

        keep_samples = self.nsamp - overlap

        self.data = np.concatenate(
            (
                self.data[:keep_samples],
                record.data,
            )
        )

        return True

    def _sample_window(self, starttime, endtime, snap_tolerance):
        """
        Return the indices of the first and last samples in a window.

        Parameters
        ----------
        starttime, endtime : Date, str, float or None
            Window bounds (see :meth:`cut`).
        snap_tolerance : float
            Snapping tolerance, as a fraction of the sampling interval.

        Returns
        -------
        tuple of int or None
            ``(first, last)`` sample indices, both included, or ``None``
            when no sample falls within the window.
        """
        if self.nsamp == 0:
            return None

        relative_start, relative_end = self._relative_time_window(
            starttime=starttime,
            endtime=endtime,
        )

        index_start = max(
            int(np.ceil(relative_start / self.delta - snap_tolerance)),
            0,
        )
        index_end = min(
            int(np.floor(relative_end / self.delta + snap_tolerance)),
            self.nsamp - 1,
        )

        if index_start > index_end:
            return None

        return index_start, index_end

    def cut(
        self,
        starttime: Any = None,
        endtime: Any = None,
        inplace: bool = True,
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> Optional[Record]:
        """
        Cut the record to the samples within a time window.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the window. Numerical values are interpreted as
            seconds from record start. Default is the record start.
        endtime : Date, str or float, optional
            End of the window, included. Numerical values are
            interpreted as seconds from record start. Default is the
            record end.
        inplace : bool, optional
            If ``True``, modify the current record. If ``False``, return
            an independent Record. Default is ``True``.
        snap_tolerance : float, optional
            Times closer than this to a sample, as a fraction of one
            sampling interval, are treated as coinciding with it.
            Default is :data:`SNAP_TOLERANCE`.

        Returns
        -------
        Record or None
            When ``inplace=False``, return the cut record. When
            ``inplace=True``, return ``None``.

        Raises
        ------
        TypeError
            If ``inplace`` is not a boolean or ``snap_tolerance`` is not a
            real number.
        ValueError
            If the snap_tolerance is invalid, the sampling interval is
            missing, the window is reversed, or no sample of the record
            falls within the window.

        Notes
        -----
        The retained samples are those whose times lie within the
        window, bounds included: the first sample at or after
        ``starttime`` and the last sample at or before ``endtime``. The
        window does not need to be aligned with the sampling grid; the
        start time of the result is the time of its first sample.
        """
        if not isinstance(inplace, bool):
            raise TypeError(
                "inplace must be a boolean"
            )

        if (
            not isinstance(snap_tolerance, Real)
            or isinstance(snap_tolerance, (bool, np.bool_))
        ):
            raise TypeError(
                "snap_tolerance must be a real number"
            )

        snap_tolerance = float(snap_tolerance)

        if not np.isfinite(snap_tolerance) or snap_tolerance < 0.0:
            raise ValueError(
                "snap_tolerance must be finite and non-negative"
            )

        if self.nsamp == 0:
            if inplace:
                return None

            return self.copy()

        if self.delta is None:
            raise ValueError(
                "sampling interval is required to cut a record"
            )

        window = self._sample_window(starttime, endtime, snap_tolerance)

        if window is None:
            raise ValueError(
                "no sample of the record falls within the requested "
                "time window"
            )

        index_start, index_end = window

        target = self if inplace else self.copy()

        # The final requested sample is included.
        target.data = self.data[index_start:index_end + 1]

        if self.time is not None:
            target.time = (
                self.time
                + index_start * self.delta
            )

        if inplace:
            return None

        return target

    def extract(
        self,
        starttime: Any = None,
        endtime: Any = None,
        snap_tolerance: float = SNAP_TOLERANCE,
    ) -> Record:
        """
        Return an independent record with the samples of a time window.

        See :meth:`cut`, of which this is the non-destructive variant.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the window.
        endtime : Date, str or float, optional
            End of the window, included.
        snap_tolerance : float, optional
            Times closer than this to a sample, as a fraction of one
            sampling interval, are treated as coinciding with it.

        Returns
        -------
        Record
            Extracted record.
        """
        return self.cut(
            starttime=starttime,
            endtime=endtime,
            inplace=False,
            snap_tolerance=snap_tolerance,
        )

    def remove_mean(self) -> None:
        """
        Remove the arithmetic mean from the waveform data.

        Notes
        -----
        The operation is performed in place. Empty records are left
        unchanged.
        """
        if self.nsamp == 0:
            return

        self.data = self.data - np.mean(self.data)

    def filter(
        self,
        highpass: Optional[float] = None,
        lowpass: Optional[float] = None,
        order: int = 4,
        minphase: bool = False,
    ) -> None:
        """
        Apply a Butterworth filter to the waveform data.

        Parameters
        ----------
        highpass : float, optional
            High-pass corner frequency in hertz.
        lowpass : float, optional
            Low-pass corner frequency in hertz.
        order : int, optional
            Butterworth filter order. Default is 4.
        minphase : bool, optional
            If ``False``, apply forward-backward zero-phase filtering.
            If ``True``, apply a single-pass causal filter. Default is
            ``False``.

        Raises
        ------
        ValueError
            If no corner frequency is provided, if the filter parameters
            are invalid, or if the record is too short for zero-phase
            filtering.

        Notes
        -----
        The operation is performed in place.

        Despite the historical parameter name, ``minphase=True`` applies
        a causal single-pass Butterworth filter. It does not construct a
        true minimum-phase equivalent filter.
        """
        if self.nsamp == 0:
            return

        if self.delta is None:
            raise ValueError(
                "sampling interval is required for filtering"
            )

        if highpass is None and lowpass is None:
            raise ValueError(
                "at least one corner frequency must be provided"
            )

        if not isinstance(order, Integral) or isinstance(
            order, (bool, np.bool_)
        ):
            raise TypeError("order must be an integer")

        order = int(order)

        if order <= 0:
            raise ValueError("order must be strictly positive")

        if not isinstance(minphase, bool):
            raise TypeError("minphase must be a boolean")

        nyquist = 0.5 * self.rate

        def validate_corner(value, name):
            if value is None:
                return None

            if (
                not isinstance(value, Real)
                or isinstance(value, (bool, np.bool_))
            ):
                raise TypeError(
                    "{0} must be a real number".format(name)
                )

            value = float(value)

            if not np.isfinite(value):
                raise ValueError(
                    "{0} must be finite".format(name)
                )

            if value <= 0.0:
                raise ValueError(
                    "{0} must be strictly positive".format(name)
                )

            if value >= nyquist:
                raise ValueError(
                    "{0} must be lower than the Nyquist frequency "
                    "({1} Hz)".format(name, nyquist)
                )

            return value

        highpass = validate_corner(highpass, "highpass")
        lowpass = validate_corner(lowpass, "lowpass")

        if (
            highpass is not None
            and lowpass is not None
            and highpass >= lowpass
        ):
            raise ValueError(
                "highpass must be lower than lowpass"
            )

        if highpass is not None and lowpass is not None:
            cutoff = (highpass, lowpass)
            filter_type = "bandpass"
        elif highpass is not None:
            cutoff = highpass
            filter_type = "highpass"
        else:
            cutoff = lowpass
            filter_type = "lowpass"

        sos = signal.butter(
            order,
            cutoff,
            btype=filter_type,
            fs=self.rate,
            output="sos",
        )

        if minphase:
            self.data = signal.sosfilt(
                sos,
                self.data,
            )
            return

        try:
            self.data = signal.sosfiltfilt(
                sos,
                self.data,
            )
        except ValueError as error:
            raise ValueError(
                "record is too short for zero-phase filtering"
            ) from error

    def resample(
        self,
        new_delta: Real,
        method: str = "interp",
        antialias: bool = True,
    ) -> None:
        """
        Resample the waveform to a new sampling interval.

        Parameters
        ----------
        new_delta : float
            New sampling interval in seconds.
        method : {"scipy", "fft", "interp"}, optional
            Resampling method:

            - ``"scipy"`` or ``"fft"`` uses Fourier-domain resampling
              through :func:`scipy.signal.resample`. The ratio between
              the new and the old sampling interval must be a fraction
              with a denominator up to 1000 (e.g. 2, 1/2, 3/2).
            - ``"interp"`` uses linear interpolation.

            Default is ``"interp"``.
        antialias : bool, optional
            Apply a low-pass anti-alias filter before interpolation when
            downsampling. Fourier-domain resampling already performs the
            corresponding spectral truncation. Default is ``True``.

        Raises
        ------
        ValueError
            If the sampling interval or method is invalid, or if the
            record is too short for the requested anti-alias filter.

        Notes
        -----
        The operation is performed in place. The start time is preserved.

        The output samples cover the original record without
        extrapolation:

        ``new_nsamp = floor(duration / new_delta) + 1``

        so the final sample lies at or before the original final
        sample, by less than one new sampling interval.

        Fourier-domain resampling treats the record as periodic, so
        records that do not start and end smoothly (e.g. not tapered)
        show oscillations near their ends. Before the transform the
        record is extended with zeros by fewer samples than the
        denominator of the interval ratio, so that its total length is
        a whole number of new sampling intervals.
        """
        if (
            not isinstance(new_delta, Real)
            or isinstance(new_delta, (bool, np.bool_))
        ):
            raise TypeError(
                "new_delta must be a real number"
            )

        new_delta = float(new_delta)

        if not np.isfinite(new_delta):
            raise ValueError(
                "new_delta must be finite"
            )

        if new_delta <= 0.0:
            raise ValueError(
                "new_delta must be strictly positive"
            )

        if not isinstance(antialias, bool):
            raise TypeError(
                "antialias must be a boolean"
            )

        method = str(method).lower()

        if method not in ("scipy", "fft", "interp"):
            raise ValueError(
                "method must be 'scipy', 'fft', or 'interp'"
            )

        if self.delta is None:
            raise ValueError(
                "sampling interval is required for resampling"
            )

        if self.nsamp == 0:
            self.delta = new_delta
            return

        if np.isclose(
            self.delta,
            new_delta,
            rtol=1e-12,
            atol=0.0,
        ):
            self.delta = new_delta
            return

        if self.nsamp == 1:
            self.delta = new_delta
            return

        old_delta = self.delta
        old_duration = self.duration

        new_nsamp = int(
            np.floor(old_duration / new_delta + SNAP_TOLERANCE)
        ) + 1

        if method in ("scipy", "fft"):
            self.data = _fft_resample(
                self.data,
                old_delta,
                new_delta,
                new_nsamp,
            )

        else:
            data = self.data

            if antialias and new_delta > old_delta:
                new_nyquist = 0.5 / new_delta

                # Keep the cutoff slightly below the new Nyquist
                # frequency to provide a finite transition band.
                cutoff = 0.9 * new_nyquist

                sos = signal.butter(
                    4,
                    cutoff,
                    btype="lowpass",
                    fs=1.0 / old_delta,
                    output="sos",
                )

                try:
                    data = signal.sosfiltfilt(
                        sos,
                        data,
                    )
                except ValueError as error:
                    raise ValueError(
                        "record is too short for anti-alias filtering"
                    ) from error

            old_time = np.arange(
                self.nsamp,
                dtype=float,
            ) * old_delta

            new_time = np.arange(
                new_nsamp,
                dtype=float,
            ) * new_delta

            # Guard against floating-point rounding of the final time.
            new_time = np.minimum(
                new_time,
                old_duration,
            )

            if np.iscomplexobj(data):
                real = np.interp(
                    new_time,
                    old_time,
                    data.real,
                )
                imag = np.interp(
                    new_time,
                    old_time,
                    data.imag,
                )
                self.data = real + 1j * imag
            else:
                self.data = np.interp(
                    new_time,
                    old_time,
                    data,
                )

        self.delta = new_delta

    def taper(self, time: float = 0.1) -> None:
        """
        Apply a symmetric Tukey taper to the waveform.

        Parameters
        ----------
        time : float, optional
            Taper duration at each end of the record, in seconds.
            A negative value applies a full cosine taper. Default is
            0.1 seconds.

        Notes
        -----
        The operation is performed in place. If the requested taper
        duration exceeds half the record duration, a full cosine taper is
        applied.
        """
        if self.nsamp == 0:
            return

        if self.delta is None:
            raise ValueError(
                "sampling interval is required for tapering"
            )

        if (
            not isinstance(time, Real)
            or isinstance(time, (bool, np.bool_))
        ):
            raise TypeError(
                "time must be a real number"
            )

        time = float(time)

        if not np.isfinite(time):
            raise ValueError(
                "time must be finite"
            )

        if time < 0.0:
            alpha = 1.0
        elif time == 0.0:
            return
        else:
            total_window = self.nsamp * self.delta
            alpha = min(
                2.0 * time / total_window,
                1.0,
            )

        self.data = self.data * tukey_window(
            self.nsamp,
            alpha,
        )

    def zero_padding(self, time: Optional[float] = None) -> None:
        """
        Append zero-valued samples to the waveform.

        Parameters
        ----------
        time : float, optional
            Padding duration in seconds. If omitted, append the same
            number of samples currently contained in the record.

        Raises
        ------
        ValueError
            If the requested duration does not correspond to an integer
            number of samples.

        Notes
        -----
        The operation is performed in place.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required for zero padding"
            )

        if time is None:
            zero_count = self.nsamp

        else:
            if (
                not isinstance(time, Real)
                or isinstance(time, (bool, np.bool_))
            ):
                raise TypeError(
                    "time must be a real number"
                )

            time = float(time)

            if not np.isfinite(time):
                raise ValueError(
                    "time must be finite"
                )

            if time < 0.0:
                raise ValueError(
                    "time must be non-negative"
                )

            sample_count = time / self.delta
            zero_count = int(round(sample_count))

            if not np.isclose(
                sample_count,
                zero_count,
                rtol=0.0,
                atol=SNAP_TOLERANCE,
            ):
                raise ValueError(
                    "padding duration is not aligned with the "
                    "sampling grid"
                )

        if zero_count == 0:
            return

        zeros = np.zeros(
            zero_count,
            dtype=self.data.dtype,
        )

        self.data = np.concatenate(
            (self.data, zeros)
        )

    def shift(self, time: float, padding: bool = True) -> None:
        """
        Shift the waveform relative to its fixed time window.

        Parameters
        ----------
        time : float
            Time shift in seconds.
        padding : bool, optional
            If ``True``, temporarily append zeros to reduce circular
            wrap-around. Default is ``True``.

        Notes
        -----
        The operation is performed in place. The record start time is not
        changed: waveform samples are shifted within the existing time
        window.

        The Fourier-domain implementation supports fractional-sample
        shifts.
        """
        if self.nsamp == 0:
            return

        if self.delta is None:
            raise ValueError(
                "sampling interval is required for shifting"
            )

        if (
            not isinstance(time, Real)
            or isinstance(time, (bool, np.bool_))
        ):
            raise TypeError(
                "time must be a real number"
            )

        time = float(time)

        if not np.isfinite(time):
            raise ValueError(
                "time must be finite"
            )

        if not isinstance(padding, bool):
            raise TypeError(
                "padding must be a boolean"
            )

        original_nsamp = self.nsamp

        if padding:
            data = np.concatenate(
                (
                    self.data,
                    np.zeros(
                        original_nsamp,
                        dtype=self.data.dtype,
                    ),
                )
            )
        else:
            data = self.data

        if np.iscomplexobj(data):
            # The shift is linear: real and imaginary parts are shifted
            # separately by the real-valued Fourier implementation.
            shifted = (
                fourier.shift_time(data.real, self.delta, time)
                + 1j * fourier.shift_time(data.imag, self.delta, time)
            )
        else:
            shifted = fourier.shift_time(data, self.delta, time)

        self.data = shifted[:original_nsamp]

    def integrate(self, method: str = "fft") -> None:
        """
        Integrate the waveform in time.

        Parameters
        ----------
        method : {"fft", "cum"}, optional
            Integration method:

            - ``"fft"`` performs periodic Fourier-domain integration.
            - ``"cum"`` performs cumulative trapezoidal integration.

            Default is ``"fft"``.

        Notes
        -----
        The operation is performed in place.

        Fourier-domain integration removes the zero-frequency component,
        because a constant integration term cannot be determined from the
        waveform alone.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required for integration"
            )

        method = str(method).lower()

        if method not in ("fft", "cum"):
            raise ValueError(
                "method must be 'fft' or 'cum'"
            )

        if self.nsamp == 0:
            return

        if self.nsamp == 1:
            self.data = np.zeros(
                1,
                dtype=np.result_type(
                    self.data.dtype,
                    float,
                ),
            )
            return

        if method == "cum":
            self.data = integrate.cumulative_trapezoid(
                self.data,
                dx=self.delta,
                initial=0.0,
            )
            return

        self.data = _fft_derivative(
            self.data,
            order=-1,
            period=self.nsamp * self.delta,
        )

    def differentiate(self, method: str = "fft") -> None:
        """
        Differentiate the waveform in time.

        Parameters
        ----------
        method : {"fft", "grad"}, optional
            Differentiation method:

            - ``"fft"`` performs periodic Fourier-domain
              differentiation.
            - ``"grad"`` uses finite differences through
              :func:`numpy.gradient`.

            Default is ``"fft"``.

        Notes
        -----
        The operation is performed in place.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required for differentiation"
            )

        method = str(method).lower()

        if method not in ("fft", "grad"):
            raise ValueError(
                "method must be 'fft' or 'grad'"
            )

        if self.nsamp == 0:
            return

        if self.nsamp == 1:
            self.data = np.zeros(
                1,
                dtype=np.result_type(
                    self.data.dtype,
                    float,
                ),
            )
            return

        if method == "grad":
            self.data = np.gradient(
                self.data,
                self.delta,
            )
            return

        self.data = _fft_derivative(
            self.data,
            order=1,
            period=self.nsamp * self.delta,
        )

    @staticmethod
    def _as_waveform_data(value, name="data"):
        """
        Convert an input object to a one-dimensional waveform array.

        Parameters
        ----------
        value : Record or array-like
            Input waveform or kernel.
        name : str, optional
            Name used in error messages.

        Returns
        -------
        numpy.ndarray
            Independent one-dimensional array.

        Raises
        ------
        TypeError
            If the input cannot be converted to an array.
        ValueError
            If the resulting array is not one-dimensional or is empty.
        """
        if isinstance(value, Record):
            array = value.data
        else:
            try:
                array = np.asarray(value)
            except Exception as error:
                raise TypeError(
                    "{0} must be a Record or an array-like "
                    "object".format(name)
                ) from error

        if array.ndim == 0:
            array = array.reshape(1)

        if array.ndim != 1:
            raise ValueError(
                "{0} must be one-dimensional".format(name)
            )

        if array.size == 0:
            raise ValueError(
                "{0} cannot be empty".format(name)
            )

        return array.copy()

    @staticmethod
    def _validate_signal_mode(mode):
        """
        Validate a SciPy convolution or correlation mode.
        """
        if mode not in ("full", "same", "valid"):
            raise ValueError(
                "mode must be 'full', 'same', or 'valid'"
            )

    @staticmethod
    def _validate_signal_method(method):
        """
        Validate a SciPy convolution or correlation method.
        """
        if method not in ("auto", "direct", "fft"):
            raise ValueError(
                "method must be 'auto', 'direct', or 'fft'"
            )

    def _replace_from_record(self, record):
        """
        Replace data and metadata using another Record.

        Parameters
        ----------
        record : Record
            Source record.

        Notes
        -----
        The header is copied and associated with the current Record,
        preventing stale parent references.
        """
        if not isinstance(record, Record):
            raise TypeError(
                "result must be a Record instance"
            )

        self.head = record.head.copy(parent=self)
        self.data = record.data.copy()

    def to_spectrum(self) -> Spectrum:
        """
        Convert the record to a frequency-domain Spectrum.

        Returns
        -------
        Spectrum
            Spectrum generated from the current record.

        Notes
        -----
        The conversion is implemented by
        :meth:`shakelab.signals.spectrum.Spectrum.from_record`; the
        spectrum receives a copy of the record header.
        """
        # Spectrum sits above Record in the import hierarchy: the
        # conversion is delegated with a local import.
        from shakelab.signals.spectrum import Spectrum

        return Spectrum.from_record(self)

    def from_spectrum(self, spectrum: Spectrum) -> None:
        """
        Replace the record with the inverse transform of a spectrum.

        Parameters
        ----------
        spectrum : object
            Object exposing a callable ``to_record()`` method.

        Raises
        ------
        TypeError
            If the supplied object does not expose ``to_record()`` or if
            the conversion does not return a Record.

        Notes
        -----
        The operation is performed in place. The converted header is
        copied so that its parent always refers to the current Record.
        """
        to_record = getattr(spectrum, "to_record", None)

        if not callable(to_record):
            raise TypeError(
                "spectrum must expose a callable to_record() method"
            )

        record = to_record()

        if not isinstance(record, Record):
            raise TypeError(
                "spectrum.to_record() must return a Record"
            )

        self._replace_from_record(record)

    def convolve(
        self,
        data: Any,
        mode: str = "full",
        method: str = "auto",
    ) -> None:
        """
        Convolve the waveform with another sequence.

        Parameters
        ----------
        data : Record or array-like
            Convolution kernel.
        mode : {"full", "same", "valid"}, optional
            Output-size mode accepted by :func:`scipy.signal.convolve`.
            Default is ``"full"``.
        method : {"auto", "direct", "fft"}, optional
            Convolution method. Default is ``"auto"``.

        Raises
        ------
        ValueError
            If the kernel is invalid or, when supplied as a Record, has
            an incompatible sampling interval.

        Notes
        -----
        The operation is performed in place.

        The record start time is preserved. A kernel is interpreted as a
        discrete filter whose first coefficient has zero lag. Its own
        start time, when supplied as a Record, is therefore ignored.
        """
        mode = str(mode).lower()
        method = str(method).lower()

        self._validate_signal_mode(mode)
        self._validate_signal_method(method)

        if self.nsamp == 0:
            return

        if isinstance(data, Record):
            if (
                self.delta is not None
                and data.delta is not None
                and not np.isclose(
                    self.delta,
                    data.delta,
                    rtol=1e-9,
                    atol=0.0,
                )
            ):
                raise ValueError(
                    "records have incompatible sampling intervals"
                )

        kernel = self._as_waveform_data(
            data,
            name="convolution kernel",
        )

        self.data = signal.convolve(
            self.data,
            kernel,
            mode=mode,
            method=method,
        )

    def deconvolve(self, data: Any) -> np.ndarray:
        """
        Deconvolve a discrete sequence from the waveform.

        Parameters
        ----------
        data : Record or array-like
            Divisor sequence representing the convolution kernel.

        Returns
        -------
        numpy.ndarray
            Remainder of the polynomial deconvolution.

        Raises
        ------
        ValueError
            If the divisor is invalid, longer than the waveform, or has
            an incompatible sampling interval.

        Notes
        -----
        The quotient replaces the waveform data in place. The remainder
        is returned so that the caller can assess the reconstruction
        error.

        This method performs exact discrete polynomial deconvolution
        through :func:`scipy.signal.deconvolve`; it is not a regularized
        spectral inverse.
        """
        if self.nsamp == 0:
            return np.array([], dtype=float)

        if isinstance(data, Record):
            if (
                self.delta is not None
                and data.delta is not None
                and not np.isclose(
                    self.delta,
                    data.delta,
                    rtol=1e-9,
                    atol=0.0,
                )
            ):
                raise ValueError(
                    "records have incompatible sampling intervals"
                )

        divisor = self._as_waveform_data(
            data,
            name="deconvolution divisor",
        )

        if divisor.size > self.nsamp:
            raise ValueError(
                "deconvolution divisor cannot be longer than the "
                "waveform"
            )

        if np.all(divisor == 0):
            raise ValueError(
                "deconvolution divisor cannot contain only zeros"
            )

        quotient, remainder = signal.deconvolve(
            self.data,
            divisor,
        )

        self.data = quotient

        return remainder

    def correlate(
        self,
        record: Record,
        mode: str = "full",
        method: str = "auto",
    ) -> None:
        """
        Cross-correlate the waveform with another sequence.

        Parameters
        ----------
        record : Record or array-like
            Reference waveform.
        mode : {"full", "same", "valid"}, optional
            Output-size mode accepted by :func:`scipy.signal.correlate`.
            Default is ``"full"``.
        method : {"auto", "direct", "fft"}, optional
            Correlation method. Default is ``"auto"``.

        Raises
        ------
        ValueError
            If the reference waveform is invalid or has an incompatible
            sampling interval.

        Notes
        -----
        The operation is performed in place.

        The output is a lag-domain sequence rather than an ordinary
        absolute-time waveform. For backward compatibility, the original
        start time is preserved. Lag coordinates should be calculated
        explicitly from the input lengths and sampling interval.
        """
        mode = str(mode).lower()
        method = str(method).lower()

        self._validate_signal_mode(mode)
        self._validate_signal_method(method)

        if self.nsamp == 0:
            return

        if isinstance(record, Record):
            if (
                self.delta is not None
                and record.delta is not None
                and not np.isclose(
                    self.delta,
                    record.delta,
                    rtol=1e-9,
                    atol=0.0,
                )
            ):
                raise ValueError(
                    "records have incompatible sampling intervals"
                )

        reference = self._as_waveform_data(
            record,
            name="correlation reference",
        )

        self.data = signal.correlate(
            self.data,
            reference,
            mode=mode,
            method=method,
        )

    def convolve_response(self, resp: Any, stages: Any = None) -> None:
        """
        Apply an instrument response to the waveform.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response object or response container. The applicable
            response is selected by source identifier and start time
            (see :func:`shakelab.signals.metadata.select_response`).
        stages : str or sequence of str, optional
            Stage types to include, e.g. ``("gain", "paz")``; None
            (default) for the complete response.

        Notes
        -----
        The operation is performed in place, in the frequency domain
        (see :meth:`shakelab.signals.spectrum.Spectrum.convolve_response`).
        The units become the output units of the response, if declared.
        """
        spectrum = self.to_spectrum()
        spectrum.convolve_response(resp, stages=stages)
        self._replace_from_record(spectrum.to_record())

    def deconvolve_response(
        self,
        resp: Any,
        waterlevel: float = 60.0,
        method: str = "water",
        stages: Any = None,
    ) -> None:
        """
        Remove an instrument response from the waveform.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response object or response container. The applicable
            response is selected by source identifier and start time
            (see :func:`shakelab.signals.metadata.select_response`).
        waterlevel : float, optional
            Water level in dB below the maximum of the response
            amplitude.
        method : {"water", "sharp", "smooth"}, optional
            Regularization of the inverse (see
            :func:`shakelab.signals.fourier.inverse_spectrum`).
        stages : str or sequence of str, optional
            Stage types to include, e.g. ``("gain", "paz")``; None
            (default) for the complete response.

        Notes
        -----
        The operation is performed in place, in the frequency domain
        (see
        :meth:`shakelab.signals.spectrum.Spectrum.deconvolve_response`).
        The units become the input units of the response, if declared.
        The record should be detrended and tapered beforehand.
        """
        spectrum = self.to_spectrum()
        spectrum.deconvolve_response(
            resp, waterlevel=waterlevel, method=method, stages=stages
        )
        self._replace_from_record(spectrum.to_record())

    def remove_sensitivity(self, resp: Any) -> None:
        """
        Divide the waveform by the overall sensitivity, in place.

        Parameters
        ----------
        resp : StageGain, StageSet, ChannelEpoch or MetadataCollection
            Response object or response container (see
            :func:`shakelab.signals.metadata.select_sensitivity`).

        Notes
        -----
        A quick correction for instruments with a flat response in the
        band of interest; no frequency-domain processing is involved.
        The units become the input units of the sensitivity, if
        declared.
        """
        sens = metadata.select_sensitivity(
            resp, sid=self.sid, time=self.time
        )
        self.data = self.data / sens.value
        if sens.input_units:
            self.head.units = sens.input_units


class Stream:
    """
    Collection of waveform records belonging to one source.

    A Stream groups one or more :class:`Record` objects sharing the same
    source identifier. Records may be contiguous, separated by gaps, or
    partially overlapping.

    Parameters
    ----------
    sid : str, optional
        Source identifier associated with the stream.

    Notes
    -----
    Records are stored in insertion order. Methods that require temporal
    ordering either sort explicitly or validate the current order.

    A Stream does not require all records to have an event identifier,
    but non-``None`` event identifiers should be unique when records are
    accessed by ``eid``.
    """

    def __init__(self, sid: Optional[str] = None) -> None:
        self.sid = sid
        self.record = self._new_record_container()

    def _new_record_container(self, records=None):
        """
        Create a new record container.

        Parameters
        ----------
        records : iterable of Record, optional
            Initial records.

        Returns
        -------
        list
            Container initialized with the provided records.
        """
        if records is None:
            records = []

        return list(records)

    def _replace_records(self, records):
        """
        Replace the current record container.

        Parameters
        ----------
        records : iterable of Record
            New records.
        """
        self.record = self._new_record_container(records)

    def __len__(self) -> int:
        """
        Return the number of records.
        """
        return len(self.record)

    def __iter__(self) -> Iterator[Record]:
        """
        Iterate over the records in insertion order.
        """
        return iter(self.record)

    def __getitem__(self, key: Any) -> Any:
        """
        Return one record, or a stream of records.

        Parameters
        ----------
        key : int, slice or str
            An integer selects a record by position and a string by
            event identifier. A slice returns a new :class:`Stream`
            (also for a circular stream) with the same source
            identifier and the selected records.

        Returns
        -------
        Record or Stream
            Selected record, or stream of the selected records.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        KeyError
            If a requested event identifier is not found or is
            ambiguous.

        Notes
        -----
        As for Python lists, a slice contains the stored records, not
        copies: changing them changes this stream. Use
        :meth:`Stream.copy` on the result for an independent selection.
        """
        if isinstance(key, str):
            return self.get_record(key)

        if isinstance(key, slice):
            stream = Stream(sid=self.sid)
            stream._replace_records(list(self.record)[key])
            return stream

        if isinstance(key, Integral) and not isinstance(
            key, (bool, np.bool_)
        ):
            return self.record[key]

        raise TypeError(
            "stream indices must be integers, slices, or event "
            "identifiers"
        )

    def __repr__(self) -> str:
        return f"Stream(sid={self.sid!r}, records={len(self.record)})"

    def __str__(self) -> str:
        if not self.record:
            return "empty stream: source identifier {0}".format(
                self.sid
            )

        lines = [
            "source identifier: {0}".format(self.sid),
            "number of records: {0}".format(len(self)),
        ]

        for index, record in enumerate(self.record):
            description = str(record).rstrip()

            lines.append(
                "Record {0}: {1}".format(
                    index,
                    description,
                )
            )

        return "\n".join(lines)

    @property
    def eid(self) -> list:
        """
        Event identifiers associated with the records.

        Returns
        -------
        list
            Event identifiers in record order. Missing identifiers are
            represented by ``None``.
        """
        return [
            record.eid
            for record in self.record
        ]

    @property
    def nsamp(self) -> int:
        """
        Total number of samples stored in the stream.

        Returns
        -------
        int
            Sum of the sample counts of all records.
        """
        return sum(
            record.nsamp
            for record in self.record
        )

    @property
    def starttime(self) -> Optional[Date]:
        """
        Earliest record start time.

        Returns
        -------
        Date or None
            Earliest available start time, or ``None`` for an empty
            stream.

        Raises
        ------
        ValueError
            If one or more records have no start time.
        """
        if not self.record:
            return None

        times = [
            record.starttime
            for record in self.record
        ]

        if any(time is None for time in times):
            raise ValueError(
                "all records must have a start time"
            )

        return min(times)

    @property
    def endtime(self) -> Optional[Date]:
        """
        Latest record end time.

        Returns
        -------
        Date or None
            Latest available end time, or ``None`` for an empty stream.

        Raises
        ------
        ValueError
            If one or more records have no valid time axis.
        """
        if not self.record:
            return None

        times = [
            record.endtime
            for record in self.record
        ]

        if any(time is None for time in times):
            raise ValueError(
                "all records must have a valid time axis"
            )

        return max(times)

    def _record_indices(self, eid):
        """
        Return all record indices matching an event identifier.

        Parameters
        ----------
        eid : str
            Event identifier to find.

        Returns
        -------
        list of int
            Matching record indices.
        """
        return [
            index
            for index, record in enumerate(self.record)
            if record.eid == eid
        ]

    def get_record(self, eid: str) -> Record:
        """
        Return the record associated with an event identifier.

        Parameters
        ----------
        eid : str
            Event identifier.

        Returns
        -------
        Record
            Matching record.

        Raises
        ------
        TypeError
            If ``eid`` is not a string.
        KeyError
            If no record matches the identifier or if the identifier is
            associated with multiple records.
        """
        if not isinstance(eid, str):
            raise TypeError(
                "eid must be a string"
            )

        indices = self._record_indices(eid)

        if not indices:
            raise KeyError(
                "event identifier not found: {0}".format(eid)
            )

        if len(indices) > 1:
            raise KeyError(
                "event identifier is not unique: {0}".format(eid)
            )

        return self.record[indices[0]]

    def info(self) -> None:
        """
        Print a human-readable description of the stream.
        """
        print(self)

    def copy(self) -> Stream:
        """
        Return an independent copy of the stream.

        Returns
        -------
        Stream
            Stream containing independent copies of all records.
        """
        stream = self.__class__(
            sid=deepcopy(self.sid),
        )

        stream._replace_records(
            record.copy()
            for record in self.record
        )

        return stream

    def append(
        self,
        record: Record,
        enforce: bool = False,
        fillvalue: float = 0.0,
        merge_tolerance: Optional[float] = None,
    ) -> None:
        """
        Append a record to the stream.

        Parameters
        ----------
        record : Record
            Record to append.
        enforce : bool, optional
            If ``False``, merge the record with the last stream record
            only when they are directly contiguous. If ``True``, allow
            integer-sample gaps and overlaps to be resolved by
            :meth:`Record.append`. Default is ``False``.
        fillvalue : scalar, optional
            Value used to fill gaps when ``enforce=True``.
        merge_tolerance : float, optional
            Maximum distance from an integer sample offset, expressed as
            a fraction of one sampling interval. Default is
            :data:`MERGE_TOLERANCE`.

        Raises
        ------
        TypeError
            If ``record`` is not a Record or ``enforce`` is not a
            boolean.
        ValueError
            If the record source identifier is incompatible with the
            stream or if record merging fails.

        Notes
        -----
        The stream stores a copy of the record, so the record passed
        by the caller is never modified, neither when it is stored nor
        when later records are merged with it.

        Records are preserved in insertion order. The new record is
        merged only with the final record currently stored in the
        stream. To guarantee meaningful temporal merging, the stream
        should be sorted before appending out-of-order data.
        """
        if not isinstance(record, Record):
            raise TypeError(
                "record must be a Record instance"
            )

        if not isinstance(enforce, bool):
            raise TypeError(
                "enforce must be a boolean"
            )

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

        record = record.copy()

        if record.sid is None:
            record.sid = self.sid

        if not self.record:
            self.record.append(record)
            return

        merged = self.record[-1].append(
            record,
            enforce=enforce,
            fillvalue=fillvalue,
            merge_tolerance=merge_tolerance,
        )

        if not merged:
            self.record.append(record)

    def remove(self, key: Any) -> Record:
        """
        Remove a record from the stream.

        Parameters
        ----------
        key : int or str
            Integer index or event identifier.

        Returns
        -------
        Record
            Removed record.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        IndexError
            If the integer index is out of range.
        KeyError
            If the event identifier is not found or is ambiguous.
        """
        if isinstance(key, str):
            indices = self._record_indices(key)

            if not indices:
                raise KeyError(
                    "event identifier not found: {0}".format(key)
                )

            if len(indices) > 1:
                raise KeyError(
                    "event identifier is not unique: {0}".format(key)
                )

            key = indices[0]

        elif isinstance(key, (bool, np.bool_)) or not isinstance(
            key, Integral
        ):
            raise TypeError(
                "record key must be an integer or an event identifier"
            )

        records = list(self.record)
        removed = records.pop(key)
        self._replace_records(records)

        return removed

    def sort(self, reverse: bool = False) -> None:
        """
        Sort records by start time.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, sort records from latest to earliest. Default is
            ``False``.

        Raises
        ------
        TypeError
            If ``reverse`` is not a boolean.
        ValueError
            If one or more records have no start time.

        Notes
        -----
        Sorting is performed in place.
        """
        if not isinstance(reverse, bool):
            raise TypeError(
                "reverse must be a boolean"
            )

        if any(
            record.starttime is None
            for record in self.record
        ):
            raise ValueError(
                "all records must have a start time"
            )

        records = sorted(
            self.record,
            key=lambda record: record.starttime,
            reverse=reverse,
        )

        self._replace_records(records)

    def extract(
        self,
        starttime: Any = None,
        endtime: Any = None,
        eid: Optional[str] = None,
        merge: bool = False,
        enforce: bool = False,
        fillvalue: float = 0.0,
        snap_tolerance: float = SNAP_TOLERANCE,
        merge_tolerance: Optional[float] = None,
    ) -> Any:
        """
        Extract the samples of the records within a time window.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the requested time window. Numerical values are
            seconds from the start of the stream (its earliest record).
        endtime : Date, str or float, optional
            End of the requested time window, included. Numerical
            values are seconds from the start of the stream.
        eid : str, optional
            Restrict extraction to the record associated with this event
            identifier; numerical times are then relative to that
            record.
        merge : bool, optional
            If ``False``, return a new Stream containing all extracted
            record segments. If ``True``, merge the selected segments
            into one Record. Default is ``False``.
        enforce : bool, optional
            When ``merge=True``, allow integer-sample gaps and overlaps
            to be resolved by :meth:`Record.append`. Default is
            ``False``.
        fillvalue : scalar, optional
            Value used to fill gaps when ``merge=True`` and
            ``enforce=True``.
        snap_tolerance : float, optional
            Times closer than this to a sample, as a fraction of one
            sampling interval, coincide with it. Default is
            :data:`SNAP_TOLERANCE`.
        merge_tolerance : float, optional
            Tolerance for joining the extracted segments when
            ``merge=True``. Default is :data:`MERGE_TOLERANCE`.

        Returns
        -------
        Stream, Record or None
            A new Stream of independent segments when ``merge=False``.
            A merged Record when ``merge=True`` and at least one record
            has samples in the window; ``None`` when ``merge=True`` and
            no record has samples in the window.

            When ``eid`` is provided, an extracted Record is returned.

        Raises
        ------
        TypeError
            If ``merge`` or ``enforce`` is not a boolean.
        KeyError
            If ``eid`` is missing or ambiguous.
        ValueError
            If selected records cannot be merged under the requested
            policy, or if the window is reversed.

        Notes
        -----
        Records without samples in the window are skipped.
        """
        if not isinstance(merge, bool):
            raise TypeError(
                "merge must be a boolean"
            )

        if not isinstance(enforce, bool):
            raise TypeError(
                "enforce must be a boolean"
            )

        if eid is not None:
            record = self.get_record(eid)

            return record.extract(
                starttime=starttime,
                endtime=endtime,
                snap_tolerance=snap_tolerance,
            )

        starttime = self._absolute_time(starttime, "starttime")
        endtime = self._absolute_time(endtime, "endtime")

        if (
            starttime is not None
            and endtime is not None
            and endtime < starttime
        ):
            raise ValueError("endtime precedes starttime")

        selected = []

        for record in self.record:
            if record._sample_window(
                starttime,
                endtime,
                snap_tolerance,
            ) is None:
                continue

            selected.append(
                record.extract(
                    starttime=starttime,
                    endtime=endtime,
                    snap_tolerance=snap_tolerance,
                )
            )

        if not merge:
            stream = Stream(sid=deepcopy(self.sid))
            stream._replace_records(selected)

            return stream

        if not selected:
            return None

        output = selected[0]

        for record in selected[1:]:
            merged = output.append(
                record,
                enforce=enforce,
                fillvalue=fillvalue,
                merge_tolerance=merge_tolerance,
            )

            if not merged:
                raise ValueError(
                    "selected records are not contiguous; use "
                    "enforce=True or merge=False"
                )

        return output

    def _absolute_time(self, value, name):
        """
        Convert a window bound to an absolute Date.

        Numerical values are seconds from the start of the stream;
        ``None`` is returned unchanged.
        """
        if value is None or isinstance(value, Date):
            return value

        if isinstance(value, str):
            return Date(value)

        if isinstance(value, Real) and not isinstance(
            value, (bool, np.bool_)
        ):
            value = float(value)

            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")

            if self.starttime is None:
                raise ValueError(
                    "relative times require a stream with records"
                )

            return self.starttime + value

        raise TypeError(
            f"{name} must be a Date, string, or relative time in "
            "seconds"
        )

    # ------------------------------------------------------------------------
    # Processing applied to every record
    # ------------------------------------------------------------------------

    def remove_mean(self) -> None:
        """
        Remove the mean from every record (see :meth:`Record.remove_mean`).
        """
        for record in self.record:
            record.remove_mean()

    def filter(
        self,
        highpass: Optional[float] = None,
        lowpass: Optional[float] = None,
        order: int = 4,
        minphase: bool = False,
    ) -> None:
        """
        Filter every record (see :meth:`Record.filter`).
        """
        for record in self.record:
            record.filter(
                highpass=highpass,
                lowpass=lowpass,
                order=order,
                minphase=minphase,
            )

    def resample(
        self,
        new_delta: Real,
        method: str = "interp",
        antialias: bool = True,
    ) -> None:
        """
        Resample every record (see :meth:`Record.resample`).
        """
        for record in self.record:
            record.resample(
                new_delta,
                method=method,
                antialias=antialias,
            )

    def taper(self, time: float = 0.1) -> None:
        """
        Taper every record (see :meth:`Record.taper`).
        """
        for record in self.record:
            record.taper(time=time)

    def shift(self, time: float, padding: bool = True) -> None:
        """
        Shift the samples of every record (see :meth:`Record.shift`).
        """
        for record in self.record:
            record.shift(time, padding=padding)

    def integrate(self, method: str = "fft") -> None:
        """
        Integrate every record (see :meth:`Record.integrate`).
        """
        for record in self.record:
            record.integrate(method=method)

    def differentiate(self, method: str = "fft") -> None:
        """
        Differentiate every record (see :meth:`Record.differentiate`).
        """
        for record in self.record:
            record.differentiate(method=method)

    def convolve_response(self, resp: Any, stages: Any = None) -> None:
        """
        Apply an instrumental response to all records.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response object or response container.
        stages
            See :meth:`Record.convolve_response`.

        Notes
        -----
        The operation is performed in place. Response selection is
        delegated to :meth:`Record.convolve_response`.
        """
        for record in self.record:
            record.convolve_response(resp, stages=stages)

    def deconvolve_response(
        self,
        resp: Any,
        waterlevel: float = 60.0,
        method: str = "water",
        stages: Any = None,
    ) -> None:
        """
        Remove an instrumental response from all records.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response object or response container.
        waterlevel, method, stages
            See :meth:`Record.deconvolve_response`.

        Notes
        -----
        The operation is performed in place. Response selection is
        delegated to :meth:`Record.deconvolve_response`.
        """
        for record in self.record:
            record.deconvolve_response(
                resp, waterlevel=waterlevel, method=method, stages=stages
            )

    def remove_sensitivity(self, resp: Any) -> None:
        """
        Divide all records by the overall sensitivity, in place.

        See :meth:`Record.remove_sensitivity`.
        """
        for record in self.record:
            record.remove_sensitivity(resp)


class StreamCollection:
    """
    Collection of waveform streams.

    A StreamCollection groups :class:`Stream` objects by source
    identifier. Depending on ``max_duration``, newly created streams are
    either standard Stream objects or time-limited CircularStream
    objects.

    Parameters
    ----------
    max_duration : float, optional
        Maximum duration, in seconds, retained by each CircularStream.
        If ``None``, standard Stream objects are used.
    """

    def __init__(self, max_duration: Optional[Real] = None) -> None:
        if max_duration is not None:
            if (
                isinstance(max_duration, bool)
                or not isinstance(max_duration, Real)
            ):
                raise TypeError(
                    "max_duration must be a real number or None"
                )

            max_duration = float(max_duration)

            if not np.isfinite(max_duration):
                raise ValueError(
                    "max_duration must be finite"
                )

            if max_duration <= 0.0:
                raise ValueError(
                    "max_duration must be strictly positive"
                )

        self.max_duration = max_duration
        self.stream = []

    def __len__(self) -> int:
        """
        Return the number of streams.
        """
        return len(self.stream)

    def __iter__(self) -> Iterator[Stream]:
        """
        Iterate over streams in insertion order.
        """
        return iter(self.stream)

    def __getitem__(self, key: Any) -> Any:
        """
        Return one stream, or a collection of streams.

        Parameters
        ----------
        key : int, slice or str
            An integer selects a stream by position and a string by
            source identifier. A slice returns a new
            :class:`StreamCollection` with the selected streams.

        Returns
        -------
        Stream or StreamCollection
            Selected stream, or collection of the selected streams.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        KeyError
            If a source identifier is not found or is ambiguous.

        Notes
        -----
        As for Python lists, a slice contains the stored streams, not
        copies. Use :meth:`StreamCollection.copy` on the result for an
        independent selection.
        """
        if isinstance(key, str):
            return self.get_stream(key)

        if isinstance(key, slice):
            collection = self.__class__(max_duration=self.max_duration)
            collection.stream = self.stream[key]
            return collection

        if isinstance(key, Integral) and not isinstance(
            key, (bool, np.bool_)
        ):
            return self.stream[key]

        raise TypeError(
            "collection indices must be integers, slices, or "
            "source identifiers"
        )

    def __repr__(self) -> str:
        return f"StreamCollection(streams={len(self.stream)})"

    def __str__(self) -> str:
        if not self.stream:
            return "empty stream collection"

        lines = [
            "number of streams: {0}".format(len(self)),
        ]

        for index, stream in enumerate(self.stream):
            description = str(stream).rstrip()

            lines.append(
                "Stream {0}: {1}".format(
                    index,
                    description,
                )
            )

        return "\n".join(lines)

    @property
    def sid(self) -> list:
        """
        Source identifiers associated with the streams.

        Returns
        -------
        list
            Source identifiers in stream order.
        """
        return [
            stream.sid
            for stream in self.stream
        ]

    @property
    def nsamp(self) -> int:
        """
        Total number of samples stored in the collection.

        Returns
        -------
        int
            Sum of the sample counts of all streams.
        """
        return sum(
            stream.nsamp
            for stream in self.stream
        )

    def _stream_indices(self, sid):
        """
        Return all stream indices matching a source identifier.

        Parameters
        ----------
        sid : str
            Source identifier to find.

        Returns
        -------
        list of int
            Matching stream indices.
        """
        return [
            index
            for index, stream in enumerate(self.stream)
            if stream.sid == sid
        ]

    def get_stream(self, sid: str) -> Stream:
        """
        Return the stream associated with a source identifier.

        Parameters
        ----------
        sid : str
            Source identifier.

        Returns
        -------
        Stream
            Matching stream.

        Raises
        ------
        TypeError
            If ``sid`` is not a string.
        KeyError
            If the source identifier is not found or is ambiguous.
        """
        if not isinstance(sid, str):
            raise TypeError(
                "sid must be a string"
            )

        indices = self._stream_indices(sid)

        if not indices:
            raise KeyError(
                "source identifier not found: {0}".format(sid)
            )

        if len(indices) > 1:
            raise KeyError(
                "source identifier is not unique: {0}".format(sid)
            )

        return self.stream[indices[0]]

    def _new_stream(self, sid):
        """
        Create a stream suitable for this collection.

        Parameters
        ----------
        sid : str
            Source identifier assigned to the new stream.

        Returns
        -------
        Stream
            Standard or circular stream, depending on the collection
            configuration.
        """
        if self.max_duration is None:
            return Stream(sid=sid)

        # Local import prevents a circular dependency because
        # circular.py imports Stream from this module.
        from shakelab.signals.circular import CircularStream

        return CircularStream(
            sid=sid,
            max_duration=self.max_duration,
        )

    def _get_or_create_stream(self, sid):
        """
        Return an existing stream or create a new one.

        Parameters
        ----------
        sid : str
            Source identifier.

        Returns
        -------
        Stream
            Existing or newly created stream.
        """
        indices = self._stream_indices(sid)

        if len(indices) > 1:
            raise KeyError(
                "source identifier is not unique: {0}".format(sid)
            )

        if indices:
            return self.stream[indices[0]]

        stream = self._new_stream(sid)
        self.stream.append(stream)

        return stream

    def info(self) -> None:
        """
        Print a human-readable description of the collection.
        """
        print(self)

    def copy(self) -> StreamCollection:
        """
        Return an independent copy of the collection.

        Returns
        -------
        StreamCollection
            Collection containing independent copies of all streams and
            records.
        """
        collection = self.__class__(
            max_duration=self.max_duration,
        )

        collection.stream = [
            stream.copy()
            for stream in self.stream
        ]

        return collection

    def append(
        self,
        data: Any,
        enforce: bool = False,
        fillvalue: float = 0.0,
        merge_tolerance: Optional[float] = None,
    ) -> None:
        """
        Append a record or stream to the collection.

        Parameters
        ----------
        data : Record or Stream
            Record or stream to append.
        enforce : bool, optional
            Policy passed to the destination stream append method.
            Default is ``False``.
        fillvalue : scalar, optional
            Value used to fill gaps by standard Stream objects when
            ``enforce=True``.
        merge_tolerance : float, optional
            Tolerance for joining contiguous records, as a fraction of
            the sampling interval. Default is :data:`MERGE_TOLERANCE`.

        Raises
        ------
        TypeError
            If ``data`` is neither a Record nor a Stream, or if
            ``enforce`` is not a boolean.
        ValueError
            If a source identifier is missing or incompatible.

        Notes
        -----
        The destination stream stores copies of the appended records,
        so the caller's records are never modified.
        """
        if not isinstance(data, (Record, Stream)):
            raise TypeError(
                "data must be a Record or Stream instance"
            )

        if not isinstance(enforce, bool):
            raise TypeError(
                "enforce must be a boolean"
            )

        sid = data.sid

        if sid is None:
            raise ValueError(
                "source identifier is required to append data"
            )

        target = self._get_or_create_stream(sid)

        if isinstance(data, Record):
            target.append(
                data,
                enforce=enforce,
                fillvalue=fillvalue,
                merge_tolerance=merge_tolerance,
            )
            return

        for record in data:
            if (
                record.sid is not None
                and record.sid != sid
            ):
                raise ValueError(
                    "record source identifier does not match its stream"
                )

            target.append(
                record,
                enforce=enforce,
                fillvalue=fillvalue,
                merge_tolerance=merge_tolerance,
            )

    def remove(self, key: Any) -> Stream:
        """
        Remove a stream from the collection.

        Parameters
        ----------
        key : int or str
            Integer index or source identifier.

        Returns
        -------
        Stream
            Removed stream.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        IndexError
            If the integer index is out of range.
        KeyError
            If the source identifier is not found or is ambiguous.
        """
        if isinstance(key, str):
            indices = self._stream_indices(key)

            if not indices:
                raise KeyError(
                    "source identifier not found: {0}".format(key)
                )

            if len(indices) > 1:
                raise KeyError(
                    "source identifier is not unique: {0}".format(key)
                )

            key = indices[0]

        elif isinstance(key, (bool, np.bool_)) or not isinstance(
            key, Integral
        ):
            raise TypeError(
                "stream key must be an integer or source identifier"
            )

        return self.stream.pop(key)

    def sort(self, sort_records: bool = True, reverse: bool = False) -> None:
        """
        Sort streams by source identifier.

        Parameters
        ----------
        sort_records : bool, optional
            If ``True``, also sort records within each stream. Default is
            ``True``.
        reverse : bool, optional
            If ``True``, sort streams from descending to ascending source
            identifier. The same ordering is requested for records.
            Default is ``False``.

        Raises
        ------
        TypeError
            If either argument is not a boolean.

        Notes
        -----
        CircularStream objects cannot be reverse-sorted internally.
        Therefore, using both ``sort_records=True`` and ``reverse=True``
        on a circular collection raises the error defined by
        CircularStream.sort().
        """
        if not isinstance(sort_records, bool):
            raise TypeError(
                "sort_records must be a boolean"
            )

        if not isinstance(reverse, bool):
            raise TypeError(
                "reverse must be a boolean"
            )

        def stream_key(stream):
            if stream.sid is None:
                return ""

            return str(stream.sid)

        self.stream.sort(
            key=stream_key,
            reverse=reverse,
        )

        if sort_records:
            for stream in self.stream:
                stream.sort(reverse=reverse)

    def extract(
        self,
        sid: Optional[str] = None,
        starttime: Any = None,
        endtime: Any = None,
        eid: Optional[str] = None,
        merge: bool = False,
        enforce: bool = False,
        fillvalue: float = 0.0,
        snap_tolerance: float = SNAP_TOLERANCE,
        merge_tolerance: Optional[float] = None,
    ) -> Any:
        """
        Extract waveform data from one stream or from all streams.

        Parameters
        ----------
        sid : str, optional
            Source identifier of the stream to extract. If omitted, the
            same extraction is applied to every stream.
        starttime : Date, str or float, optional
            Start of the requested time window. Numerical values are
            seconds from the start of each stream.
        endtime : Date, str or float, optional
            End of the requested time window, included.
        eid : str, optional
            Restrict extraction to one event identifier. Without
            ``sid``, streams without that identifier are skipped.
        merge : bool, optional
            If ``True``, merge selected records into one Record per
            stream. Default is ``False``.
        enforce : bool, optional
            Allow gaps and overlaps to be resolved when merging.
        fillvalue : scalar, optional
            Value used to fill gaps when merging with
            ``enforce=True``.
        snap_tolerance : float, optional
            Times closer than this to a sample, as a fraction of one
            sampling interval, coincide with it. Default is
            :data:`SNAP_TOLERANCE`.
        merge_tolerance : float, optional
            Tolerance for joining segments when ``merge=True``. Default
            is :data:`MERGE_TOLERANCE`.

        Returns
        -------
        Stream, Record, StreamCollection or None
            With ``sid``, the result of :meth:`Stream.extract`. Without
            ``sid``, a new collection (of standard streams) holding the
            extracted data of every stream that has samples in the
            window.
        """
        options = dict(
            starttime=starttime,
            endtime=endtime,
            eid=eid,
            merge=merge,
            enforce=enforce,
            fillvalue=fillvalue,
            snap_tolerance=snap_tolerance,
            merge_tolerance=merge_tolerance,
        )

        if sid is not None:
            return self.get_stream(sid).extract(**options)

        collection = StreamCollection()

        for stream in self.stream:
            if eid is not None and not stream._record_indices(eid):
                continue

            result = stream.extract(**options)

            if result is None or len(result) == 0:
                continue

            collection.append(result)

        return collection

    def merge(
        self,
        collection: StreamCollection,
        enforce: bool = False,
        fillvalue: float = 0.0,
        merge_tolerance: Optional[float] = None,
    ) -> None:
        """
        Append all streams from another collection.

        Parameters
        ----------
        collection : StreamCollection
            Collection to merge into the current collection.
        enforce : bool, optional
            Append policy passed to destination streams.
        fillvalue : scalar, optional
            Gap fill value used by standard streams.
        merge_tolerance : float, optional
            Tolerance for joining contiguous records, as a fraction of
            the sampling interval. Default is :data:`MERGE_TOLERANCE`.

        Raises
        ------
        TypeError
            If ``collection`` is not a StreamCollection.

        Notes
        -----
        Streams and records are not copied automatically.
        """
        if not isinstance(collection, StreamCollection):
            raise TypeError(
                "collection must be a StreamCollection instance"
            )

        for stream in collection:
            self.append(
                stream,
                enforce=enforce,
                fillvalue=fillvalue,
                merge_tolerance=merge_tolerance,
            )

    # ------------------------------------------------------------------------
    # Processing applied to every stream
    # ------------------------------------------------------------------------

    def remove_mean(self) -> None:
        """
        Remove the mean from every record (see :meth:`Record.remove_mean`).
        """
        for stream in self.stream:
            stream.remove_mean()

    def filter(
        self,
        highpass: Optional[float] = None,
        lowpass: Optional[float] = None,
        order: int = 4,
        minphase: bool = False,
    ) -> None:
        """
        Filter every record (see :meth:`Record.filter`).
        """
        for stream in self.stream:
            stream.filter(
                highpass=highpass,
                lowpass=lowpass,
                order=order,
                minphase=minphase,
            )

    def resample(
        self,
        new_delta: Real,
        method: str = "interp",
        antialias: bool = True,
    ) -> None:
        """
        Resample every record (see :meth:`Record.resample`).
        """
        for stream in self.stream:
            stream.resample(
                new_delta,
                method=method,
                antialias=antialias,
            )

    def taper(self, time: float = 0.1) -> None:
        """
        Taper every record (see :meth:`Record.taper`).
        """
        for stream in self.stream:
            stream.taper(time=time)

    def shift(self, time: float, padding: bool = True) -> None:
        """
        Shift the samples of every record (see :meth:`Record.shift`).
        """
        for stream in self.stream:
            stream.shift(time, padding=padding)

    def integrate(self, method: str = "fft") -> None:
        """
        Integrate every record (see :meth:`Record.integrate`).
        """
        for stream in self.stream:
            stream.integrate(method=method)

    def differentiate(self, method: str = "fft") -> None:
        """
        Differentiate every record (see :meth:`Record.differentiate`).
        """
        for stream in self.stream:
            stream.differentiate(method=method)

    def convolve_response(self, resp: Any, stages: Any = None) -> None:
        """
        Apply an instrumental response to all streams.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response container passed to each stream.
        stages
            See :meth:`Record.convolve_response`.

        Notes
        -----
        The operation is performed in place.
        """
        for stream in self.stream:
            stream.convolve_response(resp, stages=stages)

    def deconvolve_response(
        self,
        resp: Any,
        waterlevel: float = 60.0,
        method: str = "water",
        stages: Any = None,
    ) -> None:
        """
        Remove an instrumental response from all streams.

        Parameters
        ----------
        resp : StageResponse, StageSet, ChannelEpoch or MetadataCollection
            Response container passed to each stream.
        waterlevel, method, stages
            See :meth:`Record.deconvolve_response`.

        Notes
        -----
        The operation is performed in place.
        """
        for stream in self.stream:
            stream.deconvolve_response(
                resp, waterlevel=waterlevel, method=method, stages=stages
            )

    def remove_sensitivity(self, resp: Any) -> None:
        """
        Divide all streams by the overall sensitivity, in place.

        See :meth:`Record.remove_sensitivity`.
        """
        for stream in self.stream:
            stream.remove_sensitivity(resp)

    def __mul__(self, value: Any) -> StreamCollection:
        """
        Return a response-convolved copy of the collection.

        ``value`` may be any response container accepted by
        :meth:`Record.convolve_response`.
        """
        if not isinstance(value, _response_types()):
            raise TypeError(
                "unsupported operand type for StreamCollection "
                "multiplication: {0}".format(
                    type(value).__name__
                )
            )

        collection = self.copy()
        collection.convolve_response(value)

        return collection

    def __truediv__(self, value: Any) -> StreamCollection:
        """
        Return a response-deconvolved copy of the collection.

        ``value`` may be any response container accepted by
        :meth:`Record.deconvolve_response`.
        """
        if not isinstance(value, _response_types()):
            raise TypeError(
                "unsupported operand type for StreamCollection "
                "division: {0}".format(
                    type(value).__name__
                )
            )

        collection = self.copy()
        collection.deconvolve_response(value)

        return collection

    def read(
        self,
        file_path: Any = None,
        list_file: Any = None,
        format: Optional[str] = None,
        byte_order: Optional[str] = None,
        **format_options: Any,
    ) -> StreamCollection:
        """
        Read waveform data into the collection.

        Parameters
        ----------
        file_path : str, path-like or sequence, optional
            Waveform file path, wildcard pattern, directory, or sequence of
            waveform paths. Mutually exclusive with ``list_file``.
        list_file : str or path-like, optional
            Text file containing one waveform path per line. Mutually
            exclusive with ``file_path``. Relative paths are resolved
            against the directory containing the list file.
        format : str, optional
            Input waveform format. If omitted, the format is inferred
            independently from the extension of each input file.
        byte_order : str or None, optional
            Byte order used by formats supporting explicit selection.
            ``None`` delegates byte-order handling to the selected
            format-specific reader.
        **format_options
            Additional keyword arguments forwarded to
            :func:`shakelab.signals.io.reader`.

            For MiniSEED input, ``mseed_backend`` may be:

            - ``"auto"`` to use the configured default backend;
            - ``"libmseed"`` to use the libmseed-backed implementation;
            - ``"python"`` to use the pure-Python implementation.

            Other options are forwarded to the selected format-specific
            reader.

        Returns
        -------
        StreamCollection
            The current collection, updated with the loaded waveform data.

        Raises
        ------
        TypeError
            If an input source or a format-specific option has an
            unsupported type.
        ValueError
            If neither input source is provided, both sources are provided,
            or the waveform format cannot be determined.
        NotImplementedError
            If the selected waveform format is known but not implemented.

        Notes
        -----
        Loaded records are appended to the current collection. Existing
        streams and records are preserved.

        Input expansion, format detection, MiniSEED backend selection,
        byte-order handling, and final stream sorting are delegated entirely
        to :func:`shakelab.signals.io.reader`.

        Examples
        --------
        Read a single file using format inference:

        >>> collection.read("waveforms.mseed")

        Read all supported files in a directory:

        >>> collection.read("/path/to/waveforms")

        Read paths listed in a text file:

        >>> collection.read(list_file="waveforms.txt")

        Force the pure-Python MiniSEED backend:

        >>> collection.read(
        ...     "waveforms.mseed",
        ...     mseed_backend="python",
        ... )
        """
        io.reader(
            file_path=file_path,
            list_file=list_file,
            stream_collection=self,
            format=format,
            byte_order=byte_order,
            **format_options,
        )

        return self

    def write(
        self,
        file_path: Any,
        format: Optional[str] = None,
        byte_order: Optional[str] = None,
        owrite: bool = False,
        **format_options: Any,
    ) -> None:
        """
        Write the collection to disk.

        Parameters
        ----------
        file_path : str or path-like
            Output file path for single-file formats, or output directory for
            multi-file formats such as SAC.
        format : str, optional
            Output waveform format. If omitted, the format is inferred from
            ``file_path``.
        byte_order : str or None, optional
            Byte order used by formats supporting explicit selection.
            ``None`` delegates byte-order handling to the selected
            format-specific writer.
        owrite : bool, optional
            If ``True``, existing output files may be overwritten. Default
            is ``False``.
        **format_options
            Additional keyword arguments forwarded to
            :func:`shakelab.signals.io.writer`.
    
            For MiniSEED output, ``mseed_backend`` may be:
    
            - ``"auto"`` to use the configured default backend;
            - ``"libmseed"`` to use the libmseed-backed implementation;
            - ``"python"`` to use the pure-Python implementation.
    
            MiniSEED-specific options such as ``encoding`` and ``reclen`` may
            also be supplied when supported by the selected backend.

        Raises
        ------
        TypeError
            If a format-specific option has an unsupported type.
        ValueError
            If the output format cannot be inferred or is not recognized.
        FileExistsError
            If an output file already exists and ``owrite`` is ``False``.
        NotImplementedError
            If the selected output format is known but not implemented.

        Notes
        -----
        Format selection, output validation, MiniSEED backend selection, and
        conversion to the representation required by each backend are
        delegated entirely to :func:`shakelab.signals.io.writer`.

        MiniSEED is written to a single output file. SAC output uses
        ``file_path`` as a directory and creates one SAC file for each
        record.

        Examples
        --------
        Write MiniSEED using format inference:

        >>> collection.write("waveforms.mseed")

        Write MiniSEED with the pure-Python backend:

        >>> collection.write(
        ...     "waveforms.mseed",
        ...     mseed_backend="python",
        ...     encoding=11,
        ...     reclen=4096,
        ... )

        Write one SAC file per record:

        >>> collection.write(
        ...     "sac_output",
        ...     format="sac",
        ...     byte_order="be",
        ...     owrite=True,
        ... )
        """
        io.writer(
            input_data=self,
            file_path=file_path,
            format=format,
            byte_order=byte_order,
            owrite=owrite,
            **format_options,
        )


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def tukey_window(N: int, alpha: float = 0.5) -> np.ndarray:
    """
    Return a Tukey (tapered cosine) window.

    Parameters
    ----------
    N : int
        Number of points.
    alpha : float, optional
        Fraction of the window inside the cosine tapered region; 0
        gives a rectangular window, 1 a Hann window.

    Returns
    -------
    numpy.ndarray
        Window values; equivalent to
        :func:`scipy.signal.windows.tukey`.
    """
    return windows.tukey(N, alpha)


def _response_types():
    """Return the response classes accepted by Record operators."""
    return (
        response.StageResponse,
        response.StageSet,
        metadata.ChannelEpoch,
        metadata.MetadataCollection,
    )


def _fft_derivative(data, order, period):
    """
    Differentiate (order 1) or integrate (order -1) periodic data.

    The zero-frequency component is removed and, for an even number of
    samples, the Nyquist component is set to zero, as in the former
    ``scipy.fftpack.diff``. Complex data are processed part by part.
    """
    data = np.asarray(data)

    if np.iscomplexobj(data):
        return (
            _fft_derivative(data.real, order, period)
            + 1j * _fft_derivative(data.imag, order, period)
        )

    nsamp = data.size
    spectrum = np.fft.rfft(data)
    omega = 2.0 * np.pi * np.arange(spectrum.size) / period

    factor = np.zeros(spectrum.size, dtype=complex)
    factor[1:] = (1j * omega[1:]) ** order

    if nsamp % 2 == 0:
        factor[-1] = 0.0

    return np.fft.irfft(spectrum * factor, nsamp)


def _fft_resample(data, old_delta, new_delta, new_nsamp):
    """
    Resample data in the Fourier domain to a new sampling interval.

    Parameters
    ----------
    data : numpy.ndarray
        Samples at interval ``old_delta``.
    old_delta, new_delta : float
        Old and new sampling intervals.
    new_nsamp : int
        Number of output samples to return.

    Returns
    -------
    numpy.ndarray
        Resampled data.

    Raises
    ------
    ValueError
        If the interval ratio is not a simple fraction.
    """
    ratio = Fraction(new_delta / old_delta).limit_denominator(1000)

    if not np.isclose(
        float(ratio),
        new_delta / old_delta,
        rtol=1e-9,
        atol=0.0,
    ):
        raise ValueError(
            "Fourier resampling requires the ratio of the sampling "
            "intervals to be a simple fraction; use method='interp'"
        )

    # new_delta / old_delta = down / up: a block of `down` old samples
    # spans `up` new samples.
    down = ratio.numerator
    up = ratio.denominator

    nsamp = data.size
    padded = -(-nsamp // down) * down

    if padded > nsamp:
        data = np.concatenate(
            (data, np.zeros(padded - nsamp, dtype=data.dtype))
        )

    resampled = signal.resample(data, padded // down * up)

    return resampled[:new_nsamp]
