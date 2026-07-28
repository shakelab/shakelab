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
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Core waveform classes and processing utilities.

This module provides the fundamental data structures used throughout
ShakeLab for uniformly sampled time-series analysis.

The main classes are:

- :class:`Header`, storing waveform metadata.
- :class:`Record`, representing a continuous waveform segment together
  with its metadata.
- :class:`Stream`, representing a collection of temporally ordered
  records belonging to the same source.
- :class:`CircularStream`, implementing a time-limited circular buffer.
- :class:`StreamCollection`, grouping multiple streams.

The module also provides basic signal-processing utilities including
filtering, resampling, tapering, Fourier-domain operations, temporal
windowing, convolution, correlation, and response correction.

Waveform objects are designed to be lightweight, largely independent of
external libraries, and suitable for scientific processing pipelines.
"""

from copy import deepcopy
from numbers import Number, Real

import numpy as np
from scipy import fftpack, integrate, signal

from shakelab.signals import fourier
from shakelab.signals import response
from shakelab.signals import io
from shakelab.libutils.timeN import Date
from shakelab.libutils.geodetic import WgsPoint


class Header:
    """
    Metadata associated with a waveform record.

    A Header stores the sampling interval, sampling rate, start time,
    source identifier, event identifier, physical units, location,
    response information, and arbitrary metadata associated with a
    :class:`Record`.

    The sampling interval and sampling rate are kept synchronized:

    ``rate = 1 / delta``

    Parameters
    ----------
    delta : float, optional
        Sampling interval in seconds. It must be finite and strictly
        positive.
    time : Date or str, optional
        Start time of the associated record. Strings are converted to
        :class:`Date`.
    location : WgsPoint or sequence of float, optional
        Geographic location. A two-element sequence is interpreted as
        longitude and latitude and converted to :class:`WgsPoint`.
    sid : str, optional
        Source identifier, normally expressed as an FDSN source code.
    eid : str, optional
        Event identifier.
    units : str, optional
        Physical units of the waveform samples.
    parent : object, optional
        Object owning the header, normally a :class:`Record`.

    Notes
    -----
    The parent reference is intentionally excluded from :meth:`copy`
    unless explicitly provided. This prevents copying a Header from
    recursively duplicating its owning Record.
    """

    def __init__(
        self,
        delta=None,
        time=None,
        location=None,
        sid=None,
        eid=None,
        units=None,
        parent=None,
    ):
        self._delta = None
        self._rate = None
        self._time = None
        self._location = None
        self._parent = parent

        self.sid = sid
        self.eid = eid
        self.units = units

        self.response = None
        self.meta = {}

        if delta is not None:
            self.delta = delta

        if time is not None:
            self.time = time

        if location is not None:
            self.location = location

    def __str__(self):
        lines = [
            "sampling interval: {0}".format(self.delta),
            "sampling rate: {0}".format(self.rate),
            "starting time: {0}".format(self.time),
            "source identifier: {0}".format(self.sid),
            "event identifier: {0}".format(self.eid),
            "location: {0}".format(self.location),
            "units: {0}".format(self.units),
        ]

        return "\n".join(lines)

    @staticmethod
    def _validate_sampling_value(value, name):
        """
        Validate a sampling interval or sampling rate.

        Parameters
        ----------
        value : Real
            Value to validate.
        name : str
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

        return value

    @property
    def delta(self):
        """
        Sampling interval in seconds.
        """
        return self._delta

    @delta.setter
    def delta(self, value):
        if value is None:
            self._delta = None
            self._rate = None
            return

        value = self._validate_sampling_value(
            value,
            "delta",
        )

        self._delta = value
        self._rate = 1.0 / value

    @property
    def rate(self):
        """
        Sampling rate in samples per second.
        """
        return self._rate

    @rate.setter
    def rate(self, value):
        if value is None:
            self._rate = None
            self._delta = None
            return

        value = self._validate_sampling_value(
            value,
            "rate",
        )

        self._rate = value
        self._delta = 1.0 / value

    @property
    def time(self):
        """
        Start time of the associated record.
        """
        return self._time

    @time.setter
    def time(self, value):
        if value is None:
            self._time = None
            return

        if isinstance(value, Date):
            self._time = value
            return

        if isinstance(value, str):
            self._time = Date(value)
            return

        raise TypeError(
            "time must be a Date instance or a string"
        )

    @property
    def location(self):
        """
        Geographic location of the associated record.
        """
        return self._location

    @location.setter
    def location(self, value):
        if value is None:
            self._location = None
            return

        if isinstance(value, WgsPoint):
            self._location = value
            return

        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise ValueError(
                    "location must contain longitude and latitude"
                )

            self._location = WgsPoint(
                value[0],
                value[1],
            )
            return

        raise TypeError(
            "location must be a WgsPoint or a two-element sequence"
        )

    @property
    def nsamp(self):
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

    def info(self):
        """
        Print a human-readable description of the header.
        """
        print(self)

    def copy(self, parent=None):
        """
        Return an independent copy of the header.

        Parameters
        ----------
        parent : object, optional
            Parent to associate with the copied header. By default, the
            copied header has no parent.

        Returns
        -------
        Header
            Independent header copy.

        Notes
        -----
        The parent reference is not copied automatically because it
        normally points back to the owning Record. Copying that reference
        with ``deepcopy`` could duplicate the complete object graph or
        leave the copied Header associated with the wrong Record.
        """
        header = Header(
            delta=self.delta,
            time=deepcopy(self.time),
            location=deepcopy(self.location),
            sid=deepcopy(self.sid),
            eid=deepcopy(self.eid),
            units=deepcopy(self.units),
            parent=parent,
        )

        header.response = deepcopy(self.response)
        header.meta = deepcopy(self.meta)

        return header


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
    time : Date or str, optional
        Start time of the record.
    location : WgsPoint or sequence of float, optional
        Geographic location associated with the record.
    sid : str, optional
        Source identifier.
    eid : str, optional
        Event identifier.

    Notes
    -----
    Numerical operators return independent Record objects and never
    modify their operands. Signal-processing methods generally operate
    in place unless explicitly documented otherwise.
    """

    def __init__(
        self,
        data=None,
        delta=None,
        time=None,
        location=None,
        sid=None,
        eid=None,
    ):
        self.head = Header(
            delta=delta,
            time=time,
            location=location,
            sid=sid,
            eid=eid,
            parent=self,
        )

        self._data = np.array([], dtype=float)

        if data is not None:
            self.data = data

    def __len__(self):
        """
        Return the number of samples.
        """
        return self.nsamp

    def __getitem__(self, item):
        """
        Return one or more waveform samples.
        """
        return self.data[item]

    def __str__(self):
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
    def data(self):
        """
        One-dimensional waveform data.
        """
        return self._data

    @data.setter
    def data(self, value):
        array = np.asarray(value)

        if array.ndim == 0:
            array = array.reshape(1)

        if array.ndim != 1:
            raise ValueError(
                "record data must be one-dimensional"
            )

        self._data = array.copy()

    @property
    def sid(self):
        """
        Source identifier.
        """
        return self.head.sid

    @sid.setter
    def sid(self, value):
        self.head.sid = value

    @property
    def eid(self):
        """
        Event identifier.
        """
        return self.head.eid

    @eid.setter
    def eid(self, value):
        self.head.eid = value

    @property
    def nsamp(self):
        """
        Number of waveform samples.
        """
        return self.data.size

    @property
    def delta(self):
        """
        Sampling interval in seconds.
        """
        return self.head.delta

    @delta.setter
    def delta(self, value):
        self.head.delta = value

    @property
    def rate(self):
        """
        Sampling rate in samples per second.
        """
        return self.head.rate

    @rate.setter
    def rate(self, value):
        self.head.rate = value

    @property
    def time(self):
        """
        Start time of the record.
        """
        return self.head.time

    @time.setter
    def time(self, value):
        self.head.time = value

    @property
    def peak(self):
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

    def copy(self):
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

    def info(self):
        """
        Print a human-readable description of the record.
        """
        print(self)

    def __add__(self, value):
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

    def __radd__(self, value):
        """
        Return the result of scalar addition.
        """
        return self.__add__(value)

    def __sub__(self, value):
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

    def __rsub__(self, value):
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

    def __mul__(self, value):
        """
        Return a scaled or response-convolved record.

        Parameters
        ----------
        value : Number or StageResponse
            Numerical scale factor or instrumental response.

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

        if isinstance(value, response.StageResponse):
            record.convolve_response(value)
            return record

        raise TypeError(
            "unsupported operand type for Record multiplication: "
            "{0}".format(type(value).__name__)
        )

    def __rmul__(self, value):
        """
        Return the result of scalar multiplication.
        """
        return self.__mul__(value)

    def __truediv__(self, value):
        """
        Return a scaled or response-deconvolved record.

        Parameters
        ----------
        value : Number or StageResponse
            Numerical divisor or instrumental response.

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

        if isinstance(value, response.StageResponse):
            record.deconvolve_response(value)
            return record

        raise TypeError(
            "unsupported operand type for Record division: {0}".format(
                type(value).__name__
            )
        )

    @property
    def starttime(self):
        """
        Time of the first sample.

        Returns
        -------
        Date or None
            Record start time.
        """
        return self.head.time

    @property
    def endtime(self):
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
    def duration(self):
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

    def _sample_position(self, time, tolerance=None):
        """
        Convert a relative time to a floating-point sample position.

        Parameters
        ----------
        time : float
            Time in seconds relative to record start.
        tolerance : float, optional
            Maximum accepted distance from the nearest sample, expressed
            as a fraction of the sampling interval. The default is
            ``1e-6``.

        Returns
        -------
        int
            Index of the nearest sample.

        Raises
        ------
        ValueError
            If the requested time is not aligned with the sampling grid.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required for index conversion"
            )

        if tolerance is None:
            tolerance = 1e-6

        position = float(time) / self.delta
        index = int(round(position))

        if not np.isclose(
            position,
            index,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "time does not coincide with the sampling grid"
            )

        return index

    def overlaps(self, starttime=None, endtime=None):
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

    def time_axis(self, reference="relative", shift=0.0):
        """
        Generate the record time axis.

        Parameters
        ----------
        reference : {"relative", "absolute", "r", "a"}, optional
            Time reference. Relative time starts from zero. Absolute time
            is expressed as seconds in the internal Date reference system.
        shift : float, optional
            Additional shift in seconds.

        Returns
        -------
        numpy.ndarray
            Time coordinate for each waveform sample.

        Raises
        ------
        ValueError
            If the reference is unsupported or absolute time is requested
            without a start time.
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

        axis = np.arange(
            self.nsamp,
            dtype=float,
        ) * self.delta

        if reference in ("relative", "r"):
            return axis + shift

        if reference in ("absolute", "a"):
            if self.time is None:
                raise ValueError(
                    "start time is required for an absolute time axis"
                )

            return axis + self.time.to_seconds() + shift

        raise ValueError(
            "reference must be 'relative' or 'absolute'"
        )

    @property
    def taxis(self):
        """
        Relative time axis in seconds.
        """
        return self.time_axis()

    def append(
        self,
        record,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
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
        tolerance : float, optional
            Maximum distance from an integer sample offset, expressed as
            a fraction of one sampling interval.

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

        if (
            not isinstance(tolerance, Real)
            or isinstance(tolerance, (bool, np.bool_))
        ):
            raise TypeError("tolerance must be a real number")

        tolerance = float(tolerance)

        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError(
                "tolerance must be finite and non-negative"
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
            atol=tolerance,
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

    def cut(
        self,
        starttime=None,
        endtime=None,
        inplace=True,
        tolerance=1e-6,
    ):
        """
        Cut the record to a sample-aligned time window.

        Parameters
        ----------
        starttime : Date, str or float, optional
            First sample to retain. Numerical values are interpreted as
            seconds from record start.
        endtime : Date, str or float, optional
            Final sample to retain. Numerical values are interpreted as
            seconds from record start.
        inplace : bool, optional
            If ``True``, modify the current record. If ``False``, return
            an independent Record. Default is ``True``.
        tolerance : float, optional
            Maximum distance from the nearest sampling-grid point,
            expressed as a fraction of one sampling interval. Default is
            ``1e-6``.

        Returns
        -------
        Record or None
            When ``inplace=False``, return the cut record. When
            ``inplace=True``, return ``None``.

        Raises
        ------
        TypeError
            If ``inplace`` is not a boolean or ``tolerance`` is not a
            real number.
        ValueError
            If the tolerance is invalid, the sampling interval is
            missing, the requested interval is invalid, does not
            intersect the record, or is not aligned with the sampling
            grid.
        """
        if not isinstance(inplace, bool):
            raise TypeError(
                "inplace must be a boolean"
            )

        if (
            not isinstance(tolerance, Real)
            or isinstance(tolerance, (bool, np.bool_))
        ):
            raise TypeError(
                "tolerance must be a real number"
            )

        tolerance = float(tolerance)

        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError(
                "tolerance must be finite and non-negative"
            )

        if self.nsamp == 0:
            if inplace:
                return None

            return self.copy()

        if self.delta is None:
            raise ValueError(
                "sampling interval is required to cut a record"
            )

        relative_start, relative_end = self._relative_time_window(
            starttime=starttime,
            endtime=endtime,
        )

        if (
            relative_end < 0.0
            or relative_start > self.duration
        ):
            raise ValueError(
                "requested time window does not intersect the record"
            )

        relative_start = max(
            relative_start,
            0.0,
        )
        relative_end = min(
            relative_end,
            self.duration,
        )

        index_start = self._sample_position(
            relative_start,
            tolerance=tolerance,
        )
        index_end = self._sample_position(
            relative_end,
            tolerance=tolerance,
        )

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
        starttime=None,
        endtime=None,
        tolerance=1e-6,
    ):
        """
        Return an independent record for a sample-aligned time window.

        Parameters
        ----------
        starttime : Date, str or float, optional
            First sample to retain.
        endtime : Date, str or float, optional
            Final sample to retain.
        tolerance : float, optional
            Maximum distance from the nearest sampling-grid point,
            expressed as a fraction of one sampling interval.

        Returns
        -------
        Record
            Extracted record.
        """
        return self.cut(
            starttime=starttime,
            endtime=endtime,
            inplace=False,
            tolerance=tolerance,
        )

    def remove_mean(self):
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
        highpass=None,
        lowpass=None,
        order=4,
        minphase=False,
    ):
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

        if not isinstance(order, int) or isinstance(order, bool):
            raise TypeError("order must be an integer")

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
        new_delta,
        method="interp",
        antialias=True,
    ):
        """
        Resample the waveform to a new sampling interval.

        Parameters
        ----------
        new_delta : float
            New sampling interval in seconds.
        method : {"scipy", "fft", "interp"}, optional
            Resampling method:

            - ``"scipy"`` or ``"fft"`` uses Fourier-domain resampling
              through :func:`scipy.signal.resample`.
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

        The number of output samples is chosen to preserve the original
        first-to-last-sample duration as closely as possible:

        ``new_nsamp = round(duration / new_delta) + 1``

        If the original duration is not an exact multiple of
        ``new_delta``, the final sample time changes by at most half a new
        sampling interval.
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

        new_nsamp = max(
            1,
            int(round(old_duration / new_delta)) + 1,
        )

        if method in ("scipy", "fft"):
            self.data = signal.resample(
                self.data,
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

            # The rounding used to determine new_nsamp can place the
            # final target point marginally beyond the original record.
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

    def taper(self, time=0.1):
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

    def zero_padding(self, time=None):
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
                atol=1e-6,
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

    def shift(self, time, padding=True):
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
        original_is_real = np.isrealobj(self.data)

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

        shifted = fourier.shift_time(
            data,
            self.delta,
            time,
        )

        shifted = shifted[:original_nsamp]

        if original_is_real:
            shifted = np.real_if_close(
                shifted,
                tol=1000,
            )

        self.data = shifted

    def integrate(self, method="fft"):
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

        self.data = fftpack.diff(
            self.data,
            order=-1,
            period=self.nsamp * self.delta,
        )

    def differentiate(self, method="fft"):
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

        self.data = fftpack.diff(
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

    def to_spectrum(self):
        """
        Convert the record to a frequency-domain Spectrum.

        Returns
        -------
        Spectrum
            Spectrum generated from the current record.

        Notes
        -----
        The exact Record-Spectrum metadata contract will be reviewed in
        the dedicated refactoring of the Fourier module.
        """
        if self.delta is None:
            raise ValueError(
                "sampling interval is required to compute a spectrum"
            )

        return fourier.Spectrum(self)

    def from_spectrum(self, spectrum):
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
        data,
        mode="full",
        method="auto",
    ):
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

    def deconvolve(self, data):
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
        record,
        mode="full",
        method="auto",
    ):
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

    def _resolve_response(self, resp):
        """
        Resolve a response container for the current stream and time.

        Parameters
        ----------
        resp : response object
            Stage response, stage set, stream response, or response
            collection.

        Returns
        -------
        StageResponse or StageSet
            Response object directly applicable to the current Record.

        Raises
        ------
        TypeError
            If the response type is unsupported.
        ValueError
            If stream id or time metadata required for response selection
            is missing.
        KeyError
            If the stream or epoch is not present in the response
            container.
        """
        if isinstance(
            resp,
            (
                response.StageResponse,
                response.StageSet,
            ),
        ):
            return resp

        if isinstance(resp, response.StreamResponse):
            if self.time is None:
                raise ValueError(
                    "record start time is required to select a "
                    "response epoch"
                )

            try:
                return resp[self.time]
            except (KeyError, IndexError) as error:
                raise KeyError(
                    "no response epoch is available for the record "
                    "start time"
                ) from error

        if isinstance(resp, response.ResponseCollection):
            if self.sid is None:
                raise ValueError(
                    "source identifier is required to select a response"
                )

            if self.sid not in resp.sid:
                raise KeyError(
                    "source identifier not found in response "
                    "collection: {0}".format(self.sid)
                )

            stream_response = resp[self.sid]

            return self._resolve_response(stream_response)

        raise TypeError(
            "unsupported response object: {0}".format(
                type(resp).__name__
            )
        )

    def convolve_response(self, resp):
        """
        Apply an instrumental response to the waveform.

        Parameters
        ----------
        resp : StageResponse, StageSet, StreamResponse or ResponseCollection
            Response object or response container.

        Raises
        ------
        TypeError
            If the response object is unsupported or does not return a
            Record.
        KeyError
            If no response is available for the record source or epoch.

        Notes
        -----
        The operation is performed in place. The complete Record returned
        by the response implementation is adopted, including any updated
        units or response metadata.
        """
        resolved = self._resolve_response(resp)

        convolve_record = getattr(
            resolved,
            "convolve_record",
            None,
        )

        if not callable(convolve_record):
            raise TypeError(
                "resolved response does not expose convolve_record()"
            )

        corrected = convolve_record(self)

        self._replace_from_record(corrected)

    def deconvolve_response(self, resp):
        """
        Remove an instrumental response from the waveform.

        Parameters
        ----------
        resp : StageResponse, StageSet, StreamResponse or ResponseCollection
            Response object or response container.

        Raises
        ------
        TypeError
            If the response object is unsupported or does not return a
            Record.
        KeyError
            If no response is available for the record source or epoch.

        Notes
        -----
        The operation is performed in place. The complete Record returned
        by the response implementation is adopted, including any updated
        units or response metadata.
        """
        resolved = self._resolve_response(resp)

        deconvolve_record = getattr(
            resolved,
            "deconvolve_record",
            None,
        )

        if not callable(deconvolve_record):
            raise TypeError(
                "resolved response does not expose "
                "deconvolve_record()"
            )

        corrected = deconvolve_record(self)

        self._replace_from_record(corrected)


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

    def __init__(self, sid=None):
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

    def __len__(self):
        """
        Return the number of records.
        """
        return len(self.record)

    def __iter__(self):
        """
        Iterate over the records in insertion order.
        """
        return iter(self.record)

    def __getitem__(self, key):
        """
        Return one or more records.

        Parameters
        ----------
        key : int, slice or str
            Integer and slice keys access records by position. A string
            key accesses a record by event identifier.

        Returns
        -------
        Record or list of Record
            Selected record or records.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        KeyError
            If a requested event identifier is not found or is
            ambiguous.
        """
        if isinstance(key, str):
            return self.get_record(key)

        if isinstance(key, (int, slice)):
            if isinstance(key, int):
                return self.record[key]

            if isinstance(key, slice):
                return list(self.record)[key]

        raise TypeError(
            "stream indices must be integers, slices, or event "
            "identifiers"
        )

    def __str__(self):
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
    def eid(self):
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
    def nsamp(self):
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
    def starttime(self):
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
    def endtime(self):
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

    def get_record(self, eid):
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

    def info(self):
        """
        Print a human-readable description of the stream.
        """
        print(self)

    def copy(self):
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
        record,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
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
        tolerance : float, optional
            Maximum distance from an integer sample offset, expressed as
            a fraction of one sampling interval.

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
        Records are preserved in insertion order.

        The new record is merged only with the final record currently
        stored in the stream. To guarantee meaningful temporal merging,
        the stream should be sorted before appending out-of-order data.
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

        if record.sid is None and self.sid is not None:
            record.sid = self.sid

        if not self.record:
            self.record.append(record)
            return

        merged = self.record[-1].append(
            record,
            enforce=enforce,
            fillvalue=fillvalue,
            tolerance=tolerance,
        )

        if not merged:
            self.record.append(record)

    def remove(self, key):
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

        elif isinstance(key, bool) or not isinstance(key, int):
            raise TypeError(
                "record key must be an integer or an event identifier"
            )

        records = list(self.record)
        removed = records.pop(key)
        self._replace_records(records)

        return removed

    def sort(self, reverse=False):
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
        starttime=None,
        endtime=None,
        eid=None,
        merge=False,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
        """
        Extract records intersecting a time window.

        Parameters
        ----------
        starttime : Date, str or float, optional
            Start of the requested time window.
        endtime : Date, str or float, optional
            End of the requested time window.
        eid : str, optional
            Restrict extraction to the record associated with this event
            identifier.
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
        tolerance : float, optional
            Sampling-grid tolerance passed to Record temporal methods.

        Returns
        -------
        Stream, Record or None
            A new Stream when ``merge=False``. A merged Record when
            ``merge=True`` and at least one record intersects the
            requested window. ``None`` is returned when ``merge=True``
            and no record intersects the window.

            When ``eid`` is provided, an extracted Record is returned.

        Raises
        ------
        TypeError
            If ``merge`` or ``enforce`` is not a boolean.
        KeyError
            If ``eid`` is missing or ambiguous.
        ValueError
            If selected records cannot be merged under the requested
            policy.
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
                tolerance=tolerance,
            )

        selected = []

        for record in self.record:
            if record.overlaps(
                starttime=starttime,
                endtime=endtime,
            ):
                selected.append(
                    record.extract(
                        starttime=starttime,
                        endtime=endtime,
                        tolerance=tolerance,
                    )
                )

        if not merge:
            stream = Stream(
                sid=deepcopy(self.sid),
            )
            stream.record = selected

            return stream

        if not selected:
            return None

        output = selected[0]

        for record in selected[1:]:
            merged = output.append(
                record,
                enforce=enforce,
                fillvalue=fillvalue,
                tolerance=tolerance,
            )

            if not merged:
                raise ValueError(
                    "selected records are not contiguous; use "
                    "enforce=True or merge=False"
                )

        return output

    def convolve_response(self, resp):
        """
        Apply an instrumental response to all records.

        Parameters
        ----------
        resp : StageResponse, StageSet, StreamResponse or ResponseCollection
            Response object or response container.

        Notes
        -----
        The operation is performed in place. Response selection is
        delegated to :meth:`Record.convolve_response`.
        """
        for record in self.record:
            record.convolve_response(resp)

    def deconvolve_response(self, resp):
        """
        Remove an instrumental response from all records.

        Parameters
        ----------
        resp : StageResponse, StageSet, StreamResponse or ResponseCollection
            Response object or response container.

        Notes
        -----
        The operation is performed in place. Response selection is
        delegated to :meth:`Record.deconvolve_response`.
        """
        for record in self.record:
            record.deconvolve_response(resp)


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

    def __init__(self, max_duration=None):
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

    def __len__(self):
        """
        Return the number of streams.
        """
        return len(self.stream)

    def __iter__(self):
        """
        Iterate over streams in insertion order.
        """
        return iter(self.stream)

    def __getitem__(self, key):
        """
        Return one or more streams.

        Parameters
        ----------
        key : int, slice or str
            Integer and slice keys access streams by position. A string
            key accesses a stream by source identifier.

        Returns
        -------
        Stream or list of Stream
            Selected stream or streams.

        Raises
        ------
        TypeError
            If the key type is unsupported.
        KeyError
            If a source identifier is not found or is ambiguous.
        """
        if isinstance(key, str):
            return self.get_stream(key)

        if isinstance(key, bool):
            raise TypeError(
                "collection indices must be integers, slices, or "
                "source identifiers"
            )

        if isinstance(key, int):
            return self.stream[key]

        if isinstance(key, slice):
            return self.stream[key]

        raise TypeError(
            "collection indices must be integers, slices, or "
            "source identifiers"
        )

    def __str__(self):
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
    def sid(self):
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
    def nsamp(self):
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

    def get_stream(self, sid):
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

    def info(self):
        """
        Print a human-readable description of the collection.
        """
        print(self)

    def copy(self):
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
        data,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
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
        tolerance : float, optional
            Sampling-grid tolerance passed to the destination stream.

        Raises
        ------
        TypeError
            If ``data`` is neither a Record nor a Stream, or if
            ``enforce`` is not a boolean.
        ValueError
            If a source identifier is missing or incompatible.

        Notes
        -----
        Appended records are not copied automatically. Use
        :meth:`copy` before appending when independent objects are
        required.
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
                tolerance=tolerance,
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
                tolerance=tolerance,
            )

    def remove(self, key):
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

        elif isinstance(key, bool) or not isinstance(key, int):
            raise TypeError(
                "stream key must be an integer or source identifier"
            )

        return self.stream.pop(key)

    def sort(
        self,
        sort_records=True,
        reverse=False,
    ):
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
        sid,
        starttime=None,
        endtime=None,
        eid=None,
        merge=False,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
        """
        Extract waveform data from one stream.

        Parameters
        ----------
        sid : str
            Source identifier of the stream to extract.
        starttime : Date, str or float, optional
            Start of the requested time window.
        endtime : Date, str or float, optional
            End of the requested time window.
        eid : str, optional
            Restrict extraction to one event identifier.
        merge : bool, optional
            If ``True``, merge selected records into one Record.
            Default is ``False``.
        enforce : bool, optional
            Allow gaps and overlaps to be resolved when merging.
        fillvalue : scalar, optional
            Value used to fill gaps when merging with
            ``enforce=True``.
        tolerance : float, optional
            Sampling-grid tolerance.

        Returns
        -------
        Stream, Record or None
            Result returned by :meth:`Stream.extract`.
        """
        stream = self.get_stream(sid)

        return stream.extract(
            starttime=starttime,
            endtime=endtime,
            eid=eid,
            merge=merge,
            enforce=enforce,
            fillvalue=fillvalue,
            tolerance=tolerance,
        )

    def merge(
        self,
        collection,
        enforce=False,
        fillvalue=0.0,
        tolerance=1e-6,
    ):
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
        tolerance : float, optional
            Sampling-grid tolerance.

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
                tolerance=tolerance,
            )

    def convolve_response(self, resp):
        """
        Apply an instrumental response to all streams.

        Parameters
        ----------
        resp : ResponseCollection
            Response container passed to each stream.

        Notes
        -----
        The operation is performed in place.
        """
        for stream in self.stream:
            stream.convolve_response(resp)

    def deconvolve_response(self, resp):
        """
        Remove an instrumental response from all streams.

        Parameters
        ----------
        resp : ResponseCollection
            Response container passed to each stream.

        Notes
        -----
        The operation is performed in place.
        """
        for stream in self.stream:
            stream.deconvolve_response(resp)

    def __mul__(self, value):
        """
        Return a response-convolved copy of the collection.
        """
        if not isinstance(value, response.ResponseCollection):
            raise TypeError(
                "unsupported operand type for StreamCollection "
                "multiplication: {0}".format(
                    type(value).__name__
                )
            )

        collection = self.copy()
        collection.convolve_response(value)

        return collection

    def __truediv__(self, value):
        """
        Return a response-deconvolved copy of the collection.
        """
        if not isinstance(value, response.ResponseCollection):
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
        file_path=None,
        list_file=None,
        format=None,
        byte_order=None,
        **format_options,
    ):
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
        file_path,
        format=None,
        byte_order=None,
        owrite=False,
        **format_options,
    ):
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


def tukey_window(N, alpha=0.5):
    """
    Replicates scipy.signal.tukey(N, alpha)
    """
    if alpha <= 0:
        return np.ones(N)
    elif alpha >= 1:
        return np.hanning(N)

    x = np.linspace(0, 1, N)
    w = np.ones(x.shape)

    # first condition 0 <= x < alpha/2
    mask = x < alpha / 2
    w[mask] = 0.5 * (1 + np.cos(
        2 * np.pi / alpha * (x[mask] - alpha / 2)))

    # second condition 1 - alpha/2 < x <= 1
    mask = x > (1 - alpha / 2)
    w[mask] = 0.5 * (1 + np.cos(
        2 * np.pi / alpha * (x[mask] - 1 + alpha / 2)))

    return w
