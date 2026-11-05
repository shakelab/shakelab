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
TDMS reader.

TDMS is the binary format of National Instruments (LabVIEW, DAQ
systems), also used by DAS interrogators such as the Silixa iDAS. A
file is a sequence of segments, each with a lead-in, the description
of the objects (file, groups, channels) and their properties, and the
raw data of the channels, contiguous or interleaved (NI, "TDMS file
format internal structure").

Reading
-------
- Segments that reuse or extend the object list of the previous one,
  repeated chunks, contiguous and interleaved data in both byte
  orders; numeric channels (integers, floats, also with units,
  complex, booleans), timestamps (``datetime64[ns]``, UTC) and text.
- A last segment with less data than declared (an interrupted
  acquisition) gives a warning and the data that are present, as
  npTDMS.
- NI scales (``NI_Scale[n]_*``: linear, polynomial, sums and
  differences of scales) are applied unless ``scaled=False``; other
  scales (thermocouples, strain gauges, ...) give a warning and the raw
  values.
- DAQmx raw data and multi-dimensional channels are not supported.
- Files are mapped in memory, so the memory used is about the size of
  the data read.

ShakeLab conversion
-------------------
:func:`tdmsread` reads the numeric channels of one group as records.
The sampling rate, start time and units are taken from the channel
(LabVIEW waveform properties ``wf_increment``, ``wf_start_time``,
``wf_start_offset``, ``unit_string``), then from the group and the
file (e.g. ``SamplingFrequency[Hz]`` and ``GPSTimeStamp`` of iDAS
files), then from a file name with ``_UTC_YYYYMMDD_HHMMSS.fff``. The
properties of the channel are kept in ``record.head.extra["tdms"]``.

:func:`tdms_to_matrix` builds a DAS matrix from one or more contiguous
files. It is provisional: the DAS containers of ShakeLab are still to
be designed.
"""

from __future__ import annotations

import glob
import mmap
import os
import re
import struct
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import StreamCollection

__all__ = [
    "DasData",
    "TdmsChannel",
    "TdmsFile",
    "TdmsGroup",
    "TdmsWarning",
    "is_tdms",
    "tdms_to_matrix",
    "tdmsread",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

TDS_TAG = b"TDSm"
LEAD_IN_SIZE = 28

TOC_METADATA = 1 << 1
TOC_NEW_OBJECT_LIST = 1 << 2
TOC_RAW_DATA = 1 << 3
TOC_INTERLEAVED = 1 << 5
TOC_BIG_ENDIAN = 1 << 6
TOC_DAQMX_RAW_DATA = 1 << 7

RAW_NO_DATA = 0xFFFFFFFF
RAW_MATCHES_PREVIOUS = 0x00000000
RAW_DAQMX = (0x00001269, 0x00001369)
UINT64_MAX = 0xFFFFFFFFFFFFFFFF

TYPE_STRING = 0x20
TYPE_TIMESTAMP = 0x44

# Seconds from the LabVIEW epoch (1904-01-01) to the Unix epoch.
LABVIEW_EPOCH_OFFSET = 2082844800

TDMS_NUMPY_TYPES = {
    0x01: np.dtype("i1"),
    0x02: np.dtype("i2"),
    0x03: np.dtype("i4"),
    0x04: np.dtype("i8"),
    0x05: np.dtype("u1"),
    0x06: np.dtype("u2"),
    0x07: np.dtype("u4"),
    0x08: np.dtype("u8"),
    0x09: np.dtype("f4"),
    0x0A: np.dtype("f8"),
    0x19: np.dtype("f4"),
    0x1A: np.dtype("f8"),
    0x21: np.dtype("?"),
    0x08000C: np.dtype("c8"),
    0x10000D: np.dtype("c16"),
}

# Property names, in order of preference.
RATE_PROPERTIES = (
    "SamplingFrequency[Hz]", "sampling_rate", "Sampling Rate",
    "SampleRate", "sample_rate", "fs",
)
START_PROPERTIES = (
    "ISO8601 Timestamp", "GPSTimeStamp", "wf_start_time", "start_time",
    "StartTime",
)
SPACING_PROPERTIES = (
    "SpatialResolution[m]", "SpatialResolution", "channel_spacing",
    "Channel Spacing",
)
GAUGE_PROPERTIES = (
    "GaugeLength", "GaugeLength[m]", "gauge_length", "Gauge Length",
)
UNITS_PROPERTIES = ("unit_string", "NI_UnitDescription")

RAW_INPUT_SOURCE = 0xFFFFFFFF
_SCALE_TYPE = re.compile(r"NI_Scale\[(\d+)\]_Scale_Type")


class TdmsWarning(UserWarning):
    """Warning about TDMS data that are incomplete or not scaled."""


# ----------------------------------------------------------------------------
# TDMS objects
# ----------------------------------------------------------------------------


@dataclass
class TdmsChannel:
    """
    TDMS channel.

    Attributes
    ----------
    name : str
        Channel name.
    path : str
        TDMS object path (``/'group'/'channel'``).
    properties : dict
        Channel properties; timestamps as :class:`Date`.
    data : numpy.ndarray
        Values: numbers (scaled if requested), ``datetime64[ns]`` or
        strings.
    """

    name: str
    path: str
    properties: dict
    data: np.ndarray


@dataclass
class TdmsGroup:
    """
    TDMS group.

    Attributes
    ----------
    name : str
        Group name.
    properties : dict
        Group properties.
    channels : dict
        Channels by name, in file order.
    """

    name: str
    properties: dict = field(default_factory=dict)
    channels: dict = field(default_factory=dict)


@dataclass
class TdmsFile:
    """
    Content of a TDMS file.

    Attributes
    ----------
    path : str or None
        File path, None if read from memory.
    properties : dict
        File properties.
    groups : dict
        Groups by name, in file order.
    """

    path: Optional[str] = None
    properties: dict = field(default_factory=dict)
    groups: dict = field(default_factory=dict)

    @classmethod
    def read(cls, source: Any, scaled: bool = True) -> TdmsFile:
        """
        Read a TDMS file.

        Parameters
        ----------
        source : str, path-like, bytes or binary file object
            TDMS file, its content, or a file opened in binary mode.
        scaled : bool, optional
            If True (default), apply the NI scales of the channels.

        Raises
        ------
        ValueError
            If the file is not a valid TDMS file.
        NotImplementedError
            If the file holds DAQmx raw data or multi-dimensional
            channels.
        """
        if isinstance(source, (str, os.PathLike)):
            path = os.fspath(source)
            with open(path, "rb") as fid:
                if os.fstat(fid.fileno()).st_size == 0:
                    raise ValueError(f"empty file: {path}")
                buffer = mmap.mmap(fid.fileno(), 0, access=mmap.ACCESS_READ)
            try:
                return cls.from_bytes(buffer, scaled, path=path)
            finally:
                try:
                    buffer.close()
                except BufferError:
                    # Views held by an exception being raised: the map
                    # is closed when they are released.
                    pass
        if isinstance(source, (bytes, bytearray, memoryview)):
            return cls.from_bytes(source, scaled)
        if hasattr(source, "read"):
            content = source.read()
            if not isinstance(content, (bytes, bytearray)):
                raise TypeError("the file must be opened in binary mode")
            return cls.from_bytes(content, scaled)
        raise TypeError(
            "source must be a path, bytes or a binary file object"
        )

    @classmethod
    def from_bytes(
        cls,
        buffer: Any,
        scaled: bool = True,
        path: Optional[str] = None,
    ) -> TdmsFile:
        """
        Decode the content of a TDMS file.

        Parameters
        ----------
        buffer : bytes-like
            Content of the file.
        scaled : bool, optional
            If True (default), apply the NI scales of the channels.
        path : str, optional
            File path, kept in :attr:`path`.

        Raises
        ------
        ValueError
            If the content is not a valid TDMS file.
        """
        properties, data = _parse(memoryview(buffer).cast("B"))
        tdms = cls(path=path)
        tdms.properties = properties.get("/", {})
        for obj_path, props in properties.items():
            parts = _parse_path(obj_path)
            if len(parts) == 1:
                tdms._group(parts[0]).properties = props
        for obj_path in _channel_paths(properties, data):
            group_name, channel_name = _parse_path(obj_path)
            group = tdms._group(group_name)
            props = properties.get(obj_path, {})
            values = data.get(obj_path)
            values = (
                np.array([]) if not values
                else values[0] if len(values) == 1
                else np.concatenate(values)
            )
            if scaled and values.dtype.kind in "iufb":
                values = _apply_scaling(
                    values, obj_path, props, group.properties,
                    tdms.properties,
                )
            group.channels[channel_name] = TdmsChannel(
                channel_name, obj_path, props, values
            )
        return tdms

    def info(self) -> None:
        """Print the groups and channels of the file."""
        print(f"TDMS file: {self.path}")
        print(f"File properties: {len(self.properties)}")
        for group in self.groups.values():
            channels = list(group.channels.values())
            print(f"Group {group.name!r}: {len(group.properties)} "
                  f"properties, {len(channels)} channels")
            for channel in channels[:10]:
                print(f"  {channel.name}: {channel.data.size} values, "
                      f"{channel.data.dtype}")
            if len(channels) > 10:
                print(f"  ... ({len(channels) - 10} more)")

    def _group(self, name: str) -> TdmsGroup:
        """Group of a name, created if needed."""
        if name not in self.groups:
            self.groups[name] = TdmsGroup(name)
        return self.groups[name]


@dataclass
class DasData:
    """
    DAS matrix (provisional, until the DAS containers are designed).

    Attributes
    ----------
    data : numpy.ndarray
        Values with shape (channels, samples).
    time : numpy.ndarray or None
        Times of the samples from the start time, in seconds.
    distance : numpy.ndarray or None
        Distances of the channels along the fibre, in metres.
    sampling_rate, channel_spacing, gauge_length : float or None
        Hz and metres.
    start_time : Date or None
        Time of the first sample.
    properties : dict
        File and group properties of the first file, with the group
        name, channel names and source files.
    """

    data: np.ndarray
    time: Optional[np.ndarray]
    distance: Optional[np.ndarray]
    sampling_rate: Optional[float]
    channel_spacing: Optional[float]
    gauge_length: Optional[float]
    start_time: Optional[Date]
    properties: dict


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def tdmsread(
    source: Any,
    stream_collection: Optional[StreamCollection] = None,
    group: Optional[str] = None,
    scaled: bool = True,
    rate: Optional[float] = None,
    start: Any = None,
    network: str = "",
    station: str = "",
    location: str = "",
    channel_prefix: str = "",
    dtype: Any = None,
) -> StreamCollection:
    """
    Read the numeric channels of a TDMS group as records.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        TDMS file, its content, or a file opened in binary mode.
    stream_collection : StreamCollection, optional
        Collection to which the records are appended; a new one if
        omitted.
    group : str, optional
        Group to read. If None, the only group, or the group most
        likely to hold DAS data (numeric channels of equal length).
    scaled : bool, optional
        If True (default), apply the NI scales of the channels.
    rate : float, optional
        Sampling rate in Hz, for all channels; from the properties if
        None.
    start : Date or str, optional
        Time of the first sample, for all channels; from the
        properties or the file name if None.
    network, station, location : str, optional
        Codes of the source identifiers; the channel code is
        ``channel_prefix`` followed by the channel name.
    channel_prefix : str, optional
        Prefix of the channel codes.
    dtype : data-type, optional
        Type of the samples; the type of the file if None.

    Returns
    -------
    StreamCollection
        One record per numeric channel.

    Raises
    ------
    KeyError
        If the requested group does not exist.
    ValueError
        If no group can be chosen or the sampling rate of a channel is
        unknown.
    """
    from shakelab.signals.base import Record, StreamCollection

    tdms = TdmsFile.read(source, scaled)
    selected = _select_group(tdms, group)
    file_props = {**tdms.properties, **selected.properties}

    if stream_collection is None:
        stream_collection = StreamCollection()
    index = {stream.sid: stream for stream in stream_collection.stream}
    for name in sorted(selected.channels, key=_channel_sort_key):
        channel = selected.channels[name]
        if channel.data.dtype.kind not in "iufbc" or not channel.data.size:
            continue
        props = channel.properties
        channel_rate = rate or _channel_rate(props, file_props)
        if channel_rate is None:
            raise ValueError(
                f"the sampling rate of channel {name!r} is unknown: give "
                "rate"
            )
        channel_start = start if start is not None else _channel_start(
            props, file_props, tdms.path
        )
        data = channel.data if dtype is None else channel.data.astype(dtype)
        record = Record(
            data=data,
            delta=1.0 / float(channel_rate),
            time=_to_date(channel_start),
            sid=f"{network}.{station}.{location}.{channel_prefix}{name}",
            units=_get_property(props, UNITS_PROPERTIES),
        )
        record.head.extra["tdms"] = {
            "group": selected.name,
            "channel": name,
            "properties": dict(props),
        }
        _append_record(stream_collection, record, index)
    return stream_collection


def tdms_to_matrix(
    source: Any,
    group: Optional[str] = None,
    scaled: bool = True,
    rate: Optional[float] = None,
    channel_spacing: Optional[float] = None,
    gauge_length: Optional[float] = None,
    start: Any = None,
) -> DasData:
    """
    Build a DAS matrix from one or more TDMS files (provisional).

    Parameters
    ----------
    source : str, path-like or sequence
        TDMS file, wildcard pattern, or sequence of paths and patterns;
        several files are taken in alphabetical order and joined in
        time.
    group : str, optional
        Group to read in each file; chosen as in :func:`tdmsread` if
        None.
    scaled : bool, optional
        If True (default), apply the NI scales of the channels.
    rate : float, optional
        Sampling rate in Hz; from the properties if None.
    channel_spacing : float, optional
        Channel spacing in metres; from the properties if None (for
        iDAS files multiplied by ``Fibre Length Multiplier``).
    gauge_length : float, optional
        Gauge length in metres; from the properties if None.
    start : Date or str, optional
        Time of the first sample of the first file; from the
        properties or the file name if None.

    Returns
    -------
    DasData

    Raises
    ------
    ValueError
        If the files differ in channels, type, sampling, spacing or
        gauge length, or are not contiguous in time.
    """
    paths = _resolve_paths(source)
    blocks = []
    reference = None
    previous_end = None

    for k, path in enumerate(paths):
        tdms = TdmsFile.read(path, scaled)
        selected = _select_group(tdms, group)
        names = sorted(selected.channels, key=_channel_sort_key)
        if not names:
            raise ValueError(f"the selected group has no channels: {path}")
        arrays = [selected.channels[name].data for name in names]
        if len({a.size for a in arrays}) != 1:
            raise ValueError(f"DAS channels of different length: {path}")
        block = np.vstack(arrays)
        props = {**tdms.properties, **selected.properties}

        file_rate = rate or _get_property(props, RATE_PROPERTIES)
        if file_rate is None and "wf_increment" in props:
            file_rate = 1.0 / float(props["wf_increment"])
        if file_rate is None:
            raise ValueError(f"the sampling rate is unknown: {path}")
        layout = {
            "names": names,
            "dtype": block.dtype,
            "rate": float(file_rate),
            "spacing": _float_or_none(
                channel_spacing if channel_spacing is not None
                else _channel_spacing(props)
            ),
            "gauge": _float_or_none(
                gauge_length if gauge_length is not None
                else _get_property(props, GAUGE_PROPERTIES)
            ),
        }
        if k == 0 and start is not None:
            file_start = _to_date(start)
        else:
            file_start = _get_property(props, START_PROPERTIES)
            if file_start is None:
                file_start = _file_name_time(path)
            file_start = _to_date(file_start)

        if k == 0:
            reference = layout
            first = (props, selected.name, file_start)
            if len(paths) > 1 and file_start is None:
                raise ValueError(
                    "the start time is needed to join several files"
                )
        else:
            _check_layout(path, layout, reference)
            if file_start is None:
                raise ValueError(f"the start time is unknown: {path}")
            offset = (file_start.to_unix_ns() - previous_end) / 1e9
            if abs(offset) > max(1e-6, 1e-3 / reference["rate"]):
                raise ValueError(
                    f"the files are not contiguous in time ({offset:+g} "
                    f"s): {path}"
                )
        blocks.append(block)
        previous_end = file_start.to_unix_ns() + round(
            block.shape[1] / reference["rate"] * 1e9
        ) if file_start is not None else None

    data = blocks[0] if len(blocks) == 1 else np.concatenate(blocks, axis=1)
    props, group_name, start_time = first
    properties = dict(props)
    properties["tdms_group"] = group_name
    properties["channel_names"] = reference["names"]
    properties["source_files"] = [str(p) for p in paths]
    spacing = reference["spacing"]
    return DasData(
        data=data,
        time=np.arange(data.shape[1]) / reference["rate"],
        distance=(
            None if spacing is None else np.arange(data.shape[0]) * spacing
        ),
        sampling_rate=reference["rate"],
        channel_spacing=spacing,
        gauge_length=reference["gauge"],
        start_time=start_time,
        properties=properties,
    )


def is_tdms(head: bytes) -> bool:
    """
    Tell whether the beginning of a file is a TDMS segment.

    Parameters
    ----------
    head : bytes
        First bytes of the file (at least 4).

    Returns
    -------
    bool
        True if the file starts with the tag of a TDMS segment
        (``TDSm``); index files (``.tdms_index``, tag ``TDSh``) are not
        data files.

    Examples
    --------
    >>> with open("measurement.tdms", "rb") as fid:
    ...     is_tdms(fid.read(4))
    True
    """
    if not isinstance(head, (bytes, bytearray, memoryview)):
        raise TypeError("head must be bytes")
    return bytes(head[:4]) == TDS_TAG


# ----------------------------------------------------------------------------
# Binary parser
# ----------------------------------------------------------------------------


@dataclass
class _Object:
    """Description of a TDMS object in one segment."""

    path: str
    has_data: bool = False
    data_type: Optional[int] = None
    number_values: int = 0
    data_size: int = 0


class _Cursor:
    """Sequential reader of the metadata of a segment."""

    def __init__(self, buffer: memoryview, offset: int, endian: str,
                 end: int) -> None:
        self.buffer = buffer
        self.offset = offset
        self.endian = endian
        self.end = end

    def unpack(self, fmt: str) -> tuple:
        size = struct.calcsize(fmt)
        if self.offset + size > self.end:
            raise ValueError("TDMS metadata extend beyond the segment")
        values = struct.unpack_from(self.endian + fmt, self.buffer,
                                    self.offset)
        self.offset += size
        return values

    def string(self) -> str:
        (size,) = self.unpack("I")
        if self.offset + size > self.end:
            raise ValueError("TDMS string extends beyond the segment")
        text = bytes(self.buffer[self.offset:self.offset + size])
        self.offset += size
        return text.decode("utf-8", errors="replace").rstrip("\x00")

    def value(self, data_type: int) -> Any:
        if data_type == TYPE_STRING:
            return self.string()
        if data_type == TYPE_TIMESTAMP:
            if self.endian == "<":
                fraction, seconds = self.unpack("Qq")
            else:
                seconds, fraction = self.unpack("qQ")
            return Date.from_unix_ns(_labview_ns(seconds, fraction))
        if data_type not in TDMS_NUMPY_TYPES:
            raise NotImplementedError(
                f"TDMS property type {data_type:#x} is not supported"
            )
        dtype = TDMS_NUMPY_TYPES[data_type].newbyteorder(self.endian)
        if self.offset + dtype.itemsize > self.end:
            raise ValueError("TDMS metadata extend beyond the segment")
        value = np.frombuffer(
            self.buffer, dtype=dtype, count=1, offset=self.offset
        )[0].item()
        self.offset += dtype.itemsize
        return value


def _parse(buffer: memoryview) -> tuple:
    """Properties and data chunks of all the objects of a file."""
    size = len(buffer)
    properties: dict = {}
    data: dict = {}
    known: dict = {}
    previous: Optional[list] = None
    position = 0

    while position < size:
        if position + LEAD_IN_SIZE > size:
            _warn(f"incomplete segment lead-in at byte {position} ignored")
            break
        if bytes(buffer[position:position + 4]) != TDS_TAG:
            raise ValueError(
                f"not a TDMS segment at byte {position}: "
                f"{bytes(buffer[position:position + 4])!r}"
            )
        (toc,) = struct.unpack_from("<I", buffer, position + 4)
        endian = ">" if toc & TOC_BIG_ENDIAN else "<"
        _, next_offset, raw_offset = struct.unpack_from(
            endian + "IQQ", buffer, position + 8
        )
        meta = position + LEAD_IN_SIZE
        data_position = meta + raw_offset
        incomplete = next_offset == UINT64_MAX
        end = size if incomplete else meta + next_offset
        if end > size:
            end, incomplete = size, True
        if data_position > end:
            _warn(f"incomplete segment at byte {position} ignored")
            break
        if toc & TOC_DAQMX_RAW_DATA:
            raise NotImplementedError("DAQmx raw data are not supported")

        if toc & TOC_METADATA:
            objects = _read_metadata(
                _Cursor(buffer, meta, endian, data_position), toc,
                previous, known, properties,
            )
        elif previous is None:
            raise ValueError(
                f"segment at byte {position} has no metadata and no "
                "previous segment"
            )
        else:
            objects = previous
        if toc & TOC_RAW_DATA:
            _read_raw_data(
                buffer, data_position, end, objects, endian,
                bool(toc & TOC_INTERLEAVED), incomplete, data,
            )
        for obj in objects:
            known[obj.path] = obj
        previous = objects
        position = end
    return properties, data


def _read_metadata(
    cursor: _Cursor,
    toc: int,
    previous: Optional[list],
    known: dict,
    properties: dict,
) -> list:
    """Object list of a segment with metadata."""
    if toc & TOC_NEW_OBJECT_LIST or previous is None:
        objects: list = []
    else:
        objects = [_Object(**vars(obj)) for obj in previous]
    positions = {obj.path: k for k, obj in enumerate(objects)}

    (count,) = cursor.unpack("I")
    for _ in range(count):
        path = cursor.string()
        (index,) = cursor.unpack("I")
        if path in positions:
            obj = objects[positions[path]]
        else:
            obj = (
                _Object(**vars(known[path])) if path in known
                else _Object(path)
            )
            positions[path] = len(objects)
            objects.append(obj)

        if index == RAW_NO_DATA:
            obj.has_data = False
        elif index == RAW_MATCHES_PREVIOUS:
            if obj.data_type is None:
                raise ValueError(
                    f"{path} reuses a raw data index never defined"
                )
            obj.has_data = True
        elif index in RAW_DAQMX:
            raise NotImplementedError("DAQmx raw data are not supported")
        else:
            _read_raw_index(cursor, obj, index)

        (nprops,) = cursor.unpack("I")
        if nprops:
            props = properties.setdefault(path, {})
            for _ in range(nprops):
                name = cursor.string()
                (data_type,) = cursor.unpack("I")
                props[name] = cursor.value(data_type)
        else:
            properties.setdefault(path, {})
    return objects


def _read_raw_index(cursor: _Cursor, obj: _Object, length: int) -> None:
    """Raw data index of an object."""
    start = cursor.offset
    data_type, dimension, number_values = cursor.unpack("IIQ")
    if dimension != 1:
        raise NotImplementedError(
            f"{obj.path}: only one-dimensional channels are supported"
        )
    obj.has_data = True
    obj.data_type = data_type
    obj.number_values = int(number_values)
    if data_type == TYPE_STRING:
        (obj.data_size,) = cursor.unpack("Q")
    else:
        obj.data_size = obj.number_values * _item_size(data_type, obj.path)
    cursor.offset = max(cursor.offset, start + length - 4)


def _read_raw_data(
    buffer: memoryview,
    start: int,
    end: int,
    objects: list,
    endian: str,
    interleaved: bool,
    incomplete: bool,
    data: dict,
) -> None:
    """Raw data of a segment, by chunks."""
    channels = [obj for obj in objects if obj.has_data]
    chunk = sum(obj.data_size for obj in channels)
    total = end - start
    if chunk == 0 or total <= 0:
        return
    nfull, remainder = divmod(total, chunk)
    lengths = {obj.path: obj.number_values for obj in channels}
    final = None
    if remainder:
        _warn(f"segment data at byte {start} are {total} bytes, not a "
              f"multiple of the chunk size {chunk}: last chunk partial")
        if any(obj.data_type == TYPE_STRING for obj in channels):
            final = None
        elif interleaved or not incomplete:
            final = {
                obj.path: obj.number_values * remainder // chunk
                for obj in channels
            }
        else:
            final, left = {}, remainder
            for obj in channels:
                size = obj.data_size
                final[obj.path] = (
                    obj.number_values if left > size
                    else left // _item_size(obj.data_type, obj.path)
                )
                left -= min(left, size)

    if interleaved:
        _read_interleaved(buffer, start, channels, endian, nfull, final,
                          data)
        return
    offset = start
    for k in range(nfull + (final is not None)):
        sizes = lengths if k < nfull else final
        for obj in channels:
            count = sizes.get(obj.path, 0)
            values = _decode(buffer, offset, obj, count, endian)
            if values is not None:
                data.setdefault(obj.path, []).append(values)
            offset += (
                obj.data_size if obj.data_type == TYPE_STRING
                else count * _item_size(obj.data_type, obj.path)
            )


def _read_interleaved(
    buffer: memoryview,
    start: int,
    channels: list,
    endian: str,
    nfull: int,
    final: Optional[dict],
    data: dict,
) -> None:
    """Interleaved data: one value of each channel per row."""
    if any(obj.data_type == TYPE_STRING for obj in channels):
        raise ValueError("interleaved text channels are not valid TDMS")
    if len({obj.number_values for obj in channels}) != 1:
        raise NotImplementedError(
            "interleaved channels of different lengths are not supported"
        )
    rows = channels[0].number_values * nfull
    if final:
        rows += min(final.values())
    row = np.dtype([
        (f"f{k}", _raw_dtype(obj.data_type, endian, obj.path))
        for k, obj in enumerate(channels)
    ])
    table = np.frombuffer(buffer, dtype=row, count=rows, offset=start)
    for k, obj in enumerate(channels):
        values = table[f"f{k}"]
        data.setdefault(obj.path, []).append(
            _timestamps(values) if obj.data_type == TYPE_TIMESTAMP
            else values.astype(values.dtype.newbyteorder("="))
        )


def _decode(
    buffer: memoryview,
    offset: int,
    obj: _Object,
    count: int,
    endian: str,
) -> Optional[np.ndarray]:
    """Values of one object in a contiguous chunk."""
    if obj.data_type == TYPE_STRING:
        if count != obj.number_values:
            return None
        ends = np.frombuffer(buffer, dtype=endian + "u4", count=count,
                             offset=offset)
        base = offset + 4 * count
        starts = np.concatenate(([0], ends[:-1]))
        return np.array([
            bytes(buffer[base + a:base + b]).decode("utf-8", "replace")
            for a, b in zip(starts, ends)
        ])
    values = np.frombuffer(
        buffer, dtype=_raw_dtype(obj.data_type, endian, obj.path),
        count=count, offset=offset,
    )
    if obj.data_type == TYPE_TIMESTAMP:
        return _timestamps(values)
    return values.astype(values.dtype.newbyteorder("="))


def _raw_dtype(data_type: int, endian: str, path: str) -> np.dtype:
    """NumPy type of the raw values of a channel."""
    if data_type == TYPE_TIMESTAMP:
        if endian == "<":
            return np.dtype([("fraction", "<u8"), ("seconds", "<i8")])
        return np.dtype([("seconds", ">i8"), ("fraction", ">u8")])
    try:
        return TDMS_NUMPY_TYPES[data_type].newbyteorder(endian)
    except KeyError:
        raise NotImplementedError(
            f"{path}: TDMS data type {data_type:#x} is not supported"
        ) from None


def _item_size(data_type: int, path: str) -> int:
    """Size in bytes of one value."""
    return _raw_dtype(data_type, "<", path).itemsize


def _timestamps(values: np.ndarray) -> np.ndarray:
    """``datetime64[ns]`` of raw TDMS timestamps."""
    seconds = values["seconds"].astype(np.int64) - LABVIEW_EPOCH_OFFSET
    fraction = values["fraction"].astype(np.float64) / 2.0**64
    ns = seconds * 10**9 + np.round(fraction * 1e9).astype(np.int64)
    return ns.astype("datetime64[ns]")


def _labview_ns(seconds: int, fraction: int) -> int:
    """Unix nanoseconds of a TDMS timestamp."""
    return (seconds - LABVIEW_EPOCH_OFFSET) * 10**9 + (
        (fraction * 10**9 + 2**63) >> 64
    )


def _channel_paths(properties: dict, data: dict) -> list:
    """Paths of the channels, in file order."""
    paths = list(dict.fromkeys(list(properties) + list(data)))
    return [p for p in paths if len(_parse_path(p)) == 2]


def _parse_path(path: str) -> list:
    """Components of a TDMS object path (``/'group'/'channel'``)."""
    if path == "/":
        return []
    parts = []
    k = 0
    while k < len(path):
        if path[k] != "/" or k + 1 >= len(path) or path[k + 1] != "'":
            raise ValueError(f"invalid TDMS path: {path}")
        k += 2
        chars = []
        while k < len(path):
            if path[k] == "'":
                if path[k + 1:k + 2] == "'":
                    chars.append("'")
                    k += 2
                    continue
                k += 1
                break
            chars.append(path[k])
            k += 1
        parts.append("".join(chars))
    return parts


# ----------------------------------------------------------------------------
# NI scales
# ----------------------------------------------------------------------------


def _apply_scaling(
    values: np.ndarray,
    path: str,
    channel: dict,
    group: dict,
    file: dict,
) -> np.ndarray:
    """Values with the NI scales of the channel (or group, or file)."""
    for props in (channel, group, file):
        count = _number_of_scales(props)
        if not count:
            continue
        if props.get("NI_Scaling_Status", "unscaled") == "scaled":
            return values
        try:
            return _scale(count - 1, values, props, {})
        except _UnsupportedScale as err:
            _warn(f"{path}: NI scale {err} not supported, raw values kept")
            return values
    return values


class _UnsupportedScale(Exception):
    """A scale type that is not implemented."""


def _number_of_scales(props: dict) -> int:
    """Number of NI scales in a set of properties."""
    if "NI_Number_Of_Scales" in props:
        return int(props["NI_Number_Of_Scales"])
    indices = [
        int(m.group(1)) for m in map(_SCALE_TYPE.match, props) if m
    ]
    return max(indices) + 1 if indices else 0


def _scale(index: int, raw: np.ndarray, props: dict, cache: dict) -> Any:
    """Values of scale ``index`` (raw values for the raw source)."""
    if index == RAW_INPUT_SOURCE:
        return raw
    if index in cache:
        return cache[index]
    prefix = f"NI_Scale[{index}]_"
    kind = props.get(prefix + "Scale_Type")

    def source(name: str = "") -> Any:
        key = f"{prefix}{kind}_{name}Input_Source"
        return _scale(int(props.get(key, RAW_INPUT_SOURCE)), raw, props,
                      cache)

    if kind == "Linear":
        result = (
            source().astype(np.float64) * props[prefix + "Linear_Slope"]
            + props[prefix + "Linear_Y_Intercept"]
        )
    elif kind == "Polynomial":
        size = int(props.get(prefix + "Polynomial_Coefficients_Size", 4))
        coefficients = [
            props[f"{prefix}Polynomial_Coefficients[{k}]"]
            for k in range(size)
        ]
        result = np.polynomial.polynomial.polyval(
            source().astype(np.float64), coefficients
        )
    elif kind == "Add":
        result = source("Left_Operand_") + source("Right_Operand_")
    elif kind == "Subtract":
        result = source("Right_Operand_") - source("Left_Operand_")
    else:
        raise _UnsupportedScale(f"{index} ({kind})")
    cache[index] = result
    return result


# ----------------------------------------------------------------------------
# ShakeLab helpers
# ----------------------------------------------------------------------------


def _select_group(tdms: TdmsFile, name: Optional[str]) -> TdmsGroup:
    """Requested group, the only one, or the most DAS-like one."""
    if name is not None:
        if name not in tdms.groups:
            raise KeyError(f"TDMS group not found: {name!r}")
        return tdms.groups[name]
    if len(tdms.groups) == 1:
        return next(iter(tdms.groups.values()))
    scored = [
        (score, group) for group in tdms.groups.values()
        if (score := _das_score(group)) is not None
    ]
    if not scored:
        raise ValueError(
            "cannot choose a TDMS group: give one of "
            f"{list(tdms.groups)}"
        )
    return max(scored, key=lambda item: item[0])[1]


def _das_score(group: TdmsGroup) -> Optional[float]:
    """Likeness of a group to DAS data (None if not numeric)."""
    channels = list(group.channels.values())
    if not channels:
        return None
    if len({c.data.size for c in channels}) != 1 or not all(
        c.data.dtype.kind in "iufc" for c in channels
    ):
        return None
    numeric = sum(c.name.isdigit() for c in channels)
    score = 10.0 + min(len(channels), 10000) / 100.0
    score += numeric / len(channels)
    if group.name.lower() in ("measurement", "measurements", "data"):
        score += 5.0
    return score


def _channel_rate(props: dict, file_props: dict) -> Optional[float]:
    """Sampling rate of a channel."""
    if props.get("wf_increment"):
        return 1.0 / float(props["wf_increment"])
    value = _get_property(props, RATE_PROPERTIES)
    if value is None:
        value = _get_property(file_props, RATE_PROPERTIES)
    if value is None and file_props.get("wf_increment"):
        value = 1.0 / float(file_props["wf_increment"])
    return None if value is None else float(value)


def _channel_start(
    props: dict, file_props: dict, path: Optional[str]
) -> Any:
    """Start time of a channel."""
    if "wf_start_time" in props:
        start = _to_date(props["wf_start_time"])
        offset = float(props.get("wf_start_offset", 0.0))
        return Date.from_unix_ns(start.to_unix_ns() + round(offset * 1e9))
    value = _get_property(props, START_PROPERTIES)
    if value is None:
        value = _get_property(file_props, START_PROPERTIES)
    if value is None and path is not None:
        value = _file_name_time(path)
    return value


def _channel_spacing(props: dict) -> Optional[float]:
    """Channel spacing, corrected by the fibre length multiplier."""
    spacing = _get_property(props, SPACING_PROPERTIES)
    if spacing is None:
        return None
    multiplier = _get_property(
        props, ("Fibre Length Multiplier", "FibreLengthMultiplier")
    )
    return float(spacing) * (1.0 if multiplier is None
                             else float(multiplier))


def _get_property(props: dict, names: tuple) -> Any:
    """First property present among the names."""
    for name in names:
        if name in props:
            return props[name]
    return None


def _file_name_time(path: Any) -> Optional[str]:
    """Time in a file name ``..._UTC_YYYYMMDD_HHMMSS.fff...``."""
    match = re.search(
        r"_UTC_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}(\.\d+)?)",
        Path(path).name,
    )
    if match is None:
        return None
    y, mo, d, h, mi, s = match.groups()[:6]
    return f"{y}-{mo}-{d}T{h}:{mi}:{s}Z"


def _to_date(value: Any) -> Optional[Date]:
    """Date of a time value (Date, ISO string or None)."""
    if value is None or isinstance(value, Date):
        return value
    return Date(str(value))


def _float_or_none(value: Any) -> Optional[float]:
    """Float of a value, None if None."""
    return None if value is None else float(value)


def _channel_sort_key(name: str) -> tuple:
    """Numeric channel names in numeric order, then the others."""
    return (0, int(name), "") if name.isdigit() else (1, 0, name)


def _check_layout(path: Any, layout: dict, reference: dict) -> None:
    """Compatibility of a file with the first one of a DAS matrix."""
    if layout["names"] != reference["names"]:
        raise ValueError(f"the channels differ from the first file: {path}")
    if layout["dtype"] != reference["dtype"]:
        raise ValueError(f"the data type differs from the first file: "
                         f"{path}")
    for key, label in (("rate", "sampling rate"),
                       ("spacing", "channel spacing"),
                       ("gauge", "gauge length")):
        a, b = layout[key], reference[key]
        if (a is None) != (b is None) or (
            a is not None and not np.isclose(a, b, rtol=1e-9, atol=0)
        ):
            raise ValueError(f"the {label} differs from the first file: "
                             f"{path}")


def _resolve_paths(source: Any) -> list:
    """Sorted paths of a path, pattern or sequence of them."""
    sources = list(source) if isinstance(source, (list, tuple)) \
        else [source]
    if not sources:
        raise ValueError("no TDMS files given")
    paths = set()
    for item in sources:
        if not isinstance(item, (str, os.PathLike)):
            raise TypeError("TDMS sources must be paths or patterns")
        pattern = os.path.expanduser(os.fspath(item))
        if glob.has_magic(pattern):
            matches = [p for p in glob.glob(pattern) if os.path.isfile(p)]
            if not matches:
                raise FileNotFoundError(f"no TDMS files match {item}")
            paths.update(matches)
        elif not os.path.isfile(pattern):
            raise FileNotFoundError(f"TDMS file not found: {item}")
        else:
            paths.add(pattern)
    return sorted(paths)


def _append_record(stream_collection: Any, record: Any, index: dict) -> None:
    """
    Append a record using an index of the streams by identifier.

    ``StreamCollection.append`` looks streams up one by one, which is
    quadratic for the thousands of channels of DAS files.
    """
    sid = record.sid
    if sid in index:
        index[sid].append(record)
        return
    stream = stream_collection._new_stream(sid)
    stream.append(record)
    stream_collection.stream.append(stream)
    index[sid] = stream


def _warn(message: str) -> None:
    """Warning about the content of a file."""
    warnings.warn(message, TdmsWarning, stacklevel=4)
