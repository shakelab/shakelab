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
Lightweight TDMS reader for ShakeLab.

This module implements a native reader for numeric TDMS files, with a focus
on DAS-oriented waveform data. It provides three levels of access:

- direct TDMS object loading;
- DAS matrix-oriented loading;
- conversion to ShakeLab StreamCollection objects.

Only standard numeric TDMS channels are currently supported. DAQmx raw data,
multi-dimensional arrays and full TDMS writing are intentionally not covered.
"""

from __future__ import annotations

import struct
from copy import copy
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection


TDS_TAG = b"TDSm"

TOC_METADATA = 1 << 1
TOC_NEW_OBJECT_LIST = 1 << 2
TOC_RAW_DATA = 1 << 3
TOC_INTERLEAVED = 1 << 5
TOC_BIG_ENDIAN = 1 << 6
TOC_DAQMX_RAW_DATA = 1 << 7

RAW_NO_DATA = 0xFFFFFFFF
RAW_MATCHES_PREVIOUS = 0x00000000
UINT64_MAX = 0xFFFFFFFFFFFFFFFF

NI_EPOCH = Date("1904-01-01T00:00:00Z")


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
    0x21: np.dtype("?"),
}


@dataclass
class _TdmsObject:
    """Internal TDMS object metadata valid for one segment."""

    path: str
    has_data: bool = False
    data_type: int | None = None
    number_values: int = 0
    data_size: int = 0
    properties: dict = field(default_factory=dict)


@dataclass
class TdmsChannel:
    """TDMS channel loaded in memory."""

    name: str
    path: str
    properties: dict
    data: np.ndarray


@dataclass
class TdmsGroup:
    """TDMS group containing channels and metadata."""

    name: str
    properties: dict = field(default_factory=dict)
    channels: dict[str, TdmsChannel] = field(default_factory=dict)


@dataclass
class TdmsFile:
    """TDMS file loaded in memory."""

    path: str
    properties: dict = field(default_factory=dict)
    groups: dict[str, TdmsGroup] = field(default_factory=dict)


@dataclass
class DasData:
    """DAS-oriented matrix representation."""

    data: np.ndarray
    time: np.ndarray | None
    distance: np.ndarray | None
    sampling_rate: float | None
    channel_spacing: float | None
    gauge_length: float | None
    start_time: Date | None
    properties: dict


@dataclass
class _Segment:
    """Internal TDMS segment descriptor."""

    position: int
    toc: int
    version: int
    data_position: int
    next_position: int
    endian: str
    ordered_objects: list[_TdmsObject]
    num_chunks: int = 0


# ============================================================================
# Public API
# ============================================================================

def tdms_read(path: str | Path) -> TdmsFile:
    """Read a TDMS file into memory.

    Parameters
    ----------
    path : str or pathlib.Path
        Input TDMS file.

    Returns
    -------
    TdmsFile
        TDMS file represented as file, group and channel objects.
    """
    path = Path(path)
    file_size = path.stat().st_size

    previous_objects = {}
    previous_segment_objects = None
    properties_by_path = {}
    data_by_path = {}

    with path.open("rb") as fobj:
        while True:
            position = fobj.tell()
            tag = fobj.read(4)

            if tag == b"":
                break

            if tag != TDS_TAG:
                raise ValueError(
                    f"Invalid TDMS segment tag at byte {position}: {tag!r}"
                )

            segment = _read_lead_in(fobj, position, file_size)

            if segment.toc & TOC_DAQMX_RAW_DATA:
                raise NotImplementedError("DAQmx TDMS data are unsupported")

            segment.ordered_objects = _read_segment_objects(
                fobj,
                segment,
                previous_objects,
                previous_segment_objects,
                properties_by_path,
            )

            _calculate_segment_chunks(segment)

            if segment.toc & TOC_RAW_DATA:
                _read_segment_data(fobj, segment, data_by_path)

            for obj in segment.ordered_objects:
                previous_objects[obj.path] = obj

            previous_segment_objects = segment.ordered_objects
            fobj.seek(segment.next_position)

    return _build_tdms_file(path, properties_by_path, data_by_path)


def tdms_das_read(path: str | Path,
                  group_name: str | None = None,
                  sampling_rate: float | None = None,
                  channel_spacing: float | None = None,
                  gauge_length: float | None = None,
                  start_time=None) -> DasData:
    """Read a TDMS file as a DAS matrix.

    Parameters
    ----------
    path : str or pathlib.Path
        Input TDMS file.
    group_name : str, optional
        TDMS group to read. If None, the most DAS-like group is selected.
    sampling_rate : float, optional
        Sampling rate in Hz. If None, metadata are inspected.
    channel_spacing : float, optional
        Channel spacing in metres. If None, metadata are inspected. For
        Silixa/iDAS files, the fibre-length multiplier is applied when found.
    gauge_length : float, optional
        DAS gauge length in metres. If None, metadata are inspected.
    start_time : str, Date or datetime, optional
        Absolute start time. If None, metadata and file name are inspected.

    Returns
    -------
    DasData
        DAS matrix and related axes and metadata.
    """
    tdms = tdms_read(path)
    group = _select_tdms_das_group(tdms, preferred=group_name)

    names = sorted(group.channels, key=_channel_sort_key)
    arrays = [group.channels[name].data for name in names]

    sizes = {arr.size for arr in arrays}
    if len(sizes) != 1:
        raise ValueError("DAS channels do not have the same length")

    data = np.vstack(arrays)
    properties = _merge_properties(tdms.properties, group.properties)

    if sampling_rate is None:
        sampling_rate = _guess_sampling_rate(properties)

    if channel_spacing is None:
        channel_spacing = _guess_channel_spacing(properties)

    if gauge_length is None:
        gauge_length = _guess_gauge_length(properties)

    if start_time is None:
        start_time = _guess_tdms_start_time(path, properties)

    start_time = _to_shakelab_date(start_time)

    time = None
    if sampling_rate:
        time = np.arange(data.shape[1], dtype=float)
        time /= float(sampling_rate)

    distance = None
    if channel_spacing:
        distance = np.arange(data.shape[0], dtype=float)
        distance *= float(channel_spacing)

    properties["tdms_group"] = group.name
    properties["channel_names"] = names

    return DasData(
        data=data,
        time=time,
        distance=distance,
        sampling_rate=sampling_rate,
        channel_spacing=channel_spacing,
        gauge_length=gauge_length,
        start_time=start_time,
        properties=properties,
    )


def tdms_stream_read(input_data_source,
                     stream_collection=None,
                     group_name=None,
                     sampling_rate=None,
                     start_time=None,
                     network="XX",
                     station="DAS",
                     location="",
                     channel_prefix="",
                     dtype=None,
                     fast_append=True) -> StreamCollection:
    """Read TDMS file(s) as a ShakeLab StreamCollection.

    Parameters
    ----------
    input_data_source : str, pathlib.Path, list or tuple
        TDMS file path, or list of TDMS file paths.
    stream_collection : StreamCollection, optional
        Existing StreamCollection where new records are appended.
        If None, a new StreamCollection is created.
    group_name : str, optional
        TDMS group to convert. If None, the most DAS-like group is selected.
    sampling_rate : float, optional
        Sampling rate in Hz. If None, metadata are inspected.
    start_time : str, Date or datetime, optional
        Record start time. If None, metadata and file name are inspected.
    network : str, optional
        Network code used to build the ShakeLab stream id.
    station : str, optional
        Station code used to build the ShakeLab stream id.
    location : str, optional
        Location code used to build the ShakeLab stream id.
    channel_prefix : str, optional
        Prefix added before the TDMS channel name.
    dtype : numpy dtype, optional
        Optional output data type. If None, the native dtype is preserved.
    fast_append : bool, optional
        If True, an optimized bulk-loading append strategy is used.

    Returns
    -------
    StreamCollection
        StreamCollection containing the converted TDMS channels.

    Notes
    -----
    This function is intended for compatibility with existing ShakeLab tools.
    For matrix-oriented DAS processing, use ``tdms_das_read`` instead.
    """
    if stream_collection is None:
        stream_collection = StreamCollection()

    if isinstance(input_data_source, (list, tuple)):
        for path in input_data_source:
            tdms_stream_read(
                path,
                stream_collection=stream_collection,
                group_name=group_name,
                sampling_rate=sampling_rate,
                start_time=start_time,
                network=network,
                station=station,
                location=location,
                channel_prefix=channel_prefix,
                dtype=dtype,
                fast_append=fast_append,
            )
        return stream_collection

    tdms = tdms_read(input_data_source)
    group = _select_tdms_das_group(tdms, preferred=group_name)
    properties = _merge_properties(tdms.properties, group.properties)

    if sampling_rate is None:
        sampling_rate = _guess_sampling_rate(properties)

    if sampling_rate is None:
        raise ValueError("sampling_rate is required or must exist in TDMS")

    delta = 1.0 / float(sampling_rate)

    if start_time is None:
        start_time = _guess_tdms_start_time(input_data_source, properties)

    start_time = _to_shakelab_date(start_time)
    names = sorted(group.channels, key=_channel_sort_key)

    stream_index = None
    if fast_append:
        stream_index = {
            stream.sid: stream
            for stream in stream_collection.stream
        }

    for name in names:
        channel = group.channels[name]

        record = Record()
        record.head.sid = _format_tdms_sid(
            network,
            station,
            location,
            channel_prefix + str(name),
        )
        record.head.delta = delta

        if start_time is not None:
            record.head.time = start_time

        if dtype is None:
            record.data = np.asarray(channel.data)
        else:
            record.data = np.asarray(channel.data, dtype=dtype)

        if fast_append:
            stream_index = _append_record_fast(
                stream_collection,
                record,
                index=stream_index,
            )
        else:
            stream_collection.append(record)

    return stream_collection


def tdms_info(tdms: TdmsFile) -> None:
    """Print a compact summary of a TDMS file."""
    print(f"TDMS file: {tdms.path}")
    print(f"File properties: {len(tdms.properties)}")
    tdms_list_groups(tdms)


def tdms_list_groups(tdms: TdmsFile) -> None:
    """Print TDMS groups and basic channel information."""
    for group_name, group in tdms.groups.items():
        channels = list(group.channels.values())
        sizes = [channel.data.size for channel in channels]
        dtypes = [str(channel.data.dtype) for channel in channels]

        print(f"\nGroup: {group_name}")
        print(f"  properties: {len(group.properties)}")
        print(f"  channels: {len(group.channels)}")

        if sizes:
            print(f"  sample sizes: {sorted(set(sizes))[:5]}")
            print(f"  dtypes: {sorted(set(dtypes))}")

        for index, channel in enumerate(channels):
            print(
                f"  {channel.name}: "
                f"shape={channel.data.shape}, "
                f"dtype={channel.data.dtype}"
            )
            if index >= 9:
                print("  ...")
                break


# ============================================================================
# TDMS binary parser
# ============================================================================

def _read_lead_in(fobj, position: int, file_size: int) -> _Segment:
    """Read a TDMS segment lead-in."""
    toc = _read_u32(fobj, "<")
    endian = ">" if toc & TOC_BIG_ENDIAN else "<"

    version = _read_u32(fobj, endian)
    next_offset = _read_u64(fobj, endian)
    raw_offset = _read_u64(fobj, endian)

    metadata_position = fobj.tell()
    data_position = metadata_position + raw_offset

    if next_offset == UINT64_MAX:
        next_position = file_size
    else:
        next_position = metadata_position + next_offset

    if next_position > file_size:
        next_position = file_size

    if data_position > next_position:
        raise EOFError("Incomplete TDMS segment metadata")

    return _Segment(
        position=position,
        toc=toc,
        version=version,
        data_position=data_position,
        next_position=next_position,
        endian=endian,
        ordered_objects=[],
    )


def _read_segment_objects(fobj,
                          segment: _Segment,
                          previous_objects: dict[str, _TdmsObject],
                          previous_segment_objects: list[_TdmsObject] | None,
                          properties_by_path: dict) -> list[_TdmsObject]:
    """Read or reuse the object list for a TDMS segment."""
    if not segment.toc & TOC_METADATA:
        if previous_segment_objects is None:
            raise ValueError("Segment has no metadata and no previous segment")
        return previous_segment_objects

    if segment.toc & TOC_NEW_OBJECT_LIST or previous_segment_objects is None:
        ordered_objects = []
        existing = {}
    else:
        ordered_objects = previous_segment_objects[:]
        existing = {
            obj.path: (index, obj)
            for index, obj in enumerate(ordered_objects)
        }

    n_objects = _read_u32(fobj, segment.endian)

    for _ in range(n_objects):
        path = _read_string(fobj, segment.endian)
        raw_header = _read_u32(fobj, segment.endian)

        index, existing_obj = existing.get(path, (None, None))

        if existing_obj is not None:
            obj = _update_existing_object(
                fobj,
                existing_obj,
                raw_header,
                segment.endian,
            )
            ordered_objects[index] = obj

        elif path in previous_objects:
            obj = _update_existing_object(
                fobj,
                previous_objects[path],
                raw_header,
                segment.endian,
            )
            ordered_objects.append(obj)

        else:
            obj = _TdmsObject(path=path)

            if raw_header == RAW_MATCHES_PREVIOUS:
                raise ValueError(
                    f"Object {path} reuses raw index but has no previous "
                    f"definition"
                )

            if raw_header != RAW_NO_DATA:
                _read_raw_index(fobj, obj, raw_header, segment.endian)

            ordered_objects.append(obj)

        props = _read_properties(fobj, segment.endian)

        if props:
            properties_by_path.setdefault(path, {}).update(props)
            obj.properties.update(props)

    return ordered_objects


def _update_existing_object(fobj,
                            obj: _TdmsObject,
                            raw_header: int,
                            endian: str) -> _TdmsObject:
    """Update an existing TDMS object definition."""
    new_obj = copy(obj)

    if raw_header == RAW_NO_DATA:
        new_obj.has_data = False
        return new_obj

    if raw_header == RAW_MATCHES_PREVIOUS:
        new_obj.has_data = True
        return new_obj

    _read_raw_index(fobj, new_obj, raw_header, endian)
    return new_obj


def _read_raw_index(fobj,
                    obj: _TdmsObject,
                    raw_header: int,
                    endian: str) -> None:
    """Read a standard TDMS raw data index."""
    if raw_header < 20:
        raise ValueError(f"Invalid raw data index length: {raw_header}")

    payload_size = raw_header - 4
    start = fobj.tell()

    data_type = _read_u32(fobj, endian)
    dimension = _read_u32(fobj, endian)
    number_values = _read_u64(fobj, endian)

    if dimension != 1:
        raise NotImplementedError("Only 1-D TDMS channel arrays supported")

    obj.has_data = True
    obj.data_type = data_type
    obj.number_values = int(number_values)

    if data_type == 0x20:
        obj.data_size = int(_read_u64(fobj, endian))
    else:
        dtype = _numpy_dtype(data_type, endian)
        obj.data_size = int(number_values) * dtype.itemsize

    consumed = fobj.tell() - start
    remaining = payload_size - consumed

    if remaining < 0:
        raise ValueError(
            f"Invalid raw index parsing: consumed {consumed} bytes, "
            f"but payload size is {payload_size}"
        )

    if remaining > 0:
        fobj.read(remaining)


def _calculate_segment_chunks(segment: _Segment) -> None:
    """Compute the number of raw-data chunks in a segment."""
    chunk_size = sum(
        obj.data_size for obj in segment.ordered_objects if obj.has_data
    )
    raw_size = segment.next_position - segment.data_position

    if chunk_size == 0:
        segment.num_chunks = 0
        return

    segment.num_chunks = raw_size // chunk_size

    if raw_size % chunk_size:
        segment.num_chunks += 1


def _read_segment_data(fobj,
                       segment: _Segment,
                       data_by_path: dict[str, list[np.ndarray]]) -> None:
    """Read numeric raw data from one TDMS segment."""
    data_objects = [
        obj for obj in segment.ordered_objects
        if obj.has_data and obj.data_type != 0x20
    ]

    if not data_objects or segment.num_chunks == 0:
        return

    fobj.seek(segment.data_position)

    if segment.toc & TOC_INTERLEAVED:
        _read_interleaved_data(fobj, segment, data_objects, data_by_path)
    else:
        _read_contiguous_data(fobj, segment, data_objects, data_by_path)


def _read_contiguous_data(
        fobj,
        segment: _Segment,
        data_objects: list[_TdmsObject],
        data_by_path: dict[str, list[np.ndarray]]) -> None:
    """Read contiguous TDMS raw data."""
    for _ in range(segment.num_chunks):
        for obj in data_objects:
            dtype = _numpy_dtype(obj.data_type, segment.endian)
            data = np.fromfile(fobj, dtype=dtype, count=obj.number_values)
            data_by_path.setdefault(obj.path, []).append(data)


def _read_interleaved_data(
        fobj,
        segment: _Segment,
        data_objects: list[_TdmsObject],
        data_by_path: dict[str, list[np.ndarray]]) -> None:
    """Read interleaved TDMS raw data."""
    lengths = {obj.number_values for obj in data_objects}

    if len(lengths) != 1:
        raise NotImplementedError(
            "Interleaved TDMS data with different channel lengths are "
            "not supported"
        )

    data_types = {obj.data_type for obj in data_objects}

    if len(data_types) == 1:
        _read_interleaved_homogeneous_data(
            fobj,
            segment,
            data_objects,
            data_by_path,
        )
        return

    _read_interleaved_heterogeneous_data(
        fobj,
        segment,
        data_objects,
        data_by_path,
    )


def _read_interleaved_homogeneous_data(
        fobj,
        segment: _Segment,
        data_objects: list[_TdmsObject],
        data_by_path: dict[str, list[np.ndarray]]) -> None:
    """Read interleaved data with homogeneous channel type."""
    n_channels = len(data_objects)
    n_samples = data_objects[0].number_values * segment.num_chunks
    dtype = _numpy_dtype(data_objects[0].data_type, segment.endian)

    raw = np.fromfile(
        fobj,
        dtype=dtype,
        count=n_samples * n_channels,
    )

    if raw.size != n_samples * n_channels:
        raise EOFError("Incomplete interleaved TDMS raw data")

    data = raw.reshape(n_samples, n_channels).T.copy()

    for index, obj in enumerate(data_objects):
        data_by_path.setdefault(obj.path, []).append(data[index])


def _read_interleaved_heterogeneous_data(
        fobj,
        segment: _Segment,
        data_objects: list[_TdmsObject],
        data_by_path: dict[str, list[np.ndarray]]) -> None:
    """Read interleaved data with heterogeneous channel types."""
    n_samples = data_objects[0].number_values * segment.num_chunks
    row_width = sum(
        _numpy_dtype(obj.data_type, segment.endian).itemsize
        for obj in data_objects
    )

    raw = np.fromfile(fobj, dtype=np.uint8, count=n_samples * row_width)

    if raw.size != n_samples * row_width:
        raise EOFError("Incomplete interleaved TDMS raw data")

    raw = raw.reshape((-1, row_width))
    col0 = 0

    for obj in data_objects:
        dtype = _numpy_dtype(obj.data_type, segment.endian)
        width = dtype.itemsize
        cols = raw[:, col0:col0 + width].ravel()
        values = cols.view(dtype).copy()
        data_by_path.setdefault(obj.path, []).append(values)
        col0 += width


def _read_properties(fobj, endian: str) -> dict:
    """Read TDMS object properties."""
    properties = {}
    n_properties = _read_u32(fobj, endian)

    for _ in range(n_properties):
        name = _read_string(fobj, endian)
        data_type = _read_u32(fobj, endian)
        properties[name] = _read_value(fobj, data_type, endian)

    return properties


def _read_value(fobj, data_type: int, endian: str):
    """Read one scalar TDMS property value."""
    if data_type == 0x01:
        return struct.unpack(endian + "b", fobj.read(1))[0]
    if data_type == 0x02:
        return struct.unpack(endian + "h", fobj.read(2))[0]
    if data_type == 0x03:
        return struct.unpack(endian + "i", fobj.read(4))[0]
    if data_type == 0x04:
        return struct.unpack(endian + "q", fobj.read(8))[0]
    if data_type == 0x05:
        return struct.unpack(endian + "B", fobj.read(1))[0]
    if data_type == 0x06:
        return struct.unpack(endian + "H", fobj.read(2))[0]
    if data_type == 0x07:
        return _read_u32(fobj, endian)
    if data_type == 0x08:
        return _read_u64(fobj, endian)
    if data_type == 0x09:
        return struct.unpack(endian + "f", fobj.read(4))[0]
    if data_type == 0x0A:
        return struct.unpack(endian + "d", fobj.read(8))[0]
    if data_type == 0x20:
        return _read_string(fobj, endian)
    if data_type == 0x21:
        return bool(struct.unpack(endian + "B", fobj.read(1))[0])
    if data_type == 0x44:
        fractions = _read_u64(fobj, endian)
        seconds = struct.unpack(endian + "q", fobj.read(8))[0]
        return NI_EPOCH + seconds + fractions / 2**64

    raise NotImplementedError(f"Unsupported TDMS property type: {data_type}")


def _read_string(fobj, endian: str) -> str:
    """Read a TDMS string."""
    nbytes = _read_u32(fobj, endian)
    data = fobj.read(nbytes)

    if len(data) != nbytes:
        raise EOFError("Unexpected end of file while reading TDMS string")

    return data.decode("utf-8").rstrip("\x00")


def _read_u32(fobj, endian: str) -> int:
    """Read an unsigned 32-bit integer."""
    data = fobj.read(4)

    if len(data) != 4:
        raise EOFError("Unexpected end of file while reading uint32")

    return struct.unpack(endian + "L", data)[0]


def _read_u64(fobj, endian: str) -> int:
    """Read an unsigned 64-bit integer."""
    data = fobj.read(8)

    if len(data) != 8:
        raise EOFError("Unexpected end of file while reading uint64")

    return struct.unpack(endian + "Q", data)[0]


def _numpy_dtype(data_type: int, endian: str) -> np.dtype:
    """Return NumPy dtype for a TDMS numeric type."""
    try:
        return TDMS_NUMPY_TYPES[data_type].newbyteorder(endian)
    except KeyError as exc:
        raise NotImplementedError(
            f"Unsupported TDMS numeric data type: {data_type}"
        ) from exc


# ============================================================================
# Object construction and selection helpers
# ============================================================================

def _build_tdms_file(path: Path,
                     properties_by_path: dict,
                     data_by_path: dict) -> TdmsFile:
    """Build the public TDMS object hierarchy."""
    tdms = TdmsFile(path=str(path))
    tdms.properties.update(properties_by_path.get("/", {}))

    for obj_path, arrays in data_by_path.items():
        group_name, channel_name = _parse_channel_path(obj_path)

        group = tdms.groups.setdefault(group_name, TdmsGroup(group_name))
        group_path = _make_group_path(group_name)

        group.properties.update(properties_by_path.get(group_path, {}))

        data = np.concatenate(arrays) if len(arrays) > 1 else arrays[0]

        group.channels[channel_name] = TdmsChannel(
            name=channel_name,
            path=obj_path,
            properties=properties_by_path.get(obj_path, {}),
            data=data,
        )

    return tdms


def _select_tdms_das_group(tdms: TdmsFile,
                           preferred: str | None = None) -> TdmsGroup:
    """Select the most likely DAS data group from a TDMS file."""
    if preferred is not None:
        try:
            return tdms.groups[preferred]
        except KeyError as exc:
            raise KeyError(f"TDMS group not found: {preferred}") from exc

    if len(tdms.groups) == 1:
        return next(iter(tdms.groups.values()))

    candidates = []

    for group_name, group in tdms.groups.items():
        score = _score_tdms_das_group(group_name, group)

        if score is not None:
            candidates.append((score, group_name, group))

    if not candidates:
        raise ValueError("No suitable DAS-like TDMS group found")

    candidates.sort(reverse=True, key=lambda item: item[0])
    return candidates[0][2]


def _score_tdms_das_group(group_name: str,
                          group: TdmsGroup) -> float | None:
    """Return a DAS-likeness score for a TDMS group."""
    channels = list(group.channels.values())

    if not channels:
        return None

    sizes = [channel.data.size for channel in channels]
    dtypes = [channel.data.dtype for channel in channels]

    same_size = len(set(sizes)) == 1
    numeric = all(np.issubdtype(dtype, np.number) for dtype in dtypes)

    if not same_size or not numeric:
        return None

    n_channels = len(channels)
    numeric_names = 0

    for channel in channels:
        try:
            int(channel.name)
            numeric_names += 1
        except ValueError:
            pass

    score = 0.0
    score += 10.0
    score += min(n_channels, 10000) / 100.0
    score += numeric_names / max(n_channels, 1)

    if group_name.lower() in ("measurement", "measurements", "data"):
        score += 5.0

    return score


def _merge_properties(*items: dict) -> dict:
    """Merge property dictionaries from left to right."""
    properties = {}

    for item in items:
        if item:
            properties.update(item)

    return properties


# ============================================================================
# TDMS paths, metadata and ShakeLab helpers
# ============================================================================

def _parse_channel_path(path: str) -> tuple[str, str]:
    """Parse /'group'/'channel' TDMS paths."""
    parts = _parse_path(path)

    if len(parts) != 2:
        raise ValueError(f"Object path is not a channel path: {path}")

    return parts[0], parts[1]


def _parse_path(path: str) -> list[str]:
    """Parse TDMS object path components."""
    if path == "/":
        return []

    parts = []
    index = 0

    while index < len(path):
        if path[index] != "/":
            raise ValueError(f"Invalid TDMS path: {path}")

        index += 1

        if index >= len(path) or path[index] != "'":
            raise ValueError(f"Invalid TDMS path: {path}")

        index += 1
        chars = []

        while index < len(path):
            char = path[index]

            if char == "'":
                if index + 1 < len(path) and path[index + 1] == "'":
                    chars.append("'")
                    index += 2
                    continue

                index += 1
                break

            chars.append(char)
            index += 1

        parts.append("".join(chars))

    return parts


def _make_group_path(group_name: str) -> str:
    """Return the TDMS path for a group name."""
    safe = group_name.replace("'", "''")
    return f"/'{safe}'"


def _channel_sort_key(name: str):
    """Return a robust sorting key for TDMS channel names."""
    try:
        return int(name)
    except ValueError:
        return name


def _guess_sampling_rate(properties: dict) -> float | None:
    """Try to infer sampling rate from common TDMS properties."""
    value = _get_property(
        properties,
        (
            "SamplingFrequency[Hz]",
            "sampling_rate",
            "Sampling Rate",
            "SampleRate",
            "sample_rate",
            "fs",
        ),
    )

    if value is not None:
        return float(value)

    if "wf_increment" in properties:
        return 1.0 / float(properties["wf_increment"])

    return None


def _guess_channel_spacing(properties: dict) -> float | None:
    """Try to infer corrected DAS channel spacing in metres."""
    spacing = _get_property(
        properties,
        (
            "SpatialResolution[m]",
            "SpatialResolution",
            "channel_spacing",
            "Channel Spacing",
        ),
    )

    if spacing is None:
        return None

    spacing = float(spacing)
    multiplier = _get_property(
        properties,
        (
            "Fibre Length Multiplier",
            "FibreLengthMultiplier",
        ),
    )

    if multiplier is not None:
        spacing *= float(multiplier)

    return spacing


def _guess_gauge_length(properties: dict) -> float | None:
    """Try to infer DAS gauge length in metres."""
    value = _get_property(
        properties,
        (
            "GaugeLength",
            "GaugeLength[m]",
            "gauge_length",
            "Gauge Length",
        ),
    )

    if value is None:
        return None

    return float(value)


def _get_property(properties: dict, names: tuple[str, ...]):
    """Return the first property matching one of the given names."""
    for name in names:
        if name in properties:
            return properties[name]

    return None


def _guess_tdms_start_time(path, properties=None):
    """Infer start time from TDMS properties or file name."""
    if properties is not None:
        for key in (
            "ISO8601 Timestamp",
            "GPSTimeStamp",
            "wf_start_time",
            "start_time",
            "StartTime",
        ):
            if key in properties:
                return properties[key]

    name = Path(path).stem

    if "_UTC_" not in name:
        return None

    try:
        text = name.split("_UTC_", 1)[1]
        date_text, time_text = text.split("_", 1)

        year = date_text[0:4]
        month = date_text[4:6]
        day = date_text[6:8]

        hour = time_text[0:2]
        minute = time_text[2:4]
        second = time_text[4:]

        return f"{year}-{month}-{day}T{hour}:{minute}:{second}"

    except Exception:
        return None


def _to_shakelab_date(value):
    """Convert common time representations to a ShakeLab Date."""
    if value is None:
        return None

    if isinstance(value, Date):
        return value

    return Date(str(value))


def _format_tdms_sid(network, station, location, channel):
    """Build a ShakeLab stream id from TDMS channel information."""
    return f"{network}.{station}.{location}.{channel}"


def _append_record_fast(stream_collection, record, index=None):
    """Append a Record to a StreamCollection using a SID index.

    Parameters
    ----------
    stream_collection : StreamCollection
        Target stream collection.
    record : Record
        Record to append.
    index : dict, optional
        Mapping between stream id and stream object.

    Returns
    -------
    dict
        Updated stream index.
    """
    if index is None:
        index = {
            stream.sid: stream
            for stream in stream_collection.stream
        }

    sid = record.sid

    if sid in index:
        index[sid].append(record)
        return index

    if stream_collection.max_duration:
        stream = stream_collection.stream_class(
            id=sid,
            max_duration=stream_collection.max_duration,
        )
    else:
        stream = stream_collection.stream_class(id=sid)

    stream.append(record)
    stream_collection.stream.append(stream)
    index[sid] = stream

    return index