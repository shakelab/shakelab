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
Unit tests for :mod:`shakelab.signals.circular`.

The tests feed the circular buffer with consecutive 5-second packets
cut from the IV.MODE record and check the time limit, the protection
of the caller's records, overlaps, selections and the refusal of
block-by-block processing.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, Stream, StreamCollection
from shakelab.signals.io import reader

with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    from shakelab.signals.circular import CircularStream

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def vertical():
    """Vertical record of IV.MODE (100 Hz, 180 s)."""
    return reader(str(DATA / "MODE.mseed"))["IV.MODE..HNZ"][0]


def _packets(record, count=6, length=5):
    """Consecutive packets of ``length`` seconds (500 samples each)."""
    return [
        record.extract(i * length, i * length + length - record.delta)
        for i in range(count)
    ]


# ----------------------------------------------------------------------------
# Buffer
# ----------------------------------------------------------------------------


def test_circular_buffer_keeps_latest_window(vertical):
    """Only the last max_duration seconds are retained."""
    stream = CircularStream(sid="IV.MODE..HNZ", max_duration=10.0)

    for packet in _packets(vertical):
        stream.append(packet)

    assert float(stream.endtime - stream.starttime) == pytest.approx(10.0)
    assert stream.endtime == Date("2012-05-20T02:03:52.99")
    assert stream.nsamp == 1001


def test_circular_buffer_does_not_modify_caller_records(vertical):
    """Trimming the buffer leaves the appended records unchanged."""
    packets = _packets(vertical)
    stream = CircularStream(max_duration=10.0)

    for packet in packets:
        stream.append(packet)

    assert packets[0].nsamp == 500
    assert packets[0].time == Date("2012-05-20T02:03:23")
    assert stream.sid == "IV.MODE..HNZ"


def test_circular_buffer_merged_segment_matches_original(vertical):
    """Merging the blocks restores the continuous signal."""
    stream = CircularStream(max_duration=10.0)

    for packet in _packets(vertical):
        stream.append(packet)

    segment = stream.extract(merge=True)

    assert isinstance(segment, Record)
    np.testing.assert_array_equal(
        segment.data,
        vertical.extract(segment.time, stream.endtime).data,
    )


def test_circular_buffer_overlaps_and_order(vertical):
    """Overlaps need enforce=True; older records are rejected."""
    stream = CircularStream(max_duration=60.0)

    for packet in _packets(vertical, count=3):
        stream.append(packet)

    with pytest.raises(ValueError):
        stream.append(vertical.extract(13, 16))

    stream.append(vertical.extract(13, 16), enforce=True)
    assert stream.endtime == Date("2012-05-20T02:03:39")
    assert stream[-2].endtime == Date("2012-05-20T02:03:35.99")

    with pytest.raises(ValueError):
        stream.append(vertical.extract(0, 1))


def test_circular_buffer_gaps_are_preserved(vertical):
    """A gap leaves separate blocks."""
    stream = CircularStream(max_duration=60.0)
    stream.append(vertical.extract(0, 4.99))
    stream.append(vertical.extract(10, 14.99))

    assert len(stream) == 2

    with pytest.raises(ValueError):
        stream.extract(merge=True)


# ----------------------------------------------------------------------------
# Selections and processing
# ----------------------------------------------------------------------------


def test_circular_selections_are_standard_streams(vertical):
    """Slices and extractions are standard, processable streams."""
    stream = CircularStream(max_duration=60.0)

    for packet in _packets(vertical, count=3):
        stream.append(packet)

    selection = stream[0:2]
    window = stream.extract(2, 8)

    assert type(selection) is Stream
    assert type(window) is Stream
    assert len(selection) == 2
    assert type(stream.copy()) is CircularStream
    assert stream.copy().max_duration == 60.0


@pytest.mark.parametrize(
    "method, args",
    [
        ("filter", (0.1, 1.0)),
        ("taper", ()),
        ("resample", (0.02,)),
        ("remove_mean", ()),
        ("integrate", ()),
        ("differentiate", ()),
        ("shift", (0.5,)),
    ],
)
def test_circular_processing_is_refused(vertical, method, args):
    """Block-by-block processing raises NotImplementedError."""
    stream = CircularStream(max_duration=60.0)
    stream.append(vertical.extract(0, 4.99))

    with pytest.raises(NotImplementedError):
        getattr(stream, method)(*args)


def test_circular_validation():
    """Durations, record types and reverse sorting are validated."""
    with pytest.raises(ValueError):
        CircularStream(max_duration=0.0)
    with pytest.raises(TypeError):
        CircularStream(max_duration=True)
    with pytest.raises(TypeError):
        CircularStream().append([1, 2])
    with pytest.raises(ValueError):
        CircularStream().sort(reverse=True)
    with pytest.raises(ValueError):
        CircularStream().append(Record(data=[1.0]))


def test_collection_with_max_duration_uses_circular_streams(vertical):
    """A collection with max_duration creates circular streams."""
    collection = StreamCollection(max_duration=10.0)

    for packet in _packets(vertical):
        collection.append(packet)

    stream = collection["IV.MODE..HNZ"]

    assert isinstance(stream, CircularStream)
    assert stream.max_duration == 10.0
    assert stream.nsamp == 1001

    with pytest.raises(NotImplementedError):
        collection.filter(0.1, 1.0)
