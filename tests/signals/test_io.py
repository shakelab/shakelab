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
Unit tests for :mod:`shakelab.signals.io`.

They cover the reading of channel metadata (:func:`read_metadata`) and
their assignment while reading waveforms (``metadata`` in
:func:`reader` and :meth:`StreamCollection.read`). Real data come from
``tests/signals/data`` (see the README in that folder).
"""

from __future__ import annotations

import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.io import read_metadata, reader
from shakelab.signals.metadata import MetadataCollection, MetadataWarning

DATA = Path(__file__).parent / "data"
MSEED = str(DATA / "MODE.mseed")
SXML = str(DATA / "MODE.xml")
OTHER = str(DATA / "TT_T01.xml")


# ----------------------------------------------------------------------------
# read_metadata
# ----------------------------------------------------------------------------


def test_read_metadata_single_file():
    """One file gives the same content as the collection constructor."""
    meta = read_metadata(SXML)

    assert isinstance(meta, MetadataCollection)
    assert sorted(meta.sid) == sorted(MetadataCollection(SXML).sid)
    assert len(meta) == 3


def test_read_metadata_list():
    """The files of a list are joined in one collection."""
    meta = read_metadata([SXML, Path(OTHER)])

    assert len(meta) == 3 + 30
    assert len(meta.stations()) == 2


def test_read_metadata_pattern_and_directory(tmp_path):
    """Wildcard patterns and directories are expanded as in reader."""
    for path in (SXML, OTHER):
        shutil.copy(path, tmp_path)

    by_pattern = read_metadata(str(tmp_path / "*.xml"))
    by_directory = read_metadata(tmp_path)

    assert len(by_pattern) == len(by_directory) == 33


def test_read_metadata_repeated_file():
    """Epochs equal to those already read are ignored, silently."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        meta = read_metadata([SXML, SXML])

    assert len(meta) == 3


def test_read_metadata_ftype(tmp_path):
    """The file type can be given explicitly for unusual names."""
    path = tmp_path / "mode.station"
    shutil.copy(SXML, path)

    with pytest.raises(ValueError, match="Unknown file type"):
        read_metadata(path)

    assert len(read_metadata(path, ftype="xml")) == 3


def test_read_metadata_type_error():
    """Sources other than paths are rejected."""
    with pytest.raises(TypeError):
        read_metadata(42)

    with pytest.raises(TypeError):
        read_metadata(b"<FDSNStationXML/>")


# ----------------------------------------------------------------------------
# reader(..., metadata=...)
# ----------------------------------------------------------------------------


def test_reader_metadata_from_file():
    """Records read with a metadata file carry their channel epoch."""
    collection = reader(MSEED, metadata=SXML)
    meta = MetadataCollection(SXML)

    for stream in collection:
        record = stream[0]
        assert not record.head.meta.empty
        assert record.head.meta == meta.get(record.sid, record.time)

    assert collection["IV.MODE..HNE"][0].orientation.azimuth == 90.0


def test_reader_metadata_same_as_set_metadata():
    """Reading with metadata equals reading and then assigning them."""
    meta = MetadataCollection(SXML)

    direct = reader(MSEED, metadata=meta)
    later = reader(MSEED)
    later.set_metadata(meta)

    raw = direct["IV.MODE..HNZ"][0].copy()
    raw.deconvolve_response()
    expected = later["IV.MODE..HNZ"][0].copy()
    expected.deconvolve_response()

    np.testing.assert_allclose(raw.data, expected.data)


def test_reader_metadata_list_of_files():
    """A list of metadata files is read with read_metadata."""
    collection = reader(MSEED, metadata=[SXML, OTHER])

    assert all(not stream[0].head.meta.empty for stream in collection)


def test_reader_missing_channels_one_warning():
    """Channels without metadata are reported in a single warning."""
    with pytest.warns(MetadataWarning) as record:
        collection = reader(MSEED, metadata=OTHER)

    messages = [
        str(item.message) for item in record
        if issubclass(item.category, MetadataWarning)
    ]
    assert len(messages) == 1
    assert "3 channel(s)" in messages[0]
    for cha in ("HNE", "HNN", "HNZ"):
        assert f"IV.MODE..{cha}" in messages[0]

    assert all(stream[0].head.meta.empty for stream in collection)


def test_reader_existing_collection_untouched():
    """Records already in the collection keep their metadata."""
    collection = StreamCollection()
    collection.append(
        Record(data=[1.0, 2.0], delta=0.01, time="2026-01-01T00:00:00",
               sid="XX.STA..HHZ")
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", MetadataWarning)
        returned = reader(MSEED, stream_collection=collection,
                          metadata=SXML)

    assert returned is collection
    assert collection["XX.STA..HHZ"][0].head.meta.empty
    assert not collection["IV.MODE..HNZ"][0].head.meta.empty


def test_collection_read_metadata():
    """StreamCollection.read forwards the metadata to reader."""
    collection = StreamCollection()
    returned = collection.read(MSEED, metadata=SXML)

    assert returned is collection
    assert collection["IV.MODE..HNN"][0].orientation.azimuth == 0.0


def test_reader_metadata_checked_first(tmp_path):
    """Invalid metadata are reported before reading the data."""
    missing = str(tmp_path / "missing.mseed")

    with pytest.raises(TypeError, match="metadata must be"):
        reader(missing, metadata=42)

    with pytest.raises(FileNotFoundError):
        reader(MSEED, metadata=str(tmp_path / "missing.xml"))
