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
Tests of :mod:`shakelab.signals.libio.sac`.

The files in ``data/sac/obspy`` come from the test data of ObsPy
(``seism.sac`` is the example distributed with SAC); the expected
values follow the SAC file format description.
"""

from __future__ import annotations

import io as _io
import struct
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.libio.sac import (
    FOOTER_FIELDS,
    HEADER_SIZE,
    Sac,
    sacread,
    sacwrite,
)
from shakelab.signals.metadata import (
    ChannelEpoch,
    MetadataCollection,
    MetadataWarning,
)

DATA = Path(__file__).parent / "data" / "sac" / "obspy"


def _round_trip(record, byte_order="be"):
    """Write and read a record."""
    buffer = _io.BytesIO()
    sacwrite(buffer, record, byte_order=byte_order)
    return sacread(buffer.getvalue())


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


def test_reference_time_plus_begin():
    """The first sample is at the reference time plus B."""
    record = sacread(DATA / "seism.sac")

    # Reference 1981-088 10:38:14.000, B = 9.46 (single precision).
    assert record.time - Date("1981-03-29T10:38:14") == pytest.approx(
        9.46, abs=1e-6
    )
    assert record.nsamp == 1000
    assert record.delta == pytest.approx(0.01)


def test_milliseconds():
    """NZMSEC is in milliseconds."""
    record = sacread(DATA / "LMOW.BHE.SAC")
    assert record.time == Date("2001-04-10T00:23:00.465")


def test_byte_orders():
    """Big- and little-endian files give the same record."""
    little = sacread(DATA / "test.sac")
    big = sacread(DATA / "test.sac.swap")
    orders = [
        Sac.from_bytes((DATA / name).read_bytes()).byte_order
        for name in ("test.sac", "test.sac.swap")
    ]

    assert orders == ["le", "be"]
    assert big.time == little.time
    np.testing.assert_array_equal(big.data, little.data)


@pytest.mark.parametrize(
    "name, sid",
    [("seism.sac", ".CDV..Q"), ("LMOW.BHE.SAC", ".LMOW..BHE"),
     ("null_terminated.sac", "GD.PIN1..LYE"),
     ("dis.G.SCZ.__.BHE_short", "G.SCZ..BHE")],
)
def test_source_identifier(name, sid):
    """Undefined codes are empty; null-terminated text is cut."""
    assert sacread(DATA / name).sid == sid


def test_text_written_as_integers():
    """Text fields holding the undefined integer are undefined."""
    sac = Sac.from_bytes((DATA / "test_encode.sac").read_bytes())
    assert sac.head["KSTNM"] is None and sac.head["KEVNM"] is None


def test_channel_metadata():
    """Coordinates, burial depth and orientation go to head.meta."""
    record = sacread(DATA / "dis.G.SCZ.__.BHE_short")
    meta = record.head.meta

    assert meta.location.longitude == pytest.approx(-121.403, abs=1e-5)
    assert meta.location.latitude == pytest.approx(36.598, abs=1e-5)
    assert meta.location.elevation == pytest.approx(261.0)
    assert meta.depth == 0.0
    # CMPAZ 90, CMPINC 90 (horizontal east): SEED dip 0.
    assert meta.orientation.azimuth == pytest.approx(90.0)
    assert meta.orientation.dip == pytest.approx(0.0)


def test_vertical_orientation():
    """CMPINC 0 (up) is SEED dip -90."""
    assert sacread(DATA / "seism.sac").orientation.dip == -90.0


def test_extra_header_fields():
    """Event, picks and other fields are kept; markers are dates."""
    extra = sacread(DATA / "seism.sac").head.extra["sac"]

    assert extra["O"] == Date("1981-03-29T10:38:14")
    assert abs(extra["A"] - Date("1981-03-29T10:38:24.47")) < 1e-5
    assert extra["EVLA"] == pytest.approx(48.0, abs=1e-4)
    assert extra["IDEP"] == 50
    assert "B" not in extra and "STLA" not in extra


@pytest.mark.parametrize("name", ["testxy.sac", "seism.sac"])
def test_not_a_time_series(name, tmp_path):
    """Alphanumeric files and non time series give clear errors."""
    if name == "testxy.sac":
        with pytest.raises(ValueError, match="not a binary SAC file"):
            sacread(DATA / name)
        return

    sac = Sac.from_bytes((DATA / name).read_bytes())
    sac.head["IFTYPE"] = 4  # x-y data
    sac.data2 = np.arange(sac.data.size, dtype=np.float32)
    content = sac.to_bytes()

    raw = Sac.from_bytes(content)
    np.testing.assert_array_equal(raw.data2, sac.data2)
    with pytest.raises(ValueError, match="not an evenly sampled"):
        sacread(content)


def test_header_version_7():
    """The double-precision footer takes precedence."""
    sac = Sac.from_bytes((DATA / "seism.sac").read_bytes())
    content = bytearray(sac.to_bytes())
    content[4 * 76:4 * 77] = struct.pack(">i", 7)
    footer = [sac.head[name] if sac.head[name] is not None else -12345.0
              for name in FOOTER_FIELDS]
    footer[FOOTER_FIELDS.index("B")] = 9.4600001234
    content += struct.pack(f">{len(footer)}d", *footer)

    record = sacread(bytes(content))

    assert record.time - Date("1981-03-29T10:38:14") == pytest.approx(
        9.4600001234, abs=1e-9
    )


def test_truncated_file():
    """A file shorter than declared is rejected."""
    content = (DATA / "seism.sac").read_bytes()
    with pytest.raises(ValueError, match="shorter than declared"):
        sacread(content[:HEADER_SIZE + 100])


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("byte_order", ["be", "le"])
def test_round_trip_microseconds(byte_order):
    """Start times keep the microseconds (reference + B)."""
    record = Record(data=np.arange(100, dtype=np.float32), delta=0.005,
                    time="2026-03-01T10:20:30.123456", sid="XX.STA.00.HHZ")

    restored = _round_trip(record, byte_order)

    assert restored.time == record.time
    assert restored.sid == "XX.STA.00.HHZ"
    np.testing.assert_array_equal(restored.data, record.data)


def test_round_trip_keeps_header():
    """Picks, event and metadata survive a read-write cycle."""
    record = sacread(DATA / "seism.sac")
    restored = _round_trip(record)
    before = record.head.extra["sac"]
    after = restored.head.extra["sac"]

    assert restored.time == record.time
    assert restored.head.meta == record.head.meta
    assert before.keys() == after.keys()
    for name, value in before.items():
        if isinstance(value, Date):
            assert abs(after[name] - value) < 2e-6
        elif isinstance(value, float):
            assert after[name] == pytest.approx(value, rel=1e-6)
        else:
            assert after[name] == value


def test_metadata_written():
    """The channel metadata of a record are written to the header."""
    record = Record(data=np.zeros(10), delta=0.01,
                    time="2026-01-01T00:00:00", sid="XX.STA..HHN")
    record.head.meta = ChannelEpoch(
        "XX.STA..HHN", location=(13.25, 46.05, 120.0), depth=0.015,
        orientation=(10.0, 0.0),
    )
    sac = Sac.from_bytes(_io_bytes(record))

    assert sac.head["STLO"] == pytest.approx(13.25)
    assert sac.head["STLA"] == pytest.approx(46.05)
    assert sac.head["STEL"] == pytest.approx(120.0)
    assert sac.head["STDP"] == pytest.approx(15.0)
    assert sac.head["CMPAZ"] == pytest.approx(10.0)
    assert sac.head["CMPINC"] == pytest.approx(90.0)
    assert sac.head["IZTYPE"] == 9


def _io_bytes(record):
    """SAC content of a record."""
    buffer = _io.BytesIO()
    sacwrite(buffer, record)
    return buffer.getvalue()


def test_no_identifier_or_time():
    """Records without identifier or time can be written."""
    restored = _round_trip(Record(data=[1.0, 2.0], delta=1.0))
    assert restored.sid is None and restored.time is None


def test_existing_file(tmp_path):
    """Files are not overwritten unless requested."""
    path = tmp_path / "x.sac"
    record = Record(data=[1.0], delta=1.0, time="2026-01-01T00:00:00")
    sacwrite(path, record)

    with pytest.raises(FileExistsError):
        sacwrite(path, record)
    sacwrite(path, record, owrite=True)


def test_sac_read_sources():
    """``Sac.read`` accepts paths, bytes and binary file objects."""
    path = DATA / "seism.sac"
    content = path.read_bytes()
    expected = Sac.from_bytes(content)

    for source in (path, str(path), content, _io.BytesIO(content)):
        sac = Sac.read(source)
        assert sac.head == expected.head
        np.testing.assert_array_equal(sac.data, expected.data)
    with pytest.raises(TypeError, match="binary mode"):
        Sac.read(_io.StringIO("text"))


@pytest.mark.parametrize("byte_order", ["be", "le"])
def test_sac_write_keeps_non_time_series(byte_order, tmp_path):
    """``Sac.write`` writes what ``sacwrite`` cannot: x-y data."""
    sac = Sac.read(DATA / "seism.sac")
    sac.head["IFTYPE"] = 4  # x-y data
    sac.data2 = np.arange(sac.data.size, dtype=np.float32)
    path = tmp_path / "xy.sac"
    sac.write(path, byte_order)

    again = Sac.read(path)
    assert again.byte_order == byte_order
    assert again.head["IFTYPE"] == 4
    np.testing.assert_array_equal(again.data2, sac.data2)
    with pytest.raises(FileExistsError):
        sac.write(path)
    sac.write(path, owrite=True)
    with pytest.raises(TypeError, match="binary file object"):
        sac.write(42)


# ----------------------------------------------------------------------------
# Metadata from files and from collections
# ----------------------------------------------------------------------------


def test_consistent_metadata_without_warning():
    """Partial metadata from SAC agree with a complete epoch."""
    record = sacread(DATA / "null_terminated.sac")
    meta = record.head.meta
    complete = ChannelEpoch(
        "GD.PIN1..LYE", location=meta.location, depth=meta.depth,
        orientation=meta.orientation, sample_rate=1.0 / record.delta,
    )
    collection = MetadataCollection()
    collection.append(complete)

    with warnings.catch_warnings():
        warnings.simplefilter("error", MetadataWarning)
        record.set_metadata(collection)

    assert record.head.meta.sample_rate == complete.sample_rate


def test_conflicting_metadata_warn():
    """Different coordinates are reported with the attribute name."""
    record = sacread(DATA / "null_terminated.sac")
    other = ChannelEpoch("GD.PIN1..LYE", location=(0.0, 0.0))

    with pytest.warns(MetadataWarning, match="different location"):
        record.set_metadata(other)


def test_reader_dispatch():
    """io.reader recognizes SAC files by extension (any case)."""
    from shakelab.signals.io import reader

    collection = reader(str(DATA / "LMOW.BHE.SAC"))
    assert collection.sid == [".LMOW..BHE"]
