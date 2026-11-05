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
Tests of :mod:`shakelab.signals.libio.tdms`.

The files in ``data/tdms`` come from npTDMS (example files, and its
test scenarios with the expected data), were written with the npTDMS
writer, or are an excerpt of an iDAS acquisition (see the README of
the data folder).
"""

from __future__ import annotations

import io as _io
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.io import reader
from shakelab.signals.libio.tdms import (
    TdmsFile,
    TdmsWarning,
    tdms_to_matrix,
    tdmsread,
)

DATA = Path(__file__).parent / "data" / "tdms"
SCENARIOS = DATA / "scenarios"
GENERATED = DATA / "generated"
DAS = DATA / "das"
DAS_FILES = sorted(DAS.glob("*.tdms"))

# Scenarios whose last chunk has less data than declared.
PARTIAL = {
    "less_data_than_expected", "less_data_than_expected_interleaved",
    "incomplete_last_segment", "incomplete_last_row_of_interleaved_data",
    "single_segment_with_interleaved_data_of_different_width",
}


def _expected():
    """Expected data of the scenarios, by scenario."""
    expected = {}
    with np.load(SCENARIOS / "expected.npz") as npz:
        for key in npz.files:
            name, group, channel = key.split("/")
            expected.setdefault(name, {})[(group, channel)] = npz[key]
    return expected


EXPECTED = _expected()


# ----------------------------------------------------------------------------
# Binary format
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_scenarios_as_nptdms(name):
    """Every scenario of npTDMS gives its expected data."""
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore" if name in PARTIAL else "error", TdmsWarning
        )
        tdms = TdmsFile.read(SCENARIOS / f"{name}.tdms")

    for (group, channel), expected in EXPECTED[name].items():
        if expected.size == 0:
            continue
        data = tdms.groups[group].channels[channel].data
        assert data.dtype.kind == expected.dtype.kind
        if expected.dtype.kind in "fc":
            np.testing.assert_allclose(data, expected)
        else:
            np.testing.assert_array_equal(data, expected)


def test_partial_chunk_warns():
    """A last chunk with less data than declared gives a warning."""
    with pytest.warns(TdmsWarning, match="last chunk partial"):
        TdmsFile.read(SCENARIOS / "less_data_than_expected.tdms")


def test_big_endian_file():
    """Big-endian data and timestamps (example file of npTDMS)."""
    tdms = TdmsFile.read(DATA / "nptdms" / "big_endian.tdms")
    channel = tdms.groups["Measured Data"].channels["Amplitude sweep"]

    stamp = channel.properties["NI_ExpTimeStamp"].to_unix_ns()
    assert abs(stamp - Date("2018-11-13T23:04:49.403585Z").to_unix_ns()) \
        < 1000
    data = channel.data
    assert data.size == 3500 and data.dtype == np.float64


def test_example_files():
    """Example files of npTDMS: groups, channels and properties."""
    tdms = TdmsFile.read(DATA / "nptdms" / "Digital_Input.tdms")
    assert len(tdms.groups) == 3
    assert all(
        list(group.channels) == ["Dev1_port3_line7 - line 0"]
        for group in tdms.groups.values()
    )

    rec = tdmsread(DATA / "nptdms" / "raw_timestamps.tdms")[0][0]
    assert rec.nsamp == 128 and rec.delta == 0.001
    assert abs(rec.time.to_unix_ns()
               - Date("2024-01-24T01:48:43.068614Z").to_unix_ns()) < 1000


def test_daqmx_not_supported():
    """DAQmx raw data are reported as not supported."""
    with pytest.raises(NotImplementedError, match="DAQmx"):
        TdmsFile.read(DATA / "nptdms" / "raw1.tdms")


def test_text_and_numeric_channels():
    """Text channels do not shift the numeric ones."""
    group = TdmsFile.read(GENERATED / "mixed.tdms").groups["g"]

    np.testing.assert_array_equal(
        group.channels["a"].data, [0, 1, 2, 3, 4, 10, 11, 12, 13, 14]
    )
    np.testing.assert_array_equal(
        group.channels["b"].data, [1, 2, 3, 4, 5, 2, 4, 6, 8, 10]
    )
    assert list(group.channels["s"].data) == ["x", "yy", "zzz"] * 2


def test_truncated_file():
    """An interrupted file gives the complete data and a warning."""
    with pytest.warns(TdmsWarning):
        group = TdmsFile.read(GENERATED / "truncated.tdms").groups["g"]
    assert group.channels["a"].data.size == 5


def test_sources():
    """Paths, bytes and binary file objects are accepted."""
    path = GENERATED / "mixed.tdms"
    content = path.read_bytes()

    for source in (path, str(path), content, _io.BytesIO(content)):
        tdms = TdmsFile.read(source)
        assert tdms.groups["g"].channels["b"].data[-1] == 10
    assert TdmsFile.read(content).path is None
    with pytest.raises(TypeError, match="binary mode"):
        TdmsFile.read(_io.StringIO("x"))
    with pytest.raises(TypeError):
        TdmsFile.read(42)


@pytest.mark.parametrize(
    "content, error",
    [(b"", ValueError), (b"TDSx" + bytes(24), ValueError)],
)
def test_invalid_files(content, error, tmp_path):
    """Empty files and wrong segment tags are rejected."""
    path = tmp_path / "x.tdms"
    path.write_bytes(content)
    with pytest.raises(error):
        TdmsFile.read(path)


def test_info(capsys):
    """The summary lists groups and channels."""
    TdmsFile.read(GENERATED / "mixed.tdms").info()
    out = capsys.readouterr().out
    assert "Group 'g'" in out and "b: 10 values, int32" in out


# ----------------------------------------------------------------------------
# NI scales
# ----------------------------------------------------------------------------


def test_scales():
    """Linear, polynomial, chained and subtracted scales."""
    with pytest.warns(TdmsWarning, match="Thermocouple"):
        channels = TdmsFile.read(GENERATED / "scaled.tdms") \
            .groups["g"].channels
    raw = np.arange(6.0)

    np.testing.assert_allclose(channels["poly"].data,
                               1 + 2 * raw + 0.5 * raw**2)
    np.testing.assert_allclose(channels["chain"].data,
                               -1 + 10 * (2 * raw + 1))
    np.testing.assert_allclose(channels["subtract"].data,
                               (raw + 3) - 2 * raw)
    np.testing.assert_array_equal(channels["thermocouple"].data, raw)
    np.testing.assert_array_equal(channels["already"].data, raw)
    np.testing.assert_allclose(channels["group_scale"].data, 0.5 * raw)


def test_raw_values():
    """``scaled=False`` keeps the raw values."""
    channels = TdmsFile.read(GENERATED / "scaled.tdms", scaled=False) \
        .groups["g"].channels
    assert channels["poly"].data.dtype == np.int16
    np.testing.assert_array_equal(channels["chain"].data, np.arange(6))


# ----------------------------------------------------------------------------
# ShakeLab conversion
# ----------------------------------------------------------------------------


def test_waveform_channels():
    """LabVIEW waveform properties give rate, start time and units."""
    collection = tdmsread(GENERATED / "waveform.tdms", network="XX",
                          station="NI")
    acc = collection["XX.NI..acc"][0]
    vel = collection["XX.NI..vel"][0]

    assert acc.delta == 0.005 and vel.delta == 0.01
    assert acc.time == Date("2026-11-05T10:00:00.623456Z")
    assert vel.time == Date("2026-11-05T10:00:00.123456Z")
    assert acc.head.units == "g" and vel.head.units == "m/s"
    assert acc.head.extra["tdms"]["group"] == "Untitled"
    assert acc.head.extra["tdms"]["properties"]["wf_samples"] == 100


def test_arguments_take_precedence():
    """Rate and start time given as arguments, codes and type."""
    rec = tdmsread(
        GENERATED / "waveform.tdms", rate=50, start="2026-01-01T00:00:00",
        location="00", channel_prefix="C", dtype=np.float32,
    )["..00.Cvel"][0]

    assert rec.delta == 0.02
    assert rec.time == Date("2026-01-01T00:00:00Z")
    assert rec.data.dtype == np.float32


def test_unknown_rate_and_group():
    """A channel without rate and a missing group are errors."""
    with pytest.raises(ValueError, match="sampling rate of channel"):
        tdmsread(GENERATED / "mixed.tdms")
    with pytest.raises(KeyError, match="nothere"):
        tdmsread(GENERATED / "mixed.tdms", group="nothere")


def test_text_channels_skipped():
    """Only numeric channels become records."""
    collection = tdmsread(GENERATED / "mixed.tdms", rate=1)
    assert sorted(st.sid for st in collection) == ["...a", "...b"]


def test_das_files():
    """iDAS excerpt: properties, start time from GPS, two files joined."""
    collection = tdmsread(DAS_FILES[0], network="XX", station="DAS")
    rec = collection["XX.DAS..0"][0]

    assert len(collection) == 16
    assert rec.delta == 0.005
    assert rec.time == Date("2026-05-04T06:00:39.600619Z")
    assert rec.data.dtype == np.int16 and rec.nsamp == 200

    tdmsread(DAS_FILES[1], stream_collection=collection, network="XX",
             station="DAS")
    stream = collection["XX.DAS..15"]
    assert len(stream) == 1 and stream[0].nsamp == 400


def test_io_reader():
    """``io.reader`` reads TDMS files and patterns."""
    collection = reader(str(DAS / "*.tdms"))
    assert len(collection) == 16
    assert collection["...3"][0].nsamp == 400


# ----------------------------------------------------------------------------
# DAS matrix
# ----------------------------------------------------------------------------


def test_das_matrix():
    """Contiguous files joined in one matrix with the iDAS geometry."""
    das = tdms_to_matrix(str(DAS / "*.tdms"))

    assert das.data.shape == (16, 400)
    assert das.start_time == Date("2026-05-04T06:00:39.600619Z")
    assert das.sampling_rate == 200.0
    assert das.channel_spacing == pytest.approx(2.0 * 1.0209523838714072)
    assert das.gauge_length == 10.0
    assert das.time[-1] == pytest.approx(399 / 200.0)
    assert das.properties["source_files"] == [str(p) for p in DAS_FILES]

    single = tdms_to_matrix(DAS_FILES[1])
    np.testing.assert_array_equal(single.data, das.data[:, 200:])


def test_das_matrix_not_contiguous(tmp_path):
    """A second file that does not follow the first is rejected."""
    copy = tmp_path / "copy.tdms"
    copy.write_bytes(DAS_FILES[0].read_bytes())
    with pytest.raises(ValueError, match="not contiguous"):
        tdms_to_matrix([DAS_FILES[0], copy])


def test_das_matrix_missing_files(tmp_path):
    """Missing files and empty patterns are reported."""
    with pytest.raises(FileNotFoundError):
        tdms_to_matrix(tmp_path / "x.tdms")
    with pytest.raises(FileNotFoundError):
        tdms_to_matrix(str(tmp_path / "*.tdms"))
    with pytest.raises(ValueError):
        tdms_to_matrix([])
