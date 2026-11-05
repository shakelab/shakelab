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
Tests of :mod:`shakelab.signals.libio.dyna`.

The files in ``data/ESM/ASCII`` were downloaded from the ESM event-data
web service (ESM version 3.0, CC BY 4.0); they write
``DATABASE_VERSION`` and ``HEADER_FORMAT`` on one line.
"""

from __future__ import annotations

import io as _io
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.io import reader
from shakelab.signals.libio.dyna import (
    HEADER_FIELDS,
    Dyna,
    DynaWarning,
    dynaread,
    dynawrite,
)

ESM = Path(__file__).resolve().parents[1] / "data" / "ESM" / "ASCII"
ACC = ESM / "IV.MODE..HNZ.D.IT-2012-0008.ACC.MP.ASC"
VEL = ESM / "IV.MODE..HNZ.D.IT-2012-0008.VEL.MP.ASC"
DIS = ESM / "IV.MODE..HNZ.D.IT-2012-0008.DIS.MP.ASC"
RAW = ESM / "IV.MODE..HNZ.D.IT-2012-0008.ACC.CV.ASC"
SA = ESM / "IV.MODE..HNZ.D.IT-2012-0008.SA.MP.ASC"
SD = ESM / "IV.MODE..HNZ.D.IT-2012-0008.S.MP.ASC"


def _written(record, **options):
    """Text of a record written to DYNA."""
    buffer = _io.StringIO()
    dynawrite(buffer, record, **options)
    return buffer.getvalue()


def _header_lines(text):
    """Header lines of a DYNA text."""
    return [line for line in text.splitlines() if ":" in line]


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


def test_read_esm_file_keeps_all_samples():
    """
    All samples are read, starting from the first one.

    The ESM 3.0 header has 63 lines (two fields on one line); the old
    reader assumed 64 and lost the first sample.
    """
    dyna = Dyna.read(ACC)

    assert dyna.header["NDATA"] == 40000
    assert dyna.data.size == 40000
    assert dyna.data.dtype == np.float64
    assert list(dyna.data[:3]) == [0.0, 0.000006, 0.000013]
    assert dyna.header["LOW_CUT_FREQUENCY_HZ"] == 0.03
    assert dyna.header["HIGH_CUT_FREQUENCY_HZ"] == 40.0
    assert dyna.is_timeseries


def test_fields_on_one_line_are_separated():
    """``DATABASE_VERSION: HEADER_FORMAT: DYNA 1.2`` gives two fields."""
    dyna = Dyna.read(ACC)

    assert dyna.header["DATABASE_VERSION"] is None
    assert dyna.header["HEADER_FORMAT"] == "DYNA 1.2"
    assert len(dyna.header) == 64


def test_header_values_are_typed():
    """Empty and ``None`` values are None; numbers are converted."""
    h = Dyna.read(ACC).header

    assert h["EVENT_NAME"] is None
    assert h["MAGNITUDE_L"] is None
    assert h["MAGNITUDE_W"] == 6.1
    assert h["STATION_ELEVATION_M"] == 41.0
    assert h["FILTER_ORDER"] == 2
    assert h["PGA_CM/S^2"] == -32.733094
    assert h["TIME_PGA_S"] == 44.385
    assert h["EVENT_DATE_YYYYMMDD"] == "20120520"
    assert "PGV_CM/S" not in h


def test_dynaread_record():
    """The record carries identifier, start time, sampling and units."""
    record = dynaread(ACC)

    assert record.sid == "IV.MODE..HNZ"
    assert record.time == Date("2012-05-20T02:03:23Z")
    assert record.delta == 0.005
    assert record.nsamp == 40000
    assert record.head.units == "CM/S**2"
    assert record.head.eid == "IT-2012-0008"
    assert np.abs(record.data).max() == 32.733094


def test_dynaread_channel_metadata():
    """Station coordinates and sensor depth become channel metadata."""
    record = dynaread(ACC)

    loc = record.location
    assert (loc.longitude, loc.latitude, loc.elevation) == (
        10.949167, 44.629722, 41.0
    )
    assert record.head.meta.depth == 0.0
    assert record.head.meta.sid == "IV.MODE..HNZ"


def test_dynaread_extra_fields():
    """The other header fields are kept in ``head.extra["dyna"]``."""
    extra = dynaread(ACC).head.extra["dyna"]

    assert extra["MAGNITUDE_W"] == 6.1
    assert extra["VS30_M/S"] == 204.0
    assert extra["PROCESSING"] == "manual (Paolucci et al., 2011)"
    assert extra["DATA_LICENSE"].startswith("CC-BY4_0")
    for name in ("NETWORK", "UNITS", "NDATA", "PGA_CM/S^2",
                 "STATION_LATITUDE_DEGREE", "DATA_TYPE"):
        assert name not in extra


@pytest.mark.parametrize(
    "path, units, rate",
    [(VEL, "CM/S", 200.0), (DIS, "CM", 200.0), (RAW, "CM/S**2", 100.0)],
)
def test_dynaread_quantities(path, units, rate):
    """Velocity, displacement and unprocessed files are read."""
    record = dynaread(path)

    assert record.head.units == units
    assert record.rate == pytest.approx(rate)
    assert record.nsamp == Dyna.read(path).header["NDATA"]


def test_unprocessed_record_starts_earlier():
    """The converted (CV) record has its own start time and length."""
    record = dynaread(RAW)

    assert record.time == Date("2012-05-20T02:02:53Z")
    assert record.nsamp == 29801
    assert record.head.extra["dyna"]["PROCESSING"] == "converted to cm/s^2"


def test_milliseconds_of_start_time():
    """The milliseconds of the first sample are kept."""
    text = ACC.read_text().replace("020323.000", "020323.250", 1)

    assert dynaread(text).time == Date("2012-05-20T02:03:23.25Z")


def test_read_sources():
    """Paths, strings, content, bytes and file objects are accepted."""
    expected = dynaread(ACC).data
    text = ACC.read_text()

    for source in (str(ACC), text, text.encode(), _io.StringIO(text),
                   _io.BytesIO(text.encode())):
        assert np.array_equal(dynaread(source).data, expected)


def test_ndata_mismatch_warns():
    """A number of samples different from NDATA gives a warning."""
    text = ACC.read_text().replace("NDATA: 40000", "NDATA: 40001")

    with pytest.warns(DynaWarning, match="NDATA"):
        Dyna.read(text)


def test_invalid_numeric_value_warns():
    """A non-numeric value of a numeric field is ignored with a warning."""
    text = ACC.read_text().replace("VS30_M/S: 204", "VS30_M/S: abc")

    with pytest.warns(DynaWarning, match="VS30"):
        dyna = Dyna.read(text)
    assert dyna.header["VS30_M/S"] is None


def test_not_a_dyna_file():
    """Text without DYNA header is rejected."""
    with pytest.raises(ValueError, match="not a DYNA file"):
        Dyna.read("1.0\n2.0\n")


def test_invalid_data():
    """Non-numeric data are rejected."""
    text = ACC.read_text() + "abc\n"

    with pytest.raises(ValueError, match="invalid DYNA data"):
        Dyna.read(text)


def test_reader_selects_dyna_by_extension():
    """``io.reader`` reads ``.ASC`` files with the DYNA reader."""
    collection = reader(str(ACC))
    record = collection["IV.MODE..HNZ"][0]

    assert record.nsamp == 40000
    assert record.head.units == "CM/S**2"


# ----------------------------------------------------------------------------
# Spectra
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "path, data_type, units, first",
    [
        (SA, "ACCELERATION RESPONSE SPECTRUM", "cm/s^2", 32.740391),
        (SD, "DISPLACEMENT RESPONSE SPECTRUM", "cm", 0.000083),
    ],
)
def test_response_spectra(path, data_type, units, first):
    """Response spectra: 105 periods from 0.01 to 10 s."""
    dyna = Dyna.read(path)

    assert not dyna.is_timeseries
    assert dyna.header["DATA_TYPE"] == data_type
    assert dyna.header["UNITS"] == units
    assert dyna.abscissa.size == dyna.data.size == 105
    assert dyna.abscissa[0] == 0.01
    assert dyna.abscissa[-1] == 10.0
    assert dyna.data[0] == first
    assert dyna.start_ns is None


def test_dynaread_rejects_spectra():
    """``dynaread`` reads waveforms only."""
    with pytest.raises(ValueError, match="RESPONSE SPECTRUM"):
        dynaread(SA)


def test_spectrum_text_round_trip():
    """A spectrum written with the ESM format gives the same file."""
    original = SA.read_text()
    text = Dyna.read(SA).to_text("%12.6f")

    assert text.splitlines()[64:] == original.splitlines()[63:]
    assert len(_header_lines(text)) == 64
    again = Dyna.read(text)
    assert again.header == Dyna.read(SA).header


def test_fourier_spectrum_two_columns():
    """A Fourier amplitude spectrum: frequencies and amplitudes."""
    header = Dyna.read(ACC).header
    header.update(
        DATA_TYPE="FOURIER AMPLITUDE SPECTRUM", NDATA=3,
        DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS=None,
    )
    text = Dyna(header, [0.0, 1.5, 2.0], [0.0, 0.5, 1.0]).to_text()

    dyna = Dyna.read(text)
    assert list(dyna.abscissa) == [0.0, 0.5, 1.0]
    assert list(dyna.data) == [0.0, 1.5, 2.0]


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


def test_header_text_matches_esm():
    """Rewritten header: ESM fields and order, same values."""
    original = _header_lines(ACC.read_text())
    text = Dyna.read(ACC).to_text("%f")
    written = _header_lines(text)

    names = [line.split(":")[0] for line in original]
    k = names.index("DATABASE_VERSION")
    names.insert(k + 1, "HEADER_FORMAT")
    assert len(written) == 64
    assert [line.split(":")[0] for line in written] == names
    assert Dyna.read(text).header == Dyna.read(ACC).header
    assert text.splitlines()[64:] == ACC.read_text().splitlines()[63:]


def test_round_trip_record():
    """Read, write and read again: same data and metadata."""
    record = dynaread(ACC)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        again = dynaread(_written(record))

    assert np.array_equal(again.data, record.data)
    assert again.sid == record.sid
    assert again.time == record.time
    assert again.delta == record.delta
    assert again.head.units == record.head.units
    assert again.head.eid == record.head.eid
    assert again.location == record.location
    assert again.head.meta.depth == record.head.meta.depth
    assert again.head.extra["dyna"] == record.head.extra["dyna"]


@pytest.mark.parametrize("path", [ACC, VEL, DIS, RAW])
def test_round_trip_is_lossless(path):
    """The default format keeps all the digits of the ESM files."""
    record = dynaread(path)

    assert np.array_equal(dynaread(_written(record)).data, record.data)


def test_computed_fields():
    """Number of samples, duration and peak come from the data."""
    record = dynaread(VEL)
    record.data = record.data * 2.0
    h = Dyna.read(_written(record)).header

    k = int(np.argmax(np.abs(record.data)))
    assert h["NDATA"] == 40000
    assert h["DURATION_S"] == 200.0
    assert h["PGV_CM/S"] == pytest.approx(record.data[k])
    assert h["TIME_PGV_S"] == pytest.approx(k * 0.005)
    assert h["UNITS"] == "cm/s"
    assert h["DATA_TYPE"] == "VELOCITY"
    assert "PGA_CM/S^2" not in h


def test_start_time_written_to_the_millisecond():
    """The start time is written with its milliseconds."""
    record = dynaread(ACC)
    record.head.time = Date("2012-05-20T02:03:23.2504Z")
    text = _written(record)

    assert "DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS: 20120520_020323.250" \
        in text
    assert "DATE_TIME_FIRST_SAMPLE_PRECISION: milliseconds" in text
    assert dynaread(text).time == Date("2012-05-20T02:03:23.25Z")


def test_si_units_are_converted():
    """Records in m/s**2 are written in cm/s**2."""
    record = dynaread(ACC)
    si = record.copy()
    si.data = record.data / 100.0
    si.head.units = "M/S**2"

    again = dynaread(_written(si))
    assert again.head.units == "CM/S**2"
    assert np.allclose(again.data, record.data, rtol=0, atol=1e-9)


@pytest.mark.parametrize("units", [None, "COUNTS", "M/S**3"])
def test_unsupported_units(units):
    """Records without units or in other units are rejected."""
    record = dynaread(ACC)
    record.head.units = units

    with pytest.raises(ValueError, match="not supported"):
        _written(record)


def test_write_new_record():
    """A record built in ShakeLab is written with its metadata."""
    record = Record(
        data=np.array([0.0, 1.0, -3.0, 2.0]), delta=0.01,
        time="2020-01-01T00:00:00.5Z",
    )
    record.sid = "XX.STA.00.HNZ"
    record.location = (12.0, 45.0, 100.0)
    record.head.units = "M/S"
    text = _written(record)

    assert len(_header_lines(text)) == 64
    h = Dyna.read(text).header
    assert h["LOCATION"] == "00"
    assert h["STATION_LATITUDE_DEGREE"] == 45.0
    assert h["UNITS"] == "cm/s"
    assert h["PGV_CM/S"] == -300.0
    assert h["TIME_PGV_S"] == 0.02
    assert h["HEADER_FORMAT"] == "DYNA 1.2"
    again = dynaread(text)
    assert again.sid == "XX.STA.00.HNZ"
    assert list(again.data) == [0.0, 100.0, -300.0, 200.0]


def test_write_file(tmp_path):
    """Files are not overwritten unless requested."""
    record = dynaread(ACC)
    path = tmp_path / "out.ASC"
    dynawrite(path, record)

    with pytest.raises(FileExistsError):
        dynawrite(path, record)
    dynawrite(path, record, owrite=True)
    assert dynaread(path).nsamp == 40000


def test_header_must_be_a_dict():
    """The header of a Dyna is a dictionary."""
    with pytest.raises(TypeError):
        Dyna(header=[("NDATA", 1)])


def test_header_fields():
    """All known fields, with the three peak pairs."""
    assert len(HEADER_FIELDS) == 68
    assert HEADER_FIELDS[0] == "EVENT_NAME"
    assert HEADER_FIELDS[-1] == "USER5"
