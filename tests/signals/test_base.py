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
Unit tests for :mod:`shakelab.signals.base`.

The tests are organized by class (``Header``, ``Record``, ``Stream``,
``StreamCollection``). They fix the behavior of the public API used in
the ShakeLab demos and check the corrected defects. Real data come from
``tests/signals/data`` (station IV.MODE, Emilia earthquake of 20 May
2012; see the README in that folder).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from scipy.signal import windows

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.libutils.timeN import Date
from shakelab.signals.base import (
    GRID_TOLERANCE,
    MERGE_TOLERANCE,
    Header,
    Record,
    tukey_window,
)
from shakelab.signals.io import reader
from shakelab.signals.response import (
    ResponseCollection,
    StageGain,
    StagePoleZero,
)

DATA = Path(__file__).parent / "data"

START = "2026-01-01T00:00:00"


@pytest.fixture(scope="module")
def mode():
    """Streams of station IV.MODE (three components, 100 Hz)."""
    return reader(str(DATA / "MODE.mseed"))


@pytest.fixture(scope="module")
def mode_response():
    """Instrument responses of station IV.MODE."""
    collection = ResponseCollection()
    collection.read(str(DATA / "MODE.xml"))
    return collection


def _sine(nsamp=100, delta=0.01, frequency=1.0):
    """Record with a sine sampled at ``delta``, starting at START."""
    t = np.arange(nsamp) * delta
    return Record(
        data=np.sin(2.0 * np.pi * frequency * t),
        delta=delta,
        time=START,
    )


# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------


def test_header_initialization():
    """All constructor arguments are stored and converted."""
    header = Header(
        delta=0.01,
        time="2026-01-01T00:00:00",
        location=(13.0, 46.0),
        sid="IV.TEST..HHZ",
        eid="event_001",
        units="m/s",
    )

    assert header.delta == 0.01
    assert header.rate == 100.0
    assert header.sid == "IV.TEST..HHZ"
    assert header.eid == "event_001"
    assert header.units == "m/s"
    assert header.time == Date("2026-01-01T00:00:00")
    assert header.location == WgsPoint(longitude=13.0, latitude=46.0)
    assert header.response is None
    assert header.meta == {}


def test_header_defaults():
    """An empty header has no sampling, time or location."""
    header = Header()

    assert header.delta is None
    assert header.rate is None
    assert header.time is None
    assert header.location is None
    assert header.nsamp is None


def test_header_rate_and_delta_are_synchronized():
    """Setting the rate updates the interval and vice versa."""
    header = Header()

    header.rate = 50.0
    assert header.delta == 0.02

    header.delta = 0.005
    assert header.rate == 200.0

    header.delta = None
    assert header.rate is None


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
@pytest.mark.parametrize("name", ["delta", "rate"])
def test_header_invalid_sampling_values(name, value):
    """Sampling values must be finite and strictly positive."""
    with pytest.raises(ValueError):
        setattr(Header(), name, value)


@pytest.mark.parametrize("value", ["0.01", True, [0.01]])
def test_header_sampling_type_errors(value):
    """Sampling values must be real numbers, booleans excluded."""
    with pytest.raises(TypeError):
        Header(delta=value)


def test_header_accepts_numpy_sampling_values():
    """NumPy scalars are accepted and stored as float."""
    header = Header()
    header.rate = np.int64(100)

    assert header.delta == 0.01
    assert type(header.rate) is float


@pytest.mark.parametrize(
    "value",
    [
        "2012-05-20T02:03:23Z",
        Date("2012-05-20T02:03:23Z"),
        [2012, 5, 20, 2, 3, 23],
        Date("2012-05-20T02:03:23Z").to_seconds(),
    ],
)
def test_header_time_accepts_date_formats(value):
    """The start time accepts every format accepted by Date."""
    header = Header(time=value)

    assert header.time == Date("2012-05-20T02:03:23Z")


def test_header_time_is_copied():
    """Changing an assigned Date afterwards does not move the header."""
    date = Date("2012-05-20T02:03:23Z")
    header = Header(time=date)

    date.year = 2000

    assert header.time == Date("2012-05-20T02:03:23Z")
    assert header.time is not date


def test_header_invalid_time():
    """Unsupported or malformed times are rejected."""
    with pytest.raises(TypeError):
        Header(time=object())
    with pytest.raises(ValueError):
        Header(time="not a date")


@pytest.mark.parametrize(
    "value",
    [
        (13.0, 46.0, 250.0),
        [13.0, 46.0, 250.0],
        np.array([13.0, 46.0, 250.0]),
    ],
)
def test_header_location_sequences(value):
    """Locations are (longitude, latitude[, elevation]) sequences."""
    header = Header(location=value)

    assert header.location == WgsPoint(
        longitude=13.0, latitude=46.0, elevation=250.0
    )


def test_header_location_point_is_kept():
    """A WgsPoint is stored as is."""
    point = WgsPoint(longitude=13.0, latitude=46.0)

    assert Header(location=point).location is point


def test_header_invalid_location():
    """Wrong lengths, strings and invalid latitudes are rejected."""
    with pytest.raises(ValueError):
        Header(location=(13.0,))
    with pytest.raises(TypeError):
        Header(location="13.0, 46.0")
    with pytest.raises(TypeError):
        Header(location=13.0)
    with pytest.raises(ValueError):
        Header(location=(13.0, 95.0))


def test_header_nsamp_from_parent():
    """nsamp is read from the parent record."""
    record = Record(data=np.zeros(25), delta=0.01)

    assert record.head.nsamp == 25


def test_header_copy_is_independent():
    """A copy shares no mutable state and drops the parent."""
    parent = object()
    header = Header(
        delta=0.01,
        time="2026-01-01T00:00:00",
        sid="IV.TEST..HHZ",
        parent=parent,
    )
    header.meta["quality"] = "D"

    copied = header.copy()

    assert copied is not header
    assert copied._parent is None
    assert copied.delta == header.delta
    assert copied.sid == header.sid
    assert copied.time == header.time
    assert copied.time is not header.time
    assert copied.meta == header.meta
    assert copied.meta is not header.meta

    other = object()
    assert header.copy(parent=other)._parent is other


def test_header_text_representations(capsys):
    """repr is compact; str and info list all fields."""
    header = Header(delta=0.01, sid="IV.TEST..HHZ")

    assert repr(header).startswith("Header(sid='IV.TEST..HHZ'")
    assert "sampling rate: 100.0" in str(header)

    header.info()
    assert "source identifier: IV.TEST..HHZ" in capsys.readouterr().out


# ----------------------------------------------------------------------------
# Record: construction and properties
# ----------------------------------------------------------------------------


def test_record_construction():
    """Data become a copied one-dimensional array; metadata go to head."""
    values = [1, 2, 3]
    record = Record(
        data=values,
        delta=0.5,
        time=START,
        sid="IV.TEST..HHZ",
        eid="ev1",
        units="m/s",
    )
    values[0] = 99

    assert record.data.tolist() == [1, 2, 3]
    assert record.data.dtype.kind == "i"
    assert record.head.units == "m/s"
    assert record.head._parent is record
    assert (record.sid, record.eid) == ("IV.TEST..HHZ", "ev1")
    assert (record.delta, record.rate) == (0.5, 2.0)
    assert record.time == Date(START)


def test_record_scalar_and_invalid_data():
    """A scalar becomes one sample; multidimensional data are rejected."""
    assert Record(data=5.0).data.tolist() == [5.0]

    with pytest.raises(ValueError):
        Record(data=np.zeros((2, 2)))


def test_record_empty():
    """An empty record has no samples, peak or end time."""
    record = Record()

    assert len(record) == 0
    assert record.peak is None
    assert record.duration == 0.0
    assert record.endtime is None


def test_record_time_properties():
    """Duration spans first to last sample; endtime is the last sample."""
    record = Record(data=np.arange(11), delta=0.1, time=START)

    assert len(record) == record.nsamp == 11
    assert record.duration == pytest.approx(1.0)
    assert record.starttime == Date(START)
    assert record.endtime == Date("2026-01-01T00:00:01")
    assert record[2:4].tolist() == [2, 3]
    assert record.peak == 10


def test_record_copy_is_independent():
    """A copy shares neither data nor header."""
    record = Record(data=[1.0, 2.0], delta=1.0, time=START)
    copied = record.copy()

    copied.data[0] = 9.0
    copied.head.sid = "X"

    assert record.data[0] == 1.0
    assert record.sid is None
    assert copied.head._parent is copied


def test_record_text(capsys):
    """str and info describe start, duration, rate and size."""
    record = _sine()

    assert "duration 0.99s" in str(record)
    record.info()
    assert "100 samples" in capsys.readouterr().out


# ----------------------------------------------------------------------------
# Record: operators
# ----------------------------------------------------------------------------


def test_record_scalar_operators():
    """Scalar arithmetic returns new records and leaves operands intact."""
    record = Record(data=[1, 2, 3], delta=1.0)

    assert (record + 1).data.tolist() == [2, 3, 4]
    assert (1 + record).data.tolist() == [2, 3, 4]
    assert (record - 1).data.tolist() == [0, 1, 2]
    assert (2 - record).data.tolist() == [1, 0, -1]
    assert (record * 2).data.tolist() == [2, 4, 6]
    assert (2 * record).data.tolist() == [2, 4, 6]
    assert (record / 2).data.tolist() == [0.5, 1.0, 1.5]
    assert record.data.tolist() == [1, 2, 3]


def test_record_operator_errors():
    """Unsupported operands and division by zero raise."""
    record = Record(data=[1, 2, 3], delta=1.0)

    with pytest.raises(TypeError):
        record + "a"
    with pytest.raises(TypeError):
        record * [1, 2, 3]
    with pytest.raises(ZeroDivisionError):
        record / 0


def test_record_addition_of_records_concatenates():
    """Adding contiguous records concatenates; a gap is an error."""
    first = Record(data=[1, 2, 3], delta=1.0, time=START)
    second = Record(data=[4, 5], delta=1.0, time="2026-01-01T00:00:03")

    assert (first + second).data.tolist() == [1, 2, 3, 4, 5]
    assert first.data.tolist() == [1, 2, 3]

    late = Record(data=[4, 5], delta=1.0, time="2026-01-01T00:00:05")
    with pytest.raises(ValueError):
        first + late


# ----------------------------------------------------------------------------
# Record: instrument response
# ----------------------------------------------------------------------------


def test_record_response_correction_variants(mode, mode_response):
    """All ways of removing the response give the same record."""
    raw = mode["IV.MODE..HNZ"][0]

    expected = raw.copy()
    expected.deconvolve_response(mode_response)

    for corrected in (
        raw / mode_response,
        raw / mode_response["IV.MODE..HNZ"],
        (mode / mode_response)["IV.MODE..HNZ"][0],
    ):
        np.testing.assert_allclose(corrected.data, expected.data)


def test_record_response_correction_by_stages(mode, mode_response):
    """Stage-by-stage division matches the full response correction."""
    raw = mode["IV.MODE..HNZ"][0]
    gain_1 = StageGain({
        "description": None,
        "sensitivity": 1.0197,
        "frequency": 0.2,
        "stage_number": 1,
    })
    poles_zeros = StagePoleZero({
        "description": None,
        "input_units": "m/s**2",
        "output_units": "V",
        "pz_type": None,
        "normalization_factor": 246.0e11,
        "normalization_frequency": 0.2,
        "poles": [-981 - 1009j, -981 + 1009j, -3290 - 1263j, -3290 + 1263j],
        "zeros": [],
        "stage_number": 1,
    })
    gain_2 = StageGain({
        "description": None,
        "sensitivity": 314573,
        "frequency": 0,
        "stage_number": 2,
    })

    by_stages = raw / gain_1 / poles_zeros / gain_2

    assert by_stages.peak == pytest.approx(
        (raw / mode_response).peak, rel=1e-6
    )


def test_record_response_peak_matches_sensitivity(mode, mode_response):
    """
    The corrected PGA agrees with the total sensitivity.

    The accelerometer response is flat well above the signal band, so
    dividing by the sensitivity of the StationXML (320770 counts per
    m/s^2) must give nearly the same peak as the full correction.
    """
    raw = mode["IV.MODE..HNZ"][0]
    corrected = raw / mode_response

    assert corrected.peak == pytest.approx(raw.peak / 320770.0, rel=0.01)


# ----------------------------------------------------------------------------
# Record: time axis
# ----------------------------------------------------------------------------


def test_record_time_axis_references():
    """Relative, absolute and date-referenced axes in seconds."""
    record = Record(data=np.zeros(3), delta=0.5, time=START)

    np.testing.assert_allclose(record.time_axis(), [0.0, 0.5, 1.0])
    np.testing.assert_allclose(record.taxis, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(
        record.time_axis(shift=1.0), [1.0, 1.5, 2.0]
    )
    np.testing.assert_allclose(
        record.time_axis("absolute"),
        Date(START).to_seconds() + np.array([0.0, 0.5, 1.0]),
    )
    np.testing.assert_allclose(
        record.time_axis("2025-12-31T23:59:58"), [2.0, 2.5, 3.0]
    )
    np.testing.assert_allclose(
        record.time_axis(Date("2026-01-01T00:00:00.25")),
        [-0.25, 0.25, 0.75],
    )


def test_record_time_axis_errors():
    """Missing metadata and unknown references are rejected."""
    with pytest.raises(ValueError):
        Record(data=[1.0]).time_axis()
    with pytest.raises(ValueError):
        Record(data=[1.0], delta=1.0).time_axis("absolute")
    with pytest.raises(ValueError):
        Record(data=[1.0], delta=1.0, time=START).time_axis("bogus")


# ----------------------------------------------------------------------------
# Record: append
# ----------------------------------------------------------------------------


def test_merge_tolerance_constants():
    """Tolerances are defined once, as fractions of the interval."""
    assert MERGE_TOLERANCE == 0.5
    assert GRID_TOLERANCE == 1e-6


def test_record_append_contiguous_and_gaps():
    """Contiguous records are joined; gaps only with enforce=True."""
    record = Record(data=[1, 2, 3], delta=1.0, time=START)
    late = Record(data=[4, 5], delta=1.0, time="2026-01-01T00:00:05")

    assert record.append(late) is False
    assert record.data.tolist() == [1, 2, 3]

    assert record.append(late, enforce=True) is True
    assert record.data.tolist() == [1, 2, 3, 0, 0, 4, 5]


def test_record_append_overlap_with_enforce():
    """The appended record wins from its start time onward."""
    record = Record(data=[1, 2, 3], delta=1.0, time=START)
    overlap = Record(data=[9, 9], delta=1.0, time="2026-01-01T00:00:01")

    assert record.append(overlap, enforce=True) is True
    assert record.data.tolist() == [1, 9, 9]


def test_record_append_absorbs_time_rounding(mode):
    """
    Start times off by a fraction of a sample are still contiguous.

    MiniSEED 2 stores times with 0.1 ms resolution, so the start of a
    record can differ from the expected time by up to 50 us (0.005
    samples at 100 Hz).
    """
    raw = mode["IV.MODE..HNZ"][0]
    first = raw.extract(0.0, 10.0)
    second = raw.extract(10.01, 20.0)
    second.time = second.time + 4e-5

    assert first.append(second) is True
    np.testing.assert_array_equal(first.data, raw.extract(0.0, 20.0).data)


def test_record_append_incompatible():
    """Different identifiers or sampling intervals are rejected."""
    record = Record(data=[1], delta=1.0, time=START, sid="A")

    with pytest.raises(ValueError):
        record.append(
            Record(data=[2], delta=1.0, time="2026-01-01T00:00:01", sid="B")
        )
    with pytest.raises(ValueError):
        record.append(
            Record(data=[2], delta=0.5, time="2026-01-01T00:00:01")
        )
    with pytest.raises(TypeError):
        record.append([2])


def test_record_append_to_empty_record():
    """An empty record adopts the appended record."""
    record = Record()
    other = Record(data=[1, 2], delta=1.0, time=START, sid="A")

    assert record.append(other) is True
    assert record.sid == "A"
    assert record.data.tolist() == [1, 2]


# ----------------------------------------------------------------------------
# Record: cut and extract
# ----------------------------------------------------------------------------


def test_record_extract_on_grid(mode):
    """Bounds on samples are both included."""
    raw = mode["IV.MODE..HNZ"][0]
    window = raw.extract("2012-05-20T02:04:00", "2012-05-20T02:04:10")

    assert window.nsamp == 1001
    assert window.time == Date("2012-05-20T02:04:00")
    assert raw.nsamp == 18001


def test_record_extract_off_grid(mode):
    """Off-grid bounds keep the samples inside the window."""
    raw = mode["IV.MODE..HNZ"][0]
    window = raw.extract(
        "2012-05-20T02:04:00.003", "2012-05-20T02:04:10.005"
    )

    assert window.nsamp == 1000
    assert window.time == Date("2012-05-20T02:04:00.01")


def test_record_cut_relative_and_partial_windows():
    """Numbers are seconds from start; windows are clipped to data."""
    record = Record(data=np.arange(10), delta=1.0, time=START)

    assert record.extract(2, 4).data.tolist() == [2, 3, 4]
    assert record.extract(-5, 1).data.tolist() == [0, 1]
    assert record.extract(8.5).data.tolist() == [9]

    record.cut(2, 4)
    assert record.data.tolist() == [2, 3, 4]
    assert record.time == Date("2026-01-01T00:00:02")


def test_record_cut_errors():
    """Empty or reversed windows are rejected."""
    record = Record(data=np.arange(10), delta=1.0, time=START)

    with pytest.raises(ValueError):
        record.extract(2.2, 2.8)
    with pytest.raises(ValueError):
        record.extract(5, 3)
    with pytest.raises(ValueError):
        record.extract(20, 30)


# ----------------------------------------------------------------------------
# Record: processing
# ----------------------------------------------------------------------------


def test_record_filter():
    """A low-pass filter removes a high frequency and keeps a low one."""
    t = np.arange(2000) * 0.01
    low = np.sin(2.0 * np.pi * 0.5 * t)
    record = Record(
        data=low + np.sin(2.0 * np.pi * 20.0 * t), delta=0.01
    )

    record.filter(lowpass=5.0, order=np.int64(4))

    middle = slice(500, 1500)
    np.testing.assert_allclose(record.data[middle], low[middle], atol=0.01)


def test_record_filter_errors():
    """Invalid corners and parameters are rejected."""
    record = _sine()

    with pytest.raises(ValueError):
        record.filter()
    with pytest.raises(ValueError):
        record.filter(lowpass=50.0)
    with pytest.raises(ValueError):
        record.filter(10.0, 5.0)
    with pytest.raises(TypeError):
        record.filter(lowpass=5.0, order=2.0)


@pytest.mark.parametrize("new_delta", [0.02, 0.005])
def test_record_resample_fft(new_delta):
    """
    Fourier resampling of a periodic band-limited signal is exact.

    One sine period over the record is periodic and band-limited, so
    the resampled values must equal the sine at the new times.
    """
    record = _sine()
    record.resample(new_delta, method="fft")
    t = np.arange(record.nsamp) * new_delta

    assert record.delta == new_delta
    assert record.time == Date(START)
    np.testing.assert_allclose(
        record.data, np.sin(2.0 * np.pi * t), atol=1e-12
    )


def test_record_resample_interp_does_not_extrapolate():
    """The last resampled sample lies within the original record."""
    record = _sine()
    duration = record.duration
    record.resample(0.02, method="interp", antialias=False)
    t = np.arange(record.nsamp) * 0.02

    assert t[-1] <= duration
    np.testing.assert_allclose(
        record.data, np.sin(2.0 * np.pi * t), atol=1e-3
    )


def test_record_resample_errors():
    """Invalid intervals and non-fractional FFT ratios are rejected."""
    with pytest.raises(ValueError):
        _sine().resample(0.0)
    with pytest.raises(ValueError):
        _sine().resample(0.01, method="spline")
    with pytest.raises(ValueError):
        _sine().resample(0.01 * np.pi, method="fft")


def test_record_taper_and_window():
    """The taper zeroes the ends and keeps the center."""
    record = _sine()
    original = record.data.copy()
    record.taper(0.1)

    assert record.data[0] == pytest.approx(0.0)
    assert record.data[-1] == pytest.approx(0.0)
    assert record.data[50] == original[50]
    np.testing.assert_allclose(
        tukey_window(101, 0.3), windows.tukey(101, 0.3)
    )


def test_record_zero_padding():
    """Padding appends zeros by default doubling the length."""
    record = Record(data=[1, 2, 3], delta=0.5)
    record.zero_padding()
    assert record.data.tolist() == [1, 2, 3, 0, 0, 0]

    record = Record(data=[1, 2, 3], delta=0.5)
    record.zero_padding(1.0)
    assert record.data.tolist() == [1, 2, 3, 0, 0]

    with pytest.raises(ValueError):
        Record(data=[1, 2, 3], delta=0.5).zero_padding(0.3)


def test_record_shift():
    """An integer-sample shift moves samples within the window."""
    record = Record(data=np.r_[np.zeros(5), 1.0, np.zeros(10)], delta=1.0)
    record.shift(2.0)

    expected = np.r_[np.zeros(7), 1.0, np.zeros(8)]
    np.testing.assert_allclose(record.data, expected, atol=1e-12)
    assert record.data.dtype.kind == "f"


def test_record_shift_complex_data():
    """Complex data keep their imaginary part."""
    n = np.arange(64)
    record = Record(data=np.exp(2j * np.pi * n / 64), delta=1.0)
    record.shift(1.0, padding=False)

    np.testing.assert_allclose(
        record.data, np.exp(2j * np.pi * (n - 1) / 64), atol=1e-12
    )


def test_record_integrate_and_differentiate():
    """Fourier integration and differentiation of a periodic sine."""
    t = np.arange(100) * 0.01
    integral = _sine()
    integral.integrate()
    derivative = _sine()
    derivative.differentiate()

    np.testing.assert_allclose(
        integral.data, -np.cos(2.0 * np.pi * t) / (2.0 * np.pi), atol=1e-12
    )
    np.testing.assert_allclose(
        derivative.data, 2.0 * np.pi * np.cos(2.0 * np.pi * t), atol=1e-9
    )

    cumulative = Record(data=np.ones(5), delta=0.5)
    cumulative.integrate(method="cum")
    np.testing.assert_allclose(cumulative.data, [0.0, 0.5, 1.0, 1.5, 2.0])


def test_record_convolve_deconvolve_correlate():
    """Discrete convolution, deconvolution and correlation."""
    record = Record(data=[1.0, 2.0, 3.0], delta=1.0)
    record.convolve([1.0, 1.0])
    assert record.data.tolist() == [1.0, 3.0, 5.0, 3.0]

    remainder = record.deconvolve([1.0, 1.0])
    assert record.data.tolist() == [1.0, 2.0, 3.0]
    assert not np.any(remainder)

    record.correlate([0.0, 1.0, 0.0])
    assert record.data.tolist() == [0.0, 1.0, 2.0, 3.0, 0.0]


def test_record_spectrum_round_trip():
    """Record to Spectrum and back restores data and metadata."""
    record = _sine()
    restored = Record()
    restored.from_spectrum(record.to_spectrum())

    np.testing.assert_allclose(restored.data, record.data, atol=1e-12)
    assert restored.delta == record.delta
    assert restored.time == record.time
    assert restored.head._parent is restored
