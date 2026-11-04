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
    SNAP_TOLERANCE,
    MERGE_TOLERANCE,
    Header,
    Record,
    Stream,
    StreamCollection,
    tukey_window,
)
from shakelab.signals.io import reader
from shakelab.signals.metadata import (
    ChannelEpoch,
    MetadataCollection,
    MetadataWarning,
)
from shakelab.signals.response import (
    StageGain,
    StagePoleZero,
    StageSet,
)

DATA = Path(__file__).parent / "data"

START = "2026-01-01T00:00:00"


@pytest.fixture(scope="module")
def mode():
    """Streams of station IV.MODE (three components, 100 Hz)."""
    return reader(str(DATA / "MODE.mseed"))


@pytest.fixture(scope="module")
def mode_response():
    """Channel metadata (with responses) of station IV.MODE."""
    collection = MetadataCollection()
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
    assert isinstance(header.meta, ChannelEpoch)
    assert header.meta.sid == "IV.TEST..HHZ"
    assert header.meta.empty
    assert header.extra == {}


def test_header_defaults():
    """An empty header has no sampling, time or identifier."""
    header = Header()

    assert header.delta is None
    assert header.rate is None
    assert header.time is None
    assert header.sid is None
    assert header.meta.empty
    assert header.nsamp is None


def test_header_meta_and_sid():
    """The identifier lives in the channel epoch, stored by reference."""
    epoch = ChannelEpoch("IV.TEST..HHZ", location=(13.0, 46.0))
    header = Header(meta=epoch)

    assert header.meta is epoch
    assert header.sid == "IV.TEST..HHZ"

    header.sid = "IV.TEST..HHN"
    assert header.meta is not epoch
    assert header.meta.location == epoch.location
    assert epoch.sid == "IV.TEST..HHZ"

    header.meta = None
    assert header.meta.empty
    assert header.sid == "IV.TEST..HHN"

    with pytest.raises(ValueError):
        Header(sid="IV.OTHER..HHZ", meta=epoch)
    with pytest.raises(TypeError):
        Header(meta={})
    with pytest.raises(ValueError):
        Header(sid="TEST")


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


def test_header_nsamp_from_parent():
    """nsamp is read from the parent record."""
    record = Record(data=np.zeros(25), delta=0.01)

    assert record.head.nsamp == 25


def test_header_copy_is_independent():
    """A copy shares only the channel epoch and drops the parent."""
    parent = object()
    header = Header(
        delta=0.01,
        time="2026-01-01T00:00:00",
        sid="IV.TEST..HHZ",
        parent=parent,
    )
    header.extra["quality"] = "D"

    copied = header.copy()

    assert copied is not header
    assert copied._parent is None
    assert copied.delta == header.delta
    assert copied.sid == header.sid
    assert copied.time == header.time
    assert copied.time is not header.time
    assert copied.extra == header.extra
    assert copied.extra is not header.extra
    assert copied.meta is header.meta

    copied.sid = "IV.TEST..HHN"
    assert header.sid == "IV.TEST..HHZ"

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
    copied.head.sid = "XX.X..HHZ"

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
        raw / mode_response.get("IV.MODE..HNZ"),
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
        "normalization_factor": 2.45957e13,
        "normalization_frequency": 1.0,
        "poles": [-981 - 1009j, -981 + 1009j, -3290 - 1263j, -3290 + 1263j],
        "zeros": [],
        "stage_number": 1,
    })
    gain_2 = StageGain({
        "description": None,
        "sensitivity": 315000,
        "frequency": 0,
        "stage_number": 2,
    })

    by_stages = raw / gain_1 / poles_zeros / gain_2

    analog = raw.copy()
    analog.deconvolve_response(mode_response, stages=("gain", "paz"))

    assert by_stages.peak == pytest.approx(analog.peak, rel=1e-6)


def test_record_response_peak_matches_sensitivity(mode, mode_response):
    """
    The corrected PGA agrees with the total sensitivity.

    The accelerometer response is flat well above the signal band, so
    dividing by the sensitivity of the StationXML (321707 counts per
    m/s^2) must give nearly the same peak as the correction of the
    analog stages. The complete correction, which also inverts the
    anti-alias filters near the Nyquist frequency, stays within 2 %.
    """
    raw = mode["IV.MODE..HNZ"][0]
    simple = raw.copy()
    simple.remove_sensitivity(mode_response)
    analog = raw.copy()
    analog.deconvolve_response(mode_response, stages=("gain", "paz"))
    full = raw / mode_response

    assert simple.peak == pytest.approx(raw.peak / 321707.0)
    assert analog.peak == pytest.approx(simple.peak, rel=0.01)
    assert full.peak == pytest.approx(simple.peak, rel=0.02)


# ----------------------------------------------------------------------------
# Record: channel metadata
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        (13.0, 46.0, 250.0),
        [13.0, 46.0, 250.0],
        np.array([13.0, 46.0, 250.0]),
    ],
)
def test_record_location_sequences(value):
    """Locations are (longitude, latitude[, elevation]) sequences."""
    record = Record(location=value)

    assert record.location == WgsPoint(
        longitude=13.0, latitude=46.0, elevation=250.0
    )
    assert record.head.meta.location == record.location


def test_record_location_point_is_kept():
    """A WgsPoint is stored as is."""
    point = WgsPoint(longitude=13.0, latitude=46.0)

    assert Record(location=point).location is point


def test_record_invalid_location():
    """Wrong lengths, strings and invalid latitudes are rejected."""
    with pytest.raises(ValueError):
        Record(location=(13.0,))
    with pytest.raises(TypeError):
        Record(location=13.0)
    with pytest.raises(ValueError):
        Record(location=(13.0, 95.0))


def test_record_metadata_attributes_replace_shared_epoch():
    """Assigning a record attribute never changes the copies."""
    record = Record(data=[1.0, 2.0], delta=1.0, time=START,
                    sid="XX.STA..HHZ", location=(13.0, 46.0))
    copied = record.copy()
    assert copied.head.meta is record.head.meta

    copied.location = (14.0, 47.0)
    copied.orientation = (90.0, 0.0)
    copied.response = StageSet(sample_rate=1.0)

    assert record.location == WgsPoint(13.0, 46.0)
    assert record.orientation is None
    assert record.response is None
    assert copied.orientation.azimuth == 90.0
    assert copied.response.sample_rate == 1.0


def test_record_set_metadata_from_collection(mode, mode_response):
    """The epoch valid at the start is stored as a copy."""
    raw = mode["IV.MODE..HNZ"][0].copy()
    raw.set_metadata(mode_response)
    stored = mode_response.get("IV.MODE..HNZ", raw.time)

    assert raw.head.meta == stored
    assert raw.head.meta is not stored
    assert raw.location == stored.location
    assert raw.orientation.dip == -90.0
    assert raw.response.sensitivity.value == pytest.approx(321707.0)


def test_record_set_metadata_from_epoch():
    """An epoch must belong to the channel and be valid at the start."""
    record = Record(data=[1.0], delta=1.0, time=START, sid="XX.STA..HHZ")

    record.set_metadata(ChannelEpoch(location=(13.0, 46.0)))
    assert record.sid == "XX.STA..HHZ"
    assert record.location == WgsPoint(13.0, 46.0)

    with pytest.raises(ValueError):
        record.set_metadata(ChannelEpoch("XX.OTHER..HHZ"))
    with pytest.raises(ValueError):
        record.set_metadata(
            ChannelEpoch("XX.STA..HHZ", starttime="2027-01-01")
        )
    with pytest.raises(TypeError):
        record.set_metadata({"location": (13.0, 46.0)})

    anonymous = Record(data=[1.0], delta=1.0, time=START)
    anonymous.set_metadata(ChannelEpoch("XX.STA..HHZ"))
    assert anonymous.sid == "XX.STA..HHZ"


def test_record_set_metadata_missing_channel(mode_response):
    """A record needs an identifier present in the collection."""
    with pytest.raises(KeyError):
        Record(data=[1.0], delta=0.01, time=START,
               sid="XX.STA..HHZ").set_metadata(mode_response)
    with pytest.raises(ValueError):
        Record(data=[1.0], delta=0.01, time=START).set_metadata(
            mode_response
        )


def test_record_set_metadata_warnings():
    """Replaced metadata, other rates and epoch ends are reported."""
    def record(**kwargs):
        return Record(data=np.zeros(10), delta=1.0, time=START,
                      sid="XX.STA..HHZ", **kwargs)

    with pytest.warns(MetadataWarning, match="replaced"):
        record(location=(13.0, 46.0)).set_metadata(
            ChannelEpoch("XX.STA..HHZ", location=(1.0, 1.0))
        )
    with pytest.warns(MetadataWarning, match="sampling rate"):
        record().set_metadata(ChannelEpoch("XX.STA..HHZ", sample_rate=100))
    with pytest.warns(MetadataWarning, match="beyond the end"):
        record().set_metadata(
            ChannelEpoch("XX.STA..HHZ", endtime="2026-01-01T00:00:05")
        )


def test_stream_set_metadata_shares_one_copy(mode, mode_response):
    """The records of one epoch share a single copy of it."""
    stream = mode["IV.MODE..HNZ"].copy()
    stream.append(
        stream[0].extract(100.0, 110.0)
    )
    assert len(stream) == 2

    stream.set_metadata(mode_response)

    assert stream[0].head.meta is stream[1].head.meta
    assert stream[0].head.meta is not mode_response.get(
        "IV.MODE..HNZ", stream[0].time
    )


def test_collection_set_metadata(mode, mode_response):
    """Channels without metadata are left unchanged, with a warning."""
    collection = mode.copy()
    collection.append(Record(data=[1.0], delta=0.01, time=START,
                             sid="XX.STA..HHZ"))

    collection.append(Record(data=[1.0], delta=0.01, time=START,
                             sid="XX.STA..HHN"))

    with pytest.warns(MetadataWarning, match="XX.STA..HHZ") as caught:
        collection.set_metadata(mode_response)

    # One warning for all the channels without metadata.
    assert len(caught) == 1
    assert "2 channel(s)" in str(caught[0].message)
    assert "XX.STA..HHN" in str(caught[0].message)
    assert collection["IV.MODE..HNE"][0].orientation.azimuth == 90.0
    assert collection["XX.STA..HHZ"][0].head.meta.empty


def test_record_response_from_metadata(mode, mode_response):
    """Without argument the response of the metadata is used."""
    raw = mode["IV.MODE..HNZ"][0].copy()
    expected = raw / mode_response

    raw.set_metadata(mode_response)
    corrected = raw.copy()
    corrected.deconvolve_response()
    np.testing.assert_allclose(corrected.data, expected.data)

    simple = raw.copy()
    simple.remove_sensitivity()
    assert simple.peak == pytest.approx(raw.peak / 321707.0)

    with pytest.raises(ValueError, match="no instrument response"):
        mode["IV.MODE..HNZ"][0].copy().deconvolve_response()


def test_record_response_twice_is_rejected(mode, mode_response):
    """Declared units prevent a second correction."""
    raw = mode["IV.MODE..HNZ"][0].copy()
    raw.set_metadata(mode_response)
    raw.deconvolve_response()
    assert raw.head.units.lower() == "m/s**2"

    with pytest.raises(ValueError, match="already corrected"):
        raw.deconvolve_response()
    with pytest.raises(ValueError, match="already corrected"):
        raw.remove_sensitivity()

    raw.convolve_response()
    assert raw.head.units.lower() in ("count", "counts")


def test_record_append_metadata():
    """Empty metadata join any; different metadata are rejected."""
    first = Record(data=[1.0, 2.0], delta=1.0, time=START,
                   sid="XX.STA..HHZ")
    second = Record(data=[3.0], delta=1.0, time="2026-01-01T00:00:02",
                    sid="XX.STA..HHZ", location=(13.0, 46.0))

    assert first.append(second) is True
    assert first.location == WgsPoint(13.0, 46.0)

    third = Record(data=[4.0], delta=1.0, time="2026-01-01T00:00:03",
                   sid="XX.STA..HHZ", location=(14.0, 46.0))
    with pytest.raises(ValueError, match="channel metadata"):
        first.append(third)
    assert first.nsamp == 3


def test_stream_keeps_epochs_separate():
    """Contiguous records of different epochs remain separate."""
    stream = Stream()
    stream.append(Record(data=[1.0, 2.0], delta=1.0, time=START,
                         sid="XX.STA..HHZ", location=(13.0, 46.0)))
    stream.append(Record(data=[3.0], delta=1.0, time="2026-01-01T00:00:02",
                         sid="XX.STA..HHZ", location=(14.0, 46.0)))

    assert len(stream) == 2
    assert stream[1].location == WgsPoint(14.0, 46.0)


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
    assert SNAP_TOLERANCE == 1e-6


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
    record = Record(data=[1], delta=1.0, time=START, sid="XX.A..HHZ")

    with pytest.raises(ValueError):
        record.append(
            Record(data=[2], delta=1.0, time="2026-01-01T00:00:01",
                   sid="XX.B..HHZ")
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
    other = Record(data=[1, 2], delta=1.0, time=START, sid="XX.A..HHZ")

    assert record.append(other) is True
    assert record.sid == "XX.A..HHZ"
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


# ----------------------------------------------------------------------------
# Stream
# ----------------------------------------------------------------------------


def _two_record_stream():
    """Stream with records at 0-9 s and 20-29 s (1 Hz)."""
    stream = Stream("XX.STA..HHZ")
    stream.append(Record(data=np.arange(10), delta=1.0, time=START))
    stream.append(
        Record(
            data=np.arange(100, 110),
            delta=1.0,
            time="2026-01-01T00:00:20",
        )
    )
    return stream


def test_stream_append_demo_does_not_modify_caller_records():
    """
    Reproduces demo 00_basic_use: contiguous records are merged.

    The records passed by the caller keep their original length and
    source identifier.
    """
    pulse = np.zeros(10001)
    pulse[5001] = 1.0
    rec0 = Record(data=pulse, delta=0.001, time="2024-06-20T12:00:00",
                  sid="XX.DIR01..HHZ")
    rec1 = Record(data=pulse, delta=0.001, time="2024-06-20T12:00:10.001",
                  sid="XX.DIR01..HHZ")
    late = Record(data=pulse, delta=0.001, time="2024-06-20T15:00:00",
                  sid="XX.DIR01..HHZ")

    stream = Stream("XX.DIR01..HHZ")
    stream.append(rec0)
    stream.append(rec1)
    stream.append(late)

    assert len(stream) == 2
    assert stream[0].nsamp == 20002
    assert stream[0] is not rec0
    assert rec0.nsamp == rec1.nsamp == 10001


def test_stream_append_sets_sid_on_stored_copy():
    """A record without identifier gets the stream sid; the original not."""
    record = Record(data=[1, 2], delta=1.0, time=START)
    stream = Stream("XX.STA..HHZ")
    stream.append(record)

    assert stream[0].sid == "XX.STA..HHZ"
    assert record.sid is None


def test_stream_append_rejects_other_sources():
    """Records of another source cannot be appended."""
    stream = Stream("XX.A..HHZ")

    with pytest.raises(ValueError):
        stream.append(Record(data=[1], delta=1.0, time=START, sid="XX.B..HHZ"))
    with pytest.raises(TypeError):
        stream.append([1, 2])


def test_stream_properties_and_text(capsys):
    """Size, times, identifiers and text representations."""
    stream = _two_record_stream()

    assert len(stream) == 2
    assert stream.nsamp == 20
    assert stream.starttime == Date(START)
    assert stream.endtime == Date("2026-01-01T00:00:29")
    assert stream.eid == [None, None]
    assert repr(stream) == "Stream(sid='XX.STA..HHZ', records=2)"
    assert "number of records: 2" in str(stream)

    stream.info()
    assert "Record 1" in capsys.readouterr().out
    assert Stream().starttime is None


def test_stream_indexing():
    """Integers (NumPy included) select records; booleans are rejected."""
    stream = _two_record_stream()

    assert stream[0].data[0] == 0
    assert stream[np.int64(1)].data[0] == 100
    assert stream[-1].data[0] == 100

    with pytest.raises(TypeError):
        stream[True]
    with pytest.raises(TypeError):
        stream[1.0]


def test_stream_slice_returns_stream_of_stored_records():
    """A slice is a Stream with the same sid and the stored records."""
    stream = _two_record_stream()
    selection = stream[0:1]

    assert isinstance(selection, Stream)
    assert selection.sid == stream.sid
    assert len(selection) == 1
    assert selection[0] is stream[0]
    assert selection.copy()[0] is not stream[0]


def test_stream_event_identifiers():
    """Records can be accessed and removed by event identifier."""
    stream = Stream("XX.A..HHZ")
    stream.append(Record(data=[1], delta=1.0, time=START, eid="ev1"))
    stream.append(
        Record(data=[2], delta=1.0, time="2026-01-01T01:00:00", eid="ev2")
    )

    assert stream["ev2"].data[0] == 2
    assert stream.remove("ev1").data[0] == 1
    assert len(stream) == 1

    with pytest.raises(KeyError):
        stream["missing"]


def test_stream_remove_and_sort():
    """Records are removed by position and sorted by start time."""
    stream = _two_record_stream()
    stream.sort(reverse=True)
    assert stream[0].data[0] == 100

    stream.sort()
    assert stream[0].data[0] == 0

    removed = stream.remove(np.int64(0))
    assert removed.data[0] == 0
    assert len(stream) == 1

    with pytest.raises(TypeError):
        stream.remove(True)


def test_stream_copy_is_independent():
    """A copied stream holds independent records."""
    stream = _two_record_stream()
    copied = stream.copy()
    copied[0].data[0] = 99

    assert stream[0].data[0] == 0
    assert copied.sid == stream.sid


def test_stream_extract_relative_to_stream_start():
    """Numerical times are seconds from the start of the stream."""
    stream = _two_record_stream()

    early = stream.extract(2, 4)
    late = stream.extract(22, 24)

    assert [r.data.tolist() for r in early] == [[2, 3, 4]]
    assert [r.data.tolist() for r in late] == [[102, 103, 104]]
    assert late[0].time == Date("2026-01-01T00:00:22")


def test_stream_extract_across_records_and_merge():
    """Segments of several records can be merged when requested."""
    stream = _two_record_stream()
    window = stream.extract(8, 21)

    assert [r.data.tolist() for r in window] == [[8, 9], [100, 101]]

    with pytest.raises(ValueError):
        stream.extract(8, 21, merge=True)

    merged = stream.extract(8, 21, merge=True, enforce=True)
    assert merged.data.tolist() == [8, 9] + [0] * 10 + [100, 101]


def test_stream_extract_skips_records_without_samples():
    """Windows between samples or between records return nothing."""
    stream = _two_record_stream()

    assert len(stream.extract(2.2, 2.8)) == 0
    assert len(stream.extract(12, 18)) == 0
    assert stream.extract(12, 18, merge=True) is None

    with pytest.raises(ValueError):
        stream.extract(5, 3)


def test_stream_extract_absolute_times(mode):
    """Absolute windows work on real data."""
    stream = mode["IV.MODE..HNZ"]
    window = stream.extract(
        "2012-05-20T02:03:33", "2012-05-20T02:03:43"
    )

    assert window[0].nsamp == 1001
    assert window[0].time == Date("2012-05-20T02:03:33")
    assert stream.extract(10, 20)[0].time == window[0].time


def test_stream_processing_applies_to_every_record(mode):
    """Stream processing methods match the record methods."""
    stream = mode["IV.MODE..HNZ"].copy()
    record = stream[0].copy()

    stream.remove_mean()
    stream.filter(0.1, 20.0)
    stream.taper(1.0)
    record.remove_mean()
    record.filter(0.1, 20.0)
    record.taper(1.0)
    np.testing.assert_array_equal(stream[0].data, record.data)

    stream.resample(0.02)
    stream.integrate(method="cum")
    stream.differentiate(method="grad")
    stream.shift(0.5)
    assert stream[0].delta == 0.02
    assert stream[0].nsamp == 9001


def test_stream_response_correction(mode, mode_response):
    """Stream response removal matches the record correction."""
    stream = mode["IV.MODE..HNZ"].copy()
    expected = stream[0] / mode_response

    stream.deconvolve_response(mode_response)

    np.testing.assert_allclose(stream[0].data, expected.data)


def test_record_stream_response_of_another_channel_is_rejected(
    mode, mode_response
):
    """A stream response is applied only to records of its own channel."""
    east = mode["IV.MODE..HNE"][0]

    with pytest.raises(ValueError):
        east / mode_response.get("IV.MODE..HNZ")


# ----------------------------------------------------------------------------
# StreamCollection
# ----------------------------------------------------------------------------


def test_collection_read_and_access(mode, capsys):
    """Reading MODE.mseed gives three streams accessible by sid."""
    assert len(mode) == 3
    assert mode.sid == ["IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ"]
    assert mode.nsamp == 3 * 18001
    assert mode["IV.MODE..HNZ"] is mode.get_stream("IV.MODE..HNZ")
    assert mode[np.int64(2)] is mode["IV.MODE..HNZ"]
    assert repr(mode) == "StreamCollection(streams=3)"

    mode.info()
    assert "number of streams: 3" in capsys.readouterr().out

    with pytest.raises(TypeError):
        mode[True]
    with pytest.raises(KeyError):
        mode["IV.MODE..XXX"]


def test_collection_read_method(tmp_path):
    """StreamCollection.read appends the content of a file."""
    collection = StreamCollection()
    returned = collection.read(str(DATA / "MODE.mseed"))

    assert returned is collection
    assert len(collection) == 3


def test_collection_slice_returns_collection_of_stored_streams(mode):
    """A slice is a collection holding the stored streams."""
    selection = mode[0:2]

    assert isinstance(selection, StreamCollection)
    assert len(selection) == 2
    assert selection[0] is mode[0]


def test_collection_append_copies_records():
    """Appended records are copied and grouped by source identifier."""
    collection = StreamCollection()
    record = Record(data=[1, 2], delta=1.0, time=START, sid="XX.A..HHZ")

    collection.append(record)
    collection.append(
        Record(data=[3], delta=1.0, time="2026-01-01T00:00:02",
               sid="XX.A..HHZ")
    )
    collection.append(Record(data=[5], delta=1.0, time=START, sid="XX.B..HHZ"))

    assert collection.sid == ["XX.A..HHZ", "XX.B..HHZ"]
    assert collection["XX.A..HHZ"][0].data.tolist() == [1, 2, 3]
    assert record.data.tolist() == [1, 2]

    with pytest.raises(ValueError):
        collection.append(Record(data=[1], delta=1.0, time=START))
    with pytest.raises(TypeError):
        collection.append([1, 2])


def test_collection_remove_sort_copy_merge():
    """Streams are removed, sorted, copied and merged."""
    collection = StreamCollection()
    for sid in ("XX.C..HHZ", "XX.A..HHZ", "XX.B..HHZ"):
        collection.append(Record(data=[1], delta=1.0, time=START, sid=sid))

    collection.sort()
    assert collection.sid == ["XX.A..HHZ", "XX.B..HHZ", "XX.C..HHZ"]

    copied = collection.copy()
    copied["XX.A..HHZ"][0].data[0] = 9
    assert collection["XX.A..HHZ"][0].data[0] == 1

    assert collection.remove(np.int64(0)).sid == "XX.A..HHZ"
    assert collection.remove("XX.C..HHZ").sid == "XX.C..HHZ"
    assert collection.sid == ["XX.B..HHZ"]

    other = StreamCollection()
    other.append(Record(data=[2], delta=1.0, time=START, sid="XX.D..HHZ"))
    collection.merge(other)
    assert collection.sid == ["XX.B..HHZ", "XX.D..HHZ"]


def test_collection_extract_one_or_all_streams(mode):
    """Extraction from one stream, or from every stream."""
    one = mode.extract("IV.MODE..HNZ", 10, 20)
    assert isinstance(one, Stream)
    assert one[0].nsamp == 1001

    every = mode.extract(starttime=10, endtime=20)
    assert isinstance(every, StreamCollection)
    assert every.sid == mode.sid
    assert all(stream[0].nsamp == 1001 for stream in every)

    merged = mode.extract(starttime=10, endtime=20, merge=True)
    assert [len(stream) for stream in merged] == [1, 1, 1]

    empty = mode.extract(
        starttime="2013-01-01T00:00:00", endtime="2013-01-01T00:01:00"
    )
    assert len(empty) == 0


def test_collection_response_correction(mode, mode_response):
    """Collection response removal matches the record correction."""
    corrected = mode / mode_response
    expected = mode["IV.MODE..HNN"][0] / mode_response

    np.testing.assert_allclose(
        corrected["IV.MODE..HNN"][0].data, expected.data
    )
    assert mode["IV.MODE..HNN"][0].data.dtype.kind != "c"

    with pytest.raises(ValueError):
        mode / mode_response.get("IV.MODE..HNZ")
    with pytest.raises(TypeError):
        mode / 2


def test_collection_processing_applies_to_every_stream(mode):
    """Collection processing matches the record methods."""
    collection = mode.copy()
    record = mode["IV.MODE..HNE"][0].copy()

    collection.filter(0.1, 20.0)
    collection.taper(1.0)
    record.filter(0.1, 20.0)
    record.taper(1.0)

    np.testing.assert_array_equal(
        collection["IV.MODE..HNE"][0].data, record.data
    )

    collection.remove_mean()
    collection.resample(0.02)
    collection.integrate(method="cum")
    collection.differentiate(method="grad")
    collection.shift(0.5)
    assert all(stream[0].delta == 0.02 for stream in collection)


def test_collection_max_duration_validation():
    """max_duration must be a positive finite number or None."""
    assert StreamCollection().max_duration is None

    with pytest.raises(TypeError):
        StreamCollection(max_duration="10")
    with pytest.raises(ValueError):
        StreamCollection(max_duration=-1.0)
