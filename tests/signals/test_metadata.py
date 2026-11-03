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
Tests for :mod:`shakelab.signals.metadata`.

Data: ``MODE.xml`` (INGV, one station, three channels) and
``TT_T01.xml`` (dummy station with 12 channels and 30 epochs, written
with ObsPy; see the README of the data folder).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.signals.metadata import (
    ChannelEpoch,
    Equipment,
    MetadataCollection,
    MetadataWarning,
    Orientation,
    StationEpoch,
    select_response,
    select_sensitivity,
)
from shakelab.signals.response import (
    StageCoefficients,
    StageGain,
    StageSet,
)
from shakelab.signals.stationxml import read_stationxml

DATA = Path(__file__).parent / "data"
MODE = DATA / "MODE.xml"
TT = DATA / "TT_T01.xml"
SID = "IV.MODE..HNZ"
TIME = "2012-05-20T02:03:23"


@pytest.fixture(scope="module")
def mode():
    """Metadata of station IV.MODE."""
    return MetadataCollection(str(MODE))


@pytest.fixture(scope="module")
def tt():
    """Metadata of the dummy station TT.T01."""
    return MetadataCollection(str(TT))


def _epoch(sid="XX.TEST..HHZ", start=None, end=None, gain=1000.0):
    """Channel epoch with a single gain stage."""
    stage_set = StageSet(sample_rate=100.0)
    stage_set.append(StageGain(sensitivity=gain, frequency=1.0,
                               stage_number=1))
    return ChannelEpoch(sid, starttime=start, endtime=end,
                        sample_rate=100.0, response=stage_set)


# ----------------------------------------------------------------------------
# Value objects
# ----------------------------------------------------------------------------


def test_orientation():
    """Azimuth is reduced to [0, 360); the vector is (E, N, U)."""
    assert Orientation(360.0, 0.0).azimuth == 0.0
    assert Orientation(-90.0, 0.0).azimuth == 270.0

    np.testing.assert_allclose(Orientation(0.0, -90.0).vector, [0, 0, 1],
                               atol=1e-15)
    np.testing.assert_allclose(Orientation(90.0, 0.0).vector, [1, 0, 0],
                               atol=1e-15)
    np.testing.assert_allclose(Orientation(0.0, 0.0).vector, [0, 1, 0],
                               atol=1e-15)
    np.testing.assert_allclose(Orientation(0.0, 90.0).vector, [0, 0, -1],
                               atol=1e-15)


def test_orientation_errors():
    """Invalid orientations are rejected."""
    with pytest.raises(ValueError):
        Orientation(0.0, 91.0)
    with pytest.raises(ValueError):
        Orientation(float("nan"), 0.0)
    with pytest.raises(TypeError):
        Orientation("north", 0.0)


def test_equipment_text():
    """Equipment prints manufacturer and model, or the description."""
    assert str(Equipment(manufacturer="Kinemetrics",
                         model="Episensor")) == "Kinemetrics Episensor"
    assert str(Equipment(description="sensor X")) == "sensor X"
    assert str(Equipment()) == ""


def test_channel_epoch_validation():
    """Codes, coordinates, orientation and rates are checked."""
    epoch = ChannelEpoch("XX.TEST..HHZ", location=(13.0, 46.0, 100.0),
                         orientation=(0.0, -90.0), depth=0.01)
    assert epoch.location == WgsPoint(13.0, 46.0, 100.0)
    assert epoch.orientation == Orientation(0.0, -90.0)
    assert (epoch.network, epoch.station_code, epoch.location_code,
            epoch.channel) == ("XX", "XX.TEST", "", "HHZ")

    with pytest.raises(ValueError):
        ChannelEpoch("XX.TEST.HHZ")
    with pytest.raises(ValueError):
        StationEpoch("XX.TEST.00")
    with pytest.raises(TypeError):
        epoch.orientation = 10.0
    with pytest.raises(ValueError):
        epoch.location = (13.0,)
    with pytest.raises(ValueError):
        epoch.sample_rate = -1.0
    with pytest.raises(TypeError):
        epoch.response = "not a response"


def test_epoch_interval():
    """Start included, end excluded, open bounds."""
    epoch = _epoch(start="2020-01-01T00:00:00", end="2021-01-01T00:00:00")
    assert epoch.match("2020-01-01T00:00:00")
    assert not epoch.match("2021-01-01T00:00:00")
    assert _epoch().match("1900-01-01T00:00:00")

    later = _epoch(start="2021-01-01T00:00:00")
    assert not epoch.overlaps(later)
    assert later.overlaps(_epoch())


# ----------------------------------------------------------------------------
# StationXML
# ----------------------------------------------------------------------------


def test_read_mode(mode):
    """Station, coordinates, orientation and equipment of IV.MODE."""
    assert mode.sid == ["IV.MODE..HNZ", "IV.MODE..HNN", "IV.MODE..HNE"]
    assert len(mode) == 3

    station, = mode.stations()
    assert station.code == "IV.MODE"
    assert station.site == "Modena"
    assert station.location == WgsPoint(10.949167, 44.629722, 41.0)

    epoch = mode.get(SID, TIME)
    assert epoch.station is station
    assert epoch.location == station.location
    assert epoch.depth == 0.0
    assert epoch.orientation == Orientation(0.0, -90.0)
    assert epoch.sample_rate == 100.0
    assert str(epoch.sensor) == "Kinemetrics Episensor-ES-T"
    assert epoch.datalogger.model == "GAIA2"
    assert epoch.response.sensitivity.value == 321707.0

    assert mode.get("IV.MODE..HNE", TIME).orientation.azimuth == 90.0


def test_read_tt_epochs(tt):
    """Twelve channels, thirty epochs, depth converted to km."""
    assert len(tt.sid) == 12
    assert len(tt) == 30
    assert len(tt["TT.T01.01.HHE"]) == 6

    epoch = tt.get("TT.T01.01.HHE", "2019-10-08T18:30:00")
    assert str(epoch.starttime) == "2019-10-08T12:30:00.000000Z"
    assert epoch.endtime is None
    assert epoch.depth == pytest.approx(0.005)
    assert epoch.orientation.azimuth == 93.0
    assert epoch.sample_rate == 200.0

    early = tt.get("TT.T01.01.HHN", "2012-06-01T00:00:00")
    assert early.orientation.azimuth == 3.0
    assert tt.validate() == []


def test_read_ftypes():
    """StationXML is recognized by name or by an explicit type."""
    for ftype in (None, "xml", "sxml"):
        assert MetadataCollection(str(MODE), ftype=ftype).sid == [
            "IV.MODE..HNZ", "IV.MODE..HNN", "IV.MODE..HNE"
        ]
    with pytest.raises(ValueError):
        MetadataCollection(str(MODE), ftype="txt")


def test_channel_without_coordinates_takes_station():
    """A channel without coordinates inherits those of its station."""
    document = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" '
        'schemaVersion="1.1"><Network code="XX">'
        '<Station code="TEST" startDate="2020-01-01T00:00:00">'
        "<Latitude>46.0</Latitude><Longitude>13.0</Longitude>"
        "<Elevation>200</Elevation>"
        '<Channel code="HHZ" locationCode="00" '
        'startDate="2020-01-01T00:00:00"><Depth>15</Depth></Channel>'
        '<Channel code="HHE" locationCode="00" '
        'startDate="2020-01-01T00:00:00"><Azimuth>90</Azimuth>'
        "</Channel>"
        "</Station></Network></FDSNStationXML>"
    )
    collection = read_stationxml(document)
    vertical = collection.get("XX.TEST.00.HHZ")
    assert vertical.location == WgsPoint(13.0, 46.0, 200.0)
    assert vertical.depth == pytest.approx(0.015)
    assert vertical.orientation is None
    assert vertical.response is None
    # Azimuth without dip gives no orientation.
    assert collection.get("XX.TEST.00.HHE").orientation is None


# ----------------------------------------------------------------------------
# Collection
# ----------------------------------------------------------------------------


def test_get_and_errors(tt):
    """Epochs are selected by identifier and time."""
    with pytest.raises(KeyError):
        tt.get("XX.NONE..HHZ", TIME)
    with pytest.raises(KeyError):
        tt.get("TT.T01.01.HHE", "2000-01-01T00:00:00")
    with pytest.raises(ValueError):
        tt.get("TT.T01.01.HHE")

    # Gap between two epochs.
    with pytest.raises(KeyError):
        tt.get("TT.T01.01.HHE", "2012-01-25T12:00:00")


def test_select(tt):
    """Wildcards and time give a sub-collection of copies."""
    strong = tt.select("TT.T01.02.HN?")
    assert sorted(strong.sid) == [
        "TT.T01.02.HNE", "TT.T01.02.HNN", "TT.T01.02.HNZ"
    ]
    assert len(strong) == 6

    now = tt.select(time="2020-01-01T00:00:00")
    assert len(now) == 6
    assert all(epoch.endtime is None for epoch in now)

    epoch = strong.get("TT.T01.02.HNZ", "2017-01-01T00:00:00")
    assert epoch is not tt.get("TT.T01.02.HNZ", "2017-01-01T00:00:00")
    assert epoch.station is strong.stations()[0]


def test_stations(tt, mode):
    """Stations are shared by their channels and filtered by code."""
    merged = MetadataCollection()
    merged.extend(tt)
    merged.extend(mode)

    assert [s.code for s in merged.stations()] == ["TT.T01", "IV.MODE"]
    assert [s.code for s in merged.stations("IV.*")] == ["IV.MODE"]
    assert merged.stations(time="2000-01-01T00:00:00") == []

    stations = {id(epoch.station) for epoch in merged.select("TT.*")}
    assert len(stations) == 1


def test_append_stores_copies(mode):
    """Appended epochs are copies; equal epochs are ignored."""
    collection = MetadataCollection()
    original = mode.get(SID, TIME)
    collection.append(original)
    collection.append(original)

    assert len(collection) == 1
    stored = collection.get(SID, TIME)
    assert stored is not original
    assert stored.response is not original.response
    assert stored.station is not original.station
    assert stored == original

    with pytest.raises(TypeError):
        collection.append(StageSet())


def test_overlapping_epochs_warn():
    """Different overlapping epochs are kept, with a warning."""
    collection = MetadataCollection()
    collection.append(_epoch(start="2020-01-01T00:00:00"))
    with pytest.warns(MetadataWarning):
        collection.append(_epoch(start="2021-01-01T00:00:00", gain=2.0))

    assert len(collection) == 2
    assert len(collection.validate()) == 1

    with pytest.warns(MetadataWarning):
        epoch = collection.get("XX.TEST..HHZ", "2022-01-01T00:00:00")
    assert epoch.response.stage[0].sensitivity == 1000.0


def test_epochs_sorted_by_start():
    """Epochs of a channel are kept in order of start time."""
    collection = MetadataCollection()
    collection.append(_epoch(start="2021-01-01T00:00:00"))
    collection.append(_epoch(end="2021-01-01T00:00:00", gain=2.0))
    first, second = collection["XX.TEST..HHZ"]
    assert first.starttime is None
    assert str(second.starttime) == "2021-01-01T00:00:00.000000Z"


def test_validate_sample_rate():
    """A declared rate different from the response is reported."""
    epoch = _epoch()
    epoch.sample_rate = 50.0
    collection = MetadataCollection()
    collection.append(epoch)
    issues = collection.validate()
    assert len(issues) == 1
    assert "sampling rate" in issues[0]


def test_json_round_trip(tmp_path, tt, mode):
    """JSON keeps stations, epochs, metadata and responses."""
    collection = MetadataCollection()
    collection.extend(tt)
    collection.extend(mode)
    collection.get(SID, TIME).extra["quality"] = "D"

    path = tmp_path / "metadata.json"
    collection.write(str(path))
    loaded = MetadataCollection(str(path))

    assert loaded.sid == collection.sid
    assert len(loaded) == len(collection)
    assert len(loaded.stations()) == 2
    for one, two in zip(loaded, collection):
        assert one == two

    epoch = loaded.get(SID, TIME)
    assert epoch.extra == {"quality": "D"}
    assert epoch.station is loaded.stations("IV.MODE")[0]
    assert epoch.response.stage[4].input_sample_rate == 800.0

    freq = np.array([0.5, 1.0, 10.0])
    np.testing.assert_allclose(
        epoch.response.transfer_function(freq),
        mode.get(SID, TIME).response.transfer_function(freq),
    )


def test_json_coefficients_stage(tmp_path):
    """A JSON 'coefficients' stage is read as digital coefficients."""
    path = tmp_path / "metadata.json"
    path.write_text(
        '{"channels": [{"sid": "XX.TEST..HHZ", "response": {"stages": '
        '[{"stage_number": 1, "stage_type": "coefficients", '
        '"data": {"numerator": [0.5, 0.5]}}]}}]}'
    )
    loaded = MetadataCollection(str(path))
    stage = loaded.get("XX.TEST..HHZ").response.stage[0]
    assert isinstance(stage, StageCoefficients)


def test_info_runs(tt, capsys):
    """The summaries print without errors."""
    tt.info()
    tt.get("TT.T01.01.HHE", "2019-10-08T18:30:00").info()
    output = capsys.readouterr().out
    assert "TT.T01 (Test)" in output
    assert "azimuth 93.0" in output


def test_reading_is_silent():
    """Complete documents are read without warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        MetadataCollection(str(MODE))
        MetadataCollection(str(TT))


# ----------------------------------------------------------------------------
# Response selection
# ----------------------------------------------------------------------------


def test_select_response(mode):
    """Responses are selected from collections and epochs."""
    epoch = mode.get(SID, TIME)
    expected = epoch.response

    assert select_response(mode, SID, TIME) is expected
    assert select_response(epoch, SID, TIME) is expected
    assert select_response(epoch) is expected
    assert select_response(expected) is expected

    stage = expected.stage[0]
    assert select_response(stage, "any", None) is stage


def test_select_response_errors(mode):
    """Invalid selections raise explicit errors."""
    epoch = mode.get(SID, TIME)

    with pytest.raises(ValueError):
        select_response(epoch, "IV.MODE..HNE", TIME)
    with pytest.raises(ValueError):
        select_response(epoch, SID, "2000-01-01T00:00:00")
    with pytest.raises(ValueError):
        select_response(mode, None, TIME)
    with pytest.raises(ValueError):
        select_response(ChannelEpoch(SID), SID, TIME)
    with pytest.raises(KeyError):
        select_response(mode, "XX.NONE..HHZ", TIME)
    with pytest.raises(TypeError):
        select_response(1.0, SID, TIME)


def test_select_sensitivity(mode):
    """The overall sensitivity is selected like the response."""
    sens = select_sensitivity(mode, SID, TIME)
    assert sens.value == 321707.0
    assert sens.input_units == "m/s**2"

    gain = select_sensitivity(StageGain(sensitivity=5.0, frequency=1.0))
    assert gain.value == 5.0

    with pytest.raises(ValueError):
        select_sensitivity(StageSet())
