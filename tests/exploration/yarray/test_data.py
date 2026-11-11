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
Unit tests for :mod:`shakelab.exploration.yarray.data`.

The data model is checked on synthetic records (projections onto
rotated components, time-axis and unit checks, coordinate
covariances) and the builders on synthetic collections with
geographic positions (ENU coordinates against
``geodeticN.geo_to_enu``, time windows, gaps, selections, metadata
orientations) and on an excerpt of the SESAME array simulation.
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.data import (
    ArrayData,
    ArrayDataWarning,
    StationData,
)
from shakelab.libutils.geodeticN import MetricFrame, WgsPoint, geo_to_enu
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.io import reader

DATA = Path(__file__).parent / "data" / "sesame"
START = "2026-01-01T00:00:00"
DELTA = 0.01


def _record(sid, data=None, nsamp=500, start=START, delta=DELTA, **meta):
    """Record with identifier and optional metadata."""
    if data is None:
        data = np.random.default_rng(len(sid)).standard_normal(nsamp)

    record = Record(data=data, delta=delta, time=start, sid=sid)

    for name, value in meta.items():
        setattr(record, name, value)

    return record


def _station_records(station, location=None, start=START, nsamp=500,
                     codes="ZNE", instrument="HH"):
    """Three records of one station."""
    return [
        _record(f"XX.{station}..{instrument}{code}", start=start,
                nsamp=nsamp, location=location)
        if location is not None
        else _record(f"XX.{station}..{instrument}{code}", start=start,
                     nsamp=nsamp)
        for code in codes
    ]


def _collection(records):
    """StreamCollection with the given records."""
    collection = StreamCollection()

    for record in records:
        collection.append(record)

    return collection


LOCATIONS = {
    "A01": (13.2500, 46.0600, 120.0),
    "A02": (13.2520, 46.0600, 118.0),
    "A03": (13.2510, 46.0615, 125.0),
}


def _array_collection(**kwargs):
    """Collection of three 3-C stations with geographic locations."""
    records = []

    for station, location in LOCATIONS.items():
        records += _station_records(station, location=location, **kwargs)

    return _collection(records)


# ----------------------------------------------------------------------------
# StationData
# ----------------------------------------------------------------------------


def test_station_orientation_from_channel_codes():
    """Without metadata Z, N and E give the canonical directions."""
    station = StationData("XX.A..", _station_records("A"), (0.0, 0.0))

    np.testing.assert_array_equal(
        station.orientation, [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    )
    np.testing.assert_array_equal(station.coordinates, [0.0, 0.0, 0.0])
    assert station.ncomp == 3 and station.nsamp == 500


def test_station_orientation_from_metadata():
    """Orientations in the metadata (azimuth, dip) are used."""
    records = [
        _record("XX.A..HH1", orientation=(30.0, 0.0)),
        _record("XX.A..HH2", orientation=(120.0, 0.0)),
        _record("XX.A..HHZ", orientation=(0.0, -90.0)),
    ]
    station = StationData("XX.A..", records, (0.0, 0.0))
    angle = math.radians(30.0)

    np.testing.assert_allclose(
        station.orientation[0], [math.sin(angle), math.cos(angle), 0.0],
        atol=1e-15,
    )
    np.testing.assert_allclose(station.orientation[2], [0, 0, 1], atol=1e-15)


def test_station_unknown_orientation():
    """Channels 1 and 2 without metadata have no orientation."""
    with pytest.raises(ValueError, match="orientation of channel"):
        StationData("XX.A..", [_record("XX.A..HH1")], (0.0, 0.0))


def test_rotated_components_reconstruct_east_north():
    """
    Horizontal components rotated by 30 degrees give back E and N.

    The rotated records are the projections of known east and north
    motions on the azimuths 30 and 120 degrees.
    """
    rng = np.random.default_rng(1)
    east, north, up = rng.standard_normal((3, 400))
    angles = np.radians([30.0, 120.0])
    records = [
        _record(f"XX.A..HH{code}",
                data=east * math.sin(a) + north * math.cos(a),
                orientation=(math.degrees(a), 0.0))
        for code, a in zip("12", angles)
    ] + [_record("XX.A..HHZ", data=up)]
    station = StationData("XX.A..", records, (0.0, 0.0))

    np.testing.assert_allclose(station.canonical(), [east, north, up],
                               atol=1e-12)
    np.testing.assert_allclose(station.radial(90.0), east, atol=1e-12)
    np.testing.assert_allclose(station.radial(0.0), north, atol=1e-12)


def test_transverse_convention():
    """The transverse direction is the radial one plus 90 degrees."""
    station = StationData("XX.A..", _station_records("A"), (0.0, 0.0))

    for azimuth in (0.0, 35.0, 200.0):
        np.testing.assert_allclose(
            station.transverse(azimuth), station.radial(azimuth + 90.0),
            atol=1e-12,
        )


def test_vertical_station_does_not_support_horizontal():
    """A vertical sensor cannot reconstruct horizontal motion."""
    station = StationData("XX.A..", [_record("XX.A..HHZ")], (0.0, 0.0))

    assert station.supports((0.0, 0.0, 1.0))
    assert not station.supports((1.0, 0.0, 0.0))

    with pytest.raises(ValueError, match="cannot be reconstructed"):
        station.radial(10.0)


@pytest.mark.parametrize(
    "value, expected",
    [
        (0.25, np.diag([0.25, 0.25, 0.0])),
        ([[1.0, 0.1], [0.1, 2.0]], [[1, 0.1, 0], [0.1, 2, 0], [0, 0, 0]]),
        (np.diag([1.0, 2.0, 3.0]), np.diag([1.0, 2.0, 3.0])),
    ],
)
def test_coordinate_covariance(value, expected):
    """Scalars, 2 x 2 and 3 x 3 covariances are accepted."""
    station = StationData(
        "XX.A..", _station_records("A"), (0.0, 0.0),
        coordinate_covariance=value,
    )

    np.testing.assert_allclose(station.coordinate_covariance, expected)


@pytest.mark.parametrize(
    "value", [[[1.0, 2.0], [0.0, 1.0]], [[-1.0, 0.0], [0.0, 1.0]], [1, 2]]
)
def test_invalid_coordinate_covariance(value):
    """Asymmetric, negative or malformed covariances are rejected."""
    with pytest.raises(ValueError):
        StationData("XX.A..", _station_records("A"), (0.0, 0.0),
                    coordinate_covariance=value)


def test_station_record_checks():
    """Records must be one to three, finite, on one time axis."""
    records = _station_records("A")

    with pytest.raises(ValueError, match="one to three"):
        StationData("XX.A..", records + records[:1], (0.0, 0.0))
    with pytest.raises(ValueError, match="different time axes"):
        StationData("XX.A..", [records[0], _record("XX.A..HHN",
                                                   start=Date(START) + 0.01)],
                    (0.0, 0.0))
    with pytest.raises(ValueError, match="non-finite"):
        StationData("XX.A..", [_record("XX.A..HHZ", data=[1.0, np.nan])],
                    (0.0, 0.0))
    with pytest.raises(TypeError):
        StationData("XX.A..", [np.zeros(10)], (0.0, 0.0))


def test_station_units():
    """Records in different units are rejected."""
    first, second = _station_records("A")[:2]
    first.head.units = "m/s"

    assert StationData("XX.A..", [first, second], (0.0, 0.0)).units == "m/s"

    second.head.units = "m/s**2"

    with pytest.raises(ValueError, match="different units"):
        StationData("XX.A..", [first, second], (0.0, 0.0))


def test_station_is_immutable():
    """Attributes and arrays of a station cannot be changed."""
    station = StationData("XX.A..", _station_records("A"), (1.0, 2.0))

    with pytest.raises(AttributeError):
        station.coordinates = (0.0, 0.0)
    with pytest.raises(ValueError):
        station.coordinates[0] = 5.0


# ----------------------------------------------------------------------------
# ArrayData
# ----------------------------------------------------------------------------


def _array(nsta=3):
    """Array of 3-C stations in a local frame."""
    return ArrayData([
        StationData(f"XX.A{i}..", _station_records(f"A{i}"), (10.0 * i, 0))
        for i in range(nsta)
    ])


def test_array_properties():
    """Shapes, time axis and access by identifier."""
    array = _array()

    assert array.nsta == 3 and len(array) == 3
    assert array.sids == ("XX.A0..", "XX.A1..", "XX.A2..")
    assert array.nsamp == 500 and array.delta == DELTA
    assert array.starttime == Date(START)
    assert array.duration == pytest.approx(4.99)
    assert array["XX.A1.."] is array[1]
    assert array.coordinates.shape == (3, 3)
    assert array.canonical().shape == (3, 3, 500)
    assert array.coordinate_covariances is None
    assert array.frame is None

    with pytest.raises(KeyError):
        array["XX.B.."]


def test_array_covariances_with_missing_stations():
    """Stations without covariance get zeros if another has one."""
    stations = list(_array())
    stations[1] = StationData(
        stations[1].sid, stations[1].records, (10.0, 0.0),
        coordinate_covariance=0.01,
    )
    covariances = ArrayData(stations).coordinate_covariances

    assert covariances.shape == (3, 3, 3)
    assert covariances[0].sum() == 0.0
    assert covariances[1, 0, 0] == 0.01


def test_array_checks():
    """Duplicate identifiers and different time axes are rejected."""
    station = StationData("XX.A..", _station_records("A"), (0.0, 0.0))
    shifted = StationData(
        "XX.B..", _station_records("B", start=Date(START) + 0.004),
        (5.0, 0.0),
    )

    with pytest.raises(ValueError, match="unique"):
        ArrayData([station, station])
    with pytest.raises(ValueError, match="no resampling"):
        ArrayData([station, shifted])
    with pytest.raises(ValueError):
        ArrayData([])
    with pytest.raises(TypeError):
        ArrayData([station], frame=(13.0, 46.0))


def test_array_projection_requires_support():
    """Horizontal motion needs horizontal components at every station."""
    stations = list(_array(2)) + [
        StationData("XX.V..", [_record("XX.V..HHZ")], (0.0, 5.0))
    ]
    array = ArrayData(stations)

    assert array.vertical().shape == (3, 500)
    assert array.unsupported_stations((1.0, 0.0, 0.0)) == ("XX.V..",)

    with pytest.raises(ValueError, match="XX.V.."):
        array.horizontal()


# ----------------------------------------------------------------------------
# Builders
# ----------------------------------------------------------------------------


def test_from_collection_geographic_coordinates():
    """
    Positions from the metadata are ENU coordinates in a centred frame.

    The reference is the mean longitude, latitude and elevation of the
    stations; the coordinates equal ``geo_to_enu`` in that frame.
    """
    array = ArrayData.from_collection(_array_collection())
    mean = np.mean(list(LOCATIONS.values()), axis=0)

    np.testing.assert_allclose(array.frame.ref_geo, mean)

    for station in array:
        expected = geo_to_enu(LOCATIONS[station.sid.split(".")[1]],
                              array.frame)
        np.testing.assert_allclose(station.coordinates, expected)

    # 0.002 degrees of longitude at 46.06 N: about 154.6 m.
    assert array["XX.A02."].coordinates[0] - array["XX.A01."].coordinates[
        0
    ] == pytest.approx(154.6, abs=0.5)


def test_from_collection_reference_point():
    """A reference point fixes the origin of the frame."""
    point = WgsPoint(13.25, 46.06, 120.0)
    array = ArrayData.from_collection(_array_collection(), reference=point)

    assert array.frame == MetricFrame((13.25, 46.06, 120.0), "enu")
    np.testing.assert_allclose(array["XX.A01."].coordinates, 0.0, atol=1e-6)


def test_from_collection_explicit_coordinates():
    """Explicit local coordinates replace the metadata."""
    coordinates = {"XX.A01.": (0, 0), "XX.A02.": (10, 0), "XX.A03.": (0, 10)}
    array = ArrayData.from_collection(
        _array_collection(), coordinates=coordinates
    )

    assert array.frame is None
    np.testing.assert_array_equal(array["XX.A03."].coordinates, [0, 10, 0])

    with pytest.raises(ValueError, match="missing"):
        ArrayData.from_collection(
            _array_collection(), coordinates={"XX.A01.": (0, 0)}
        )


def test_from_collection_missing_location():
    """Without metadata or coordinates the positions are unknown."""
    collection = _collection(_station_records("A01"))

    with pytest.raises(ValueError, match="unknown"):
        ArrayData.from_collection(collection)


def test_from_collection_common_window():
    """The default window is the interval covered by all channels."""
    records = (
        _station_records("A01", location=LOCATIONS["A01"])
        + _station_records("A02", location=LOCATIONS["A02"],
                           start=Date(START) + 1.0)
    )
    array = ArrayData.from_collection(_collection(records))

    assert array.starttime == Date(START) + 1.0
    assert array.nsamp == 400
    np.testing.assert_array_equal(
        array["XX.A01."].records[0].data, records[0].data[100:]
    )


def test_from_collection_explicit_window():
    """An explicit window is cut with both bounds included."""
    array = ArrayData.from_collection(
        _array_collection(),
        starttime=Date(START) + 1.0,
        endtime=Date(START) + 2.0,
    )

    assert array.nsamp == 101
    assert array.starttime == Date(START) + 1.0


def test_from_collection_does_not_modify_input():
    """The records of the array are new records."""
    collection = _array_collection()
    array = ArrayData.from_collection(collection, endtime=Date(START) + 1.0)
    array[0].records[0].data[:] = 0.0

    assert collection.records()[0].nsamp == 500
    assert np.any(collection.records()[0].data != 0.0)


def _gap_collection():
    """Collection in which station A02 has a gap in channel HHN."""
    records = _station_records("A01", location=LOCATIONS["A01"])
    records += _station_records("A03", location=LOCATIONS["A03"])
    z, n, e = _station_records("A02", location=LOCATIONS["A02"])
    first = _record("XX.A02..HHN", nsamp=200, location=LOCATIONS["A02"])
    second = _record("XX.A02..HHN", nsamp=200, start=Date(START) + 3.0,
                     location=LOCATIONS["A02"])
    return _collection(records + [z, e, first, second])


def test_from_collection_drops_stations_with_gaps():
    """Stations with gaps are excluded with a warning."""
    with pytest.warns(ArrayDataWarning, match="XX.A02."):
        array = ArrayData.from_collection(_gap_collection())

    assert array.sids == ("XX.A01.", "XX.A03.")
    assert array.excluded[0][0] == "XX.A02."
    assert "XX.A02..HHN" in array.excluded[0][1]


def test_from_collection_raises_on_gaps():
    """With gaps='raise' an incomplete station is an error."""
    with pytest.raises(ValueError, match="XX.A02..HHN"):
        ArrayData.from_collection(_gap_collection(), gaps="raise")

    with pytest.raises(ValueError):
        ArrayData.from_collection(_gap_collection(), gaps="fill")


def test_from_collection_selection():
    """Channels of different instruments must be selected."""
    records = _station_records("A01", location=LOCATIONS["A01"])
    records += _station_records("A01", location=LOCATIONS["A01"],
                                instrument="HN")
    records += _station_records("A02", location=LOCATIONS["A02"])
    collection = _collection(records)

    with pytest.raises(ValueError, match="select"):
        ArrayData.from_collection(collection)

    array = ArrayData.from_collection(collection, select="XX.*..HH?")
    assert array.sids == ("XX.A01.", "XX.A02.")

    vertical = ArrayData.from_collection(collection, select=["*.*..HHZ"])
    assert all(station.ncomp == 1 for station in vertical)

    with pytest.raises(ValueError, match="no channel"):
        ArrayData.from_collection(collection, select="YY.*")


def test_from_collection_orientations_and_covariance():
    """Per-station orientations and covariances are applied."""
    rotation = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    array = ArrayData.from_collection(
        _array_collection(),
        orientations={"XX.A02.": rotation},
        coordinate_covariance={"XX.A01.": 0.04},
    )

    np.testing.assert_array_equal(array["XX.A02."].orientation, rotation)
    assert array["XX.A01."].coordinate_covariance[0, 0] == 0.04
    assert array["XX.A03."].coordinate_covariance is None

    common = ArrayData.from_collection(_array_collection(),
                                       coordinate_covariance=0.01)
    assert all(st.coordinate_covariance[1, 1] == 0.01 for st in common)


def test_from_collection_type_checks():
    """Inputs of the wrong type are rejected."""
    with pytest.raises(TypeError):
        ArrayData.from_collection(_station_records("A"))
    with pytest.raises(TypeError):
        ArrayData.from_collection(_array_collection(), coordinates=[(0, 0)])


def test_from_records_equals_from_collection():
    """Records and collections give the same array."""
    records = []

    for station, location in LOCATIONS.items():
        records += _station_records(station, location=location)

    first = ArrayData.from_records(records)
    second = ArrayData.from_collection(_collection(records))

    assert first.sids == second.sids
    np.testing.assert_array_equal(first.coordinates, second.coordinates)
    np.testing.assert_array_equal(first.canonical(), second.canonical())

    with pytest.raises(ValueError, match="identifier"):
        ArrayData.from_records([Record(data=[1.0], delta=1.0, time=START)])
    with pytest.raises(TypeError):
        ArrayData.from_records([np.zeros(3)])


def test_sesame_excerpt():
    """
    Three stations of the SESAME simulation with local coordinates.

    The SAC files keep East and North (m) in USER7 and USER8; the
    vertical motion of the array is the Z record of each station.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        collection = reader(str(DATA / "*.sac"))

    coordinates = {}

    for stream in collection:
        sac = stream.records()[0].head.extra["sac"]
        coordinates[stream.sid.rsplit(".", 1)[0]] = (
            sac["USER7"], sac["USER8"], sac["USER9"]
        )

    array = ArrayData.from_collection(collection, coordinates=coordinates)

    assert array.sids == (".S1001.", ".S1002.", ".S1003.")
    assert array.nsamp == 6858
    assert array.rate == pytest.approx(1.0 / 0.00875, rel=1e-6)
    np.testing.assert_array_equal(array[".S1002."].coordinates,
                                  [2028.0, 2000.0, 0.0])
    np.testing.assert_array_equal(
        array.vertical()[0],
        collection[".S1001..Z"].records()[0].data,
    )
