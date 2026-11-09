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
Tests of :mod:`shakelab.signals.stationtext`.

The text format of the FDSN station service, at channel and station
level, read from the files of ``tests/signals/data/stationtext``
directly, through :class:`MetadataCollection` and through
:func:`metadata_reader`.
"""

from __future__ import annotations

import os
import shutil

import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.io import metadata_reader
from shakelab.signals.metadata import MetadataCollection
from shakelab.signals.stationtext import is_station_text, read_station_text

DATA = os.path.join(os.path.dirname(__file__), "data", "stationtext")


def test_channels():
    """Coordinates, depth in km, orientation, rate and sensitivity."""
    collection = read_station_text(os.path.join(DATA, "channels.txt"))
    assert collection.sid == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ", "IV.MODE.00.HHZ"
    ]

    (hnz,) = collection["IV.MODE..HNZ"]
    assert hnz.location.longitude == pytest.approx(10.949167)
    assert hnz.location.latitude == pytest.approx(44.629722)
    assert hnz.location.elevation == pytest.approx(41.0)
    assert hnz.orientation.dip == -90.0
    assert hnz.sample_rate == 100.0
    assert hnz.sensor.description == "Kinemetrics Episensor-ES-T"
    assert hnz.starttime == Date("2008-04-07T08:00:00")
    assert hnz.endtime is None
    assert hnz.station is None

    sensitivity = hnz.response.sensitivity
    assert sensitivity.value == 321707.0
    assert sensitivity.frequency == 10.625
    assert sensitivity.input_units == "m/s**2"
    assert hnz.response.stage == []

    (hhz,) = collection["IV.MODE.00.HHZ"]
    assert hhz.depth == pytest.approx(0.0025)
    assert hhz.endtime == Date("2020-01-01T00:00:00")


def test_stations():
    """Station epochs with site names."""
    collection = read_station_text(os.path.join(DATA, "stations.txt"))
    assert len(collection) == 0
    stations = collection.stations()
    assert [station.code for station in stations] == ["IV.MODE", "IV.MOXX"]
    assert stations[0].site == "Modena"
    assert stations[1].endtime == Date("2015-01-01T00:00:00")


def test_content_and_empty_fields():
    """Bytes are accepted; empty numeric fields become None."""
    text = (
        b"#Network|Station|Location|Channel|Latitude|Longitude|"
        b"Elevation|Depth|Azimuth|Dip|SensorDescription|Scale|ScaleFreq|"
        b"ScaleUnits|SampleRate|StartTime|EndTime\n"
        b"XX|TEST||BHZ|||||||||||||\n"
    )
    (epoch,) = read_station_text(text)
    assert epoch.sid == "XX.TEST..BHZ"
    assert epoch.location is None
    assert epoch.orientation is None
    assert epoch.response is None


def test_wrong_lines():
    """Lines with other numbers of fields are rejected."""
    with pytest.raises(ValueError, match="line 2"):
        read_station_text("#header\nIV|MODE|x\n")


# ----------------------------------------------------------------------------
# Recognition and metadata readers
# ----------------------------------------------------------------------------


def test_is_station_text(tmp_path):
    """Header lines of the station and channel levels."""
    assert is_station_text(os.path.join(DATA, "channels.txt"))
    assert is_station_text(os.path.join(DATA, "stations.txt"))
    assert is_station_text(
        b"#Network|Station|Latitude|Longitude|Elevation|SiteName|"
        b"StartTime|EndTime\n"
    )
    assert not is_station_text(b"#Network|Station|Other\n")
    assert not is_station_text(b"Network|Station|Latitude|Longitude|"
                               b"Elevation|SiteName|StartTime|EndTime\n")
    assert not is_station_text(str(tmp_path / "missing.txt"))
    assert not is_station_text(3)


def test_metadata_collection_reads_text(tmp_path):
    """Recognized from the content; explicit format for any name."""
    collection = MetadataCollection(os.path.join(DATA, "channels.txt"))
    assert len(collection) == 4

    renamed = tmp_path / "channels.dat"
    shutil.copy(os.path.join(DATA, "channels.txt"), renamed)
    collection = MetadataCollection(str(renamed))
    assert "IV.MODE..HNZ" in collection

    collection = MetadataCollection()
    with open(os.path.join(DATA, "stations.txt"), "rb") as fobj:
        collection.read(fobj.read())
    assert len(collection.stations()) == 2

    collection = MetadataCollection(str(renamed), format="fdsntext")
    assert len(collection) == 4


def test_metadata_reader_reads_text(tmp_path):
    """Text files found in a directory are read, others skipped."""
    folder = tmp_path / "meta"
    folder.mkdir()
    for name in ("channels.txt", "stations.txt"):
        shutil.copy(os.path.join(DATA, name), folder / name)
    (folder / "notes.txt").write_text("not metadata\n")

    with pytest.warns(UserWarning, match="notes.txt"):
        collection = metadata_reader(str(folder))
    assert len(collection) == 4
    assert [s.code for s in collection.stations()] == ["IV.MODE", "IV.MOXX"]
