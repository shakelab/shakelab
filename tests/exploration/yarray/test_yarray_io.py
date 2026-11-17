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
Unit tests for :mod:`shakelab.exploration.yarray.io`.

Array files are read from the SESAME excerpt (also with explicit
orientations, reader options and errors); result files are written,
appended and read back, and other files are rejected.
"""

from __future__ import annotations

import json
import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
)
from shakelab.exploration.yarray.config import BeamformingConfig
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.io import (
    read_array_file,
    read_result,
    read_result_header,
    write_result,
    write_table,
)
from shakelab.libutils.timeN import Date

DATA = Path(__file__).parent / "data"
ARRAY = DATA / "sesame_array.json"


@pytest.fixture(autouse=True)
def _ignore_sac_warnings():
    """The SESAME files have local coordinates in STLA/STLO."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="STLO")
        yield


def _edited(tmp_path, change):
    """Copy of the array file with the data, modified by ``change``."""
    shutil.copytree(DATA / "sesame", tmp_path / "sesame")
    data = json.loads(ARRAY.read_text())
    change(data)
    path = tmp_path / "array.json"
    path.write_text(json.dumps(data))
    return path


def read_array(path):
    """Array of the stations of an array file."""
    return ArrayData(read_array_file(path)[0])


def test_read_array():
    """Stations, coordinates, components Z, N, E and data."""
    array = read_array(ARRAY)

    assert array.sids == ("S1001", "S1002", "S1003", "S1006", "S1017",
                          "S1019", "S1021", "S1023", "S1028")
    assert array.frame is None
    assert array.nsamp == 6858
    np.testing.assert_array_equal(array["S1002"].coordinates,
                                  [2028.0, 2000.0, 0.0])
    np.testing.assert_array_equal(
        array["S1001"].orientation, [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    )


def test_explicit_orientation_and_reader(tmp_path):
    """
    Extended components: orientation, reader options, covariance.

    The north file declared as a component at azimuth 45 degrees gives
    an orientation (sin 45, cos 45, 0), placed after Z and E.
    """
    def change(data):
        station = data["stations"][0]
        station["components"]["H"] = {
            "file": station["components"].pop("N"),
            "orientation": [0.70710678, 0.70710678, 0.0],
            "reader": {"format": "sac"},
        }
        station["coordinate_covariance"] = [[0.01, 0.0], [0.0, 0.01]]
        station["coordinates"] = [1.0, 2.0]
        data["reader"] = {"format": "sac", "options": {}}

    array = read_array(_edited(tmp_path, change))
    station = array["S1001"]

    np.testing.assert_allclose(station.orientation[2],
                               [2 ** -0.5, 2 ** -0.5, 0.0])
    assert station.coordinate_covariance[0, 0] == 0.01
    np.testing.assert_array_equal(station.coordinates, [1.0, 2.0, 0.0])


@pytest.mark.parametrize(
    "change, error, match",
    [
        (lambda d: d["stations"][0]["components"].update(H="x.sac"),
         ValueError, "1 to 3"),
        (lambda d: d["stations"][0]["components"].update(
            Z={"file": d["stations"][0]["components"]["Z"]},
            Q=d["stations"][0]["components"].pop("N")),
         ValueError, "orientation"),
        (lambda d: d["stations"][0]["components"].update(Z="missing.sac"),
         FileNotFoundError, "missing"),
        (lambda d: d.update(version=2), ValueError, "version"),
        (lambda d: d["coordinate_system"].update(unit="km"),
         ValueError, "metres"),
        (lambda d: d.update(stations=[]), ValueError, "non-empty"),
        (lambda d: d["stations"][0].update(coordinates="here"),
         TypeError, "coordinates"),
    ],
)
def test_invalid_array_files(tmp_path, change, error, match):
    """Invalid array files give clear errors."""
    with pytest.raises(error, match=match):
        read_array(_edited(tmp_path, change))


def test_array_file_without_coordinates(tmp_path):
    """A station of the array file may have no position."""
    def change(data):
        del data["stations"][0]["coordinates"]

    array = read_array(_edited(tmp_path, change))

    assert not array["S1001"].has_position
    assert array["S1002"].has_position
    assert not array.has_positions


def test_reference_of_the_array_file(tmp_path):
    """The reference of the file is the origin of the local frame."""
    path = _edited(tmp_path, lambda d: d.update(reference=[13.25, 46.06]))
    stations, reference = read_array_file(path)
    array = ArrayData(path)

    assert reference == (13.25, 46.06)
    assert array.frame.ref_geo == pytest.approx((13.25, 46.06, 0.0))
    assert read_array_file(ARRAY)[1] is None

    other = ArrayData(reference=(13.0, 46.0))

    with pytest.raises(ValueError, match="differs"):
        other.read(path)

    with pytest.raises(ValueError, match="longitude"):
        read_array_file(_edited(tmp_path / "bad",
                                lambda d: d.update(reference="here")))


def test_invalid_json(tmp_path):
    """Syntax errors report their position."""
    path = tmp_path / "array.json"
    path.write_text("{ ")

    with pytest.raises(ValueError, match="line 1"):
        read_array(path)


def _result():
    """A small result on the SESAME excerpt."""
    config = BeamformingConfig(
        frequencies=[5.0], grid_min=100.0, grid_max=2000.0, grid_size=20,
        azimuth_size=36, method="bartlett", window_length=10,
        windows_per_block=20,
    )
    return Beamformer(read_array(ARRAY), config).process()


def test_result_lines(tmp_path):
    """
    Header and lines are written, appended and read back.

    The header is checked when appending; arrays go to a .npz file per
    line, with a number added when the name exists.
    """
    path = tmp_path / "out" / "result.jsonl"
    header = {"kind": "test", "config": {"a": 1}, "units": None}
    times = {1.0: (Date("2026-11-16T10:15:00.5"),)}
    write_result(header, {"x": [1.0, np.nan]}, path,
                 arrays={"v": np.arange(3)}, times=times)

    with pytest.raises(FileExistsError,
                       match=r"\(use owrite=True or append=True\)"):
        write_result(header, {"x": []}, path)

    write_result(header, {"x": [np.inf]}, path, times=times,
                 arrays={"v": np.ones(2)}, append=True)
    read, entries = read_result(path)

    assert read == header
    assert read_result_header(path) == header
    assert [entry["x"] for entry in entries] == [[1.0, None], [None]]
    np.testing.assert_array_equal(entries[0]["arrays"]["v"], np.arange(3))
    assert entries[1]["arrays"]["v"].tolist() == [1.0, 1.0]
    assert sorted(item.name for item in path.with_name(
        "result_arrays").iterdir()) == ["20261116T101500.npz",
                                        "20261116T101500_1.npz"]
    assert read_result(path, arrays=False)[1][0]["arrays"] is None

    with pytest.raises(ValueError, match="holds test results"):
        write_result({**header, "kind": "other"}, {}, path, append=True)
    with pytest.raises(ValueError, match="configuration"):
        write_result({**header, "units": "m"}, {}, path, append=True)

    write_result(header, {"x": []}, path, owrite=True)
    assert len(read_result(path)[1]) == 1

    joined = BeamformingResult()
    joined.append([_result(), _result()])
    joined.write(tmp_path / "bf.jsonl")
    assert BeamformingResult(tmp_path / "bf.jsonl").nblocks == joined.nblocks


def test_read_results_rejects_other_files(tmp_path):
    """Other files, other versions, bad lines or no analysis."""
    path = tmp_path / "other.jsonl"

    with pytest.raises(ValueError, match="jsonl"):
        write_result({}, {}, tmp_path / "result.JSON")

    path.write_text('{"a": 1}\n')

    with pytest.raises(ValueError, match="not a YARRAY result"):
        read_result(path)

    path.write_text('{"format": "shakelab.yarray.result", "version": 9}\n')

    with pytest.raises(ValueError, match="version 9"):
        read_result(path)

    path.write_text('{"format": "shakelab.yarray.result", "version": 2}\n')

    with pytest.raises(ValueError, match="no analysis"):
        read_result(path)

    with path.open("a") as stream:
        stream.write("[1, 2\n")

    with pytest.raises(ValueError, match="line 2"):
        read_result(path)


def test_tables(tmp_path):
    """CSV tables: missing values empty, times in ISO 8601."""
    path = tmp_path / "table.csv"
    write_table(path, {"t": [Date("2026-01-01")], "v": [np.nan],
                       "n": [3]})

    assert path.read_text().splitlines() == [
        "t,v,n", "2026-01-01T00:00:00.000000000Z,,3"
    ]

    with pytest.raises(FileExistsError):
        write_table(path, {"a": []})
