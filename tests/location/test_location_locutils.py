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
"""Tests for the semicolon-separated station list reader."""

from __future__ import annotations

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.location.locutils import load_station_list


def test_station_list_reads_coordinates_and_skips_bad_rows(tmp_path):
    """Keep station codes and longitude, latitude and signed elevation."""
    path = tmp_path / "stations.csv"
    path.write_text(
        "code;network;latitude;longitude;elevation;comment\n"
        " AAA ;OX;46,25;13,125;175,5;ridge\n"
        "BBB;IV;44.5;11.25;-5;borehole\n"
        "\n"
        ";OX;46;13;0;missing code\n"
        "CCC;OX;46;13\n"
        "DDD;OX;bad;13;0\n",
        encoding="ascii",
    )
    stations = load_station_list(path)
    assert set(stations) == {"AAA", "BBB"}
    assert stations["AAA"].station_id == "AAA"
    assert stations["AAA"].location == WgsPoint(13.125, 46.25, 175.5)
    assert stations["BBB"].location == WgsPoint(11.25, 44.5, -5)


def test_station_list_with_only_header(tmp_path):
    """Return an empty mapping when the file has no station rows."""
    path = tmp_path / "stations.csv"
    path.write_text("code;network;latitude;longitude;elevation\n",
                    encoding="ascii")
    assert load_station_list(path) == {}
