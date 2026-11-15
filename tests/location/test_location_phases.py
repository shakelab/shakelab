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
"""Regression checks for the phase files supplied with the examples."""

from __future__ import annotations

from pathlib import Path

import pytest

from shakelab.libutils.timeN import Date
from shakelab.location.phases import (
    PhasePick,
    load_picks,
    normalize_phases,
    read_hpl_picks,
    read_pun_origin,
)


DATA = Path(__file__).resolve().parents[2] / "examples/location/data"
PICKING = DATA / "Picking"
EVENT = PICKING / "20240327211937"


@pytest.mark.parametrize(
    ("extension", "count"), [(".obs", 168), (".dat", 167), (".xml", 168)],
)
def test_load_socchieve_picks(extension, count):
    """Dispatch by extension and decode the recorded arrival times."""
    picks = load_picks(EVENT.with_suffix(extension))
    assert len(picks) == count
    assert all(isinstance(pick.time, Date) for pick in picks)
    arrivals = {(pick.station, pick.phase[0]): pick.time for pick in picks}
    assert arrivals["INVI", "P"] == Date("2024-03-27T21:19:39.70Z")
    assert arrivals["PURA", "S"] == Date("2024-03-27T21:19:41.24Z")


def test_hpl_origin_and_primary_picks():
    """Read the catalog origin and primary HPL records with blank fields."""
    origin, picks = read_hpl_picks(EVENT.with_suffix(".hpl"))
    pun = read_pun_origin(EVENT.with_suffix(".pun"))
    assert origin.time == pun.time == Date("2024-03-27T21:19:37.37Z")
    assert len(picks) == 93
    assert load_picks(EVENT.with_suffix(".hpl")) == picks
    assert pun.depth_km == pytest.approx(9.39)
    assert pun.longitude == pytest.approx(12.8076666667)
    assert pun.latitude == pytest.approx(46.3588333333)

    by_station = {pick.station: pick for pick in picks}
    assert by_station["MOV"].weight == 2
    assert by_station["MOV"].polarity is None
    assert by_station["POLC"].phase == "P"
    assert by_station["POLC"].onset is None


def test_quakeml_pick_metadata_and_time_offset(tmp_path):
    """Retain waveform identifiers and normalize time offsets to UTC."""
    path = tmp_path / "picks.XML"
    path.write_text(
        """<eventParameters xmlns="urn:example:quakeml">
        <pick>
            <time><value>2024-03-28T00:19:39.700+03:00</value></time>
            <waveformID networkCode="OX" stationCode="ABC"
                        locationCode="00" channelCode="HHZ"/>
            <phaseHint>Pg</phaseHint>
            <evaluationMode>manual</evaluationMode>
            <evaluationStatus>reviewed</evaluationStatus>
        </pick>
        <pick><waveformID stationCode="UNTIMED"/></pick>
        <pick><waveformID networkCode="OX"/></pick>
        </eventParameters>""",
        encoding="ascii",
    )
    pick, untimed = load_picks(path)
    assert (pick.network, pick.station, pick.location, pick.channel) == (
        "OX", "ABC", "00", "HHZ",
    )
    assert pick.phase == "Pg"
    assert pick.time == Date("2024-03-27T21:19:39.70Z")
    assert pick.picker == "manual"
    assert pick.extra["evaluation_status"] == "reviewed"
    assert untimed.station == "UNTIMED"
    assert untimed.time is None


def test_normalize_phase_families_in_place():
    """Collapse P and S families while retaining other phase labels."""
    picks = [PhasePick("STA", phase, None)
             for phase in (" Pg ", "sn", "Lg", "")]
    assert normalize_phases(picks) is picks
    assert [pick.phase for pick in picks] == ["P", "S", "Lg", ""]


def test_pick_reader_rejects_unknown_extension(tmp_path):
    """Report unsupported formats before opening the file."""
    with pytest.raises(ValueError, match="Unsupported pick file extension"):
        load_picks(tmp_path / "picks.txt")
