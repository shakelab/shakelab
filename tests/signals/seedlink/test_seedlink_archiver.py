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
Tests of :mod:`shakelab.signals.seedlink.archiver`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.seedlink import archiver
from shakelab.signals.seedlink.archiver import Archive

TIME = Date("2026-11-05T08:17:42.5").to_unix_ns()


@pytest.mark.parametrize("rotate, seconds", [
    (600, 600), (2.5, 2.5), ("30s", 30), ("10min", 600), (" 1 h", 3600),
    ("1d", 86400), ("1.5h", 5400),
])
def test_rotate(rotate, seconds):
    assert archiver._rotate_ns(rotate) == round(seconds * 1e9)


@pytest.mark.parametrize("rotate", ["10m", "record", "h", "-1h", 0, -5])
def test_rotate_invalid(rotate):
    with pytest.raises(ValueError):
        archiver._rotate_ns(rotate)


def test_rotate_type():
    with pytest.raises(TypeError):
        archiver._rotate_ns(True)
    with pytest.raises(TypeError):
        Archive(5)


def test_paths():
    sid = "IV.MODE..HNZ"
    assert Archive("data/today.mseed").path(sid, TIME) == Path(
        "data/today.mseed"
    )
    assert Archive("d/{net}.{sta}.{loc}.{cha}.mseed").path(sid, TIME) == (
        Path("d/IV.MODE..HNZ.mseed")
    )
    hourly = Archive("d/{sta}.{time:%Y%m%dT%H%M}.mseed", rotate="1h")
    assert hourly.path(sid, TIME) == Path("d/MODE.20261105T0800.mseed")
    ten = Archive("d/{cha}_{time:%H%M%S}_{year}{jday}", rotate=600)
    assert ten.path(sid, TIME) == Path("d/HNZ_081000_2026309")


def test_sds_path():
    sds = Archive("archive", layout="sds")
    assert sds.path("OX.CAE.00.HHZ", TIME) == Path(
        "archive/2026/OX/CAE/HHZ.D/OX.CAE.00.HHZ.D.2026.309"
    )
    early = Date("2026-01-05T00:00:00").to_unix_ns()
    assert sds.path("OX.CAE..HHZ", early).name == "OX.CAE..HHZ.D.2026.005"


@pytest.mark.parametrize("path, kwargs", [
    ("d/{station}.mseed", {}),
    ("d/{time:%H}.mseed", {}),
    ("d/{net", {}),
    ("archive", {"layout": "sds", "rotate": "1d"}),
    ("archive", {"layout": "SDS"}),
])
def test_invalid_archive(path, kwargs):
    with pytest.raises(ValueError):
        Archive(path, **kwargs)


def test_write_and_file_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(archiver, "MAX_OPEN_FILES", 2)
    archive = Archive(tmp_path / "x" / "{sta}.mseed")
    for index, station in enumerate(["A", "B", "C", "A"]):
        archive.write(bytes([index]) * 4, f"N.{station}..HHZ", TIME)
    assert len(archive._files) == 2
    archive.close()
    assert archive._files == {}
    assert (tmp_path / "x" / "A.mseed").read_bytes() == b"\0" * 4 + b"\3" * 4
    assert (tmp_path / "x" / "C.mseed").read_bytes() == b"\2" * 4
