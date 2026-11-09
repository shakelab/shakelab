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
Tests of the command-line client :mod:`shakelab.signals.fdsnws.cli`,
against the fake server of ``fdsn_server.py``.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from shakelab.signals import io
from shakelab.signals.fdsnws.cli import main

from fdsn_server import SIGNALS_DATA, read_file

WINDOW = ["--start", "2012-05-20T02:04:00", "--end", "2012-05-20T02:05:30"]


@pytest.fixture(autouse=True)
def mseed_backend(monkeypatch):
    """Records read with the default backend of signals.io."""
    try:
        import cymseed3  # noqa: F401
    except ImportError:
        monkeypatch.setattr(io, "USE_LIBMSEED", False)


def test_waveforms_raw(server, tmp_path, capsys):
    """Records written as received, with the StationXML."""
    output = tmp_path / "mode.mseed"
    xml = tmp_path / "mode.xml"
    code = main(["waveforms", server.url, "IV.MODE..HN?", *WINDOW,
                 "-o", str(output), "--metadata", str(xml)])
    assert code == 0
    raw = read_file(SIGNALS_DATA, "MODE.mseed")
    assert output.read_bytes() == raw
    assert xml.read_bytes() == read_file(SIGNALS_DATA, "MODE.xml")
    assert "bytes written" in capsys.readouterr().err

    assert main(["waveforms", server.url, "IV.MODE", *WINDOW,
                 "-o", str(output)]) == 1
    assert "exists" in capsys.readouterr().err


def test_waveforms_corrected(server, tmp_path):
    """Corrected records written with signals.io, one SAC per channel."""
    template = str(tmp_path / "{net}.{sta}.{loc}.{cha}.sac")
    code = main(["waveforms", server.url, "IV.MODE..HN[ZN]", *WINDOW,
                 "--correct", "--pad", "30", "--lowpass", "40",
                 "-o", template])
    assert code == 0
    assert sorted(os.listdir(tmp_path)) == [
        "IV.MODE..HNN.sac", "IV.MODE..HNZ.sac"
    ]
    (record,) = io.reader(str(tmp_path / "IV.MODE..HNZ.sac")).records()
    assert record.nsamp == 9001
    assert np.abs(record.data).max() == pytest.approx(0.3273, rel=1e-2)


def test_waveforms_format_needs_processing(server, tmp_path, capsys):
    """A format other than MiniSEED needs a correction or a filter."""
    code = main(["waveforms", server.url, "IV.MODE", *WINDOW,
                 "-o", str(tmp_path / "x.sac"), "--format", "sac"])
    assert code == 1
    assert "--format" in capsys.readouterr().err


def test_waveforms_without_data(server, tmp_path, capsys):
    """No data: exit status 1 and no file."""
    server.replies["/fdsnws/dataselect/1/query"] = lambda req: (204, b"", {})
    code = main(["waveforms", server.url, "IV.MODE", *WINDOW,
                 "-o", str(tmp_path / "x.mseed")])
    assert code == 1
    assert "no data" in capsys.readouterr().err
    assert os.listdir(tmp_path) == []


def test_metadata_table(server, capsys):
    """One line per channel epoch, or per station."""
    assert main(["metadata", server.url, "IV.MODE..HNZ"]) == 0
    out, err = capsys.readouterr()
    assert out.startswith("IV.MODE..HNZ")
    assert "10.9492" in out and "44.6297" in out and "100 Hz" in out
    assert "1 channel epoch(s)" in err

    assert main(["metadata", server.url, "IV", "--level", "station",
                 "--text"]) == 0
    out, err = capsys.readouterr()
    assert [line.split()[0] for line in out.splitlines()] == [
        "IV.MODE", "IV.MOXX"
    ]


def test_metadata_to_file(server, tmp_path):
    """The reply written to a file."""
    output = tmp_path / "channels.txt"
    assert main(["metadata", server.url, "IV.MODE", "--text",
                 "-o", str(output)]) == 0
    assert output.read_text().startswith("#Network")


def test_events(server, tmp_path, capsys):
    """Table of the events, or QuakeML written to a file."""
    assert main(["events", server.url, "--minmag", "5.5",
                 "--region", "9", "13", "43.5", "46"]) == 0
    out, err = capsys.readouterr()
    lines = out.splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("2012-05-20T02:03:52")
    assert "5.9 Mw" in lines[0]
    query = server.requests[-1].query
    assert query["minmagnitude"] == "5.5"
    assert query["minlongitude"] == "9"

    output = tmp_path / "events.xml"
    assert main(["events", server.url, "-o", str(output)]) == 0
    assert b"quakeml" in output.read_bytes()


def test_availability(server, capsys):
    """One line per span."""
    assert main(["availability", server.url, "IV.MODE..HNZ", *WINDOW,
                 "--merge", "quality"]) == 0
    out, err = capsys.readouterr()
    assert len(out.splitlines()) == 2
    assert server.requests[-1].query["merge"] == "quality"


def test_datacenters(capsys):
    """The built-in table, optionally by service and with addresses."""
    assert main(["datacenters", "--service", "event", "--urls"]) == 0
    out, err = capsys.readouterr()
    assert "USGS" in out
    assert "https://earthquake.usgs.gov/fdsnws/event/1/" in out
    assert "BGR " not in out


def test_errors(capsys):
    """Errors: message on standard error and exit status 1."""
    assert main(["metadata", "NOWHERE"]) == 1
    assert "unknown data center" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(["waveforms", "INGV"])
