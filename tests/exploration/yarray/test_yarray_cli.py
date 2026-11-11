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
Unit tests for :mod:`shakelab.exploration.yarray.cli`.

The command is run in-process on the SESAME excerpt: processing from
an array file, options replacing the configuration file, the section
``graphics``, saved results and errors.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from shakelab.exploration.yarray.cli import main
from shakelab.exploration.yarray.config import BeamformingConfig
from shakelab.exploration.yarray.io import read_results

DATA = Path(__file__).parent / "data"
ARRAY = str(DATA / "sesame_array.json")
OPTIONS = [
    "--input", ARRAY, "--fmin", "3", "--fmax", "8", "--nfreq", "3",
    "--grid-min", "100", "--grid-max", "2000", "--grid-size", "30",
    "--azimuth-size", "36", "--window-length", "10",
    "--windows-per-block", "20",
]


def _config_file(tmp_path, **extra):
    """Configuration file with a graphics section."""
    data = {
        "frequencies": {"values": [3.0, 5.0]},
        "grid": {"minimum": 100, "maximum": 2000, "size": 30,
                 "azimuth_size": 36},
        "window": {"length": 10, "windows_per_block": 20},
        "method": {"name": "capon"},
        "output": {"mode": "both"},
        **extra,
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(data))
    return path


def test_process_from_options(tmp_path, capsys):
    """Options alone define the configuration; result written."""
    output = tmp_path / "result.pkl"
    written = tmp_path / "used.json"

    assert main(OPTIONS + ["-o", str(output), "--write-config",
                           str(written)]) == 0

    out = capsys.readouterr().out
    assert "ArrayData: 9 stations" in out
    assert "BeamformingResult: capon, vertical" in out
    assert len(read_results(output)) == 1
    assert BeamformingConfig.read(written).frequencies.size == 3


def test_options_replace_the_configuration(tmp_path):
    """Command-line values take precedence over the file."""
    config = _config_file(tmp_path)
    written = tmp_path / "used.json"

    assert main(["--input", ARRAY, "--config", str(config), "--method",
                 "bartlett", "--frequencies", "4", "--write-config",
                 str(written)]) == 0

    used = BeamformingConfig.read(written)
    assert used.method.value == "bartlett"
    assert list(used.frequencies.values) == [4.0]
    assert used.window.length == 10.0


def test_append_and_plot_saved_result(tmp_path):
    """Results accumulate; a saved result is plotted without data."""
    pytest.importorskip("matplotlib").use("Agg")
    config = _config_file(tmp_path)
    output = tmp_path / "result.pkl"

    for _ in range(2):
        assert main(["--input", ARRAY, "--config", str(config), "-o",
                     str(output), "--append"]) == 0

    assert len(read_results(output)) == 2
    assert main(["--result", str(output), "--plot-spectrum", "--frequency",
                 "5", "--spectrum-view", "azimuth", "--save-figures",
                 str(tmp_path / "figures")]) == 0
    assert (tmp_path / "figures" / "spectrum_f5_b0_azimuth_power.png"
            ).is_file()
    assert main(["--result", str(output), "--plot-density",
                 "--density-coordinate", "log_velocity", "--save-figures",
                 str(tmp_path / "figures")]) == 0
    assert (tmp_path / "figures" / "density_vertical.png").is_file()


def test_graphics_section(tmp_path):
    """Figures requested by the configuration file."""
    pytest.importorskip("matplotlib").use("Agg")
    figures = tmp_path / "figures"
    config = _config_file(tmp_path, graphics={
        "dispersion": True, "ranks": [1],
        "save": {"directory": str(figures), "format": "png", "dpi": 50},
    })

    assert main(["--input", ARRAY, "--config", str(config)]) == 0
    assert (figures / "dispersion_all.png").is_file()


@pytest.mark.parametrize(
    "arguments, message",
    [
        (["--input", ARRAY, "--fmin", "3"], "go together"),
        (["--input", ARRAY], "frequencies and grid"),
        (["--input", "missing.json", "--frequencies", "3", "--grid-min",
          "100", "--grid-max", "200", "--grid-size", "5"], "missing"),
        (OPTIONS + ["--method", "bartlett", "--components",
                    "rayleigh_joint"], "Capon"),
        (["--data", str(DATA / "sesame" / "*.sac"), "--frequencies", "3",
          "--grid-min", "100", "--grid-max", "200", "--grid-size", "5"],
         "unknown"),
        (["--result", "r.pkl", "--output", "x.pkl"], "needs data"),
        (OPTIONS + ["--plot-spectrum"], "--frequency"),
    ],
)
def test_errors(arguments, message, capsys):
    """Errors are reported on one line with exit status 1."""
    assert main(arguments) == 1
    assert message in capsys.readouterr().err


def test_existing_output(tmp_path, capsys):
    """Existing result files are protected unless --owrite."""
    output = tmp_path / "result.pkl"
    output.write_bytes(b"")

    assert main(OPTIONS + ["-o", str(output)]) == 1
    assert "exists" in capsys.readouterr().err
    assert main(OPTIONS + ["-o", str(output), "--owrite"]) == 0


def test_response_command(tmp_path, capsys):
    """Response of a coordinate file, velocity limits and figures."""
    pytest.importorskip("matplotlib").use("Agg")
    stations = tmp_path / "stations.txt"
    stations.write_text(
        "# name x y\nA 0 0\nB 30 0\nC 0 40\nD 25 35 1\nE,-20,10\n"
        "F 10 -30\nG -35 -20\n"
    )
    figures = tmp_path / "figures"

    assert main(["response", "--coordinates", str(stations), "--fmin", "1",
                 "--fmax", "10", "--nfreq", "3", "--search-max", "0.03",
                 "--plot-all", "--save-figures", str(figures)]) == 0
    output = capsys.readouterr().out
    assert "usable wavenumbers" in output and "minimum velocity" in output
    assert sorted(path.stem for path in figures.iterdir()) == [
        "array_geometry", "array_response", "resolution", "search_domain",
        "velocity_limits",
    ]


def test_response_command_from_array(capsys):
    """Response of an array file; errors of the options."""
    assert main(["response", "--input", ARRAY, "--grid-max", "0.3",
                 "--grid-size", "61"]) == 0
    assert "9 stations" in capsys.readouterr().out
    assert main(["response", "--input", ARRAY, "--plot-search-domain"]) == 1
    assert "--search-max" in capsys.readouterr().err
    assert main(["response", "--input", ARRAY, "--grid-max", "0.3"]) == 1
    assert "go together" in capsys.readouterr().err
