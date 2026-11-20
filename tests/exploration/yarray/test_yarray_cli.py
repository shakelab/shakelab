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

The command is run in-process: the beamforming on the SESAME excerpt
(processing from an array file, options replacing the configuration
file, the section ``graphics``, saved results and errors), the FDD and
the enhanced FDD on a shear frame written to SAC files, the response
of coordinate and array files.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.signal import cont2discrete, lfilter

from shakelab.exploration.yarray.cli import main
from shakelab.exploration.yarray.config import BeamformingConfig, FDDConfig
from shakelab.exploration.yarray.fdd import FDDResult
from shakelab.exploration.yarray.io import read_result
from shakelab.signals.base import Record
from shakelab.signals.io import writer

DATA = Path(__file__).parent / "data"
ARRAY = str(DATA / "sesame_array.json")
OPTIONS = [
    "beamforming", "--input", ARRAY, "--fmin", "3", "--fmax", "8",
    "--nfreq", "3",
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
    output = tmp_path / "result.jsonl"
    written = tmp_path / "used.json"

    assert main(OPTIONS + ["-o", str(output), "--write-config",
                           str(written)]) == 0

    out = capsys.readouterr().out
    assert "ArrayData: 9 stations" in out
    assert "BeamformingResult: capon, vertical" in out
    assert len(read_result(output)[1]) == 1
    assert BeamformingConfig(written).frequencies.size == 3


def test_options_replace_the_configuration(tmp_path):
    """Command-line values take precedence over the file."""
    config = _config_file(tmp_path)
    written = tmp_path / "used.json"

    assert main(["beamforming", "--input", ARRAY, "--config", str(config),
                 "--method", "bartlett", "--frequencies", "4",
                 "--write-config", str(written)]) == 0

    used = BeamformingConfig(written)
    assert used.method.value == "bartlett"
    assert list(used.frequencies.values) == [4.0]
    assert used.window.length == 10.0


def test_append_and_plot_saved_result(tmp_path):
    """Results accumulate; a saved result is plotted without data."""
    pytest.importorskip("matplotlib").use("Agg")
    config = _config_file(tmp_path)
    output = tmp_path / "result.jsonl"

    for _ in range(2):
        assert main(["beamforming", "--input", ARRAY, "--config",
                     str(config), "-o", str(output), "--append"]) == 0

    assert len(read_result(output)[1]) == 2
    assert main(["beamforming", "--result", str(output), "--plot-spectrum",
                 "--frequency", "5", "--spectrum-view", "azimuth",
                 "--save-figures", str(tmp_path / "figures")]) == 0
    assert (tmp_path / "figures" / "spectrum_f5_b0_azimuth_power.png"
            ).is_file()
    assert main(["beamforming", "--result", str(output), "--plot-density",
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

    assert main(["beamforming", "--input", ARRAY, "--config",
                 str(config)]) == 0
    assert (figures / "dispersion_all.png").is_file()


@pytest.mark.parametrize(
    "arguments, message",
    [
        (["beamforming", "--input", ARRAY, "--fmin", "3"],
         "minimum and maximum"),
        (["beamforming", "--input", ARRAY], "frequencies and grid"),
        (["beamforming", "--input", "missing.json", "--frequencies", "3",
          "--grid-min", "100", "--grid-max", "200", "--grid-size", "5"],
         "missing"),
        (OPTIONS + ["--method", "bartlett", "--components",
                    "rayleigh_joint"], "Capon"),
        (["beamforming", "--data", str(DATA / "sesame" / "*.sac"),
          "--frequencies", "3", "--grid-min", "100", "--grid-max", "200",
          "--grid-size", "5"], "unknown"),
        (["beamforming", "--result", "r.jsonl", "--output", "x.jsonl"],
         "needs data"),
        (OPTIONS + ["--plot-spectrum"], "--frequency"),
    ],
)
def test_errors(arguments, message, capsys):
    """Errors are reported on one line with exit status 1."""
    assert main(arguments) == 1
    assert message in capsys.readouterr().err


def test_csv_and_light_output(tmp_path):
    """--csv writes the maxima, also of a saved result; --no-arrays."""
    config = _config_file(tmp_path)
    output = tmp_path / "result.jsonl"
    table = tmp_path / "maxima.csv"

    assert main(["beamforming", "--input", ARRAY, "--config", str(config),
                 "-o", str(output), "--no-arrays", "--csv", str(table)]) == 0
    assert read_result(output)[1][0]["arrays"] is None
    assert table.read_text().startswith("block,time,frequency")

    again = tmp_path / "again.csv"
    assert main(["beamforming", "--result", str(output), "--csv",
                 str(again)]) == 0
    assert again.read_text() == table.read_text()


def test_existing_output(tmp_path, capsys):
    """Existing result files are protected unless --owrite."""
    output = tmp_path / "result.jsonl"
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
    assert main(["response", "--input", ARRAY, "--grid-max", "0.3"]) == 0
    assert "extent +/-0.3" in capsys.readouterr().out
    assert main(["response", "--input", ARRAY, "--grid-size", "60"]) == 1
    assert "odd" in capsys.readouterr().err


def test_missing_command(capsys):
    """The subcommand is required."""
    with pytest.raises(SystemExit) as error:
        main([])

    assert error.value.code == 2
    assert "command" in capsys.readouterr().err


# ----------------------------------------------------------------------------
# Frequency domain decomposition
# ----------------------------------------------------------------------------

DELTA = 0.02
FDD_OPTIONS = [
    "--fmin", "0.5", "--fmax", "10", "--window-type", "fixed",
    "--window-length", "20", "--windows-per-block", "10",
    "--components", "east",
]


def _frame_files(folder, duration=1200.0, seed=0):
    """
    SAC files of a four-storey shear frame under white noise.

    As in the tests of the FDD: uniform storeys, natural frequencies
    from 1.5 to 8.1 Hz, 2 % damping, exact modal superposition; one
    east component per floor, without positions. Return the natural
    frequencies.
    """
    stiffness = 2.0 * np.eye(4) - np.eye(4, k=1) - np.eye(4, k=-1)
    stiffness[-1, -1] = 1.0
    omega2, shapes = eigh(stiffness)
    omega2 *= (2.0 * np.pi * 1.5) ** 2 / omega2[0]
    rng = np.random.default_rng(seed)
    forces = shapes.T @ rng.standard_normal((4, int(duration / DELTA)))
    motion = np.zeros(forces.shape)

    for k, omega in enumerate(np.sqrt(omega2)):
        num, den, _ = cont2discrete(
            ([1.0], [1.0, 0.04 * omega, omega ** 2]), DELTA, "zoh"
        )
        motion += np.outer(shapes[:, k], lfilter(num[0], den, forces[k]))

    for floor, series in enumerate(motion):
        record = Record(data=series, delta=DELTA,
                        time="2026-01-01T00:00:00",
                        sid=f"XX.F{floor:02d}..HNE")
        writer(record, folder / f"F{floor:02d}.sac")

    return np.sqrt(omega2) / (2.0 * np.pi)


@pytest.fixture(scope="module")
def frame(tmp_path_factory):
    """Pattern of the SAC files of the frame and its frequencies."""
    folder = tmp_path_factory.mktemp("frame")
    frequency = _frame_files(folder)
    return str(folder / "*.sac"), frequency


def _rows(path):
    """Rows of a CSV file."""
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_fdd_command(tmp_path, frame, capsys):
    """Processing, result and configuration files, peaks and modes."""
    pattern, frequency = frame
    output = tmp_path / "fdd.jsonl"
    written = tmp_path / "used.json"
    modes = tmp_path / "modes.csv"

    assert main(["fdd", "--data", pattern, *FDD_OPTIONS, "-o", str(output),
                 "--write-config", str(written), "--modes-csv",
                 str(modes)]) == 0

    out = capsys.readouterr().out
    assert "FDDResult(east, 4 channels" in out
    assert out.count("Mode at") == 4
    rows = _rows(modes)
    assert np.allclose([float(row["frequency"]) for row in rows],
                       frequency, rtol=0.03)
    assert {row["method"] for row in rows} == {"fdd"}
    assert all(row["damping"] == "" for row in rows)
    assert all(float(row["mac"]) > 0.99 for row in rows)
    assert "XX.F03. east std" in rows[0]
    assert FDDConfig(written).components.value == "east"
    # 0.5 to 10 Hz every 0.05 Hz (20 s windows), both bounds included
    assert len(FDDResult(output).singular_values()[0]) == 191


def test_fdd_from_saved_result(tmp_path, frame, capsys):
    """Peaks in a band, given modes and tables from a saved result."""
    pattern, frequency = frame
    output = tmp_path / "fdd.jsonl"
    table = tmp_path / "vectors.csv"
    modes = tmp_path / "modes.csv"

    assert main(["fdd", "--data", pattern, *FDD_OPTIONS, "-o",
                 str(output), "--no-arrays"]) == 0
    capsys.readouterr()

    assert main(["fdd", "--result", str(output), "--peak-range", "4", "7",
                 "--prominence", "3", "--peak-count", "2", "--csv",
                 str(table), "--modes-csv", str(modes)]) == 0
    assert capsys.readouterr().out.count("Mode at") == 2
    assert np.allclose([float(row["frequency"]) for row in _rows(modes)],
                       frequency[1:3], rtol=0.03)
    assert table.read_text().startswith("frequency,block,time,value")

    assert main(["fdd", "--result", str(output), "--modes", "4.3", "1.5",
                 "--modes-csv", str(modes), "--owrite"]) == 0
    assert np.allclose([float(row["frequency"]) for row in _rows(modes)],
                       [1.5, 4.3])


def test_efdd_command(tmp_path, frame, capsys):
    """Enhanced FDD of the peaks: damping in the table, figures."""
    pytest.importorskip("matplotlib").use("Agg")
    pattern, frequency = frame
    modes = tmp_path / "modes.csv"
    figures = tmp_path / "figures"

    assert main(["fdd", "--data", pattern, *FDD_OPTIONS, "--efdd",
                 "--mac", "0.8", "--cycles", "100", "--levels", "0.1", "0.9",
                 "--modes-csv", str(modes), "--plot-singular-values",
                 "--plot-mode-shape", "--along", "index", "--plot-mac",
                 "--plot-efdd", "--save-figures", str(figures)]) == 0

    assert "EFDD: peak" in capsys.readouterr().out
    rows = _rows(modes)
    assert {row["method"] for row in rows} == {"efdd"}
    assert np.allclose([float(row["frequency"]) for row in rows],
                       frequency, rtol=0.01)
    assert all(0.015 < float(row["damping"]) < 0.03 for row in rows)

    names = sorted(path.stem for path in figures.iterdir())
    assert names[:4] == [f"efdd_f{float(row['frequency']):.4g}"
                         for row in rows]
    assert names[4] == "mac" and names[-1] == "singular_values"
    assert len(names) == 10


def test_fdd_graphics_section(tmp_path, frame):
    """Figures requested by the configuration file of the FDD."""
    pytest.importorskip("matplotlib").use("Agg")
    pattern, _ = frame
    figures = tmp_path / "figures"
    config = tmp_path / "fdd.json"
    config.write_text(json.dumps({
        "frequencies": {"minimum": 0.5, "maximum": 10.0},
        "window": {"type": "fixed", "length": 20.0,
                   "windows_per_block": 10},
        "components": "east",
        "graphics": {"singular_values": {"enabled": True, "count": 2},
                     "save": {"directory": str(figures), "dpi": 50}},
    }))

    assert main(["fdd", "--data", pattern, "--config", str(config)]) == 0
    assert [path.name for path in figures.iterdir()] == [
        "singular_values.png"
    ]


@pytest.mark.parametrize(
    "arguments, message",
    [
        (["--fmin", "0.5"], "minimum and maximum"),
        ([], "needs frequencies"),
        (FDD_OPTIONS + ["--cycles", "100", "--duration", "60", "--efdd"],
         "not both"),
        (FDD_OPTIONS + ["--plot-efdd"], "needs --efdd"),
        (FDD_OPTIONS + ["--modes-csv", "m.csv", "--peak-range", "9.8",
                        "10"], "no modes"),
    ],
)
def test_fdd_errors(arguments, message, frame, capsys):
    """Errors of the FDD on one line with exit status 1."""
    assert main(["fdd", "--data", frame[0], *arguments]) == 1
    assert message in capsys.readouterr().err


def test_fdd_result_errors(capsys):
    """Options that need data are refused with a saved result."""
    for extra in (["--efdd"], ["--output", "x.jsonl"]):
        assert main(["fdd", "--result", "r.jsonl", *extra]) == 1
        assert "needs data" in capsys.readouterr().err
