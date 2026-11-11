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
Array processing from the command line (``slab-yarray``).

Every command of this example is run through the entry point of the
command (:func:`shakelab.exploration.yarray.cli.main`), with the same
arguments as in a shell, on the nine SESAME stations; files and figures
go to a temporary folder. In a shell, after ``pip install -e .``, the
same commands are::

    slab-yarray --input sesame_array_9.json --config capon.json \\
        -o capon.pkl --plot-dispersion --save-figures figures
    slab-yarray --input sesame_array_9.json --config capon.json \\
        --method music --fmin 4 --fmax 8 --nfreq 5 -o music.pkl
    slab-yarray --input sesame_array_9.json --config capon.json \\
        --components rayleigh_joint --windows-per-block 40 \\
        --write-config rtbf.json -o rtbf.pkl
    slab-yarray --input day_2.json --config capon.json \\
        -o daily.pkl --append
    slab-yarray --result capon.pkl --plot-spectrum --frequency 6 \\
        --spectrum-view azimuth --save-figures figures
    slab-yarray --result capon.pkl --plot-density \\
        --density-coordinate log_velocity --save-figures figures
    slab-yarray response --input sesame_array_9.json --fmin 2 --fmax 20 \\
        --search-max 0.1 --plot-all --save-figures figures
    slab-yarray response --coordinates layout.txt --plot-response

(``python -m shakelab.exploration.yarray`` is equivalent). Waveforms
with metadata are read directly with ``--data "data/*.mseed"
--metadata array.xml --select "XX.*..HH?"``; the SESAME files have no
geographic positions, so they are read through the array file.

Run from the root of the repository::

    python examples/exploration/yarray_command_line.py
"""

from __future__ import annotations

import json
import shlex
import tempfile
import warnings
from pathlib import Path

from shakelab.exploration.yarray.cli import main
from shakelab.exploration.yarray.io import read_results

DATA = Path("tests/exploration/yarray/data").resolve()
ARRAY = DATA / "sesame_array_9.json"

warnings.filterwarnings("ignore", message="STLO")


def run(command: str) -> None:
    """Run one command line and stop on errors."""
    print(f"\n$ slab-yarray {command}")
    status = main(shlex.split(command))

    if status != 0:
        raise RuntimeError(f"exit status {status}")


with tempfile.TemporaryDirectory() as folder:
    work = Path(folder)
    figures = work / "figures"

    # 1. Configuration file: processing sections and figures (the
    #    section "graphics" is read by the command line only).
    config = {
        "frequencies": {"minimum": 3.0, "maximum": 10.0, "size": 10},
        "grid": {"minimum": 100.0, "maximum": 1500.0, "size": 60,
                 "azimuth_size": 72},
        "window": {"type": "fixed", "length": 2.0,
                   "windows_per_block": 30},
        "method": {"name": "capon"},
        "output": {"mode": "both", "maxima": 2},
        "graphics": {"dispersion": True, "ranks": [1]},
    }
    (work / "capon.json").write_text(json.dumps(config, indent=2))

    # 2. Processing with the configuration; the figures of the section
    #    "graphics" are saved.
    run(f"--input {ARRAY} --config {work / 'capon.json'} "
        f"-o {work / 'capon.pkl'} --save-figures {figures}")

    # 3. Options of the command line replace the values of the file.
    run(f"--input {ARRAY} --config {work / 'capon.json'} "
        f"--method music --fmin 4 --fmax 8 --nfreq 5 "
        f"-o {work / 'music.pkl'}")

    # 4. Rayleigh three-component beamforming; the configuration used
    #    (file and options) is written for later runs.
    run(f"--input {ARRAY} --config {work / 'capon.json'} "
        f"--components rayleigh_joint --windows-per-block 40 "
        f"--ellipticity-angle-size 35 --maxima 1 "
        f"--write-config {work / 'rtbf.json'} -o {work / 'rtbf.pkl'} "
        f"--plot-ellipticity --save-figures {figures}")
    print(json.loads((work / "rtbf.json").read_text())["components"])

    # 5. Successive analyses appended to one file, as for data
    #    processed day by day (here twice the same record, with the
    #    vertical and the radial components).
    for components in ("vertical", "radial"):
        run(f"--input {ARRAY} --config {work / 'capon.json'} "
            f"--components {components} -o {work / 'daily.pkl'} --append")

    print(f"{len(read_results(work / 'daily.pkl'))} results in the file")

    # 6. A saved result is only plotted: one spectrum and the density
    #    of the maxima.
    run(f"--result {work / 'capon.pkl'} --plot-spectrum --frequency 6 "
        f"--spectrum-view azimuth --save-figures {figures}")
    run(f"--result {work / 'capon.pkl'} --plot-density "
        f"--density-coordinate log_velocity --save-figures {figures}")

    # 7. Geometric response of the array: limits, velocities at the
    #    frequencies of interest, a search domain and all the figures.
    run(f"response --input {ARRAY} --fmin 2 --fmax 20 --nfreq 5 "
        f"--search-max 0.1 --plot-all --save-figures {figures}")

    # 8. The same for a layout still to be deployed: a text file of
    #    names and local coordinates.
    (work / "layout.txt").write_text(
        "# name  east  north (m)\n"
        "C 0 0\nA1 15 0\nA2 -7.5 13\nA3 -7.5 -13\n"
        "B1 0 40\nB2 34.6 -20\nB3 -34.6 -20\n"
    )
    run(f"response --coordinates {work / 'layout.txt'} "
        f"--resolution-level 0.7 --plot-response "
        f"--save-figures {figures / 'layout'}")

    print("\nfigures: " + ", ".join(sorted(
        str(path.relative_to(figures)) for path in figures.rglob("*.png")
    )))
