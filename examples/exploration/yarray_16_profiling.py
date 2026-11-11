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
Processing times.

The cost of the main choices of a configuration on the SESAME test
data. Times depend on the machine: compare them with each other.

Run (from any folder)::

    python examples/exploration/yarray_16_profiling.py
"""

import time
import warnings
from pathlib import Path

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

array = ArrayData(DATA / "sesame_array.json")
base = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=50,
    grid_min=100.0, grid_max=2000.0, grid_size=100, azimuth_size=180,
    window_type="fixed", window_length=2.0, windows_per_block=40,
)
cases = {
    "bartlett": {"method": "bartlett"},
    "capon": {},
    "music": {"method": "music"},
    "radial": {"components": "radial"},
    "three components": {"components": "three_component_separate"},
    "rfft nearest": {"transform": "rfft"},
    "rfft interp": {"transform": "rfft", "fft_selection": "interp"},
    "200 freq, direct": {"nfreq": 200},
    "200 freq, rfft": {"nfreq": 200, "transform": "rfft"},
    "rtbf, 19 ellipticities": {"components": "rayleigh_joint", "nfreq": 10,
                               "ellipticity_angle_size": 19},
    "rtbf, 61 ellipticities": {"components": "rayleigh_joint", "nfreq": 10},
}

for name, changes in cases.items():
    config = base.copy()
    config.update(**changes)
    start = time.perf_counter()
    Beamformer(array, config).process()
    print(f"{name:24s} {time.perf_counter() - start:6.2f} s")

# Repeated analyses (as in real time), with and without the cache.
for cache in [False, True]:
    beamformer = Beamformer(array, base, cache=cache)
    start = time.perf_counter()

    for _ in range(5):
        beamformer.process()

    print(f"5 analyses, cache={cache!s:5s}   "
          f"{time.perf_counter() - start:6.2f} s")
