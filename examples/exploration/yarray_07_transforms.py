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
Direct transform and RFFT.

The direct transform works at any frequency and with windows of a
number of cycles; the RFFT (fixed windows only) takes the nearest FFT
bins or interpolates, and is faster with many frequencies.

Run (from any folder)::

    python examples/exploration/yarray_07_transforms.py
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

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=40,
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
)
beamformer = Beamformer(ArrayData(DATA / "sesame_array.json"), config)

for transform, selection in [("direct", "nearest"), ("rfft", "nearest"),
                             ("rfft", "interp")]:
    beamformer.config.update(transform=transform, fft_selection=selection)
    start = time.perf_counter()
    result = beamformer.process()
    print(f"{transform:6s} {selection:8s} {len(result.frequencies):3d} "
          f"frequencies in {time.perf_counter() - start:.2f} s")

# RFFT with a range without nfreq: every bin (0.5 Hz with 2 s windows).
beamformer.config.update(frequencies={"minimum": 3.0, "maximum": 10.0})
print(beamformer.process().frequencies)

# Windows of 10 cycles of each frequency (direct transform).
beamformer.config.update(transform="direct", nfreq=20, window_type="cycles",
                         window_length=10)
print(len(beamformer.process().maxima))
