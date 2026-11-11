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

The Fourier coefficients of the windows are computed in two ways:

- direct: the transform at each requested frequency (any frequency,
  windows of fixed duration or of a number of cycles);
- RFFT: one FFT per window (fixed windows only), taking the nearest
  bin (the frequencies become the bins, repeated bins are processed
  once) or interpolating at the requested frequencies.

Run (from any folder)::

    python examples/exploration/yarray_07_transforms.py
"""

import time
import warnings
from pathlib import Path

from shakelab.exploration import yarray as ya

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array.json")

config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=40),
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
)

for transform, selection in [("direct", "nearest"), ("rfft", "nearest"),
                             ("rfft", "interp")]:
    start = time.perf_counter()
    result = ya.beamform(array, config.replace(transform=transform,
                                               fft_selection=selection))
    elapsed = time.perf_counter() - start
    print(f"{transform:6s} {selection:8s} {len(result.frequencies):3d} "
          f"frequencies in {elapsed:.2f} s; first: "
          f"{result.frequencies[:3].round(3)}")

# With 2 s windows the FFT bins are 0.5 Hz apart: the nearest
# selection keeps 15 distinct frequencies between 3 and 10 Hz.

# With the RFFT a range without size takes every bin in the range.
bins = config.replace(
    transform="rfft",
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0),
)
print(ya.beamform(array, bins).frequencies)

# Windows of a number of cycles need the direct transform.
cycles = config.replace(window=ya.WindowConfig(type="cycles", length=10,
                                               windows_per_block=30))
print(len(ya.beamform(array, cycles).maxima))
