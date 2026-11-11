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

Measures the cost of the main choices of a configuration on the SESAME
test data. Times depend on the machine: compare them with each other.

Run (from any folder)::

    python examples/exploration/yarray_16_profiling.py
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
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=50),
    grid=ya.SearchGrid(100.0, 2000.0, 100, azimuth_size=180),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=40),
)

# 1. Estimators.
for method in ["bartlett", "capon", "music"]:
    start = time.perf_counter()
    ya.beamform(array, config.replace(method=method))
    print(f"{method:24s} {time.perf_counter() - start:6.2f} s")

# 2. Components.
for components in ["vertical", "radial", "three_component_separate"]:
    start = time.perf_counter()
    ya.beamform(array, config.replace(components=components))
    print(f"{components:24s} {time.perf_counter() - start:6.2f} s")

# 3. Transforms, with 50 and 200 frequencies.
for size in [50, 200]:
    frequencies = ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=size)

    for transform, selection in [("direct", "nearest"),
                                 ("rfft", "nearest"), ("rfft", "interp")]:
        start = time.perf_counter()
        ya.beamform(array, config.replace(
            frequencies=frequencies, transform=transform,
            fft_selection=selection,
        ))
        print(f"{size} frequencies, {transform} {selection:8s} "
              f"{time.perf_counter() - start:6.2f} s")

# 4. RTBF: the most expensive analysis (a system of dimension 18 for
#    each frequency, block, azimuth and ellipticity).
rtbf = config.replace(
    components="rayleigh_joint",
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=10),
)
for size in [19, 61]:
    start = time.perf_counter()
    ya.beamform(array, rtbf.replace(
        ellipticity_grid=ya.EllipticityGrid(-85.0, 85.0, size)
    ))
    print(f"rtbf, {size} ellipticities    "
          f"{time.perf_counter() - start:6.2f} s")

# 5. Repeated analyses of the same array, with and without the cache
#    of the steering vectors (as in real time).
for cache in [False, True]:
    beamformer = ya.Beamformer(config, cache=cache)
    start = time.perf_counter()
    for _ in range(5):
        beamformer.process(array)
    print(f"5 analyses, cache={cache!s:5s}   "
          f"{time.perf_counter() - start:6.2f} s")
