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
Processing times of the array analysis.

Times the main choices of a configuration on the nine SESAME stations
(three components, 60 s), so that their cost can be weighed:

A. estimators (Bartlett, Capon, MUSIC) and component modes;
B. Fourier transform: direct transform at each frequency versus the
   RFFT of each window (nearest bin or interpolation), as the number
   of frequencies grows, with the same fixed windows. With windows of
   2 s the bins are 0.5 Hz apart: the nearest selection keeps only the
   distinct bins (15 between 3 and 10 Hz), so its time stops growing,
   while the interpolation keeps every frequency;
C. Rayleigh three-component beamforming (RTBF) with the three
   transforms and the size of the ellipticity grid;
D. repeated analyses (real time) with and without the cache of the
   steering vectors;
E. planning: geometric response and Monte Carlo experiments.

Each case is run ``--repeat`` times and the best time is kept. The
default sizes are small, so that the example runs in a few seconds;
``--full`` uses sizes closer to real work (50 to 200 frequencies, 180
azimuths, 100 radial values; about a minute)::

    python examples/exploration/yarray_profiling.py --full --repeat 3

Run from the root of the repository::

    python examples/exploration/yarray_profiling.py
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np

from shakelab.exploration.yarray import (
    Beamformer,
    BeamformingConfig,
    EllipticityGrid,
    FrequencyConfig,
    MonteCarloConfig,
    PerformanceConfig,
    PlaneWaveSource,
    SearchGrid,
    WindowConfig,
    analyze_response,
    array_response,
    beamform,
    evaluate_performance,
    read_array,
)

warnings.filterwarnings("ignore", message="STLO")

parser = argparse.ArgumentParser(description="YArray processing times")
parser.add_argument("--full", action="store_true",
                    help="sizes closer to real work")
parser.add_argument("--repeat", type=int, default=1,
                    help="runs of each case (the best is kept)")
arguments, _ = parser.parse_known_args()

FULL = arguments.full
FREQUENCIES = (50, 100, 200) if FULL else (10, 20)
GRID = SearchGrid(100.0, 2000.0, 100 if FULL else 40,
                  azimuth_size=180 if FULL else 60)

array = read_array("tests/exploration/yarray/data/sesame_array_9.json")
rows = []


def measure(group: str, name: str, func) -> float:
    """Best time of ``func()`` over the repetitions; record a row."""
    best = np.inf

    for _ in range(max(arguments.repeat, 1)):
        start = time.perf_counter()
        func()
        best = min(best, time.perf_counter() - start)

    rows.append((group, name, best))
    print(f"{group}  {name:44s} {best:8.3f} s")
    return best


def config(size: int = 10, **options) -> BeamformingConfig:
    """Fixed windows of 2 s, Capon, ``size`` frequencies 3-10 Hz."""
    options.setdefault("window", WindowConfig(type="fixed", length=2.0,
                                              windows_per_block=40))
    return BeamformingConfig(
        frequencies=FrequencyConfig(minimum=3.0, maximum=10.0, size=size),
        grid=GRID,
        **options,
    )


# A. Estimators and component modes.
for method in ("bartlett", "capon", "music"):
    measure("A", f"{method}, vertical",
            lambda: beamform(array, config(method=method)))

for components in ("radial", "transverse", "three_component_separate"):
    measure("A", f"capon, {components}",
            lambda: beamform(array, config(components=components)))

# B. Direct transform and RFFT, as the frequencies grow.
for size in FREQUENCIES:
    direct = measure("B", f"{size} frequencies, direct",
                     lambda: beamform(array, config(size)))
    nearest = measure("B", f"{size} frequencies, rfft nearest",
                      lambda: beamform(array, config(size, transform="rfft")))
    measure("B", f"{size} frequencies, rfft interp",
            lambda: beamform(array, config(size, transform="rfft",
                                           fft_selection="interp")))
    kept = len(beamform(array, config(size, transform="rfft")).frequencies)
    print(f"   direct / rfft nearest: {direct / nearest:.2f} "
          f"({kept} distinct bins)")

# C. RTBF: transforms and ellipticity grid.
for transform in ("direct", "rfft"):
    measure("C", f"rtbf, {transform}, 61 ellipticities",
            lambda: beamform(array, config(components="rayleigh_joint",
                                           transform=transform)))

measure("C", "rtbf, direct, 19 ellipticities",
        lambda: beamform(array, config(
            components="rayleigh_joint",
            ellipticity_grid=EllipticityGrid(-85.0, 85.0, 19),
        )))

# D. Repeated analyses of the same array (as the monitor does).
windows = 10 if FULL else 4
for cache in (False, True):
    beamformer = Beamformer(config(), cache=cache)
    measure("D", f"{windows} analyses, cache={cache}",
            lambda: [beamformer.process(array) for _ in range(windows)])

# E. Planning.
measure("E", "geometric response and analysis",
        lambda: analyze_response(array_response(array)))
response = array_response(array)
experiment = PerformanceConfig(
    snr_db=0.0,
    grid=SearchGrid(0.0, 0.24, 41, quantity="wavenumber",
                    azimuth_size=72),
    monte_carlo=MonteCarloConfig(trials=200 if FULL else 50,
                                 realizations=20, seed=1),
)
measure("E", f"{experiment.monte_carlo.trials} Monte Carlo trials, capon",
        lambda: evaluate_performance(response, PlaneWaveSource(0.15, 60.0),
                                     experiment))

total = sum(row[2] for row in rows)
print(f"\n{len(rows)} cases, {total:.1f} s in total "
      f"({'full' if FULL else 'quick'} sizes, best of {arguments.repeat})")
