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
Density of the spectral maxima.

The maxima of many blocks (a long record, the analyses of a monitor)
scatter around the dispersion curves. Their density at each frequency,
a mixture of Gaussian kernels, shows the curves; its summary gives
mode, median and a quantile interval, its peaks the candidate modes.

Run (from any folder)::

    python examples/exploration/yarray_10_density.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration import yarray as ya
from shakelab.exploration.yarray.statistics import maximum_density
from shakelab.graphics.yarray import plot_density

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array_9.json")

# Short blocks with large overlap: many blocks, many maxima.
config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=15),
    grid=ya.SearchGrid(100.0, 1500.0, 80, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=15, block_overlap=0.8),
    maxima=3,
)
result = ya.beamform(array, config)

# 1. Density in log-velocity of the maxima of all ranks, weighted by
#    their power relative to the strongest of their spectrum.
density = maximum_density(
    result,
    coordinate="log_velocity",       # or "wavenumber", "log_slowness"
    weighting="relative_power",      # or "uniform", "power"
)
density.info()

# 2. Mode, median and the 16-84 % interval at each frequency (in the
#    coordinate of the density: here the logarithm of the velocity).
summary = density.summary()
for frequency, mode in zip(summary.frequencies, np.exp(summary.mode)):
    print(f"{frequency:5.2f} Hz: mode {mode:5.0f} m/s")

# 3. Up to two peaks per frequency (candidate modes).
peaks = density.peaks(2, minimum_relative_density=0.2)
print(np.exp(peaks.values[:3]))

# 4. Only the strongest maxima, uniform weights, kernels as wide as two
#    steps of the grid.
strongest = maximum_density(result, ranks=[1], bandwidth_factor=2.0,
                            coordinate="log_velocity")

plot_density(density, peaks=peaks, summary=summary)
plt.show()
