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

The maxima of many blocks scatter around the dispersion curves; their
density at each frequency shows the curves, its summary gives mode,
median and a quantile interval, its peaks the candidate modes.

Run (from any folder)::

    python examples/exploration/yarray_10_density.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_density

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=15,
    grid_min=100.0, grid_max=1500.0, grid_size=80, azimuth_size=72,
    window_type="fixed", window_length=2.0,
    windows_per_block=15, block_overlap=0.8,   # many blocks
    maxima=3,
)
result = Beamformer(ArrayData(DATA / "sesame_array.json"), config).process()

density = result.density(
    coordinate="log_velocity",          # or "wavenumber", "log_slowness"
    weighting="relative_power",         # or "uniform", "power"
)
density.info()

summary = density.summary()
peaks = density.peaks(2, minimum_relative_density=0.2)

for frequency, mode in zip(summary.frequencies, np.exp(summary.mode)):
    print(f"{frequency:5.2f} Hz: mode {mode:5.0f} m/s")

plot_density(density, peaks=peaks, summary=summary)
plt.show()
