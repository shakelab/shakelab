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
Frequency-wavenumber analysis of a small three-component array.

Nine stations of the SESAME ambient-vibration simulation (60 s, local
coordinates in an array file) are processed in two ways:

1. vertical Capon beamforming, with fixed windows and the RFFT;
2. Rayleigh three-component beamforming (RTBF), which also gives the
   signed ellipticity.

The same configuration is then applied window by window as in real
time, with an :class:`ArrayMonitor` fed directly (the SeedLink client
calls it in the same way); the results are joined and the density of
their maxima shows the dispersion curve.

Run from the root of the repository::

    python examples/exploration/yarray_sesame.py
"""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    Beamformer,
    BeamformingConfig,
    EllipticityGrid,
    FrequencyConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.io import read_array
from shakelab.exploration.yarray.monitor import ArrayMonitor
from shakelab.exploration.yarray.statistics import maximum_density
from shakelab.graphics.yarray import (
    plot_density,
    plot_dispersion,
    plot_ellipticity,
    plot_spectrum,
)
from shakelab.signals.base import StreamCollection

ARRAY = "tests/exploration/yarray/data/sesame_array_9.json"

# The SAC files keep the local coordinates also in STLA/STLO.
warnings.filterwarnings("ignore", message="STLO")

# 1. Read the array: coordinates and files of the components.
array = read_array(ARRAY)
array.info()

# 2. Vertical Capon beamforming with fixed windows of 2 s.
vertical = BeamformingConfig(
    frequencies=FrequencyConfig(minimum=3.0, maximum=12.0, size=20),
    grid=SearchGrid(100.0, 2000.0, 80, azimuth_size=90),
    window=WindowConfig(type="fixed", length=2.0, windows_per_block=20),
    method="capon",
    transform="rfft",
    output="both",
    maxima=2,
)
result = Beamformer(vertical).process(array)
result.info()

# 3. Rayleigh three-component beamforming (signed ellipticity), on a
#    coarser grid: each frequency, block and azimuth needs a matrix of
#    dimension 18 and steering vectors for every ellipticity, the most
#    expensive analysis (see yarray_profiling.py).
rtbf = vertical.replace(
    components="rayleigh_joint", maxima=1,
    frequencies=FrequencyConfig(minimum=3.0, maximum=12.0, size=10),
    grid=SearchGrid(100.0, 2000.0, 40, azimuth_size=48),
    ellipticity_grid=EllipticityGrid(-85.0, 85.0, 31),
)
joint = Beamformer(rtbf).process(array)
joint.info()

figure, axes = plt.subplots(1, 3, figsize=(16, 5))
plot_dispersion(result, ax=axes[0])
plot_spectrum(result, 6.0, block=0, ax=axes[1])
plot_ellipticity(joint, ax=axes[2])
figure.tight_layout()

# 4. The same analysis window by window, as in real time: the monitor
#    receives the data of the array (here the records of the file) and
#    keeps the results, which can be queried and joined.
#    The records keep the identifiers of the SAC files (".S1001..Z").
collection = StreamCollection()

for station in array:
    for record in station.records:
        collection.append(record.copy())

#    Windows of 20 s hold 19 windows of 2 s: blocks of 10 windows.
short = vertical.replace(
    window=vertical.window.replace(windows_per_block=10),
)
monitor = ArrayMonitor(
    Beamformer(short, cache=True), window=20.0, min_stations=3,
    coordinates={f".{st.sid}.": st.coordinates for st in array},
    background=False,
)
start = array.starttime

for offset in (0.0, 20.0, 40.0):
    monitor.process(collection, start + offset, start + offset + 19.99)

print(monitor.status())
joined = monitor.result()
joined.info()

# 5. Density of the maxima of all the windows, with its peaks.
density = maximum_density(joined, coordinate="log_slowness")
plot_density(density, peaks=density.peaks(3, minimum_relative_density=0.2))

plt.show()
