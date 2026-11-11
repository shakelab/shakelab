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
A frequency-wavenumber analysis.

Vertical Capon beamforming of the nine SESAME stations: the spectral
maxima give phase velocity and azimuth at each frequency and for each
covariance block.

Run (from any folder)::

    python examples/exploration/yarray_03_beamforming.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_dispersion, plot_spectrum

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

array = ArrayData(DATA / "sesame_array.json")
config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=12,
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    spectra=True, maxima=2,
)

beamformer = Beamformer(array, config)
result = beamformer.process()                 # window common to all stations
result.info()

# A part of the record only.
first = beamformer.process("2003-01-01T00:01:45", "2003-01-01T00:02:30")

for maximum in result.maxima_at(rank=1, frequency=5.0):
    print(f"{maximum.time}  {maximum.velocity:6.0f} m/s  "
          f"azimuth {maximum.azimuth:5.1f}  power {maximum.power:.3g}")

result.dispersion().info()

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_dispersion(result, ax=axes[0])
plot_spectrum(result, 5.0, block=0, ax=axes[1])
figure.tight_layout()
plt.show()
