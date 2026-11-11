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
covariance block (here a few blocks over 60 s).

Run (from any folder)::

    python examples/exploration/yarray_03_beamforming.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_dispersion, plot_spectrum

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array.json")

config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
    method="capon",
    output="both",
    maxima=2,
)

# 1. Process: a function for a single analysis...
result = ya.beamform(array, config)

#    ...or a Beamformer, for repeated analyses with the same
#    configuration (cache=True keeps the steering vectors).
beamformer = ya.Beamformer(config, cache=True)
result = beamformer.process(array)
result.info()

# 2. The maxima: rank 1 is the strongest of each spectrum.
for maximum in result.maxima_at(rank=1, frequency=5.0):
    print(f"{maximum.time}  {maximum.frequency:.2f} Hz  "
          f"{maximum.velocity:6.0f} m/s  azimuth {maximum.azimuth:5.1f}  "
          f"power {maximum.power:.3g}")

print(result.power_units)    # (unit)**2/Hz: data not corrected

# 3. Figures: maxima versus frequency and one spectrum.
figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_dispersion(result, ax=axes[0])
plot_spectrum(result, 5.0, block=0, ax=axes[1])
figure.tight_layout()
plt.show()
