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
Bartlett, Capon and MUSIC.

The same data and grid processed with the three estimators:

- Bartlett: conventional beamforming, robust, resolution of the array;
- Capon: minimum variance, sharper maxima (needs more windows per
  block than stations);
- MUSIC: noise-subspace pseudo-spectrum, sharpest, without physical
  units.

Run (from any folder)::

    python examples/exploration/yarray_04_methods.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_spectrum

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array_9.json")

config = ya.BeamformingConfig(
    frequencies=[6.0],
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
    output="both",
)

figure, axes = plt.subplots(1, 3, figsize=(17, 5))

for axis, method in zip(axes, ["bartlett", "capon", "music"]):
    result = ya.beamform(array, config.replace(method=method))
    best = result.maxima_at(rank=1)[0]
    print(f"{method:9s} {best.velocity:6.0f} m/s, azimuth "
          f"{best.azimuth:5.1f}, power {best.power:.3g} "
          f"{result.power_units}")
    plot_spectrum(result, 6.0, block=0, ax=axis)
    axis.set_title(method)

# Capon options: relative diagonal loading (stability) and the check
# on the number of windows (strict=True by default).
capon = config.replace(method="capon", diagonal_loading=0.01)

# MUSIC option: number of sources excluded from the noise subspace.
music = config.replace(method="music", music_sources=2)

figure.tight_layout()
plt.show()
