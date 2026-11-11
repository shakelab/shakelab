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
Rayleigh three-component beamforming (RTBF).

Radial and vertical coefficients are processed together (a matrix of
dimension 2 x stations) with steering vectors ``[-i e a, a]``: the
spectral maxima give the phase velocity of Rayleigh waves and their
signed ellipticity ``e`` (radial over vertical amplitude; the sign
gives the sense of the particle motion). Wathelet et al. (2018).

Run (from any folder)::

    python examples/exploration/yarray_06_rtbf.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_dispersion, plot_ellipticity

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array_9.json")

config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=10),
    grid=ya.SearchGrid(100.0, 1500.0, 40, azimuth_size=48),
    # Capon needs more windows per block than the matrix dimension:
    # 2 x 9 stations = 18.
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=40),
    method="capon",                    # RTBF works with Capon only
    components="rayleigh_joint",
    # Ellipticity searched in its angle atan(e), from -85 to 85 degrees.
    ellipticity_grid=ya.EllipticityGrid(-85.0, 85.0, 35),
)

result = ya.beamform(array, config)
result.info()

for maximum in result.maxima_at(rank=1):
    print(f"{maximum.frequency:5.2f} Hz  {maximum.velocity:6.0f} m/s  "
          f"ellipticity {maximum.ellipticity:+.2f} "
          f"({maximum.ellipticity_angle:+.0f} degrees)")

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_dispersion(result, ax=axes[0])
plot_ellipticity(result, ax=axes[1])
figure.tight_layout()
plt.show()
