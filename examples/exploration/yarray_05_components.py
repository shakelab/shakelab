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
Horizontal components.

For each azimuth of the grid the horizontal motion is projected on
the propagation direction (radial: Rayleigh waves) and normal to it
(transverse: Love waves).

Run (from any folder)::

    python examples/exploration/yarray_05_components.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_dispersion

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=12,
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    components="three_component_separate",   # or radial, transverse
)
result = Beamformer(ArrayData(DATA / "sesame_array.json"), config).process()
print(result.components)

figure, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)

for axis, component in zip(axes, result.components):
    result.dispersion(component=component).info()
    plot_dispersion(result, component=component, ax=axis)

figure.tight_layout()
plt.show()
