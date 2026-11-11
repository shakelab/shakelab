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

For each tested azimuth the east and north records are projected on
the radial direction (along the propagation) and on the transverse one
(radial rotated by 90 degrees clockwise). Radial and vertical motion
carry Rayleigh waves, transverse motion carries Love waves.

Run (from any folder)::

    python examples/exploration/yarray_05_components.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_dispersion

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array_9.json")

config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
    method="capon",
)

# 1. One component at a time.
radial = ya.beamform(array, config.replace(components="radial"))
transverse = ya.beamform(array, config.replace(components="transverse"))

# 2. The three components in one run, on the same grid.
three = ya.beamform(
    array, config.replace(components="three_component_separate")
)
print(three.components)          # ('vertical', 'radial', 'transverse')

for component in three.components:
    maxima = three.maxima_at(rank=1, component=component, frequency=5.0)
    print(component, [round(item.velocity) for item in maxima])

# 3. One figure per component.
figure, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)

for axis, component in zip(axes, three.components):
    plot_dispersion(three, component=component, ax=axis)

figure.tight_layout()
plt.show()
