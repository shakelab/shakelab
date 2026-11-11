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
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_dispersion

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=12,
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    components="three_component_separate",   # or radial, transverse
)
result = Beamformer(ArrayData("array.json"), config).process()

for component in result.components:      # vertical, radial, transverse
    result.dispersion(component=component).info()

plot_dispersion(result, component="radial")
plt.show()
