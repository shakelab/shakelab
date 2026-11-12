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

Radial and vertical motion processed together give the phase velocity
of Rayleigh waves and their signed ellipticity (radial over vertical
amplitude; positive for prograde, negative for retrograde particle
motion). Wathelet et al. (2018).
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_dispersion, plot_ellipticity

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=10,
    grid_min=100.0, grid_max=1500.0, grid_size=40, azimuth_size=48,
    window_type="fixed", window_length=2.0,
    windows_per_block=40,                 # more than 2 x stations
    components="rayleigh_joint",          # Capon only
    ellipticity_angle_size=35,            # atan(e) from -85 to 85 degrees
)
result = Beamformer(ArrayData("array.json"), config).process()
result.dispersion().info()

# Density of the maxima in ellipticity angle, e.g. of one mode
density = result.density(
    coordinate="ellipticity_angle",
    weighting="relative_power",
    velocity=(150.0, 400.0),              # default: all the maxima
)
density.info()

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_dispersion(result, ax=axes[0])
plot_ellipticity(result, ax=axes[1])
plt.show()
