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

The same data processed with the three estimators: Bartlett is robust,
Capon gives sharper maxima (more windows per block than stations
needed), MUSIC the sharpest, without physical units.
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_spectrum

config = BeamformingConfig(
    frequencies=[6.0],
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    spectra=True,
)
beamformer = Beamformer(ArrayData("array.json"), config)
figure, axes = plt.subplots(1, 3, figsize=(17, 5))

for axis, method in zip(axes, ["bartlett", "capon", "music"]):
    beamformer.config.method = method
    result = beamformer.process()
    plot_spectrum(result, 6.0, ax=axis)

# Options: diagonal loading of Capon, sources excluded by MUSIC.
beamformer.config.update(method="capon", diagonal_loading=0.01)
beamformer.config.update(method="music", music_sources=2)

plt.show()
