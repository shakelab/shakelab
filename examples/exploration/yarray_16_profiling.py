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
Processing times.

The cost of the main choices of a configuration. Times depend on the
machine and on the data: compare them with each other.
"""

import time

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)

array = ArrayData("array.json")
config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=50,
    grid_min=100.0, grid_max=2000.0, grid_size=100, azimuth_size=180,
    window_type="fixed", window_length=2.0, windows_per_block=40,
)
cases = [
    {"method": "bartlett"},
    {"method": "music"},
    {"components": "three_component_separate"},
    {"transform": "rfft"},
    {"transform": "rfft", "fft_selection": "interp"},
    {"components": "rayleigh_joint", "ellipticity_angle_size": 19},
]

for changes in cases:
    start = time.perf_counter()
    Beamformer(array, config.copy(**changes)).process()
    print(changes, time.perf_counter() - start)

# Repeated analyses (as in real time): the cache keeps the steering
# vectors while coordinates and grid do not change.
beamformer = Beamformer(array, config, cache=True)
