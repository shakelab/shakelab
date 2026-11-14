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
The configuration.

Every parameter of the analysis has a default and the name of the
corresponding command-line option (``--grid-max`` is ``grid_max``).
Only frequencies and search grid are required.
"""

from shakelab.exploration.yarray import BeamformingConfig, SearchGrid

config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=12,              # log-spaced frequencies
    grid_min=100.0, grid_max=1500.0, grid_size=60,  # velocities (m/s)
    azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    method="capon",
    spectra=True,                               # keep the whole spectra
    maxima=2,
)
config.info()

# Changes are validated together; parameters are also attributes.
config.update(method="music", music_sources=2)
config.window_length = 4.0
bartlett = config.copy(method="bartlett")

# Groups of parameters can be given as objects or dictionaries.
config.frequencies = [2.0, 3.0, 5.0, 8.0]
config.grid = SearchGrid(0.0, 0.25, 51, quantity="wavenumber")
config.window = {"type": "cycles", "length": 20}

# Files: written, read, and read with changes (as on the command line).
config.write("capon.json")
config = BeamformingConfig("capon.json", method="bartlett")
