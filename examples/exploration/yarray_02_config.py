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
The processing configuration.

BeamformingConfig collects every parameter of the analysis. It is
immutable: ``replace`` returns a modified copy. The same configuration
is written to and read from JSON, and used by the command line.

Run (from any folder; files are written to ./output)::

    python examples/exploration/yarray_02_config.py
"""

from pathlib import Path

from shakelab.exploration import yarray as ya

OUTPUT = Path("output")             # files written by the example
OUTPUT.mkdir(exist_ok=True)

# 1. A configuration from Python. Only frequencies and grid are needed;
#    every other parameter has a default.
config = ya.BeamformingConfig(
    # 12 frequencies from 3 to 10 Hz, logarithmically spaced
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
    # 60 velocities from 100 to 1500 m/s (log), 72 azimuths
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    # windows of 2 s, 30 windows per covariance block
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
    method="capon",
    output="both",              # keep maxima and full spectra
    maxima=2,                   # two maxima per spectrum
)
print(config)

# 2. Other ways of giving frequencies and grid.
explicit = ya.FrequencyConfig(values=[2.0, 3.0, 5.0, 8.0])
cycles = ya.WindowConfig(type="cycles", length=20)    # 20 cycles
slowness = ya.SearchGrid(0.0005, 0.01, 50, quantity="slowness")
wavenumber = ya.SearchGrid(0.0, 0.25, 51, quantity="wavenumber")

# 3. Modified copies.
music = config.replace(method="music")
longer = config.replace(window=config.window.replace(length=4.0))

# 4. JSON: write, read, and the equivalent dictionary.
config.write(OUTPUT / "capon.json", owrite=True)
again = ya.BeamformingConfig.read(OUTPUT / "capon.json")
print(again.to_dict())

# 5. Invalid combinations are refused when the configuration is built.
try:
    config.replace(components="rayleigh_joint", method="bartlett")
except ValueError as error:
    print(f"refused: {error}")
