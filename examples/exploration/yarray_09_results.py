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
Working with the results.

A BeamformingResult holds maxima, spectra (if requested), configuration
and stations of one or more analyses. Results of successive windows
are joined, written to files and read back.

Run (from any folder; files are written to ./output)::

    python examples/exploration/yarray_09_results.py
"""

import warnings
from pathlib import Path

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
    BeamformingResult,
)

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

array = ArrayData(DATA / "sesame_array.json")
config = BeamformingConfig(
    fmin=3.0, fmax=10.0, nfreq=12,
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=20,
    spectra=True, maxima=3,
)
beamformer = Beamformer(array, config)
result = beamformer.process()

# Maxima, spectra, dispersion curve.
best = result.maxima_at(rank=1, frequency=6.0)[0]
print(best.velocity, best.slowness, best.azimuth, best.time)
spectrum = result.spectrum(6.0)
print(spectrum.power.shape)          # (blocks, azimuths, radial values)
result.dispersion().info()

# Successive windows joined into one result.
middle = array.starttime + 30.0
joined = beamformer.process(endtime=middle - array.delta)
joined.append(beamformer.process(starttime=middle))
joined.info()
print(joined.stations_at(6.0, block=0))

# Files: one result, or results accumulated in one file.
result.write("output/result.pkl", owrite=True)
joined.write("output/daily.pkl", owrite=True)
beamformer.process(starttime=middle).write("output/daily.pkl", append=True)

daily = BeamformingResult("output/daily.pkl")
print(daily.nblocks, len(daily.maxima))
