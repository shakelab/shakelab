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
"""

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
    BeamformingResult,
)

array = ArrayData("array.json")
beamformer = Beamformer(array, BeamformingConfig("capon.json", maxima=3,
                                                 spectra=True))
result = beamformer.process()

# Maxima, spectra, dispersion curve.
best = result.maxima_at(rank=1, frequency=6.0)[0]
print(best.velocity, best.slowness, best.azimuth, best.time)
spectrum = result.spectrum(6.0)       # power (blocks, azimuths, radii)
result.dispersion().info()

# Successive windows joined into one result.
joined = beamformer.process("2026-11-01T10:00", "2026-11-01T10:10")
joined.append(beamformer.process("2026-11-01T10:10", "2026-11-01T10:20"))
joined.info()
print(joined.stations_at(6.0, block=0))

# Files: one result, or results accumulated in one file.
result.write("result.pkl")
beamformer.process().write("daily.pkl", append=True)
daily = BeamformingResult("daily.pkl")
