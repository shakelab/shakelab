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
Frequency domain decomposition (FDD).

The cross-spectral matrix of all the channels is decomposed at each
frequency: the peaks of the first singular value are the resonances
(modes of a structure, two-dimensional resonances of a valley) and the
first singular vector is the shape of each mode.
"""

from shakelab.exploration.yarray import (
    ArrayData,
    FDD,
    FDDConfig,
    FDDResult,
    mac,
)

array = ArrayData("array.json")

config = FDDConfig(
    fmin=0.1, fmax=1.0, nfreq=100,
    frequency_scale="linear",           # default: "log"
    window_type="cycles",               # or "fixed" (seconds)
    window_length=20,
    windows_per_block=20,
    subfrequencies=5, bandwidth=0.05,   # default: 1, 0 (no averaging)
    normalization="coherence",          # default: "none"
    components="transverse",            # or "vertical", "radial", ...
    azimuth=30.0,                       # e.g. the axis of a valley
)
config.info()

fdd = FDD(array, config)
result = fdd.process()
result.info()

# Singular values: geometric mean and logarithmic spread over blocks.
frequency, mean, log_std = result.singular_values()

# Resonances and their shapes.
peaks = result.peaks(fmin=0.2, fmax=0.8)      # prominence 2 by default
first = result.mode(peaks[0])
first.info()
second = result.mode(peaks[1])
print(mac(first.shape, second.shape))         # 1: same shape, 0: orthogonal

# Files: results of successive windows accumulated in one file.
result.write("fdd.pkl")
fdd.process("2026-11-01T10:00", "2026-11-01T11:00").write(
    "fdd.pkl", append=True)
joined = FDDResult("fdd.pkl")
