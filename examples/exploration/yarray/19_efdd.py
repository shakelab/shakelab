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
Enhanced FDD: natural frequency and damping of the modes of a building.

The FDD gives the resonances as peaks of the first singular value.
The enhanced FDD takes, around each peak, the lines whose shape matches
the mode (the bell), transforms them to the free decay of the mode and
measures its period and logarithmic decrement.
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import ArrayData, FDD, FDDConfig
from shakelab.graphics.yarray import plot_efdd
from shakelab.signals.io import reader

records = reader("building/*.mseed", metadata="building.xml")
array = ArrayData(records, select="XX.B*..HN?")   # one station per floor

config = FDDConfig(
    fmin=0.5, fmax=15.0,
    window_type="fixed", window_length=20.0,      # seconds
    windows_per_block=10,
    components="east",                            # or "north", ...
)

fdd = FDD(array, config)
result = fdd.process()
peaks = result.peaks()

# One mode per peak, from a local spectrum of each mode.
modes = fdd.efdd(
    peaks,
    mac=0.8,                # lines of the bell: MAC with the peak shape
    cycles=100,             # windows of 100 cycles (or duration=60.0)
    levels=(0.1, 0.9),      # extrema of the correlation used
)

for mode in modes:
    mode.info()             # frequency, damping, shape

print(modes[0].efdd.decrement, modes[0].efdd.damped_frequency)

# Other figures: plot_mode_shape(modes[0], along="z"), plot_mac(modes).
plot_efdd(modes[0])
plt.show()
