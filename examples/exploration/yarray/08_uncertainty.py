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
Uncertainty of the beam power.

``statistical_uncertainty="empirical"`` (Bartlett and Capon): standard
error of the power from the spread of the windows of a block.
``geometric_uncertainty="expected_coherence"`` (Bartlett): power
expected with uncertain station coordinates.
"""

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_spectrum
from shakelab.signals.io import metadata_reader, reader

config = BeamformingConfig(
    frequencies=[6.0],
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    spectra=True, statistical_uncertainty="empirical",
)
result = Beamformer(ArrayData("array.json"), config).process()
best = result.maxima_at(rank=1)[0]
print(best.power, best.power_standard_error)

# Coordinates with a variance of 1 m**2 (one value, or a 2 x 2 or
# 3 x 3 covariance matrix, also per station).
collection = reader("data/*.mseed", metadata=metadata_reader("array.xml"))
uncertain = ArrayData(collection, coordinate_covariance=1.0)
config.update(method="bartlett", statistical_uncertainty="none",
              geometric_uncertainty="expected_coherence")
expected = Beamformer(uncertain, config).process()

plot_spectrum(result, 6.0, quantity="standard_error")
plt.show()
