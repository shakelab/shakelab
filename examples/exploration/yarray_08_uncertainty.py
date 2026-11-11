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

Run (from any folder)::

    python examples/exploration/yarray_08_uncertainty.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_spectrum
from shakelab.signals.io import reader

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

config = BeamformingConfig(
    frequencies=[6.0],
    grid_min=100.0, grid_max=1500.0, grid_size=60, azimuth_size=72,
    window_type="fixed", window_length=2.0, windows_per_block=30,
    spectra=True, statistical_uncertainty="empirical",
)
result = Beamformer(ArrayData(DATA / "sesame_array.json"), config).process()
best = result.maxima_at(rank=1)[0]
print(f"power {best.power:.3g} +/- {best.power_standard_error:.2g}")

# Coordinates with a variance of 1 m**2 (a dictionary gives one value,
# or a 3 x 3 covariance, per station).
collection = reader(str(DATA / "sesame" / "*.sac"))
coordinates = {
    stream.sid.rsplit(".", 1)[0]: (
        stream.records()[0].head.extra["sac"]["USER7"],
        stream.records()[0].head.extra["sac"]["USER8"],
    )
    for stream in collection
}
uncertain = ArrayData(collection, coordinates=coordinates,
                      coordinate_covariance=1.0)
config.update(method="bartlett", statistical_uncertainty="none",
              geometric_uncertainty="expected_coherence")
print(Beamformer(uncertain, config).process().maxima_at(rank=1)[0].velocity)

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_spectrum(result, 6.0, ax=axes[0])
plot_spectrum(result, 6.0, quantity="standard_error", ax=axes[1])
figure.tight_layout()
plt.show()
