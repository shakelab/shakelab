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

Two options of the configuration:

- statistical_uncertainty="empirical" (Bartlett and Capon): standard
  error of the power at every grid point, from the spread of the
  windows of the block (corrected for their overlap);
- geometric_uncertainty="expected_coherence" (Bartlett): the power
  expected when the station coordinates are uncertain (Gaussian
  errors with a given covariance).

Both are uncertainties of the power, not of the position of a maximum.

Run (from any folder)::

    python examples/exploration/yarray_08_uncertainty.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration import yarray as ya
from shakelab.graphics.yarray import plot_spectrum
from shakelab.signals.io import reader

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")


config = ya.BeamformingConfig(
    frequencies=[6.0],
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=30),
    method="capon",
    output="both",
)

# 1. Empirical standard error of the power.
array = ya.read_array(DATA / "sesame_array_9.json")
result = ya.beamform(array,
                     config.replace(statistical_uncertainty="empirical"))
best = result.maxima_at(rank=1)[0]
print(f"power {best.power:.3g} +/- {best.power_standard_error:.2g}")

# 2. Uncertain coordinates: a variance of 1 m**2 (1 m of standard
#    deviation) for every station; a dictionary gives one value (or a
#    3 x 3 covariance) per station.
collection = reader(str(DATA / "sesame" / "*.sac"))
coordinates = {
    stream.sid.rsplit(".", 1)[0]: (
        stream.records()[0].head.extra["sac"]["USER7"],
        stream.records()[0].head.extra["sac"]["USER8"],
    )
    for stream in collection
}
uncertain = ya.ArrayData.from_collection(
    collection, coordinates=coordinates, coordinate_covariance=1.0,
)
expected = ya.beamform(uncertain, config.replace(
    method="bartlett", geometric_uncertainty="expected_coherence",
))
print(expected.maxima_at(rank=1)[0].velocity)

# 3. Power and its standard error.
figure, axes = plt.subplots(1, 2, figsize=(13, 5))
plot_spectrum(result, 6.0, ax=axes[0])
plot_spectrum(result, 6.0, quantity="standard_error", ax=axes[1])
figure.tight_layout()
plt.show()
