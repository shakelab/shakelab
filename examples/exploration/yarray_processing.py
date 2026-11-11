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
Array processing as a library: every option of the beamforming.

Nine three-component stations of the SESAME ambient-vibration
simulation (60 s) are processed with the Python interface of
:mod:`shakelab.exploration.yarray`:

1. three ways of building the array (array file; waveforms read by
   ShakeLab with a table of coordinates; records already in memory);
2. a configuration built in Python, saved to JSON, read back and
   modified with ``replace``;
3. the three estimators (Bartlett, Capon, MUSIC) on the vertical
   component;
4. the component modes: radial, transverse, the three components
   separately and the Rayleigh three-component beamforming (RTBF) with
   the signed ellipticity;
5. the uncertainties: empirical standard error of the power and
   expected power over uncertain station coordinates;
6. the result: maxima, spectra, units, joining, files and the density
   of the maxima.

The data are not corrected for the instrument, so the power is in
(unit)**2/Hz; every step works the same with physical units.

Run from the root of the repository::

    python examples/exploration/yarray_processing.py
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration.yarray import (
    ArrayData,
    BeamformingConfig,
    BeamformingResult,
    EllipticityGrid,
    FrequencyConfig,
    SearchGrid,
    WindowConfig,
    beamform,
    read_array,
    read_results,
    write_result,
)
from shakelab.exploration.yarray.statistics import maximum_density
from shakelab.graphics.yarray import (
    plot_density,
    plot_dispersion,
    plot_ellipticity,
    plot_spectrum,
)
from shakelab.signals.io import reader

DATA = Path("tests/exploration/yarray/data")

# The SAC files keep the local coordinates also in STLA/STLO.
warnings.filterwarnings("ignore", message="STLO")

# ----------------------------------------------------------------------------
# 1. The array
# ----------------------------------------------------------------------------

# a) Array file: local coordinates and one waveform file per component.
array = read_array(DATA / "sesame_array_9.json")
array.info()

# b) Waveforms read by ShakeLab and a table of coordinates (metres, in
#    a local frame). With geographic positions in the metadata (e.g.
#    StationXML) the coordinates are computed in a local ENU frame and
#    the table is not needed.
collection = reader(str(DATA / "sesame" / "*.sac"))
table = {
    stream.sid.rsplit(".", 1)[0]: (
        stream.records()[0].head.extra["sac"]["USER7"],
        stream.records()[0].head.extra["sac"]["USER8"],
    )
    for stream in collection
}
from_files = ArrayData.from_collection(collection, coordinates=table,
                                       select="*.S10*..*")
print(from_files.sids)

# c) Records already in memory (here those of five stations); the
#    keys of the coordinates are the NET.STA.LOC of the records.
from_records = ArrayData.from_records(
    [record for station in array.stations[:5]
     for record in station.records],
    coordinates=table,
)
print(from_records.sids)

# ----------------------------------------------------------------------------
# 2. The configuration
# ----------------------------------------------------------------------------

config = BeamformingConfig(
    frequencies=FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
    grid=SearchGrid(100.0, 1500.0, 60, azimuth_size=72),  # velocity, log
    window=WindowConfig(type="fixed", length=2.0, windows_per_block=30),
    method="capon",
    output="both",              # maxima and full spectra
    maxima=2,                   # two maxima per block
)

with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "capon.json"
    config.write(path)          # the same file is used by slab-yarray
    config = BeamformingConfig.read(path)

print(config.to_dict()["window"])

# ----------------------------------------------------------------------------
# 3. Estimators on the vertical component
# ----------------------------------------------------------------------------

results = {
    method: beamform(array, config.replace(method=method))
    for method in ("bartlett", "capon", "music")
}

for method, result in results.items():
    best = result.maxima_at(rank=1, frequency=5.0)[0]
    print(f"{method:9s} 5 Hz: {best.velocity:6.0f} m/s, azimuth "
          f"{best.azimuth:5.1f}, power {best.power:.3g} "
          f"{result.power_units or ''}")

# ----------------------------------------------------------------------------
# 4. Component modes
# ----------------------------------------------------------------------------

# More windows per block: the radial-vertical matrix of the RTBF has
# dimension 2 x 9 and Capon needs more windows than that.
three = beamform(array, config.replace(
    components="three_component_separate", maxima=1,
))
print(three.components)

rtbf = beamform(array, config.replace(
    components="rayleigh_joint",
    maxima=1,
    window=config.window.replace(windows_per_block=40),
    ellipticity_grid=EllipticityGrid(-85.0, 85.0, 35),
))

for maximum in rtbf.maxima_at(rank=1)[:4]:
    print(f"{maximum.frequency:5.2f} Hz {maximum.velocity:6.0f} m/s "
          f"ellipticity {maximum.ellipticity:+.2f} "
          f"({maximum.ellipticity_angle:+.0f} degrees)")

# ----------------------------------------------------------------------------
# 5. Uncertainties
# ----------------------------------------------------------------------------

# Empirical standard error of the power (Bartlett and Capon).
empirical = beamform(array, config.replace(
    statistical_uncertainty="empirical",
))
best = empirical.maxima_at(rank=1)[0]
print(f"power {best.power:.3g} +/- {best.power_standard_error:.2g}")

# Expected power over uncertain coordinates (Bartlett): 1 m standard
# deviation for every station (variance in square metres).
uncertain = ArrayData.from_collection(
    collection, coordinates=table, coordinate_covariance=1.0,
)
expected = beamform(uncertain, config.replace(
    method="bartlett", geometric_uncertainty="expected_coherence",
))
print(expected.maxima_at(rank=1)[0])

# ----------------------------------------------------------------------------
# 6. Results
# ----------------------------------------------------------------------------

result = results["capon"]
result.info()
spectrum = result.spectrum(6.0)       # nearest analysis frequency
print(spectrum.power.shape)           # (blocks, azimuths, radial values)

#    Results of consecutive data join into one (blocks renumbered);
#    a file can collect the results of successive analyses. Two halves
#    of 30 s hold 29 windows of 2 s: blocks of 20 windows.
halves = config.replace(window=config.window.replace(windows_per_block=20))
first, second = (
    beamform(ArrayData([
        type(station)(station.sid,
                      [record.extract(start, end)
                       for record in station.records],
                      station.coordinates)
        for station in array
    ]), halves)
    for start, end in ((None, array.starttime + 29.99),
                       (array.starttime + 30.0, None))
)

with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "results.pkl"
    write_result(first, path)
    write_result(second, path, append=True)
    joined = BeamformingResult.concatenate(read_results(path))

joined.info()

#    Density of the maxima: dispersion curve and its modes.
density = maximum_density(joined, coordinate="log_velocity")
density.info()
summary = density.summary()      # in the coordinate: log(velocity)

for frequency, mode in zip(summary.frequencies, np.exp(summary.mode)):
    print(f"{frequency:5.2f} Hz: mode {mode:6.0f} m/s")

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

figure, axes = plt.subplots(2, 3, figsize=(17, 10))

for axis, (method, item) in zip(axes[0], results.items()):
    plot_spectrum(item, 6.0, ax=axis)
    axis.set_title(f"{method}, 6 Hz, block 0")

plot_dispersion(three, component="radial", ax=axes[1, 0])
plot_ellipticity(rtbf, ax=axes[1, 1])
plot_density(density, peaks=density.peaks(2), ax=axes[1, 2])
figure.tight_layout()
plt.show()
