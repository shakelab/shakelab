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

A BeamformingResult holds the maxima, the spectra (if requested), the
configuration and the stations. Results of consecutive data are
joined, saved to files and read back.

Run (from any folder; files are written to ./output)::

    python examples/exploration/yarray_09_results.py
"""

import warnings
from pathlib import Path

from shakelab.exploration import yarray as ya

OUTPUT = Path("output")             # files written by the example
OUTPUT.mkdir(exist_ok=True)

# Test data of the repository (the example runs from any folder).
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "exploration" / "yarray" / "data"

warnings.filterwarnings("ignore", message="STLO")

array = ya.read_array(DATA / "sesame_array.json")

config = ya.BeamformingConfig(
    frequencies=ya.FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
    grid=ya.SearchGrid(100.0, 1500.0, 60, azimuth_size=72),
    window=ya.WindowConfig(type="fixed", length=2.0,
                           windows_per_block=20),
    output="both",
    maxima=3,
)

result = ya.beamform(array, config)
result.info()

# 1. Selection of the maxima by rank, component and frequency.
maxima = result.maxima_at(rank=1, frequency=6.0)
first = maxima[0]
print(first.velocity, first.slowness, first.kx, first.ky, first.time)

# 2. The spectrum of a frequency: power of every block on the grid.
spectrum = result.spectrum(6.0)
print(spectrum.power.shape)        # (blocks, azimuths, radial values)
print(spectrum.azimuths[:4], spectrum.grid_values[:4])
print(spectrum.times)              # centres of the blocks

# 3. What the result knows about itself.
print(result.config.method, result.sids, result.coordinates.shape)
print(result.frequencies)          # analysis frequencies

# 4. Results of consecutive data join into one (blocks renumbered).
middle = array.starttime + 30.0
first_half = ya.ArrayData([
    ya.StationData(station.sid,
                   [record.extract(None, middle - 0.01)
                    for record in station.records],
                   station.coordinates)
    for station in array
])
second_half = ya.ArrayData([
    ya.StationData(station.sid,
                   [record.extract(middle, None)
                    for record in station.records],
                   station.coordinates)
    for station in array
])
first = ya.beamform(first_half, config)
second = ya.beamform(second_half, config)
joined = ya.BeamformingResult.concatenate([first, second])
joined.info()

# 5. Files (pickle; only for trusted files): one result, or successive
#    results appended to the same file.
ya.write_result(result, OUTPUT / "result.pkl", owrite=True)
ya.write_result(first, OUTPUT / "daily.pkl", owrite=True)
ya.write_result(second, OUTPUT / "daily.pkl", append=True)

result = ya.read_result(OUTPUT / "result.pkl")
results = ya.read_results(OUTPUT / "daily.pkl")              # list of results
joined = ya.read_result(OUTPUT / "daily.pkl", concatenate=True)
print(len(results), len(joined.maxima))
