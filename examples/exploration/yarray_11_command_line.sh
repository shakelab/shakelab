#!/bin/sh
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
#
# Array processing from the command line.
#
# slab-yarray is installed with ShakeLab (pip install -e .);
# "python -m shakelab.exploration.yarray" is the same command.
# Run from the root of the repository:
#
#     sh examples/exploration/yarray_11_command_line.sh

set -e

ARRAY=tests/exploration/yarray/data/sesame_array.json
mkdir -p output

# 1. A configuration file (written here with the shell; the section
#    "graphics" is used only by the command line).
cat > output/capon.json << 'END'
{
  "frequencies": {"minimum": 3.0, "maximum": 10.0, "size": 12},
  "grid": {"minimum": 100.0, "maximum": 1500.0, "size": 60,
           "azimuth_size": 72},
  "window": {"type": "fixed", "length": 2.0, "windows_per_block": 30},
  "method": {"name": "capon"},
  "output": {"spectra": true, "maxima": 2},
  "graphics": {"dispersion": true, "ranks": [1]}
}
END

# 2. Process an array file with the configuration; save the result and
#    the figures of the section "graphics".
slab-yarray --input $ARRAY --config output/capon.json \
    -o output/capon.pkl --save-figures output/figures --owrite

# 3. Options replace the values of the file: MUSIC on 5 frequencies.
slab-yarray --input $ARRAY --config output/capon.json \
    --method music --fmin 4 --fmax 8 --nfreq 5 -o output/music.pkl --owrite

# 4. Rayleigh three-component beamforming; write the configuration
#    actually used (file plus options) for later runs.
slab-yarray --input $ARRAY --config output/capon.json \
    --components rayleigh_joint --windows-per-block 40 \
    --ellipticity-angle-size 35 --maxima 1 \
    --write-config output/rtbf.json -o output/rtbf.pkl --owrite \
    --plot-ellipticity --save-figures output/figures

# 5. Successive analyses (e.g. one file per day) appended to one file.
slab-yarray --input $ARRAY --config output/capon.json -o output/daily.pkl \
    --append
slab-yarray --input $ARRAY --config output/capon.json -o output/daily.pkl \
    --append

# 6. Waveform files with their metadata (positions and orientations),
#    selecting channels and a time window:
#
#    slab-yarray --data "data/*.mseed" --metadata array.xml \
#        --select "XX.*..HH?" --start 2026-11-01T00:00 \
#        --end 2026-11-01T01:00 --config output/capon.json -o hour.pkl

# 7. Plot a saved result without processing: one spectrum, then the
#    density of the maxima; --show opens the figures on screen.
slab-yarray --result output/capon.pkl --plot-spectrum --frequency 6 \
    --spectrum-view azimuth --save-figures output/figures
slab-yarray --result output/capon.pkl --plot-density \
    --density-coordinate log_velocity --save-figures output/figures

# 8. Geometric response of the array: resolution and aliasing limits,
#    velocity limits from 2 to 20 Hz, a search domain and all figures.
slab-yarray response --input $ARRAY --fmin 2 --fmax 20 --nfreq 5 \
    --search-max 0.1 --plot-all --save-figures output/figures

# 9. The same for a layout still to be deployed: a text file with one
#    station per line, "name east north" in metres.
cat > output/layout.txt << 'END'
C    0.0    0.0
A1  15.0    0.0
A2  -7.5   13.0
A3  -7.5  -13.0
B1   0.0   40.0
B2  34.6  -20.0
B3 -34.6  -20.0
END
slab-yarray response --coordinates output/layout.txt --plot-response \
    --save-figures output/figures/layout

# Every option: slab-yarray --help, slab-yarray response --help
