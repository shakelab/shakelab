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
Operational assessment of an array.

Where will the f-k analysis be reliable? The performance experiment is
repeated over a grid of frequencies, velocities and azimuths, and each
point is marked usable when it meets explicit criteria; otherwise it
gets the factor that limits it. Recommendations summarize the result.

Run (from any folder)::

    python examples/exploration/yarray_15_operational.py
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    ArrayResponse,
    PerformanceConfig,
)
from shakelab.graphics.yarray import plot_operational

DATA = Path(__file__).resolve().parents[2] / "tests/exploration/yarray/data"
warnings.filterwarnings("ignore", message="STLO")  # local SAC coordinates

response = ArrayResponse(ArrayData(DATA / "sesame_array.json"))
config = PerformanceConfig(snr_db=-3.0, method="capon",
                           grid_min=0.0, grid_max=0.24, grid_size=41,
                           azimuth_size=72, trials=50, realizations=20,
                           seed=1)

assessment = response.assess(config,
                             frequencies=[3, 5, 8, 12],
                             velocities=[150, 250, 400, 600, 1000],
                             azimuths=[0, 60, 120],  # default: 12
                             min_local=0.8, max_velocity_rmse=0.15)
assessment.info()                        # usable velocities, advice

plot_operational(assessment)
plt.show()
