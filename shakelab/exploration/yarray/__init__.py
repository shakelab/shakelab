# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Seismic-array processing (YArray).

Frequency-wavenumber analysis of synchronized seismic arrays, usable
from Python, from the command line and within the callbacks of the
SeedLink client: array data (:mod:`.data`), configuration
(:mod:`.config`), estimators (:mod:`.estimators`), beamforming
(:mod:`.beamforming`), statistics of the maxima (:mod:`.statistics`),
input and output (:mod:`.io`), real-time
monitor (:mod:`.monitor`) and command line (:mod:`.cli`, installed as
``slab-yarray``). Figures are in :mod:`shakelab.graphics.yarray`. See
``README.md`` in this folder.
"""

from __future__ import annotations

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
    BeamformingWarning,
    FrequencySpectrum,
    SpectralMaximum,
    beamform,
)
from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    EllipticityGrid,
    FrequencyConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.data import (
    ArrayData,
    ArrayDataWarning,
    StationData,
)
from shakelab.exploration.yarray.io import (
    read_array,
    read_result,
    read_results,
    write_result,
)
from shakelab.exploration.yarray.monitor import ArrayMonitor

__all__ = [
    "ArrayData",
    "ArrayDataWarning",
    "ArrayMonitor",
    "Beamformer",
    "BeamformingConfig",
    "BeamformingResult",
    "BeamformingWarning",
    "EllipticityGrid",
    "FrequencyConfig",
    "FrequencySpectrum",
    "SearchGrid",
    "SpectralMaximum",
    "StationData",
    "WindowConfig",
    "beamform",
    "read_array",
    "read_result",
    "read_results",
    "write_result",
]
