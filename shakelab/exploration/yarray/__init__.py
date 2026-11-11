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
geometric response and limits (:mod:`.response`), Monte Carlo
performance of the estimators (:mod:`.performance`), operational
assessment for survey planning (:mod:`.operational`), file formats
(:mod:`.io`), real-time monitor (:mod:`.monitor`) and command line
(:mod:`.cli`, installed as ``slab-yarray``). Figures are in
:mod:`shakelab.graphics.yarray`. See ``README.md`` in this folder.

.. code-block:: python

    from shakelab.exploration.yarray import (
        ArrayData,
        ArrayResponse,
        Beamformer,
    )

    array = ArrayData("array.json")
    result = Beamformer(array, "capon.json").process()
    result.dispersion().info()

    ArrayResponse(array).analyze().info()
"""

from __future__ import annotations

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
    BeamformingWarning,
    DispersionCurve,
    FrequencySpectrum,
    SpectralMaximum,
)
from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    EllipticityGrid,
    FrequencyConfig,
    PerformanceConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.data import (
    ArrayData,
    ArrayDataWarning,
    StationData,
)
from shakelab.exploration.yarray.monitor import ArrayMonitor
from shakelab.exploration.yarray.operational import OperationalAssessment
from shakelab.exploration.yarray.performance import (
    PerformanceResult,
    PerformanceSweepResult,
    PlaneWaveSource,
)
from shakelab.exploration.yarray.response import (
    ArrayGeometry,
    ArrayResponse,
    ResponseAnalysis,
    ResponseGrid,
)

__all__ = [
    "ArrayData",
    "ArrayDataWarning",
    "ArrayGeometry",
    "ArrayMonitor",
    "ArrayResponse",
    "Beamformer",
    "BeamformingConfig",
    "BeamformingResult",
    "BeamformingWarning",
    "DispersionCurve",
    "EllipticityGrid",
    "FrequencyConfig",
    "FrequencySpectrum",
    "OperationalAssessment",
    "PerformanceConfig",
    "PerformanceResult",
    "PerformanceSweepResult",
    "PlaneWaveSource",
    "ResponseAnalysis",
    "ResponseGrid",
    "SearchGrid",
    "SpectralMaximum",
    "StationData",
    "WindowConfig",
]
