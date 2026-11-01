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
Regression tests for :meth:`shakelab.site.soil.Model1D.soil_response`.

The method builds a :class:`~shakelab.signals.spectrum.Spectrum` from
the transfer function of a layered model; it failed before the
spectrum constructor accepted spectral values.
"""

from __future__ import annotations

import numpy as np

from shakelab.signals.spectrum import Spectrum
from shakelab.site.soil import Model1D


def test_single_layer_resonance():
    """
    The SH resonance of a layer over a half-space is at Vs / (4 H).

    A 20 m layer with Vs = 300 m/s has its fundamental frequency at
    3.75 Hz (quarter-wavelength resonance).
    """
    model = Model1D()
    model.add_layer({"hl": 20.0, "vp": 600.0, "vs": 300.0, "dn": 1900.0,
                     "qp": 50.0, "qs": 25.0})
    model.add_layer({"hl": 0.0, "vp": 2000.0, "vs": 1000.0, "dn": 2200.0,
                     "qp": 100.0, "qs": 50.0})

    spectrum = model.soil_response(0.01, 1024)
    peak = spectrum.frequency[np.argmax(spectrum.amplitude)]

    assert isinstance(spectrum, Spectrum)
    assert spectrum.nsamp == 1024
    assert abs(peak - 3.75) < 0.15
