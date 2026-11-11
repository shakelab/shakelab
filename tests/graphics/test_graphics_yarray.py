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
Unit tests for :mod:`shakelab.graphics.yarray`.

The figures are drawn with the non-interactive backend on results of
the SESAME excerpt; the tests check the plotted data and the errors,
not the appearance. Skipped without matplotlib.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from shakelab.exploration.yarray.beamforming import beamform  # noqa: E402
from shakelab.exploration.yarray.config import (  # noqa: E402
    BeamformingConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.io import read_array  # noqa: E402
from shakelab.graphics.yarray import (  # noqa: E402
    plot_dispersion,
    plot_ellipticity,
    plot_spectrum,
)

ARRAY = (
    Path(__file__).parents[1] / "exploration" / "yarray" / "data"
    / "sesame_array.json"
)


@pytest.fixture(scope="module")
def results():
    """Vertical and Rayleigh joint results with spectra."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        array = read_array(ARRAY)

    config = BeamformingConfig(
        [3.0, 5.0], SearchGrid(100.0, 2000.0, 30, azimuth_size=36),
        window=WindowConfig(length=10, windows_per_block=12),
        output="both", maxima=2, statistical_uncertainty="empirical",
    )
    vertical = beamform(array, config)
    joint = beamform(array, config.replace(
        components="rayleigh_joint", statistical_uncertainty="none",
    ))
    yield vertical, joint
    plt.close("all")


def test_plot_dispersion(results):
    """One series per rank, velocities of the maxima."""
    vertical, _ = results
    figure, axis = plot_dispersion(vertical)
    points = np.concatenate([item.get_offsets() for item in axis.collections])

    assert len(axis.collections) == 2
    assert len(points) == len(vertical.maxima)
    assert axis.get_yscale() == "log"

    _, axis = plot_dispersion(vertical, ranks=[1], log_velocity=False)
    assert len(axis.collections) == 1 and axis.get_yscale() == "linear"


def test_plot_ellipticity(results):
    """Ellipticity or its angle of the joint processing."""
    vertical, joint = results
    _, axis = plot_ellipticity(joint, ranks=[1], angle=True)
    values = axis.collections[0].get_offsets()[:, 1]

    np.testing.assert_allclose(
        sorted(values),
        sorted(item.ellipticity_angle for item in joint.maxima_at(rank=1)),
    )

    with pytest.raises(ValueError, match="ellipticity"):
        plot_ellipticity(vertical)


@pytest.mark.parametrize("view", ["kxky", "azimuth"])
@pytest.mark.parametrize("quantity", ["power", "standard_error"])
def test_plot_spectrum(results, view, quantity):
    """Spectra in both views, normalized power and standard error."""
    vertical, _ = results
    figure, axis = plot_spectrum(vertical, 5.0, block=0, view=view,
                                 quantity=quantity)
    image = axis.collections[0].get_array()

    if quantity == "power":
        assert np.nanmax(image) == pytest.approx(1.0)

    assert len(figure.axes) == 2


def test_plot_spectrum_ellipticity_and_units(results):
    """Ellipticity of the joint processing; power in units."""
    vertical, joint = results
    figure, _ = plot_spectrum(joint, 3.0, quantity="ellipticity")
    assert figure.axes[1].get_ylabel() == "Signed ellipticity"

    figure, _ = plot_spectrum(vertical, 3.0, normalize=False)
    assert "(unit)**2/Hz" in figure.axes[1].get_ylabel()


def test_plot_errors(results):
    """Invalid selections are rejected."""
    vertical, joint = results

    with pytest.raises(ValueError, match="block"):
        plot_spectrum(vertical, 5.0, block=99)
    with pytest.raises(ValueError, match="view"):
        plot_spectrum(vertical, 5.0, view="polar")
    with pytest.raises(ValueError, match="standard error"):
        plot_spectrum(joint, 5.0, quantity="standard_error")
    with pytest.raises(ValueError, match="not processed"):
        plot_dispersion(vertical, component="radial")
    with pytest.raises(TypeError):
        plot_dispersion(vertical.maxima)
