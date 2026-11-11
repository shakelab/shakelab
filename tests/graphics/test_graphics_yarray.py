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

from shakelab.exploration.yarray.beamforming import (  # noqa: E402
    Beamformer,
)
from shakelab.exploration.yarray.config import (  # noqa: E402
    BeamformingConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.data import ArrayData  # noqa: E402
from shakelab.exploration.yarray.operational import (  # noqa: E402
    OperationalConfig,
    assess_array,
)
from shakelab.exploration.yarray.performance import (  # noqa: E402
    EstimatorConfig,
    MonteCarloConfig,
    PerformanceConfig,
    PlaneWaveSource,
    evaluate_performance,
    performance_sweep,
)
from shakelab.exploration.yarray.response import (  # noqa: E402
    analyze_response,
    analyze_search_domain,
    array_response,
)
from shakelab.exploration.yarray.statistics import (  # noqa: E402
    maximum_density,
)
from shakelab.graphics.yarray import (  # noqa: E402
    plot_array_geometry,
    plot_array_response,
    plot_density,
    plot_dispersion,
    plot_ellipticity,
    plot_operational,
    plot_performance,
    plot_performance_sweep,
    plot_resolution,
    plot_search_domain,
    plot_spectrum,
    plot_velocity_limits,
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
        array = ArrayData(ARRAY)

    config = BeamformingConfig(
        frequencies=[3.0, 5.0],
        grid=SearchGrid(100.0, 2000.0, 30, azimuth_size=36),
        window=WindowConfig(length=10, windows_per_block=20),
        spectra=True, maxima=2, statistical_uncertainty="empirical",
    )
    beamformer = Beamformer(array, config)
    vertical = beamformer.process()
    beamformer.config.update(
        components="rayleigh_joint", statistical_uncertainty="none",
    )
    joint = beamformer.process()
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


@pytest.mark.parametrize("coordinate", ["log_slowness", "wavenumber"])
def test_plot_density(results, coordinate):
    """
    Relative density with peaks and summary.

    Logarithmic coordinates are shown as velocity on a log axis.
    """
    vertical, _ = results
    density = maximum_density(vertical, coordinate=coordinate)
    figure, axis = plot_density(density, peaks=density.peaks(2),
                                summary=density.summary(), scale="log")

    if coordinate == "log_slowness":
        assert axis.get_yscale() == "log"
        assert "velocity" in axis.get_ylabel()
    else:
        assert "wavenumber" in axis.get_ylabel()

    assert np.nanmax(axis.collections[0].get_array()) == pytest.approx(1.0)

    with pytest.raises(ValueError, match="scale"):
        plot_density(density, scale="sqrt")
    with pytest.raises(TypeError):
        plot_density(vertical)


# ----------------------------------------------------------------------------
# Response, performance and operational assessment
# ----------------------------------------------------------------------------


SEVEN = np.array([[0, 0], [30, 0], [0, 40], [25, 35], [-20, 10],
                  [10, -30], [-35, -20]], dtype=float)


@pytest.fixture(autouse=True)
def _close_figures():
    """Close the figures of each test."""
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def analysis():
    """Response analysis of a seven-station array."""
    yield analyze_response(array_response(SEVEN))
    plt.close("all")


@pytest.fixture(scope="module")
def performance(analysis):
    """Experiment, sweep and assessment on the seven-station array."""
    config = PerformanceConfig(
        snr_db=-3.0,
        grid=SearchGrid(0.0, 0.2, 21, quantity="wavenumber",
                        azimuth_size=36),
        estimator=EstimatorConfig("capon"),
        monte_carlo=MonteCarloConfig(20, 20, seed=1),
    )
    response = analysis.response
    source = PlaneWaveSource(0.15, 60.0)
    result = evaluate_performance(response, source, config)
    sweep = performance_sweep(response, source, config, [-5.0, 5.0],
                              [10, 20])
    assessment = assess_array(response, config, OperationalConfig(
        [3.0, 10.0], [200.0, 800.0], azimuths=[0.0, 90.0]
    ))
    yield result, sweep, assessment
    plt.close("all")


def test_plot_array_geometry(analysis):
    """Stations, centroid, labels and principal axes."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        array = ArrayData(ARRAY)

    _, axis = plot_array_geometry(array, labels=True)
    assert len(axis.collections[0].get_offsets()) == 9
    assert [text.get_text() for text in axis.texts] == list(array.sids)
    _, axis = plot_array_geometry(analysis.response, principal_axes=False,
                                  centered=True)
    np.testing.assert_allclose(axis.collections[0].get_offsets(),
                               analysis.response.geometry
                               .centered_coordinates)

    with pytest.raises(ValueError, match="ArrayData"):
        plot_array_geometry(analysis.response, labels=True)
    with pytest.raises(TypeError):
        plot_array_geometry(SEVEN)


def test_plot_array_response(analysis):
    """Response in decibels with lobe contours and peaks."""
    _, axis = plot_array_response(analysis, peaks=3)
    image = axis.collections[0].get_array()
    assert image.max() == pytest.approx(0.0) and image.min() >= -40.0
    assert len([line for line in axis.lines if line.get_marker() == "x"]) \
        == 3
    _, axis = plot_array_response(analysis.response, scale="linear")
    assert axis.collections[0].get_array().max() == pytest.approx(1.0)

    with pytest.raises(ValueError, match="scale"):
        plot_array_response(analysis, scale="sqrt")
    with pytest.raises(ValueError, match="floor"):
        plot_array_response(analysis, floor=3.0)


def test_plot_resolution_and_limits(analysis):
    """Directional radii, velocity limits and search domain."""
    _, axis = plot_resolution(analysis)
    np.testing.assert_allclose(axis.lines[0].get_ydata()[:-1],
                               analysis.main_lobe.radii)
    _, axis = plot_resolution(analysis.main_lobe, quantity="separation",
                              polar=True)
    assert axis.name == "polar"
    limits = analysis.velocity_limits([1.0, 2.0, 5.0])
    _, axis = plot_velocity_limits(limits)
    np.testing.assert_allclose(axis.lines[0].get_ydata(), limits.minimum)
    domain = analyze_search_domain(analysis.response, 0.03)
    _, axis = plot_search_domain(analysis, domain)
    assert "safe" in axis.get_title()

    with pytest.raises(ValueError, match="quantity"):
        plot_resolution(analysis, quantity="area")
    with pytest.raises(ValueError, match="polar"):
        plot_resolution(analysis, polar=True, ax=axis)


def test_plot_performance(performance):
    """Estimates by outcome, sweep and envelope."""
    result, sweep, assessment = performance
    from matplotlib.collections import PathCollection

    _, axis = plot_performance(result)
    points = np.concatenate([
        item.get_offsets() for item in axis.collections
        if isinstance(item, PathCollection)
    ])
    assert len(points) == result.samples.size
    _, axis = plot_performance_sweep(sweep, "alias_probability")
    assert axis.collections[0].get_array().shape == (2, 2)
    _, axis = plot_operational(assessment)
    np.testing.assert_allclose(
        axis.collections[0].get_array(),
        assessment.conservative_local_probability.T,
    )

    with pytest.raises(ValueError, match="quantity"):
        plot_performance_sweep(sweep, "power")
    with pytest.raises(ValueError, match="quantity"):
        plot_operational(assessment, "power")
    with pytest.raises(TypeError):
        plot_performance(sweep)
