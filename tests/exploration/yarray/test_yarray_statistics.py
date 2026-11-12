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
Unit tests for :mod:`shakelab.exploration.yarray.statistics`.

Densities are checked on results with known maxima: normalization,
kernel widths (grid resolution, Silverman's rule, fixed), mode and
quantiles of a single Gaussian kernel, weights, and the peaks and
prominences of two modes against ``scipy.signal.peak_prominences``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import norm

from shakelab.exploration.yarray.beamforming import (
    BeamformingResult,
    SpectralMaximum,
)
from shakelab.exploration.yarray.config import BeamformingConfig, SearchGrid
from shakelab.exploration.yarray.statistics import (
    DensityCoordinate,
    MaximumDensity,
    density_peaks,
    maximum_density,
)
from shakelab.libutils.timeN import Date

GRID = SearchGrid(100.0, 1000.0, 200, scale="linear")
CONFIG = BeamformingConfig(frequencies=[4.0], grid=GRID)
START = Date("2026-01-01T00:00:00")


def _maximum(velocity, frequency=4.0, block=0, rank=1, power=1.0,
             component="vertical", angle=None):
    """Spectral maximum of a velocity (and ellipticity angle)."""
    return SpectralMaximum(
        block=block, time=START, frequency=frequency, component=component,
        rank=rank, wavenumber=2 * math.pi * frequency / velocity,
        azimuth=0.0, power=power,
        ellipticity=None if angle is None else math.tan(math.radians(angle)),
        ellipticity_angle=angle,
    )


def _rayleigh(velocity, angle, **kwargs):
    """Maximum of the Rayleigh joint processing."""
    return _maximum(velocity, component="rayleigh_joint", angle=angle,
                    **kwargs)


def _result(maxima, config=CONFIG):
    """Result holding the given maxima."""
    frequencies = sorted({item.frequency for item in maxima})
    return BeamformingResult._from_parts(
        maxima=maxima, spectra=None, config=config,
        block_times={f: (START,) for f in frequencies}, sids=(),
        coordinates=np.empty((0, 3)), units=None,
    )


def test_single_maximum_is_one_gaussian_kernel():
    """
    One maximum gives a Gaussian of the grid spacing at its value.

    In log-velocity the spacing of a linear velocity grid of step
    dv is dv / v; mode and median are at log(v), the 16% and 84%
    quantiles at log(v) -+ 0.994 sigma.
    """
    density = maximum_density(_result([_maximum(400.0)]),
                              coordinate="log_velocity", size=2001)
    sigma = (900.0 / 199.0) / 400.0
    summary = density.summary()

    assert density.bandwidths[0] == pytest.approx(sigma, rel=0.01)
    assert np.trapezoid(density.density[0], density.values) == (
        pytest.approx(1.0)
    )
    assert summary.mode[0] == pytest.approx(math.log(400.0), abs=1e-3)
    assert summary.median[0] == pytest.approx(math.log(400.0), abs=1e-4)
    assert summary.upper[0] - summary.median[0] == pytest.approx(
        norm.ppf(0.84) * sigma, rel=0.01
    )
    np.testing.assert_allclose(density.velocity, np.exp(density.values))


def test_coordinates_and_counts():
    """Wavenumber and log-slowness grids; counts by frequency."""
    maxima = [_maximum(400.0), _maximum(410.0), _maximum(300.0, 6.0)]
    result = _result(
        maxima, BeamformingConfig(frequencies=[4.0, 6.0], grid=GRID)
    )

    density = maximum_density(result)
    assert density.coordinate is DensityCoordinate.WAVENUMBER
    np.testing.assert_array_equal(density.frequencies, [4.0, 6.0])
    np.testing.assert_array_equal(density.counts, [2, 1])

    slowness = maximum_density(result, coordinate="log_slowness")
    assert slowness.summary().mode[1] == pytest.approx(
        math.log(1.0 / 300.0), abs=0.01
    )

    with pytest.raises(ValueError, match="depends on f"):
        density.velocity


def test_bandwidth_modes():
    """Fixed width; Silverman's rule (Silverman, 1986, eq. 3.31)."""
    values = np.array([380.0, 390.0, 400.0, 405.0, 430.0])
    result = _result([_maximum(v) for v in values])

    fixed = maximum_density(result, coordinate="log_velocity",
                            bandwidth_mode="fixed", bandwidth=0.05)
    assert fixed.bandwidths[0] == 0.05

    adaptive = maximum_density(result, coordinate="log_velocity",
                               bandwidth_mode="adaptive")
    # 0.9 min(standard deviation, IQR / 1.349) n**-0.2, with the
    # quartiles of the weighted empirical distribution function.
    logs = np.log(values)
    q25, q75 = np.interp([0.25, 0.75], np.arange(1, 6) / 5, logs)
    scale = min(np.std(logs), (q75 - q25) / 1.349)
    assert adaptive.bandwidths[0] == pytest.approx(
        0.9 * scale * len(values) ** -0.2, rel=1e-12
    )

    wider = maximum_density(result, bandwidth_factor=2.0)
    assert wider.bandwidths[0] == pytest.approx(
        2 * maximum_density(result).bandwidths[0]
    )


def test_weights_and_ranks():
    """Relative power weights within each block; rank selection."""
    maxima = [
        _maximum(400.0, block=0, rank=1, power=10.0),
        _maximum(700.0, block=0, rank=2, power=1.0),
        _maximum(400.0, block=1, rank=1, power=0.1),
        _maximum(700.0, block=1, rank=2, power=0.1),
    ]
    result = _result(maxima)

    uniform = maximum_density(result, coordinate="log_velocity").peaks(2)
    relative = maximum_density(result, coordinate="log_velocity",
                               weighting="relative_power").peaks(2)
    first = maximum_density(result, coordinate="log_velocity", ranks=[1])

    # Uniform: two equal modes; relative power: (1 + 1) against
    # (0.1 + 1), so the second mode has 0.55 of the first.
    assert uniform.relative_density[0, 1] == pytest.approx(1.0, rel=1e-3)
    assert relative.relative_density[0, 1] == pytest.approx(0.55, rel=0.01)
    assert first.counts[0] == 2


def test_peaks_and_prominence():
    """
    Two modes: heights and prominences as scipy.signal.

    The weaker mode sits on the flank of the stronger one: its
    prominence is measured from the saddle between them (smaller than
    its height), not from the lowest density on its side (as in the
    August 2026 implementation, which gave the height).
    """
    maxima = [_maximum(400.0)] * 4 + [_maximum(330.0)]
    density = maximum_density(_result(maxima), coordinate="log_velocity",
                              bandwidth_mode="fixed", bandwidth=0.06,
                              size=512)
    peaks = density.peaks(number=3)

    assert peaks.count[0] == 2
    np.testing.assert_allclose(np.exp(peaks.values[0, :2]), [400.0, 330.0],
                               rtol=0.02)
    assert peaks.prominence[0, 0] == pytest.approx(1.0)
    assert peaks.prominence[0, 1] < 0.8 * peaks.relative_density[0, 1]

    strong = density_peaks(density, minimum_prominence=0.2)
    assert strong.count[0] == 1
    tall = density_peaks(density, minimum_relative_density=0.9)
    assert tall.count[0] == 1


def test_relative_density_and_info(capsys):
    """Per-frequency normalization with a floor; summary text."""
    density = maximum_density(_result([_maximum(400.0)]))

    relative = density.relative_density(floor=0.01)
    assert relative.max() == 1.0 and relative.min() == 0.01

    density.info()
    assert "1 maxima" in capsys.readouterr().out

    with pytest.raises(ValueError):
        density.relative_density(floor=1.0)


def test_ellipticity_angle_density():
    """
    Density of the ellipticity angle: kernels of the step of the
    ellipticity grid (170 / 60 degrees by default), grid from -90 to
    90 degrees, signed ellipticity of the grid.
    """
    maxima = [_rayleigh(300.0, -30.0, block=b) for b in range(3)]
    maxima.append(_rayleigh(300.0, 40.0, block=3))
    density = maximum_density(_result(maxima), coordinate="ellipticity_angle",
                              size=721)

    assert density.coordinate is DensityCoordinate.ELLIPTICITY_ANGLE
    assert density.values[[0, -1]].tolist() == [-90.0, 90.0]
    assert density.bandwidths[0] == pytest.approx(170.0 / 60.0)
    assert np.trapezoid(density.density[0], density.values) == (
        pytest.approx(1.0)
    )
    assert density.summary().mode[0] == pytest.approx(-30.0, abs=0.3)
    peaks = density.peaks(2)
    np.testing.assert_allclose(peaks.values[0], [-30.0, 40.0], atol=0.3)
    np.testing.assert_allclose(density.ellipticity,
                               np.tan(np.radians(density.values)))

    with pytest.raises(ValueError, match="ellipticity angle"):
        density.velocity
    with pytest.raises(ValueError, match="not an ellipticity"):
        maximum_density(_result(maxima)).ellipticity


def test_ellipticity_angle_bandwidth_of_exact_angles():
    """
    The angular processing finds the angle without a grid: the kernels
    take the step of the default ellipticity grid, whatever the
    configuration holds.
    """
    config = BeamformingConfig(
        frequencies=[4.0], grid=GRID, components="rayleigh_angular",
        ellipticity_angle_size=11,
    )
    maxima = [_maximum(300.0, component="rayleigh_angular", angle=-20.0)]
    density = maximum_density(_result(maxima, config),
                              coordinate="ellipticity_angle")

    assert density.bandwidths[0] == pytest.approx(170.0 / 60.0)
    joint = maximum_density(
        _result([_rayleigh(300.0, -20.0)],
                CONFIG.copy(components="rayleigh_joint",
                            ellipticity_angle_size=11)),
        coordinate="ellipticity_angle",
    )
    assert joint.bandwidths[0] == pytest.approx(17.0)


def test_ellipticity_angle_is_periodic():
    """-90 and 90 degrees are the same motion: the kernels wrap."""
    density = maximum_density(_result([_rayleigh(300.0, 88.0)]),
                              coordinate="ellipticity_angle", size=361)
    row = density.density[0]

    assert row[0] == pytest.approx(row[-1])
    # The mass beyond 90 degrees reappears above -90 degrees.
    assert row[density.values < -85.0].max() > 0.1 * row.max()
    assert np.trapezoid(row, density.values) == pytest.approx(1.0, abs=1e-3)


def test_ellipticity_angle_needs_rayleigh_maxima():
    """Maxima without ellipticity cannot give an angle density."""
    with pytest.raises(ValueError, match="rayleigh_joint"):
        maximum_density(_result([_maximum(400.0)]),
                        coordinate="ellipticity_angle")


def test_velocity_selection(capsys):
    """Maxima selected by phase velocity, e.g. one mode."""
    maxima = [_rayleigh(200.0, -40.0), _rayleigh(210.0, -42.0, block=1),
              _rayleigh(600.0, 20.0, block=2)]
    result = _result(maxima)

    slow = maximum_density(result, coordinate="ellipticity_angle",
                           velocity=(None, 300.0))
    assert slow.counts.tolist() == [2]
    assert slow.velocity_range == (None, 300.0)
    assert slow.summary().mode[0] < 0.0
    fast = result.density(coordinate="ellipticity_angle",
                          velocity=(500, None))
    assert fast.counts.tolist() == [1]
    assert fast.summary().mode[0] == pytest.approx(20.0, abs=1.0)
    assert maximum_density(result).velocity_range is None
    fast.info()
    assert "velocity 500 to - m/s" in capsys.readouterr().out

    with pytest.raises(ValueError, match="no maximum"):
        maximum_density(result, velocity=(300.0, 500.0))
    with pytest.raises(ValueError, match="below"):
        maximum_density(result, velocity=(300.0, 200.0))
    with pytest.raises(ValueError, match="positive"):
        maximum_density(result, velocity=(0.0, 200.0))
    with pytest.raises(TypeError, match="pair"):
        maximum_density(result, velocity=300.0)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"coordinate": "velocity"}, "coordinate"),
        ({"bandwidth_mode": "fixed"}, "bandwidth"),
        ({"bandwidth": 0.1}, "bandwidth"),
        ({"size": 8}, "size"),
        ({"ranks": [3]}, "no maximum"),
        ({"component": "radial"}, "component"),
        ({"bandwidth_factor": 0.0}, "bandwidth_factor"),
    ],
)
def test_invalid_options(kwargs, match):
    """Invalid options and empty selections are rejected."""
    with pytest.raises(ValueError, match=match):
        maximum_density(_result([_maximum(400.0)]), **kwargs)


def test_ambiguous_component_and_types():
    """Several components need a choice; wrong inputs fail."""
    result = _result([_maximum(400.0), _maximum(400.0, component="radial")])

    with pytest.raises(ValueError, match="choose"):
        maximum_density(result)

    assert maximum_density(result, component="radial").component == (
        "radial"
    )

    with pytest.raises(TypeError):
        maximum_density(result.maxima)
    with pytest.raises(TypeError):
        density_peaks("density")
    with pytest.raises(ValueError):
        maximum_density(_result([_maximum(400.0)])).summary((0.6, 0.9))

    assert isinstance(maximum_density(result, component="vertical"),
                      MaximumDensity)
