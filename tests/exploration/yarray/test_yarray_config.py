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
Unit tests for :mod:`shakelab.exploration.yarray.config`.

The configuration classes are checked for validation, conversion of
strings to enumerations, incompatible options, JSON round trips and
the reading of the configuration files of the YArray examples of
August 2026.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    BeamformingMethod,
    ComponentMode,
    EllipticityGrid,
    FrequencyConfig,
    GridQuantity,
    GridScale,
    SearchGrid,
    TransformMethod,
    WindowConfig,
    WindowType,
)

# Configuration of example 10 of the YArray suite of August 2026.
EXAMPLE = {
    "version": 1,
    "frequencies": {"minimum": 2.0, "maximum": 12.0, "size": 100,
                    "scale": "log"},
    "spectral": {"transform": "direct", "fft_selection": "nearest",
                 "fft_padding": None},
    "method": {"name": "capon", "diagonal_loading": 0.001,
               "music_sources": 1},
    "components": "vertical",
    "grid": {"quantity": "velocity", "scale": "log", "minimum": 100.0,
             "maximum": 2000.0, "size": 100, "azimuth_size": 180},
    "window": {"type": "fixed", "length": 5.0, "taper": 0.1,
               "overlap": 0.5, "windows_per_block": 40,
               "block_overlap": 0.5},
    "ellipticity": {"minimum_angle": -85.0, "maximum_angle": 85.0,
                    "size": 61},
    "output": {"mode": "both", "maxima": 3, "neighborhood": [3, 3],
               "minimum_relative_power": 0.0},
    "uncertainty": {"statistical": "none", "geometric": "none"},
    "strict": True,
    "graphics": {"dispersion": True},
}

GRID = SearchGrid(100.0, 2000.0, 50)


# ----------------------------------------------------------------------------
# Sections
# ----------------------------------------------------------------------------


def test_frequency_config():
    """Explicit values and log or linear ranges."""
    np.testing.assert_allclose(
        FrequencyConfig(minimum=1.0, maximum=100.0, size=3).targets(),
        [1.0, 10.0, 100.0],
    )
    np.testing.assert_allclose(
        FrequencyConfig(minimum=1.0, maximum=3.0, size=3,
                        scale="linear").targets(),
        [1.0, 2.0, 3.0],
    )
    assert FrequencyConfig(values=[2.0, 3.0]).is_explicit

    with pytest.raises(ValueError, match="FFT bin"):
        FrequencyConfig(minimum=1.0, maximum=3.0).targets()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"values": [2.0, 1.0]},
        {"values": [0.0, 1.0]},
        {"values": [1.0], "minimum": 1.0},
        {"minimum": 1.0},
        {"minimum": 3.0, "maximum": 1.0},
        {"minimum": 1.0, "maximum": 3.0, "size": 1},
        {"minimum": 1.0, "maximum": 3.0, "scale": "cubic"},
    ],
)
def test_invalid_frequency_config(kwargs):
    """Mixed, decreasing or incomplete definitions are rejected."""
    with pytest.raises(ValueError):
        FrequencyConfig(**kwargs)


def test_window_config():
    """Strings become enumerations; sizes in samples."""
    window = WindowConfig(type="fixed", length=5.0)

    assert window.type is WindowType.FIXED
    assert window.size(0.01) == 500
    assert WindowConfig(length=20).size(0.01, 4.0) == 500

    with pytest.raises(ValueError, match="frequency"):
        WindowConfig().size(0.01)

    for kwargs in ({"overlap": 1.0}, {"taper": 1.5}, {"length": 0.0},
                   {"windows_per_block": 0}):
        with pytest.raises(ValueError):
            WindowConfig(**kwargs)


def test_search_grid():
    """Radial values, wavenumbers and wavevectors."""
    grid = SearchGrid(100.0, 400.0, 3, azimuth_size=4)

    np.testing.assert_allclose(grid.values, [100.0, 200.0, 400.0])
    np.testing.assert_allclose(grid.azimuth_values, [0, 90, 180, 270])
    np.testing.assert_allclose(
        grid.wavenumbers(2.0), 4 * np.pi / grid.values
    )
    kx, ky = grid.wavevectors(2.0)
    assert kx.shape == (4, 3)
    np.testing.assert_allclose(kx[1], grid.wavenumbers(2.0))
    np.testing.assert_allclose(ky[1], 0.0, atol=1e-15)

    slowness = SearchGrid(0.001, 0.01, 2, quantity="slowness")
    np.testing.assert_allclose(
        slowness.wavenumbers(1.0), 2 * np.pi * np.array([0.001, 0.01])
    )


def test_search_grid_explicit_azimuths():
    """Explicit azimuths are reduced to [0, 360) and kept in JSON."""
    grid = SearchGrid(0.0, 0.1, 5, quantity="wavenumber", scale="linear",
                      azimuths=[10.0, 370.5])

    np.testing.assert_allclose(grid.azimuth_values, [10.0, 10.5])
    assert grid.azimuth_size == 2
    assert grid.to_dict()["azimuths"] == [10.0, 10.5]

    with pytest.raises(ValueError):
        SearchGrid(0.0, 0.1, 5, quantity="wavenumber", scale="log")
    with pytest.raises(ValueError):
        SearchGrid(100.0, 2000.0, 5, azimuths=[10.0, 370.0])
    with pytest.raises(ValueError):
        SearchGrid(0.0, 2000.0, 5)


def test_search_grid_default_scale():
    """Default spacing: log for velocities, linear for wavenumbers."""
    assert SearchGrid(100.0, 2000.0, 5).scale is GridScale.LOG
    assert SearchGrid(0.001, 0.01, 5, quantity="slowness").scale is \
        GridScale.LOG
    grid = SearchGrid(0.0, 0.1, 5, quantity="wavenumber")
    assert grid.scale is GridScale.LINEAR
    assert grid.to_dict()["scale"] == "linear"


def test_ellipticity_grid():
    """Angles strictly inside (-90, 90); values are their tangents."""
    grid = EllipticityGrid(-45.0, 45.0, 3)

    np.testing.assert_allclose(grid.values, [-1.0, 0.0, 1.0], atol=1e-15)

    with pytest.raises(ValueError):
        EllipticityGrid(-90.0, 45.0)


# ----------------------------------------------------------------------------
# BeamformingConfig
# ----------------------------------------------------------------------------


def test_strings_and_defaults():
    """Strings are converted; frequencies become a FrequencyConfig."""
    config = BeamformingConfig([2.0, 4.0], GRID, method="Bartlett",
                               components="radial")

    assert config.method is BeamformingMethod.BARTLETT
    assert config.components is ComponentMode.RADIAL
    assert config.frequencies.is_explicit
    assert config.window.windows_per_block == 40
    assert config.peak_neighborhood == (3, 3)

    with pytest.raises(ValueError, match="'capon'"):
        BeamformingConfig([2.0], GRID, method="beam")


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"transform": "rfft"}, "fixed duration"),
        ({"components": "rayleigh_joint", "method": "bartlett"}, "Capon"),
        ({"statistical_uncertainty": "empirical", "method": "music"},
         "Bartlett and Capon"),
        ({"geometric_uncertainty": "expected_coherence"}, "Bartlett"),
        ({"diagonal_loading": -1.0}, "negative"),
        ({"minimum_relative_power": 2.0}, r"\[0, 1\]"),
        ({"peak_neighborhood": (3,)}, "two"),
    ],
)
def test_incompatible_options(kwargs, match):
    """Unsupported combinations are rejected with a clear message."""
    with pytest.raises(ValueError, match=match):
        BeamformingConfig([2.0], GRID, **kwargs)


def test_direct_transform_needs_sampled_frequencies():
    """A range without size (FFT bins) needs the RFFT."""
    frequencies = FrequencyConfig(minimum=1.0, maximum=5.0)

    with pytest.raises(ValueError, match="direct"):
        BeamformingConfig(frequencies, GRID)

    config = BeamformingConfig(
        frequencies, GRID, transform="rfft",
        window=WindowConfig(type="fixed", length=4.0),
    )
    assert config.transform is TransformMethod.RFFT


def test_wrong_section_types():
    """Sections must be the configuration classes."""
    with pytest.raises(TypeError):
        BeamformingConfig([2.0], {"minimum": 100})
    with pytest.raises(TypeError):
        BeamformingConfig([2.0], GRID, strict="yes")


def test_replace():
    """replace returns a validated copy."""
    config = BeamformingConfig([2.0], GRID)
    music = config.replace(method="music", music_sources=2)

    assert music.method is BeamformingMethod.MUSIC
    assert config.method is BeamformingMethod.CAPON

    window = config.window.replace(length=10.0)
    assert window.length == 10.0

    with pytest.raises(ValueError):
        config.replace(components="rayleigh_joint", method="bartlett")


def test_reads_august_example():
    """
    The configurations of the August examples are read unchanged.

    Unknown sections (``graphics``) are ignored.
    """
    config = BeamformingConfig.from_dict(EXAMPLE)

    assert config.window.type is WindowType.FIXED
    assert config.grid.quantity is GridQuantity.VELOCITY
    assert config.maxima == 3
    assert config.frequencies.size == 100

    expected = {key: value for key, value in EXAMPLE.items()
                if key != "graphics"}
    assert config.to_dict() == expected


def test_minimal_json_uses_defaults():
    """Only frequencies and grid are required."""
    config = BeamformingConfig.from_dict({
        "frequencies": {"values": [2.0, 3.0]},
        "grid": {"minimum": 100, "maximum": 1000, "size": 10},
        "method": "bartlett",
    })

    assert config.method is BeamformingMethod.BARTLETT
    assert config.window.length == 20.0

    with pytest.raises(ValueError, match="version"):
        BeamformingConfig.from_dict({**EXAMPLE, "version": 2})
    with pytest.raises(ValueError, match="frequencies and grid"):
        BeamformingConfig.from_dict({"grid": EXAMPLE["grid"]})
    with pytest.raises(TypeError):
        BeamformingConfig.from_dict({**EXAMPLE, "window": [1, 2]})


def test_read_and_write(tmp_path):
    """Files round trip; existing files are protected."""
    path = tmp_path / "config.json"
    config = BeamformingConfig(
        FrequencyConfig(minimum=1.0, maximum=8.0, size=12),
        SearchGrid(0.0, 0.05, 40, quantity="wavenumber", scale="linear"),
        components="rayleigh_joint",
        window=WindowConfig(length=15, windows_per_block=60),
    )
    config.write(path)

    assert BeamformingConfig.read(path).to_dict() == config.to_dict()

    with pytest.raises(FileExistsError):
        config.write(path)

    config.write(path, owrite=True)
    path.write_text("{ invalid", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid JSON"):
        BeamformingConfig.read(path)

    assert json.loads(json.dumps(config.to_dict())) == config.to_dict()
