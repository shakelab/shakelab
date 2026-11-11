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


def _config(**kwargs):
    """Configuration with explicit frequencies and the test grid."""
    return BeamformingConfig(frequencies=[2.0], grid=GRID, **kwargs)


def test_strings_and_defaults():
    """Strings are converted; frequencies become a FrequencyConfig."""
    config = BeamformingConfig(frequencies=[2.0, 4.0], grid=GRID,
                               method="Bartlett", components="radial")

    assert config.method is BeamformingMethod.BARTLETT
    assert config.components is ComponentMode.RADIAL
    assert config.frequencies.is_explicit
    assert config.window.windows_per_block == 40
    assert config.peak_neighborhood == (3, 3)
    assert config.complete

    with pytest.raises(ValueError, match="'capon'"):
        _config(method="beam")


def test_single_parameters_equal_objects():
    """Single parameters build the same objects as the classes."""
    single = BeamformingConfig(
        fmin=3.0, fmax=10.0, nfreq=12, grid_min=100.0, grid_max=2000.0,
        grid_size=50, window_type="fixed", window_length=2.0,
        windows_per_block=30, ellipticity_angle_size=31,
    )
    objects = BeamformingConfig(
        frequencies=FrequencyConfig(minimum=3.0, maximum=10.0, size=12),
        grid=GRID,
        window={"type": "fixed", "length": 2.0, "windows_per_block": 30},
        ellipticity_grid=EllipticityGrid(size=31),
    )

    assert single == objects
    assert single.fmin == 3.0 and single.window_length == 2.0
    assert single.grid_size == 50 and single.azimuth_size == 180
    assert BeamformingConfig().fmin is None


def test_single_parameters_apply_on_objects():
    """Single parameters replace fields of the objects and the file."""
    config = BeamformingConfig(
        window=WindowConfig(type="fixed", length=5.0), window_length=4.0
    )
    assert config.window.to_dict() == WindowConfig(
        type="fixed", length=4.0
    ).to_dict()

    config = BeamformingConfig(EXAMPLE, method="music", nfreq=20)
    assert config.method is BeamformingMethod.MUSIC
    assert config.frequencies.size == 20
    assert config.fmin == 2.0


def test_update_and_attributes():
    """update and assignment validate all the values together."""
    config = _config()
    config.update(method="music", music_sources=2)
    config.grid_max = 3000.0

    assert config.method is BeamformingMethod.MUSIC
    assert config.grid.maximum == 3000.0

    with pytest.raises(ValueError):
        config.update(components="rayleigh_joint", method="bartlett")

    assert config.components is ComponentMode.VERTICAL

    config.update(components="rayleigh_joint", method="capon")
    assert config.components is ComponentMode.RAYLEIGH_JOINT

    with pytest.raises(TypeError, match="unknown"):
        config.update(grid_maximum=1.0)
    with pytest.raises(AttributeError):
        config.speed = 1.0
    with pytest.raises(AttributeError):
        config.window.length = 1.0


def test_frequency_and_grid_changes():
    """Ranges replace values; derived defaults follow the new values."""
    config = _config()
    config.update(fmin=1.0, fmax=4.0, nfreq=4)
    assert not config.frequencies.is_explicit

    config.frequencies = [3.0, 5.0]
    assert config.frequencies.is_explicit

    config.update(grid_quantity="wavenumber", grid_min=0.0, grid_max=0.1)
    assert config.grid.scale is GridScale.LINEAR

    config.azimuths = [0.0, 90.0]
    assert config.azimuth_size == 2

    config.azimuth_size = 36
    assert config.azimuths is None and config.azimuth_size == 36


def test_copy_and_equality():
    """Copies are independent; equality compares the values."""
    config = _config()
    copy = config.copy()
    copy.method = "bartlett"

    assert config.method is BeamformingMethod.CAPON
    assert copy != config
    copy.method = "capon"
    assert copy == config


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
        _config(**kwargs)


def test_direct_transform_needs_sampled_frequencies():
    """A range without size (FFT bins) needs the RFFT."""
    frequencies = FrequencyConfig(minimum=1.0, maximum=5.0)

    with pytest.raises(ValueError, match="direct"):
        BeamformingConfig(frequencies=frequencies, grid=GRID)

    config = BeamformingConfig(
        frequencies=frequencies, grid=GRID, transform="rfft",
        window_type="fixed", window_length=4.0,
    )
    assert config.transform is TransformMethod.RFFT


def test_wrong_types():
    """Objects of the wrong type and non-boolean flags are rejected."""
    with pytest.raises(TypeError):
        BeamformingConfig(grid=[100, 2000, 10])
    with pytest.raises(TypeError):
        _config(strict="yes")


def test_reads_august_example():
    """
    The configurations of the August examples are read unchanged.

    Unknown sections (``graphics``) are ignored.
    """
    config = BeamformingConfig(EXAMPLE)

    assert config.window.type is WindowType.FIXED
    assert config.grid.quantity is GridQuantity.VELOCITY
    assert config.maxima == 3
    assert config.frequencies.size == 100

    assert config.spectra is True
    expected = {key: value for key, value in EXAMPLE.items()
                if key != "graphics"}
    expected["output"] = {"spectra": True, "maxima": 3,
                          "neighborhood": [3, 3],
                          "minimum_relative_power": 0.0}
    assert config.to_dict() == expected

    for mode, spectra in (("maxima", False), ("spectra", True)):
        data = {**EXAMPLE, "output": {"mode": mode}}
        assert BeamformingConfig(data).spectra is spectra

    with pytest.raises(TypeError, match="boolean"):
        BeamformingConfig(spectra="yes")


def test_minimal_json_uses_defaults():
    """Missing fields take the defaults; frequencies and grid are needed."""
    config = BeamformingConfig({
        "frequencies": {"values": [2.0, 3.0]},
        "grid": {"minimum": 100, "maximum": 1000, "size": 10},
        "method": "bartlett",
    })

    assert config.method is BeamformingMethod.BARTLETT
    assert config.window.length == 20.0

    with pytest.raises(ValueError, match="version"):
        BeamformingConfig({**EXAMPLE, "version": 2})
    with pytest.raises(TypeError):
        BeamformingConfig({**EXAMPLE, "window": [1, 2]})

    partial = BeamformingConfig({"grid": EXAMPLE["grid"]})
    assert not partial.complete
    assert "frequencies" not in partial.to_dict()

    with pytest.raises(ValueError, match="needs frequencies"):
        partial.check_complete()


def test_read_and_write(tmp_path, capsys):
    """Files round trip; existing files are protected."""
    path = tmp_path / "config.json"
    config = BeamformingConfig(
        frequencies=FrequencyConfig(minimum=1.0, maximum=8.0, size=12),
        grid=SearchGrid(0.0, 0.05, 40, quantity="wavenumber"),
        components="rayleigh_joint",
        window=WindowConfig(length=15, windows_per_block=60),
    )
    config.write(path)

    assert BeamformingConfig(path) == config

    other = _config(method="music")
    other.read(path)
    assert other == config

    with pytest.raises(FileExistsError):
        config.write(path)

    config.write(path, owrite=True)
    path.write_text("{ invalid", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid JSON"):
        BeamformingConfig(path)

    assert json.loads(json.dumps(config.to_dict())) == config.to_dict()

    config.info()
    assert '"rayleigh_joint"' in capsys.readouterr().out
