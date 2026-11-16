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
Unit tests for :mod:`shakelab.exploration.yarray.fdd`.

The reference is a four-storey shear frame excited by white noise at
every floor, simulated exactly by modal superposition (each mode a
single-degree-of-freedom oscillator discretized with a zero-order
hold): the peaks of the first singular value are its natural
frequencies and the modes its mode shapes.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.signal import cont2discrete, lfilter

from shakelab.exploration.yarray.beamforming import BeamformingResult
from shakelab.exploration.yarray.config import (
    FDDComponents,
    FDDConfig,
    FDDNormalization,
)
from shakelab.exploration.yarray.data import ArrayData, StationData
from shakelab.exploration.yarray.fdd import (
    FDD,
    FDDResult,
    FDDWarning,
    Mode,
    mac,
)
from shakelab.signals.base import Record

DELTA = 0.02
START = "2026-01-01T00:00:00"
STOREYS = 4
DAMPING = 0.02


def _shear_frame(first=1.5):
    """Frequencies (Hz) and mass-normalized shapes of a uniform frame."""
    stiffness = (
        2.0 * np.eye(STOREYS) - np.eye(STOREYS, k=1) - np.eye(STOREYS, k=-1)
    )
    stiffness[-1, -1] = 1.0
    omega2, shapes = eigh(stiffness)
    scale = (2.0 * np.pi * first) ** 2 / omega2[0]
    return np.sqrt(omega2 * scale) / (2.0 * np.pi), shapes


def _response(duration=1200.0, seed=0, noise=0.01):
    """Floor displacements ``(STOREYS, nsamp)`` under white noise."""
    frequency, shapes = _shear_frame()
    rng = np.random.default_rng(seed)
    nsamp = int(round(duration / DELTA))
    forces = rng.standard_normal((STOREYS, nsamp))
    modal_forces = shapes.T @ forces
    motion = np.zeros((STOREYS, nsamp))

    for k in range(STOREYS):
        omega = 2.0 * np.pi * frequency[k]
        num, den, _ = cont2discrete(
            ([1.0], [1.0, 2.0 * DAMPING * omega, omega ** 2]), DELTA, "zoh"
        )
        motion += np.outer(shapes[:, k], lfilter(num[0], den,
                                                 modal_forces[k]))

    motion += noise * motion.std() * rng.standard_normal(motion.shape)
    return motion


def _array(motion, code="E"):
    """ArrayData of the floors (vertical line, 3 m storeys)."""
    stations = [
        StationData(
            f"XX.F{i:02d}.",
            [Record(data=series, delta=DELTA, time=START,
                    sid=f"XX.F{i:02d}..HN{code}")],
            (0.0, 0.0, 3.0 * (i + 1)),
        )
        for i, series in enumerate(motion)
    ]
    return ArrayData(stations)


def _config(**changes):
    """Fixed 20 s windows on the FFT frequencies from 0.5 to 10 Hz."""
    parameters = {
        "fmin": 0.5, "fmax": 10.0, "window_type": "fixed",
        "window_length": 20.0, "window_overlap": 0.5,
        "windows_per_block": 10, "block_overlap": 0.5,
        "components": "east",
    }
    parameters.update(changes)
    return FDDConfig(**parameters)


@pytest.fixture(scope="module")
def frame_result():
    """FDD of the shear frame (shared: do not modify)."""
    return FDD(_array(_response()), _config()).process()


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------


def test_config_defaults_and_json(tmp_path):
    """Defaults, enumerations from strings and a JSON round trip."""
    config = FDDConfig(fmin=0.1, fmax=1.0, nfreq=10, normalization=None)

    assert config.normalization is FDDNormalization.NONE
    assert config.components is FDDComponents.VERTICAL
    assert config.vectors == 3 and config.subfrequencies == 1

    config.update(components="transverse", azimuth=390.0,
                  normalization="coherence", subfrequencies=5,
                  bandwidth=0.1)
    assert config.azimuth == pytest.approx(30.0)

    config.write(tmp_path / "fdd.json")
    assert FDDConfig(tmp_path / "fdd.json") == config
    assert "FDDConfig" in config.summary()


@pytest.mark.parametrize(
    "changes, error",
    [
        ({"components": "radial"}, ValueError),
        ({"subfrequencies": 3}, ValueError),
        ({"bandwidth": 0.2}, ValueError),
        ({"bandwidth": -0.1, "subfrequencies": 2}, ValueError),
        ({"vectors": 0}, ValueError),
        ({"vectors": 1.5}, TypeError),
        ({"normalization": "pearson"}, ValueError),
        ({"unknown": 1}, TypeError),
    ],
)
def test_config_invalid(changes, error):
    """Missing azimuth, inconsistent sub-frequencies, wrong values."""
    with pytest.raises(error):
        FDDConfig(fmin=0.1, fmax=1.0, nfreq=10, **changes)


# ----------------------------------------------------------------------------
# Shear frame
# ----------------------------------------------------------------------------


def test_frame_frequencies(frame_result):
    """The peaks of the first singular value are the natural frequencies."""
    expected, _ = _shear_frame()
    peaks = frame_result.peaks()

    assert peaks.size == STOREYS
    # FFT lines of the 20 s windows every 0.05 Hz: the noise can move
    # a peak by one line (checked with other seeds).
    np.testing.assert_allclose(peaks, expected, atol=0.075)


def test_frame_mode_shapes(frame_result):
    """The modes are the shapes of the frame, stable among the blocks."""
    expected, shapes = _shear_frame()

    for k, frequency in enumerate(expected):
        mode = frame_result.mode(frequency)

        assert isinstance(mode, Mode)
        assert mac(mode.shape, shapes[:, k]) > 0.99
        assert mode.mac > 0.95
        assert np.linalg.norm(mode.shape) == pytest.approx(1.0)
        assert mode.shape[np.argmax(np.abs(mode.shape))] > 0.0
        assert mode.blocks == frame_result.values[mode.frequency].shape[0]
        assert mode.damping is None and mode.method == "fdd"


def test_frame_singular_values(frame_result):
    """One dominant singular value at each mode; statistics shapes."""
    expected, _ = _shear_frame()
    frequency, mean, log_std = frame_result.singular_values()

    assert mean.shape == log_std.shape == (frequency.size, STOREYS)
    assert np.all(np.diff(mean, axis=1) <= 0.0)

    for value in expected:
        row = np.argmin(np.abs(frequency - value))
        assert mean[row, 0] > 20.0 * mean[row, 1]

    assert np.all(log_std > 0.0)


def test_frame_channels_and_times(frame_result):
    """Channels, coordinates and block times of the result."""
    assert frame_result.nchannels == STOREYS
    assert frame_result.channels[0] == "XX.F00. east"
    np.testing.assert_allclose(frame_result.coordinates[:, 2],
                               [3.0, 6.0, 9.0, 12.0])

    times = frame_result.block_times[float(frame_result.frequencies[0])]
    # 20 s windows, half overlap, blocks of 10 windows, half overlap.
    assert len(times) == 22
    assert float(times[1] - times[0]) == pytest.approx(50.0)
    assert "4 channels" in frame_result.summary()


# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------


def test_coherence_removes_channel_gains():
    """With coherence the vectors do not depend on the channel gains."""
    motion = _response(duration=300.0, seed=1)
    gains = np.array([1.0, 50.0, 0.1, 3.0])[:, None]
    config = _config(normalization="coherence", fmax=3.0)
    first = FDD(_array(motion), config).process()
    second = FDD(_array(gains * motion), config).process()

    for frequency in first.frequencies:
        np.testing.assert_allclose(
            np.abs(second.vectors[frequency][:, :, 0]),
            np.abs(first.vectors[frequency][:, :, 0]), atol=1e-8,
        )


def test_transverse_is_east_for_azimuth_zero():
    """Transverse to a northward azimuth is the east component."""
    array = _array(_response(duration=200.0, seed=2))
    east = FDD(array, _config(fmax=3.0)).process()
    transverse = FDD(
        array, _config(fmax=3.0, components="transverse", azimuth=0.0)
    ).process()

    for frequency in east.frequencies:
        np.testing.assert_allclose(transverse.values[frequency],
                                   east.values[frequency])


def test_cycles_and_subfrequencies():
    """Windows in cycles and sub-frequencies; short data warn."""
    array = _array(_response(duration=120.0, seed=3))
    config = FDDConfig(
        fmin=0.2, fmax=5.0, nfreq=8, window_length=10,
        windows_per_block=5, subfrequencies=3, bandwidth=0.05,
        components="east", vectors=2,
    )

    with pytest.warns(FDDWarning):
        result = FDD(array, config).process()

    assert result.values[float(result.frequencies[0])].shape[0] == 0
    last = float(result.frequencies[-1])
    assert result.vectors[last].shape[1:] == (STOREYS, 2)

    with pytest.raises(ValueError):
        result.mode(float(result.frequencies[0]))
    with pytest.raises(ValueError):
        result.mode(last, index=2)
    with pytest.raises(ValueError):
        FDD(array, FDDConfig(fmin=0.5, fmax=5.0)).process()


def test_processor_errors():
    """Missing array or configuration, wrong types."""
    with pytest.raises(ValueError):
        FDD(config=_config()).process()
    with pytest.raises(ValueError):
        FDD(_array(_response(duration=60.0))).process()
    with pytest.raises(TypeError):
        FDD(array="array")
    with pytest.raises(TypeError):
        FDD(config={"fmin": 1.0})


def test_without_positions():
    """
    The decomposition needs no positions.

    The coordinates of the result are NaN for the stations without a
    position, and such results can still be joined.
    """
    motion = _response(duration=300.0)
    located = _array(motion)
    stations = list(located)
    stations[2] = StationData(stations[2].sid, stations[2].records)
    array = ArrayData(stations)
    result = FDD(array, _config()).process()
    reference = FDD(located, _config()).process()

    assert not array.has_positions
    assert np.isnan(result.coordinates[2]).all()
    np.testing.assert_array_equal(result.coordinates[[0, 1, 3]],
                                  reference.coordinates[[0, 1, 3]])
    np.testing.assert_allclose(
        result.values[result.frequencies[5]],
        reference.values[reference.frequencies[5]],
    )

    joined = result.copy()
    joined.append(result)

    assert joined.nblocks == 2 * result.nblocks
    assert np.isnan(result.mode(1.5).coordinates[2]).all()


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


def test_result_files_and_append(tmp_path, frame_result):
    """Write, append, read back joined; incompatible results fail."""
    path = tmp_path / "fdd.pkl"
    frame_result.write(path)
    frame_result.write(path, append=True)
    joined = FDDResult(path)

    assert joined.nblocks == 2 * frame_result.nblocks
    np.testing.assert_array_equal(joined.frequencies,
                                  frame_result.frequencies)
    assert joined.mode(1.5).blocks == 2 * frame_result.mode(1.5).blocks

    copy = frame_result.copy()
    assert copy.nblocks == frame_result.nblocks

    with pytest.raises(FileExistsError):
        frame_result.write(path)

    other = FDD(
        _array(_response(duration=120.0, seed=4)),
        _config(fmax=3.0),
    ).process()

    with pytest.raises(ValueError):
        other.write(path, append=True)
    with pytest.raises(ValueError):
        copy.append(other)
    with pytest.raises(ValueError):
        BeamformingResult(path)
    with pytest.raises(ValueError):
        FDDResult().write(tmp_path / "empty.pkl")


def test_peaks_options(frame_result):
    """Range, count and singular value index of the peaks."""
    peaks = frame_result.peaks(fmin=3.0, fmax=7.0)

    assert np.all((peaks >= 3.0) & (peaks <= 7.0))
    np.testing.assert_allclose(frame_result.peaks(count=1), 1.5)
    # Without prominence, small fluctuations are peaks too.
    assert frame_result.peaks(prominence=1.0).size > STOREYS

    with pytest.raises(ValueError):
        frame_result.peaks(index=STOREYS)
    with pytest.raises(ValueError):
        frame_result.peaks(prominence=0.5)
    with pytest.raises(ValueError):
        FDDResult().peaks()


def test_mac():
    """Modal assurance criterion of proportional and orthogonal shapes."""
    assert mac([1.0, 2.0], [-2.0, -4.0]) == pytest.approx(1.0)
    assert mac([1.0, 0.0], [0.0, 1.0]) == pytest.approx(0.0)
    assert mac([1.0, 1j], [1j, -1.0]) == pytest.approx(1.0)

    with pytest.raises(ValueError):
        mac([1.0, 2.0], [1.0])
    with pytest.raises(ValueError):
        mac([0.0, 0.0], [1.0, 1.0])


def test_no_warnings_on_complete_data(frame_result):
    """The frame analysis fills every frequency."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        FDD(_array(_response(duration=300.0, seed=5)),
            _config(fmax=3.0)).process()

    assert all(
        values.shape[0] for values in frame_result.values.values()
    )
