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

import csv
import dataclasses
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
    modes_to_csv,
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
# Enhanced FDD
# ----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def frame_modes():
    """Enhanced FDD of the four modes of the shear frame (shared)."""
    fdd = FDD(_array(_response()), _config())
    return fdd, fdd.efdd(_shear_frame()[0])


def test_efdd_of_the_shear_frame(frame_modes):
    """
    Natural frequencies, damping and shapes of the shear frame.

    With 20 minutes of data and 2% damping the estimates of six
    realizations scatter by about 10% in damping and 0.3% in frequency
    around the true values (diary, 2026-11-16): the bounds are four
    times wider than that scatter.
    """
    frequency, shapes = _shear_frame()
    _, modes = frame_modes

    for k, mode in enumerate(modes):
        assert mode.method == "efdd"
        assert mode.frequency == pytest.approx(frequency[k], rel=0.005)
        assert mode.damping == pytest.approx(DAMPING, abs=0.006)
        assert mac(mode.shape, shapes[:, k]) > 0.999
        assert mode.mac > 0.95

        details = mode.efdd
        assert details.duration == pytest.approx(
            100.0 / frequency[k], abs=DELTA
        )
        assert details.bell[np.argmin(np.abs(
            details.frequency - details.peak_frequency
        ))]
        assert details.used.sum() >= 3
        assert np.all(np.abs(details.extrema[details.used, 1]) <= 0.9)
        assert details.correlation[0] == 1.0
        assert mode.damping == pytest.approx(
            details.decrement / np.hypot(2.0 * np.pi, details.decrement)
        )
        assert mode.frequency > details.damped_frequency

    assert "EFDD" in modes[0].summary()


def test_efdd_options_and_errors(frame_modes):
    """Fixed duration, levels and invalid options."""
    fdd, _ = frame_modes
    mode = fdd.efdd(4.3, duration=40.0, levels=(0.2, 0.8))[0]

    assert mode.efdd.duration == pytest.approx(40.0)
    assert mode.frequency == pytest.approx(_shear_frame()[0][1], rel=0.01)

    with pytest.raises(TypeError):
        fdd.efdd(4.3, cycles=50, duration=40.0)
    for options in ({"mac": 0.0}, {"levels": (0.9, 0.1)}, {"cycles": -1},
                    {"duration": 2000.0}):
        with pytest.raises(ValueError):
            fdd.efdd(4.3, **options)
    with pytest.raises(ValueError, match="Nyquist"):
        fdd.efdd(30.0)
    with pytest.raises(ValueError):
        FDD(config=_config()).efdd(4.3)


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


def test_result_files_and_append(tmp_path, frame_result):
    """Write, append, read back joined; incompatible results fail."""
    path = tmp_path / "fdd.jsonl"
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
        FDDResult().write(tmp_path / "empty.jsonl")


def test_light_result_files(tmp_path, frame_result):
    """
    Without its arrays a file gives the first value and vector.

    The first singular value and vector of every block are in the
    JSON Lines file, exactly; all the values and the kept vectors are
    in the .npz files. Peaks and the first mode are the same; light
    and complete results do not join.
    """
    path = tmp_path / "fdd.jsonl"
    frame_result.write(path)
    complete = FDDResult(path)
    light = FDDResult(path, arrays=False)
    frequency = frame_result.frequencies[20]

    np.testing.assert_array_equal(complete.values[frequency],
                                  frame_result.values[frequency])
    np.testing.assert_array_equal(complete.vectors[frequency],
                                  frame_result.vectors[frequency])
    np.testing.assert_array_equal(light.values[frequency][:, 0],
                                  frame_result.values[frequency][:, 0])
    np.testing.assert_array_equal(light.vectors[frequency][:, :, 0],
                                  frame_result.vectors[frequency][:, :, 0])
    np.testing.assert_array_equal(light.peaks(), frame_result.peaks())
    # Same vectors, stored in arrays of different layout: the sums of
    # NumPy may differ in the last digit on some platforms.
    assert light.mode(1.5).mac == pytest.approx(frame_result.mode(1.5).mac,
                                                rel=1e-12)
    assert light.singular_values()[1].shape[1] == 1
    assert light.block_times == frame_result.block_times

    with pytest.raises(ValueError):
        light.mode(1.5, index=1)
    with pytest.raises(ValueError):
        light.peaks(index=1)
    with pytest.raises(ValueError, match="numbers of singular"):
        light.copy().append(complete)

    frame_result.write(tmp_path / "only.jsonl", arrays=False)

    assert not (tmp_path / "only_arrays").exists()
    assert FDDResult(tmp_path / "only.jsonl").values[frequency].shape[1] == 1

    for item in (tmp_path / "fdd_arrays").iterdir():
        item.unlink()

    with pytest.raises(FileNotFoundError, match="arrays=False"):
        FDDResult(path)


def test_csv_of_the_first_vector(tmp_path, frame_result):
    """One row per frequency and block, real and imaginary parts."""
    table = tmp_path / "fdd.csv"
    frame_result.to_csv(table)
    rows = table.read_text().splitlines()
    header = rows[0].split(",")

    assert header[:4] == ["frequency", "block", "time", "value"]
    assert header[4:6] == ["XX.F00. east real", "XX.F00. east imag"]
    assert len(header) == 4 + 2 * STOREYS
    assert len(rows) == frame_result.nblocks + 1


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


def test_modes_to_csv(tmp_path, frame_modes, frame_result):
    """Table of modes of the FDD and of the enhanced FDD; errors."""
    _, modes = frame_modes
    path = tmp_path / "modes.csv"
    modes_to_csv(path, modes)

    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert list(rows[0])[:8] == [
        "frequency", "method", "index", "value", "value_log_std", "mac",
        "blocks", "damping",
    ]
    assert list(rows[0])[8:10] == ["XX.F00. east", "XX.F00. east std"]
    assert [float(row["damping"]) for row in rows] == pytest.approx(
        [mode.damping for mode in modes]
    )
    assert [float(row["XX.F03. east"]) for row in rows] == pytest.approx(
        [mode.shape[3] for mode in modes]
    )

    with pytest.raises(FileExistsError):
        modes_to_csv(path, modes)

    plain = [frame_result.mode(value) for value in frame_result.peaks()]
    modes_to_csv(path, plain, owrite=True)

    with path.open(encoding="utf-8", newline="") as stream:
        assert all(row["damping"] == "" for row in csv.DictReader(stream))

    other = dataclasses.replace(plain[0], stations=("A",) * STOREYS)

    with pytest.raises(ValueError):
        modes_to_csv(path, [plain[1], other], owrite=True)
    with pytest.raises(ValueError):
        modes_to_csv(path, [], owrite=True)
    with pytest.raises(TypeError):
        modes_to_csv(path, [1.0], owrite=True)


def test_no_warnings_on_complete_data(frame_result):
    """The frame analysis fills every frequency."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        FDD(_array(_response(duration=300.0, seed=5)),
            _config(fmax=3.0)).process()

    assert all(
        values.shape[0] for values in frame_result.values.values()
    )
