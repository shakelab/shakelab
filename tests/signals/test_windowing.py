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
Unit tests for :mod:`shakelab.signals.windowing`.

Window layouts are checked against explicit loops, tapers against
``scipy.signal.get_window``, and the effective number of realizations
against a Monte Carlo estimate of the variance of averaged periodograms
of white noise (Percival and Walden, 1993, Section 6.17).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import signal

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.windowing import (
    effective_realizations,
    segments,
    split_record,
    taper,
    window_size,
    window_starts,
)

START = "2026-01-01T00:00:00"


# ----------------------------------------------------------------------------
# window_size
# ----------------------------------------------------------------------------


def test_window_size_from_duration_and_cycles():
    """Durations and cycles give the same number of samples."""
    assert window_size(0.01, duration=5.0) == 500
    assert window_size(0.01, cycles=20, frequency=4.0) == 500
    assert window_size(0.004, duration=1.0) == 250


def test_window_size_rounds_down_inexact_durations():
    """A duration between two multiples of delta is rounded down."""
    assert window_size(0.01, duration=0.0259) == 2
    assert window_size(0.01, cycles=1, frequency=3.0) == 33


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"duration": 1.0, "cycles": 2.0, "frequency": 1.0},
        {"cycles": 2.0},
        {"duration": 1.0, "frequency": 1.0},
        {"duration": 0.015},
        {"duration": -1.0},
    ],
)
def test_window_size_invalid_arguments(kwargs):
    """Inconsistent or too short windows are rejected."""
    with pytest.raises(ValueError):
        window_size(0.01, **kwargs)


def test_window_size_types():
    """Non-numerical arguments raise TypeError."""
    with pytest.raises(TypeError):
        window_size("0.01", duration=1.0)
    with pytest.raises(TypeError):
        window_size(0.01, duration=True)


# ----------------------------------------------------------------------------
# window_starts and segments
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "nsamp, size, overlap, expected",
    [
        (10, 4, None, [0, 4]),
        (10, 4, 0.5, [0, 2, 4, 6]),
        (10, 4, 0.9, [0, 1, 2, 3, 4, 5, 6]),
        (4, 4, 0.0, [0]),
        (3, 4, 0.0, []),
        (0, 1, 0.0, []),
    ],
)
def test_window_starts(nsamp, size, overlap, expected):
    """Only complete windows are returned."""
    starts = window_starts(nsamp, size, overlap=overlap)

    assert starts.dtype == np.int64
    np.testing.assert_array_equal(starts, expected)


def test_window_starts_with_step():
    """The step can be given in samples."""
    np.testing.assert_array_equal(
        window_starts(10, 4, step=3), [0, 3, 6]
    )


def test_window_starts_also_gives_blocks():
    """Blocks of windows use the same layout as windows of samples."""
    windows = window_starts(1000, 100, overlap=0.5)
    blocks = window_starts(windows.size, 8, overlap=0.5)

    assert windows.size == 19
    np.testing.assert_array_equal(blocks, [0, 4, 8])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"overlap": 1.0},
        {"overlap": -0.1},
        {"overlap": 0.5, "step": 2},
        {"step": 0},
    ],
)
def test_window_starts_invalid(kwargs):
    """Invalid overlaps and steps are rejected."""
    with pytest.raises(ValueError):
        window_starts(10, 4, **kwargs)


def test_window_starts_types():
    """Non-integer counts raise TypeError."""
    with pytest.raises(TypeError):
        window_starts(10.0, 4)
    with pytest.raises(TypeError):
        window_starts(10, 4, overlap="0.5")


@pytest.mark.parametrize("axis", [0, 1, -1])
def test_segments_match_explicit_loop(axis):
    """Segments equal the slices given by window_starts."""
    data = np.arange(3 * 50, dtype=float).reshape(3, 50)
    data = np.moveaxis(data, 1, axis)
    starts = window_starts(50, 8, overlap=0.25)

    result = segments(data, 8, overlap=0.25, axis=axis)
    expected = np.stack([
        np.moveaxis(np.take(data, range(i, i + 8), axis=axis), axis, -1)
        for i in starts
    ])

    np.testing.assert_array_equal(result, expected)


def test_segments_are_read_only_views():
    """Segments share memory with the data and cannot be modified."""
    data = np.arange(100.0)
    result = segments(data, 10, step=5)

    assert np.shares_memory(result, data)
    assert not result.flags.writeable


def test_segments_without_complete_windows():
    """An array shorter than one window gives zero windows."""
    result = segments(np.zeros((2, 5)), 8)

    assert result.shape == (0, 2, 8)


def test_segments_invalid_axis():
    """Out-of-range axes and scalars are rejected."""
    with pytest.raises(ValueError):
        segments(np.zeros(10), 4, axis=1)
    with pytest.raises(ValueError):
        segments(3.0, 1)
    with pytest.raises(TypeError):
        segments(np.zeros(10), 4, axis=0.0)


# ----------------------------------------------------------------------------
# Tapers
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("window", ["hann", ("tukey", 0.1), "boxcar"])
def test_taper_matches_get_window(window):
    """Named windows are the periodic windows of SciPy."""
    np.testing.assert_allclose(
        taper(64, window), signal.get_window(window, 64, fftbins=True)
    )


def test_taper_normalizations():
    """The rms and sum normalizations are exact."""
    rms = taper(100, ("tukey", 0.2), norm="rms")
    total = taper(100, ("tukey", 0.2), norm="sum")

    assert np.mean(rms ** 2) == pytest.approx(1.0, rel=1e-12)
    assert np.sum(total) == pytest.approx(1.0, rel=1e-12)


def test_taper_sum_norm_calibrates_amplitude():
    """With norm='sum' a sinusoid of amplitude A has a DFT of A / 2."""
    size = 200
    values = taper(size, "hann", norm="sum")
    time = np.arange(size)
    sine = 3.0 * np.cos(2.0 * np.pi * 10 * time / size)

    assert abs(np.fft.rfft(values * sine)[10]) == pytest.approx(1.5)


def test_taper_from_array():
    """An explicit array is accepted and checked."""
    np.testing.assert_array_equal(taper(3, [1.0, 2.0, 1.0]), [1, 2, 1])

    with pytest.raises(ValueError):
        taper(4, [1.0, 2.0, 1.0])
    with pytest.raises(ValueError):
        taper(3, [0.0, 0.0, 0.0], norm="rms")
    with pytest.raises(ValueError):
        taper(3, "hann", norm="peak")


# ----------------------------------------------------------------------------
# Effective realizations
# ----------------------------------------------------------------------------


def test_effective_realizations_limits():
    """Non-overlapping or single windows are independent."""
    values = taper(64, "hann")

    assert effective_realizations(values, 64, 10) == pytest.approx(10.0)
    assert effective_realizations(values, 32, 1) == 1.0


def test_effective_realizations_hann_half_overlap():
    """Analytical value for a Hann taper with 50% overlap."""
    # For a periodic Hann taper shifted by half its length rho = 1/6.
    count = 20
    expected = count / (1.0 + 2.0 * (1.0 - 1.0 / count) / 36.0)

    assert effective_realizations(
        taper(64, "hann"), 32, count
    ) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "window, overlap, count",
    [("hann", 0.5, 20), (("tukey", 0.1), 0.5, 20), ("hann", 0.75, 10)],
)
def test_effective_realizations_against_monte_carlo(window, overlap, count):
    """
    The formula predicts the variance of averaged periodograms.

    For white noise each periodogram is exponentially distributed, so
    that the effective number equals mean**2 / variance of the average
    (Percival and Walden, 1993, Section 6.17). With 20000 trials the
    Monte Carlo estimate has a relative error of about 1%: the
    tolerance is 4%.
    """
    rng = np.random.default_rng(12345)
    size = 64
    step = int(size * (1.0 - overlap))
    values = taper(size, window, norm="rms")
    nsamp = (count - 1) * step + size

    noise = rng.standard_normal((20000, nsamp))
    windows = segments(noise, size, step=step)
    coefficients = np.fft.rfft(windows * values, axis=-1)[..., 10]
    average = np.mean(np.abs(coefficients) ** 2, axis=0)
    estimate = np.mean(average) ** 2 / np.var(average)

    assert effective_realizations(values, step, count) == pytest.approx(
        estimate, rel=0.04
    )


def test_effective_realizations_invalid():
    """Invalid tapers and counts are rejected."""
    with pytest.raises(ValueError):
        effective_realizations([1.0], 1, 2)
    with pytest.raises(ValueError):
        effective_realizations(np.zeros(8), 4, 2)
    with pytest.raises(ValueError):
        effective_realizations(np.ones(8), 4, 0)
    with pytest.raises(TypeError):
        effective_realizations(np.ones(8), 4.0, 2)


# ----------------------------------------------------------------------------
# split_record
# ----------------------------------------------------------------------------


def test_split_record():
    """Windows are independent records with the right start times."""
    record = Record(
        data=np.arange(1000.0), delta=0.01, time=START,
        sid="XX.STA..HHZ", units="m/s",
    )

    pieces = split_record(record, 2.0, overlap=0.5)

    assert len(pieces) == 9
    assert pieces[3].starttime == Date(START) + 3.0
    assert pieces[3].sid == "XX.STA..HHZ"
    assert pieces[3].head.units == "m/s"
    np.testing.assert_array_equal(pieces[3].data, np.arange(300, 500))

    pieces[0].data[0] = -1.0
    assert record.data[0] == 0.0


def test_split_record_short_or_invalid():
    """Short records give no windows; invalid inputs are rejected."""
    record = Record(data=np.zeros(10), delta=0.1, time=START)

    assert split_record(record, 2.0) == []

    with pytest.raises(TypeError):
        split_record(np.zeros(10), 1.0)
    with pytest.raises(ValueError):
        split_record(Record(data=np.zeros(10)), 1.0)
