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
Unit tests for :mod:`shakelab.signals.spectrogram`.

The transforms are checked against known signals (amplitude and phase
of sinusoids, localisation of impulses, ridge of a chirp, power
spectral density of white noise), against SciPy's ``ShortTimeFFT`` and
against the defining property of the Stockwell transform (the sum over
time is the Fourier spectrum). The Spectrogram class is checked for
the metadata contract, options, cutting and edge effects.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import signal

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.spectrogram import (
    Spectrogram,
    cwt_coi,
    cwt_transform,
    morlet_wavelet,
    power_to_db,
    st_transform,
    wft_transform,
)

START = "2026-01-01T00:00:00"
DELTA = 0.01


def _cosine(amplitude=2.5, freq=5.0, phase=0.3, nsamp=4000):
    """Samples of a cosine."""
    time = np.arange(nsamp) * DELTA
    return amplitude * np.cos(2.0 * np.pi * freq * time + phase)


def _record(data, units="m/s**2"):
    """Record with metadata."""
    return Record(
        data=data, delta=DELTA, time=START, sid="XX.STA..HNZ",
        units=units,
    )


# ----------------------------------------------------------------------------
# Calibration
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("freq", [1.0, 5.0, 20.0])
def test_amplitude_scaling_of_a_sinusoid(freq):
    """A sinusoid of amplitude A gives A with every transform."""
    data = _cosine(freq=freq)
    middle = data.size // 2

    _, cwt = cwt_transform(data, DELTA, [freq])
    _, st = st_transform(data, DELTA, [freq])
    wfreq, wtime, wft = wft_transform(data, DELTA, 400)
    row = np.argmin(np.abs(wfreq - freq))

    assert abs(cwt[0, middle]) == pytest.approx(2.5, rel=1e-3)
    assert abs(st[0, middle]) == pytest.approx(2.5, rel=1e-6)
    assert abs(wft[row, wtime.size // 2]) == pytest.approx(2.5, rel=1e-6)


def test_zero_frequency_is_the_mean():
    """The zero frequency gives the mean (one-sided factor 1)."""
    data = np.full(1000, 3.0)

    _, st = st_transform(data, DELTA, [0.0])
    _, _, wft = wft_transform(data, DELTA, 100)

    np.testing.assert_allclose(np.abs(st[0]), 3.0)
    assert abs(wft[0, 5]) == pytest.approx(3.0)


def test_phase_is_common_to_the_transforms():
    """The phase is that of the signal at the time of the coefficient."""
    data = _cosine(amplitude=1.0, freq=5.0, phase=0.7)
    freq, time, wft = wft_transform(data, DELTA, 400)
    column = time.size // 2
    sample = int(round(time[column] / DELTA))
    expected = np.angle(np.exp(1j * (2 * np.pi * 5.0 * time[column] + 0.7)))

    _, cwt = cwt_transform(data, DELTA, [5.0])
    _, st = st_transform(data, DELTA, [5.0])

    assert np.angle(wft[np.argmin(np.abs(freq - 5.0)), column]) == (
        pytest.approx(expected)
    )
    assert np.angle(cwt[0, sample]) == pytest.approx(expected, abs=1e-6)
    assert np.angle(st[0, sample]) == pytest.approx(expected, abs=1e-6)


def test_psd_scaling_of_white_noise():
    """With scaling="psd" the power is the one-sided PSD (as Welch)."""
    rng = np.random.default_rng(0)
    data = 2.0 * rng.standard_normal(100000)

    _, _, psd = wft_transform(
        data, DELTA, 512, scaling="psd", power=True
    )

    # Variance 4 over 50 Hz: 0.08 units**2/Hz.
    assert psd[5:-5].mean() == pytest.approx(0.08, rel=0.02)


def test_energy_scaling_of_the_wavelet_transform():
    """With scaling="energy" the modulus grows as 1 / sqrt(f)."""
    ratio = []

    for freq in (2.0, 8.0):
        _, cwt = cwt_transform(
            _cosine(freq=freq), DELTA, [freq], scaling="energy"
        )
        ratio.append(abs(cwt[0, 2000]))

    assert ratio[0] / ratio[1] == pytest.approx(2.0, rel=1e-3)


# ----------------------------------------------------------------------------
# Windowed Fourier transform
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("length, step, nfft", [
    (128, 32, None), (100, 25, None), (128, 40, 256),
])
def test_wft_matches_scipy_short_time_fft(length, step, nfft):
    """The centred transform equals SciPy's ShortTimeFFT."""
    data = np.random.default_rng(1).standard_normal(1000)

    _, time, wft = wft_transform(
        data, DELTA, length, step=step, nfft=nfft, scaling=None
    )

    sft = signal.ShortTimeFFT(
        signal.get_window("hann", length), hop=step, fs=1 / DELTA,
        mfft=nfft, scale_to=None,
    )
    reference = sft.stft(data)
    columns = np.arange(time.size) - sft.p_min

    np.testing.assert_allclose(wft, reference[:, columns], atol=1e-10)
    np.testing.assert_allclose(time, sft.t(data.size)[columns])


def test_wft_windows_cover_the_data():
    """Centred windows reach the last sample; others stay inside."""
    data = np.ones(1000)

    _, time, _ = wft_transform(data, DELTA, 256, step=100)
    assert time[0] == 0.0
    assert time[-1] == pytest.approx(9.0)

    _, time, spectrum = wft_transform(
        data, DELTA, 256, step=100, centred=False
    )
    assert time[0] == pytest.approx(1.28)
    assert spectrum.shape == (129, 8)

    with pytest.raises(ValueError):
        wft_transform(data[:100], DELTA, 256, centred=False)


def test_wft_windows():
    """Windows by name, tuple or array; arrays of the wrong length fail."""
    data = np.random.default_rng(2).standard_normal(1000)

    by_tuple = wft_transform(data, DELTA, 128, window=("tukey", 0.25))[2]
    by_array = wft_transform(
        data, DELTA, 128, window=signal.get_window(("tukey", 0.25), 128)
    )[2]

    np.testing.assert_allclose(by_tuple, by_array)

    with pytest.raises(ValueError):
        wft_transform(data, DELTA, 128, window=np.ones(64))


# ----------------------------------------------------------------------------
# Wavelet and Stockwell transforms
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("nsamp", [1000, 1001])
def test_impulse_is_localised(nsamp):
    """The modulus peaks at the impulse, for even and odd lengths."""
    data = np.zeros(nsamp)
    data[300] = 1.0

    _, cwt = cwt_transform(data, DELTA, [2.0, 10.0, 30.0])
    _, st = st_transform(data, DELTA, [2.0, 10.0, 30.0])

    np.testing.assert_array_equal(np.argmax(np.abs(cwt), axis=1), 300)
    np.testing.assert_array_equal(np.argmax(np.abs(st), axis=1), 300)


def test_chirp_ridge():
    """The ridge of the wavelet transform follows a linear chirp."""
    time = np.arange(4000) * DELTA
    data = signal.chirp(time, 1.0, time[-1], 20.0)
    freq = np.linspace(0.5, 25.0, 200)

    _, cwt = cwt_transform(data, DELTA, freq)

    for sample in (1000, 2000, 3000):
        expected = 1.0 + 19.0 * time[sample] / time[-1]
        ridge = freq[np.argmax(np.abs(cwt[:, sample]))]
        assert ridge == pytest.approx(expected, rel=0.03)


def test_stockwell_sums_to_the_fourier_spectrum():
    """Summed over time, the Stockwell transform is the DFT."""
    data = np.random.default_rng(3).standard_normal(256)

    freq, st = st_transform(data, DELTA, scaling=None)

    np.testing.assert_allclose(freq, np.fft.rfftfreq(256, DELTA))
    np.testing.assert_allclose(st.sum(axis=1), np.fft.rfft(data))


def test_stockwell_frequencies_are_dft_frequencies():
    """Requested frequencies become the nearest DFT frequencies."""
    freq, st = st_transform(np.zeros(1000), DELTA, [1.04, 1.06, 1.1, 50])

    np.testing.assert_allclose(freq, [1.0, 1.1, 50.0])
    assert st.shape == (3, 1000)

    with pytest.raises(ValueError):
        st_transform(np.zeros(1000), DELTA, [60.0])


@pytest.mark.parametrize("axis", [0, 1, 2, -1])
def test_axis_layout(axis):
    """Frequency first, time axis in its original position."""
    data = np.random.default_rng(4).standard_normal((3, 5, 512))
    moved = np.moveaxis(data, 2, axis)
    position = axis % 3 + 1

    cwt = cwt_transform(moved, DELTA, [5.0, 10.0], axis=axis)[1]
    st = st_transform(moved, DELTA, [5.0, 10.0], axis=axis)[1]
    wft = wft_transform(moved, DELTA, 64, axis=axis)[2]

    expected_shape = list(moved.shape)
    assert cwt.shape == (2,) + moved.shape
    assert st.shape == (2,) + moved.shape
    expected_shape[axis] = 16
    assert wft.shape == (33,) + tuple(expected_shape)

    single = cwt_transform(data[1, 2], DELTA, [5.0, 10.0])[1]
    np.testing.assert_allclose(
        np.moveaxis(cwt, position, -1)[:, 1, 2], single
    )


def test_morlet_wavelet_norms():
    """Unit energy, or an envelope that sums to one."""
    energy = morlet_wavelet(101, 8.0)
    amplitude = morlet_wavelet(101, 8.0, norm="amplitude")

    assert np.sum(np.abs(energy) ** 2) == pytest.approx(1.0)
    assert np.sum(np.abs(amplitude)) == pytest.approx(1.0)
    assert np.argmax(np.abs(energy)) == 50

    with pytest.raises(ValueError):
        morlet_wavelet(101, 8.0, norm="peak")


def test_cone_of_influence():
    """sqrt(2) times the wavelet standard deviation, in seconds."""
    np.testing.assert_allclose(
        cwt_coi([1.0, 2.0]), np.sqrt(2) * 6 / (2 * np.pi) / [1.0, 2.0]
    )


@pytest.mark.parametrize("call", [
    lambda: cwt_transform(np.zeros(100), 0.0, [1.0]),
    lambda: cwt_transform(np.zeros(100), DELTA, [0.0]),
    lambda: cwt_transform(np.zeros(100), DELTA, [60.0]),
    lambda: cwt_transform(np.zeros(100), DELTA, [1.0], axis=1),
    lambda: cwt_transform(np.zeros(100), DELTA, [1.0], scaling="psd"),
    lambda: st_transform(np.zeros(100), DELTA, [-1.0]),
    lambda: st_transform(np.float64(1.0), DELTA),
    lambda: wft_transform(np.zeros(100), DELTA, 0),
    lambda: wft_transform(np.zeros(100), DELTA, 10, step=0),
    lambda: wft_transform(np.zeros(100), DELTA, 10, nfft=8),
])
def test_invalid_arguments(call):
    """Out-of-range arguments raise ValueError."""
    with pytest.raises(ValueError):
        call()


def test_invalid_types():
    """Arguments of the wrong type raise TypeError."""
    with pytest.raises(TypeError):
        wft_transform(np.zeros(100), "0.01", 10)
    with pytest.raises(TypeError):
        wft_transform(np.zeros(100), DELTA, 10.0)
    with pytest.raises(TypeError):
        st_transform(np.zeros(100), DELTA, axis=1.0)


def test_power_to_db():
    """Decibels relative to a reference, with an optional floor."""
    np.testing.assert_allclose(
        power_to_db([0.001, 1.0, 100.0], floor=30), [-10.0, 0.0, 20.0]
    )
    np.testing.assert_allclose(
        power_to_db([1.0, 100.0], reference=100.0), [-20.0, 0.0]
    )
    assert power_to_db([0.0])[0] == -np.inf

    with pytest.raises(ValueError):
        power_to_db([1.0], floor=0)


# ----------------------------------------------------------------------------
# Spectrogram
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["wft", "cwt", "st"])
def test_spectrogram_from_record(method):
    """The spectrogram keeps the metadata and is calibrated."""
    record = _record(_cosine(freq=5.0))

    spec = record.spectrogram(method)

    assert isinstance(spec, Spectrogram)
    assert spec.method == method
    assert spec.sid == "XX.STA..HNZ"
    assert spec.units == "m/s**2"
    assert spec.starttime == Date(START)
    assert spec.data.shape == (spec.nfreq, spec.ntime)
    assert spec.frequency.shape == (spec.nfreq,)
    assert spec.time.shape == (spec.ntime,)
    assert spec.coi.shape == (spec.nfreq,)

    row = np.argmin(np.abs(spec.frequency - 5.0))
    assert spec.amplitude[row, spec.ntime // 2] == pytest.approx(
        2.5, rel=0.05
    )

    # The header is a copy.
    spec.head.sid = "XX.OTHER..HNZ"
    assert record.sid == "XX.STA..HNZ"


def test_spectrogram_options():
    """Options in seconds and frequency grids."""
    record = _record(_cosine())

    spec = record.spectrogram("wft", window=2.0, overlap=0.75, fmax=20)
    assert spec.parameters["window_length"] == 200
    assert spec.parameters["step"] == 50
    assert spec.time[1] == pytest.approx(0.5)
    assert spec.frequency.max() <= 20
    np.testing.assert_allclose(spec.coi, 1.0)

    spec = record.spectrogram(
        "cwt", fmin=1, fmax=10, nfreq=5, spacing="linear"
    )
    np.testing.assert_allclose(spec.frequency, [1, 3.25, 5.5, 7.75, 10])
    assert spec.ntime == record.nsamp

    spec = record.spectrogram("cwt", fmin=1, fmax=16, nfreq=5)
    np.testing.assert_allclose(spec.frequency, [1, 2, 4, 8, 16])

    spec = record.spectrogram("st", fmin=0, fmax=1, nfreq=3)
    assert spec.frequency[0] == 0.0


def test_spectrogram_invalid_options():
    """Options of another method, unknown methods and bad records."""
    record = _record(_cosine())

    with pytest.raises(TypeError):
        record.spectrogram("cwt", window=2.0)
    with pytest.raises(TypeError):
        record.spectrogram("wft", nfreq=10)
    with pytest.raises(TypeError):
        Spectrogram.from_record(np.zeros(10))
    with pytest.raises(ValueError):
        record.spectrogram("fft")
    with pytest.raises(ValueError):
        record.spectrogram("cwt", fmin=10, fmax=1)
    with pytest.raises(ValueError):
        record.spectrogram("wft", overlap=1.0)
    with pytest.raises(ValueError):
        Record(data=np.zeros(10)).spectrogram()


def test_spectrogram_amplitude_power_db():
    """Derived quantities of the coefficients."""
    spec = Spectrogram(
        [[3 + 4j, 0.0]], frequency=[1.0], time=[0.0, 1.0]
    )

    np.testing.assert_allclose(spec.amplitude, [[5.0, 0.0]])
    np.testing.assert_allclose(spec.power, [[25.0, 0.0]])
    np.testing.assert_allclose(spec.phase[0, 0], np.arctan2(4, 3))
    np.testing.assert_allclose(
        spec.db(floor=20), [[10 * np.log10(25), 10 * np.log10(0.25)]]
    )

    with pytest.raises(ValueError):
        Spectrogram([[1.0, 2.0]], frequency=[1.0], time=[0.0])


def test_spectrogram_cut():
    """Cut by time (Date, string or seconds) and by frequency."""
    record = _record(_cosine())
    spec = record.spectrogram("cwt", fmin=1, fmax=16, nfreq=5)

    cut = spec.cut(
        Date(START) + 10.0, "2026-01-01T00:00:20", fmin=2, fmax=8,
        inplace=False,
    )
    np.testing.assert_allclose(cut.frequency, [2, 4, 8])
    assert cut.time[0] == pytest.approx(10.0)
    assert cut.time[-1] == pytest.approx(20.0)
    assert cut.data.shape == (3, 1001)
    assert cut.coi.shape == (3,)
    assert spec.ntime == record.nsamp

    assert spec.cut(5.0, 6.0) is None
    assert spec.time[0] == pytest.approx(5.0)

    with pytest.raises(ValueError):
        spec.cut(fmin=100)
    with pytest.raises(ValueError):
        spec.cut(6.0, 5.0)
    with pytest.raises(TypeError):
        spec.cut([1.0])


def test_spectrogram_cone_of_influence():
    """The mask marks the edges of the record, also after a cut."""
    spec = _record(_cosine()).spectrogram("cwt", fmin=1, fmax=4, nfreq=3)
    mask = spec.coi_mask

    assert mask.shape == spec.data.shape
    assert mask[:, 0].all() and mask[:, -1].all()
    assert not mask[:, spec.ntime // 2].any()

    # Lower frequencies have wider edges.
    assert mask[0].sum() > mask[-1].sum()

    cut = spec.cut(30.0, 39.99, inplace=False)
    assert cut.coi_mask[:, -1].all()
    assert not cut.coi_mask[:, 0].any()


def test_spectrogram_copy_and_repr():
    """Copies are independent; text descriptions name the method."""
    spec = _record(_cosine()).spectrogram("st", nfreq=4)
    copy = spec.copy()
    copy.data[:] = 0

    assert spec.amplitude.max() > 0
    assert "method='st'" in repr(spec)
    assert "method st" in str(spec)
    assert Spectrogram().nfreq == 0
