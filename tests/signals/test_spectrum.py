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
Unit tests for :mod:`shakelab.signals.spectrum`.

The tests cover the construction from spectral values and from time
series, the round trip with Record for even and odd lengths, the
metadata contract, arithmetic, filtering, interpolation, inversion and
smoothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record
from shakelab.signals.spectrum import Spectrum

START = "2026-01-01T00:00:00"


def _noise(nsamp=1000, delta=0.01, seed=0):
    """Record of white noise with metadata."""
    rng = np.random.default_rng(seed)
    return Record(
        data=rng.standard_normal(nsamp),
        delta=delta,
        time=START,
        sid="XX.STA..HHZ",
        units="m/s",
    )


# ----------------------------------------------------------------------------
# Construction
# ----------------------------------------------------------------------------


def test_spectrum_from_values():
    """The constructor takes spectral values, e.g. a transfer function."""
    spectrum = Spectrum(np.ones(6), delta=0.1, nsamp=10)

    assert spectrum.nfreq == len(spectrum) == 6
    assert spectrum.nsamp == 10
    assert spectrum.dfreq == pytest.approx(1.0)
    np.testing.assert_allclose(spectrum.frequency, [0, 1, 2, 3, 4, 5])
    assert spectrum.data.dtype.kind == "c"
    assert Spectrum(np.ones(6), delta=0.1).nsamp == 10
    assert Spectrum(np.ones(6), delta=0.1, nsamp=11).nsamp == 11


def test_spectrum_length_mismatch():
    """nsamp must match the number of spectral values."""
    with pytest.raises(ValueError):
        Spectrum(np.ones(6), delta=0.1, nsamp=20)
    with pytest.raises(ValueError):
        Spectrum(np.ones((2, 3)), delta=0.1)


def test_spectrum_from_samples():
    """The spectrum of a sine has amplitude T/2 at its frequency."""
    delta = 0.01
    t = np.arange(1000) * delta
    spectrum = Spectrum.from_samples(np.sin(2 * np.pi * 5.0 * t), delta)
    k = np.argmin(np.abs(spectrum.frequency - 5.0))

    assert spectrum.amplitude[k] == pytest.approx(5.0)
    assert spectrum.delta == delta
    assert spectrum.nsamp == 1000


@pytest.mark.parametrize("nsamp", [1000, 1001])
def test_record_round_trip(nsamp):
    """Record -> Spectrum -> Record restores data and metadata."""
    record = _noise(nsamp)
    spectrum = record.to_spectrum()
    restored = spectrum.to_record()

    np.testing.assert_allclose(restored.data, record.data, atol=1e-12)
    assert spectrum.nsamp == restored.nsamp == nsamp
    assert restored.time == Date(START)
    assert restored.sid == record.sid
    assert restored.head.units == "m/s"
    assert restored.head._parent is restored

    other = Record()
    other.from_spectrum(spectrum)
    np.testing.assert_allclose(other.data, record.data, atol=1e-12)


def test_from_record_copies_metadata_without_parent():
    """The spectrum header is a copy, not linked to the record."""
    record = _noise()
    spectrum = Spectrum.from_record(record)

    assert spectrum.head is not record.head
    assert spectrum.head._parent is None
    spectrum.head.sid = "other"
    assert record.sid == "XX.STA..HHZ"


def test_from_record_errors():
    """Only Records with a sampling interval can be transformed."""
    with pytest.raises(TypeError):
        Spectrum.from_record(np.zeros(10))
    with pytest.raises(ValueError):
        Record(data=np.zeros(10)).to_spectrum()


def test_copy_and_text(capsys):
    """Copies are independent; text representations are compact."""
    spectrum = _noise().to_spectrum()
    copied = spectrum.copy()
    copied.data[1] = 0.0
    copied.head.sid = "other"

    assert spectrum.data[1] != 0.0
    assert spectrum.head.sid == "XX.STA..HHZ"
    assert repr(spectrum).startswith("Spectrum(sid='XX.STA..HHZ'")

    spectrum.info()
    assert "501 frequencies" in capsys.readouterr().out


# ----------------------------------------------------------------------------
# Arithmetic
# ----------------------------------------------------------------------------


def test_spectrum_arithmetic():
    """Products and ratios with scalars and compatible spectra."""
    spectrum = _noise().to_spectrum()

    np.testing.assert_allclose((spectrum * 2).data, 2 * spectrum.data)
    np.testing.assert_allclose((2 * spectrum).data, 2 * spectrum.data)
    np.testing.assert_allclose(
        (spectrum / np.float32(2)).data, spectrum.data / 2
    )
    ratio = spectrum / spectrum
    np.testing.assert_allclose(ratio.data, 1.0)


def test_spectrum_arithmetic_errors():
    """Incompatible axes, zero and unsupported operands are rejected."""
    spectrum = _noise(1000).to_spectrum()

    with pytest.raises(ValueError):
        spectrum / _noise(500).to_spectrum()
    with pytest.raises(ZeroDivisionError):
        spectrum / 0
    with pytest.raises(TypeError):
        spectrum * "2"


# ----------------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------------


def test_spectrum_filter_butterworth():
    """The Butterworth filter halves the power at the corner."""
    spectrum = Spectrum(np.ones(501), delta=0.01, nsamp=1000)
    spectrum.filter(lowpass=10.0, order=3)
    k = np.argmin(np.abs(spectrum.frequency - 10.0))

    assert spectrum.amplitude[k] == pytest.approx(1 / np.sqrt(2))
    assert spectrum.amplitude[0] == pytest.approx(1.0)


def test_spectrum_filter_sharp_and_errors():
    """The sharp filter zeroes the stop band; bad options raise."""
    spectrum = Spectrum(np.ones(501), delta=0.01, nsamp=1000)
    spectrum.filter(highpass=5.0, lowpass=10.0, method="sharp")
    freq = spectrum.frequency
    band = (freq >= 5.0) & (freq <= 10.0)

    assert np.all(spectrum.data[band] == 1.0)
    assert np.all(spectrum.data[~band] == 0.0)

    with pytest.raises(ValueError):
        spectrum.filter()
    with pytest.raises(ValueError):
        spectrum.filter(lowpass=1.0, method="chebyshev")


def test_spectrum_filter_matches_record_band():
    """Filtering in frequency removes the components outside the band."""
    t = np.arange(4000) * 0.01
    low = np.sin(2 * np.pi * 0.5 * t)
    record = Record(data=low + np.sin(2 * np.pi * 20.0 * t), delta=0.01)

    spectrum = record.to_spectrum()
    spectrum.filter(lowpass=5.0, order=8)
    filtered = spectrum.to_samples()

    np.testing.assert_allclose(filtered, low, atol=1e-3)


def test_spectrum_interpolate():
    """Linear interpolation of real and imaginary parts."""
    spectrum = Spectrum(
        np.array([0, 1 + 1j, 2 + 4j]), delta=0.25, nsamp=4
    )

    np.testing.assert_allclose(
        spectrum.interpolate([0.5, 1.5]), [0.5 + 0.5j, 1.5 + 2.5j]
    )
    with pytest.raises(ValueError):
        spectrum.interpolate([3.0])


def test_spectrum_invert():
    """The regularized inverse of a flat spectrum is its reciprocal."""
    spectrum = Spectrum(np.full(6, 2.0 + 0j), delta=0.1, nsamp=10)
    spectrum.invert(waterlevel=120, method="sharp")

    np.testing.assert_allclose(spectrum.data, 0.5)


@pytest.mark.parametrize("method", ["gaussian", "konno-ohmachi"])
def test_spectrum_logsmooth(method):
    """Smoothing reduces the scatter of a noise spectrum."""
    spectrum = _noise(4096).to_spectrum()
    original = spectrum.copy()

    spectrum.logsmooth(method=method)

    band = slice(400, 2000)
    assert np.std(np.log(spectrum.amplitude[band])) < 0.5 * np.std(
        np.log(original.amplitude[band])
    )
    assert spectrum.data[0] == pytest.approx(original.data[0])


def test_spectrum_logsmooth_amplitude_only():
    """With phase=False the phase is kept."""
    spectrum = _noise(1024).to_spectrum()
    original = spectrum.copy()

    spectrum.logsmooth(phase=False)

    np.testing.assert_allclose(
        np.angle(spectrum.data[1:]), np.angle(original.data[1:]),
        atol=1e-12,
    )
