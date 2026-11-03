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
Unit tests for :mod:`shakelab.signals.response`.

The tests cover reading StationXML and JSON files, the selection of the
response by channel and epoch, the transfer functions of the stages
(analog, digital at their own sampling rate, tables, polynomials) and
of stage sets, units, sensitivity, validation of the metadata and the
application of responses to spectra and records (water level, round
trip, delay of the digital filters).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from shakelab.signals.base import Record
from shakelab.signals.metadata import (
    ChannelEpoch,
    MetadataCollection,
    select_sensitivity,
)
from shakelab.signals.response import (
    ResponseWarning,
    Sensitivity,
    StageCoefficients,
    StageGain,
    StagePoleZero,
    StagePolynomial,
    StageResponseList,
    StageSet,
    StageUnsupported,
    expand_fir_symmetry,
    same_units,
)
from shakelab.signals.spectrum import Spectrum

DATA = Path(__file__).parent / "data"

SID = "IV.MODE..HNZ"
TIME = "2012-05-20T02:03:23"


@pytest.fixture(scope="module")
def mode_metadata():
    """Channel metadata (with responses) of station IV.MODE."""
    return MetadataCollection(str(DATA / "MODE.xml"))


def _velocimeter(sid="XX.TEST..HHZ"):
    """
    Channel epoch of a 1 Hz velocimeter (counts per m/s).

    Two zeros at the origin and two poles at 1 Hz (damping 0.707),
    normalized at 5 Hz, followed by a gain of 1000.
    """
    corner = 2.0 * np.pi * 1.0
    poles = corner * np.array([-0.707 + 0.707j, -0.707 - 0.707j])
    paz = StagePoleZero(
        input_units="m/s",
        output_units="V",
        normalization_factor=1.0,
        poles=poles,
        zeros=[0j, 0j],
        stage_number=1,
    )
    stage_set = StageSet()
    stage_set.append(paz)
    stage_set.append(StageGain(sensitivity=1000.0, stage_number=2))

    return ChannelEpoch(sid, starttime="2000-01-01T00:00:00",
                        response=stage_set)


def _record(sid="XX.TEST..HHZ", nsamp=2000, delta=0.01):
    """Band-limited signal (2-10 Hz sines with a taper)."""
    t = np.arange(nsamp) * delta
    data = (
        np.sin(2.0 * np.pi * 3.0 * t)
        + 0.5 * np.sin(2.0 * np.pi * 7.0 * t)
    ) * np.hanning(nsamp)
    return Record(data=data, delta=delta, time="2026-01-01T00:00:00",
                  sid=sid)


# ----------------------------------------------------------------------------
# Transfer functions
# ----------------------------------------------------------------------------


def test_mode_transfer_function(mode_metadata):
    """The response matches the declared sensitivity and units."""
    stage_set = mode_metadata.get(SID, TIME).response

    # Analog part: flat, equal to the product of the stage gains.
    freq = np.array([0.2, 1.0, 5.0])
    analog = stage_set.transfer_function(freq, stages=("gain", "paz"))
    np.testing.assert_allclose(np.abs(analog), 1.0197 * 315000, rtol=1e-3)

    # Complete response: the sensitivity is declared at 10.625 Hz.
    full = stage_set.transfer_function(np.array([10.625]))
    np.testing.assert_allclose(np.abs(full), 321707.0, rtol=1e-4)

    assert stage_set.units() == ("m/s**2", "count")
    assert stage_set.units(("gain", "paz")) == ("m/s**2", "V")
    assert stage_set.sensitivity.value == 321707.0
    assert stage_set.sample_rate == 100.0


def test_mode_digital_stages(mode_metadata):
    """The FIR stages carry their own sampling rate and delay."""
    stage_set = mode_metadata.get(SID, TIME).response
    firs = [s for s in stage_set.stage if isinstance(s, StageCoefficients)]

    assert [s.input_sample_rate for s in firs] == [800.0, 800.0, 200.0]
    assert [s.decimation_factor for s in firs] == [1, 4, 2]
    assert [s.correction for s in firs] == [0.0, 0.029375, 0.3125]
    assert [s.numerator.size for s in firs] == [1, 48, 126]

    # Above the passband of the last FIR (100 Hz -> anti-alias cut).
    full = stage_set.transfer_function(np.array([49.9]))
    analog = stage_set.transfer_function(
        np.array([49.9]), stages=("gain", "paz")
    )
    assert np.abs(full[0]) < 0.05 * np.abs(analog[0])


def test_gain_transfer_function():
    """A gain stage is constant and real."""
    stage = StageGain(sensitivity=2.5)
    values = stage.transfer_function(np.array([0.0, 1.0, 10.0]))

    np.testing.assert_array_equal(values, 2.5 + 0j)


def test_paz_without_zeros_is_not_modified():
    """Missing zeros are treated as empty without changing the stage."""
    stage = StagePoleZero(poles=[-1.0 + 0j], normalization_factor=1.0)
    values = stage.transfer_function(np.array([0.0]))

    assert stage.zeros is None
    np.testing.assert_allclose(values, 1.0)


def test_paz_hertz_matches_radians():
    """Poles and zeros in Hz give the same response as in rad/s."""
    pole_hz = -1.0 + 0.5j
    freq = np.array([0.1, 1.0, 3.0])

    hertz = StagePoleZero(
        pz_type="LAPLACE (HERTZ)",
        poles=[pole_hz],
        normalization_factor=1.0,
    )
    radians = StagePoleZero(
        pz_type="LAPLACE (RADIANS/SECOND)",
        poles=[2.0 * np.pi * pole_hz],
        normalization_factor=2.0 * np.pi,
    )

    np.testing.assert_allclose(
        hertz.transfer_function(freq), radians.transfer_function(freq)
    )


def test_digital_paz_matches_fir():
    """``0.5 (z + 1) / z`` equals the FIR ``[0.5, 0.5]``."""
    freq = np.array([0.0, 10.0, 25.0, 50.0])
    paz = StagePoleZero(
        pz_type="DIGITAL (Z-TRANSFORM)",
        normalization_factor=0.5,
        zeros=[-1.0 + 0j],
        poles=[0j],
        input_sample_rate=100.0,
    )
    fir = StageCoefficients(numerator=[0.5, 0.5], input_sample_rate=100.0)

    np.testing.assert_allclose(
        paz.transfer_function(freq), fir.transfer_function(freq),
        atol=1e-12,
    )

    paz.input_sample_rate = None
    with pytest.raises(ValueError):
        paz.transfer_function(freq)


def test_fir_transfer_function():
    """A two-point average passes zero frequency and stops Nyquist."""
    stage = StageCoefficients(numerator=[0.5, 0.5])

    with pytest.raises(ValueError):
        stage.transfer_function(np.array([0.0]))

    values = stage.transfer_function(np.array([0.0, 50.0]), delta=0.01)
    np.testing.assert_allclose(np.abs(values), [1.0, 0.0], atol=1e-12)


def test_fir_uses_its_own_sampling_rate():
    """The stage rate takes precedence over the data sampling."""
    stage = StageCoefficients(numerator=[0.5, 0.5], input_sample_rate=200.0)
    values = stage.transfer_function(np.array([0.0, 100.0]), delta=0.01)

    np.testing.assert_allclose(np.abs(values), [1.0, 0.0], atol=1e-12)


def test_fir_time_correction():
    """A correction equal to the delay of a pure delay gives H = 1."""
    stage = StageCoefficients(
        numerator=[0.0, 1.0], input_sample_rate=100.0, delay=0.01,
    )
    freq = np.array([1.0, 10.0, 30.0])
    np.testing.assert_allclose(
        stage.transfer_function(freq),
        np.exp(-2j * np.pi * freq * 0.01),
    )

    stage.correction = 0.01
    np.testing.assert_allclose(stage.transfer_function(freq), 1.0)


def test_response_list_interpolation():
    """Amplitude is interpolated in log-log, phase linearly."""
    stage = StageResponseList(
        frequency=[10.0, 1.0],
        amplitude=[100.0, 1.0],
        phase=[90.0, 0.0],
    )
    freq = np.array([1.0, np.sqrt(10.0), 10.0])
    values = stage.transfer_function(freq)

    np.testing.assert_allclose(np.abs(values), [1.0, 10.0, 100.0])
    np.testing.assert_allclose(
        np.degrees(np.angle(values)), 90.0 * (freq - 1.0) / 9.0
    )

    with pytest.warns(ResponseWarning):
        outside = stage.transfer_function(np.array([0.5, 20.0]))
    np.testing.assert_array_equal(outside, 0.0)


def test_polynomial_stage():
    """
    The polynomial gives the input from the output (SEED blockette 62).

    Example of the SEED manual: a pressure transducer with
    ``p = 600 + 100 * V`` mbar has a forward gain of 0.01 V/mbar.
    """
    linear = StagePolynomial(
        coefficients=[600.0, 100.0], input_units="mbar", output_units="V"
    )

    np.testing.assert_allclose(linear.evaluate([0.0, 5.0]), [600.0, 1100.0])
    np.testing.assert_allclose(
        linear.transfer_function(np.array([0.01, 0.1])), 0.01
    )

    offset = StagePolynomial(coefficients=[1.0])
    with pytest.raises(ValueError):
        offset.transfer_function(np.array([1.0]))

    nonlinear = StagePolynomial(coefficients=[0.0, 1.0, 0.1])
    with pytest.raises(ValueError):
        nonlinear.transfer_function(np.array([1.0]))


def test_unsupported_stage():
    """An unsupported stage cannot be evaluated, but can be excluded."""
    stage_set = StageSet()
    stage_set.append(StageGain(sensitivity=2.0))
    stage_set.append(StageUnsupported(element="Coefficients (OTHER)"))

    with pytest.raises(NotImplementedError):
        stage_set.transfer_function(np.array([1.0]))

    values = stage_set.transfer_function(np.array([1.0]), stages="gain")
    np.testing.assert_allclose(values, 2.0)


def test_stage_selection():
    """The complete response is the default; subsets are selectable."""
    stage_set = StageSet()
    stage_set.append(StageGain(sensitivity=2.0))
    stage_set.append(StageCoefficients(numerator=[0.0, 1.0]))
    freq = np.array([1.0, 10.0])

    np.testing.assert_allclose(
        stage_set.transfer_function(freq, 0.01),
        2.0 * np.exp(-2j * np.pi * freq * 0.01),
    )
    np.testing.assert_allclose(
        stage_set.transfer_function(freq, 0.01, stages=["gain"]), 2.0
    )

    with pytest.raises(ValueError):
        stage_set.transfer_function(freq, 0.01, stages=("other",))


def test_stage_set_units():
    """Units come from the first and last stages that declare them."""
    stage_set = _velocimeter().response

    assert stage_set.input_units == "m/s"
    assert stage_set.output_units == "V"

    stage_set.append(StageCoefficients(input_units="V", output_units="counts"))
    assert stage_set.output_units == "counts"
    assert stage_set.units(("gain", "paz")) == ("m/s", "V")

    empty = StageSet()
    assert empty.units() == (None, None)
    empty.sensitivity = Sensitivity(1.0, 1.0, "m/s", "count")
    assert empty.units() == ("m/s", "count")


def test_same_units():
    """Unit names are compared without case and count/counts."""
    assert same_units("M/S**2", "m/s**2")
    assert same_units("COUNTS", "count")
    assert not same_units("m/s", "m/s**2")


def test_select_sensitivity(mode_metadata):
    """The overall sensitivity is selected like the response."""
    sens = select_sensitivity(mode_metadata, SID, TIME)
    assert (sens.value, sens.input_units) == (321707.0, "m/s**2")

    gain = select_sensitivity(StageGain(sensitivity=5.0, frequency=1.0))
    assert (gain.value, gain.frequency) == (5.0, 1.0)

    with pytest.raises(ValueError):
        select_sensitivity(StageSet())

    with pytest.raises(TypeError):
        select_sensitivity(StagePoleZero(poles=[-1 + 0j]))


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------


def test_validate_mode(mode_metadata):
    """The metadata of MODE are consistent."""
    assert mode_metadata.validate() == []
    assert mode_metadata.validate(rtol=1e-5) == []


def test_validate_asymmetric_fir():
    """A FIR whose delay implies linear phase must be symmetric."""
    fir = StageCoefficients(
        numerator=[0.25, 0.5, 0.25], input_sample_rate=100.0,
        delay=0.01, stage_number=3,
    )
    assert fir.validate() == []

    fir.numerator = [0.25, 0.5, 0.025]
    issues = fir.validate()
    assert len(issues) == 1
    assert "stage 3" in issues[0] and "not symmetric" in issues[0]

    # Without the linear-phase delay, asymmetry is legitimate.
    fir.delay = 0.0
    assert fir.validate() == []


def test_validate_paz_normalization():
    """``|A0 * H(fn)|`` must be 1 at the normalization frequency."""
    paz = StagePoleZero(
        normalization_factor=10.0, normalization_frequency=0.0,
        poles=[-10.0 + 0j], stage_number=1,
    )
    assert paz.validate() == []

    paz.normalization_factor = 10.1
    assert "poles-and-zeros" in paz.validate()[0]
    assert paz.validate(rtol=0.02) == []


def test_validate_consistent_response():
    """A consistent response has no issues."""
    stage_set = _consistent()
    assert stage_set.validate() == []


def test_validate_reports_inconsistencies():
    """Gains, units and decimation chain are checked."""
    stage_set = _consistent()
    stage_set.sensitivity.value = 2000.0
    stage_set.stage[2].input_units = "mV"
    stage_set.sample_rate = 50.0

    issues = " | ".join(stage_set.validate())

    assert "product of stage gains" in issues
    assert "total response amplitude" in issues
    assert "input units 'mV'" in issues
    assert "decimation chain ends at 100.0 Hz" in issues


def _consistent():
    """Accelerometer, digitizer and FIR with consistent metadata."""
    stage_set = StageSet(sample_rate=100.0)
    stage_set.append(StagePoleZero(
        input_units="m/s**2", output_units="V",
        normalization_factor=1e6, normalization_frequency=0.0,
        poles=[-1e6 + 0j], stage_number=1,
    ))
    stage_set.append(StageGain(sensitivity=10.0, stage_number=1))
    stage_set.append(StageCoefficients(
        input_units="V", output_units="count", stage_number=2,
        input_sample_rate=200.0, decimation_factor=2, delay=0.0025,
        numerator=[0.25, 0.5, 0.25],
    ))
    stage_set.append(StageGain(sensitivity=100.0, stage_number=2))
    stage_set.sensitivity = Sensitivity(1000.0, 0.0, "m/s**2", "count")
    return stage_set


# ----------------------------------------------------------------------------
# Application to data
# ----------------------------------------------------------------------------


def test_spectrum_response_round_trip():
    """Convolution followed by deconvolution restores the spectrum."""
    record = _record()
    spectrum = record.to_spectrum()
    original = spectrum.data.copy()

    spectrum.convolve_response(_velocimeter())
    spectrum.deconvolve_response(_velocimeter(), waterlevel=100.0)

    band = (spectrum.frequency > 0.5) & (spectrum.frequency < 40.0)
    np.testing.assert_allclose(
        spectrum.data[band], original[band], rtol=1e-6, atol=1e-9
    )


def test_record_response_round_trip():
    """Convolution and deconvolution of a record restore the signal."""
    record = _record()
    corrected = record.copy()

    corrected.convolve_response(_velocimeter())
    corrected.deconvolve_response(_velocimeter())

    np.testing.assert_allclose(corrected.data, record.data, atol=1e-3)
    assert corrected.head.units == "m/s"


def test_waterlevel_limits_amplification():
    """The inverse amplitude is limited by the water level."""
    spectrum = Spectrum(np.ones(1001, dtype=complex), delta=0.01, nsamp=2000)
    spectrum.head.sid = "XX.TEST..HHZ"
    spectrum.head.time = "2026-01-01T00:00:00"

    response = _velocimeter()
    stage_set = response.response
    peak = np.abs(stage_set.transfer_function(spectrum.frequency)).max()

    for waterlevel in (20.0, 40.0):
        corrected = spectrum.copy()
        corrected.deconvolve_response(response, waterlevel=waterlevel)
        limit = 10.0 ** (waterlevel / 20.0) / peak

        assert np.abs(corrected.data).max() == pytest.approx(limit)
        assert corrected.data[0] == 0.0


def test_response_of_another_channel_is_rejected():
    """A stream response is not applied to another channel."""
    record = _record(sid="XX.OTHER..HHZ")

    with pytest.raises(ValueError):
        record.deconvolve_response(_velocimeter())


def test_mode_correction_units(mode_metadata):
    """Correcting the accelerometer gives m/s**2."""
    record = Record(
        data=np.random.default_rng(0).standard_normal(1000),
        delta=0.01,
        time=TIME,
        sid=SID,
    )
    record.deconvolve_response(mode_metadata, waterlevel=40.0,
                               method="smooth")

    assert record.head.units == "m/s**2"


def test_remove_sensitivity(mode_metadata):
    """Division by the sensitivity, without frequency processing."""
    record = Record(data=np.arange(5.0) * 321707.0, delta=0.01,
                    time=TIME, sid=SID)
    record.remove_sensitivity(mode_metadata)

    np.testing.assert_allclose(record.data, np.arange(5.0))
    assert record.head.units == "m/s**2"

    spectrum = Spectrum(np.full(3, 321707.0 + 0j), delta=0.01, nsamp=4)
    spectrum.head.sid = SID
    spectrum.head.time = TIME
    spectrum.remove_sensitivity(mode_metadata)
    np.testing.assert_allclose(spectrum.data, 1.0)


def _lag(mode, stages_a, stages_b):
    """Lag (samples) between two corrections of a random signal."""
    rng = np.random.default_rng(3)
    data = np.convolve(rng.standard_normal(3000), np.hanning(9), "same")
    data *= np.hanning(data.size)
    record = Record(data=data, delta=0.01, time=TIME, sid=SID)

    first = record.copy()
    first.deconvolve_response(mode, stages=stages_a)
    second = record.copy()
    second.deconvolve_response(mode, stages=stages_b)

    lags = np.arange(-100, 101)
    xcorr = [np.dot(second.data[100:-100],
                    np.roll(first.data, lag)[100:-100]) for lag in lags]
    return int(lags[int(np.argmax(xcorr))])


def test_fir_delay_is_cancelled_by_time_correction(mode_metadata):
    """
    The FIR delays are cancelled by the declared time correction.

    MODE declares delays of 0.029375 and 0.3125 s, corrected in the
    time stamps by the datalogger: the complete correction does not
    shift the signal. Without the correction it would advance it by
    about 34 samples at 100 Hz.
    """
    analog = ("gain", "paz")
    assert _lag(mode_metadata, None, analog) == 0

    collection = mode_metadata.select(SID)
    for stage in collection.get(SID, TIME).response.stage:
        if getattr(stage, "correction", None):
            stage.correction = 0.0
    assert _lag(collection, None, analog) == 34


@pytest.mark.parametrize(
    "symmetry, expected",
    [
        ("NONE", [0.1, 0.4, 0.5]),
        ("odd", [0.1, 0.4, 0.5, 0.4, 0.1]),
        ("EVEN", [0.1, 0.4, 0.5, 0.5, 0.4, 0.1]),
        (None, [0.1, 0.4, 0.5]),
    ],
)
def test_expand_fir_symmetry(symmetry, expected):
    """Symmetric FIR filters are expanded from their first half."""
    assert expand_fir_symmetry([0.1, 0.4, 0.5], symmetry) == expected


def test_expand_fir_symmetry_unknown():
    """An unknown symmetry raises ValueError."""
    with pytest.raises(ValueError):
        expand_fir_symmetry([1.0], "HALF")


def test_paz_without_poles():
    """Poles and zeros may be absent: the stage is a pure gain."""
    stage = StagePoleZero(normalization_factor=3.0)
    np.testing.assert_allclose(stage.transfer_function([0.0, 1.0]), 3.0)


def _velocity_stage_set(gain_frequency):
    """Velocity sensor normalized at 1 Hz, gain declared elsewhere."""
    poles = [-0.037 + 0.037j, -0.037 - 0.037j]
    paz = StagePoleZero(
        pz_type="LAPLACE (RADIANS/SECOND)", poles=poles,
        zeros=[0j, 0j], normalization_frequency=1.0, stage_number=1,
    )
    paz.normalization_factor = 1.0 / abs(paz.transfer_function([1.0])[0])
    stage_set = StageSet()
    stage_set.append([
        paz,
        StageGain(sensitivity=1500.0, frequency=gain_frequency,
                  stage_number=1),
    ])
    return stage_set


def test_stage_normalized_at_gain_frequency():
    """
    The stage amplitude equals the gain at the gain frequency.

    SEED blockette 58: the gain applies at its frequency. With the
    poles and zeros normalized at 1 Hz and the gain declared at
    0.02 Hz, the stage is rescaled by ``1 / |T(0.02)|`` (as evalresp).
    """
    stage_set = _velocity_stage_set(0.02)
    paz = stage_set.stage[0]
    np.testing.assert_allclose(
        abs(stage_set.transfer_function([0.02])[0]), 1500.0, rtol=1e-12
    )
    # The selected transfer element alone has unit amplitude there.
    np.testing.assert_allclose(
        abs(stage_set.transfer_function([0.02], stages="paz")[0]), 1.0,
        rtol=1e-12,
    )
    # Shape unchanged: ratio to the declared stage is constant.
    frequency = np.array([0.01, 0.1, 1.0, 10.0])
    ratio = stage_set.transfer_function(frequency) / (
        1500.0 * paz.transfer_function(frequency)
    )
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-12)
    assert abs(ratio[0] - 1.0) > 1e-3


def test_stage_consistent_gain_unchanged():
    """With the gain at the normalization frequency nothing changes."""
    stage_set = _velocity_stage_set(1.0)
    frequency = np.array([0.01, 0.1, 1.0, 10.0])
    np.testing.assert_allclose(
        stage_set.transfer_function(frequency),
        1500.0 * stage_set.stage[0].transfer_function(frequency),
        rtol=1e-12,
    )


def test_stage_gain_at_zero_amplitude_warns():
    """A gain at a frequency where the stage is zero is not used."""
    stage_set = _velocity_stage_set(0.0)
    with pytest.warns(ResponseWarning, match="not normalized"):
        value = stage_set.transfer_function([1.0])
    np.testing.assert_allclose(abs(value[0]), 1500.0, rtol=1e-12)


def test_digital_stage_normalized_at_zero_frequency():
    """A digital filter with gain at 0 Hz gets unit gain in DC."""
    stage_set = StageSet()
    stage_set.append([
        StageCoefficients(numerator=[0.2, 0.4, 0.2], input_sample_rate=100.0,
                          stage_number=1),
        StageGain(sensitivity=1.0, frequency=0.0, stage_number=1),
    ])
    np.testing.assert_allclose(stage_set.transfer_function([0.0]), 1.0)
