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
response by channel and epoch, the transfer functions of the stages and
of stage sets, units, and the application of responses to spectra and
records (water level, round trip).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from shakelab.signals.base import Record
from shakelab.signals.response import (
    ResponseCollection,
    StageFIR,
    StageGain,
    StagePoleZero,
    StagePolynomial,
    StageSet,
    StreamResponse,
    select_response,
)
from shakelab.signals.spectrum import Spectrum

DATA = Path(__file__).parent / "data"

SID = "IV.MODE..HNZ"
TIME = "2012-05-20T02:03:23"


@pytest.fixture(scope="module")
def mode_response():
    """Instrument responses of station IV.MODE."""
    return ResponseCollection(str(DATA / "MODE.xml"))


def _velocimeter(sid="XX.TEST..HHZ"):
    """
    Stream response of a 1 Hz velocimeter (counts per m/s).

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
    stage_set = StageSet("2000-01-01T00:00:00", None)
    stage_set.append(paz)
    stage_set.append(StageGain(sensitivity=1000.0, stage_number=2))

    stream_response = StreamResponse(sid)
    stream_response.append(stage_set)
    return stream_response


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
# Reading and selection
# ----------------------------------------------------------------------------


def test_read_stationxml(mode_response):
    """StationXML gives one stream response per channel."""
    assert mode_response.sid == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ"
    ]
    assert len(mode_response) == 3

    for ftype in ("xml", "sxml"):
        collection = ResponseCollection(str(DATA / "MODE.xml"), ftype=ftype)
        assert collection.sid == mode_response.sid


def test_read_unknown_type_is_rejected(tmp_path):
    """Unknown file types raise ValueError."""
    with pytest.raises(ValueError):
        ResponseCollection(str(tmp_path / "file.txt"))


def test_missing_channel_and_epoch(mode_response):
    """A missing channel or epoch raises KeyError."""
    with pytest.raises(KeyError):
        mode_response["XX.NONE..HHZ"]

    with pytest.raises(KeyError):
        mode_response[SID].extract("1990-01-01T00:00:00")

    stage_set = mode_response.extract(SID, TIME)
    assert isinstance(stage_set, StageSet)


def test_select_response(mode_response):
    """Responses are selected by channel and time."""
    expected = mode_response[SID].extract(TIME)

    assert select_response(mode_response, SID, TIME) is expected
    assert select_response(mode_response[SID], SID, TIME) is expected
    assert select_response(expected) is expected

    stage = expected.stage[0]
    assert select_response(stage, "any", None) is stage


def test_select_response_errors(mode_response):
    """Invalid selections raise explicit errors."""
    with pytest.raises(ValueError):
        select_response(mode_response[SID], "IV.MODE..HNE", TIME)

    with pytest.raises(ValueError):
        select_response(mode_response[SID], SID, None)

    with pytest.raises(ValueError):
        select_response(mode_response, None, TIME)

    with pytest.raises(KeyError):
        select_response(mode_response, "XX.NONE..HHZ", TIME)

    with pytest.raises(TypeError):
        select_response(1.0, SID, TIME)


def test_append_merges_epochs():
    """Stream responses with the same identifier are merged."""
    collection = ResponseCollection()
    collection.append(_velocimeter())
    collection.append(_velocimeter())

    assert len(collection) == 1
    assert len(collection["XX.TEST..HHZ"]) == 2

    with pytest.raises(ValueError):
        collection.append(StageSet())


def test_json_round_trip(tmp_path, mode_response):
    """Writing and reading JSON keeps stages and epochs."""
    collection = ResponseCollection()
    for resp in mode_response:
        collection.append(resp)

    late = _velocimeter(SID).stageset[0]
    late.starttime = None
    early = mode_response[SID].stageset[0].copy()
    early.starttime = None
    early.endtime = mode_response[SID].stageset[0].starttime
    collection[SID].append(early)

    path = tmp_path / "response.json"
    collection.write(str(path))
    loaded = ResponseCollection(str(path))

    assert loaded.sid == collection.sid
    assert len(loaded[SID]) == 2

    freq = np.array([0.5, 1.0, 10.0])
    np.testing.assert_allclose(
        loaded.extract(SID, TIME).transfer_function(freq),
        mode_response.extract(SID, TIME).transfer_function(freq),
    )


def test_json_coefficients_stage(tmp_path):
    """A JSON 'coefficients' stage is read as a FIR stage."""
    path = tmp_path / "response.json"
    path.write_text(
        '[{"sid": "XX.TEST..HHZ", "starttime": null, "endtime": null,'
        ' "stages": [{"stage_number": 1, "stage_type": "coefficients",'
        ' "data": {"numerator": [0.5, 0.5]}}]}]'
    )
    loaded = ResponseCollection(str(path))

    assert isinstance(loaded.extract("XX.TEST..HHZ").stage[0], StageFIR)


# ----------------------------------------------------------------------------
# Transfer functions
# ----------------------------------------------------------------------------


def test_mode_transfer_function(mode_response):
    """The total response matches the declared sensitivity."""
    stage_set = mode_response.extract(SID, TIME)
    total = stage_set.transfer_function(np.array([0.2, 1.0, 5.0]))

    np.testing.assert_allclose(np.abs(total), 320770.0, rtol=1e-3)
    assert stage_set.input_units == "m/s**2"


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


def test_fir_transfer_function():
    """A two-point average passes zero frequency and stops Nyquist."""
    stage = StageFIR(numerator=[0.5, 0.5])

    with pytest.raises(ValueError):
        stage.transfer_function(np.array([0.0]))

    values = stage.transfer_function(np.array([0.0, 50.0]), delta=0.01)
    np.testing.assert_allclose(np.abs(values), [1.0, 0.0], atol=1e-12)


def test_polynomial_transfer_function():
    """The polynomial stage evaluates a rational function of s."""
    stage = StagePolynomial(numerator=[1.0], denominator=[1.0, 1.0])
    freq = np.array([1.0])
    s = 2j * np.pi * freq

    np.testing.assert_allclose(stage.transfer_function(freq), 1 / (1 + s))


def test_inactive_stages_are_ignored():
    """Only active stage types contribute to a stage set."""
    stage_set = StageSet()
    stage_set.append(StageGain(sensitivity=2.0))
    stage_set.append(StageFIR(numerator=[0.0, 1.0]))

    values = stage_set.transfer_function(np.array([1.0, 10.0]), 0.01)
    np.testing.assert_allclose(values, 2.0)


def test_stage_set_units():
    """Units come from the first and the last stage, if declared."""
    stage_set = _velocimeter().stageset[0]

    assert stage_set.input_units == "m/s"
    assert stage_set.output_units is None

    stage_set.stage[-1] = StagePoleZero(
        input_units="V", output_units="counts", poles=[]
    )
    assert stage_set.output_units == "counts"

    assert StageSet().input_units is None


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
    stage_set = response.stageset[0]
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


def test_mode_correction_units(mode_response):
    """Correcting the accelerometer gives m/s**2."""
    record = Record(
        data=np.random.default_rng(0).standard_normal(1000),
        delta=0.01,
        time=TIME,
        sid=SID,
    )
    record.deconvolve_response(mode_response, waterlevel=40.0,
                               method="smooth")

    assert record.head.units == "m/s**2"
