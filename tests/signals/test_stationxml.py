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
Unit tests for :mod:`shakelab.signals.stationxml`.

The tests read the StationXML of station IV.MODE and small synthetic
documents covering every response element: symmetric FIR filters,
poles and zeros of each type, digital and analog coefficients, response
lists, polynomials, decimation and time correction, channels without
response and unknown elements.
"""

from __future__ import annotations

import io
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.signals.response import (
    ResponseWarning,
    StageCoefficients,
    StagePoleZero,
    StagePolynomial,
    StageResponseList,
    StageUnsupported,
)
from shakelab.signals.stationxml import (
    iter_channels,
    read_stationxml,
    stationxml_to_dict,
)
from shakelab.libutils.xml import read_xml

DATA = Path(__file__).parent / "data"
MODE = DATA / "MODE.xml"
SID = "XX.TEST.00.HHZ"
TIME = "2020-06-01T00:00:00"


def _document(response, channel_extra="", second_epoch=False):
    """StationXML with one channel whose ``<Response>`` is given."""
    epoch = (
        '<Channel code="HHZ" locationCode="00" '
        'startDate="{start}" {end}>'
        "<SampleRate>100</SampleRate>"
        f"{channel_extra}"
        f"{response}"
        "</Channel>"
    )
    channels = epoch.format(
        start="2020-01-01T00:00:00", end='endDate="2021-01-01T00:00:00"'
    )
    if second_epoch:
        channels += epoch.format(start="2021-01-01T00:00:00", end="")
    return (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" '
        'schemaVersion="1.1">'
        '<Network code="XX"><Station code="TEST">'
        f"{channels}"
        "</Station></Network></FDSNStationXML>"
    )


def _stages(response):
    """Read a document and return the stages of its first epoch."""
    collection = read_stationxml(_document(response))
    return collection.extract(SID, TIME).stage


def _stage(content, number=1, gain=1.0):
    """A ``<Stage>`` with one element and a gain."""
    return (
        f'<Stage number="{number}">{content}'
        f"<StageGain><Value>{gain}</Value><Frequency>1</Frequency>"
        "</StageGain></Stage>"
    )


UNITS = (
    "<InputUnits><Name>m/s</Name></InputUnits>"
    "<OutputUnits><Name>V</Name></OutputUnits>"
)
DIGITAL_UNITS = (
    "<InputUnits><Name>count</Name></InputUnits>"
    "<OutputUnits><Name>count</Name></OutputUnits>"
)


def _fir(symmetry, coefficients):
    """A ``<FIR>`` element with ``NumeratorCoefficient`` values."""
    values = "".join(
        f'<NumeratorCoefficient i="{i}">{c}</NumeratorCoefficient>'
        for i, c in enumerate(coefficients)
    )
    return (
        f"<FIR>{DIGITAL_UNITS}<Symmetry>{symmetry}</Symmetry>{values}"
        "</FIR>"
    )


# ----------------------------------------------------------------------------
# MODE
# ----------------------------------------------------------------------------


def test_read_mode():
    """The IV.MODE accelerometer is read with all its stages."""
    collection = read_stationxml(MODE)

    assert collection.sid == [
        "IV.MODE..HNZ", "IV.MODE..HNN", "IV.MODE..HNE"
    ]
    stage_set = collection.extract("IV.MODE..HNZ", "2012-05-20T02:03:23")

    assert [s.stage_type for s in stage_set.stage] == [
        "paz", "gain", "coefficients", "gain", "coefficients", "gain",
        "coefficients", "gain",
    ]
    assert [s.stage_number for s in stage_set.stage] == [
        1, 1, 2, 2, 3, 3, 4, 4
    ]
    assert stage_set.sensitivity.value == 321707.0
    assert stage_set.sensitivity.frequency == 10.625
    assert stage_set.sample_rate == 100.0
    assert stage_set.units() == ("m/s**2", "count")

    fir = stage_set.stage[4]
    assert fir.input_sample_rate == 800.0
    assert fir.decimation_factor == 4
    assert fir.delay == pytest.approx(0.029375)
    assert fir.correction == pytest.approx(0.029375)


def test_read_sources():
    """Path, text, bytes and file objects give the same responses."""
    text = MODE.read_text(encoding="utf-8")
    for source in (str(MODE), MODE, text, text.encode("utf-8"),
                   io.BytesIO(text.encode("utf-8"))):
        assert read_stationxml(source).sid == [
            "IV.MODE..HNZ", "IV.MODE..HNN", "IV.MODE..HNE"
        ]


def test_not_stationxml_is_rejected():
    """Other XML documents raise ValueError."""
    with pytest.raises(ValueError):
        read_stationxml("<quakeml/>")


def test_iter_channels_and_dict():
    """Channels are found by code; the dictionary keeps strings."""
    root = read_xml(_document(""))
    assert [sid for sid, _ in iter_channels(root)] == [SID]

    content = stationxml_to_dict(_document(""))
    assert content["tag"] == "FDSNStationXML"


# ----------------------------------------------------------------------------
# Elements
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "symmetry, expected",
    [
        ("NONE", [0.1, 0.4, 0.5]),
        ("EVEN", [0.1, 0.4, 0.5, 0.5, 0.4, 0.1]),
        ("ODD", [0.1, 0.4, 0.5, 0.4, 0.1]),
    ],
)
def test_fir_symmetry(symmetry, expected):
    """Symmetric FIR filters are expanded (FDSN definitions)."""
    stages = _stages(
        "<Response>" + _stage(_fir(symmetry, [0.1, 0.4, 0.5])) +
        "</Response>"
    )

    assert isinstance(stages[0], StageCoefficients)
    np.testing.assert_allclose(stages[0].numerator, expected)


def test_fir_unknown_symmetry():
    """An unknown symmetry raises ValueError."""
    with pytest.raises(ValueError):
        _stages("<Response>" + _stage(_fir("HALF", [1.0])) + "</Response>")


@pytest.mark.parametrize(
    "pz_type",
    ["LAPLACE (RADIANS/SECOND)", "LAPLACE (HERTZ)", "DIGITAL (Z-TRANSFORM)"],
)
def test_poles_zeros(pz_type):
    """Poles, zeros, normalization and type are read."""
    content = (
        f"<PolesZeros>{UNITS}"
        f"<PzTransferFunctionType>{pz_type}</PzTransferFunctionType>"
        "<NormalizationFactor>2.5</NormalizationFactor>"
        "<NormalizationFrequency>1</NormalizationFrequency>"
        '<Zero number="0"><Real>0</Real><Imaginary>0</Imaginary></Zero>'
        '<Pole number="0"><Real>-1</Real><Imaginary>2</Imaginary></Pole>'
        '<Pole number="1"><Real>-1</Real><Imaginary>-2</Imaginary></Pole>'
        "</PolesZeros>"
    )
    paz = _stages("<Response>" + _stage(content) + "</Response>")[0]

    assert isinstance(paz, StagePoleZero)
    assert paz.pz_type == pz_type
    assert paz.normalization_factor == 2.5
    np.testing.assert_allclose(paz.poles, [-1 + 2j, -1 - 2j])
    np.testing.assert_allclose(paz.zeros, [0j])
    assert (paz.input_units, paz.output_units) == ("m/s", "V")


def test_decimation_and_correction():
    """Decimation parameters are attached to digital stages."""
    content = (
        f"<Coefficients>{DIGITAL_UNITS}"
        "<CfTransferFunctionType>DIGITAL</CfTransferFunctionType>"
        "<Numerator>0.0</Numerator><Numerator>1.0</Numerator>"
        "</Coefficients>"
        "<Decimation><InputSampleRate>100</InputSampleRate>"
        "<Factor>1</Factor><Offset>0</Offset><Delay>0.01</Delay>"
        "<Correction>0.01</Correction></Decimation>"
    )
    fir = _stages("<Response>" + _stage(content) + "</Response>")[0]

    assert (fir.input_sample_rate, fir.decimation_factor) == (100.0, 1)
    assert (fir.delay, fir.correction) == (0.01, 0.01)
    np.testing.assert_allclose(
        fir.transfer_function(np.array([1.0, 20.0])), 1.0
    )


def test_empty_coefficients_is_unit_response():
    """An analog-to-digital converter stage has H = 1."""
    content = (
        "<Coefficients><InputUnits><Name>V</Name></InputUnits>"
        "<OutputUnits><Name>count</Name></OutputUnits>"
        "<CfTransferFunctionType>DIGITAL</CfTransferFunctionType>"
        "</Coefficients>"
    )
    stages = _stages("<Response>" + _stage(content, gain=1e6) +
                     "</Response>")

    np.testing.assert_allclose(
        stages[0].transfer_function(np.array([1.0]), delta=0.01), 1.0
    )
    assert stages[1].sensitivity == 1e6


@pytest.mark.parametrize(
    "cf_type, scale",
    [("ANALOG (RADIANS/SECOND)", 1.0), ("ANALOG (HERTZ)", 2.0 * np.pi)],
)
def test_analog_coefficients(cf_type, scale):
    """Analog coefficients are ascending powers of s (SEED eq. 7)."""
    # H(s) = 1 / (1 + s/w0), with w0 = 10 in the units of the type.
    content = (
        f"<Coefficients>{UNITS}"
        f"<CfTransferFunctionType>{cf_type}</CfTransferFunctionType>"
        "<Numerator>1</Numerator>"
        "<Denominator>1</Denominator><Denominator>0.1</Denominator>"
        "</Coefficients>"
    )
    stage = _stages("<Response>" + _stage(content) + "</Response>")[0]
    freq = np.array([0.5, 1.0, 5.0])
    s = 2j * np.pi * freq / scale

    assert stage.cf_type == cf_type
    np.testing.assert_allclose(
        stage.transfer_function(freq), 1.0 / (1.0 + 0.1 * s)
    )


def test_unknown_coefficients_type_is_unsupported():
    """An unknown coefficients type gives an unsupported stage."""
    content = (
        f"<Coefficients>{UNITS}"
        "<CfTransferFunctionType>OTHER</CfTransferFunctionType>"
        "</Coefficients>"
    )
    with pytest.warns(ResponseWarning):
        stages = _stages("<Response>" + _stage(content) + "</Response>")

    assert isinstance(stages[0], StageUnsupported)
    assert stages[0].output_units == "V"
    with pytest.raises(NotImplementedError):
        stages[0].transfer_function(np.array([1.0]))


def test_response_list():
    """Response tables are read in order."""
    rows = "".join(
        f"<ResponseListElement><Frequency>{f}</Frequency>"
        f"<Amplitude>{a}</Amplitude><Phase>{p}</Phase>"
        "</ResponseListElement>"
        for f, a, p in [(1, 1, 0), (10, 100, 90)]
    )
    content = f"<ResponseList>{UNITS}{rows}</ResponseList>"
    stage = _stages("<Response>" + _stage(content) + "</Response>")[0]

    assert isinstance(stage, StageResponseList)
    np.testing.assert_allclose(stage.frequency, [1.0, 10.0])
    np.testing.assert_allclose(stage.phase, [0.0, 90.0])


def test_polynomials():
    """Stage and instrument polynomials are read as non-linear stages."""
    polynomial = (
        "<{tag}>"
        "<InputUnits><Name>degC</Name></InputUnits>"
        "<OutputUnits><Name>V</Name></OutputUnits>"
        "<ApproximationType>MACLAURIN</ApproximationType>"
        "<FrequencyLowerBound>0</FrequencyLowerBound>"
        "<FrequencyUpperBound>0</FrequencyUpperBound>"
        "<ApproximationLowerBound>-20</ApproximationLowerBound>"
        "<ApproximationUpperBound>60</ApproximationUpperBound>"
        "<MaximumError>0.01</MaximumError>"
        '<Coefficient number="1">0.5</Coefficient>'
        '<Coefficient number="0">1.0</Coefficient>'
        '<Coefficient number="2">0.1</Coefficient>'
        "</{tag}>"
    )
    response = (
        "<Response>"
        + polynomial.format(tag="InstrumentPolynomial")
        + _stage(polynomial.format(tag="Polynomial"))
        + "</Response>"
    )
    collection = read_stationxml(_document(response))
    stage_set = collection.extract(SID, TIME)

    for stage in (stage_set.instrument_polynomial, stage_set.stage[0]):
        assert isinstance(stage, StagePolynomial)
        np.testing.assert_allclose(stage.coefficients, [1.0, 0.5, 0.1])
        assert stage.approximation_upper_bound == 60.0
        np.testing.assert_allclose(stage.evaluate([2.0]), [2.4])

    with pytest.raises(ValueError):
        stage_set.transfer_function(np.array([1.0]))


def test_sensitivity_and_epochs():
    """Each epoch is a stage set with its own sensitivity."""
    response = (
        "<Response><InstrumentSensitivity><Value>1000</Value>"
        "<Frequency>1</Frequency>"
        "<InputUnits><Name>m/s</Name></InputUnits>"
        "<OutputUnits><Name>count</Name></OutputUnits>"
        "</InstrumentSensitivity></Response>"
    )
    collection = read_stationxml(_document(response, second_epoch=True))

    assert len(collection) == 1
    assert len(collection[SID]) == 2
    late = collection.extract(SID, "2022-01-01T00:00:00")
    assert late.sensitivity.value == 1000.0
    assert late.units() == ("m/s", "count")
    assert late.endtime is None


def test_channel_without_response_and_unknown_elements():
    """Missing responses and unknown elements are reported."""
    with pytest.warns(ResponseWarning):
        assert len(read_stationxml(_document(""))) == 0

    with pytest.warns(ResponseWarning):
        stages = _stages(
            "<Response>" + _stage("<Mystery/>", gain=3.0) + "</Response>"
        )
    assert [s.stage_type for s in stages] == ["gain"]


def test_mode_reads_without_warnings():
    """A complete document is read silently."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_stationxml(MODE)
