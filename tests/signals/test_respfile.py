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
Unit tests for :mod:`shakelab.signals.respfile`.

The tests compare the RESP and StationXML versions of the same nominal
response (EarthScope NRL, sensor and datalogger of IV.MODE), compare
the response of IU.ANMO.00.BHZ with the values computed by evalresp,
and cover each blockette with small synthetic RESP texts.
"""

from __future__ import annotations

import io
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.signals.response import (
    ResponseCollection,
    ResponseWarning,
    StageCoefficients,
    StageGain,
    StagePoleZero,
    StagePolynomial,
    StageResponseList,
    StageUnsupported,
)
from shakelab.signals.respfile import parse_resp_time, read_resp
from shakelab.signals.stationxml import read_stationxml

DATA = Path(__file__).parent / "data" / "RESP"
NRL_RESP = DATA / "NRL_MODE_HNZ.resp"
NRL_XML = DATA / "NRL_MODE_HNZ.xml"
ANMO = DATA / "RESP.IU.ANMO.00.BHZ"
SID = "XX.TEST.00.HHZ"
SENSITIVITY = 8.6473e8  # stage 0 of the second epoch


def _resp(*blocks, location="00", start="2020,001", end="No Ending Time"):
    """RESP text of one channel epoch with the given blockette lines."""
    header = [
        "#  synthetic RESP",
        "B050F03     Station:     TEST",
        "B050F16     Network:     XX",
        f"B052F03     Location:    {location}",
        "B052F04     Channel:     HHZ",
        f"B052F22     Start date:  {start}",
        f"B052F23     End date:    {end}",
    ]
    return "\n".join(header + [line for block in blocks for line in block])


def _stages(*blocks, **kwargs):
    """Stage set read from a synthetic RESP text."""
    collection = read_resp(_resp(*blocks, **kwargs))
    return collection[SID].stageset[0]


def _paz(code="A", number=1, poles=((-1.0, 1.0), (-1.0, -1.0)),
         zeros=((0.0, 0.0),)):
    lines = [
        f"B053F03     Transfer function type:    {code}",
        f"B053F04     Stage sequence number:     {number}",
        "B053F05     Response in units lookup:  M/S - Velocity in m/s",
        "B053F06     Response out units lookup: V - Volts",
        "B053F07     A0 normalization factor:   2.0",
        "B053F08     Normalization frequency:   1.0",
        f"B053F09     Number of zeroes:          {len(zeros)}",
        f"B053F14     Number of poles:           {len(poles)}",
    ]
    lines += [
        f"B053F10-13  {i}  {re}  {im}  0.0  0.0"
        for i, (re, im) in enumerate(zeros)
    ]
    lines += [
        f"B053F15-18  {i}  {re}  {im}  0.0  0.0"
        for i, (re, im) in enumerate(poles)
    ]
    return lines


def _gain(number, value, frequency="1.0 HZ"):
    return [
        f"B058F03     Stage sequence number:     {number}",
        f"B058F04     Gain:                      {value}",
        f"B058F05     Frequency of gain:         {frequency}",
        "B058F06     Number of calibrations:    0",
    ]


def _decimation(number, rate, factor=1, delay=0.0, correction=0.0):
    return [
        f"B057F03     Stage sequence number:     {number}",
        f"B057F04     Input sample rate (HZ):    {rate}",
        f"B057F05     Decimation factor:         {factor:05d}",
        "B057F06     Decimation offset:         00000",
        f"B057F07     Estimated delay (seconds): {delay}",
        f"B057F08     Correction applied (seconds): {correction}",
    ]


def _coefficients(number, numerators, denominators=(), code="D",
                  units=("COUNTS", "COUNTS")):
    lines = [
        f"B054F03     Transfer function type:    {code}",
        f"B054F04     Stage sequence number:     {number}",
        f"B054F05     Response in units lookup:  {units[0]} - In",
        f"B054F06     Response out units lookup: {units[1]} - Out",
        f"B054F07     Number of numerators:      {len(numerators)}",
        f"B054F10     Number of denominators:    {len(denominators)}",
    ]
    lines += [f"B054F08-09  {i}  {c}  0.0" for i, c in enumerate(numerators)]
    lines += [
        f"B054F11-12  {i}  {c}  0.0" for i, c in enumerate(denominators)
    ]
    return lines


def _fir(number, symmetry, coefficients):
    lines = [
        f"B061F03     Stage sequence number:     {number}",
        "B061F04     Response Name:             FIR_TEST",
        f"B061F05     Symmetry Code:             {symmetry}",
        "B061F06     Response in units lookup:  COUNTS - Digital counts",
        "B061F07     Response out units lookup: COUNTS - Digital counts",
        f"B061F08     Number of Coefficients:    {len(coefficients)}",
    ]
    lines += [f"B061F09     {i}  {c}" for i, c in enumerate(coefficients)]
    return lines


def _polynomial(number, coefficients=(0.0, 0.5), frequency_units="B"):
    lines = [
        "B062F03     Transfer function type:    P",
        f"B062F04     Stage sequence number:     {number}",
        "B062F05     Response in units lookup:  PA - Pressure in Pascals",
        "B062F06     Response out units lookup: V - Volts",
        "B062F07     Polynomial Approximation Type: M",
        f"B062F08     Valid Frequency Units:     {frequency_units}",
        "B062F09     Lower Valid Frequency Bound: 0.0",
        "B062F10     Upper Valid Frequency Bound: 6.283185307179586",
        "B062F11     Lower Bound of Approximation: 600.0",
        "B062F12     Upper Bound of Approximation: 1100.0",
        "B062F13     Maximum Absolute Error:    0.0",
        f"B062F14     Number of Coefficients:    {len(coefficients)}",
    ]
    lines += [
        f"B062F15-16  {i}  {c}  0.0" for i, c in enumerate(coefficients)
    ]
    return lines


@pytest.fixture(scope="module")
def anmo():
    """Response of IU.ANMO.00.BHZ (two epochs)."""
    return read_resp(ANMO)["IU.ANMO.00.BHZ"]


# ----------------------------------------------------------------------------
# Real files
# ----------------------------------------------------------------------------


def test_nrl_resp_equals_stationxml():
    """RESP and StationXML of the same NRL response give the same model."""
    resp = read_resp(NRL_RESP)
    sxml = read_stationxml(NRL_XML)
    assert resp.sid == sxml.sid == ["IV.MODE.00.HNZ"]

    first = resp[0].stageset[0]
    second = sxml[0].stageset[0]
    assert first.starttime == second.starttime
    assert first.endtime is None and second.endtime is None
    assert first.sample_rate == second.sample_rate == 100.0
    assert first.sensitivity.to_dict() == second.sensitivity.to_dict()
    assert [type(s) for s in first.stage] == [type(s) for s in second.stage]

    for one, two in zip(first.stage, second.stage):
        for key in one._KEYMAP:
            if key == "description":
                continue
            a, b = one[key], two[key]
            if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
                np.testing.assert_allclose(a, b, rtol=1e-12)
            else:
                assert a == b, (one.stage_number, key)

    frequency = np.logspace(-2, np.log10(45.0), 50)
    np.testing.assert_allclose(
        first.transfer_function(frequency),
        second.transfer_function(frequency),
        rtol=1e-12,
    )
    assert first.validate() == []


def test_anmo_epochs(anmo):
    """Two epochs with truncated SEED times; the second is open."""
    assert len(anmo.stageset) == 2
    first, second = anmo.stageset
    assert str(first.starttime) == "1998-10-26T20:00:00.000000Z"
    assert str(first.endtime) == "2000-10-19T16:00:00.000000Z"
    assert second.starttime == first.endtime
    assert second.endtime is None


def test_anmo_stages(anmo):
    """Zeros, gains, decimation and overall sensitivity are read."""
    stage_set = anmo.extract("2015-01-01T00:00:00")
    paz = stage_set.stage[0]
    assert isinstance(paz, StagePoleZero)
    assert paz.pz_type == "LAPLACE (RADIANS/SECOND)"
    np.testing.assert_array_equal(paz.zeros, [0j, 0j])
    assert len(paz.poles) == 5
    assert (paz.input_units, paz.output_units) == ("M/S", "V")

    gains = [s for s in stage_set.stage if isinstance(s, StageGain)]
    assert gains[0].frequency == 0.02

    assert stage_set.sensitivity.value == SENSITIVITY
    assert stage_set.sensitivity.frequency == 0.02
    assert stage_set.sensitivity.input_units == "M/S"
    assert stage_set.sensitivity.output_units == "COUNTS"
    assert stage_set.sample_rate == 20.0


def test_anmo_against_evalresp(anmo):
    """
    Amplitude as computed by evalresp; phase up to the residual delay.

    evalresp evaluates symmetric FIR filters with zero phase, i.e. it
    assumes their delay fully corrected. ShakeLab applies the declared
    corrections (SEED blockette 57, field 8), which here leave a delay
    equal to the sum of ``(N - 1) / 2 / fs - correction``.
    """
    amp = np.loadtxt(DATA / "AMP.IU.ANMO.00.BHZ")
    phase = np.loadtxt(DATA / "PHASE.IU.ANMO.00.BHZ")
    stage_set = anmo.extract("2015-01-01T00:00:00")
    response = stage_set.transfer_function(amp[:, 0])

    np.testing.assert_allclose(np.abs(response), amp[:, 1], rtol=1e-4)

    residual = sum(
        (np.size(s.numerator) - 1) / 2.0 / s.input_sample_rate
        - (s.correction or 0.0)
        for s in stage_set.stage
        if isinstance(s, StageCoefficients) and s.input_sample_rate
    )
    expected = np.deg2rad(phase[:, 1]) - 2 * np.pi * amp[:, 0] * residual
    difference = np.angle(response * np.exp(-1j * expected))
    np.testing.assert_allclose(difference, 0.0, atol=1e-3)


@pytest.mark.parametrize("name", ["RESP.IU.ANMO.00.BHZ", "NRL_MODE_HNZ.resp"])
def test_collection_infers_resp(name):
    """ResponseCollection recognizes RESP files by their name."""
    collection = ResponseCollection(str(DATA / name))
    assert len(collection) == 1


def test_collection_explicit_resp(tmp_path):
    """ftype="resp" reads files with any name."""
    path = tmp_path / "response.txt"
    path.write_text(NRL_RESP.read_text())
    collection = ResponseCollection(str(path), ftype="resp")
    assert collection.sid == ["IV.MODE.00.HNZ"]

    other = tmp_path / "mode.resp"
    other.write_text(NRL_RESP.read_text())
    assert ResponseCollection(str(other)).sid == ["IV.MODE.00.HNZ"]


# ----------------------------------------------------------------------------
# Sources
# ----------------------------------------------------------------------------


def test_sources():
    """Paths, text, bytes and open files give the same result."""
    text = NRL_RESP.read_text()
    sources = [
        NRL_RESP,
        str(NRL_RESP),
        text,
        text.encode(),
        io.StringIO(text),
        io.BytesIO(text.encode()),
    ]
    for source in sources:
        collection = read_resp(source)
        assert collection.sid == ["IV.MODE.00.HNZ"]
        assert len(collection[0].stageset[0].stage) == 8


def test_missing_file():
    """A missing path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_resp("/nonexistent/RESP.XX.NONE..HHZ")


def test_bad_source_type():
    """Unsupported sources raise TypeError."""
    with pytest.raises(TypeError):
        read_resp(3.14)


def test_no_response():
    """Text without channel responses raises ValueError."""
    with pytest.raises(ValueError):
        read_resp("# nothing\n# here\n")


# ----------------------------------------------------------------------------
# Times and header
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2019,123,04:05:06.789", "2019-05-03T04:05:06.789000Z"),
        ("2019,123,04:05:06", "2019-05-03T04:05:06.000000Z"),
        ("2019,123,04:05", "2019-05-03T04:05:00.000000Z"),
        ("1998,299,20", "1998-10-26T20:00:00.000000Z"),
        ("1999,013", "1999-01-13T00:00:00.000000Z"),
        ("2001", "2001-01-01T00:00:00.000000Z"),
    ],
)
def test_parse_resp_time(text, expected):
    """SEED times may be truncated after any subfield."""
    assert str(parse_resp_time(text)) == expected


@pytest.mark.parametrize("text", [None, "", "No Ending Time"])
def test_parse_resp_time_open(text):
    """Missing end times mean an open interval."""
    assert parse_resp_time(text) is None


@pytest.mark.parametrize("text", ["2019-05-03", "yesterday", "2019,12a"])
def test_parse_resp_time_invalid(text):
    """Times not in SEED format raise ValueError."""
    with pytest.raises(ValueError):
        parse_resp_time(text)


def test_empty_location():
    """rdseed writes "??" for an empty location code."""
    collection = read_resp(_resp(_paz(), location="??"))
    assert collection.sid == ["XX.TEST..HHZ"]


def test_epochs_and_channels():
    """Each channel header starts a new epoch; headers may be partial."""
    text = "\n".join([
        _resp(_paz(), _gain(1, 10.0), start="2020,001", end="2021,001"),
        "B052F03     Location:    00",
        "B052F04     Channel:     HHZ",
        "B052F22     Start date:  2021,001",
        "B052F23     End date:    No Ending Time",
        *_paz(), *_gain(1, 20.0),
        "B050F03     Station:     OTHER",
        "B050F16     Network:     XX",
        "B052F03     Location:    ",
        "B052F04     Channel:     HHN",
        "B052F22     Start date:  2020,001",
        *_paz(),
    ])
    collection = read_resp(text)
    assert collection.sid == [SID, "XX.OTHER..HHN"]
    first, second = collection[SID].stageset
    assert first.stage[1].sensitivity == 10.0
    assert second.stage[1].sensitivity == 20.0
    assert second.endtime is None


# ----------------------------------------------------------------------------
# Blockettes
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "code, pz_type",
    [
        ("A [Laplace Transform (Rad/sec)]", "LAPLACE (RADIANS/SECOND)"),
        ("A - Laplace transform analog response, in rad/sec",
         "LAPLACE (RADIANS/SECOND)"),
        ("B", "LAPLACE (HERTZ)"),
        ("D", "DIGITAL (Z-TRANSFORM)"),
    ],
)
def test_paz_types(code, pz_type):
    """Blockette 53: transfer function types A, B and D."""
    stage = _stages(_paz(code=code)).stage[0]
    assert isinstance(stage, StagePoleZero)
    assert stage.pz_type == pz_type
    assert stage.normalization_factor == 2.0
    np.testing.assert_array_equal(stage.zeros, [0j])
    np.testing.assert_array_equal(stage.poles, [-1 + 1j, -1 - 1j])


def test_paz_composite_unsupported():
    """Blockette 53 of composite type is kept as unsupported."""
    with pytest.warns(ResponseWarning, match="not supported"):
        stage = _stages(_paz(code="C")).stage[0]
    assert isinstance(stage, StageUnsupported)
    assert (stage.input_units, stage.output_units) == ("M/S", "V")
    with pytest.raises(NotImplementedError):
        stage.transfer_function([1.0])


def test_paz_without_poles():
    """Poles and zeros may be absent: the stage is a pure gain."""
    stage = _stages(_paz(poles=(), zeros=())).stage[0]
    np.testing.assert_allclose(stage.transfer_function([0.1, 1.0]), 2.0)


def test_gain_frequency_with_units():
    """Blockette 58: frequencies may carry the "HZ" suffix."""
    stage_set = _stages(_paz(), _gain(1, 1500.0, "2.000000E-02 HZ"))
    gain = stage_set.stage[1]
    assert isinstance(gain, StageGain)
    assert gain.stage_number == 1
    assert (gain.sensitivity, gain.frequency) == (1500.0, 0.02)


def test_overall_sensitivity():
    """Blockette 58 of stage 0 is the overall sensitivity."""
    stage_set = _stages(
        _paz(), _gain(1, 1500.0),
        _coefficients(2, [], units=("V", "COUNTS")),
        _decimation(2, 100.0), _gain(2, 1000.0, "0.0"),
        _gain(0, 1.5e6, "1.0"),
    )
    assert len(stage_set.stage) == 4
    sens = stage_set.sensitivity
    assert (sens.value, sens.frequency) == (1.5e6, 1.0)
    assert (sens.input_units, sens.output_units) == ("M/S", "COUNTS")


def test_coefficients_digital():
    """Blockette 54: digital filter with decimation and correction."""
    stage_set = _stages(
        _coefficients(1, [0.25, 0.5, 0.25]),
        _decimation(1, 200.0, factor=2, delay=0.005, correction=0.005),
    )
    stage = stage_set.stage[0]
    assert isinstance(stage, StageCoefficients)
    assert stage.cf_type == "DIGITAL"
    np.testing.assert_array_equal(stage.numerator, [0.25, 0.5, 0.25])
    np.testing.assert_array_equal(stage.denominator, [1.0])
    assert stage.input_sample_rate == 200.0
    assert stage.decimation_factor == 2
    assert stage.decimation_offset == 0
    assert stage.delay == stage.correction == 0.005
    assert stage_set.sample_rate == 100.0
    # Symmetric filter with corrected delay: zero phase.
    response = stage.transfer_function([10.0, 40.0])
    np.testing.assert_allclose(np.angle(response), 0.0, atol=1e-12)


def test_coefficients_iir_and_analog():
    """Blockette 54: denominators and analog types."""
    stage = _stages(
        _coefficients(1, [1.0], [1.0, -0.5]), _decimation(1, 100.0)
    ).stage[0]
    np.testing.assert_array_equal(stage.denominator, [1.0, -0.5])
    np.testing.assert_allclose(stage.transfer_function([0.0]), 2.0)

    stage = _stages(_coefficients(1, [0.0, 1.0], code="B")).stage[0]
    assert stage.cf_type == "ANALOG (HERTZ)"
    np.testing.assert_allclose(stage.transfer_function([2.0]), 2.0j)


def test_coefficients_without_numerators():
    """An analog-to-digital converter has a unit response."""
    stage = _stages(
        _coefficients(1, [], units=("V", "COUNTS")), _decimation(1, 100.0)
    ).stage[0]
    np.testing.assert_array_equal(stage.numerator, [1.0])
    assert (stage.input_units, stage.output_units) == ("V", "COUNTS")


def test_coefficients_continuation():
    """A blockette continued in a second one of the same stage."""
    first = _coefficients(1, [0.1, 0.2])
    second = _coefficients(1, [0.3])
    stage_set = _stages(first, second, _decimation(1, 100.0))
    assert len(stage_set.stage) == 1
    np.testing.assert_array_equal(stage_set.stage[0].numerator,
                                  [0.1, 0.2, 0.3])


@pytest.mark.parametrize(
    "symmetry, expected",
    [
        ("A", [0.1, 0.4, 0.5]),
        ("B", [0.1, 0.4, 0.5, 0.4, 0.1]),
        ("C", [0.1, 0.4, 0.5, 0.5, 0.4, 0.1]),
    ],
)
def test_fir_symmetry(symmetry, expected):
    """Blockette 61: symmetry codes A (none), B (odd), C (even)."""
    stage = _stages(_fir(1, symmetry, [0.1, 0.4, 0.5])).stage[0]
    assert isinstance(stage, StageCoefficients)
    assert stage.description == "FIR_TEST"
    np.testing.assert_allclose(stage.numerator, expected)


def test_fir_unknown_symmetry():
    """An unknown symmetry code raises ValueError."""
    with pytest.raises(ValueError):
        _stages(_fir(1, "X", [1.0]))


@pytest.mark.parametrize("with_index", [True, False])
def test_response_list(with_index):
    """Blockette 55, with or without the row index."""
    rows = [(1.0, 2.0, 0.0, -90.0, 0.0), (10.0, 4.0, 0.0, -180.0, 0.0)]
    lines = [
        "B055F03     Stage sequence number:     1",
        "B055F04     Response in units lookup:  NM - Displacement",
        "B055F05     Response out units lookup: COUNTS - Digital counts",
        "B055F06     Number of responses:       2",
    ] + [
        "B055F07-11  " + (f"{i}  " if with_index else "")
        + "  ".join(str(v) for v in row)
        for i, row in enumerate(rows)
    ]
    stage = _stages(lines).stage[0]
    assert isinstance(stage, StageResponseList)
    np.testing.assert_array_equal(stage.frequency, [1.0, 10.0])
    np.testing.assert_array_equal(stage.amplitude, [2.0, 4.0])
    np.testing.assert_array_equal(stage.phase, [-90.0, -180.0])
    assert stage.input_units == "NM"


def test_polynomial_stage():
    """Blockette 62 in the cascade; frequency bounds in rad/s."""
    stage = _stages(_polynomial(1, frequency_units="A")).stage[0]
    assert isinstance(stage, StagePolynomial)
    assert stage.approximation_type == "MACLAURIN"
    np.testing.assert_array_equal(stage.coefficients, [0.0, 0.5])
    assert stage.frequency_upper_bound == pytest.approx(1.0)
    assert stage.approximation_lower_bound == 600.0
    assert (stage.input_units, stage.output_units) == ("PA", "V")


def test_instrument_polynomial():
    """Blockette 62 of stage 0 is the overall polynomial."""
    stage_set = _stages(_polynomial(1), _polynomial(0, (600.0, 2.0)))
    assert len(stage_set.stage) == 1
    np.testing.assert_array_equal(
        stage_set.instrument_polynomial.coefficients, [600.0, 2.0]
    )
    assert stage_set.sample_rate is None


@pytest.mark.parametrize("btype", ["056", "060"])
def test_unsupported_blockettes(btype):
    """Generic responses and response references are unsupported."""
    number_field = "F04" if btype == "060" else "F03"
    lines = [
        f"B{btype}{number_field}     Stage sequence number:     1",
    ]
    if btype == "056":
        lines += [
            "B056F04     Response in units lookup:  M/S - Velocity",
            "B056F05     Response out units lookup: V - Volts",
        ]
    with pytest.warns(ResponseWarning):
        stage = _stages(lines).stage[0]
    assert isinstance(stage, StageUnsupported)
    assert btype in stage.element


def test_decimation_without_filter():
    """A decimation without response blockette is a unit filter."""
    stage_set = _stages(_paz(), _decimation(2, 100.0, factor=2))
    stage = stage_set.stage[1]
    assert isinstance(stage, StageCoefficients)
    assert stage.stage_number == 2
    assert stage_set.sample_rate == 50.0


def test_not_a_number():
    """Malformed numeric fields raise ValueError."""
    with pytest.raises(ValueError):
        _stages(_gain(1, "abc"))
    with pytest.raises(ValueError):
        _stages(_paz() + ["B053F15-18  2  x  0.0  0.0  0.0"])


def test_warnings_free_on_real_files():
    """The reference files are read without warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_resp(NRL_RESP)
        read_resp(ANMO)
