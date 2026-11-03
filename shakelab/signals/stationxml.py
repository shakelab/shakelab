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
Reader of FDSN StationXML instrument responses.

The reader converts the ``<Response>`` of every channel epoch of a
StationXML document into the response model of
:mod:`shakelab.signals.response`.

Mapping of the response elements
--------------------------------
=====================================  =================================
StationXML                             ShakeLab
=====================================  =================================
``InstrumentSensitivity``              ``StageSet.sensitivity``
``InstrumentPolynomial``               ``StageSet.instrument_polynomial``
``PolesZeros`` (any type)              ``StagePoleZero``
``Coefficients`` (any type)            ``StageCoefficients``
``FIR`` (symmetry expanded)            ``StageCoefficients``
``ResponseList``                       ``StageResponseList``
``Polynomial``                         ``StagePolynomial``
``Decimation``                         attributes of the digital stage
``StageGain``                          ``StageGain``
=====================================  =================================

Design choices
--------------
- The metadata are read as written: units are not rewritten and
  inconsistent values are not corrected (see
  :meth:`shakelab.signals.response.StageSet.validate`).
- Elements of unknown type are kept as ``StageUnsupported``, so that
  an incomplete response cannot be used without notice.
- Unknown elements and channels without response are reported with
  :class:`shakelab.signals.response.ResponseWarning`.
- :func:`iter_channels` walks the channel epochs and can be reused by
  readers of other channel metadata.

Notes
-----
FIR symmetry follows the FDSN definitions (see
:func:`shakelab.signals.response.expand_fir_symmetry`): with ``ODD``
symmetry the total number of coefficients is odd and the central
coefficient is stored once; with ``EVEN`` symmetry it is even and every
stored coefficient is mirrored.
"""

from __future__ import annotations

import warnings
from typing import Any, Iterator, Optional
from xml.etree.ElementTree import Element

from shakelab.libutils.xml import find_float, find_text, read_xml, to_dict
from shakelab.signals import response as rspm

__all__ = [
    "read_stationxml",
    "iter_channels",
    "parse_response",
    "stationxml_to_dict",
]

_CF_TYPES = ("DIGITAL", "ANALOG (RADIANS/SECOND)", "ANALOG (HERTZ)")

_TRANSFER_ELEMENTS = (
    "PolesZeros",
    "Coefficients",
    "ResponseList",
    "FIR",
    "Polynomial",
)


# ----------------------------------------------------------------------------
# Document
# ----------------------------------------------------------------------------


def read_stationxml(source: Any) -> rspm.ResponseCollection:
    """
    Read the instrument responses of a StationXML document.

    Parameters
    ----------
    source
        File path, XML text, bytes or open file (see
        :func:`shakelab.libutils.xml.read_xml`).

    Returns
    -------
    ResponseCollection
        One stream response per channel, one stage set per epoch.

    Raises
    ------
    ValueError
        If the document is not FDSN StationXML.
    """
    root = read_xml(source)
    if root.tag != "FDSNStationXML":
        raise ValueError("Not a valid FDSN StationXML document")

    collection = rspm.ResponseCollection()

    for sid, channel in iter_channels(root):
        element = channel.find("Response")
        if element is None:
            warnings.warn(
                f"{sid}: channel epoch without response, skipped",
                rspm.ResponseWarning,
                stacklevel=2,
            )
            continue

        stage_set = parse_response(
            element,
            starttime=channel.get("startDate"),
            endtime=channel.get("endDate"),
            sample_rate=find_float(channel, "SampleRate"),
        )
        stream_response = rspm.StreamResponse(sid)
        stream_response.append(stage_set)
        collection.append(stream_response)

    return collection


def iter_channels(root: Element) -> Iterator[tuple]:
    """
    Iterate over the channel epochs of a StationXML document.

    Parameters
    ----------
    root
        Root element (``FDSNStationXML``), without namespaces.

    Yields
    ------
    tuple
        ``(sid, channel)``: source identifier
        ``NET.STA.LOC.CHA`` and ``<Channel>`` element.
    """
    for network in root.iter("Network"):
        net = network.get("code", "")
        for station in network.iter("Station"):
            sta = station.get("code", "")
            for channel in station.iter("Channel"):
                loc = channel.get("locationCode", "")
                cha = channel.get("code", "")
                yield f"{net}.{sta}.{loc}.{cha}", channel


def stationxml_to_dict(source: Any) -> dict:
    """
    Return a StationXML document as nested dictionaries.

    See :func:`shakelab.libutils.xml.to_dict`; values are kept as
    strings.
    """
    return to_dict(read_xml(source))


# ----------------------------------------------------------------------------
# Response
# ----------------------------------------------------------------------------


def parse_response(
    element: Element,
    starttime: Any = None,
    endtime: Any = None,
    sample_rate: Optional[float] = None,
) -> rspm.StageSet:
    """
    Convert a ``<Response>`` element into a stage set.

    Parameters
    ----------
    element
        ``<Response>`` element.
    starttime, endtime
        Validity interval of the epoch.
    sample_rate
        Sampling rate of the channel in Hz.

    Returns
    -------
    StageSet
        Stages in order, with the overall sensitivity or polynomial.
    """
    stage_set = rspm.StageSet(starttime, endtime, sample_rate=sample_rate)

    sens = element.find("InstrumentSensitivity")
    if sens is not None:
        stage_set.sensitivity = rspm.Sensitivity(
            value=find_float(sens, "Value"),
            frequency=find_float(sens, "Frequency"),
            input_units=find_text(sens, "InputUnits/Name"),
            output_units=find_text(sens, "OutputUnits/Name"),
        )

    polynomial = element.find("InstrumentPolynomial")
    if polynomial is not None:
        stage_set.instrument_polynomial = _parse_polynomial(polynomial)

    for stage in element.findall("Stage"):
        stage_set.append(_parse_stage(stage))

    return stage_set


def _parse_stage(stage: Element) -> list:
    """Return the stages (transfer element and gain) of a ``<Stage>``."""
    number = int(stage.get("number"))
    parsed = []
    decimation = None
    transfer = None

    for child in stage:
        tag = child.tag
        if tag in _TRANSFER_ELEMENTS:
            transfer = _PARSERS[tag](child)
        elif tag == "Decimation":
            decimation = child
        elif tag == "StageGain":
            gain = rspm.StageGain(
                sensitivity=find_float(child, "Value"),
                frequency=find_float(child, "Frequency"),
            )
            gain.stage_number = number
            parsed.append(gain)
        else:
            warnings.warn(
                f"stage {number}: unrecognized element {tag!r} ignored",
                rspm.ResponseWarning,
                stacklevel=3,
            )

    if transfer is not None:
        transfer.stage_number = number
        if decimation is not None and "input_sample_rate" in (
            transfer._KEYMAP
        ):
            _set_decimation(transfer, decimation)
        if isinstance(transfer, rspm.StageUnsupported):
            warnings.warn(
                f"stage {number}: {transfer.element} is not supported; "
                "the response cannot be evaluated completely",
                rspm.ResponseWarning,
                stacklevel=3,
            )
        parsed.insert(0, transfer)

    return parsed


def _set_decimation(stage: rspm.StageResponse, element: Element) -> None:
    """Copy the ``<Decimation>`` parameters to a digital stage."""
    stage.input_sample_rate = find_float(element, "InputSampleRate")
    factor = find_text(element, "Factor")
    offset = find_text(element, "Offset")
    stage.decimation_factor = None if factor is None else int(factor)
    stage.decimation_offset = None if offset is None else int(offset)
    stage.delay = find_float(element, "Delay")
    stage.correction = find_float(element, "Correction")


def _units(element: Element) -> dict:
    """Return the input and output units of a transfer element."""
    return {
        "description": element.get("description")
        or find_text(element, "Description"),
        "input_units": find_text(element, "InputUnits/Name"),
        "output_units": find_text(element, "OutputUnits/Name"),
    }


def _values(element: Element, tag: str) -> list:
    """Return the numeric values of the sub-elements ``tag``."""
    return [
        float(item.text) for item in element.findall(tag)
        if item.text and item.text.strip()
    ]


def _complex_values(element: Element, tag: str) -> list:
    """Return the complex values of the sub-elements ``tag``."""
    return [
        complex(find_float(item, "Real", 0.0),
                find_float(item, "Imaginary", 0.0))
        for item in element.findall(tag)
    ]


def _parse_paz(element: Element) -> rspm.StagePoleZero:
    """Convert a ``<PolesZeros>`` element."""
    return rspm.StagePoleZero(
        pz_type=find_text(element, "PzTransferFunctionType"),
        normalization_factor=find_float(element, "NormalizationFactor"),
        normalization_frequency=find_float(
            element, "NormalizationFrequency"
        ),
        poles=_complex_values(element, "Pole"),
        zeros=_complex_values(element, "Zero"),
        **_units(element),
    )


def _parse_coefficients(element: Element) -> rspm.StageResponse:
    """
    Convert a ``<Coefficients>`` element.

    Without numerator coefficients the stage has a unit response (e.g.
    an analog-to-digital converter). Unknown types give an unsupported
    stage.
    """
    cf_type = (find_text(element, "CfTransferFunctionType") or "").upper()
    units = _units(element)

    if cf_type not in _CF_TYPES:
        return rspm.StageUnsupported(
            element=f"Coefficients ({cf_type or 'unknown type'})", **units
        )

    return rspm.StageCoefficients(
        cf_type=cf_type,
        numerator=_values(element, "Numerator") or [1.0],
        denominator=_values(element, "Denominator") or [1.0],
        **units,
    )


def _parse_fir(element: Element) -> rspm.StageCoefficients:
    """
    Convert a ``<FIR>`` element, expanding symmetric coefficients.

    Raises
    ------
    ValueError
        If the symmetry is unknown.
    """
    items = element.findall("NumeratorCoefficient")
    if all(item.get("i") is not None for item in items):
        items = sorted(items, key=lambda item: int(item.get("i")))
    coefficients = [
        float(item.text) for item in items
        if item.text and item.text.strip()
    ]

    coefficients = rspm.expand_fir_symmetry(
        coefficients, find_text(element, "Symmetry") or "NONE"
    )

    return rspm.StageCoefficients(
        cf_type="DIGITAL",
        numerator=coefficients or [1.0],
        denominator=[1.0],
        **_units(element),
    )


def _parse_response_list(element: Element) -> rspm.StageResponseList:
    """Convert a ``<ResponseList>`` element."""
    items = element.findall("ResponseListElement")
    return rspm.StageResponseList(
        frequency=[find_float(item, "Frequency") for item in items],
        amplitude=[find_float(item, "Amplitude") for item in items],
        phase=[find_float(item, "Phase", 0.0) for item in items],
        **_units(element),
    )


def _parse_polynomial(element: Element) -> rspm.StagePolynomial:
    """Convert a ``<Polynomial>`` or ``<InstrumentPolynomial>`` element."""
    items = element.findall("Coefficient")
    if all(item.get("number") is not None for item in items):
        items = sorted(items, key=lambda item: int(item.get("number")))
    return rspm.StagePolynomial(
        approximation_type=find_text(
            element, "ApproximationType", "MACLAURIN"
        ),
        frequency_lower_bound=find_float(element, "FrequencyLowerBound"),
        frequency_upper_bound=find_float(element, "FrequencyUpperBound"),
        approximation_lower_bound=find_float(
            element, "ApproximationLowerBound"
        ),
        approximation_upper_bound=find_float(
            element, "ApproximationUpperBound"
        ),
        maximum_error=find_float(element, "MaximumError"),
        coefficients=[float(item.text) for item in items],
        **_units(element),
    )


_PARSERS = {
    "PolesZeros": _parse_paz,
    "Coefficients": _parse_coefficients,
    "ResponseList": _parse_response_list,
    "FIR": _parse_fir,
    "Polynomial": _parse_polynomial,
}
