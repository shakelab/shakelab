# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Parsers for StationXML files to extract instrument response information.

This module provides functions to parse FDSN-compliant StationXML files and
convert the instrument response metadata into a ShakeLab ResponseCollection
object. Supported response types include poles and zeros (PAZ), FIR filters,
polynomials, and gain stages.
"""
import numpy as np
import xml.etree.ElementTree as ET
import io, os, re
from pathlib import Path

import shakelab.signals.response as rspm
from shakelab.libutils.utils import cast_value


PZTYPEMAP = ["LAPLACE (RADIANS/SECOND)",
             "LAPLACE (HERTZ)",
             "DIGITAL (Z-TRANSFORM)"]

UNITSMAP = [""]


def parse_sxml(xml):
    """
    Parse a StationXML string or file and return a ResponseCollection.

    This function reads a StationXML document (either from a file path
    or as a string) and extracts network, station, and channel
    response information.
    It builds a  ResponseCollection composed of StreamResponse objects,
    each containing one or more StageSet entries representing the instrument 
    response chain.

    Args:
        xml (str or Path): Path to a StationXML file or an XML string.

    Returns:
        ResponseCollection: A ShakeLab ResponseCollection object
        containing parsed instrument responses.
    """
    if isinstance(xml, (str, Path)):
        if not xml.lstrip().startswith("<?xml"):
            path = Path(xml)
            if path.is_file():
                with open(path, "r", encoding="utf-8") as f:
                    xml = f.read()
            else:
                raise ValueError(f"File not found: {path}")
    else:
        raise TypeError("Input must be a StationXML string or a file path.")

    xml = xml_strip(xml)
    it = ET.iterparse(io.StringIO(xml))

    # Strip namespaces for cleaner tag access
    for _, el in it:
        el.tag = re.sub(r'\{.*\}', '', el.tag)

    root = it.root

    if 'FDSNStationXML' not in root.tag:
        raise ValueError('Not a valid FDSN StationXML document')

    rcoll = rspm.ResponseCollection()

    for network in root.findall('.//Network'):
        net = network.attrib.get('code')

        for station in network.findall('.//Station'):
            sta = station.attrib.get('code')

            for channel in station.findall('.//Channel'):
                cha = channel.attrib.get('code')
                loc = channel.attrib.get('locationCode')
                start = channel.attrib.get('startDate')
                end = channel.attrib.get('endDate')

                sid = f"{net}.{sta}.{loc}.{cha}"

                if sid not in rcoll:
                    rcoll.append(rspm.StreamResponse(sid))

                stage_set = rspm.StageSet(start, end)
                stage_set.append(parse_response(channel))
                rcoll[sid].append(stage_set)

    return rcoll


def parse_response(channel):
    """
    Parse all response stages of a Channel element and return a list 
    of Stage objects.

    This function inspects each <Stage> element under the given 
    <Channel> element in the StationXML structure. It parses recognized 
    stage types (StageGain, PolesZeros, FIR, Polynomial), converting 
    them into corresponding ShakeLab Stage objects. Unsupported tags 
    such as Decimation or ResponseList are ignored.

    Args:
        channel (xml.etree.ElementTree.Element): The <Channel> XML 
            element containing one or more <Stage> definitions.

    Returns:
        list: A list of Stage objects, each representing a component 
        of the instrument response chain for the given channel.
    """
    stage_list = []
    response = channel.find('Response')

    for stage in response.findall('Stage'):
        stage_number = int(stage.get('number'))

        parsed = None

        for child in stage:
            tag = child.tag.split('}')[-1]  # Handle possible namespace

            if tag == 'StageGain':
                parsed = parse_gain(child)

            elif tag == 'PolesZeros':
                parsed = parse_paz(child)

            elif tag == 'Coefficients':
                #parsed = parse_coefficients(child)
                pass

            elif tag == 'FIR':
                parsed = parse_fir(child)

            elif tag == 'Polynomial':
                parsed = parse_polynomial(child)

            elif tag in ['Decimation', 'ResponseList']:
                continue  # Not handled as separate Stage objects

            else:
                print(f'Unrecognized stage element: {tag}')

            if parsed:
                parsed.stage_number = stage_number
                stage_list.append(parsed)
                parsed = None  # Reset for safety

    return stage_list


def parse_gain(element):
    """
    Parse a <StageGain> XML element and return a StageGain object.

    Extracts the gain value, frequency, optional description, and 
    stage number from the provided XML element, converting them into 
    a ShakeLab StageGain object used to model amplification stages 
    within an instrument response.

    Args:
        element (xml.etree.ElementTree.Element): The <StageGain> XML 
            element to be parsed.

    Returns:
        StageGain: A ShakeLab StageGain object containing gain 
        information for the corresponding stage.
    """
    data = {
        'description': None,
        'sensitivity': None,
        'frequency': None,
        'stage_number': None
    }

    data['sensitivity'] = cast_value(element.findtext('Value'), float)
    data['frequency'] = cast_value(element.findtext('Frequency'), float)
    data['description'] = element.get('description')
    data['stage_number'] = cast_value(element.get('stage_number'), int)

    return rspm.StageGain(data)


def parse_paz(element):
    """
    Parse a <PolesZeros> XML element and return a StagePoleZero object.

    Extracts poles, zeros, normalization factor and frequency, units,
    and transfer function type from the given StationXML element, and 
    returns a ShakeLab StagePoleZero object representing a stage of the 
    instrument response in the Laplace domain.

    Args:
        element (xml.etree.ElementTree.Element): The <PolesZeros> XML 
            element to be parsed.

    Returns:
        StagePoleZero: A ShakeLab StagePoleZero object representing the
        parsed response stage.
    """
    data = {
        'description': None,
        'input_units': None,
        'output_units': None,
        'pz_type': None,
        'normalization_factor': None,
        'normalization_frequency': None,
        'poles': [],
        'zeros': [],
        'stage_number': None
    }

    # Input/output units
    data['input_units'] = element.findtext('./InputUnits/Name')
    data['output_units'] = element.findtext('./OutputUnits/Name')
    data['pz_type'] = element.findtext('PzTransferFunctionType')

    # Normalization
    data['normalization_factor'] = cast_value(
        element.findtext('NormalizationFactor'), float
    )
    data['normalization_frequency'] = cast_value(
        element.findtext('NormalizationFrequency'), float
    )

    # Poles
    for pole in element.findall('.//Pole'):
        re = cast_value(pole.findtext('Real'), float)
        im = cast_value(pole.findtext('Imaginary'), float)
        data['poles'].append(complex(re, im))

    # Zeros
    for zero in element.findall('.//Zero'):
        re = cast_value(zero.findtext('Real'), float)
        im = cast_value(zero.findtext('Imaginary'), float)
        data['zeros'].append(complex(re, im))

    return rspm.StagePoleZero(data)


def parse_coefficients(element):
    """
    Parse a <Coefficients> XML element and return a StageFIR object.

    This function extracts digital filter coefficients from a StationXML
    <Coefficients> element. It parses the numerator and denominator arrays,
    as well as input/output units and the stage number. If the transfer
    function type is not 'DIGITAL', it raises an error.

    Args:
        element (xml.etree.ElementTree.Element): The <Coefficients> XML
            element to parse.

    Returns:
        StageFIR: A ShakeLab FIR response stage with the parsed coefficients.

    Raises:
        [TO IMPROVE!!]
    """
    tf_type = element.findtext('CfTransferFunctionType')

    # Common fields
    input_units = element.findtext('InputUnits/Name')
    output_units = element.findtext('OutputUnits/Name')

    # Numerator
    numerator = [
        cast_value(num.text, float)
        for num in element.findall('.//Numerator')
        if num.text is not None
    ]

    # Denominator (may be optional)
    denominator = [
        cast_value(den.text, float)
        for den in element.findall('.//Denominator')
        if den.text is not None
    ]

    if not numerator:
        numerator = [1.0]

    if not denominator:
        denominator = [1.0]

    if tf_type == 'DIGITAL':
        data = {
            'input_units': input_units,
            'output_units': output_units,
            'numerator': numerator,
            'denominator': denominator,
            'stage_number': cast_value(element.get('stage_number'), int)
        }
        return rspm.StageFIR(data)

    raise NotImplementedError(f"TF type {tf_type} not supported.")


def parse_polynomial(element):
    """
    Parse a <Polynomial> XML element and return a StagePolynomial object.

    This function extracts a polynomial response stage from a StationXML
    <Polynomial> element. It includes input/output units, numerator
    coefficients, and stage number. Denominator is always set to [1.0],
    as only feedforward coefficients are supported.

    Args:
        element (xml.etree.ElementTree.Element): The <Polynomial> XML
            element to parse.

    Returns:
        StagePolynomial: A ShakeLab polynomial response stage.
    """
    data = {
        'description': element.get('description'),
        'input_units': element.findtext('InputUnits/Name'),
        'output_units': element.findtext('OutputUnits/Name'),
        'numerator': [],
        'denominator': [1.0],
        'stage_number': cast_value(element.get('stage_number'), int)
    }

    for coeff in element.findall('Coefficient'):
        val = cast_value(coeff.text, float)
        data['numerator'].append(val)

    if not data['numerator']:
        data['numerator'] = [1.0]

    return rspm.StagePolynomial(data)


def parse_fir(element):
    """
    Parse an <FIR> XML element and return a StageFIR object.

    This function reads an FIR (Finite Impulse Response) response stage from
    a StationXML <FIR> element. It extracts the input/output units, stage
    number, and numerator coefficients. The denominator is always set to [1.0],
    as FIR filters are non-recursive.

    Args:
        element (xml.etree.ElementTree.Element): The <FIR> XML element to parse.

    Returns:
        StageFIR: A ShakeLab FIR filter response stage.
    """
    data = {
        'input_units': element.findtext('InputUnits/Name'),
        'output_units': element.findtext('OutputUnits/Name'),
        'numerator': [],
        'denominator': [1.0],  # Default for FIR
        'stage_number': cast_value(element.get('stage_number'), int)
    }

    for coeff in element.findall('.//Numerator'):
        if coeff.text is not None:
            data['numerator'].append(cast_value(coeff.text, float))

    if not data['numerator']:
        data['numerator'] = [1.0]

    return rspm.StageFIR(data)


def stationxml_to_dict(xml):
    """
    Parse StationXML file and return its contents as a dictionary.

    Args:
        xml (str): Path to StationXML file or XML content as string.

    Returns:
        dict: Dictionary representing StationXML data.
    """
    if os.path.isfile(xml):
        with open(xml, 'r') as f:
            xml = f.read()

    xml = xml_strip(xml)

    it = ET.iterparse(io.StringIO(xml))

    # Strip namespaces
    for _, el in it:
        el.tag = re.sub(r'\{.*\}', '', el.tag)

    return {it.root.tag : node_to_dict(it.root)}


def node_to_dict(node):
    """
    Recursively convert an ElementTree node to a dictionary.

    Args:
        node (ElementTree.Element): XML element to convert.

    Returns:
        dict: Dictionary representing the XML element.
    """
    xml_dict = {'attrib' : None, 'value' : None}

    xml_dict['attrib'] = node.attrib if node.attrib else None

    if node.text:
        xml_dict['value'] = _convert(node.text)
    else:
        xml_dict['value'] = []
        for elem in node:
            xml_dict['value'] += [{elem.tag : node_to_dict(elem)}]
            if elem.tail:
                print(elem.tail)

    return xml_dict


def xml_strip(xml):
    """
    Remove white spaces and new lines from an XML string.

    Args:
        xml (str): XML content as a string.

    Returns:
        str: XML content with white spaces and new lines removed.
    """
    lines = xml.split('\n')
    buffer = []
    for x in lines:
        y = x.strip()
        if y: buffer.append(y)
    buffer = ''.join(buffer)
    return buffer


def _convert(string):
    """
    Convert a string to a float if possible, otherwise return the string.

    Args:
        string (str): String to convert.

    Returns:
        float or str: Converted float or the original string.
    """
    try:
        return float(string)
    except:
        return string

