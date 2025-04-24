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
    """
    is_str_or_path = isinstance(xml, (str, Path))
    is_xml = xml.lstrip().startswith('<?xml')
    
    if is_str_or_path and not is_xml:
        path = Path(xml)
        if path.is_file():
            with open(path, 'r') as f:
                xml = f.read()

    xml = xml_strip(xml)
    it = ET.iterparse(io.StringIO(xml))

    # Remove namespaces
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
    Parse all response stages of a channel and return a list of Stage objects.
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
                #parsed = parse_fir(child)
                pass

            elif tag == 'Polynomial':
                #parsed = parse_polynomial(child)
                pass

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
    Parse a StageGain from XML and return a StageGain object.
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
    Parse a PAZ stage from StationXML and return a StagePoleZero object.
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
    Parse a digital coefficient stage and return a StageFIR object.
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
    Parse a Polynomial stage from StationXML and return a StagePolynomial
    object.
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
    Parse an FIR filter stage from StationXML and return a StageFIR object.
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

