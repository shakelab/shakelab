# ****************************************************************************
#
# Copyright (C) 2019-2023, ShakeLab Developers.
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

import shakelab.signals.response as rspm


def parse_sxml(xml):

    if os.path.isfile(xml):
        with open(xml, 'r') as f:
            xml = f.read()

    # Preliminarily clean xml
    xml = xml_strip(xml)

    it = ET.iterparse(io.StringIO(xml))

    # Strip namespaces
    for _, el in it:
        el.tag = re.sub(r'\{.*\}', '', el.tag)

    root = it.root

    if 'FDSNStationXML'not in root.tag:
        raise ValueError('Not a valid FDSN StationXML')

    rcoll = rspm.ResponseCollection()

    for network in root.findall('.//Network'):
        network_code = _getattrib(network, 'code')

        for station in network.findall('.//Station'):
            station_code = _getattrib(station, 'code')

            for channel in station.findall('.//Channel'):
                channel_code = _getattrib(channel, 'code')
                location_code = _getattrib(channel, 'locationCode')
                starttime = _getattrib(channel, 'startDate')
                endtime = _getattrib(channel, 'endDate')

                fdsn_code = '{0}.{1}.{2}.{3}'.format(network_code,
                                                     station_code,
                                                     location_code,
                                                     channel_code)

                if fdsn_code not in rcoll.sid:
                    strmr = rspm.StreamResponse(fdsn_code)
                    rcoll.append(strmr)

                srec = rspm.StageRecord(starttime, endtime)
                srec.append(parse_response(channel))

                rcoll[fdsn_code].append(srec)

    return rcoll

def _getattrib(element, key):
    """
    """
    if key in element.attrib:
        return element.attrib[key]
    else:
        return None

def parse_response(channel):
    """
    """
    stage_list = []
    response = channel.find('Response')

    for stage in response.findall(".//Stage"):
        stage_number = int(stage.attrib['number'])

        for child in stage:
            if child.tag == 'StageGain':
                stage_list.append(parse_gain(child))

            if child.tag == 'PolesZeros':
                stage_list.append(parse_polezero(child))

            stage_list[-1].stage_number = stage_number

    return stage_list

def parse_gain(element):
    """
    """
    value = element.find('Value').text
    frequency = element.find('Frequency').text

    stage = rspm.StageGain()
    stage.sensitivity = float(value)

    return stage

def parse_polezero(element):
    """
    """
    input_units = element.find('InputUnits').find('Name').text
    output_units = element.find('InputUnits').find('Name').text
    normalization_factor = element.find('NormalizationFactor').text
    normalization_frequency = element.find('NormalizationFrequency').text

    zeros = []
    for zero in element.findall('.//Zero'):
        real = float(zero.find('Real').text)
        imag = float(zero.find('Imaginary').text)
        zeros.append(real + imag*1j)

    poles = []
    for pole in element.findall('.//Pole'):
        real = float(pole.find('Real').text)
        imag = float(pole.find('Imaginary').text)
        poles.append(real + imag*1j)

    stage = rspm.StagePoleZero()
    stage.input_units = input_units
    stage.output_units = output_units
    stage.normalization_factor = float(normalization_factor)
    stage.normalization_frequency = float(normalization_frequency)
    stage.poles = np.array(poles)
    stage.zeros = np.array(zeros)

    return stage

def parse_coefficients():
    """
    """
    pass

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

