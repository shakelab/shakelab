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
Module to parse catalogues in quakeml format.
"""
import xml.etree.ElementTree as ET
import io, os, re

import shakelab.seismicity.catalogue as cat

def read_quakeml(xml):
    """
    xml can be a file or string (e.g. from webservice)
    """
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

    if 'quakeml' not in root.tag:
        raise ValueError('Not a valid QuakeML')

    edb = cat.EqDatabase()

    for event in root.findall('.//event'):

        # Find event ID
        event_id = _getattrib(event, 'publicID')
        print(event_id)

        preferred_magnitude_id = event.find('.//preferredMagnitudeID')
        if preferred_magnitude_id is not None:
            preferred_magnitude_id = preferred_magnitude_id.text

        parse_magnitude(event)
        parse_location(event)

    return event

def parse_magnitude(event):
    """
    """
    for magnitude in event.findall('.//magnitude'):

        mtype = magnitude.find('.//type').text
        value = magnitude.find('.//mag').find('.//value').text
        try:
            agency = magnitude.find('.//creationInfo').find('.//agencyID').text
        except:
            agency = None
        print(mtype, value, agency)

def parse_location(event):
    """
    """
    for origin in event.findall('.//origin'):

        latitude = origin.find('.//latitude').find('.//value').text
        longitude = origin.find('.//longitude').find('.//value').text
        try:
            agency = magnitude.find('.//creationInfo').find('.//agencyID').text
        except:
            agency = None
        print(longitude, latitude, agency)

def _getattrib(element, key):
    """
    """
    if key in element.attrib:
        return element.attrib[key]
    else:
        return None

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
