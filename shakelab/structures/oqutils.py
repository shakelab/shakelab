# ****************************************************************************
#
# Copyright (C) 2019-2020, ShakeLab Developers.
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
from shakelab.structures.fragility import (FragilityModelParametric,
                                           FragilityModelDiscrete)
import xml.etree.cElementTree as xet


def fragility_to_xml(fragility_collection, xml_file, ndl=0.1):
    """
    """
    nrml = xet.Element('nrml', {
                'xmlns': 'http://openquake.org/xmlns/nrml/0.5',
                'xmlns:gml': 'http://www.opengis.net/gml'})

    fm = xet.SubElement(nrml, 'fragilityModel', {
                'assetCategory': 'buildings',
                'id': 'buildings',
                'lossCategory': 'structural'})

    xet.SubElement(fm, 'description').text = 'Fragility model'
    xet.SubElement(fm, 'limitStates').text = 'D1 D2 D3 D4 D5'

    for m in fragility_collection.model:
        if isinstance(m, FragilityModelParametric):
            ff = xet.SubElement(fm, 'fragilityFunction', {
                        'id': '{0}'.format(m.id),
                        'format': 'continuous',
                        'shape': 'lognormal'})

            im = xet.SubElement(ff, 'imls', {
                        'imt': '{0}'.format(m.gmt),
                        'noDamageLimit': '{0}'.format(ndl),
                        'minIML': '{0}'.format(m.bounds[0]),
                        'maxIML': '{0}'.format(m.bounds[1])})

            for dsl in m.damage_state.keys():
                mean = m.damage_state[dsl][0]
                stdv = m.damage_state[dsl][1]
                par = xet.SubElement(ff, 'params', {
                            'ls': '{0}'.format(dsl),
                            'mean': '{0}'.format(mean),
                            'stddev': '{0}'.format(stdv)})

        if isinstance(m, FragilityModelDiscrete):
            ff = xet.SubElement(fm, 'fragilityFunction', {
                        'id': '{0}'.format(m.id),
                        'format': 'discrete'})

            im = xet.SubElement(ff, 'imls', {
                        'imt': '{0}'.format(m.gmt),
                        'noDamageLimit': '{0}'.format(ndl)})
            im.text = ' '.join(str(n) for n in m.gmi)

            for dsl in m.damage_state.keys():
                poe = xet.SubElement(ff, 'poes', {
                            'ls': '{0}'.format(dsl)})
                poe.text = ' '.join(str(n) for n in m.damage_state[dsl])

    indent(nrml)

    tree = xet.ElementTree(nrml)
    tree.write(xml_file, encoding='utf-8', xml_declaration=True)


def exposure_to_xml(exposure, taxonomy_tree, xml_file):
    """
    WARNING: this is specific to our case study (Friuly Region)
    and will be made more general in the future
    """
    nrml = xet.Element('nrml', {
                'xmlns': 'http://openquake.org/xmlns/nrml/0.5',
                'xmlns:gml': 'http://www.opengis.net/gml'})

    em = xet.SubElement(nrml, 'exposureModel', {
                'id': 'buildings',
                'category': 'buildings',
                'taxonomySource': 'GEM taxonomy'})

    xet.SubElement(em, 'description').text = 'Buildings'

    con = xet.SubElement(em, 'conversions')
    ctp = xet.SubElement(con, 'costTypes')

    for name in ['structural', 'nonstructural',
                 'contents', 'business_interruption']:
        xet.SubElement(ctp, 'costType', {
                    'name': name,
                    'unit': 'EUR',
                    'type': 'per_asset'})
    xet.SubElement(con, 'area', {'type': 'per_asset', 'unit': 'SQM'})
    xet.SubElement(em, 'tagNames').text = 'sezione comune'

    ass = xet.SubElement(em, 'assets')
    for li in exposure.location:
        for tax in li.taxonomy:
            tte = taxonomy_tree.get_element(tax.id)
            for bri, brw in tte.branch.items():
                nob = int(tax.number_of_buildings * brw)
                if nob > 0:
                    id = '_'.join([li.id, tax.id, bri])

                    ast = xet.SubElement(ass, 'asset', {
                            'id': id,
                            'name': li.code,
                            'area': str(li.area),
                            'number': str(nob),
                            'taxonomy': bri})
                    xet.SubElement(ast, 'location', {
                            'lon': str(li.longitude),
                            'lat': str(li.latitude)})
                    occ = xet.SubElement(ast, 'occupancies')
                    for mytime in ['day','transit','night']:
                        occ_s = xet.SubElement(occ, 'occupancy', {
                            'occupants': str(int(tax.occupants['day'])),  
                            'period': mytime})
                    cst = xet.SubElement(ast, 'costs') 
                    for mytype in ['structural', 'nonstructural',
                        'contents', 'business_interruption']:
                        cst_s = xet.SubElement(cst, 'cost', {
                            'type': mytype,
                            'value': '0'})
                    tags = xet.SubElement(ast, 'tags' , {
                            'sezione': 'all_included',  
                            'comune': li.code})

    indent(nrml)

    tree = xet.ElementTree(nrml)
    tree.write(xml_file, encoding='utf-8', xml_declaration=True)


def indent(elem, level=0):
    """
    Code from Fredrik Lundh to create pretty indented xml files
    http://effbot.org/zone/element-lib.htm#prettyprint
    """
    i = "\n" + level*"    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
