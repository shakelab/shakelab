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

import math
import re
import warnings
from copy import deepcopy

from lxml import etree
from shakelab.libutils.time import Date as shakeDate
import shakelab.signals.xmltemplate as xmlt

# Define the recognized StationXML schemas
READABLE_VERSIONS = ("1.0", "1.1", "1.2")

msg_missingreq = "Element %s lacks required FDSN attribute or tag and thus " \
                 "it cannot be read. It will not be part of the final " \
                 "inventory object."



def _get_fdsn_info(xmldoc):
    """
    Return StationXML namespace and version string if the namespace corresponds 
    to the FDSN standard, or ``None`` otherwise.
    """
    root = xmldoc.getroot()
    try:
        match = re.match(
            r'{http://www.fdsn.org/xml/station/[0-9]+}FDSNStationXML',
            root.tag)
        root_ns = root.tag[1:].partition("}")[0]
        assert match is not None
    except Exception:
        return None
    try:
        version = root.attrib["schemaVersion"]
    except KeyError:
        return None
    return root_ns, version


def read_stationxml(path_or_file_object, level='station'):
    """
    Function reading a StationXML file. Fields required by the FDSN StationXML
    standard are checked as a condition to parse each object.

    :param [str, binary] path_or_file_object: 
        File name or file like object.
    :param [str] level: 
        Level of detail to read from file. One among 'response','channel', 
        'station' or 'network'.
    :output [dict] inv: 
        A dictionary compliant to the FDSN Inventory schema.
    """
    xmldoc = etree.parse(path_or_file_object)
    root = xmldoc.getroot()
    root_ns = _get_fdsn_info(xmldoc)[0]

    def _ns(tagname):
        return "{%s}%s" % (root_ns, tagname)

    inv = xmlt.initialize_metadata(level=level)

    # schemaVersion, Source and Created fields must be present
    inv["schemaVersion"] = root.attrib.get("schemaVersion")
    inv["Source"] = root.find(_ns("Source")).text
    inv["Created"] = shakeDate(root.find(_ns("Created")).text)

    tmp = {}
    tmp["Sender"] = _tag2obj(root, _ns("Sender"), str)
    tmp["Module"] = _tag2obj(root, _ns("Module"), str)
    tmp["ModuleURI"] = _tag2obj(root, _ns("ModuleURI"), str)
    for label in ["Sender", "Module", "ModuleURI"]:
        if (tmp[label] is not None):
            inv[label] = tmp[label]

    networks = []
    for network in root.findall(_ns("Network")):
        networks.append(_read_network(network, _ns, level))
    inv["Network"] = networks
    return inv


def _read_base_node(element, outputobj, _ns):
    """
    Reads the base node structure from element (except code) and saves it in
    outputobj dictionary.
    """
    tmp = {}
    tmp["alternateCode"] = _attr2obj(element, "alternateCode", str)
    tmp["historicalCode"] = _attr2obj(element, "historicalCode", str)
    tmp["startDate"] = _attr2obj(element, "startDate", shakeDate)
    tmp["endDate"] = _attr2obj(element, "endDate", shakeDate)
    tmp["restrictedStatus"] = _attr2obj(element, "restrictedStatus", str)
    tmp["sourceID"] = _attr2obj(element, "sourceID", str)
    tmp["Description"] = _tag2obj(element, _ns("Description"), str)
    tmp["Identifier"] = _tag2list(element, "Identifier", _ns, _read_identifier)
    tmp["Comment"] = _tag2list(element, "Comment", _ns, _read_comment)
    data_availability = element.find(_ns("DataAvailability"))
    tmp["DataAvailability"] = _read_data_availability(data_availability, _ns)
    for label in ["alternateCode", "historicalCode", "startDate", "endDate",
                  "restrictedStatus", "sourceID", "Description", "Identifier", "Comment", "DataAvailability"]:
        if tmp[label] is not None:
            outputobj[label] = tmp[label]


def _read_network(net_elem, _ns, level):
    """
    Returns either a dictionary based on the FDSN Network element schema, 
    or ``None`` if a required field is missing.
    """
    code_req = net_elem.get("code")
    if code_req is None:
        # code attribute must be present
        warnings.warn((msg_missingreq % net_elem.tag), UserWarning)
        return None
    
    net = xmlt.t_req_network(level=level)
    net["code"] = code_req

    _read_base_node(net_elem, net, _ns)
    tmp = {}
    tmp["Operator"] = _tag2list(net_elem, "Operator", _ns, _read_operator)
    tmp["TotalNumberStations"] = \
        _tag2obj(net_elem, _ns("TotalNumberStations"), int)
    tmp["SelectedNumberStations"] = \
        _tag2obj(net_elem, _ns("SelectedNumberStations"), int)
    for label in ["Operator", "TotalNumberStations", "SelectedNumberStations"]:
        if tmp[label] is not None:
            net[label] = tmp[label]    
        
    if level in ('station', 'channel', 'response'):
        stations = []
        for station in net_elem.findall(_ns("Station")):
            sta = _read_station(station, _ns, level)
            if sta is not None:
                stations.append(sta)
        # Catch if all stations are skipped
        if not stations:
            warnings.warn("Level '%s' was required but 'Station' info is "
                          "not available in parsed metadata." % level, UserWarning)
            stations = None
        net["Station"] = stations
    return net


def _read_station(sta_elem, _ns, level):
    """
    Returns either a dictionary based on the FDSN Station element schema, 
    or ``None`` if a required field is missing.
    """
    station = xmlt.t_req_station(level=level)

    code_req = sta_elem.get("code")
    latitude_req = _read_floattype(sta_elem, _ns("Latitude"), unit=True,
                                   datum=True, required=True)
    longitude_req = _read_floattype(sta_elem, _ns("Longitude"), unit=True,
                                    datum=True, required=True)
    elevation_req = _read_floattype(sta_elem, _ns("Elevation"), unit=True,
                                    required=True)
    site_req = _read_site(sta_elem.find(_ns("Site")), _ns)
    if None in [code_req, latitude_req, longitude_req, elevation_req, site_req]:
        # code, Latitude, Longitude, Elevation and Site tags must be present
        warnings.warn((msg_missingreq % sta_elem.tag), UserWarning)
        return None

    station["code"] = code_req
    _read_base_node(sta_elem, station, _ns)
    station["Latitude"] = latitude_req
    station["Longitude"] = longitude_req
    station["Elevation"] = elevation_req
    station["Site"] = site_req

    tmp = {}
    # water level only for schemas > 1.0
    tmp["WaterLevel"] = _read_floattype(sta_elem, _ns("WaterLevel"), unit=True)
    tmp["Vault"] = _tag2obj(sta_elem, _ns("Vault"), str)
    tmp["Geology"] = _tag2obj(sta_elem, _ns("Geology"), str)
    tmp["Equipment"] = _tag2list(sta_elem, "Equipment", _ns, _read_equipment)
    tmp["Operator"] = _tag2list(sta_elem, "Operator", _ns, _read_operator)
    tmp["CreationDate"] = _tag2obj(sta_elem, _ns("CreationDate"), shakeDate)
    tmp["TerminationDate"] = \
        _tag2obj(sta_elem, _ns("TerminationDate"), shakeDate)
    tmp["TotalNumberChannels"] = \
        _tag2obj(sta_elem, _ns("TotalNumberChannels"), int)
    tmp["SelectedNumberChannels"] = \
        _tag2obj(sta_elem, _ns("SelectedNumberChannels"), int)
    tmp["ExternalReference"] = \
        _tag2list(sta_elem, "ExternalReference", _ns, _read_external_reference)
    for label in ["WaterLevel", "Vault", "Geology", "Equipment", "Operator", 
                  "CreationDate", "TerminationDate", "TotalNumberChannels", 
                  "SelectedNumberChannels", "ExternalReference"]:
        if tmp[label] is not None:
            station[label] = tmp[label]

    if level in ('channel', 'response'):
        channels = []
        for channel in sta_elem.findall(_ns("Channel")):
            chan = _read_channel(channel, _ns, level)
            if chan is not None:
                channels.append(chan)
        # Catch if all channels are skipped
        if not channels:
            warnings.warn("Level '%s' was required but 'Channel' info is "
                          "not available in parsed metadata." % level, UserWarning)
            channels = None        
        station["Channel"] = channels

    return station


def _read_floattype(parent, tag, unit=False, datum=False, required=False):
    """
    Utility to parse dictionaries with float values
    """
    elem = parent.find(tag)
    if elem is None:
        return None
    convert = tagtext2obj(elem, float)
    #return None if tag text is not present, when required by FDSN XML standards
    if required:
        if convert is None:
            return None
    obj = {}
    tmp = {}
    tmp["value"] = convert
    tmp["plusError"] = _attr2obj(elem, "plusError", float)
    tmp["minusError"] = _attr2obj(elem, "minusError", float)
    tmp["measurementMethod"] = _attr2obj(elem, "measurementMethod", str)
    for label in ["value", "plusError", "minusError", "measurementMethod"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    if unit:
        tmp["unit"] = _attr2obj(elem, "unit", str)
        if tmp["unit"] is not None:
            obj["unit"] = tmp["unit"]
    if datum:
        tmp["datum"] = _attr2obj(elem, "datum", str)
        if tmp["datum"] is not None:
            obj["datum"] = tmp["datum"]
    #catch if dictionary is empty
    if not obj:
        return None
    return obj


def _read_floattype_list(parent, tag, unit=False, datum=False, number=False, 
                         required=False):
    """
    Utility to parse lists of dictionaries with float values
    """
    elems = parent.findall(tag)
    if len(elems) == 0:
        return None
    objs = []
    for elem in elems:
        convert = tagtext2obj(elem, float)
        #append None if tag text is not present, when required by FDSN XML 
        # standards
        if required:
            if convert is None:
                objs.append(None)
                continue        
        obj = {}
        tmp = {}
        tmp["value"] = convert
        tmp["plusError"] = _attr2obj(elem, "plusError", float)
        tmp["minusError"] = _attr2obj(elem, "minusError", float)
        tmp["measurementMethod"] = _attr2obj(elem, "measurementMethod", str)
        for label in ["value", "plusError", "minusError", "measurementMethod"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        if unit:
            tmp["unit"] = _attr2obj(elem, "unit", str)
            if tmp["unit"] is not None:
                obj["unit"] = tmp["unit"]            
        if datum:
            tmp["datum"] = _attr2obj(elem, "datum", str)
            if tmp["datum"] is not None:
                obj["datum"] = tmp["datum"]
        if number:
            tmp["number"] = _attr2obj(elem, "number", int)
            if tmp["number"] is not None:
                obj["number"] = tmp["number"]
        #catch if dictionary is empty
        if not obj:
            continue
        objs.append(obj)

    #catch if list is empty or malformed
    if objs.count(None) == len(objs):
        return None
    return objs

#capire come implementare selezione su locationcode
def _read_channel(cha_elem, _ns, level):
    """
    Returns either a dictionary based on the FDSN Channel element schema, 
    or ``None`` if a required field is missing.
    """
    channel = xmlt.t_req_channel(level=level)

    code_req = cha_elem.get("code")
    location_code_req = cha_elem.get("locationCode")
    latitude_req = _read_floattype(cha_elem, _ns("Latitude"),
                                   unit=True, datum=True, required=True)
    longitude_req = _read_floattype(cha_elem, _ns("Longitude"),
                                    unit=True, datum=True, required=True)
    elevation_req = _read_floattype(cha_elem, _ns("Elevation"),
                                    unit=True, required=True)
    depth_req = _read_floattype(cha_elem, _ns("Depth"), unit=True,
                                required=True)
    if None in [code_req, location_code_req, latitude_req, longitude_req, 
                elevation_req, depth_req]:
        # code, locationCode, Latitude, Longitude, Elevation and Depth 
        # tags must be present
        warnings.warn((msg_missingreq % cha_elem.tag), UserWarning)
        return None

    channel["code"] = code_req
    channel["locationCode"] = location_code_req
    _read_base_node(cha_elem, channel, _ns)
    tmp = {}
    tmp["ExternalReference"] = \
        _tag2list(cha_elem, "ExternalReference", _ns, _read_external_reference)
    if tmp["ExternalReference"] is not None:
        channel["ExternalReference"] = tmp["ExternalReference"]
    channel["Longitude"] = longitude_req
    channel["Latitude"] = latitude_req
    channel["Elevation"] = elevation_req
    channel["Depth"] = depth_req
    tmp["Azimuth"] = _read_floattype(cha_elem, _ns("Azimuth"), unit=True)
    tmp["Dip"] = _read_floattype(cha_elem, _ns("Dip"), unit=True)
    # water level only for schemas > 1.0
    tmp["WaterLevel"] = _read_floattype(cha_elem, _ns("WaterLevel"), unit=True)
    tmp["Type"] = _tags2obj(cha_elem, _ns("Type"), str)
    tmp["SampleRate"] = _read_floattype(cha_elem, _ns("SampleRate"), unit=True)
    tmp["SampleRateRatio"] = _read_sample_rate_ratio(cha_elem, _ns)
    tmp["ClockDrift"] = _read_floattype(cha_elem, _ns("ClockDrift"), unit=True)
    calibration_units = cha_elem.find(_ns("CalibrationUnits"))
    tmp["CalibrationUnits"] = _read_units(calibration_units, _ns)
    sensor = cha_elem.find(_ns("Sensor"))
    tmp["Sensor"] = _read_equipment(sensor, _ns)
    pre_amplifier = cha_elem.find(_ns("PreAmplifier"))
    tmp["PreAmplifier"] = _read_equipment(pre_amplifier, _ns)
    data_logger = cha_elem.find(_ns("DataLogger"))
    tmp["DataLogger"] = _read_equipment(data_logger, _ns)
    tmp["Equipment"] = _tag2list(cha_elem, "Equipment", _ns, _read_equipment)

    for label in ["Azimuth", "Dip", "WaterLevel", "Type", "SampleRate", 
                  "SampleRateRatio", "ClockDrift", "CalibrationUnits", 
                  "Sensor", "PreAmplifier", "DataLogger", "Equipment"]:
        if tmp[label] is not None:
            channel[label] = tmp[label]

    if level == 'response':
        resp = _read_response(cha_elem, _ns)
        fdsnhint = False
        if resp is not None:
            try:
                hint1 = ("Polynomial" in resp["Stage"])
                hint2 = ("InstrumentSensitivity" in resp)
                fdsnhint = (hint1 or hint2)
            except:
                pass
        if not fdsnhint:
            # FDSN schema strongly suggests that either InstrumentSensitivity  
            # or InstrumentPolynomial be present
            warnings.warn("Warning: neither InstrumentSensitivity nor "
                          "InstrumentPolynomial are available - this is "
                          "discouraged in FDSN schema " %  UserWarning)     
        channel["Response"] = resp
    return channel


def _read_response(cha_elem, _ns):
    """
    Utility to parse FDSN Response dictionaries
    """
    obj = None
    response = cha_elem.find(_ns("Response"))
    if response is not None:
        obj = {}
        tmp = {}
        tmp["resourceId"] = _attr2obj(response, "resourceId", str) 
        tmp["InstrumentSensitivity"] = \
            _read_instrument_sensitivity(response, _ns)

        instrument_polynomial = response.find(_ns("InstrumentPolynomial"))
        tmp["InstrumentPolynomial"] = \
            _read_polynomial(instrument_polynomial, _ns)
        tmp["Stage"] = _tag2list(response, "Stage", _ns, _read_stage)
        for label in ["resourceId", "InstrumentSensitivity", 
                      "InstrumentPolynomial", "Stage"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        # catch malformed response objects
        if not obj:
            obj = None
    return obj


def _read_stage(stage_elem, _ns):
    """
    Utility to parse all FDSN Stage types; returns ``None`` if a required 
    field is missing.
    """
    number_req = _attr2obj(stage_elem, "number", int)

    if number_req is None:
        # number tag must be present
        warnings.warn((msg_missingreq % stage_elem.tag), UserWarning)
        return None

    obj = {}
    obj["number"] = number_req
    resourceid = _attr2obj(stage_elem, "resourceId", str)
    if resourceid is not None:
        obj["resourceId"] = resourceid

    stage_polynomial = stage_elem.find(_ns("Polynomial"))
    polynomial = _read_polynomial(stage_polynomial, _ns)

    if polynomial is not None:
        obj["Polynomial"] = polynomial

    else:
        gain_req = _read_stage_gain(stage_elem, _ns)
        if gain_req is None:
            # stagegain attribute must be present
            warnings.warn((msg_missingreq % stage_elem.tag), UserWarning)
            return None
        
        tmp = {}
        tmp["PolesZeros"] = _read_poleszeros(stage_elem, _ns)
        tmp["Coefficients"] = _read_coefficients(stage_elem, _ns)
        tmp["ResponseList"] = _read_resplist(stage_elem, _ns)
        tmp["FIR"] = _read_fir(stage_elem, _ns)

        label = []
        if tmp["PolesZeros"] is not None:
            label.append("PolesZeros")
        if tmp["Coefficients"] is not None:
            label.append("Coefficients")
        if tmp["ResponseList"] is not None:
            label.append("ResponseList")
        if tmp["FIR"] is not None:
            label.append("FIR")
        if len(label) != 1:
            # exactly one response attribute must be present
            warnings.warn("Element '%s' has either none or more than one "
                "response type." % stage_elem.tag, UserWarning)
            return None      
              
        obj[label[0]] = tmp[label[0]]
        decimation = _read_decimation(stage_elem, _ns)
        if decimation is not None:
            obj["Decimation"] = decimation         
        obj["StageGain"] = gain_req

    return obj


def _read_stage_gain(stage_elem, _ns):
    """
    Utility to parse FDSN Gain dictionaries; returns ``None`` if 
    a required field is missing.
    """    
    obj = None
    stage_gain = stage_elem.find(_ns("StageGain"))
    if stage_gain is not None:
        value_req = _tag2obj(stage_gain, _ns("Value"), float)
        frequency_req = _tag2obj(stage_gain, _ns("Frequency"), float)
        if None in [value_req, frequency_req]:
            # Value and Frequency tags must be present
            warnings.warn((msg_missingreq % stage_elem.tag), UserWarning)
            return None
        obj = {}
        obj["Value"] = value_req
        obj["Frequency"] = frequency_req
    return obj

def _read_poleszeros(stage_elem, _ns):
    """
    Utility to parse FDSN PolesZeros dictionaries; returns ``None`` if 
    a required field is missing.
    """   
    obj = None
    pz = stage_elem.find(_ns("PolesZeros"))
    if pz is not None:
        input_units = pz.find(_ns("InputUnits"))
        input_units_req = _read_units(input_units, _ns)
        output_units = pz.find(_ns("OutputUnits"))
        output_units_req = _read_units(output_units, _ns)
        pz_transfer_function_type_req = \
            _tag2obj(pz, _ns("PzTransferFunctionType"), str)
        normalization_factor_req = \
            _tag2obj(pz, _ns("NormalizationFactor"), float)
        normalization_frequency_req = \
            _read_floattype(pz, _ns("NormalizationFrequency"), unit=True,
                            required=True)

        if None in [input_units_req, output_units_req,  
                    pz_transfer_function_type_req, normalization_factor_req, 
                    normalization_frequency_req]:
            # InputUnits, OutputUnits, PzTransferFunctionType, 
            # NormalizationFactor and NormalizationFrequency tags 
            # must be present
            warnings.warn((msg_missingreq % pz.tag), UserWarning)
            return None

        obj = {}
        tmp = {}
        tmp["name"] = _attr2obj(pz, "name", str)
        tmp["resourceId"] = _attr2obj(pz, "resourceId", str)
        tmp["Description"] = _tag2obj(pz, _ns("Description"), str)
        for label in ["name", "resourceId", "Description"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        obj["InputUnits"] = input_units_req
        obj["OutputUnits"] = output_units_req
        obj["PzTransferFunctionType"] = pz_transfer_function_type_req
        obj["NormalizationFactor"] = normalization_factor_req
        obj["NormalizationFrequency"] = normalization_frequency_req
        tmp["Zero"] = _tag2list(pz, "Zero", _ns, _tag2pole_or_zero)
        tmp["Pole"] = _tag2list(pz, "Pole", _ns, _tag2pole_or_zero)
        for label in ["Zero", "Pole"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
    return obj

def _tag2pole_or_zero(element, _ns):
    """
    Utility to parse Pole and Zero dictionaries; returns ``None`` if 
    a required field is missing.
    """
    real_req = _read_floattype(element, _ns("Real"), required=True)
    imaginary_req = _read_floattype(element, _ns("Imaginary"), required=True)
    if None in [real_req, imaginary_req]:
        # Real and Imaginary tags must be present
        warnings.warn((msg_missingreq % element.tag), UserWarning)
        return None    
    obj = {}
    tmp = {}
    tmp["number"] = _attr2obj(element, "number", int)
    if tmp["number"] is not None:
        obj["number"] = tmp["number"]
    obj["Real"] = real_req
    obj["Imaginary"] = imaginary_req
    return obj

def _read_coefficients(stage_elem, _ns):
    """
    Utility to parse FDSN Coefficients dictionaries; returns ``None`` if 
    a required field is missing.
    """
    obj = None
    coeffs = stage_elem.find(_ns("Coefficients"))
    if coeffs is not None:
        input_units = coeffs.find(_ns("InputUnits"))
        input_units_req = _read_units(input_units, _ns)
        output_units = coeffs.find(_ns("OutputUnits"))
        output_units_req = _read_units(output_units, _ns)
        cf_transfer_function_type_req = \
            _tag2obj(coeffs, _ns("CfTransferFunctionType"), str)        
        if None in [input_units_req, output_units_req, 
                    cf_transfer_function_type_req]:
            # InputUnits, OutputUnits,CfTransferFunctionType tags must be 
            # present
            warnings.warn((msg_missingreq % coeffs.tag), UserWarning)
            return None
        
        obj = {}
        tmp = {}
        tmp["name"] = _attr2obj(coeffs, "name", str)
        tmp["resourceId"] = _attr2obj(coeffs, "resourceId", str)
        tmp["Description"] = _tag2obj(coeffs, _ns("Description"), str)
        for label in ["name", "resourceId", "Description"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        obj["InputUnits"] = input_units_req
        obj["OutputUnits"] = output_units_req
        obj["CfTransferFunctionType"] = cf_transfer_function_type_req
        tmp["Numerator"] = \
            _read_floattype_list(coeffs, _ns("Numerator"), number=True)
        tmp["Denominator"] = \
            _read_floattype_list(coeffs, _ns("Denominator"), number=True)
        for label in ["Numerator", "Denominator"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
    return obj


def _read_resplist(stage_elem, _ns):
    """
    Utility to parse FDSN ResponseList dictionaries; returns ``None`` if 
    a required field is missing.
    """
    obj = None
    resp = stage_elem.find(_ns("ResponseList"))
    if resp is not None:
        input_units = resp.find(_ns("InputUnits"))
        input_units_req = _read_units(input_units, _ns)
        output_units = resp.find(_ns("OutputUnits"))
        output_units_req = _read_units(output_units, _ns)
        if None in [input_units_req, output_units_req]:
            # InputUnits and OutputUnits tags must be present
            warnings.warn((msg_missingreq % resp.tag), UserWarning)
            return None
        
        obj = {}
        tmp = {}
        tmp["name"] = _attr2obj(resp, "name", str)
        tmp["resourceId"] = _attr2obj(resp, "resourceId", str)
        tmp["Description"] = _tag2obj(resp, _ns("Description"), str)
        for label in ["name", "resourceId", "Description"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        obj["InputUnits"] = input_units_req
        obj["OutputUnits"] = output_units_req
        tmp["ResponseListElement"] = \
            _tag2list(resp, "ResponseListElement", _ns, _read_resplist_elem)
        if tmp["ResponseListElement"] is not None:
            obj["ResponseListElement"] = tmp["ResponseListElement"]
    return obj


def _read_resplist_elem(element, _ns):
    """
    Utility to parse FDSN ResponseList elements; returns ``None`` if 
    a required field is missing.
    """
    frequency_req = _read_floattype(element, _ns("Frequency"), unit=True, 
                                    required=True)
    amplitude_req = _read_floattype(element, _ns("Amplitude"), unit=True, 
                                    required=True)
    phase_req = _read_floattype(element, _ns("Phase"), unit=True, required=True)

    if None in [frequency_req, amplitude_req, phase_req]:
        # Frequency, Amplitude and Phase tags must be present
        warnings.warn((msg_missingreq % element.tag), UserWarning)
        return None    
    obj = {}
    obj["Frequency"] = frequency_req
    obj["Amplitude"] = amplitude_req
    obj["Phase"] = phase_req
    return obj


def _read_fir(stage_elem, _ns):
    """
    Utility to parse FDSN FIR dictionaries; returns ``None`` if a required field
    is missing.
    """
    obj = None
    fir = stage_elem.find(_ns("FIR"))
    if fir is not None:
        input_units = fir.find(_ns("InputUnits"))
        input_units_req = _read_units(input_units, _ns)
        output_units = fir.find(_ns("OutputUnits"))
        output_units_req = _read_units(output_units, _ns)
        symmetry_req = _tag2obj(fir, _ns("Symmetry"), str)
        if None in [input_units_req, output_units_req, symmetry_req]:
            # InputUnits, OutputUnits and Symmetry tags must be present
            warnings.warn((msg_missingreq % fir.tag), UserWarning)
            return None

        obj = {}
        tmp = {}
        tmp["name"] = _attr2obj(fir, "name", str)
        tmp["resourceId"] = _attr2obj(fir, "resourceId", str)
        tmp["Description"] = _tag2obj(fir, _ns("Description"), str)
        for label in ["name", "resourceId", "Description"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
        obj["InputUnits"] = input_units_req
        obj["OutputUnits"] = output_units_req
        obj["Symmetry"] = symmetry_req
        tmp["NumeratorCoefficient"] = \
            _tag2list(fir, "NumeratorCoefficient", _ns, _read_numerator_coeff)
        if tmp["NumeratorCoefficient"] is not None:
            obj["NumeratorCoefficient"] = tmp["NumeratorCoefficient"]
    return obj


def _read_numerator_coeff(element, _ns):
    """
    Utility to parse FDSN NumeratorCoefficient dictionaries; returns 
    ``None`` if a required field is missing.
    """
    obj = {}
    tmp = {}
    tmp["value"] = tagtext2obj(element, float)
    tmp["i"] = _attr2obj(element, "i", int)
    for label in ["value", "i"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    #catch if numerator coeffcient is empty
    if not obj:
        return None
    return obj


def _read_decimation(stage_elem, _ns):
    """
    Utility to parse FDSN Decimation dictionaries; returns ``None`` if a 
    required field is missing.
    """
    obj = None
    decim = stage_elem.find(_ns("Decimation"))
    if decim is not None:
        sample_rate_req = _read_floattype(decim, _ns("InputSampleRate"),
                                          unit=True, required=True)
        factor_req = _tag2obj(decim, _ns("Factor"), int)
        offset_req = _tag2obj(decim, _ns("Offset"), int)
        delay_req = _read_floattype(decim, _ns("Delay"), unit=True,
                                    required=True)
        correction_req = _read_floattype(decim, _ns("Correction"), unit=True,
                                         required=True)
        if None in [sample_rate_req, factor_req, offset_req, delay_req, 
                    correction_req]:
            # InputSampleRate, Factor, Offset, Delay and Correction tags 
            # must be present
            warnings.warn((msg_missingreq % stage_elem.tag), UserWarning)
            return None
        obj = {}
        obj["InputSampleRate"] = sample_rate_req
        obj["Factor"] = factor_req
        obj["Offset"] = offset_req
        obj["Delay"] = delay_req
        obj["Correction"] = correction_req
    return obj


def _read_instrument_sensitivity(response_elem, _ns):
    """
    Utility to parse FDSN InstrumentSensitivity dictionaries; returns ``None`` 
    if a required field is missing.
    """
    obj = None
    sensitivity = response_elem.find(_ns("InstrumentSensitivity"))
    if sensitivity is not None:
        value_req = _tag2obj(sensitivity, _ns("Value"), float)
        frequency_req = _tag2obj(sensitivity, _ns("Frequency"), float)
        input_units = sensitivity.find(_ns("InputUnits"))
        input_units_req = _read_units(input_units, _ns)
        output_units = sensitivity.find(_ns("OutputUnits"))
        output_units_req = _read_units(output_units, _ns)
        if None in [value_req, frequency_req, input_units_req, 
                    output_units_req]:
            # Value, Frequency, InputUnits and OutputUnits tags must be present
            warnings.warn((msg_missingreq % response_elem.tag), UserWarning)
            return None
        tmp = {}
        obj = {}
        obj["Value"] = value_req
        obj["Frequency"] = frequency_req
        obj["InputUnits"] = input_units_req
        obj["OutputUnits"] = output_units_req
        tmp["FrequencyStart"] = \
            _tag2obj(sensitivity, _ns("FrequencyStart"), float)
        tmp["FrequencyEnd"] = \
            _tag2obj(sensitivity, _ns("FrequencyEnd"), float)
        tmp["FrequencyDBVariation"] = \
            _tag2obj(sensitivity, _ns("FrequencyDBVariation"), float)
        for label in ["FrequencyStart", "FrequencyEnd", "FrequencyDBVariation"]:
            if tmp[label] is not None:
                obj[label] = tmp[label]
    return obj


def _read_polynomial(poly_elem, _ns):
    """
    Utility to parse FDSN Stage Polynomial and InstrumentPolynomial
    dictionaries; returns ``None`` if a required field is missing.
    """
    if poly_elem is None:
        return None

    input_units = poly_elem.find(_ns("InputUnits"))
    input_units_req = _read_units(input_units, _ns)
    output_units = poly_elem.find(_ns("OutputUnits"))
    output_units_req = _read_units(output_units, _ns)
    approx_type_req = _tag2obj(poly_elem, _ns("ApproximationType"), str)
    f_lower_req = _read_floattype(poly_elem, _ns("FrequencyLowerBound"), 
                                  unit=True, required=True)
    f_upper_req = _read_floattype(poly_elem, _ns("FrequencyUpperBound"), 
                                  unit=True, required=True)
    approx_lower_req = \
        _tag2obj(poly_elem, _ns("ApproximationLowerBound"), float)
    approx_upper_req = \
        _tag2obj(poly_elem, _ns("ApproximationUpperBound"), float)
    max_err_req = _tag2obj(poly_elem, _ns("MaximumError"), float)
    coeffs_req = _read_floattype_list(poly_elem, _ns("Coefficient"),
                                      number=True, required=True)
    if None in [input_units_req, output_units_req, approx_type_req,
                f_lower_req, f_upper_req, approx_lower_req, approx_upper_req,
                max_err_req, coeffs_req]:
        # InputUnits, OutputUnits, ApproximationType, FrequencyLowerBound,
        # FrequencyUpperBound, ApproximationLowerBound, ApproximationUpperBound,
        # MaximumError and Coefficient tags must be present
        warnings.warn((msg_missingreq % poly_elem.tag), UserWarning)
        return None

    obj = {}
    tmp = {}
    tmp["name"] = _attr2obj(poly_elem, "name", str)
    tmp["resourceId"] = _attr2obj(poly_elem, "resourceId", str)
    tmp["Description"] = _tag2obj(poly_elem, _ns("Description"), str)
    for label in ["name", "resourceId", "Description"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    obj["InputUnits"] = input_units_req
    obj["OutputUnits"] = output_units_req
    obj["ApproximationType"] = approx_type_req
    obj["FrequencyLowerBound"] = f_lower_req
    obj["FrequencyUpperBound"] = f_upper_req
    obj["ApproximationLowerBound"] = approx_lower_req
    obj["ApproximationUpperBound"] = approx_upper_req
    obj["MaximumError"] = max_err_req
    obj["Coefficient"] = coeffs_req
    return obj


def _read_external_reference(ref_elem, _ns):
    """
    Utility to parse FDSN ExternalReference dictionaries; returns 
    ``None`` if a required field is missing.
    """
    uri_req = _tag2obj(ref_elem, _ns("URI"), str)
    description_req = _tag2obj(ref_elem, _ns("Description"), str)
    if None in [uri_req, description_req]:
        # URI and Description tags must be present
        warnings.warn((msg_missingreq % ref_elem.tag), UserWarning)
        return None
    obj = {}
    obj["URI"] = uri_req
    obj["Description"] = description_req
    return obj


def _read_operator(operator_elem, _ns):
    """
    Utility to parse FDSN Operator dictionaries; returns 
    ``None`` if a required field is missing.
    """
    agency_req = operator_elem.find(_ns("Agency")).text
    if agency_req is None:
        # Agency tag must be present
        warnings.warn((msg_missingreq % operator_elem.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    obj["Agency"] = agency_req
    tmp["Contact"] = _tag2list(operator_elem, "Contact", _ns, _read_person)
    tmp["Website"] = _tag2obj(operator_elem, _ns("WebSite"), str)
    for label in ["Contact", "Website"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    return obj


def _read_data_availability(avail_elem, _ns):
    """
    Utility to parse FDSN DataAvailability dictionaries; returns ``None`` if a
    required field is missing.
    """
    if avail_elem is None:
        return None
    obj = {}
    tmp = {}
    extent = avail_elem.find(_ns("Extent"))
    spans = avail_elem.findall(_ns("Span"))
    # Recovery from empty Extent tag + no spans.
    if extent is None and not spans:
        return None
    tmp["Extent"] = _read_data_availability_extent(extent, _ns) 
    tmp["Span"] = \
        _tag2list(avail_elem, "Span", _ns, _read_data_availability_span)
    for label in ["Extent", "Span"]:
        if tmp[label] is not None:
            obj[label] = tmp[label] 
    if not obj:
        return None  
    return obj


def _read_data_availability_extent(element, _ns):
    if element is None:
        return None
    end_req = element.attrib['end']
    start_req = element.attrib['start']
    if None in [end_req, start_req]:
        # end and start attributes must be present
        warnings.warn((msg_missingreq % element.tag), UserWarning)
        return None
    obj = {}
    obj["end"] = shakeDate(end_req)
    obj["start"] = shakeDate(start_req)
    return obj


def _read_data_availability_span(element, _ns):
    end_req = element.attrib['end']
    number_of_segments_req = element.attrib['numberSegments']
    start_req = element.attrib['start']
    if None in [end_req, number_of_segments_req, start_req]:
        # end, numberSegments and start attributes must be present
        warnings.warn((msg_missingreq % element.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    obj["end"] = shakeDate(end_req)
    tmp["maximumTimeTear"] = _attr2obj(element, "maximumTimeTear", float)
    if tmp["maximumTimeTear"] is not None:
        obj["maximumTimeTear"] = tmp["maximumTimeTear"]
    obj["numberSegments"] = int(number_of_segments_req)
    obj["start"] = shakeDate(start_req)
    return obj


def _read_sample_rate_ratio(cha_elem, _ns):
    """
    Utility to parse FDSN SampleRateRatio dictionaries; returns ``None`` if a
    required field is missing.
    """
    # Parse the optional sample rate ratio
    obj = None
    sample_rate_ratio = cha_elem.find(_ns("SampleRateRatio"))
    if sample_rate_ratio is not None:
        number_samples_req = \
            _tag2obj(sample_rate_ratio, _ns("NumberSamples"), int)
        number_seconds_req = \
            _tag2obj(sample_rate_ratio, _ns("NumberSeconds"), int)
        if None in [number_samples_req, number_seconds_req]:
            # NumberSamples and NumberSeconds tags must be present
            warnings.warn((msg_missingreq % sample_rate_ratio.tag), UserWarning)
            return None
        obj = {}
        obj["NumberSamples"] = number_samples_req
        obj["NumberSeconds"] = number_seconds_req
    return obj


def _read_equipment(equip_elem, _ns):
    """
    Utility to parse FDSN Sensor, PreAmplifier, DataLogger and Equipment
    dictionaries.
    """
    if equip_elem is None:
        return None
    obj = {}
    tmp = {}
    tmp["resourceId"] = _attr2obj(equip_elem, "resourceId", str)
    tmp["Type"] = _tag2obj(equip_elem, _ns("Type"), str)
    tmp["Description"] = _tag2obj(equip_elem, _ns("Description"), str)
    tmp["Manufacturer"] = _tag2obj(equip_elem, _ns("Manufacturer"), str)
    tmp["Vendor"] = _tag2obj(equip_elem, _ns("Vendor"), str)
    tmp["Model"] = _tag2obj(equip_elem, _ns("Model"), str)
    tmp["SerialNumber"] = _tag2obj(equip_elem, _ns("SerialNumber"), str)
    tmp["InstallationDate"] = \
        _tag2obj(equip_elem, _ns("InstallationDate"), shakeDate)
    tmp["RemovalDate"] = \
        _tag2obj(equip_elem, _ns("RemovalDate"), shakeDate)
    tmp["CalibrationDate"] = \
        _tags2obj(equip_elem, _ns("CalibrationDate"), shakeDate)
    for label in ["resourceId", "Type", "Description", "Manufacturer", 
                  "Vendor", "Model", "SerialNumber", "InstallationDate", 
                  "RemovalDate", "CalibrationDate"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    #catch if equipment is empty
    if not obj:
        return None
    return obj

def _read_units(unit_elem, _ns):
    """
    Utility to parse FDSN InputUnits, OutputUnits and CalibrationUnits
    dictionaries; returns ``None`` if a required field is missing.
    """
    if unit_elem is None:
        return None
    name_req = unit_elem.find(_ns("Name")).text
    if name_req is None:
        # name tag must be present
        warnings.warn((msg_missingreq % unit_elem.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    obj["Name"] = name_req
    tmp["Description"] = _tag2obj(unit_elem, _ns("Description"), str)
    if tmp["Description"] is not None:
        obj["Description"] = tmp["Description"]
    return obj


def _read_site(site_elem, _ns):
    """
    Utility to parse FDSN Site dictionaries; returns ``None`` if a required
    field is missing.
    """
    if site_elem is None:
        return None
    name_req = site_elem.find(_ns("Name")).text
    if name_req is None:
        # name tag must be present
        warnings.warn((msg_missingreq % site_elem.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    obj["Name"] = name_req
    tmp["Description"] = _tag2obj(site_elem, _ns("Description"), str)
    tmp["Town"] = _tag2obj(site_elem, _ns("Town"), str)
    tmp["County"] = _tag2obj(site_elem, _ns("County"), str)
    tmp["Region"] = _tag2obj(site_elem, _ns("Region"), str)
    tmp["Country"] = _tag2obj(site_elem, _ns("Country"), str)
    for label in ["Description", "Town", "County", "Region", "Country"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    return obj

def _read_identifier(element, _ns):
    """
    Utility to parse FDSN Identifier dictionaries.
    """
    obj = {}
    tmp = {}
    tmp["value"] = element.text
    tmp["type"] = element.get('type')
    for label in ["value", "type"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    if not obj:
        return None
    return obj


def _read_comment(comment_elem, _ns):
    """
    Utility to parse FDSN Comment dictionaries; returns ``None`` if a 
    required field is missing.
    """
    value_req = _tag2obj(comment_elem, _ns("Value"), str)
    if value_req is None:
        # Value tag must be present
        warnings.warn((msg_missingreq % comment_elem.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    obj["Value"] = value_req
    tmp["id"] = _attr2obj(comment_elem, "id", int)
    tmp["subject"] = _attr2obj(comment_elem, "subject", str)
    tmp["BeginEffectiveTime"] = \
        _tag2obj(comment_elem, _ns("BeginEffectiveTime"), shakeDate)
    tmp["EndEffectiveTime"] = \
        _tag2obj(comment_elem, _ns("EndEffectiveTime"), shakeDate)
    tmp["Author"] = _tag2list(comment_elem, "Author", _ns, _read_person)
    for label in ["id", "tysubjectpe", "BeginEffectiveTime", 
                  "EndEffectiveTime", "Author"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    return obj


def _read_person(person_elem, _ns):
    """
    Utility to parse FDSN Author and Contact dictionaries.
    """
    person = {}
    tmp = {}
    tmp["Name"] = _tags2obj(person_elem, _ns("Name"), str)
    tmp["Agency"] = _tags2obj(person_elem, _ns("Agency"), str)
    tmp["Email"] = _tags2obj(person_elem, _ns("Email"), str)
    tmp["Phone"] = _tag2list(person_elem, "Phone", _ns, _read_phone)
    for label in ["Name", "Agency", "Email", "Phone"]:
        if tmp[label] is not None:
            person[label] = tmp[label]
    #catch if contact is empty
    if not person:
        return None
    return person


def _read_phone(phone_elem, _ns):
    """
    Utility to parse FDSN Phone dictionaries; returns ``None`` if a 
    required field is missing.
    """
    area_code_req = _tag2obj(phone_elem, _ns("AreaCode"), int)
    phone_number_req = _tag2obj(phone_elem, _ns("PhoneNumber"), int)
    if None in [area_code_req, phone_number_req]:
        # AreaCode and PhoneNumber tags must be present
        warnings.warn((msg_missingreq % phone_elem.tag), UserWarning)
        return None
    obj = {}
    tmp = {}
    tmp["description"] = _attr2obj(phone_elem, "description", str)
    tmp["CountryCode"] = _tag2obj(phone_elem, _ns("CountryCode"), int)
    for label in ["description", "CountryCode"]:
        if tmp[label] is not None:
            obj[label] = tmp[label]
    obj["AreaCode"] = area_code_req
    obj["PhoneNumber"] = phone_number_req
    return obj

####DA QUI finire le spiegazioni

def tagtext2obj(element, convert):
    """ Utility to get text value from "element" tag object"""
    try:
        obj = convert(element.text)
    except Exception:
        warnings.warn(
            "'%s' could not be converted to a float. Will be skipped. Please "
            "contact to report this issue." % etree.tostring(element),
            UserWarning)
        return None
    # Catch NaNs.
    if math.isnan(obj):
        warnings.warn("Tag '%s' has a value of NaN. It will be skipped." %
                      element.getparent().tag, UserWarning)
        return None
    return obj

def _convert(text, convert):
    if convert is str:
        return _convert_str(text)
    try:
        return convert(text)
    except Exception:
        return None


def _convert_str(text):
    """ Utility to catch empty trings"""
    if not text:
        return ''
    return text

def _tag2list(element, tag, _ns, read_function):
    """ Utility to apply read_function parser to each instance in "element" 
    tag objects"""
    elements = element.findall(_ns(tag))
    if not len(elements):
        return None
    obj = [read_function(elem, _ns) for elem in elements]
    if obj.count(None) == len(obj):
        return None
    return obj


def _tag2obj(element, tag, convert):
    try:
        text = element.find(tag).text
    except Exception:
        return None
    return _convert(text, convert)


def _tags2obj(element, tag, convert):
    elements = element.findall(tag)
    if not len(elements):
        return None
    values = []
    for elem in elements:
        values.append(_convert(elem.text, convert))
        if not values:
            return None
    return values


def _attr2obj(element, attr, convert):
    attribute = element.get(attr)
    if attribute is None:
        return None
    try:
        return convert(attribute)
    except Exception:
        return None
