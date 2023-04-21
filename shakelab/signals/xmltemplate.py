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

from copy import deepcopy
import warnings

# Empty dictionary templates to store the minimum inventory information  
# required by the most recent FDSN StationXML schema (1.2)

t_req_float_unit = {'value':None}

t_req_float_unit_datum = {'value':None}

t_req_site = {'Name':None}

t_req_units = {'Name':None}

t_req_instrumsensitivity = {'Value':None, 
                            'Frequency':None, 
                            'InputUnits':t_req_units, 
                            'OutputUnits':t_req_units
                           }

t_req_response = {'InstrumentSensitivity':t_req_instrumsensitivity}


def initialize_metadata(level='station'):
    """
    Initialize minimal, empty dictionaries following the most recent FDSN 
    StationXML Inventory schema (1.2); for each metadata level, only elements 
    deemed as 'required' by the schema are initialized
    """
    try:
        inventory = t_req_inventory(level=level)
    except:
        warnings.warn("Warning: Metadata initialization failed",  UserWarning)
        return None
    return inventory


def t_req_inventory(level='station'):
    """Initialize empty FDSN Inventory template according to requested level"""
    t_inventory = {'schemaVersion':None,
                   'Source':None,
                   'Created':None,
                   'Network':[t_req_network(level=level)]
                  }
    return t_inventory


def t_req_network(level='station'):
    """Initialize empty network template according to requested level"""
    t_network = {'code':None}
    if level in ['station', 'channel', 'response']:
        t_network['Station'] = [t_req_station(level=level)]
    return t_network


def t_req_station(level='station'):
    """Initialize empty station template according to requested level"""
    t_station = {'code':None, 
                 'Latitude':deepcopy(t_req_float_unit_datum), 
                 'Longitude':deepcopy(t_req_float_unit_datum), 
                 'Elevation':deepcopy(t_req_float_unit), 
                 'Site':deepcopy(t_req_site)
                }
    if level in ['channel', 'response']:
        t_station['Channel'] = [t_req_channel(level=level)]
    return t_station


def t_req_channel(level='channel'):
    """Initialize empty channel template according to requested level"""
    t_channel = {'code':None, 
                 'locationCode':None, 
                 'Latitude':deepcopy(t_req_float_unit_datum), 
                 'Longitude':deepcopy(t_req_float_unit_datum), 
                 'Elevation':deepcopy(t_req_float_unit), 
                 'Depth':deepcopy(t_req_float_unit)
                }
    if level in ['response']:
        t_channel['Response'] = [deepcopy(t_req_response)]
    return t_channel
