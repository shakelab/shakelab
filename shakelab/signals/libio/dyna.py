# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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
'''
Functionalities to handle the DYNA (1.2) file format
(from ITACA - Italian Accelerometric Archive)
'''
import os
from shakelab.libutils.time import Date
from shakelab.signals import base

DYNA_VERSION=1.2

def dynaread(input_file, units='A'):
    '''
    '''
    dy = Dyna()
    dy.read(input_file)

    record = base.Record()
    record.head.delta = dy.header['SAMPLING_INTERVAL_S']
    record.head.time = Date()
    record.head.time.strparse(
        dy.header['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS'],
        '%Y%m%d_%H%M%S'
        )
    record.head.sid = '{0}.{1}.{2}.{3}'.format(
        dy.header['NETWORK'],
        dy.header['STATION_CODE'],
        dy.header['LOCATION'],
        dy.header['STREAM']
        )
    record.data = dy.data
    return record

def dynawrite(record, output_file):
    '''
    '''
    dy = Dyna(
        sid=record.head.sid,
        data=record.data,
        starttime=record.time,
        dt=record.delta,
        units='A')
    dy.write(output_file)

# Template header dictionary with default values
HEADER_DEFAULTS = {
    'EVENT_NAME': '',
    'EVENT_ID': '',
    'EVENT_DATE_YYYYMMDD': '',
    'EVENT_TIME_HHMMSS': '',
    'EVENT_LATITUDE_DEGREE': None,
    'EVENT_LONGITUDE_DEGREE': None,
    'EVENT_DEPTH_KM': None,
    'HYPOCENTER_REFERENCE': '',
    'MAGNITUDE_W': None,
    'MAGNITUDE_W_REFERENCE': '',
    'MAGNITUDE_L': None,
    'MAGNITUDE_L_REFERENCE': '',
    'FOCAL_MECHANISM': 'U',
    'NETWORK': '',
    'STATION_CODE': '',
    'STATION_NAME': '',
    'STATION_LATITUDE_DEGREE': None,
    'STATION_LONGITUDE_DEGREE': None,
    'STATION_ELEVATION_M': None,
    'LOCATION': '',
    'SENSOR_DEPTH_M': None,
    'VS30_M/S': None,
    'SITE_CLASSIFICATION_EC8': '',
    'MORPHOLOGIC_CLASSIFICATION': '',
    'EPICENTRAL_DISTANCE_KM': None,
    'EARTHQUAKE_BACKAZIMUTH_DEGREE': None,
    'DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS': '',
    'DATE_TIME_FIRST_SAMPLE_PRECISION': '',
    'SAMPLING_INTERVAL_S': None,
    'NDATA': None,
    'DURATION_S': None,
    'STREAM': '',
    'UNITS': 'cm/s2',
    'INSTRUMENT': '',
    'INSTRUMENT_ANALOG/DIGITAL': 'D',
    'INSTRUMENTAL_FREQUENCY_HZ': None,
    'INSTRUMENTAL_DAMPING': None,
    'FULL_SCALE_G': None,
    'N_BIT_DIGITAL_CONVERTER': None,
    'PGA_CM/S^2': None,
    'PGV_CM/S': None,
    'PGD_CM': None,
    'TIME_PGA_S': None,
    'TIME_PGV_S': None,
    'TIME_PGD_S': None,
    'BASELINE_CORRECTION': 'BASELINE NOT REMOVED',
    'FILTER_TYPE': '',
    'FILTER_ORDER': None,
    'LOW_CUT_FREQUENCY_HZ': None,
    'HIGH_CUT_FREQUENCY_HZ': None,
    'LATE/NORMAL_TRIGGERED': 'NT',
    'DATABASE_VERSION': '',
    'HEADER_FORMAT': 'DYNA 1.2',
    'DATA_TYPE': 'PGA',
    'PROCESSING': '',
    'DATA_TIMESTAMP_YYYYMMDD_HHMMSS': '',
    'DATA_LICENSE': '',
    'DATA_CITATION': '',
    'DATA_CREATOR': '',
    'ORIGINAL_DATA_MEDIATOR_CITATION': '',
    'ORIGINAL_DATA_MEDIATOR': '',
    'ORIGINAL_DATA_CREATOR_CITATION': '',
    'ORIGINAL_DATA_CREATOR': '',
    'USER1': '',
    'USER2': '',
    'USER3': '',
    'USER4': '',
    'USER5': ''
}

# Type specification for each header field
HEADER_TYPES = {
    'EVENT_LATITUDE_DEGREE': float,
    'EVENT_LONGITUDE_DEGREE': float,
    'EVENT_DEPTH_KM': float,
    'MAGNITUDE_W': float,
    'MAGNITUDE_L': float,
    'STATION_LATITUDE_DEGREE': float,
    'STATION_LONGITUDE_DEGREE': float,
    'STATION_ELEVATION_M': float,
    'SENSOR_DEPTH_M': float,
    'VS30_M/S': float,
    'EPICENTRAL_DISTANCE_KM': float,
    'EARTHQUAKE_BACKAZIMUTH_DEGREE': float,
    'SAMPLING_INTERVAL_S': float,
    'NDATA': int,
    'DURATION_S': float,
    'INSTRUMENTAL_FREQUENCY_HZ': float,
    'INSTRUMENTAL_DAMPING': float,
    'FULL_SCALE_G': float,
    'N_BIT_DIGITAL_CONVERTER': int,
    'PGX_UNITS': float,
    'TIME_PGX_S': float,
    'FILTER_ORDER': int,
    'LOW_CUT_FREQUENCY_HZ': float,
    'HIGH_CUT_FREQUENCY_HZ': float,
    'PGA_CM/S^2': float,
    'PGV_CM/S': float,
    'PGD_CM': float,
    'TIME_PGA_S': float,
    'TIME_PGV_S': float,
    'TIME_PGD_S': float
}

class Dyna:
    '''
    A class to handle ITACA ASCII (DYNA 1.2) files.

    Attributes:
    - header (dict): Contains the header information.
    - data (list): Contains the waveform data.
    - starttime (datetime): Time of the first sample.
    - dt (float): Sampling interval.
    '''
    
    def __init__(self, sid=None, starttime=None, dt=None,
                       data=None, units='A'):
        '''
        Initializes the Dyna object with optional data, starttime, and dt.
        
        Parameters:
        - data (list): List of waveform data (default: None).
        - starttime (datetime): Time of the first sample (default: None).
        - dt (float): Sampling interval in seconds (default: None).
        '''
        self.units = units
        self._header_init()
        self.data = data if data else []

        if sid:
            code = sid.split('.')
            self.header['NETWORK'] = code[0]
            self.header['STATION_CODE'] = code[1]
            self.header['LOCATION'] = code[2]
            self.header['STREAM'] = code[3]
        if starttime:
            if not isinstance(starttime, Date):
                starttime = Date(starttime)
            dstr = starttime.strformat('%Y%m%d_%H%M%S')
            self.header['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS'] = dstr
        if dt:
            self.header['SAMPLING_INTERVAL_S'] = dt
        self.header['NDATA'] = len(self.data) if data else 0

    def _header_init(self):
        '''
        '''
        self.header = HEADER_DEFAULTS.copy()

        if self.units == 'A':
            self.header['UNITS'] = 'cm/s^2'
            for key in ['PGV_CM/S', 'PGD_CM', 'TIME_PGV_S', 'TIME_PGD_S']:
                self.header.pop(key)
        elif self.units == 'V':
            self.header['UNITS'] = 'cm/s'
            for key in ['PGA_CM/S2', 'PGD_CM', 'TIME_PGA_S', 'TIME_PGD_S']:
                self.header.pop(key)
        elif self.units == 'D':
            self.header['UNITS'] = 'cm'
            for key in ['PGA_CM/S2', 'PGV_CM/S', 'TIME_PGA_S', 'TIME_PGV_S']:
                self.header.pop(key)
        else:
            raise ValueError(f'Units {units} not recognized')

    def _convert_header_value(self, key, value):
        '''
        Converts a header value to the appropriate data type based on the key.
        
        Parameters:
        - key (str): The header field key.
        - value (str): The value to be converted.
        
        Returns:
        - The value converted to the appropriate type (int, float, or str).
        '''
        if key in HEADER_DEFAULTS:
            if key in HEADER_TYPES:
                target_type = HEADER_TYPES[key]
                try:
                    return target_type(value)
                except ValueError:
                    return None  # If conversion fails, return None
            else:
                return value
        else:
            raise ValueError(f'Header field {key} not recognized')

    def read(self, input_file):
        '''
        Reads the ITACA ASCII (DYNA 1.2) file into the Dyna object.
        
        Parameters:
        - input_file (str): Path to the ITACA ASCII file.
        '''
        if not os.path.exists(input_file):
            raise FileNotFoundError(f'File not found: {input_file}')

        with open(input_file, 'r') as f:
            self.header = {}
            # Read the header: the first 64 lines
            for _ in range(64):
                line = f.readline().strip()
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    self.header[key] = self._convert_header_value(key, value)
            
            # Read the waveform data
            self.data = []
            for line in f:
                try:
                    self.data.append(float(line.strip()))
                except ValueError:
                    continue

    def write(self, output_file):
        '''
        Writes the Dyna object's header and data to an ITACA ASCII
        (DYNA 1.2) file.
        
        Parameters:
        - output_file (str): Path where the file will be written.
        '''
        with open(output_file, 'w') as f:
            # Write the header
            for key, value in self.header.items():
                if value is None or value == '':
                    f.write(f'{key}: \n')
                elif isinstance(value, float):
                    f.write(f'{key}: {value:.6f}\n')
                elif isinstance(value, int):
                    f.write(f'{key}: {value}\n')
                else:
                    f.write(f'{key}: {value}\n')
            
            # Ensure the header has 64 lines
            num_header_lines = len(self.header)
            if num_header_lines < 64:
                for _ in range(64 - num_header_lines):
                    f.write('\n')
            
            # Write the waveform data
            for datum in self.data:
                f.write(f'{datum:.6e}\n')


