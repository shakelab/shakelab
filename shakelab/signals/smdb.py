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
Coolection of methods to parse text files
"""

from shakelab.libutils.time import Date

class Itaca(object):

    def __init__(self, file=[]):
        """
        """

        # Variable initialisation
        self.head = {}
        self.data = []

        # Import ASCII file
        if file:
            self.read(file)

    def read(self, file):
        """
        """

        with open(file, 'r') as f:

            # Read header
            for i in range(0, 64):
                line = f.readline()

                key, value = line.split(':', 1)
                value = value.strip()

                cast = itaca_head_struc[key]
                if value:
                    if cast == 's':
                        self.head[key] = str(value)
                    elif cast == 'i':
                            self.head[key] = int(value)
                    elif cast == 'f':
                            self.head[key] = float(value)
                else:
                    self.head[key] = None

            # Loop over data
            for value in f:
                if value:
                    self.data.append(float(value))

    def write(self, file):
        """
        """
        pass

    def sampling_rate(self):
        """
        """
        return self.head['SAMPLING_INTERVAL_S']

    def time_date(self):
        """
        """

        date = self.head['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS']
        
        year = int(date[0:4])
        month = int(date[4:6])
        day = int(date[6:8])
        hour = int(date[9:11])
        minute = int(date[11:13])
        second = float(date[13:])

        return Date([year, month, day, hour, minute, second])


itaca_head_struc = {
    'EVENT_NAME': 's',
    'EVENT_ID': 's',
    'EVENT_DATE_YYYYMMDD': 'i',
    'EVENT_TIME_HHMMSS': 'i',
    'EVENT_LATITUDE_DEGREE': 'f',
    'EVENT_LONGITUDE_DEGREE': 'f',
    'EVENT_DEPTH_KM': 'f',
    'HYPOCENTER_REFERENCE': 's',
    'MAGNITUDE_W': 'f',
    'MAGNITUDE_W_REFERENCE': 's',
    'MAGNITUDE_L': 'f',
    'MAGNITUDE_L_REFERENCE': 's',
    'FOCAL_MECHANISM': 's',
    'NETWORK': 's',
    'STATION_CODE': 's',
    'STATION_NAME': 's',
    'STATION_LATITUDE_DEGREE': 'f',
    'STATION_LONGITUDE_DEGREE': 'f',
    'STATION_ELEVATION_M': 'f',
    'LOCATION': 's',
    'SENSOR_DEPTH_M': 'f',
    'VS30_M/S': 'f',
    'SITE_CLASSIFICATION_EC8': 's',
    'MORPHOLOGIC_CLASSIFICATION': 's',
    'EPICENTRAL_DISTANCE_KM': 'f',
    'EARTHQUAKE_BACKAZIMUTH_DEGREE': 'f',
    'DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS': 's',
    'DATE_TIME_FIRST_SAMPLE_PRECISION': 's',
    'SAMPLING_INTERVAL_S': 'f',
    'NDATA': 'i',
    'DURATION_S': 'f',
    'STREAM': 's',
    'UNITS': 's',
    'INSTRUMENT': 's',
    'INSTRUMENT_ANALOG/DIGITAL': 's',
    'INSTRUMENTAL_FREQUENCY_HZ': 's',
    'INSTRUMENTAL_DAMPING': 's',
    'FULL_SCALE_G': 's',
    'N_BIT_DIGITAL_CONVERTER': 'i',
    'PGA_CM/S^2': 'f',
    'TIME_PGA_S': 'f',
    'PGV_CM/S': 'f',
    'TIME_PGV_S': 'f',
    'PGD_CM': 'f',
    'TIME_PGD_S': 'f',
    'BASELINE_CORRECTION': 's',
    'FILTER_TYPE': 's',
    'FILTER_ORDER': 's',
    'LOW_CUT_FREQUENCY_HZ': 'f',
    'HIGH_CUT_FREQUENCY_HZ': 'f',
    'LATE/NORMAL_TRIGGERED': 's',
    'DATABASE_VERSION': 's',
    'HEADER_FORMAT': 's',
    'DATA_TYPE': 's',
    'PROCESSING': 's',
    'DATA_TIMESTAMP_YYYYMMDD_HHMMSS': 's',
    'DATA_LICENSE': 's',
    'DATA_CITATION': 's',
    'DATA_CREATOR': 's',
    'ORIGINAL_DATA_MEDIATOR_CITATION': 's',
    'ORIGINAL_DATA_MEDIATOR': 's',
    'ORIGINAL_DATA_CREATOR_CITATION': 's',
    'ORIGINAL_DATA_CREATOR': 's',
    'USER1': 's',
    'USER2': 's',
    'USER3': 's',
    'USER4': 's',
    'USER5': 's'}
