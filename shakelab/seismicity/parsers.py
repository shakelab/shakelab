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
Module containing earthquake catalogue parsers for different formats.
"""

import json

import shakelab.seismicity.catalogue as cat


def read(file_name, type='bin'):
    """
    """
    edb = EqDatabase()

    if type == 'bin':

        edb.load(file_name)

    elif type == 'json':

        with open(file_name, 'r') as f:
            data = json.load(f)['EqDatabase']
            edb.header = data['Header']

            for ee in data['Event']:
                event = Event(ee['Id'])

                for ms in ee['Magnitude']:
                    event.add_magnitude(ms)
                for ls in ee['Location']:
                    event.add_location(ls)

                edb.add_event(event)

    elif type == 'csv':
        pass

    elif type == 'isf':
        pass

    elif type == 'ndk':
        pass

    return edb

def write(edb, file_name, type='bin'):
    """
    """
    if not isinstance(edb, EqDatabase):
        raise ValueError('not a valid database')

    if type == 'bin':

        edb.dump(file_name)

    if type == 'csv':

        at = AsciiTable()

    elif type == 'json':

        data = edb.export_to_dict()

        with open(file_name, 'w') as f:
            json.dump(data, f, indent=True)


def read_ISF(isf_file, name=None, version=None, info=None):
    """
    Importer for ISC bulletin in ISF format
    """

    def create_event(line):
        
        id = line[6:14].strip()
        region = line[15:80]
        return cat.Event(id)

    def create_location_solution(line):

        loc = cat.LocationSolution()
        loc['Year'] = line[0:4].strip(' ')
        loc['Month'] = line[5:7].strip(' ')
        loc['Day'] = line[8:10].strip(' ')
        loc['Hour'] = line[10:13].strip(' ')
        loc['Minute'] = line[14:16].strip(' ')
        loc['Second'] = line[17:22].strip(' ')
        loc['Latitude'] = line[36:44].strip(' ')
        loc['Longitude'] = line[45:54].strip(' ')
        loc['Depth'] = line[71:76].strip(' ')
        loc['SecError'] = line[24:29].strip(' ')
        loc['DepError'] = line[78:82].strip(' ')
        loc['LocCode'] = line[118:127].strip(' ')
        loc['LocPrime'] = False
        return loc

    def create_magnitude_solution(line):

        mag = cat.MagnitudeSolution()
        mag['MagType'] = line[0:5].strip(' ')
        mag['MagSize'] = line[6:10].strip(' ')
        mag['MagError'] = line[11:14].strip(' ')
        mag['MagCode'] = line[20:29].strip(' ')
        mag['MagPrime'] = False
        return mag

    # Initialising database
    edb = cat.EqDatabase(name, version, info)

    with open(isf_file, 'r') as fid:

        # Bulletin header
        data_type = fid.readline().strip()
        if data_type != 'DATA_TYPE EVENT IMS1.0':
            raise ValueError('wring data type format')

        bulletin = fid.readline().strip()

        # Loop over events
        while True:

            line = fid.readline()

            if 'Event' in line:
                event = create_event(line)
                edb.event.append(event)

            elif 'Date' in line:
                # Loop over location solutions
                while True:
                    line2 = fid.readline()
                    if line2.strip():
                        if 'PRIME' in line2:
                            edb[-1].location[-1].prime = True
                        elif 'CENTROID' in line2:
                            pass
                        else:
                            loc = create_location_solution(line2)
                            edb[-1].location.add(loc)
                    else:
                        break

            elif 'Magnitude' in line:
                # Loop over magnitude solutions
                while True:
                    line2 = fid.readline()
                    if line2.strip():
                        mag = create_magnitude_solution(line2)
                        edb[-1].magnitude.add(mag)
                    else:
                        break

            if not line:
                break

    return edb
