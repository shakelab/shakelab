# ****************************************************************************
#
# Copyright (C) 2019-2022, ShakeLab Developers.
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

from shakelab.libutils.time import Date
from shakelab.signals.libio.mseed import ByteStream, MiniSeed

import requests
import json

DATA_CENTER_REGISTRY = {
    'AUSPASS' : 'http://auspass.edu.au',
    'BGR' : 'https://www.bgr.bund.de',
    'EMSC' : 'https://www.emsc-csem.org',
    'ETH' : 'https://www.ethz.ch',
    'GEOFON' : 'https://geofon.gfz-potsdam.de',
    'ICGC' : 'https://www.icgc.cat/en/terratremols',
    'IESDMC' : 'http://batsws.earth.sinica.edu.tw/fdsnws',
    'INGV' : 'http://www.ingv.it',
    'IPGP' : 'http://datacenter.ipgp.fr',
    'IRISDMC' : 'https://ds.iris.edu',
    'ISC' : 'http://www.isc.ac.uk',
    'KAGSR' : 'http://sdis.emsd.ru',
    'KOERI' : 'https://www.koeri.boun.edu.tr',
    'LMU' : 'https://www.uni-muenchen.de',
    'NCEDC' : 'https://ncedc.org',
    'NIEP' : 'https://www.infp.ro',
    'NOA' : 'http://eida.gein.noa.gr',
    'ORFEUS' : 'https://orfeus-eu.org',
    'RASPISHAKE' : 'http://raspberryshake.org',
    'RESIF' : 'https://seismology.resif.fr',
    'SCEDC' : 'https://scedc.caltech.edu',
    'UIB-NORSAR' : 'http://eida.geo.uib.no',
    'USGS' : 'https://earthquake.usgs.gov',
    'USP' : 'http://www.moho.iag.usp.br',
    'OGS-SCP3' : 'http://158.110.30.85:8080',
    'OGS-ANT' : 'http://158.110.30.202:5600'
}

defaults = {
    "starttime" : None,
    "endtime" : None,
    "network" : "*",
    "station" : "*",
    "location" : "*",
    "channel" : "*",
    "quality" : "B",
    "nodata" : "404",
    "format" : "miniseed"
}

"""
fc = FDSNClient()
fc.get_data(params={"network":"NL",
                    "station": "HGN",
                    "starttime": "2017-01-01T00:00:00",
                    "endtime": "2017-01-01T00:01:00"})
"""

class FDSNClient(object):
    """
    """
    def __init__(self, data_center='ORFEUS'):
        """
        """
        self.url = None

        if data_center in DATA_CENTER_REGISTRY.keys():
            self.url = DATA_CENTER_REGISTRY[data_center]
        else:
            if 'http' in data_center:
                self.url = data_center
            else:
                raise ValueError('Not a valid data center')

        self.data = None

    def get_station(self):
        """
        """
        pass

    def get_data(self, params=None, file_name=None, **kwargs):
        """
        """
        # Initialising / adding default parameters
        if params is None:
            params = {**defaults}
        else:
            if isinstance(params, dict):
                params = {**defaults, **params}
            elif isinstance(params, (tuple, list)):
                if '.' in params[0]:
                    net, sta, loc, chn = params[0].split(".")
                    params = {**defaults,
                              'network' : net,
                              'station' : sta,
                              'location' : loc,
                              'channel' : chn,
                              'starttime' : params[1],
                              'endtime' : params[2]} 
                else:
                    params = {**defaults,
                              'network' : params[0],
                              'station' : params[1],
                              'location' : params[2],
                              'channel' : params[3],
                              'starttime' : params[4],
                              'endtime' : params[5]}

        # Updating default parameters with keyword arguments
        for key, value in kwargs.items():
            if key in defaults.keys():
                params[key] = value

        # Check for empty fields
        for key in params:
            if params[key] == '':
                params[key] = '*'

        # Date conversion
        starttime = params['starttime']
        if isinstance(starttime, Date):
            params['starttime'] = starttime.get_date(dtype='s')

        endtime = params['endtime']
        if isinstance(endtime, Date):
            params['endtime'] = endtime.get_date(dtype='s')

        version = 1
        query = "/fdsnws/dataselect/{0}/query".format(version)

        resp = requests.get(self.url + query, params=params)

        if resp.content:

            if b'Error' in resp.content:
                print(resp.content.decode())

            else:
                if file_name is None:
                    if params['format'] == 'miniseed':
                        byte_stream = ByteStream(resp.content)
                        self.data = MiniSeed(byte_stream)
                    else:
                        raise ValueError('Output file name must be specified')
                else:
                    with open(file_name, 'wb') as f:
                        f.write(resp.content)

        else:
            print('No data available')

    def get_event(self):
        """
        """
        pass

    def get_info(self):
        """
        """
        pass

def get_fdsn_data_center_registry():
    """
    Data centers from the FDSN registry
    """
    url = "https://www.fdsn.org/ws/datacenters/1/query"

    # Temporary patch of ssl problem, must be resolved
    requests.packages.urllib3.disable_warnings()

    resp = requests.get(url, verify=False)
    data = json.loads(resp.content.decode())

    for dc in data["datacenters"]:
        print('{0:15} {1}'.format(dc['name'], dc['website']))
