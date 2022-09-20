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
    'IRIS' : 'https://service.iris.edu',
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

FDSN_VERSION = 1

STATION_DEFAULTS = {
    "network" : "*",
    "station" : "*",
    "location" : "*",
    "channel" : "*",
    "starttime" : None,
    "endtime" : None,
    "startbefore" : None,
    "startafter" : None,
    "endbefore" : None,
    "endafter" : None,
    "sensor" : "*",
    "level" : "station",
    "includerestricted" : "true",
    "includeavailability" : "false",
    "includecomments" : "true",
    "updatedafter" : "*",
    "format" : "xml",
    "matchtimeseries" : "false",
    "nodata" : "204"
}

BOX_SEARCH_DEFAULTS = {
    "minlatitude" : -90,
    "maxlatitude" : 90,
    "minlongitude" : -180,
    "maxlongitude" : 180,
}

RAD_SEARCH_DEFAULTS = {
    "latitude" : 0,
    "longitude" : 0,
    "maxradius" : "180",
    "minradius" : 0,
}

DATASELECT_DEFAULTS = {
    "network" : "*",
    "station" : "*",
    "location" : "*",
    "channel" : "*",
    "starttime" : None,
    "endtime" : None,
    "quality" : "B",
    "nodata" : "204",
    "format" : "miniseed"
}

class FDSNClient(object):
    """
    """
    def __init__(self, data_center='ORFEUS'):
        """
        """
        self.url = _init_data_center(data_center)
        self.data = None

    def query_station(self, params=None, box_bounds=None, rad_bounds=None,
                      file_name=None, **kwargs):
        """
        """
        # Initialising parameters
        params = _params_init(params, STATION_DEFAULTS)

        # Updating parameters
        params = _params_update(params, kwargs)

        # Check for non standard values
        params = _params_check(params)

        resp = _fdsn_query(self.url, 'station', params)

        print(resp.url)
        print(resp)
        print(resp.content.decode())

    def query_data(self, params=None, file_name=None, **kwargs):
        """
        """
        if isinstance(params, (tuple, list)):
            if '.' in params[0]:
                net, sta, loc, chn = params[0].split(".")
                params = {'network' : net,
                          'station' : sta,
                          'location' : loc,
                          'channel' : chn,
                          'starttime' : params[1],
                          'endtime' : params[2]} 
            else:
                params = {'network' : params[0],
                          'station' : params[1],
                          'location' : params[2],
                          'channel' : params[3],
                          'starttime' : params[4],
                          'endtime' : params[5]}

        # Initialising parameters
        params = _params_init(params, DATASELECT_DEFAULTS)

        # Updating parameters
        params = _params_update(params, kwargs)

        # Check for non standard values
        params = _params_check(params)

        # Date conversion
        starttime = params['starttime']
        if isinstance(starttime, Date):
            params['starttime'] = starttime.get_date(dtype='s')

        endtime = params['endtime']
        if isinstance(endtime, Date):
            params['endtime'] = endtime.get_date(dtype='s')

        resp = _fdsn_query(self.url, 'dataselect', params)

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

def _init_data_center(data_center):
    """
    """
    data_center_url = None

    if data_center in DATA_CENTER_REGISTRY.keys():
        data_center_url = DATA_CENTER_REGISTRY[data_center]
    else:
        if 'http' in data_center:
            data_center_url = data_center
        else:
            raise ValueError('Not a valid data center')

    return data_center_url

def _params_init(params, defaults):
    """
    Initialising default parameters
    """
    if params:
        params = {**defaults, **params}
    else:
        params = {**defaults}

def _params_update(params, update_params):
    """
    Updating parameters
    """
    if updated_params:
        params = {**params, **updated_params}

    return params

def _params_check(params):
    """
    """
    # Checking for empty fields
    params = {k: ("*" if v=="" else v) for (k,v) in params.items()}

    # Remove None entries
    params = {k:v for (k,v) in params.items() if k is not None}

    return params

def _fdsn_query(data_center_url, interface, params):
    """
    """
    query = "/fdsnws/{0}/{1}/query".format(interface, FDSN_VERSION)
    resp = requests.get(data_center_url + query, params=params)

    return resp

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
