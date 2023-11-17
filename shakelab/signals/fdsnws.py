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

from shakelab.libutils.time import Date
from shakelab.signals.base import StreamCollection
from shakelab.signals.stationxml import parse_sxml

from copy import deepcopy
import requests
import json

DATA_CENTER_REGISTRY = {
    'AUSPASS' : ' http://auspass.edu.au:8080',
    'BGR' : 'https://eida.bgr.de',
    'EMSC' : 'https://www.seismicportal.eu',
    'ETH' : 'https://eida.ethz.ch',
    'GEOFON' : 'https://geofon.gfz-potsdam.de',
    'ICGC' : 'https://ws.icgc.cat',
    'IESDMC' : 'http://batsws.earth.sinica.edu.tw',
    'INGV' : 'https://webservices.ingv.it',
    'IPGP' : 'http://ws.ipgp.fr',
    'IRISDMC' : 'https://service.iris.edu',
    'ISC' : 'http://www.isc.ac.uk',
    'KAGSR' : 'http://sdis.emsd.ru',
    'KOERI' : 'https://eida.koeri.boun.edu.tr',
    'LMU' : 'https://erde.geophysik.uni-muenchen.de',
    'NCEDC' : 'http://service.ncedc.org',
    'NIEP' : 'https://eida-sc3.infp.ro',
    'NOA' : 'http://eida.gein.noa.gr',
    'ORFEUS' : 'https://orfeus-eu.org',
    'RASPISHAKE' : 'https://fdsnws.raspberryshakedata.com',
    'RESIF' : 'https://ws.resif.fr',
    'SCEDC' : 'https://service.scedc.caltech.edu',
    'UIB-NORSAR' : 'http://eida.geo.uib.no',
    'USGS' : 'https://earthquake.usgs.gov',
    'USP' : 'http://seisrequest.iag.usp.br',
    'OGS-SCP' : 'http://158.110.30.85:8080',
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
    "level" : "response",
    "includerestricted" : "true",
    "includeavailability" : "false",
    "updateafter" : None,
    "matchtimeseries" : "false",
    "format" : "xml",
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

    def get_waveform(self, fdsn_code, starttime, endtime,
                     correct=False, file_name=None):
        """
        """
        fc = FDSNCode(fdsn_code)

        sc = self.query_data(fc.get('dict'),
                             starttime=starttime, endtime=endtime)

        if correct:
            xml = self.query_station(fc.get('dict'), level='response')
            rc = parse_sxml(xml)
            sc.deconvolve_response(rc)

        if file_name is None:
            return sc

    def query_station(self, params={}, box_bounds=None, rad_bounds=None,
                            file_name=None, **kwargs):
        """
        """
        # Initiale and update query parameters
        params = _params_update(params, STATION_DEFAULTS, **kwargs)

        # Check for non standard values
        params = _params_check(params)

        resp = _fdsn_query(self.url, 'station', params)

        if resp.content:

            if b'Error' in resp.content:
                print(resp.content.decode())

            else:
                return resp.content.decode()

        else:
            print('No station available')


    def query_data(self, params={}, file_name=None, **kwargs):
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

        # Updating parameters
        params = _params_update(params, DATASELECT_DEFAULTS, **kwargs)

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
                        sc = StreamCollection()
                        sc.read(resp.content)
                        return sc
                    else:
                        raise ValueError('Format not supported')
                else:
                    with open(file_name, 'wb') as f:
                        f.write(resp.content)

        else:
            print('No data available')

    def query_event(self):
        """
        """
        pass

    def query_info(self):
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

def _params_update(params, defaults, **kwargs):
    """
    Updating default parameters
    """
    params = {**defaults, **params}

    for key, value in kwargs.items():
        if key in defaults.keys():
            params[key] = value

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

    return {dc['name'] : dc['website'] for dc in data["datacenters"]}


class FDSNCode(object):
    """
    """
    def __init__(self, code=None):
        self.network = ''
        self.station = ''
        self.location = ''
        self.channel = ''

        if code is not None:
            self.set(code)

    def __repr__(self):
        self.get()

    def set(self, code):
        """
        """
        if isinstance(code, (list, tuple)):
            self.network = code[0]
            self.station = code[1]
            self.location = code[2]
            self.channel = code[3]

        elif isinstance(code, str):
            code = code.split('.')
            self.network = code[0]
            self.station = code[1]
            self.location = code[2]
            self.channel = code[3]

        elif isinstance(code, dict):
            self.network = code['network']
            self.station = code['station']
            self.location = code['location']
            self.channel = code['channel']

        else:
            raise TypeError('Not a supported input type')

    def get(self, dtype='str'):
        """
        """
        if dtype == 'list':
            return [self.network, self.station,
                    self.location, self.channel]

        elif dtype == 'str':
            return '{0}.{1}.{2}.{3}'.format(self.network,
                                            self.station,
                                            self.location,
                                            self.channel)

        elif dtype == 'dict':
            return {'network' : self.network,
                    'station' : self.station,
                    'location' : self.location,
                    'channel' : self.channel}

        else:
            raise TypeError('Not a supported output type')