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
"""
from copy import deepcopy

import logging
import json
import io

import requests
from requests.exceptions import RequestException, Timeout

import shakelab.signals.base as base
from shakelab.libutils.time import Date
from shakelab.signals.stationxml import parse_sxml
from shakelab.signals.io import writer

USE_LIBMSEED = True

def get_mseed_module(use_libmseed=USE_LIBMSEED):
    """
    Load the appropriate module conditionally.
    """
    if use_libmseed:
        from shakelab.signals.libio import cymseed as mseed
    else:
        from shakelab.signals.libio import mseed
    return mseed

FDSN_VERSION = 1

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

# Set up logging configuration
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class FDSNClient(object):
    """
    """
    def __init__(self, data_center='ORFEUS'):
        """
        """
        self.url = _init_data_center(data_center)

    def get_waveform(self, fdsn_code, starttime, endtime,
                     correct=False, output=None, format='mseed'):
        """
        """
        fc = FDSNCode(fdsn_code)

        sc = self.query_data(fc.get('dict'),
                             starttime=starttime, endtime=endtime)

        if sc is not None:
            if correct:
                xml = self.query_station(fc.get('dict'), level='response')
                rc = parse_sxml(xml)
                sc.deconvolve_response(rc)

            if output is None:
                return sc
            else:
                writer(output, sc, format)
                return True

    def get_response(self, fdsn_code, starttime, endtime, file_name=None):
        """
        """
        fc = FDSNCode(fdsn_code)
        xml = self.query_station(fc.get('dict'), level='response')

        if file_name is None:
            return parse_sxml(xml)
        else:
            logger.info("Writing response to file.")
            with open(file_name, 'w') as f:
                f.write(xml)

    def query_station(self, params={}, box_bounds=None, rad_bounds=None,
                            file_name=None, **kwargs):
        """
        """
        # Initiale and update query parameters
        params = _params_update(params, STATION_DEFAULTS, **kwargs)

        # Check for non standard values
        params = _params_check(params)

        # Fetch data using the _fdsn_query function
        resp_data, content_type, error = _fdsn_query(
            self.url, 'station', params
            )
        
        # Check for errors
        if error:
            print(f"Error occurred: {error}")
        else:
            if resp_data:  # Check if there is any data
                if isinstance(resp_data, bytes) and b'Error' in resp_data:
                    # Decode bytes and print the error message
                    print(resp_data.decode()) 
                else:
                    # Return the decoded content if it's not an error
                    if isinstance(resp_data, bytes):
                        return resp_data.decode()
                    else:
                        return resp_data
            else:
                print('No station available')
                return None


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


        # Convert time objects to iso8601 strings
        starttime = params['starttime']
        if isinstance(starttime, Date):
            params['starttime'] = starttime.iso8601
    
        endtime = params['endtime']
        if isinstance(endtime, Date):
            params['endtime'] = endtime.iso8601

        # Check for non standard values
        params = _params_check(params)

        resp_data, content_type, error = _fdsn_query(
            self.url, 'dataselect', params
            )
        
        if error:
            print(f"Error occurred: {error}")
            return None
        else:
            if resp_data:  # Check if there is data
                if content_type == 'bytes' and b'Error' in resp_data:
                    print(resp_data.decode())
                else:
                    if file_name is None:
                        if params['format'] == 'miniseed':
                            sc = base.StreamCollection()
        
                            mseed = get_mseed_module(USE_LIBMSEED)
        
                            # Import data from miniseed binary buffer
                            mseed.msread(resp_data, stream_collection=sc)
        
                            # Cut waveform to proper time window (TO CHECK)
                            for stream in sc:
                                for record in stream:
                                    record.cut(starttime, endtime, True)
                            return sc
                        else:
                            raise ValueError('Format not supported')
                    else:
                        with open(file_name, 'wb') as f:
                            f.write(resp_data)
        
            else:
                print('No data available')
                return None

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

def _fdsn_query(
        url: str, service: str, params: dict,
        retries: int = 3, timeout: int = 10
        ) -> tuple:
    """
    Queries the FDSN service with robust error handling and retry mechanism.
    
    Args:
        url (str): Base URL of the FDSN service.
        service (str): The specific service being queried ('station', 'event', etc.).
        params (dict): Query parameters.
        retries (int): Number of retry attempts in case of failure (default is 3).
        timeout (int): Timeout duration in seconds for each request.
    
    Returns:
        tuple: (response_data, content_type, error_message)
               - response_data: The content of the response (JSON, text, or bytes).
               - content_type: Type of the response ('json', 'text', 'bytes').
               - error_message: Error message if applicable, otherwise None.
    """
    full_url = f"{url}/fdsnws/{service}/{FDSN_VERSION}/query"
    logger.info(
        f"Querying {service} service at {full_url} with params: {params}"
        )

    for attempt in range(retries):
        try:
            response = requests.get(full_url, params=params, timeout=timeout)

            # Raises HTTPError for bad HTTP codes
            response.raise_for_status()
            
            content_type = response.headers.get('Content-Type', '').lower()
            logger.info(f"Received response with content type: {content_type}")
            
            if 'application/json' in content_type:
                return response.json(), 'json', None
            elif 'text/' in content_type:
                return response.text, 'text', None
            return response.content, 'bytes', None

        except Timeout:
            logger.warning(
                f"Timeout on attempt {attempt}/{retries} for URL: {full_url}"
                )
            if attempt == retries:
                return None, None, "Timeout occurred after multiple attempts"

        except requests.HTTPError as e:
            logger.error(
                f"HTTP error: {e} on attempt {attempt}/{retries}"
                )
            if attempt == retries:
                return None, None, f"HTTP error: {e}"

        except requests.RequestException as e:
            logger.error(
                f"Request exception: {e} on attempt {attempt}/{retries}"
                )
            if attempt == retries:
                return None, None, f"Request exception: {e}"

    logger.error(
        f"Failed to fetch data from {full_url} after {retries} attempts."
        )
    return None, None, f"Failed after {retries} attempts"


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
    _KEYMAP = {
        'network' : '',
        'station' : '',
        'location' : '',
        'channel' : ''}

    def __init__(self, code=None, **kwargs):

        # Initialise attributes to default value
        for key in self._KEYMAP:
            self._update_attribute(key, self._KEYMAP[key])

        # Update attributes as a unique block argument
        if code is not None:
            self.set(code)

        # Update attributes as individual arguments
        for key in kwargs:
            self._update_attribute(key, kwargs[key])

    def __repr__(self):
        return self.get('str')

    def __eq__(self, code):

        if isinstance(code, str):
            return self.get('str') == code

        if isinstance(code, dict):
            return self.get('dict') == code

        if isinstance(code, (list, tuple)):
            return self.get('list') == code

        else:
            return False

    def _update_attribute(self, key, value):
        """
        """
        if key in self._KEYMAP:
            exec('self.{0}={1}[0]'.format(key, [value]))

    def _get_attribute(self, key):
        """
        """
        return eval('self.{0}'.format(key))

    def set(self, code=None):
        """
        """
        if code is not None:
            if isinstance(code, str):
                self.set(dict(zip(self._KEYMAP.keys(), code.split('.'))))

            elif isinstance(code, dict):
                for key in code:
                    self._update_attribute(key, code[key])

            elif isinstance(code, (list, tuple)):
                self.set(dict(zip(self._KEYMAP.keys(), code)))

            else:
                raise TypeError('Not a supported input type')

    def get(self, dtype='str'):
        """
        """
        if dtype == 'str':
            return '{0}.{1}.{2}.{3}'.format(self.network,
                                            self.station,
                                            self.location,
                                            self.channel)

        elif dtype == 'dict':
            return {'network' : self.network,
                    'station' : self.station,
                    'location' : self.location,
                    'channel' : self.channel}

        elif dtype == 'list':
            return [self.network,
                    self.station,
                    self.location,
                    self.channel]

        else:
            raise TypeError('Not a supported output type')

