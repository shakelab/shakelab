# ****************************************************************************
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
# ****************************************************************************

import io
import json
import logging
from collections import namedtuple
from copy import deepcopy

import requests
from requests.exceptions import RequestException, Timeout

import shakelab.signals.base as base
from shakelab.libutils.time import Date
from shakelab.signals.io import writer
from shakelab.signals.stationxml import parse_sxml
from shakelab.seismicity.quakeml import parse_quakeml

USE_LIBMSEED = True
FDSN_VERSION = 1
FDSNResponse = namedtuple('FDSNResponse', ['data', 'content_type', 'error'])

# Set up logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

# Registry of known FDSN data centers
DATA_CENTER_REGISTRY = {
    'AUSPASS': 'http://auspass.edu.au:8080',
    'BGR': 'https://eida.bgr.de',
    'EMSC': 'https://www.seismicportal.eu',
    'ETH': 'https://eida.ethz.ch',
    'GEOFON': 'https://geofon.gfz-potsdam.de',
    'GEONET': 'https://www.geonet.org.nz',
    'ICGC': 'https://ws.icgc.cat',
    'IESDMC': 'http://batsws.earth.sinica.edu.tw',
    'IGN': 'http://www.ign.es',
    'INGV': 'https://webservices.ingv.it',
    'IPGP': 'http://ws.ipgp.fr',
    'IRISDMC': 'https://service.iris.edu',
    'ISC': 'http://www.isc.ac.uk',
    'KAGSR': 'http://sdis.emsd.ru',
    'KOERI': 'https://eida.koeri.boun.edu.tr',
    'LMU': 'https://erde.geophysik.uni-muenchen.de',
    'NCEDC': 'http://service.ncedc.org',
    'NIEP': 'https://eida-sc3.infp.ro',
    'NOA': 'http://eida.gein.noa.gr',
    'ORFEUS': 'https://orfeus-eu.org',
    'RASPISHAKE': 'https://fdsnws.raspberryshakedata.com',
    'RESIF': 'https://ws.resif.fr',
    'SCEDC': 'https://service.scedc.caltech.edu',
    'UIB-NORSAR': 'http://eida.geo.uib.no',
    'USGS': 'https://earthquake.usgs.gov',
    'USP': 'http://seisrequest.iag.usp.br',
    'OGS-SCP': 'http://158.110.30.85:8080',
    'OGS-ANT': 'http://158.110.30.202:5600'
}

# Default parameters for station queries
STATION_DEFAULTS = {
    "network": "*",
    "station": "*",
    "location": "*",
    "channel": "*",
    "starttime": None,
    "endtime": None,
    "startbefore": None,
    "startafter": None,
    "endbefore": None,
    "endafter": None,
    "level": "response",
    "includerestricted": "true",
    "includeavailability": "false",
    "updateafter": None,
    "matchtimeseries": "false",
    "format": "xml",
    "nodata": "204"
}

# Default parameters for waveform data queries
DATASELECT_DEFAULTS = {
    "network": "*",
    "station": "*",
    "location": "*",
    "channel": "*",
    "starttime": None,
    "endtime": None,
    "quality": "B",
    "minimumlength": None,
    "longestonly": None,
    "validate": None,
    "nodata": "204",
    "format": "miniseed"
}

# Default parameters for event data queries
EVENT_DEFAULTS = {
    "starttime": None,
    "endtime": None,
    "minlatitude": None,
    "maxlatitude": None,
    "minlongitude": None,
    "maxlongitude": None,
    "latitude": None,
    "longitude": None,
    "minradius": None,
    "maxradius": None,
    "mindepth": None,
    "maxdepth": None,
    "minmagnitude": None,
    "maxmagnitude": None,
    "magnitudetype": None,
    "eventtype": None,
    "catalog": None,
    "contributor": None,
    "eventid": None,
    "limit": None,
    "offset": None,
    "format": "xml",
    "nodata": "204"
}
 
def get_mseed_module(use_libmseed=USE_LIBMSEED):
    """
    Dynamically import the MiniSEED module depending on configuration.

    Parameters
    ----------
    use_libmseed : bool
        If True, use cymseed (C-based), else pure Python mseed reader.

    Returns
    -------
    module
        The imported MiniSEED I/O module.
    """
    if use_libmseed:
        from shakelab.signals.libio import cymseed as mseed
    else:
        from shakelab.signals.libio import mseed
    return mseed


class FDSNClient:
    """
    Client to interact with FDSN web services.
    """

    def __init__(self, data_center='ORFEUS'):
        """
        Initialize FDSN client with given data center.
        """
        self.url = self._init_data_center(data_center)

    def get_waveform(self, fdsn_code, starttime, endtime,
                     correct=False, file_name=None, format='mseed'):
        """
        Retrieve waveform data and optionally save to file.
    
        Parameters
        ----------
        fdsn_code : str
            FDSN code (network.station.location.channel).
        starttime : str
            Start time of the waveform request (ISO8601 format).
        endtime : str
            End time of the waveform request (ISO8601 format).
        correct : bool, optional
            Apply instrument response correction (default is False).
        file_name : str, optional
            If given, save the waveform to this file.
        format : str, optional
            Output file format (default is 'mseed').
    
        Returns
        -------
        StreamCollection or bool
            Waveform data as StreamCollection if file_name is None,
            otherwise True if saved to file.
        """
        fc = FDSNCode(fdsn_code)
        sc = self.query_data(fc.get('dict'),
                             starttime=starttime,
                             endtime=endtime)

        # Cut data precisely to defined time window
        for stream in sc:
            for record in stream:
                record.cut(starttime, endtime)

        if sc and correct:
            xml = self.query_station(fc.get('dict'),
                                     starttime=starttime,
                                     endtime=endtime,
                                     level='response')
            rc = parse_sxml(xml)
            sc.deconvolve_response(rc)

        if file_name:
            writer(sc, file_name, format)
            return True
        else:
            return sc

    def get_response(self, fdsn_code, starttime, endtime, file_name=None):
        """
        Retrieve station metadata (including response) and
        optionally save to file.
    
        Parameters
        ----------
        fdsn_code : str
            FDSN code (network.station.location.channel).
        starttime : str
            Start time of the metadata request (ISO8601 format).
        endtime : str
            End time of the metadata request (ISO8601 format).
        file_name : str, optional
            If given, save the station metadata to this file.
    
        Returns
        -------
        dict or bool
            Parsed station metadata as dictionary if file_name is None,
            otherwise True if saved to file.
        """
        fc = FDSNCode(fdsn_code)
        xml = self.query_station(fc.get('dict'),
                                 starttime=starttime,
                                 endtime=endtime,
                                 level='response')

        if file_name:
            logger.info("Writing response to file.")
            with open(file_name, 'w') as f:
                f.write(xml)
            return True
        else:
            return parse_sxml(xml)

    def get_catalogue(self, starttime, endtime,
                      magnitude=None, latitude=None, longitude=None,
                      depth=None, file_name=None):
        """
        Retrieve earthquake events catalogue and optionally save to file.
    
        Parameters
        ----------
        starttime : str
            Start time of the catalogue request (ISO8601 format).
        endtime : str
            End time of the catalogue request (ISO8601 format).
        magnitude : list or tuple, optional
            (min_magnitude, max_magnitude) bounds.
        latitude : list or tuple, optional
            (min_latitude, max_latitude) bounds.
        longitude : list or tuple, optional
            (min_longitude, max_longitude) bounds.
        depth : list or tuple, optional
            (min_depth, max_depth) bounds (in kilometers).
        file_name : str, optional
            If given, save the catalogue to this file.
    
        Returns
        -------
        str or bool
            XML data as string, or True if saved to file.
        """
        params = {
            'starttime': starttime,
            'endtime': endtime,
        }
    
        if magnitude is not None:
            params['minmagnitude'] = magnitude[0]
            params['maxmagnitude'] = magnitude[1]
    
        if latitude is not None:
            params['minlatitude'] = latitude[0]
            params['maxlatitude'] = latitude[1]
    
        if longitude is not None:
            params['minlongitude'] = longitude[0]
            params['maxlongitude'] = longitude[1]
    
        if depth is not None:
            if depth[0] is not None:
                params['mindepth'] = depth[0]
            if depth[1] is not None:
                params['maxdepth'] = depth[1]
    
        data = self.query_event(params)
    
        if file_name:
            logger.info(f"Writing catalogue to {file_name}")
            with open(file_name, 'w') as f:
                f.write(data)
            return True

        if data:
            return parse_quakeml(data)
        else:
            logger.warning("No catalogue data returned.")
            return None

    def query_station(self, params=None, **kwargs):
        """
        Query station metadata from the FDSN service.
        """
        params = _params_update(params or {}, STATION_DEFAULTS, **kwargs)
        params = _params_check(params)

        response = _fdsn_query(self.url, 'station', params)

        if response.error:
            logger.error(f"Error fetching station metadata: {response.error}")
            return None

        return (response.data.decode()
                if isinstance(response.data, bytes) else response.data)

    def query_data(self, params=None, file_name=None, **kwargs):
        """
        Query waveform data from the FDSN service.
        """
        params = _params_update(params or {}, DATASELECT_DEFAULTS, **kwargs)
        params = _params_check(params)

        response = _fdsn_query(self.url, 'dataselect', params)

        if response.error:
            logger.error(f"Error fetching waveform data: {response.error}")
            return None

        if file_name:
            with open(file_name, 'wb') as f:
                f.write(response.data)
            return True

        if params['format'] == 'miniseed':
            sc = base.StreamCollection()
            mseed = get_mseed_module(USE_LIBMSEED)
            mseed.msread(response.data, stream_collection=sc)
            return sc

        raise ValueError(f"Unsupported format: {params['format']}")

    def query_event(self, params=None, **kwargs):
        """
        Query earthquake event data from the FDSN service.
        """
        params = _params_update(params or {}, EVENT_DEFAULTS, **kwargs)
        params = _params_check(params)

        response = _fdsn_query(self.url, 'event', params)

        if response.error:
            logger.error(f"Error fetching event data: {response.error}")
            return None

        return (response.data.decode()
                if isinstance(response.data, bytes) else response.data)

    def query_info(self, service='station'):
        """
        Query WADL service information for the specified service.
        """
        full_url = (f"{self.url}/fdsnws/{service}/{FDSN_VERSION}/"
                    "application.wadl")

        try:
            response = requests.get(full_url, timeout=5)
            response.raise_for_status()
            return response.text

        except RequestException as e:
            logger.error(f"Error fetching service info: {e}")
            return None

    def _init_data_center(self, data_center):
        """
        Resolve and return the URL for the specified data center.
        """
        if data_center in DATA_CENTER_REGISTRY:
            return DATA_CENTER_REGISTRY[data_center]

        if 'http' in data_center:
            return data_center

        raise ValueError(f"Invalid data center: {data_center}")


class FDSNCode:
    """
    Helper class to parse and format FDSN codes.
    """
    _KEYS = ('network', 'station', 'location', 'channel')

    def __init__(self, code=None, **kwargs):
        """
        Initialize FDSNCode instance.
        """
        for key in self._KEYS:
            setattr(self, key, '')

        if code:
            self.set(code)

        for key, value in kwargs.items():
            if key in self._KEYS:
                setattr(self, key, value)

    def __repr__(self):
        return self.get('str')

    def __eq__(self, other):
        """
        Compare two FDSNCode objects.
        """
        if isinstance(other, str):
            return self.get('str') == other

        if isinstance(other, dict):
            return self.get('dict') == other

        if isinstance(other, (list, tuple)):
            return self.get('list') == other

        return False

    def set(self, code):
        """
        Set code values from a string, dict or list.
        """
        if isinstance(code, str):
            parts = code.split('.')
            for key, value in zip(self._KEYS, parts):
                setattr(self, key, value)

        elif isinstance(code, dict):
            for key, value in code.items():
                setattr(self, key, value)

        elif isinstance(code, (list, tuple)):
            for key, value in zip(self._KEYS, code):
                setattr(self, key, value)

        else:
            raise TypeError("Unsupported input type")

    def get(self, format='str'):
        """
        Return code in specified format.
        """
        if format == 'str':
            return (f"{self.network}.{self.station}."
                    f"{self.location}.{self.channel}")

        if format == 'dict':
            return {key: getattr(self, key) for key in self._KEYS}

        if format == 'list':
            return [getattr(self, key) for key in self._KEYS]

        raise TypeError("Unsupported output format")


def _params_update(params, defaults, **kwargs):
    """
    Update parameters based on defaults and additional keyword arguments.
    """
    updated = {**defaults, **params, **kwargs}
    return updated


def _params_check(params):
    """
    Clean up parameters, replacing empty string fields with wildcard '*'.
    """
    if isinstance(params.get('starttime'), Date):
        params['starttime'] = params['starttime'].iso8601

    if isinstance(params.get('endtime'), Date):
        params['endtime'] = params['endtime'].iso8601

    param_clean = {}
    for key, value in params.items():
        if value is None:
            continue
        if isinstance(value, str) and value == '':
            param_clean[key] = '*'
        else:
            param_clean[key] = value
    return param_clean


def _fdsn_query(url, service, params, retries=3, timeout=10):
    """
    Perform a robust HTTP query to the specified FDSN service.
    """
    full_url = f"{url}/fdsnws/{service}/{FDSN_VERSION}/query"

    logger.info(f"Querying {service} at {full_url} with params: {params}")

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(full_url, params=params, timeout=timeout)
            response.raise_for_status()

            content_type = response.headers.get('Content-Type', '').lower()

            if 'application/json' in content_type:
                return FDSNResponse(json.loads(response.content),
                                    'json', None)

            if 'application/xml' in content_type:
                return FDSNResponse(response.text, 'text', None)

            if 'text/' in content_type:
                return FDSNResponse(response.text, 'text', None)

            return FDSNResponse(response.content, 'bytes', None)

        except Timeout:
            logger.warning(f"Timeout attempt {attempt}/{retries} ",
                           f"on {full_url}")

        except RequestException as e:
            logger.error(f"Request error: {e} on attempt ",
                         f"{attempt}/{retries}")

    return FDSNResponse(None, None, "Failed after retries")


def get_fdsn_data_center_registry():
    """
    Fetch the current FDSN data center registry from official website.
    """
    url = "https://www.fdsn.org/ws/datacenters/1/query"
    requests.packages.urllib3.disable_warnings()

    try:
        resp = requests.get(url, verify=False)
        data = resp.json()
        return {dc['name']: dc['website'] for dc in data['datacenters']}

    except Exception as e:
        logger.error(f"Failed to fetch FDSN datacenter registry: {e}")
        return {}
