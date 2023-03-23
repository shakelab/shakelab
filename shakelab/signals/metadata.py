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
"""

from copy import deepcopy
import warnings

from shakelab.signals import base
from shakelab.signals.xmltemplate import initialize_metadata
from shakelab.signals.xmlparser import read_stationxml
from shakelab.signals.polezero import SensorResponse, paz_map
from shakelab.libutils.time import Date


class Metadata():
    """
    A class used to handle waveform metadata; all the information is currently 
    stored in the 'data' attribute. 
    TODO: implement a proper set of classes to handle different metadata parts, 
    and the instrumental response in particular.

    ...

    Attributes
    ----------
    data : dict
        General container for metadata info, may be related to more than one 
        station or channel. 
        (default: ``shakelab.signals.metadatatemplate.initialize_metadata``, 
        empty FDSN-compliant dictionary)

    Methods
    -------
    read(mdfile=None, format='stationxml', sid=None, level='station')
    select(sid)


    """

    def __init__(self, mdfile=None, format='stationxml', sid=None,
                 level='station'):
        """
        Initialize data attribute as an FDSN-compliant dictionary.

        """
        self.data = initialize_metadata(level=level)
        if mdfile is not None:
            self.read(mdfile, format=format, sid=sid, level=level)


    def read(self, mdfile=None, format='stationxml', sid=None, time=None,
             level='station'):
        """
        Populate the data attribute using information from the input mdfile, 
        up to the required level of information. If a stream identifier (sid) 
        is specified, only its pertaining metadata is kept and a specific time
        can be requested.

        NOTE: currently implemented only for StationXML format

        :param [str, binary] mdfile:
            String containing a path or open file-like object containing 
            input metadata information; optional (default is ``None``)

        :param [str] format:
            Format of the input metadata, used to choose the called parser;
            optional (default is 'stationxml')

        :param [str, dict, StreamId] sid: 
            Network.code, Station.code, Channel.locationCode and Channel.code 
            for which metadata should be loaded. Can be a ``StreamId`` object, 
            an FDSN nslc string, or a dictionary (must include keys 'network', 
            'station', 'channel' and 'location'); optional (default is ``None``)

        :param [str, Date] time: 
            The operating time for which information is required for the 
            selected stream, either a ``Date`` instance or an ISO8601 string 
            (default is ``None``)

        :param [str] level: 
            The level of detail to read from input file, either 'network', 'station', 'channel' or 'response'; optional (default is 'station')

        """
        self.data = load_metadata(mdfile, format=format, sid=sid, time=time,
                                  level=level)

    def write(self):
        pass

    def add(self):
        pass

    def remove(self):
        pass

    def select(self, sid, time=None): 
        """
        Return a new ``Metadata`` instance whose data attribute is a pruned 
        copy of the original data containing information pertaining only to the 
        specified stream identifier (sid) at the selected time.

        NOTE: in order for the 'time' option to work, self.data must include 
        information at least at 'channel' level (no automatic check is made!)

        :param [str, dict, list, tuple, StreamId] sid: 
            Network.code, Station.code, Channel.locationCode and Channel.code 
            for which metadata should be kept. 

        :param [str, Date] time: 
            The operating time for which information is required for the 
            selected stream, either a ``Date`` instance or an ISO8601 string 
            (default is ``None``)    
        """
        pruned = self.__class__()

        try:
            md_data = deepcopy(self.data)
            pruned.data = select_nslc(md_data, sid, time)

        except:
            warnings.warn("Selection of metadata for '%s' failed " % sid,
                          UserWarning)
        return pruned


    def __getitem__(self, key):
        """
        Return either a dictionary or a list of dictionaries taken from a copy 
        of the instance data. The level of information is specified by the key. 
        The search is hierarchical and child levels (e.g. 'Channel') are only extracted if the parent levels are available and unique.
    
        :param [str] key:
            Level of information to parse input metadata for. It can be either
            'Network', 'Station', 'Location' or 'Channel'. 

        :return [dict, list] md_item:
            Either a list of dictionaries, or a single dictionary (it the list would contain a single item). Note that the dictionaries follow the FDSN convention only for the corresponding level of information (they do not keep the information related to parent levels)

        """

        try:
            if key in ['Network', 'Station', 'Channel', 'Response']:
                md_item = deepcopy(self.data['Network'])
                if len(md_item) == 1:
                    md_item = md_item[0]
                    if key in ['Station', 'Channel', 'Response']:
                        md_item = md_item['Station']
                        if len(md_item) == 1:
                            md_item = md_item[0]
                            if key in ['Channel', 'Response']:
                                md_item = md_item['Channel']
                                if len(md_item) == 1:
                                    md_item = md_item[0]
                                    if key in ['Response']:
                                        md_item = md_item['Response']
        except:
            warnings.warn("An exception in slicing has occurred", UserWarning)
            return None
        return md_item


    def populate_response(self):
        """
        TODO; extract only the instrumental response from metadata and store it 
        in a useful container (either a dictionary or a dedicated Class; must 
        work for any metadata format)
        """

        paz = deepcopy(paz_map)
        pzstage = self.data['Channel'][0]['Response']['Stage'][0]['PolesZeros']
        paz['input_units'] = pzstage['InputUnits']['Name']
        paz['output_units'] = pzstage['OutputUnits']['Name']
        paz['normalization_factor'] = pzstage['NormalizationFactor']
        paz['normalization_frequency'] = \
            pzstage['NormalizationFrequency']['value']

        poles = pzstage['Pole']
        paz['poles'] = \
                [p['Real']['value']+p['Imaginary']['value']*1j for p in poles]
        zeros = pzstage['Zero']
        paz['zeros'] = \
                [z['Real']['value']+z['Imaginary']['value']*1j for z in zeros]
        return paz


def find_lst_index(lst, key, value):
    """Utility for dictionary parsing"""
    for i, dic in enumerate(lst):
        if dic[key] == value:
            return i
    raise ValueError


def find_lst_indexes(lst, key, value, time=None):
    """Utility for dictionary parsing at Channel level"""
    try:
        indexes = []
        for i, dic in enumerate(lst):
            if dic[key] == value:
                if time:
                    if not isinstance(time, Date):
                        time = Date(time)
                    st_time = dic['startDate']
                    try:
                        end_time = dic['endDate']
                    except:
                        end_time = None 
                    if (time<=st_time):
                        continue
                    if end_time and (time>=end_time):
                        continue
                indexes.append(i)
        if not indexes:
            raise ValueError
        return indexes
    except:
        raise ValueError


def load_metadata(mdfile, format='stationxml', sid=None, time=None, 
                  level='station'):
    """
    Load metadata information at the required level from an input mdfile 
    (either file or memory) into an FDSN-compliant dictionary. Information can 
    be requested for a specific stream identifier (sid) at a specific epoch
    (time).
    (cf. ``shakelab.signals.metadatatemplate.initialize_metadata``)

    NOTE: currently implemented only for StationXML format

    :param [str, binary] mdfile:
        String containing a path or open file-like object containing 
        input metadata information

    :param [str] format:
        Format of the input metadata; optional (default is 'stationxml')

    :param [str, dict, list, tuple, StreamId] sid: 
        Network.code, Station.code, Channel.locationCode and Channel.code for 
        which metadata should be loaded; optional

    :param [str, Date] time: 
        The operating time for which information is required for the selected
        stream, either a ``Date`` instance or an ISO8601 string (default is 
        ``None``)        

    :param [str] level: 
        The level of detail to read from input file, either 'network', 'station', 'channel' or 'response' (default is 'station')

    :return [dict] mdata:
        Dictionary compliant to the FDSN metadata standard

    :raises UserWarning: if format is not implemented 
    """

    if format == 'stationxml':
        try:
            mdata = read_stationxml(mdfile, level=level)
            if sid is not None:
                mdata = select_nslc(mdata, sid, time)
            return mdata
        except:
            warnings.warn("Warning: read_stationxml failed ",  UserWarning)
            pass
    else:
        warnings.warn("Warning: chosen format is not implemented yet; "
                      "available formats are: 'stationxml' ",  UserWarning)


def select_nslc(metadata, sid, time=None):
    """
    Filter a metadata dictionary according to the requested stream identifier 
    (net, sta, loc, cha). If time is specified, only Channels with 
    startDate <= time <= endTime are selected, to provide the correct instrumental response at the required epoch. Returns a dictionary.

    :param [dict] metadata:
        Dictionary compliant to the FDSN metadata standard; can be a 
        ``Metadata.data`` instance    
    
    :param [str, dict, list, tuple, StreamId] sid: 
        Network.code, Station.code, Channel.locationCode and Channel.code for 
        which metadata should be loaded. Can be a ``StreamId`` object, an FDSN 
        nslc string ('network.station.location.channel'), a list, a tuple or a 
        dictionary (must include keys 'network', 'station', 'channel' and 
        'location')
    
    :param [str, Date] time: 
        The operating time for which information is required for the selected
        stream, either a ``Date`` instance or an ISO8601 string (default is 
        ``None``)    

    :return [dict] filtered_md:
        Dictionary compliant to the FDSN metadata standard    
        
    NOTE: the request parameters must follow the nslc hyerarchical order (e.g., 
    a station can only be requested if also a network value is specified); time 
    requests only work if level is at least 'channel'
    """

    filtered_md = deepcopy(metadata)
    try:
        if not isinstance(sid, base.StreamId):
            sid = base.StreamId(sid)
        net = sid.network
        sta = sid.station
        loc = sid.location
        cha = sid.channel

    except:
        warnings.warn("Selection of metadata for '%s' failed " % sid, 
                      UserWarning)

    if not net:
        warnings.warn("No selection specified - copying input metadata", 
                      UserWarning)
    if net:
        try:
            md_nets = filtered_md['Network']
            n_index = find_lst_index(md_nets, 'code', net)
            selected_net = md_nets[n_index]
            filtered_md['Network'] = [selected_net]

            if sta:
                md_stas = selected_net['Station']
                s_index = find_lst_index(md_stas, 'code', sta)
                selected_sta = md_stas[s_index]
                filtered_md['Network'][0]['Station'] = [selected_sta]

                if cha:
                    md_chans = selected_sta['Channel']
                    c_index = find_lst_indexes(md_chans, 'code', cha, time=time)
                    selected_chan = [md_chans[i] for i in c_index]
                    filtered_md['Network'][0]['Station'][0]['Channel'] = \
                        selected_chan

                    if loc:
                        l_index = find_lst_indexes(selected_chan, 
                                                   'locationCode', loc, 
                                                   time=time)
                        selected_loc = [selected_chan[i] for i in l_index] 
                        filtered_md['Network'][0]['Station'][0]['Channel'] = \
                            selected_loc

        except:
            warnings.warn("Warning: Information for requested nslc combination" 
                          " was not found in parsed metadata",  UserWarning)
            return None 

    return filtered_md

