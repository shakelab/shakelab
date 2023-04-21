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
Seedlink client implementation
"""

import socket
from threading import Thread, Event
import xml.etree.ElementTree as et

import matplotlib.pyplot as plt

from shakelab.signals.mseed import ByteStream, Record

SL_DEFAULT_PORT = 18000
BUFFER_SIZE = 1024


class Client():

    def __init__(self, host=None, port=SL_DEFAULT_PORT):
        self._s = None
        self.host = host
        self.port = port
        self.id = []

        if self.host is not None:
            self.connect(host, port)

    def _send(self, string):
        """
        """
        self._s.send(bytes('{0}\r'.format(string), 'utf8'))

    def _recv(self, buffer_size=BUFFER_SIZE):
        """
        """
        return self._s.recv(buffer_size).decode()

    def connect(self, host, port=SL_DEFAULT_PORT):
        """
        """
        try:
            self._s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

            self._s.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 65536)
            self._s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self._s.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)

            self._s.connect((host, port))
            self._s.setblocking(1)

            self._send('HELLO')
            response = self._recv()
            print(response)

        except:
            print('Error')

    def info(self, level):
        """
        Requests an INFO packet.
        Level should be one of the following:
            ID, CAPABILITIES, STATIONS, STREAMS, GAPS, CONNECTIONS, ALL
        """
        levels = ['ID', 'CAPABILITIES', 'STATIONS', 'STREAMS',
                  'GAPS', 'CONNECTIONS', 'ALL']

        if level not in levels:
            print('Error')

        else:
            self._send('INFO ' + level)

            info = ""
            while True:
                header = self._s.recv(8).decode()

                data = b''
                while len(data) < 512:
                    data += self._s.recv(8)

                byte_stream = ByteStream(data=data, byte_order='be')
                record = Record(byte_stream)
                info += record.decode()

                # to check
                if '*' not in header: break

            self.id = []
            root = et.fromstring(info)
            for child in root:
                if child.tag == 'station':
                   self.id.append(child.attrib)
                else:
                   print('Not identified tag')

            return info


    def station(self, station_code, network_code=''):
        """
        Turns on multi-station mode, used to transfer data of
        multiple stations over a single TCP channel.
        """
        command = 'STATION {0} {1}\r'.format(station_code, network_code)
        self._s.send(bytes(command, 'utf8'))
        response = self._s.recv(BUFFER_SIZE).decode()
        print(repr(response.strip()))


    def select(self, pattern=''):
        """
        """
        command = 'SELECT {0}\r'.format(pattern)
        print(command)
        self._s.send(bytes(command, 'utf8'))
        response = self._s.recv(BUFFER_SIZE).decode()
        print(repr(response.strip()))

    def start_thread(self):
        """
        """
        event = Event()
        self._t = Thread(target=self.get_data)
        self._t.start()

    def close_thread(self):
        self._t.join()

    def get_data(self):
        """
        """
        

    def start(self):
        """
        """
        self._s.send(b'END\r')

        for i in range(1024):

            header = self._s.recv(8).decode()
            print('Received:', repr(header.strip()))

            # Ensuring that received data is exactly 512 bytes
            # In UNIX is equivalent to: s.recv(512, socket.MSG_WAITALL)
            data = b''
            while len(data) < 512:
                data += self._s.recv(8)
                if not data: break

            byte_stream = ByteStream(data=data, byte_order='be')
            record = Record(byte_stream)

            #print(record.header)
            #print(record.blockette)
            
            plt.figure(1)
            plt.plot(record.data)
            plt.show(block=False)


    def close(self):
        """
        Closes the connection to seedlink server.
        """
        try:
            self._s.send(b'BYE\r')
            self._s.close()
        except:
            print('No connection to close')




class RingBuffer():
    """
    Modified from the example of Sébastien Keim.
    """
    
    def __init__(self, size):
        self.size = size
        self.data = [None] * size
        self._i0 = 0

    def append_value(self, data):
        """
        Append an element overwriting the oldest one.
        """
        self.data[self._i0] = data
        self._i0 = (self._i0+1) % self.size

    def append_array(self, data):
        """
        """
        for d in data:
            self.append_value(d)

    def get_buffer(self):
        """
        Return list of elements in correct order
        """
        return self.data[self._i0:]+self.data[:self._i0]

