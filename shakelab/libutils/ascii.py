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
A tool to manipulate tabular data from/to ascii files.
"""

from shakelab.libutils.utils import cast_value


class AsciiTable():

    def __init__(self, header=[]):

        self.header = header
        self.database = []

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            item = self.database[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return item

    def add_element(self, data):
        """
        Add a new data element (row) to the database table.
        Data must be in header's format.
        """

        newitem = {}
        for i, key in enumerate(self.header):
            newitem[key] = data[i]
        self.database.append(newitem)

    def add_key(self, key, index=-1, data=None):
        """
        Add an header key at given position (default is last element).
        Data structure can optionally be inflated with a scalar or a
        list values of values (default is None).
        """

        # Check if key is already stored
        if key in self.header:
            print('Warning: key already in header.')
            return

        # Add key to header list (last element by default)
        if index == -1:
            index = len(self.header)
        self.header.insert(index, key)

        # Loop over data
        for i, item in enumerate(self.database):

            # Check value types
            if isinstance(data, (str, int, float)):
                element = data
            else:
                element = data[i]

            # Add element at corresponding key
            self.database[i][key] = element

    def remove_key(self, key):
        """
        Remove a given key from header and data structure.
        """

        # Remove key from header
        self.header.pop(self.header.index(key))

        # Remove elements from data
        for i, item in enumerate(self.database):
            self.database[i].pop(key)

    def rename_key(self, old_key, new_key):
        """
        Rename a given key in the header and the data structure.
        """

        # Rename header's key
        self.header[self.header.index(old_key)] = new_key

        # Rename key in data structure
        for i, item in enumerate(self.database):
            self.database[i][new_key] = self.data[i].pop(old_key)

    def extract(self, key, dtype=float):
        """
        Extract data values by key.
        Data type can be specified.
        """

        return [cast_value(item[key], dtype) for item in self.database]

    @property
    def size(self):
        """
        Return size of the database.
        """

        enum = len(self.database)
        hnum = len(self.header)

        return enum, hnum

    def concatenate(self, table):
        """
        Merge two data tables.
        Header structure must be identical.
        """

        if self.header == table.header:
            for i in range(0, table.size()[0]):
                self.data.append(table.data[i])
        else:
            print('Error: headers do not match...')

    def read(self, ascii_file, header=None, dtype=str, delimiter=',',
             skipline=0, comment='#', empty=None):
        """
        Import data from ascii file (tabular)
        """

        self.header = []
        self.database = []

        # Open input ascii file
        with open(ascii_file, 'r') as f:

            # Ignore initial line(s) if necessary
            for i in range(0, skipline):
                f.readline()

            # Import header (skip comments)
            if header is None:
                while 1:
                    line = f.readline()
                    if line[0] != comment:
                        break
                header = line.strip().split(delimiter)
            self.header = header

            # Loop over lines
            for line in f:

                # Skip comments, if any
                if line[0] != comment:
                    value = line.strip().split(delimiter)

                    # Loop over data values
                    data = []
                    for i, h in enumerate(header):

                        # Skip empty header fields
                        if h not in ['', None]:

                            # Data type(s) switch
                            if isinstance(dtype, list):
                                dtp = dtype[i]
                            else:
                                dtp = dtype

                            data.append(cast_value(value[i], dtp, empty))
                    self.add_element(data)
            return

        # Warn user if model file does not exist
        print('Error: file not found.')

    def write(self, ascii_file, header=True, delimiter=','):
        """
        Export data object into an ascii file.
        """
        with open(ascii_file, 'w') as f:

            # Write header
            if header:
                f.write(delimiter.join(self.header) + '\n')

            # Write data (loop over rows)
            for i, item in enumerate(self.database):
                data = [cast_value(item[j], str) for j in self.header]
                data = delimiter.join(data)

                if i < (self.size()[0] - 1):
                    f.write(data + '\n')
                else:
                    f.write(data)
            return

        # Warn user if model file does not exist
        print('Cannot open file.')
