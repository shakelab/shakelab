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
Module for earthquake catalogue storage and manipulation.
"""

import numpy as np
import pickle
from copy import deepcopy

from shakelab.libutils.time import Date
from shakelab.libutils.geodetic import WgsPoint
from shakelab.libutils.ascii import AsciiTable
from shakelab.libutils.utils import cast_value


# A key string is associated to an attribute name;
# for each attribute, type and default values are defined
_MAGMAP = {'MagSize' : ('size', float, None),
           'MagError' : ('error', float, None),
           'MagType' : ('type', str, None),
           'MagCode' : ('code', str, None),
           'MagPrime' : ('prime', bool, False)}

_LOCMAP = {'Year' : ('year', int, None),
           'Month' : ('month', int, 1),
           'Day' : ('day', int, 1),
           'Hour' : ('hour', int, 0),
           'Minute' : ('minute', int, 0),
           'Second' : ('second', float, 0.),
           'Latitude' : ('latitude', float, None),
           'Longitude' : ('longitude', float, None),
           'Depth' : ('depth', float, 0.),
           'SecError' : ('secerror', float, None),
           'LatError' : ('laterror', float, None),
           'LonError' : ('lonerror', float, None),
           'DepError' : ('deperror', float, None),
           'LocCode' : ('code', str, None),
           'LocPrime' : ('prime', bool, False)}


class BaseSolution(object):
    """
    Base class.
    Should not be called directly as KEYMAP is empty.
    """

    _KEYMAP = {}

    def __init__(self, data=None):

        # Initialise attributes to default value
        for key in self._KEYMAP:
            name = self._KEYMAP[key][0]
            value = self._KEYMAP[key][2]
            exec('self.{0}={1}[0]'.format(name, [value]))

        if data is not None:
            self.set(data)

    def __setitem__(self, key, value):

        if key in self._KEYMAP:
            name = self._KEYMAP[key][0]
            value = cast_value(value, self._KEYMAP[key][1])
            exec('self.{0}={1}[0]'.format(name, [value]))
        else:
            raise KeyError('{0}'.format(key))

    def __getitem__(self, key):

        if key in self._KEYMAP:
            return eval('self.{0}'.format(self._KEYMAP[key][0]))
        else:
            raise KeyError('{0}'.format(key))

    def __str__(self):
        msg = ''
        for key in self.keys:
            msg += '\t{0}: {1}\n'.format(key, self[key])
        return msg

    def set(self, data):
        """
        """
        if isinstance(data, dict):
            for (key, value) in data.items():
                self[key] = value
        else:
            raise ValueError('not a valid data format')

    def get(self):
        """
        """
        data = {}
        for key in self.keys:
            data[key] = self[key]
        return data

    @property
    def keys(self):
        return self._KEYMAP.keys()


class MagnitudeSolution(BaseSolution):
    """
    """

    _KEYMAP = _MAGMAP


class LocationSolution(BaseSolution):
    """
    """

    _KEYMAP = _LOCMAP

    @property
    def date(self):
        return Date([self.year, self.month, self.day,
                     self.hour, self.minute, self.second])

    @property
    def hypocentre(self):
        return WgsPoint(self.latitude, self.longitude, -self.depth)


class BaseProperty(object):
    """
    Base class for event property (magnitude, location, ....)
    The base class should not be used directly, as solution type
    is not yet defined.
    """

    _SOLUTION_TYPE = None

    def __init__(self, solution=None):
        self.solution = []

        if solution is not None:
            self.add_solution(solution, prime=True)

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            solution = self.solution[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return solution

    def __getitem__(self, idx):
        if idx < len(self.solution):
            return self.solution[idx]
        else:
            raise ValueError('Index larger than available solutions')

    def __len__ (self):
        return len(self.solution)

    def __str__(self):
        msg = ''
        for idx, es in enumerate(self.solution):
            msg += 'Solution # {0}:\n'.format(idx)
            msg += es.__str__()
        return msg

    def add(self, solution, prime=None):
        """
        """

        if isinstance(solution, dict):
            solution = self._SOLUTION_TYPE(solution)

        elif isinstance(solution, self._SOLUTION_TYPE):
            pass

        else:
            raise ValueError('not a valid solution format')

        if prime is not None:
            solution.prime = prime

        # Reset all prime flags if new prime is given
        if solution.prime:
            for es in self.solution:
                es.prime = False

        self.solution.append(solution)

        # If no prime is given, assuming last solution is prime
        #if not any([s.prime for s in self.solution]):
        #    self.solution[-1].prime = True

    def remove(self, idx=-1):

        if idx < len(self.solution):
            del self.solution[idx]
        else:
            raise ValueError('Index larger than available solutions')

    @property
    def prime(self):
        """
        Return the preferred (prime) magnitude solution.
        If no solution is marked as prime, last solution is provided.
        """
        if self.solution:
            prime_idx = -1
            for idx, es in enumerate(self.solution):
                if es.prime:
                    prime_idx = idx
            return self.solution[prime_idx]
        else:
            return None


class Magnitude(BaseProperty):
    """
    """
    _SOLUTION_TYPE = MagnitudeSolution


class Location(BaseProperty):
    """
    """
    _SOLUTION_TYPE = LocationSolution


class Event(object):
    """
    """
    def __init__(self, id, magnitude=None, location=None):
        self.id = id
        self.magnitude = Magnitude()
        self.location = Location()

        if magnitude is not None:
            self.add_magnitude(magnitude)

        if location is not None:
            self.add_location(location)

    def __getitem__(self, key):
        if key == 'Id':
            return self.id
        elif key == 'Magnitude':
            return self.magnitude
        elif key == 'Location':
            return self.location

    def __str__(self):
        msg = 'EVENT {0}\n'.format(self.id)
        msg += 'MAGNITUDE SOLUTIONS:\n'
        msg += self.magnitude.__str__()
        msg += 'LOCATION SOLUTIONS:\n'
        msg += self.location.__str__()
        return msg

    def add_magnitude(self, magnitude, prime=None):
        """
        """
        self.magnitude.add(magnitude, prime=prime)

    def add_location(self, location, prime=None):
        """
        """
        self.location.add(location, prime=prime)

    def remove_magnitude(idx=-1):
        """
        """
        self.magnitude.remove(idx)

    def remove_location(idx=-1):
        """
        """
        self.location.remove(idx)

    def copy(self):
        """
        """
        return deepcopy(self)


class EqDatabase(object):
    """
    """

    def __init__(self, name=None, version=None, info=None):
        self.header = {'Name': name,
                       'Version': version,
                       'Info': info}
        self.event = []

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        try:
            event = self.event[self._counter]
            self._counter += 1
        except IndexError:
            raise StopIteration
        return event

    def __getitem__(self, idx):
        return self.event[idx]

    def __str__(self):
        msg = ''
        for (key, value) in self.header.items():
            if value is not None:
                msg += '{0}: {1}\n'.format(key, cast_value(value, str, ''))

        msg += 'Number of events: {0}\n'.format(len(self))

        min, max = self.get_range('Year')
        msg += 'Year rage: ({0}, {1})\n'.format(min, max)

        min, max = self.get_range('MagSize')
        msg += 'Magnitude rage: ({0}, {1})\n'.format(min, max)

        min, max = self.get_range('Depth')
        msg += 'Depth rage: ({0}, {1})\n'.format(min, max)

        return msg

    def __len__(self):
        return len(self.event)

    def add(self, event):
        """
        """
        if not isinstance(event, Event):
            raise ValueError('not a valid event')

        # Check if the event already exists
        idx = self._get_index(event.id)
        
        if idx is not None:

            # Reset previous prime flags
            for ms in self.event[idx].magnitude:
                ms.prime = False
            for ls in self.event[idx].location:
                ls.prime = False

            # Append new solution(s)
            for ms in event.magnitude:
                self.event[idx].add_magnitude(ms)
            for ls in event.location:
                self.event[idx].add_location(ls)
        else:
            self.event.append(event)

    def remove(self, id):
        """
        """
        idx = self._get_index(id)
        if idx is not None:
            del self.event[idx]
        else:
            raise ValueError('id not found')

    def add_magnitude(self, id, magnitude, prime=None):
        """
        """
        idx = self._get_index(id)
        if idx is not None:
            self.event[idx].add_magnitude(magnitude, prime=prime)
        else:
            raise ValueError('id not found')

    def add_location(self, id, location, prime=None):
        """
        """
        idx = self._get_index(id)
        if idx is not None:
            self.event[idx].add_location(location, prime=prime)
        else:
            raise ValueError('id not found')

    def _get_index(self, id):
        """
        """
        for idx in range(len(self.event)):
            if id == self.event[idx].id:
                return idx
        return None

    def get_event(self, id):
        """
        """
        for event in self.event:
            if id == event.id:
                return event
        return None

    def filter(self, key, operator, value, delete_empty=True):
        """
        Note: this method will filter solutions in place.
        None items will be filtered out
        """
        if key in _LOCMAP:
            solution_key = 'Location'
        elif key in _MAGMAP:
            solution_key = 'Magnitude'
        else:
            raise ValueError('not a valid key')

        if operator not in ['eq', 'ne', 'gt', 'lt', 'ge', 'le']:
            raise ValueError('not a valid operator')

        ide_list = []

        for ide, event in enumerate(self.event):

            idx_list = []
            event_solution = event[solution_key]

            for idx, solution in enumerate(event_solution):
                target = solution[key]

                if target is None:
                    idx_list.append(idx)
                else:
                    if operator == 'eq' and target != value:
                        idx_list.append(idx)

                    elif operator == 'ne' and target == value:
                        idx_list.append(idx)

                    elif operator == 'gt' and target <= value:
                        idx_list.append(idx)

                    elif operator == 'lt' and target >= value:
                        idx_list.append(idx)

                    elif operator == 'ge' and target < value:
                        idx_list.append(idx)

                    elif operator == 'le' and target > value:
                        idx_list.append(idx)

            # Delete empty solutions in reverse order
            idx_list.sort(reverse=True)

            for idx in idx_list:
                del event_solution.solution[idx]

            if len(event_solution) == 0:
                ide_list.append(ide)

        if delete_empty:
            # Delete empty events in reverse order
            ide_list.sort(reverse=True)

            for ide in ide_list:
                del self.event[ide]

    def list_values(self, key, remove_empty=True, any=False):
        """
        """
        values = []

        if key.lower() == 'id':
            return [e['Id'] for e in self]

        elif key.lower() == 'date':
            return [e.location.prime.date for e in self]

        elif key in _LOCMAP:
            for e in self:
                if e.location.prime:
                    values.append(e.location.prime[key])
                else:
                    values.append(None)

        elif key in _MAGMAP:
            for e in self:
                if e.magnitude.prime:
                    values.append(e.magnitude.prime[key])
                else:
                    values.append(None)

        else:
            raise ValueError('not a valid key')

        if remove_empty:
            values = [v for v in values if v is not None]

        return values

    def get_range(self, key):
        """
        """
        values = self.list_values(key, remove_empty=True)
        return [min(values), max(values)]

    def sort(self, key='Date', reverse=False):
        """
        Sort the database events by a given key.
    
        Parameters
        ----------
        key : str
            Field name to sort by. It can be:
            - 'Id'
            - 'Date' (prime location date)
            - Any key from location (_LOCMAP)
            - Any key from magnitude (_MAGMAP)
        reverse : bool, default False
            If True, sort in descending order.
        """
        if key.lower() == 'id':
            sort_values = [e.id for e in self.event]
    
        elif key.lower() == 'date':
            # Sort by prime location date (converted to seconds)
            sort_values = [
                e.location.prime.date.to_seconds()
                if e.location.prime else None
                for e in self.event
            ]
    
        elif key in _LOCMAP:
            sort_values = [
                e.location.prime[key] if e.location.prime else None
                for e in self.event
            ]
    
        elif key in _MAGMAP:
            sort_values = [
                e.magnitude.prime[key] if e.magnitude.prime else None
                for e in self.event
            ]
    
        else:
            raise ValueError(f'Not a valid key: {key}')
    
        # Replace None with a placeholder that won't break sorting
        sort_values = [
            v if v is not None else float('-inf') for v in sort_values
        ]
    
        idx = sorted(range(len(sort_values)),
                     key=lambda i: sort_values[i],
                     reverse=reverse)
    
        self.event = [self.event[i] for i in idx]

    def sort_by_date(self, reverse=False):
        """
        Sort events by prime location date.
        Wrapper for sort_by('Date').
        """
        self.sort('Date', reverse=reverse)

    def sort_by_magnitude(self, reverse=True):
        """
        Sort events by prime magnitude size.
        Wrapper for sort('MagSize').
    
        Parameters
        ----------
        reverse : bool, default False
            If True, sort in descending order (largest first).
        """
        self.sort('MagSize', reverse=reverse)

    def append(self, db):
        """
        Append another EqDatabase to the current one, avoiding duplicate ids.
    
        Events with IDs already present in the current catalogue are
        skipped. These skipped events are collected into a new EqDatabase,
        which is returned to the caller.
    
        Parameters
        ----------
        db : EqDatabase
            The catalogue to append.
    
        Returns
        -------
        EqDatabase
            A new EqDatabase containing only the events skipped because
            their IDs were already present in this catalogue.
        """
        if not isinstance(db, EqDatabase):
            raise ValueError("Input must be an EqDatabase")
    
        duplicates_db = EqDatabase(name="Duplicate Events")
    
        for event in db:
            if self._get_index(event.id) is not None:
                duplicates_db.add(event.copy())
                continue
            self.add(event.copy())
    
        if len(duplicates_db) > 0:
            dup_ids = [ev.id for ev in duplicates_db]
            print(f"[append] Skipped {len(duplicates_db)} duplicate events")
    
        return duplicates_db

    def read(self):
        """
        """
        pass

    def write(self):
        """
        """
        pass

    def load(self, file_name):
        """
        Load the database from pickle file
        """
        with open(file_name, 'rb') as f:
            buf = pickle.load(f)
            self.header = buf.header
            self.event = buf.event
            f.close()
            return

    def dump(self, file_name):
        """
        Save the database to a pickle file
        """
        with open(file_name, 'wb') as f:
            pickle.dump(self, f, protocol=2)
            f.close()
            return

    def copy(self):
        """
        """
        return deepcopy(self)