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
Module for earthquake catalogue storage and manipulation.
"""

import numpy as np

from shakelab.libutils.time import Date
from shakelab.libutils.geodetic import WgsPoint

class MagnitudeSolution(object):
    """
    """

    KEY_LIST = ['MagSize', 'MagError', 'MagType', 'MagPrime']

    def __init__(self, size=None, error=None, type=None,
                       code=None, prime=False):

        self.import_data([size, error, type, code], prime)

    def __getitem__(self, key):

        if key == 'MagSize':
            return self.size
        elif key == 'MagError':
            return self.error
        elif key == 'MagType':
            return self.type
        elif key == 'MagCode':
            return self.code
        elif key == 'MagPrime':
            return self.prime
        else:
            raise KeyError('{0}'.format(key))

    def __setitem__(self, key, value):

        if key == 'MagSize':
            self.size = cast_value(value, float)
        elif key == 'MagError':
            self.error = cast_value(value, float)
        elif key == 'MagType':
            self.type = cast_value(value, str)
        elif key == 'MagCode':
            self.code = cast_value(value, str)
        elif key == 'MagPrime':
            self.code = cast_value(value, bool)
        else:
            raise KeyError('{0}'.format(key))

    def __str__(self):
        msg = ''
        for key in self.keys:
            msg += '\t{0}: {1}\n'.format(key, self[key])
        return msg

    def import_data(self, data, prime=False):
        """
        """
        if isinstance(data, dict):
            for (key, value) in data.items():
                self[key] = value
        elif isinstance(data, list):
            self.size = cast_value(data[0], float)
            self.error = cast_value(data[1], float)
            self.type = cast_value(data[2], str)
            self.code = cast_value(data[3], str)
        else:
            raise ValueError('not a valid data format')

        self.prime = cast_value(prime, bool)

    @property
    def keys(self):
        return self.KEY_LIST

class LocationSolution(object):
    """
    """

    KEY_LIST = ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second',
                'Latitude', 'Longitude', 'Depth',
                'SecError', 'LatError', 'LonError', 'DepError',
                'LocCode', 'LocPrime']

    def __init__(self, year=None, month=None, day=None,
                       hour=None, minute=None, second=None,
                       latitude=None, longitude=None, depth=None,
                       secerr=None, laterr=None, lonerr=None, deperr=None,
                       code=None, prime=False):

        self.date = Date()
        self.hypo = WgsPoint(None, None)

        self.import_data([year, month, day, hour, minute, second,
                          latitude, longitude, depth,
                          secerr, laterr, lonerr, deperr, code],
                          prime)

    def __getitem__(self, key):

        if key == 'Year':
            return self.date.year
        elif key == 'Month':
            return self.date.month
        elif key == 'Day':
            return self.date.day
        elif key == 'Hour':
            return self.date.hour
        elif key == 'Minute':
            return self.date.minute
        elif key == 'Second':
            return self.date.second
        elif key == 'SecError':
            return self.date.error
        elif key == 'Latitude':
            return self.hypo.latitude
        elif key == 'Longitude':
            return self.hypo.longitude
        elif key == 'Depth':
            return self.hypo.elevation
        elif key == 'LatError':
            return self.hypo.laterror
        elif key == 'LonError':
            return self.hypo.lonerror
        elif key == 'DepError':
            return self.hypo.eleerror
        elif key == 'LocCode':
            return self.code
        elif key == 'LocPrime':
            return self.prime
        else:
            raise KeyError('{0}'.format(key))

    def __setitem__(self, key, value):

        if key == 'Year':
            self.date.year = cast_value(value, int)
        elif key == 'Month':
            self.date.month = cast_value(value, int)
        elif key == 'Day':
            self.date.day = cast_value(value, int)
        elif key == 'Hour':
            self.date.hour = cast_value(value, int)
        elif key == 'Minute':
            self.date.minute = cast_value(value, int)
        elif key == 'Second':
            self.date.second = cast_value(value, float)
        elif key == 'SecError':
            self.date.error = cast_value(value, float)
        elif key == 'Latitude':
            self.hypo.latitude = cast_value(value, float)
        elif key == 'Longitude':
            self.hypo.longitude = cast_value(value, float)
        elif key == 'Depth':
            self.hypo.elevation = cast_value(value, float)
        elif key == 'LatError':
            self.hypo.laterror = cast_value(value, float)
        elif key == 'LonError':
            self.hypo.lonerror = cast_value(value, float)
        elif key == 'deperror':
            self.hypo.deperror = cast_value(value, float)
        elif key == 'LocCode':
            self.code = cast_value(value, str)
        elif key == 'LocPrime':
            self.prime = cast_value(value, bool)
        else:
            raise KeyError('{0}'.format(key))

    def __str__(self):
        msg = ''
        for key in self.keys:
            msg += '\t{0}: {1}\n'.format(key, self[key])
        return msg

    def import_data(self, data, prime=False):
        """
        """
        if isinstance(data, dict):
            for (key, value) in data.items():
                self[key] = value

        elif isinstance(data, list):
            self.date.year = cast_value(data[0], int)
            self.date.month = cast_value(data[1], int, 1)
            self.date.day = cast_value(data[2], int, 1)
            self.date.hour = cast_value(data[3], int, 0)
            self.date.minute = cast_value(data[4], int, 0)
            self.date.second = cast_value(data[5], float, 0.0)
            self.date.error = cast_value(data[6], float)
            self.hypo.latitude = cast_value(data[7], float)
            self.hypo.longitude = cast_value(data[8], float)
            self.hypo.elevation = cast_value(data[9], float)
            self.hypo.laterror = cast_value(data[10], float)
            self.hypo.lonerror = cast_value(data[11], float)
            self.hypo.eleerror = cast_value(data[12], float)
            self.code = cast_value(data[13], str)

        else:
            raise ValueError('not a valid data format')

        self.prime = cast_value(prime, bool)

    @property
    def keys(self):
        return self.KEY_LIST

class EventProperty(object):
    """
    Base class for event property (magnitude, solution, ....)
    The base class should not be used directly, as solution type
    is not yet defined.
    """

    SOLUTION_TYPE = None

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
            print('Index larger than available solutions')
            return None

    def __str__(self):
        msg = ''
        for idx, es in enumerate(self.solution):
            msg += 'Solution # {0}:\n'.format(idx)
            msg += es.__str__()
        return msg

    def add_solution(self, solution, prime=False):
        """
        """
        # Reset all prime flags if new prime is given
        if prime:
            for es in self.solution:
                es.prime = False

        if isinstance(solution, (dict, list)):
            es = self.SOLUTION_TYPE()
            es.import_data(solution, prime)
            self.solution.append(es)

        elif isinstance(solution, self.SOLUTION_TYPE):
            solution.prime = prime
            self.solution.append(solution)

        # If no prime is given, assuming last solution is prime
        if not any([s.prime for s in self.solution]):
            self.solution[-1].prime = True

    @property
    def prime(self):
        """
        Return the preferred (prime) magnitude solution.
        If no solution is marked as prime, last solution is provided.
        """
        prime_idx = -1
        for idx, es in enumerate(self.solution):
            if es.prime:
                prime_idx = idx
                return self.solution[prime_idx]
        return None


class Magnitude(EventProperty):
    """
    """
    SOLUTION_TYPE = MagnitudeSolution


class Location(EventProperty):
    """
    """
    SOLUTION_TYPE = LocationSolution


class Event(object):
    """
    """
    def __init__(self, id, magnitude=None, location=None):
        self.id = id
        self.magnitude = Magnitude()
        self.location = Location()

        if magnitude is not None:
            self.add_magnitude(magnitude, prime=True)

        if location is not None:
            self.add_location(location, prime=True)

    def add_magnitude(self, magnitude, prime=False):
        """
        """
        self.magnitude.add_solution(magnitude, prime)

    def add_location(self, location, prime=False):
        """
        """
        self.location.add_solution(location, prime)

    def remove_magnitude(idx):
        """
        """
        pass

    def remove_location(idx):
        """
        """
        pass


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
                msg += '\t{0}: {1}\n'.format(key, cast_value(value, str, ''))

        msg += '\tNumber of events: {0}'.format(len(self))
        return msg

    def __len__(self):
        return len(self.event)

    def add_event(self, event):
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

    def add_magnitude(self, id, magnitude, prime=False):
        """
        """
        idx = self._get_index(id)
        if idx is not None:
            self.event[idx].add_magnitude(magnitude, prime)
        else:
            raise ValueError('id not found')

    def add_location(self, id, location, prime=False):
        """
        """
        idx = self._get_index(id)
        if idx is not None:
            self.event[idx].add_location(location, prime)
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


def cast_value(value, dtype, default=None):
    """
    """
    if value not in [None, 'None', 'none', 'NaN', 'nan', '', []]:
        value =  dtype(value)
    else:
        value = default
    return value


def generate_synthetic_catalogue(aval, bval, mmin, mmax, duration=1.):
    """
    Must define an initial date!

    also, it is probably worth to generalize the MFD as input
    instead of aval, bval.
    """

    # Annual rate
    rate = (10.**aval)*(10.**(-bval*mmin)-10.**(-bval*mmax))
    catlen = int(rate*duration)

    # Magnitude distribution:
    # Derived as the inverse of the 
    # double-truncated cumulative GR
    # See Inverse transform sampling (also known as inversion sampling,
    # the inverse probability integral transform)

    Um = np.random.rand(catlen)
    C = 1.-(10.**(-bval*(mmax-mmin)))
    magnitude = mmin-np.log10(1.-(Um*C))/bval

    # Average interval time obtained
    # from the inverse of the cumulative
    # of the exponential distribution
    # (Poisson assumption)
    # TOCHECK!

    Ut = np.random.rand(catlen)
    dT = -np.log(1-Ut)/catlen
    T = np.cumsum(dT)

    # Seconds per year
    # (Not counting leap second years)
    secyear = duration*365.*24.*60.*60.
    timesec = T*secyear

    return timesec, magnitude
