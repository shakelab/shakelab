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
Time and date functionalities
"""
from decimal import Decimal
import time as systime

# CONSTANTS
MSEC = 60.
HSEC = 3600.
DSEC = 86400.
YDAYS = 365.


class Date(object):
    """
    Note: error is in seconds
    """

    def __init__(self, date=None, error=0., timezone='UTC'):
        """
        """
        self.year = None
        self.month = None
        self.day = None
        self.hour = None
        self.minute = None
        self.second = None
        self.error = error
        self.offset = 0.

        if timezone != 'UTC':
            # self.offset = offset_from_timezone(timezone)
            pass

        if date is not None:
            if date == 'now':
                self.set_date(utc_now())
            else:
                self.set_date(date)

    def __eq__(self, value):
        """
        """
        t0 = self.to_seconds(decimal=True)
        t1 = to_decimal(value)

        return (t0 == t1)

    def __lt__(self, value):
        """
        """
        t0 = self.to_seconds(decimal=True)
        t1 = to_decimal(value)

        return (t0 < t1)

    def __le__(self, value):
        """
        """
        t0 = self.to_seconds(decimal=True)
        t1 = to_decimal(value)

        return (t0 <= t1)

    def __gt__(self, value):
        """
        """
        t0 = self.to_seconds(decimal=True)
        t1 = to_decimal(value)

        return (t0 > t1)

    def __ge__(self, value):
        """
        """
        t0 = self.to_seconds(decimal=True)
        t1 = to_decimal(value)

        return (t0 >= t1)

    def __add__(self, value):
        """
        """
        if isinstance(value, (int, float)):
            t0 = self.to_seconds(decimal=True)
            t1 = Decimal(value)
            return Date(float(t0 + t1))

        elif isinstance(value, Date):
            t0 = self.to_seconds(decimal=True)
            t1 = value.to_seconds(decimal=True)
            return float(t0 + t1)

        else:
            print('Not a supported operation.')
            return None

    def __sub__(self, value):
        """
        """
        if isinstance(value, (int, float)):
            t0 = self.to_seconds(decimal=True)
            t1 = Decimal(value)
            return Date(float(t0 - t1))

        elif isinstance(value, Date):
            t0 = self.to_seconds(decimal=True)
            t1 = value.to_seconds(decimal=True)
            return float(t0 - t1)

        else:
            print('Not a supported operation.')
            return None

    def __repr__(self):
        """
        """
        return self.get_date(dtype='string')

    @property
    def ordinal_day(self):
        """
        """
        return to_ordinal_day(self.year, self.month, self.day)

    @ordinal_day.setter
    def ordinal_day(self, value):
        """
        """
        (self.month, self.day) = from_ordinal_day(self.year, value)

    def set_date(self, date):
        """
        Ordinal format is recognized from 
        """
        if isinstance(date, (list, tuple)):

            if len(date) == 6:
                self.year = int(date[0])
                self.month = int(date[1])
                self.day = int(date[2])
                self.hour = int(date[3])
                self.minute = int(date[4])
                self.second = float(date[5])

            elif len(date) == 5:
                self.year = int(date[0])
                self.day = int(date[1])
                self.hour = int(date[2])
                self.minute = int(date[3])
                self.second = float(date[4])
                (self.month, self.day) = from_ordinal_day(self.year, self.day)

        elif isinstance(date, dict):

            self.year = date['year']
            self.hour = date['hour']
            self.minute = date['minute']
            self.second = date['second']

            if 'month' in dict:
                self.month = date['month']
                self.day = date['day']
            else:
                # Assuming ordinal format
                (self.month, self.day) = from_ordinal_day(self.year, self.day)

        elif isinstance(date, str):

            dbuf = read_iso8601_date(date)
            self.year = dbuf[0]
            self.month = dbuf[1]
            self.day = dbuf[2]
            self.hour = dbuf[3]
            self.minute = dbuf[4]
            self.second = dbuf[5]
            self.offset = dbuf[6]

        elif isinstance(date, (int, float, Decimal)):

            self.from_seconds(date)

        else:
            raise ValueError('Format not recognized')

        self._selfcheck()

    def get_date(self, dtype='list', format='calendar'):
        """
        NOTE: ordinal format not yet implemented in output
        """
        if dtype in ('l', 'list'):

            date = [self.year, self.month, self.day,
                    self.hour, self.minute, self.second]

        elif dtype in ('s', 'str', 'string', 'iso8601'):
            """
            # ISO-8601 format
            date = '{0:04d}-'.format(self.year)
            date += '{0:02d}-'.format(self.month)
            date += '{0:02d}T'.format(self.day)
            date += '{0:02d}:'.format(self.hour)
            date += '{0:02d}:'.format(self.minute)
            date += '{0:07.4f}'.format(self.second)
            """
            date = write_iso8601_date(self.year, self.month, self.day,
                                      self.hour, self.minute, self.second)

        else:
            raise ValueError('Format not recognized')

        return date

    def _selfcheck(self):
        """
        """
        if self.year < 1:
            raise ValueError('Year must be >= 1')

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must be between 1 and 12')

        if self.day < 1 or self.day > 31:
            raise ValueError('Day must be between 1 and 31')

        if self.hour < 0 or self.hour > 23:
            raise ValueError('Hours must be between 0 and 23')

        if self.minute < 0 or self.minute > 59:
            raise ValueError('Minutes must be between 0 and 59')

        if self.second < 0 or self.second > 60:
            raise ValueError('Seconds must be between 0 and 60')

    def to_seconds(self, decimal=False):
        """
        """
        return date_to_second(self.year, self.month, self.day,
                              self.hour, self.minute, self.second,
                              decimal=decimal)

    def from_seconds(self, second):
        """
        """
        self.set_date(second_to_date(float(second)))

    @property
    def seconds(self):
        """
        """
        return self.to_seconds()

    @seconds.setter
    def seconds(self, value):
        """
        """
        self.from_seconds(value)

    def shift_time(self, time, units='s'):
        """
        """
        if units in ['s', 'second', 'seconds']:
            seconds = time
        elif units in ['m', 'minute', 'minutes']:
            seconds = time * MSEC
        elif units in ['h', 'hour', 'hours']:
            seconds = time * HSEC
        elif unit in ['d', 'day', 'days']:
            seconds = time * DSEC
        else:
            print('Time units not yet implemented')

        self.from_seconds(self.to_seconds() + seconds)

def leap_check(year):
    """
    Check if leap year.
    """
    c0 = (year % 4 == 0)
    c1 = (year % 100 != 0)
    c2 = (year % 400 == 0)

    return (c0 and c1) or c2

def leap_num(year):
    """
    Compute the number of leap years before current year
    (starting from 1 a.c.)
    """
    n0 = (year-1)//4
    n1 = (year-1)//100
    n2 = (year-1)//400

    return n0 - n1 + n2

def days_in_month(year):
    """
    """
    if leap_check(year):
        mdays = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    else:
        mdays = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

    return mdays

def date_to_second(year=1, month=1, day=1, hour=0, minute=0, second=0.0,
                decimal=False):
    """
    Convert a date to seconds.
    Not yet implemented for b.c. years.
    """

    if year < 1:
      print('Error: Year must be positive (> 1)')
      return None

    #MDAYS = [0] + days_in_month(2012)[:-1]
    MDAYS = [0] + days_in_month(year)[:-1]    # TO CHECK

    ysec = (year - 1) * YDAYS * DSEC
    ysec += leap_num(year) * DSEC
    msec = MDAYS[int(month) - 1] * DSEC
    dsec = (day - 1) * DSEC

    sec = ysec + msec + dsec + hour*HSEC+ minute*MSEC

    if decimal:
        return Decimal(sec) + Decimal(second)
    else:
        return sec + second

def from_ordinal_day(year, day):
    """
    Converts a ordinal day of the year to the corresponding month and day.

    Parameters:
        year (int): The year to consider for leap year calculations.
        day (int): The ordinal day of the year.

    Returns:
        tuple: A tuple (month, day) where 'month' is the month and 'day'
               is the day of the month.
    """
    MDAYS = days_in_month(year)

    for m, d in enumerate(MDAYS):
        if day <= d:
            month = m + 1
            if m > 0:
                day -= MDAYS[m - 1]
            break

    return (month, day)

def to_ordinal_day(year, month, day):
    """
    Converts a given month and day to the ordinal day of the year.
    
    Parameters:
        year (int): The year to consider for leap year calculations.
        month (int): The month (1 for January, 2 for February, etc.).
        day (int): The day of the month.
    
    Returns:
        int: The ordinal day of the year.
    """
    MDAYS = days_in_month(year)
    
    if month == 1:
        return day
    else:
        return MDAYS[month - 2] + day

def second_to_date(second):
    """
    Implemented using a direct search approach
    (it would be nicer to use a more elengant algorithm).
    """
    # Loop over years (accounting for leap days)
    year = 0
    while True:
        year += 1
        secy = YDAYS * DSEC
        if leap_check(year):
            secy += DSEC
        if second >= secy:
            second -= secy
        else:
            break

    # Number of total days in the last year
    day = (second // DSEC) + 1

    # Find month corresponding to day of the year
    MDAYS = days_in_month(year)

    for m, d in enumerate(MDAYS):
        if day <= d:
            month = m + 1
            if m > 0:
                day -= MDAYS[m - 1]
                second -= MDAYS[m - 1] * DSEC
            break
    second -= (day - 1) * DSEC

    # Hours
    hour = second // HSEC
    second -= hour * HSEC

    # Minutes
    minute = second // MSEC
    second -= minute * MSEC

    return [int(year), int(month), int(day),
            int(hour), int(minute), second]

def read_iso8601_date(date_str):
    """
    Parses a date string in ISO 8601 format and returns individual components.

    Args:
        date_str (str): A string representing a date in ISO 8601 format.

    Returns:
        tuple: A tuple containing the year, month, day, hour, minute, second,
               and time offset (for time zone).

    Note:
        The ISO 8601 format supported by this function is
        'YYYY-MM-DDTHH:MM:SS.sssZ' or 'YYYY-MM-DDTHH:MM:SS.sss±HH:MM',
        where 'Z' indicates UTC and '±HH:MM' represents the time offset
        from UTC.
    """
    [date_part, time_part] = date_str.split('T')

    # Extract date part
    date_len = len(date_part.split('-'))

    if date_len == 3:
        year, month, day = map(int, date_part.split('-'))

    if date_len == 2:
        year, day = map(int, date_part.split('-'))
        (month, day) = from_ordinal_day(year, day)

    time_offset = 0

    # Extract time part
    if 'Z' in time_part:
        time_part = time_part.strip('Z')

    elif '+' in time_part:
        [time_part, offset_part] = time_part.split('+')
        offsets = offset_part.split(':')
        time_offset = int(offsets[0]) * HSEC
        if len(offsets) > 1:
            time_offset = time_offset + int(offsets[1]) * MSEC
        if len(offsets) > 2:
            time_offset = time_offset + int(offsets[2])

    elif '-' in time_part:
        [time_part, offset_part] = time_part.split('-')
        offsets = offset_part.split(':')
        time_offset = - int(offsets[0]) * HSEC
        if len(offsets) > 1:
            time_offset = time_offset - int(offsets[1]) * MSEC
        if len(offsets) > 2:
            time_offset = time_offset - int(offsets[2])

    if '.' in time_part:
        [time_part, msecond] = time_part.split('.')
        mslen = len(msecond)
        msecond = int(msecond)

    else:
        mslen = 0
        msecond = 0

    hour, minute, second = map(int, time_part.split(':'))

    # Add decimal part to seconds
    second += msecond * 10**(-mslen)

    return year, month, day, hour, minute, second, time_offset

def write_iso8601_date(year, month, day, hour, minute, second, timezone='Z'):
    """
    Converts year, month, day, hour, minute, and second to ISO 8601 format.

    Args:
        year (int): The year.
        month (int): The month (1-12).
        day (int): The day of the month (1-31).
        hour (int): The hour (0-23).
        minute (int): The minute (0-59).
        second (float or int): The second (including fractions of a second).
        timezone (str) : Local time zone (Z is UTC default)

    Returns:
        str: ISO 8601 formatted string.
    """
    iso8601_str = "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:09.6f}{:s}".format(
        year, month, day, hour, minute, second, timezone
    )
    return iso8601_str

def utc_now():
    """
    Get the current UTC date and time in ISO 8601 format.

    Returns:
        str: ISO 8601 formatted string of current UTC date and time.
    """
    dbuf = systime.gmtime(systime.time())
    return write_iso8601_date(dbuf[0], dbuf[1], dbuf[2],
                              dbuf[3], dbuf[4], dbuf[5])

def to_decimal(value):
    """
    """
    if isinstance(value, Date):
        t0 = value.to_seconds(decimal=True)
    elif isinstance(value, (int, float)):
        t0 = Decimal(value)
    else:
        raise ValueError('Not a valid time format')

    return t0