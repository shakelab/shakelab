# ============================================================================
#
# Copyright (C) 2019 Valerio Poggi.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ============================================================================
"""
Time and date functionalities
"""

from datetime import datetime

class Date(object):
    """
    """

    def __init__(self, date=None):
        if date is not None:
            self.set(date)

    def set(self, date):
        """
        """
        self.year = int(date[0])
        self.month = int(date[1])
        self.day = int(date[2])
        self.hour = int(date[3])
        self.minute = int(date[4])
        self.second = float(date[5])

    def get(self):
        """
        """
        return [self.year, self.month, self.day,
                self.hour, self.minute, self.second]

    def to_second(self):
        """
        """
        return date_to_sec(self.year, self.month, self.day,
                           self.hour, self.minute, self.second)

    def from_second(self, second):
        """
        """
        self.set(sec_to_date(second))

    def shift_time(self, second, unit='s'):
        """
        """
        if unit is 'm':
            second *= 60
        if unit is 'h':
            second *= 3600
        if unit is 'd':
            second *= 86400

        self.from_second(self.to_second() + second)

    def __add__(self, value):
        """
        """
        if isinstance(value, (int, float)):
            dnew = Date()
            dnew.from_second(self.to_second() + value)
            return dnew
        else:
            print('Not a supported operation')

    def __sub__(self, value):
        """
        """
        if isinstance(value, (int, float)):
            dnew = Date()
            dnew.from_second(self.to_second() - value)
            return dnew
        elif isinstance(value, Date):
            return self.to_second() - value.to_second()
        else:
            print('Not a supported operation')

    def __repr__(self):
        """
        """
        return repr(self.get())


def date_to_sec(year, month, day, hour, minute, second):
    """
    """
    isecond = int(second // 1)
    msecond = int((second % 1) * 1e6)

    t0 = datetime(1 ,1 ,1 ,0 ,0 ,0 ,0)
    t1 = datetime(year, month, day, hour, minute, isecond, msecond)
    return (t1 - t0).total_seconds()

def sec_to_date(second=0.):
    """
    NOTE: Datetime produces rounding errors for unknown reason.
    """
    sref = date_to_sec(1970, 1, 1, 0, 0, 0.)
    d = datetime.utcfromtimestamp(second - sref)
    fsecond = d.second + (d.microsecond * 1e-6)
    return [d.year, d.month, d.day, d.hour, d.minute, fsecond]

# -----------------------------------------------------------------------------
# MY IMPLEMENTATION - UNDER DEVELOPMENT

DSEC = 24.*3600.
YDAYS = 365.

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
    """

    n0 = (year-1)//4
    n1 = (year-1)//100
    n2 = (year-1)//400

    return n0 - n1 + n2

def date_to_sec_mine(year=1, month=1, day=1, hour=0, minute=0, second=0.0):
    """
    Convert a date to seconds.
    Not yet implemented for b.c. years.
    """

    if year < 1:
      print('Error: Year must be positive (> 1)')
      return None

    if leap_check(year):
        MDAYS = [0.,31.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.]
    else:
        MDAYS = [0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.]

    ysec = (year - 1) * YDAYS * DSEC
    ysec += leap_num(year) * DSEC
    msec = MDAYS[int(month) - 1] * DSEC
    dsec = (day - 1) * DSEC

    sec = ysec + msec + dsec + hour*3600.+ minute*60. + second*1.

    return sec

def sec_to_date_mine(second):
    """
    IN PROGRESS!!!
    """
    # Number of days in 4 years (including a leap year)
    ysec = YDAYS * DSEC
    ysec4 = ysec*4 + DSEC

    nsec4 = (second // ysec4)
    if nsec4 > 0:
        sres = second - (nsec4 * ysec4)
    else:
        sres = second

    nsec = sres // ysec
    year = (nsec4 * 4) + nsec

    # Residual seconds in the last year
    if nsec > 0:
        sres = sres - (nsec * ysec)

    # Number of days left in the last year
    days = sres // DSEC

    if leap_check(year):
        MDAYS = [31.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,365.]
    else:
        MDAYS = [31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.,366.]

    for m,d in enumerate(MDAYS):
        if d > days:
            month = m
            break

    sres = sres - (MDAYS[int(month)] * DSEC)

    print(year, month, sres)

