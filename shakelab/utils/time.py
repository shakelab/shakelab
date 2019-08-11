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
"""

DSEC = 24.*3600.
YDAYS = 365.

class Date(object):
    """
    """

    year = 1
    month = 1
    day = 1
    hour = 0
    minute = 0
    second = 0.0

    def __init__(self, date_list=None):
        if date_list is not None:
            self.set(date_list)

    def set(self, date_list):
        self.year = date_list[0]
        self.month = date_list[1]
        self.day = date_list[2]
        self.hour = date_list[3]
        self.minute = date_list[4]
        self.second = date_list[5]

    def get(self):
        return [self.year, self.month, self.day,
                self.hour, self.minute, self.second]

    def to_sec(self):
        return date_to_sec(self.year, self.month, self.day,
                           self.hour, self.minute, self.second)

    def shift_time(self):
        print('to_implement')

    def __add__(self, date):
        if isinstance(date, (int, float)):
            return self.to_sec() + date
        if isinstance(date, Date):
            return self.to_sec() + date.to_sec()

    def __sub__(self, date):
        if isinstance(date, (int, float)):
            return self.to_sec() - date
        if isinstance(date, Date):
            return self.to_sec() - date.to_sec()

# TO IMPLEMENT
def sec_to_date(second):
    """
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

def date_to_sec(year=1, month=1, day=1, hour=0, minute=0, second=0.0):
    """
    Convert a date to seconds.
    Not yet implemented for b.c. years.
    """

    if year < 1:
      print('Error: Year must be > 1')
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
