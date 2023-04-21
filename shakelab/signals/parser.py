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
import numpy as np

from shakelab.signals.base import Record


def read_record(ascii_file, dt=None, skipline=0, delimiter=','):
    """
    Generic parser for tabular (ascii) data files.
    """

    data = []

    # Open input ascii file
    with open(ascii_file, 'r') as f:

        # Ignore initial line(s) if necessary
        for i in range(0, skipline):
            next(f)

        # Import data (strip an collapse white spaces)
        for line in f:
            line = " ".join(line.split())
            data.append(line.split(delimiter))

        # Split data columns into recordings
        rnum = len(data[0])
        rec_list = [0]*rnum
        for i in range(0, rnum):
            rec_list[i] = Record()
            rec_list[i].data = np.array([d[i] for d in data], dtype=float)
            if dt is not None:
                rec_list[i].dt = dt

        return rec_list

    # Warn user if model file does not exist
    print('Error: file not found.')
