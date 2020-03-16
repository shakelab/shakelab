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
Module for basic waveform analysis
"""

from shakelab.signals import mseed, sac
from scipy import signal


def import_record(file, format='sac', path=None, byte_order='le',
                     **kwargs):
    """
    """

    # Initialise an empty trace
    rec_list = []

    # Import recordings from file
    if format is 'mseed':
        ms = miniseed.MiniSeed(file, byte_order=byte_order)
        for mr in ms.record:
            rec = Record()
            rec.dt = mr.sampling_rate()
            rec.time = mr.time_date()
            rec.data = mr.data
            rec_list.append(rec)

    elif format is 'sac':
        sac = sac.Sac(file, byte_order=byte_order)
        rec = Record()
        rec.dt = sac.sampling_rate()
        rec.time = sac.time_date()
        rec.data = sac.data[0]

    elif format is 'ascii':
        raise NotImplementedError('format not yet implemented')

    elif format is 'seisan':
        raise NotImplementedError('format not yet implemented')

    elif format is 'seg2':
        raise NotImplementedError('format not yet implemented')

    elif format is 'dat':
        raise NotImplementedError('format not yet implemented')

    elif format is 'gse':
        raise NotImplementedError('format not yet implemented')

    else:
        raise NotImplementedError('format not recognized')

    return rec_list



class Record(object):
    """
    Single recording.
    """

    def __init__(self):
        self.dt = None
        self.data = []
        self.time = Date()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, sliced):
        return self.data[sliced]

    def filter(self, fl=None, fh=None, order=2):
        """
        """

        # Corner frequencies
        corners = []

        if fl:
            corners.append(2. * fl * self.dt)
            filter_type = 'high'

        if fh:
            corners.append(2. * fh * self.dt)
            filter_type = 'low'

        if fl and hl:
            filter_type = 'band'

        if len(corners) > 0:
            # Butterworth filter
            b, a = signal.butter(order, corners, btype=filter_type)

            # Filtering seismic record
            zi = signal.lfilter_zi(b, a);
            self.data,_ = signal.lfilter(b, a, self.data, zi=zi*self.data[0])

    def cut(self, start, stop):
        """
        """
        pass

