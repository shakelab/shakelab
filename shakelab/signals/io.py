# ============================================================================
#
# Copyright (C) 2019, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ============================================================================
"""

"""

from shakelab.signals.iolibs import sac
from shakelab.signals.waveform import Recording


def import_recording (file, use_path=None, file_type='sac', byte_order='le',
                      **kwargs):
    """
    """

    # Initialise an empty trace
    rec = Recording()

    # Import recording from file
    if file_type is 'sac':

        if isinstance(file, list):
            files = file
        else:
            files = [file]

        for f in files:
            s = sac.Sac()
            s.read(f, byte_order=byte_order)

            rec.dt = s.head['DELTA']
            rec.channel.append(s.data[0])

    if file_type is 'ascii':
        print('ascii')

    return rec

def export_recording(file, use_path=None, file_type='sac', **kwargs):
    print('export')
