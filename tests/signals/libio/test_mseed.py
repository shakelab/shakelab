# ****************************************************************************
#
# Copyright (C) 2019-2022, ShakeLab Developers.
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
Testing code for the miniseed library
"""

import unittest

from shakelab.signals.libio.mseed import MiniSeed


class TestMiniSeed(unittest.TestCase):

    def test_read(self):
        """
        """
        fibonacci = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233]

        for enc in [1, 3, 4, 10, 11]:
            for rlen in [256, 512, 1024, 2048, 4096]:
                for bord in ['le']:

                    fname = 'fibonacci_enc{0}'.format(enc)
                    fname += '_len{0}'.format(rlen)
                    fname += '_{0}.mseed'.format(bord)

                    ms = MiniSeed('data/' + fname, byte_order=bord)
                    data = ms.stream['AA.NOA.00.Z'][0].data
                    self.assertEqual(data, fibonacci)


if __name__ == '__main__':
    unittest.main()