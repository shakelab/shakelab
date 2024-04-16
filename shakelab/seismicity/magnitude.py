# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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
Module to handle magnitude calculation.
"""

import numpy as np

def ml_IASPEI_2011(wa, dist):
    """
    DOI: 10.2312/GFZ.NMSOP-2_DS_3.1
    """
    return np.log10(wa) + 1.11*np.log10(dist) + 0.00189*dist - 2.09

def ml_Hutton_Boore87(wa, dist):
    """
    Attenuation law based on Hutton and Boore, BSSA, 1987:
    -log (Ao) = 1.110 log(r/100) + 0.00189(r - 100) + 3.0
    """

    c0 = -18.0471
    c1 = np.log10(wa)
    c2 = 1.105 * np.log10(dist)
    c3 = 147.111 * np.log10(dist * 4.015e-5 + 1.33885)

    return (c0 + c1 + c2 + c3)
