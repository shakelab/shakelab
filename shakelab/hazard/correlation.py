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


def baker_jayaram(t1, t2):
    """
    Compute the correlation of epsilons for the NGA ground motion models.
    Adapted from the original code of Jack Baker.

    Reference:
    Baker JW, Jayaram N. (2008). Correlation of spectral acceleration values
    from NGA ground motion models. Earthquake Spectra, 24(1), 299–317

    input:
        t1, t2 - The two periods of interest. The periods may be equal,
                 and there is no restriction on which one is larger.
    output:
        rho - The predicted correlation coefficient
    """

    tmin = np.min([t1, t2])
    tmax = np.max([t1, t2])

    if tmin < 0.01:
        print('Warning: minimum period not in the validity range')
    if tmax > 10:
        print('Warning: maximum period not in the validity range')

    C1 = 1. - np.cos(np.pi/2. - np.log(tmax/np.max([tmin, 0.109])) * 0.366)

    if tmax < 0.2:
        C2 = 1. - 0.105*(1. - 1./(1. + np.exp(100.*tmax - 5.)))
        C2 *= (tmax-tmin)/(tmax-0.0099)

    if tmax < 0.109:
        C3 = C2
    else:
        C3 = C1

    C4 = C1 + 0.5*(np.sqrt(C3) - C3) * (1 + np.cos(np.pi*tmin/0.109))

    if tmax <= 0.109:
        rho = C2
    elif tmin > 0.109:
        rho = C1
    elif tmax < 0.2:
        rho = np.min([C2, C4])
    else:
        rho = C4

    return rho