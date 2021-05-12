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

import numpy as _np
from scipy.constants import g

from shakelab.modeling.gmpe.base import GMPE

class BragatoSlejko2005(GMPE):
    """
    Reference pubblication:
        Bragato P.L. and Slejko D. (2005).
        Empirical Ground-Motion Attenuation Relations for the
        Eastern Alps in the Magnitude Range 2.5–6.3.
        Bulletin of the Seismological Society of America.
        Vol. 95, No. 1, pp. 252–276,

    Distance: epicentral
    """

    REFERENCE_VELOCITY = 800.
    DISTANCE_METRIC = 'epicentral'
    MAGNITUDE_TYPE = 'ML'

    _COEFF_FILE = 'bragato_slejko_2005.json'
    _COEFF_SET = 'epicentral'

    def ground_motion(self, imt, mag, dist):

        # Extract coefficients for the given intensity measure type
        C = self.get_coefficients(imt)

        # Distance term
        r = _np.sqrt(dist**2 + C['h']**2)

        # Compute mean
        mean = C['a'] + (C['b'] + C['c']*mag)*mag
        mean += (C['d'] + C['e']*mag**3)*_np.log10(r)

        # Compute standard deviation
        stdv = C['s']

        # Convert to natural log
        mean *= _np.log(10.)
        stdv *= _np.log(10.)

        if imt == 'PGV':
            # Convert from cm/s to m/s
            mean -= _np.log(100.)

        return (mean, stdv)


class BragatoSlejko2005JB(BragatoSlejko2005):
    """
    """

    DISTANCE_METRIC = 'joyner-boore'
    _COEFF_SET = 'joyner-boore'
