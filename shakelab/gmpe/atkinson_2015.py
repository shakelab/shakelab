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

class Atkinson2015(GMPE):
    """
    Reference pubblication:
        Atkinson, G. A. (2015) Ground-Motion Prediction Equation for
        Small-to-Moderate Events at Short Hypocentral Distances,
        with Application to Induced-Seismicity Hazards.
        Bulletin of the Seismological Society of America. 105(2)
    """

    REFERENCE_VELOCITY = 760.
    DISTANCE_METRIC = 'hypocentral'
    MAGNITUDE_TYPE = 'Mw'

    _COEFF_FILE = 'atkinson_2015.json'

    def ground_motion(self, imt, mag, dist):

        # Extract coefficients for the given intensity measure type
        C = self.get_coefficients(imt)

        # Compute mean
        heff = _np.max([1., 10.**(-1.72+0.43*mag)])
        mean = C['c0'] + C['c1']*mag + C['c2']*(mag**2.)
        mean += C['c3']*_np.log10(_np.sqrt(dist**2. + heff**2.))

        # Compute standard deviation
        stdv = C['phi']

        # Convert to natural log
        mean *= _np.log(10.)
        stdv *= _np.log(10.)

        if not imt == 'PGV':
            # Convert from cm/s2 to g
            mean -= _np.log(100.*g)
        else:
            # Convert from cm/s to m/s
            mean -= _np.log(100.)

        return (mean, stdv)

