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

import numpy as _np

def source_spectrum(m0, omega, omega_corner, order=2., rho=2800.,
                              beta=3200., gamma=2, n=1):
    """
    Compute the source spectrum shape using the Brune omega-square
    model

    :param float f:
        The input frequency

    :param class event:
        The event-related parameters

    :param float gamma:
         Falloff of the displacement amplitude spectra
         above the comer frequency

    :param 

    :return float ss:
        The Brune source spectrum shape
    """

    THETA = 0.55
    F = 2.
    CSI = _np.sqrt(0.5)

    C = (THETA*F*CSI)/(4.*_np.pi*rho*beta**3)

    plateau = m0*C*omega**order

    spectrum = plateau/(1.+(omega/omega_corner)**(n*gamma))**(1./n)

    return spectrum


def magnitude_to_moment(mw):
    """
    Hanks & Kanamori equation
    """
    return 10.**((3./2.)*(mw+10.7))


def moment_to_magnitude(m0):
    """
    Hanks & Kanamori equation
    """
    return ((2./3.)*_np.log10(m0))-10.7







