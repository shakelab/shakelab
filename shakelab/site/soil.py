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
The module contains several functions to derive and manipulate average
soil parameters, such as travel-time average velocity and site kappa.
"""

import numpy as _np
import scipy.optimize as _spo


def depth_weighted_average(thickness, soil_param, depth):
    """
    Compute the weighted average of a soil property at
    arbitrary depth.

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array soil_param:
        array of soil properties (e.g. slowness, density)

    :param float depth:
        averaging depth in meters

    :return float mean_param:
        the weighted mean of the given soil property
    """

    mean_param = 0
    total_depth = 0

    for tk, sp in zip(thickness[:-1], soil_param[:-1]):
        if (tk + total_depth) < depth:
            mean_param += tk*sp/depth
        else:
            mean_param += (depth - total_depth)*sp/depth
            break
        total_depth += tk

    # Check for the half-space
    if total_depth == _np.sum(thickness[:-1]):
        mean_param += (depth - total_depth)*soil_param[-1]/depth

    return mean_param


def traveltime_velocity(thickness, s_velocity, depth=30):
    """
    The function calculates the travel-time average (harmonic mean)
    velocity at arbitrary depth (e.g. the widespread Vs30).

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param float depth:
        averaging depth in meters; if depth is not specified,
        depth is fixed to 30m

    :return float mean_velocity:
        the average velocity in m/s
    """

    # Converting velocity to slowness
    slowness = 1./s_velocity

    # Harmonic averaging is done in slowness
    mean_slowness = depth_weighted_average(thickness, slowness, depth)

    # Back to velocity
    mean_velocity = 1./mean_slowness

    return mean_velocity


def compute_site_kappa(thickness, s_velocity, s_quality, depth=[]):
    """
    This function calucalte the site attenuation parameter Kappa(0)
    for a given soil profile at arbitrary depth.

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param numpy.array s_quality:
        array of layer's shear-wave quality factors (adimensional)

    :param float depth:
        averaging depth in meters; if depth is not specified,
        the last layer interface is used instead

    :return float kappa0:
        the site attenuation parameter kappa(0) in seconds
    """

    # If depth not given, using the whole profile
    if not depth:
        depth = _np.sum(thickness)

    # Kappa vector
    layer_kappa = depth/(s_velocity*s_quality)

    # Computing the average kappa of the profile
    kappa0 = depth_weighted_average(thickness, layer_kappa, depth)

    return kappa0


def quarter_wavelength_average(thickness, s_velocity, density, frequency):
    """
    This function solves the quarter-wavelength problem (Boore 2003)
    and return the frequency-dependent average velocity and density

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param numpy.array density:
        array of layer's densities in kg/m3

    :param numpy.array frequency:
        array of frequencies in Hz for the calculation

    :return numpy.array qwl_depth:
        array of averaging depths

    :return numpy.array qwl_velocity:
        array of quarter-wavelength average velocities

    :return numpy.array qwl_density:
        array of quarter-wavelength average dencities
    """

    # Initialisation
    freq_num = len(frequency)
    slowness = 1./s_velocity

    qwl_depth = _np.zeros(freq_num)
    qwl_velocity = _np.zeros(freq_num)
    qwl_density = _np.zeros(freq_num)

    for nf in range(freq_num):

        # Upper depth bound for the search
        ubnd = _np.max(1./(4.*frequency[nf]*slowness))

        # Input arguments for the search function
        args = (thickness, slowness, frequency[nf])

        # Compute the quarter-wavelength depth
        qwl_depth[nf] = _spo.fminbound(_qwl_fit_func, 0., ubnd, args)

        # Computing average soil property at the qwl-depth
        qwl_velocity[nf] = 1./depth_weighted_average(thickness,
                                                     slowness,
                                                     qwl_depth[nf])

        # Computing average soil property at the qwl-depth
        qwl_density[nf] = depth_weighted_average(thickness,
                                                 density,
                                                 qwl_depth[nf])

    return qwl_depth, qwl_velocity, qwl_density


def _qwl_fit_func(search_depth, thickness, slowness, frequency):
    """
    Internal function to recursively search for the qwl depth.
    """

    qwl_slowness = depth_weighted_average(thickness,
                                          slowness,
                                          search_depth)

    # Misfit is computed as a simple L1 norm
    misfit = _np.abs(search_depth - (1./(4.*frequency*qwl_slowness)))

    return misfit


def gt_soil_class(vs30, code='EC8'):
    """
    Compute geotechnical soil class from a given vs30 according
    to a specified building code.

    :param float vs30:
        The travel-time average over the first 30m

    :param string code:
        The reference building code for the classification;
        default is EC8

    reuturn string gt_class:
        Label of the geotechnical soil class
    """

    if code == 'EC8':
        if vs30 >= 800.:
            gt_class = 'A'
        elif vs30 >= 360. and vs30 < 800.:
            gt_class = 'B'
        elif vs30 >= 180. and vs30 < 360.:
            gt_class = 'C'
        elif vs30 < 180.:
            gt_class = 'D'
        else:
            gt_class = None
            print('Warning: no class assigned')

    # NEHRP (BSSC 1997)
    if code == 'NEHRP':
        if vs30 >= 1500.:
            gt_class = 'A'
        elif vs30 >= 760. and vs30 < 1500.:
            gt_class = 'B'
        elif vs30 >= 360. and vs30 < 760.:
            gt_class = 'C'
        elif vs30 >= 180. and vs30 < 360.:
            gt_class = 'D'
        elif vs30 < 180.:
            gt_class = 'E'
        else:
            gt_class = None
            print('Warning: no class assigned')

    return gt_class
