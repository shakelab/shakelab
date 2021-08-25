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
Collection of methods to compute seismic amplification from a
one-dimensional soil column.
"""

import numpy as _np


def frequency_axis(fmin, fmax, fnum, log=True):
    """
    Compute a linear or logarithmic frequency axis

    :param float fmin:
        Minimum frequency

    :param float fmax:
        Maximum frequency

    :param int fnum:
        Number of frequencies

    :param boolean log:
        Switch between linear or logarithmic spacing
        (default is logarithmic)

    :return numpy.array freq:
        The frequency axis
    """

    if log:
        freq = _np.logspace(_np.log10(fmin), _np.log10(fmax), fnum)
    else:
        freq = _np.linspace(fmin, fmax, fnum)

    return freq


def impedance_amplification(top_vs, top_dn, ref_vs=[], ref_dn=[], inc_ang=0.):
    """
    This function calculates the amplification due to a single seismic
    impedance contrast as formalized in Joyner et al. (1981) and
    Boore (2013).

    Site (top) parameters can be either scalars or arrays, as for
    the case of quarter-wavelenght approximation. If no reference
    is provided, the last value of the site array is used.
    Default angle of incidence is vertical.

    :param float or numpy.array top_vs:
        topmost shear-wave velocity in m/s

    :param float or numpy.array top_dn:
        topmost density in kg/m3

    :param float or numpy.array ref_vs:
        lowermost (reference) shear-wave velocity in m/s

    :param float or numpy.array ref_dn:
        lowermost (reference) density in kg/m3

    :param float inc_ang:
        angle of incidence in degrees, relative to the vertical
        (default is vertical incidence)

    :return float or array imp_amp:
        amplification due to the seismic impedance contrast
    """

    # Setting the reference, if not provided
    if not _np.sum(ref_vs):
        ref_vs = top_vs[-1] if _np.size(top_vs) > 1. else top_vs

    if not _np.sum(ref_dn):
        ref_dn = top_dn[-1] if _np.size(top_dn) > 1. else top_dn

    # Computing square-root impedance amplification
    imp_amp = _np.sqrt((ref_dn*ref_vs)/(top_dn*top_vs))

    if inc_ang > 0.:
        # Converting incident angle from degrees to radiants
        inc_ang *= _np.pi/180.

        # Effective angle of incidence computed using Snell's law
        eff_ang = _np.arcsin((top_vs/ref_vs)*_np.sin(inc_ang))

        # Correcting for non-vertical incidence
        imp_amp *= _np.sqrt(_np.cos(inc_ang)/_np.cos(eff_ang))

    return imp_amp


def attenuation_decay(freq, kappa):
    """
    Compute the frequency-dependent attenuation decay function
    for a given site Kappa (0).

    :param float or numpy.array freq:
            array of frequencies in Hz for the calculation

    :param float kappa:
        site-specific attenutation operator (kappa0)

    :return numpy.array att_func:
        attenuation decay function (array)
    """

    # Computing attenuation function
    att_fun = _np.exp(-_np.pi*kappa*freq)

    return att_fun


def sh_transfer_function(freq, hl, vs, dn, qs=None, inc_ang=0., depth=0.):
    """
    Compute the SH-wave transfer function using Knopoff formalism
    (implicit layer matrix scheme). Calculation can be done for an
    arbitrary angle of incidence (0-90), with or without anelastic
    attenuation (qs is optional).

    It return the displacements computed at arbitrary depth.
    If depth = -1, calculation is done at each layer interface
    of the profile.

    NOTE: the implicit scheme is simple to understand and to implement,
    but is also computationally intensive and memory consuming.
    For the future, an explicit (recursive) scheme should be implemented.

    :param float or numpy.array freq:
        array of frequencies in Hz for the calculation

    :param numpy.array hl:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array vs:
        array of layer's shear-wave velocities in m/s

    :param numpy.array dn:
        array of layer's densities in kg/m3

    :param numpy.array qs:
        array of layer's shear-wave quality factors (adimensional)

    :param float inc_ang:
        angle of incidence in degrees, relative to the vertical
        (default is vertical incidence)

    :param float or numpy.array depth:
        dephts in meters at which displacements are calculated
        (default is the free surface)

    :return numpy.array dis_mat:
        matrix of displacements computed at each depth (complex)
    """

    # Precision of the complex type
    CTP = 'complex128'

    # Check for single frequency value
    if isinstance(freq, (int, float)):
        freq = _np.array([freq])

    # Model size
    lnum = len(hl)
    fnum = len(freq)

    # Variable recasting to numpy complex
    hl = _np.array(hl, dtype=CTP)
    vs = _np.array(vs, dtype=CTP)
    dn = _np.array(dn, dtype=CTP)

    # Attenuation using complex velocities
    if qs is not None:
        qs = _np.array(qs, dtype=CTP)
        vs *= ((2.*qs*1j)/(2.*qs*1j-1.))

    # Conversion to angular frequency
    angf = 2.*_np.pi*freq

    # Layer boundary depth (including free surface)
    bounds = interface_depth(hl, dtype=CTP)

    # Check for depth to calculate displacements
    if isinstance(depth, (int, float)):
        if depth < 0.:
            depth = _np.array(bounds)
        else:
            depth = _np.array([depth])
    znum = len(depth)

    # Computing angle of propagation within layers

    iD = _np.zeros(lnum, dtype=CTP)
    iM = _np.zeros((lnum, lnum), dtype=CTP)

    iD[0] = _np.sin(inc_ang)
    iM[0, -1] = 1.

    for nl in range(lnum-1):
        iM[nl+1, nl] = 1./vs[nl]
        iM[nl+1, nl+1] = -1./vs[nl+1]

    iA = _np.linalg.solve(iM, iD)
    iS = _np.arcsin(iA)

    # Elastic parameters

    # Lame Parameters : shear modulus
    mu = dn*(vs**2.)

    # Horizontal slowness
    ns = _np.cos(iS)/vs

    # Data vector initialisation

    # Layer's amplitude vector (incognita term)
    amp_vec = _np.zeros(lnum*2, dtype=CTP)

    # Layer matrix
    lay_mat = _np.zeros((lnum*2, lnum*2), dtype=CTP)

    # Input motion vector (known term)
    inp_vec = _np.zeros(lnum*2, dtype=CTP)
    inp_vec[-1] = 1.

    # Output layer's displacement matrix
    dis_mat = _np.zeros((znum, fnum), dtype=CTP)

    # Loop over frequencies

    for nf in range(fnum):

        # Reinitialise the layer matrix
        lay_mat *= 0.

        # Free surface constraints
        lay_mat[0, 0] = 1.
        lay_mat[0, 1] = -1.

        # Interface constraints
        for nl in range(lnum-1):
            row = (nl*2)+1
            col = nl*2

            exp_dsa = _np.exp(1j*angf[nf]*ns[nl]*hl[nl])
            exp_usa = _np.exp(-1j*angf[nf]*ns[nl]*hl[nl])

            # Displacement continuity conditions
            lay_mat[row, col+0] = exp_dsa
            lay_mat[row, col+1] = exp_usa
            lay_mat[row, col+2] = -1.
            lay_mat[row, col+3] = -1.

            # Stress continuity conditions
            lay_mat[row+1, col+0] = mu[nl]*ns[nl]*exp_dsa
            lay_mat[row+1, col+1] = -mu[nl]*ns[nl]*exp_usa
            lay_mat[row+1, col+2] = -mu[nl+1]*ns[nl+1]
            lay_mat[row+1, col+3] = mu[nl+1]*ns[nl+1]

        # Input motion constraints
        lay_mat[-1, -1] = 1.

        # Solving linear system of wave's amplitudes
        try:
            amp_vec = _np.linalg.solve(lay_mat, inp_vec)
        except:
            amp_vec[:] = _np.nan

        # Solving displacements at depth
        for nz in range(znum):

            # Check in which layer falls the calculation depth
            if depth[nz] <= hl[0]:
                nl = 0
                dh = depth[nz]
            elif depth[nz] > bounds[-1]:
                nl = lnum-1
                dh = depth[nz] - bounds[-1]
            else:
                # There might be a more pythonic way to do that...
                zlim = lambda x: x >= depth[nz]
                nl = list(map(zlim, bounds)).index(True) - 1
                dh = depth[nz] - bounds[nl]

            # Displacement of the up-going and down-going waves
            exp_dsa = _np.exp(1j*angf[nf]*ns[nl]*dh)
            exp_usa = _np.exp(-1j*angf[nf]*ns[nl]*dh)

            dis_dsa = amp_vec[nl*2]*exp_dsa
            dis_usa = amp_vec[nl*2+1]*exp_usa

            dis_mat[nz, nf] = dis_dsa + dis_usa

    return dis_mat


def interface_depth(hl, dtype='float64'):
    """
    Utility to calcualte the depth of the layer's interface
    (including the free surface) from a 1d thickness profile.

    :param numpy.array hl:
        array of layer's thicknesses in meters (half-space is 0.)

    :param string dtype:
        data type for variable casting (optional)

    :return numpy.array depth:
        array of interface depths in meters
    """

    depth = [sum(hl[:i]) for i in range(len(hl))]
    depth = _np.array(depth, dtype=dtype)

    return depth


def resonance_frequency(freq, spec):
    """
    Identify resonance frequencies on an amplification spectrum.

    :param float or numpy.array freq:
        array of frequencies in Hz for the calculation

    :param float or numpy.array spec:
        the amplification spectrum

    :return list resf:
        list of tuples containg each a resonance frequency
        and the corresponding amplitude
    """

    resf = []
    spec = _np.abs(spec)

    for nf in range(0, len(freq[:-2])):
        # Three-points search for local maxima
        a0 = spec[nf]
        a1 = spec[nf+1]
        a2 = spec[nf+2]

        if (a1-a0) > 0 and (a2-a1) < 0:
            resf.append((freq[nf+1], a1))

    return resf
