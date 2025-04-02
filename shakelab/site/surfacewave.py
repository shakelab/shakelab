# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
Module for surface wave computation
"""
import warnings
import numpy as np

from typing import Callable, Iterable
from scipy.optimize import root_scalar


def love_wave_dispersion(hl, vs, dn, frequencies, num_modes):
    """
    Computes Love wave dispersion curves over multiple frequencies and modes.

    Parameters:
    - hl: Array of layer thicknesses (m).
    - vs: Array of S-wave velocities (m/s).
    - dn: Array of densities (kg/m^3).
    - frequencies: Array of frequencies (Hz) to compute dispersion curves for.
    - num_modes: Number of modes to compute for each frequency.
    
    Returns:
    - phase_velocities_matrix: Matrix of phase velocities for each mode
      and frequency.
    """
    phase_velocities_matrix = np.zeros((num_modes, len(frequencies)))

    for i, freq in enumerate(frequencies):
        modes = love_wave_roots(freq, len(hl), hl, vs, dn)
        # Store the first 'num_modes' modes in the matrix
        phase_velocities_matrix[:len(modes), i] = modes[:num_modes]
        #print(f"Computed phase velocities for {freq:.2f} Hz: {modes}")

    return phase_velocities_matrix

def love_wave_roots(freq, lnum, hl, vs, dn):
    """
    Finds multiple phase velocities (modes) for the given frequency
    using multi_root.

    Parameters:
    - freq: Frequency (Hz).
    - lnum: Number of layers.
    - hl: Array of layer thicknesses (m).
    - vs: Array of S-wave velocities (m/s).
    - dn: Array of densities (kg/m^3).

    Returns:
    - A list of phase velocities corresponding to different modes for
      the given frequency.
    """

    # Angular frequency
    omega = 2 * np.pi * freq

    def search_function(vl):
        """
        """
        if isinstance(vl, np.ndarray):
            results = np.zeros(vl.shape, dtype=complex)
            for idx, v in enumerate(vl):
                results[idx] = love_wave_solution(lnum, hl, vs, dn, omega, v)
            return results
        else:
            return love_wave_solution(lnum, hl, vs, dn, omega, vl)

    # Set bounds for phase velocity (between min and max vs)
    min_vs = np.min(vs)
    max_vs = np.max(vs)

    # Use multi_root to find multiple roots (modes)
    modes = multi_root(search_function, [min_vs, max_vs], n=500)

    return sorted(modes)


def love_wave_solution(lnum, hl, vs, dn, omega, vl):
    """
    Function to compute the characteristic equation for a scalar vl.
    Python translation of the original Geopsy function from Marc Wathelet.
    Note: this version does not need complex arithmetics.
    """
    k = omega / vl  # Horizontal wavenumber (scalar)
    i0 = lnum - 1  # Starting from the bottom layer
    
    slowS = 1.0 / vs[i0]
    ks = omega * slowS
    # Ensure nu is real by taking abs before sqrt
    nu = np.sqrt(np.abs(k**2 - ks**2))
    l21 = dn[i0] * nu
    l22 = slowS**2

    # Iterating over layers from bottom to top
    for i in range(i0 - 1, -1, -1):
        slowS = 1.0 / vs[i]
        ks = omega * slowS
        # Same here, ensure valid sqrt
        nu = np.sqrt(np.abs(k**2 - ks**2))
        q = hl[i] * nu
        mui = dn[i] * (vs[i]**2)
        numui = nu * mui
        
        # Handle scalar comparisons correctly
        if k < ks:
            g11 = np.cos(q)
            sinq = np.sin(q)
            g12 = sinq / numui
            g21 = -numui * sinq
        elif np.isclose(k, ks):
            g11 = 1.0
            g12 = hl[i] / mui
            g21 = 0.0
        else:
            fac = np.exp(-2.0 * q) if q < 21.2 else 0.0
            g11 = (1.0 + fac) * 0.5
            sinq = (1.0 - fac) * 0.5
            g12 = sinq / numui
            g21 = numui * sinq
        
        l21n = l21 * g11 + l22 * g21
        l22n = l21 * g12 + l22 * g11
        maxL = max(np.abs(l21n), np.abs(l22n))
        
        if maxL > 1.0e5:
            scale_factor = 1.0e5 / maxL
            l21 = l21n * scale_factor
            l22 = l22n * scale_factor
        else:
            l21 = l21n
            l22 = l22n

    return l21

# ----------------------------------------------------------------------------
# UTILITY FUNCTIONS

def multi_root(f: Callable, bracket: Iterable[float],
               args: Iterable = (), n: int = 100) -> np.ndarray:
    """
    Find all roots of f in `bracket`, given that resolution `n` covers
    the sign change. Fine-grained root finding is performed with
    `scipy.optimize.root_scalar`.

    credits to https://gist.github.com/maxiimilian

    Parameters
    ----------
    f: Callable
        Function to be evaluated
    bracket: Sequence of two floats
        Specifies interval within which roots are searched.
    args: Iterable, optional
        Iterable passed to `f` for evaluation
    n: int
        Number of points sampled equidistantly from bracket to evaluate `f`.
        Resolution has to be high enough to cover sign changes of all roots
        but not finer than that.
        Actual roots are found using `scipy.optimize.root_scalar`.

    Returns
    -------
    roots: np.ndarray
        Array containing all unique roots that were found in `bracket`.
    """
    # Evaluate function in given bracket
    x = np.linspace(*bracket, n)
    y = f(x, *args)

    # Find where adjacent signs are not equal
    sign_changes = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]

    # Find roots around sign changes
    root_finders = (
        root_scalar(
            f=f,
            args=args,
            bracket=(x[s], x[s+1])
        )
        for s in sign_changes
    )

    roots = np.array([
        r.root if r.converged else np.nan for r in root_finders
    ])

    if np.any(np.isnan(roots)):
        warnings.warn("Not all root finders converged for estimated brackets!"
                      "Maybe increase resolution `n`.")
        roots = roots[~np.isnan(roots)]

    roots_unique = np.unique(roots)
    if len(roots_unique) != len(roots):
        warnings.warn("One root was found multiple times. "
                      "Try to increase or decrease resolution `n`"
                      "to see if this warning disappears.")

    return roots_unique
