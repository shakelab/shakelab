# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Frequency-wavenumber estimators.

Functions on arrays that evaluate the beam power of a cross-spectral
matrix ``C`` (``n x n``) at the wavevectors of a search grid, given by
their steering vectors ``a(k)``:

- Bartlett (conventional) beamforming: ``a^H C a / n**2``;
- Capon (minimum variance, high resolution): ``1 / (a^H R^-1 a)``,
  with the diagonally loaded matrix ``R``;
- MUSIC: ``1 / (a^H P_noise a)``, a pseudo-power that is large where
  ``a`` is orthogonal to the noise subspace.

The same functions serve the processing of recorded data
(:mod:`.beamforming`) and the controlled simulations of the performance
analysis, so that both use the same definitions.

Conventions
-----------
Horizontal coordinates x = East, y = North; wavevector ``(kx, ky)`` in
radians per distance unit; steering vectors
``a_i(k) = exp(-1j * (kx * x_i + ky * y_i))`` (not normalized). With a
cross-spectral matrix of power spectral densities (see
:mod:`shakelab.signals.crossspectral`), Bartlett and Capon give the
density of a plane wave at its wavevector, in the units of the data
squared per Hz.

References
----------
Capon, J. (1969). High-resolution frequency-wavenumber spectrum
analysis. Proceedings of the IEEE, 57(8), 1408-1418.

Schmidt, R. O. (1986). Multiple emitter location and signal parameter
estimation. IEEE Transactions on Antennas and Propagation, 34(3),
276-280.

Lacoss, R. T., Kelly, E. J., and Toksöz, M. N. (1969). Estimation of
seismic noise structure using arrays. Geophysics, 34(1), 21-38.
"""

from __future__ import annotations

from enum import Enum
from numbers import Integral
from typing import Any

import numpy as np

from shakelab.signals.crossspectral import cross_spectral_matrix

__all__ = [
    "BeamformingMethod",
    "bartlett_power",
    "bartlett_power_statistics",
    "beam_power",
    "capon_power",
    "capon_power_statistics",
    "check_realizations",
    "music_power",
    "regularized_matrix",
    "steering_vectors",
]


class BeamformingMethod(Enum):
    """Beam-power estimators."""

    BARTLETT = "bartlett"
    CAPON = "capon"
    MUSIC = "music"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _square_matrix(matrix: Any) -> np.ndarray:
    """Return a complex square matrix."""
    matrix = np.asarray(matrix, dtype=complex)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("the cross-spectral matrix must be square")

    return matrix


def _steering(steering: Any, size: int) -> np.ndarray:
    """Return steering vectors ``(..., size)``."""
    steering = np.asarray(steering, dtype=complex)

    if steering.ndim < 1 or steering.shape[-1] != size:
        raise ValueError(
            f"steering vectors must have {size} elements (last axis)"
        )

    return steering


def _method(method: Any) -> BeamformingMethod:
    """Return a method from a member or its name."""
    if isinstance(method, BeamformingMethod):
        return method

    try:
        return BeamformingMethod(str(method).lower())
    except ValueError:
        names = ", ".join(repr(item.value) for item in BeamformingMethod)
        raise ValueError(f"method must be one of {names}") from None


def _reciprocal(denominator: np.ndarray) -> np.ndarray:
    """Return 1 / x where x > 0, NaN elsewhere."""
    power = np.full(denominator.shape, np.nan)
    valid = denominator > 0.0
    power[valid] = 1.0 / denominator[valid]
    return power


def _check_effective(effective: float, count: int) -> float:
    """Return a valid effective number of realizations."""
    effective = float(effective)

    if not np.isfinite(effective) or not 0.0 < effective <= count:
        raise ValueError(
            "effective_realizations must be in (0, number of realizations]"
        )

    return effective


def _check_coefficients(coefficients: Any, steering: Any) -> tuple:
    """Return coefficients ``(L, n)`` (L >= 2) and steering vectors."""
    coefficients = np.asarray(coefficients, dtype=complex)

    if coefficients.ndim != 2 or coefficients.shape[0] < 2:
        raise ValueError(
            "coefficients must be (n_realizations >= 2, n_channels)"
        )

    return coefficients, _steering(steering, coefficients.shape[1])


# ----------------------------------------------------------------------------
# Matrices and steering vectors
# ----------------------------------------------------------------------------


def regularized_matrix(matrix: Any, loading: float) -> np.ndarray:
    """
    Return the Hermitian part of a matrix with diagonal loading.

    ``R = (C + C^H) / 2 + loading * trace(C) / n * I``: the loading is
    relative to the mean power on the diagonal.

    Parameters
    ----------
    matrix : array-like
        Cross-spectral matrix ``(n, n)``.
    loading : float
        Relative diagonal loading, non-negative.

    Returns
    -------
    numpy.ndarray
        Regularized matrix.

    Raises
    ------
    ValueError
        If the matrix is not square, the loading is negative or not
        finite, or the diagonal power is not positive (with loading).
    """
    matrix = _square_matrix(matrix)
    matrix = 0.5 * (matrix + matrix.conj().T)
    loading = float(loading)

    if not np.isfinite(loading) or loading < 0.0:
        raise ValueError("loading must be finite and non-negative")

    if loading == 0.0:
        return matrix

    size = matrix.shape[0]
    mean_power = float(np.real(np.trace(matrix))) / size

    if not np.isfinite(mean_power) or mean_power <= 0.0:
        raise ValueError("the matrix has no positive diagonal power")

    return matrix + loading * mean_power * np.eye(size)


def check_realizations(
    count: int,
    size: int,
    method: Any,
    strict: bool = True,
) -> None:
    """
    Check that a matrix average has enough realizations.

    Capon and MUSIC invert or decompose the matrix: with ``strict`` they
    require more realizations than its dimension (otherwise the
    matrix is singular). Bartlett has no requirement.

    Parameters
    ----------
    count : int
        Number of realizations (windows) averaged.
    size : int
        Dimension of the matrix.
    method : BeamformingMethod or str
        Estimator.
    strict : bool, optional
        Apply the requirement. Default is True.

    Raises
    ------
    ValueError
        If the requirement is not met.
    """
    method = _method(method)

    if not isinstance(count, Integral) or not isinstance(size, Integral):
        raise TypeError("count and size must be integers")

    if strict and method is not BeamformingMethod.BARTLETT and count <= size:
        raise ValueError(
            f"{method.value} needs more realizations than the matrix "
            f"dimension: {count} for a {size} x {size} matrix (average "
            "more windows per block, or strict=False with diagonal "
            "loading)"
        )


def steering_vectors(coordinates: Any, kx: Any, ky: Any) -> np.ndarray:
    """
    Return the steering vectors of plane waves.

    Parameters
    ----------
    coordinates : array-like
        Station coordinates ``(n, 2)`` or ``(n, 3)`` (east, north[,
        up]); only east and north are used.
    kx, ky : array-like
        East and north angular wavenumbers, same shape.

    Returns
    -------
    numpy.ndarray
        ``exp(-1j * (kx * x + ky * y))`` with shape ``kx.shape + (n,)``.

    Raises
    ------
    ValueError
        If the shapes are inconsistent or the values are not finite.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    kx = np.asarray(kx, dtype=float)
    ky = np.asarray(ky, dtype=float)

    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("coordinates must be (n_stations, 2 or 3)")

    if kx.shape != ky.shape:
        raise ValueError("kx and ky must have the same shape")

    if not (
        np.all(np.isfinite(coordinates[:, :2]))
        and np.all(np.isfinite(kx)) and np.all(np.isfinite(ky))
    ):
        raise ValueError("coordinates and wavenumbers must be finite")

    phase = (
        kx[..., None] * coordinates[:, 0] + ky[..., None] * coordinates[:, 1]
    )

    return np.exp(-1j * phase)


# ----------------------------------------------------------------------------
# Estimators
# ----------------------------------------------------------------------------


def bartlett_power(
    matrix: Any,
    steering: Any,
    attenuation: Any = None,
) -> np.ndarray:
    """
    Return the Bartlett beam power.

    ``P = a^H C a / n**2``. With ``attenuation``, the expected power
    over independent Gaussian errors of the station coordinates: the
    steering vectors are damped by the expected coherence ``g`` and the
    diagonal is kept, ``P = (a g)^H C (a g) / n**2 + sum_i (1 - g_i**2)
    C_ii / n**2``.

    Parameters
    ----------
    matrix : array-like
        Cross-spectral matrix ``(n, n)``.
    steering : array-like
        Steering vectors ``(..., n)``.
    attenuation : array-like, optional
        Expected coherence of each station, same shape as ``steering``.

    Returns
    -------
    numpy.ndarray
        Power with shape ``steering.shape[:-1]``, non-negative.

    Raises
    ------
    ValueError
        If the shapes are inconsistent.
    """
    matrix = _square_matrix(matrix)
    steering = _steering(steering, matrix.shape[0])

    if attenuation is None:
        projected = np.einsum(
            "...i,ij,...j->...", steering.conj(), matrix, steering,
            optimize=True,
        )
    else:
        attenuation = np.asarray(attenuation, dtype=float)

        if attenuation.shape != steering.shape:
            raise ValueError("attenuation must have the shape of steering")

        damped = steering * attenuation
        projected = np.einsum(
            "...i,ij,...j->...", damped.conj(), matrix, damped,
            optimize=True,
        ) + (1.0 - attenuation ** 2) @ np.real(np.diag(matrix))

    return np.maximum(np.real(projected) / matrix.shape[0] ** 2, 0.0)


def bartlett_power_statistics(
    coefficients: Any,
    steering: Any,
    effective_realizations: float,
    attenuation: Any = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the Bartlett power and its empirical standard error.

    The power of each realization is computed separately; the standard
    error is their standard deviation divided by the square root of the
    effective number of independent realizations.

    Parameters
    ----------
    coefficients : array-like
        Fourier coefficients ``(n_realizations, n)``, at least two.
    steering : array-like
        Steering vectors ``(..., n)``.
    effective_realizations : float
        Effective number of independent realizations.
    attenuation : array-like, optional
        Expected coherence (see :func:`bartlett_power`).

    Returns
    -------
    power, standard_error : numpy.ndarray
        Arrays with shape ``steering.shape[:-1]``.

    Raises
    ------
    ValueError
        If the shapes are inconsistent or fewer than two realizations
        are given.
    """
    coefficients, steering = _check_coefficients(coefficients, steering)
    count, size = coefficients.shape
    effective = _check_effective(effective_realizations, count)
    shape = steering.shape[:-1]
    flat = steering.reshape(-1, size)

    if attenuation is not None:
        attenuation = np.asarray(attenuation, dtype=float)

        if attenuation.shape != steering.shape:
            raise ValueError("attenuation must have the shape of steering")

        attenuation = attenuation.reshape(-1, size)
        flat = flat * attenuation

    power = np.abs(coefficients @ flat.conj().T) ** 2

    if attenuation is not None:
        power += (np.abs(coefficients) ** 2) @ (1.0 - attenuation ** 2).T

    power /= size ** 2
    error = np.sqrt(np.var(power, axis=0, ddof=1) / effective)

    return power.mean(axis=0).reshape(shape), error.reshape(shape)


def capon_power(matrix: Any, steering: Any, loading: float) -> np.ndarray:
    """
    Return the Capon beam power.

    ``P = 1 / (a^H R^-1 a)`` with the loaded matrix ``R`` (see
    :func:`regularized_matrix`).

    Parameters
    ----------
    matrix : array-like
        Cross-spectral matrix ``(n, n)``.
    steering : array-like
        Steering vectors ``(..., n)``.
    loading : float
        Relative diagonal loading.

    Returns
    -------
    numpy.ndarray
        Power with shape ``steering.shape[:-1]`` (NaN where the
        quadratic form is not positive).
    """
    matrix = regularized_matrix(matrix, loading)
    steering = _steering(steering, matrix.shape[0])
    flat = steering.reshape(-1, matrix.shape[0])
    solved = np.linalg.solve(matrix, flat.T)
    denominator = np.real(np.einsum("gi,ig->g", flat.conj(), solved))

    return _reciprocal(denominator).reshape(steering.shape[:-1])


def capon_power_statistics(
    coefficients: Any,
    steering: Any,
    effective_realizations: float,
    loading: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the Capon power and its first-order standard error.

    The quadratic form ``D = a^H R^-1 a`` is linear in the realizations
    to first order (delta method): its standard error is estimated from
    the spread of the contributions of the single realizations,
    including the dependence of the loading on the trace, and
    propagated to ``P = 1 / D`` as ``P**2 * se(D)``.

    Parameters
    ----------
    coefficients : array-like
        Fourier coefficients ``(n_realizations, n)``, at least two.
    steering : array-like
        Steering vectors ``(..., n)``.
    effective_realizations : float
        Effective number of independent realizations.
    loading : float
        Relative diagonal loading.

    Returns
    -------
    power, standard_error : numpy.ndarray
        Arrays with shape ``steering.shape[:-1]``.
    """
    coefficients, steering = _check_coefficients(coefficients, steering)
    count, size = coefficients.shape
    effective = _check_effective(effective_realizations, count)
    matrix = regularized_matrix(cross_spectral_matrix(coefficients), loading)
    flat = steering.reshape(-1, size)
    solved = np.linalg.solve(matrix, flat.T)
    denominator = np.real(np.einsum("gi,ig->g", flat.conj(), solved))
    power = _reciprocal(denominator)

    # Contribution of each realization to a^H R^-1 C R^-1 a, i.e. to
    # the first-order change of the denominator.
    single = np.abs(coefficients @ solved.conj()) ** 2

    if loading != 0.0:
        trace = np.sum(np.abs(coefficients) ** 2, axis=1)
        norm = np.sum(np.abs(solved) ** 2, axis=0)
        single += float(loading) * trace[:, None] * norm[None, :] / size

    error_denominator = np.sqrt(np.var(single, axis=0, ddof=1) / effective)
    shape = steering.shape[:-1]

    return power.reshape(shape), (power ** 2 * error_denominator).reshape(
        shape
    )


def music_power(matrix: Any, steering: Any, sources: int) -> np.ndarray:
    """
    Return the MUSIC pseudo-power.

    ``P = 1 / (a^H E E^H a)``, with ``E`` the eigenvectors of the
    ``n - sources`` smallest eigenvalues (noise subspace).

    Parameters
    ----------
    matrix : array-like
        Cross-spectral matrix ``(n, n)``.
    steering : array-like
        Steering vectors ``(..., n)``.
    sources : int
        Number of sources (signal subspace dimension), smaller than n.

    Returns
    -------
    numpy.ndarray
        Pseudo-power with shape ``steering.shape[:-1]``; it has no
        physical units.

    Raises
    ------
    ValueError
        If ``sources`` is not in ``[1, n)``.
    """
    matrix = _square_matrix(matrix)
    matrix = 0.5 * (matrix + matrix.conj().T)
    size = matrix.shape[0]
    steering = _steering(steering, size)

    if (
        not isinstance(sources, Integral)
        or isinstance(sources, bool)
        or not 1 <= sources < size
    ):
        raise ValueError("sources must be an integer in [1, n)")

    noise = np.linalg.eigh(matrix)[1][:, :size - sources]
    denominator = np.sum(np.abs(steering @ noise.conj()) ** 2, axis=-1)

    return _reciprocal(denominator)


def beam_power(
    matrix: Any,
    steering: Any,
    count: int,
    method: Any,
    strict: bool = True,
    diagonal_loading: float = 0.001,
    music_sources: int = 1,
    attenuation: Any = None,
) -> np.ndarray:
    """
    Return the beam power of the selected estimator.

    Parameters
    ----------
    matrix : array-like
        Cross-spectral matrix ``(n, n)``.
    steering : array-like
        Steering vectors ``(..., n)``.
    count : int
        Number of realizations averaged in ``matrix``.
    method : BeamformingMethod or str
        Estimator.
    strict : bool, optional
        Require more realizations than ``n`` for Capon and MUSIC.
    diagonal_loading : float, optional
        Relative diagonal loading of Capon. Default is 0.001.
    music_sources : int, optional
        Number of sources of MUSIC. Default is 1.
    attenuation : array-like, optional
        Expected coherence of the stations (Bartlett only).

    Returns
    -------
    numpy.ndarray
        Power with shape ``steering.shape[:-1]``.

    Raises
    ------
    ValueError
        If the realizations are too few, or ``attenuation`` is given
        for Capon or MUSIC.
    """
    method = _method(method)
    matrix = _square_matrix(matrix)
    check_realizations(count, matrix.shape[0], method, strict)

    if method is BeamformingMethod.BARTLETT:
        return bartlett_power(matrix, steering, attenuation)

    if attenuation is not None:
        raise ValueError("the coordinate uncertainty needs Bartlett")

    if method is BeamformingMethod.CAPON:
        return capon_power(matrix, steering, diagonal_loading)

    return music_power(matrix, steering, music_sources)
