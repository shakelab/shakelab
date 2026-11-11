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
Geometric response and ambiguity of two-dimensional arrays.

The response depends on the station positions only: no data, estimator
or processing configuration is involved. Coordinates are east (x) and
north (y); azimuths are clockwise from north; wavenumbers are angular
(radians per unit of the coordinates).

For stations at ``r_n`` the normalized array factor and the power
response (geometric ambiguity function) are::

    B(dk) = 1 / N sum_n exp(-i dk . r_n),     A(dk) = |B(dk)|**2,

and, for two wavevectors, ``|a(k1)^H a(k2)|**2 = A(k2 - k1)`` with
normalized steering vectors ``a``. The central lobe of ``A`` therefore
describes the local resolution, its remote maxima the geometric
ambiguities (aliases, grating lobes, strong sidelobes).

:func:`array_response` evaluates ``A`` on a Cartesian grid
(:class:`ResponseGrid`); :func:`analyze_response` extracts the main
lobe, the ambiguity structure and the derived limits
(:class:`ResponseAnalysis`); :func:`analyze_search_domain` checks a
circular search domain; :func:`velocity_limits` converts the
wavenumber limits to phase velocities. The full two-dimensional
response is the primary quantity: scalar limits are summaries and keep
the thresholds they were obtained with.

Coordinate uncertainty, as independent Gaussian covariances ``S_i`` of
the stations, is propagated exactly: with
``q_i = exp(-dk^T S_i dk / 2)`` the expected response is::

    E[A] = (|sum_i q_i exp(-i dk . r_i)|**2 + sum_i (1 - q_i**2)) / N**2,

where the diagonal term accounts for the same error appearing in both
factors of each auto term.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from scipy.ndimage import label, maximum_filter
from scipy.spatial.distance import pdist

from shakelab.exploration.yarray.data import ArrayData

__all__ = [
    "AmbiguityAnalysis",
    "AmbiguityPeak",
    "ArrayGeometry",
    "ArrayLimits",
    "ArrayResponse",
    "MainLobe",
    "ResponseAnalysis",
    "ResponseGrid",
    "SearchDomainAnalysis",
    "VelocityLimits",
    "analyze_response",
    "analyze_search_domain",
    "array_response",
    "velocity_limits",
]

_RESPONSE_LEVEL = 0.5
_AMBIGUITY_LEVEL = 0.5
_AZIMUTH_SIZE = 180
_PEAK_NEIGHBORHOOD = 3
_SAMPLES_PER_RESOLUTION = 12
_EXTENT_FACTOR = 1.25
_MAX_GRID_SIZE = 1201
_LEVEL_TOLERANCE = 1e-12
_NEIGHBORS = (
    (-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)
)


def _readonly(value: Any, dtype: Any = float) -> np.ndarray:
    """Return a read-only copy of an array."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _integer(value: Any, name: str, minimum: int) -> int:
    """Return an integer not smaller than ``minimum``."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")

    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")

    return int(value)


def _level(value: Any, name: str) -> float:
    """Return a relative response level in (0, 1)."""
    value = float(value)

    if not np.isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be between 0 and 1 (excluded)")

    return value


def _coordinates(value: Any) -> np.ndarray:
    """Return validated horizontal coordinates ``(n, 2)``."""
    coordinates = np.asarray(value, dtype=float)

    if coordinates.ndim != 2 or coordinates.shape[1] not in (2, 3):
        raise ValueError("coordinates must be (n_stations, 2 or 3)")

    if coordinates.shape[0] < 2:
        raise ValueError("at least two stations are required")

    coordinates = coordinates[:, :2]

    if not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates must be finite")

    if np.unique(coordinates, axis=0).shape[0] != coordinates.shape[0]:
        raise ValueError("station coordinates must be distinct")

    return coordinates


def _covariances(value: Any, count: int) -> Optional[np.ndarray]:
    """Return validated horizontal covariances ``(n, 2, 2)`` or None."""
    if value is None:
        return None

    covariance = np.asarray(value, dtype=float)

    if covariance.shape == (count, 3, 3):
        covariance = covariance[:, :2, :2]

    if covariance.shape != (count, 2, 2):
        raise ValueError(
            "coordinate covariances must be (n_stations, 2, 2) or "
            "(n_stations, 3, 3)"
        )

    if not np.all(np.isfinite(covariance)):
        raise ValueError("coordinate covariances must be finite")

    if not np.allclose(
        covariance, np.swapaxes(covariance, 1, 2), rtol=0.0, atol=1e-12
    ):
        raise ValueError("coordinate covariances must be symmetric")

    covariance = 0.5 * (covariance + np.swapaxes(covariance, 1, 2))
    eigenvalues = np.linalg.eigvalsh(covariance)
    scale = max(1.0, float(np.max(np.abs(eigenvalues))))

    if np.any(eigenvalues < -8.0 * np.finfo(float).eps * scale):
        raise ValueError("coordinate covariances must be positive "
                         "semidefinite")

    return covariance


def _format(value: Optional[float]) -> str:
    """Return a compact representation of an optional number."""
    return "not found" if value is None else f"{value:.6g}"


# ----------------------------------------------------------------------------
# Geometry and grid
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class ArrayGeometry:
    """
    Summary of the horizontal station geometry.

    Build it with :meth:`from_coordinates` or :meth:`from_array`.

    Attributes
    ----------
    coordinates : numpy.ndarray
        Horizontal coordinates ``(n, 2)`` (east, north).
    centered_coordinates : numpy.ndarray
        Coordinates relative to their centroid.
    centroid : numpy.ndarray
        Centroid of the coordinates.
    covariance : numpy.ndarray
        Population covariance ``(2, 2)`` of the centered coordinates.
    principal_values : numpy.ndarray
        Eigenvalues of ``covariance``, decreasing.
    principal_directions : numpy.ndarray
        Corresponding unit eigenvectors, by columns.
    minimum_spacing, maximum_aperture : float
        Minimum and maximum distance between two stations.
    rank : int
        Rank of the centered coordinates (1 for a linear array).

    Notes
    -----
    ``rayleigh_proxy``, ``alias_period_proxy`` and ``nyquist_proxy``
    are the classical rules of thumb, for diagnostics only: the limits
    of :func:`analyze_response` come from the actual response.
    """

    coordinates: np.ndarray
    centered_coordinates: np.ndarray
    centroid: np.ndarray
    covariance: np.ndarray
    principal_values: np.ndarray
    principal_directions: np.ndarray
    minimum_spacing: float
    maximum_aperture: float
    rank: int

    @classmethod
    def from_coordinates(cls, coordinates: Any) -> ArrayGeometry:
        """
        Return the geometry of station coordinates.

        Parameters
        ----------
        coordinates : array-like
            Coordinates ``(n, 2)`` or ``(n, 3)``; only east and north
            are used.

        Returns
        -------
        ArrayGeometry
            Geometry summary.

        Raises
        ------
        ValueError
            If there are fewer than two stations, or repeated or
            non-finite coordinates.
        """
        coordinates = _coordinates(coordinates)
        centroid = np.mean(coordinates, axis=0)
        centered = coordinates - centroid
        covariance = centered.T @ centered / coordinates.shape[0]
        values, vectors = np.linalg.eigh(covariance)
        order = np.argsort(values)[::-1]
        distances = pdist(coordinates)

        return cls(
            coordinates=_readonly(coordinates),
            centered_coordinates=_readonly(centered),
            centroid=_readonly(centroid),
            covariance=_readonly(covariance),
            principal_values=_readonly(values[order]),
            principal_directions=_readonly(vectors[:, order]),
            minimum_spacing=float(np.min(distances)),
            maximum_aperture=float(np.max(distances)),
            rank=int(np.linalg.matrix_rank(centered)),
        )

    @classmethod
    def from_array(cls, array: ArrayData) -> ArrayGeometry:
        """
        Return the geometry of the stations of an array.

        Parameters
        ----------
        array : ArrayData
            Array data (e.g. from :func:`read_array`).

        Returns
        -------
        ArrayGeometry
            Geometry summary.
        """
        if not isinstance(array, ArrayData):
            raise TypeError("array must be an ArrayData")

        return cls.from_coordinates(array.coordinates)

    @property
    def n_stations(self) -> int:
        """Number of stations."""
        return int(self.coordinates.shape[0])

    @property
    def is_linear(self) -> bool:
        """True for a linear (rank-one) geometry."""
        return self.rank == 1

    @property
    def is_two_dimensional(self) -> bool:
        """True if both horizontal wavevector components are resolved."""
        return self.rank == 2

    @property
    def rayleigh_proxy(self) -> float:
        """Classical resolution proxy ``2 pi / maximum_aperture``."""
        return float(2.0 * np.pi / self.maximum_aperture)

    @property
    def alias_period_proxy(self) -> float:
        """Classical alias period proxy ``2 pi / minimum_spacing``."""
        return float(2.0 * np.pi / self.minimum_spacing)

    @property
    def nyquist_proxy(self) -> float:
        """Classical Nyquist proxy ``pi / minimum_spacing``."""
        return float(np.pi / self.minimum_spacing)

    @property
    def anisotropy_ratio(self) -> float:
        """Ratio of the principal extents (infinite if linear)."""
        major, minor = (float(value) for value in self.principal_values)

        if minor <= 0.0:
            return float("inf")

        return float(np.sqrt(major / minor))

    def info(self) -> None:
        """Print a summary of the geometry."""
        print(
            f"ArrayGeometry: {self.n_stations} stations, rank {self.rank}, "
            f"spacing {self.minimum_spacing:.6g} to "
            f"{self.maximum_aperture:.6g}, anisotropy "
            f"{self.anisotropy_ratio:.4g}"
        )


def _geometry(source: Any) -> ArrayGeometry:
    """Return the geometry of an array, a geometry or coordinates."""
    if isinstance(source, ArrayGeometry):
        return source

    if isinstance(source, ArrayData):
        return ArrayGeometry.from_array(source)

    return ArrayGeometry.from_coordinates(source)


@dataclass(frozen=True)
class ResponseGrid:
    """
    Symmetric Cartesian grid of angular wavenumbers.

    Parameters
    ----------
    maximum : float
        Extent: both axes span ``[-maximum, maximum]``.
    size : int
        Odd number of samples per axis (at least five), so that the
        origin is a grid point.

    Raises
    ------
    ValueError
        If ``maximum`` is not positive or ``size`` is even or too
        small.
    """

    maximum: float
    size: int

    def __post_init__(self) -> None:
        maximum = float(self.maximum)

        if not np.isfinite(maximum) or maximum <= 0.0:
            raise ValueError("the grid maximum must be positive")

        size = _integer(self.size, "size", 5)

        if size % 2 == 0:
            raise ValueError("the grid size must be odd")

        object.__setattr__(self, "maximum", maximum)
        object.__setattr__(self, "size", size)

    @classmethod
    def from_geometry(
        cls,
        geometry: Any,
        extent_factor: float = _EXTENT_FACTOR,
        samples_per_resolution: int = _SAMPLES_PER_RESOLUTION,
        max_size: int = _MAX_GRID_SIZE,
    ) -> ResponseGrid:
        """
        Return a practical grid for a geometry.

        The extent is ``extent_factor * 2 pi / minimum_spacing`` and
        the step ``(2 pi / maximum_aperture) / samples_per_resolution``:
        numerical sampling choices only.

        Parameters
        ----------
        geometry : ArrayGeometry, ArrayData or array-like
            Geometry, array or station coordinates.
        extent_factor : float, optional
            Extent relative to the alias period proxy. Default is 1.25.
        samples_per_resolution : int, optional
            Samples per Rayleigh proxy (at least four). Default is 12.
        max_size : int, optional
            Largest size allowed. Default is 1201.

        Returns
        -------
        ResponseGrid
            Grid.

        Raises
        ------
        ValueError
            If the grid would exceed ``max_size`` (it is never
            coarsened silently, so that numerical and physical
            resolution are not confused): give an explicit grid or a
            larger ``max_size``.
        """
        geometry = _geometry(geometry)
        extent_factor = float(extent_factor)

        if not np.isfinite(extent_factor) or extent_factor <= 0.0:
            raise ValueError("extent_factor must be positive")

        samples = _integer(samples_per_resolution, "samples_per_resolution",
                           4)
        max_size = _integer(max_size, "max_size", 5)

        if max_size % 2 == 0:
            max_size -= 1

        maximum = extent_factor * geometry.alias_period_proxy
        half = int(np.ceil(maximum / (geometry.rayleigh_proxy / samples)))
        size = 2 * half + 1

        if size > max_size:
            raise ValueError(
                f"the automatic response grid needs {size} x {size} "
                f"samples, more than max_size={max_size}; give an "
                "explicit ResponseGrid or a larger max_size"
            )

        return cls(maximum=float(maximum), size=size)

    @property
    def values(self) -> np.ndarray:
        """Wavenumbers of each axis."""
        return np.linspace(-self.maximum, self.maximum, self.size)

    @property
    def step(self) -> float:
        """Grid step."""
        return float(2.0 * self.maximum / (self.size - 1))

    @property
    def center_index(self) -> int:
        """Index of the origin on each axis."""
        return self.size // 2

    def mesh(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the ``kx`` and ``ky`` meshes (rows are ``ky``)."""
        return np.meshgrid(self.values, self.values, indexing="xy")


# ----------------------------------------------------------------------------
# Response
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class ArrayResponse:
    """
    Geometric array response on a Cartesian grid.

    Attributes
    ----------
    geometry : ArrayGeometry
        Station geometry.
    grid : ResponseGrid
        Wavenumber grid.
    factor : numpy.ndarray
        Complex normalized array factor ``(size, size)``, rows ``ky``.
    power : numpy.ndarray
        Power response ``|factor|**2``.
    expected_power : numpy.ndarray or None
        Expected power under the coordinate covariances, if given.
    coordinate_covariances : numpy.ndarray or None
        Horizontal covariances ``(n, 2, 2)`` used for
        ``expected_power``.
    """

    geometry: ArrayGeometry
    grid: ResponseGrid
    factor: np.ndarray
    power: np.ndarray
    expected_power: Optional[np.ndarray] = None
    coordinate_covariances: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        shape = (self.grid.size, self.grid.size)

        for name in ("factor", "power", "expected_power"):
            value = getattr(self, name)

            if value is not None and value.shape != shape:
                raise ValueError(f"{name} does not match the grid")

    @property
    def kx(self) -> np.ndarray:
        """East wavenumbers (columns)."""
        return self.grid.values

    @property
    def ky(self) -> np.ndarray:
        """North wavenumbers (rows)."""
        return self.grid.values

    def selected_power(self, expected: bool = False) -> np.ndarray:
        """
        Return the nominal or the expected power response.

        Parameters
        ----------
        expected : bool, optional
            If True, the response expected under the coordinate
            uncertainties.

        Returns
        -------
        numpy.ndarray
            Power response.

        Raises
        ------
        ValueError
            If ``expected`` is True and no covariances were given.
        """
        if not expected:
            return self.power

        if self.expected_power is None:
            raise ValueError("no coordinate covariances: the expected "
                             "response is not available")

        return self.expected_power

    def info(self) -> None:
        """Print a summary of the response."""
        uncertainty = (
            ", with coordinate uncertainty"
            if self.expected_power is not None else ""
        )
        print(
            f"ArrayResponse: {self.geometry.n_stations} stations, grid "
            f"{self.grid.size} x {self.grid.size}, step "
            f"{self.grid.step:.6g}, extent +/-{self.grid.maximum:.6g}"
            f"{uncertainty}"
        )


def array_response(
    array: Any,
    grid: Optional[ResponseGrid] = None,
    coordinate_covariances: Any = None,
    extent_factor: float = _EXTENT_FACTOR,
    samples_per_resolution: int = _SAMPLES_PER_RESOLUTION,
    max_grid_size: int = _MAX_GRID_SIZE,
    chunk_size: int = 65536,
) -> ArrayResponse:
    """
    Return the geometric response of an array.

    Parameters
    ----------
    array : ArrayData, ArrayGeometry or array-like
        Array, geometry or station coordinates ``(n, 2)`` or ``(n, 3)``
        (east, north[, up]; only the horizontal ones are used).
    grid : ResponseGrid, optional
        Wavenumber grid. Default: :meth:`ResponseGrid.from_geometry`
        with the three following options.
    coordinate_covariances : array-like, optional
        Covariances ``(n, 2, 2)`` or ``(n, 3, 3)`` of the coordinates,
        for the expected response. Default: those of an ArrayData, if
        any.
    extent_factor, samples_per_resolution, max_grid_size : optional
        Options of the automatic grid.
    chunk_size : int, optional
        Grid points evaluated at once. Default is 65536.

    Returns
    -------
    ArrayResponse
        Array factor, power and, with covariances, expected power.

    Notes
    -----
    Coordinates are centered before evaluating the phases, so the
    power does not depend on a translation of the array.
    """
    if (
        coordinate_covariances is None and isinstance(array, ArrayData)
    ):
        coordinate_covariances = array.coordinate_covariances

    geometry = _geometry(array)
    covariance = _covariances(coordinate_covariances, geometry.n_stations)

    if grid is None:
        grid = ResponseGrid.from_geometry(
            geometry, extent_factor, samples_per_resolution, max_grid_size
        )

    elif not isinstance(grid, ResponseGrid):
        raise TypeError("grid must be a ResponseGrid")

    chunk_size = _integer(chunk_size, "chunk_size", 1)
    kx, ky = grid.mesh()
    points = np.column_stack((kx.ravel(), ky.ravel()))
    count = geometry.n_stations
    factor = np.empty(points.shape[0], dtype=complex)
    expected = None if covariance is None else np.empty(points.shape[0])
    centered = geometry.centered_coordinates

    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        vectors = points[start:stop]
        steering = np.exp(-1j * (vectors @ centered.T))
        factor[start:stop] = np.sum(steering, axis=1) / count

        if expected is not None:
            variance = np.einsum(
                "pi,nij,pj->pn", vectors, covariance, vectors,
                optimize=True,
            )
            attenuation = np.exp(-0.5 * np.maximum(variance, 0.0))
            expected[start:stop] = (
                np.abs(np.sum(steering * attenuation, axis=1)) ** 2
                + np.sum(1.0 - attenuation ** 2, axis=1)
            ) / count ** 2

    factor = factor.reshape(grid.size, grid.size)
    power = np.abs(factor) ** 2
    center = grid.center_index
    factor[center, center] = 1.0
    power[center, center] = 1.0

    if expected is not None:
        expected = expected.reshape(grid.size, grid.size)
        expected[center, center] = 1.0
        expected = _readonly(np.clip(expected, 0.0, 1.0))

    return ArrayResponse(
        geometry=geometry,
        grid=grid,
        factor=_readonly(factor, complex),
        power=_readonly(power),
        expected_power=expected,
        coordinate_covariances=(
            None if covariance is None else _readonly(covariance)
        ),
    )


# ----------------------------------------------------------------------------
# Analysis results
# ----------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class MainLobe:
    """
    Central lobe of the response at a given level.

    Attributes
    ----------
    level : float
        Relative power defining the lobe.
    mask : numpy.ndarray
        Grid points of the connected region containing the origin.
    area : float
        Area of the lobe (grid-cell approximation).
    azimuths : numpy.ndarray
        Axis azimuths in ``[0, 180)`` degrees.
    radii : numpy.ndarray
        Distance of the first level crossing along each azimuth.
    """

    level: float
    mask: np.ndarray
    area: float
    azimuths: np.ndarray
    radii: np.ndarray

    @property
    def equivalent_radius(self) -> float:
        """Radius of the circle with the area of the lobe."""
        return float(np.sqrt(self.area / np.pi))

    @property
    def resolution_separation(self) -> np.ndarray:
        """Full width of the lobe along each azimuth (twice the radii)."""
        return 2.0 * self.radii

    @property
    def best_radius(self) -> float:
        """Smallest radius."""
        return float(np.nanmin(self.radii))

    @property
    def worst_radius(self) -> float:
        """Largest radius."""
        return float(np.nanmax(self.radii))

    @property
    def best_resolution_separation(self) -> float:
        """Smallest full width."""
        return 2.0 * self.best_radius

    @property
    def worst_resolution_separation(self) -> float:
        """Largest full width."""
        return 2.0 * self.worst_radius


@dataclass(frozen=True)
class AmbiguityPeak:
    """
    Local maximum of the response outside the main lobe.

    Attributes
    ----------
    rank : int
        Rank by decreasing power, from 1.
    kx, ky : float
        Wavevector displacement of the maximum.
    power : float
        Relative power.

    Notes
    -----
    ``azimuth`` is that of the displacement, not a propagation
    direction.
    """

    rank: int
    kx: float
    ky: float
    power: float

    @property
    def wavenumber(self) -> float:
        """Length of the displacement."""
        return float(np.hypot(self.kx, self.ky))

    @property
    def azimuth(self) -> float:
        """Azimuth of the displacement, clockwise from north."""
        return float(np.mod(np.degrees(np.arctan2(self.kx, self.ky)), 360.0))

    @property
    def power_db(self) -> float:
        """Relative power in decibels."""
        if self.power <= 0.0:
            return float("-inf")

        return float(10.0 * np.log10(self.power))


@dataclass(frozen=True)
class AmbiguityAnalysis:
    """
    Remote ambiguity structure at a given level.

    Attributes
    ----------
    level : float
        Relative power threshold.
    peaks : tuple of AmbiguityPeak
        All local maxima outside the main lobe, by decreasing power.
    significant_peaks : tuple of AmbiguityPeak
        Maxima with power at least ``level``.
    peak_sidelobe_level : float
        Power of the strongest maximum (0 if none).
    peak_distance : float or None
        Distance of the nearest significant maximum.
    onset_distance : float or None
        Shortest distance to a remote connected region with power at
        least ``level``: it controls the search-domain limits.
    ambiguity_margin : float
        ``1 - peak_sidelobe_level``.
    remote_region_count : int
        Number of remote regions above ``level``.
    """

    level: float
    peaks: tuple
    significant_peaks: tuple
    peak_sidelobe_level: float
    peak_distance: Optional[float]
    onset_distance: Optional[float]
    ambiguity_margin: float
    remote_region_count: int

    @property
    def aliasing_radius(self) -> Optional[float]:
        """Half the onset distance (None if no remote region)."""
        if self.onset_distance is None:
            return None

        return 0.5 * self.onset_distance


@dataclass(frozen=True, eq=False)
class ArrayLimits:
    """
    Resolution and ambiguity limits derived from the response.

    Attributes
    ----------
    resolution_level, ambiguity_level : float
        Levels of the main lobe and of the remote regions.
    azimuths, main_lobe_radii, resolution_separations : numpy.ndarray
        Directional radii and full widths of the main lobe.
    resolution_radius : float
        Largest radius of the main lobe: the smallest usable
        wavenumber.
    resolution_separation : float
        Twice ``resolution_radius``.
    ambiguity_peak_distance, ambiguity_onset_distance : float or None
        See :class:`AmbiguityAnalysis`.
    aliasing_radius : float or None
        Half the onset distance: the largest wavenumber of a
        full-azimuth circular search free of ambiguity.
    rayleigh_proxy, alias_period_proxy, nyquist_proxy : float
        Classical proxies, for comparison.
    """

    resolution_level: float
    ambiguity_level: float
    azimuths: np.ndarray
    main_lobe_radii: np.ndarray
    resolution_separations: np.ndarray
    resolution_radius: float
    resolution_separation: float
    ambiguity_peak_distance: Optional[float]
    ambiguity_onset_distance: Optional[float]
    aliasing_radius: Optional[float]
    rayleigh_proxy: float
    alias_period_proxy: float
    nyquist_proxy: float

    @property
    def has_finite_upper_limit(self) -> bool:
        """True if a remote ambiguity was found."""
        return self.aliasing_radius is not None

    @property
    def usable_interval(self) -> tuple:
        """Usable wavenumbers ``(resolution_radius, aliasing_radius)``."""
        return (self.resolution_radius, self.aliasing_radius)


@dataclass(frozen=True, eq=False)
class VelocityLimits:
    """
    Phase-velocity limits at given frequencies.

    Attributes
    ----------
    frequencies : numpy.ndarray
        Frequencies in Hz.
    minimum, maximum : numpy.ndarray
        Velocities of the upper (aliasing) and of the lower
        (resolution) wavenumber limit.
    resolution_level, ambiguity_level : float
        Levels of the limits.
    """

    frequencies: np.ndarray
    minimum: np.ndarray
    maximum: np.ndarray
    resolution_level: float
    ambiguity_level: float


@dataclass(frozen=True, eq=False)
class ResponseAnalysis:
    """
    Main lobe, ambiguities and limits of a response.

    Attributes
    ----------
    response : ArrayResponse
        Analysed response.
    main_lobe : MainLobe
        Central lobe.
    ambiguity : AmbiguityAnalysis
        Remote ambiguity structure.
    limits : ArrayLimits
        Derived limits.
    expected : bool
        True if the expected response was analysed.
    """

    response: ArrayResponse
    main_lobe: MainLobe
    ambiguity: AmbiguityAnalysis
    limits: ArrayLimits
    expected: bool = False

    def velocity_limits(self, frequencies: Any) -> VelocityLimits:
        """Return the velocity limits (see :func:`velocity_limits`)."""
        return velocity_limits(frequencies, self.limits)

    def info(self) -> None:
        """Print a summary of the analysis."""
        geometry = self.response.geometry
        ambiguity = self.ambiguity
        limits = self.limits
        print(
            f"ResponseAnalysis: {geometry.n_stations} stations, spacing "
            f"{geometry.minimum_spacing:.6g} to "
            f"{geometry.maximum_aperture:.6g}, anisotropy "
            f"{geometry.anisotropy_ratio:.4g}"
            + (" (expected response)" if self.expected else "")
        )
        print(
            f"  main lobe at {self.main_lobe.level:g}: radius "
            f"{self.main_lobe.best_radius:.6g} to "
            f"{self.main_lobe.worst_radius:.6g}"
        )
        print(
            f"  ambiguity at {ambiguity.level:g}: peak sidelobe "
            f"{ambiguity.peak_sidelobe_level:.4g}, "
            f"{ambiguity.remote_region_count} remote regions, nearest "
            f"peak {_format(ambiguity.peak_distance)}, onset "
            f"{_format(ambiguity.onset_distance)}"
        )
        if limits.aliasing_radius is None:
            upper = "with no upper limit in the grid"
        elif limits.aliasing_radius <= limits.resolution_radius:
            upper = (f"to {limits.aliasing_radius:.6g}, empty: the main "
                     "lobe is wider than the alias-free domain")
        else:
            upper = f"to {limits.aliasing_radius:.6g}"
        print(
            f"  usable wavenumbers: from {limits.resolution_radius:.6g} "
            f"{upper} (proxies: Rayleigh "
            f"{limits.rayleigh_proxy:.6g}, Nyquist "
            f"{limits.nyquist_proxy:.6g})"
        )


@dataclass(frozen=True)
class SearchDomainAnalysis:
    """
    Ambiguity of a full-azimuth circular search domain.

    The differences of the wavevectors with ``|k| <= K`` fill the disk
    ``|dk| <= 2 K``; only remote regions above the level are
    considered.

    Attributes
    ----------
    maximum_wavenumber : float
        Radius ``K`` of the search domain.
    difference_radius : float
        ``2 K``.
    ambiguity_level : float
        Level of the remote regions.
    maximum_ambiguity : float
        Largest response of the reachable remote regions (0 if none):
        a squared steering-vector correlation, not a probability.
    ambiguity_margin : float
        ``1 - maximum_ambiguity``.
    is_safe : bool
        True if no remote region is reachable.
    ambiguity_peaks : tuple of AmbiguityPeak
        Reachable significant maxima.
    """

    maximum_wavenumber: float
    difference_radius: float
    ambiguity_level: float
    maximum_ambiguity: float
    ambiguity_margin: float
    is_safe: bool
    ambiguity_peaks: tuple

    def info(self) -> None:
        """Print a summary of the analysis."""
        status = "safe" if self.is_safe else "ambiguous"
        print(
            f"SearchDomainAnalysis: k <= {self.maximum_wavenumber:.6g} is "
            f"{status} at level {self.ambiguity_level:g} (largest "
            f"reachable ambiguity {self.maximum_ambiguity:.4g}, "
            f"{len(self.ambiguity_peaks)} peaks)"
        )


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------


def _main_lobe_mask(power: np.ndarray, grid: ResponseGrid,
                    level: float) -> np.ndarray:
    """Return the region above ``level`` connected to the origin."""
    components, _ = label(power >= level - _LEVEL_TOLERANCE,
                          structure=np.ones((3, 3), dtype=int))
    center = grid.center_index
    identifier = int(components[center, center])

    if identifier == 0:
        raise ValueError("the origin is below the main-lobe level")

    mask = components == identifier

    if mask[0].any() or mask[-1].any() or mask[:, 0].any() or \
            mask[:, -1].any():
        raise ValueError("the main lobe reaches the grid boundary: "
                         "increase ResponseGrid.maximum")

    return mask


def _bilinear(values: np.ndarray, grid: ResponseGrid, kx: Any,
              ky: Any) -> np.ndarray:
    """Return a bilinear interpolation of a grid."""
    x = np.clip((np.asarray(kx, dtype=float) + grid.maximum) / grid.step,
                0.0, grid.size - 1.0)
    y = np.clip((np.asarray(ky, dtype=float) + grid.maximum) / grid.step,
                0.0, grid.size - 1.0)
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = np.minimum(x0 + 1, grid.size - 1)
    y1 = np.minimum(y0 + 1, grid.size - 1)
    wx = x - x0
    wy = y - y0

    return (
        (1.0 - wx) * (1.0 - wy) * values[y0, x0]
        + wx * (1.0 - wy) * values[y0, x1]
        + (1.0 - wx) * wy * values[y1, x0]
        + wx * wy * values[y1, x1]
    )


def _radii(power: np.ndarray, grid: ResponseGrid, level: float,
           azimuths: np.ndarray) -> np.ndarray:
    """Return the first level crossing along each azimuth."""
    step = grid.step / 4.0
    radii = np.full(azimuths.shape, np.nan)

    for index, azimuth in enumerate(azimuths):
        angle = np.radians(azimuth)
        direction = np.array([np.sin(angle), np.cos(angle)])
        extent = np.divide(
            grid.maximum, np.abs(direction), out=np.full(2, np.inf),
            where=np.abs(direction) > 0.0,
        )
        radial = np.linspace(0.0, float(np.min(extent)),
                             int(np.ceil(np.min(extent) / step)) + 1)
        sampled = _bilinear(power, grid, radial * direction[0],
                            radial * direction[1])
        below = np.flatnonzero(sampled <= level)

        if below.size == 0:
            continue

        upper = int(below[0])

        if upper == 0:
            radii[index] = 0.0
            continue

        r0, r1 = radial[upper - 1], radial[upper]
        p0, p1 = sampled[upper - 1], sampled[upper]

        if np.isclose(p0, p1):
            radii[index] = r1
        else:
            fraction = float(np.clip((level - p0) / (p1 - p0), 0.0, 1.0))
            radii[index] = r0 + fraction * (r1 - r0)

    if not np.all(np.isfinite(radii)):
        raise ValueError("a main-lobe level crossing is outside the grid: "
                         "increase ResponseGrid.maximum")

    return radii


def _peak_offset(minus: float, center: float, plus: float) -> float:
    """Return the parabolic sub-sample offset of a maximum."""
    denominator = minus - 2.0 * center + plus

    if (
        not np.isfinite(denominator) or denominator >= 0.0
        or np.isclose(denominator, 0.0)
    ):
        return 0.0

    return float(np.clip(0.5 * (minus - plus) / denominator, -0.5, 0.5))


def _ambiguity_peaks(power: np.ndarray, grid: ResponseGrid,
                     lobe: np.ndarray, neighborhood: int,
                     minimum_power: float) -> tuple:
    """Return the local maxima outside the main lobe."""
    neighborhood = _integer(neighborhood, "peak_neighborhood", 3)

    if neighborhood % 2 == 0:
        raise ValueError("peak_neighborhood must be odd")

    minimum_power = float(minimum_power)

    if not np.isfinite(minimum_power) or not 0.0 <= minimum_power < 1.0:
        raise ValueError("minimum_peak_power must be in [0, 1)")

    filtered = maximum_filter(power, size=neighborhood, mode="nearest")
    candidates = (
        np.isclose(power, filtered, rtol=0.0, atol=32.0 * np.finfo(float).eps)
        & (power >= minimum_power) & ~lobe
    )
    candidates[[0, -1], :] = False
    candidates[:, [0, -1]] = False
    axis = grid.values
    peaks = []

    for row, column in np.argwhere(candidates):
        dx = _peak_offset(power[row, column - 1], power[row, column],
                          power[row, column + 1])
        dy = _peak_offset(power[row - 1, column], power[row, column],
                          power[row + 1, column])
        kx = float(axis[column] + dx * grid.step)
        ky = float(axis[row] + dy * grid.step)
        value = float(_bilinear(power, grid, np.array([kx]),
                                np.array([ky]))[0])
        peaks.append((kx, ky, value))

    peaks.sort(key=lambda item: item[2], reverse=True)

    return tuple(
        AmbiguityPeak(rank=index + 1, kx=kx, ky=ky, power=value)
        for index, (kx, ky, value) in enumerate(peaks)
    )


def _remote_mask(power: np.ndarray, grid: ResponseGrid,
                 level: float) -> tuple[np.ndarray, int]:
    """Return the remote regions above ``level`` and their number."""
    components, _ = label(power >= level - _LEVEL_TOLERANCE,
                          structure=np.ones((3, 3), dtype=int))
    center = grid.center_index
    central = int(components[center, center])

    if central == 0:
        raise ValueError("the origin is below the ambiguity level")

    remote = (components != 0) & (components != central)

    return remote, int(np.unique(components[remote]).size)


def _onset_distance(power: np.ndarray, grid: ResponseGrid,
                    remote: np.ndarray, level: float) -> Optional[float]:
    """Return the sub-grid distance of the nearest remote region."""
    indices = np.argwhere(remote)

    if indices.size == 0:
        return None

    axis = grid.values
    best = float("inf")

    for row, column in indices:
        x1, y1 = axis[column], axis[row]
        p1 = power[row, column]

        if np.hypot(x1, y1) >= best:
            continue

        for dr, dc in _NEIGHBORS:
            nr, nc = row + dr, column + dc

            if not (0 <= nr < grid.size and 0 <= nc < grid.size):
                continue

            if remote[nr, nc]:
                continue

            x0, y0 = axis[nc], axis[nr]
            p0 = power[nr, nc]

            if np.hypot(x0, y0) >= np.hypot(x1, y1) or p0 >= level or \
                    p1 - p0 <= 0.0:
                continue

            fraction = float(np.clip((level - p0) / (p1 - p0), 0.0, 1.0))
            best = min(best, float(np.hypot(x0 + fraction * (x1 - x0),
                                            y0 + fraction * (y1 - y0))))

    if np.isfinite(best):
        return best

    return float(np.min(np.hypot(axis[indices[:, 1]], axis[indices[:, 0]])))


def _two_dimensional(response: ArrayResponse) -> None:
    """Check the type and the rank of a response."""
    if not isinstance(response, ArrayResponse):
        raise TypeError("response must be an ArrayResponse")

    if not response.geometry.is_two_dimensional:
        raise ValueError(
            "two-dimensional limits cannot be determined for a linear "
            "array (the response itself is valid)"
        )


def analyze_response(
    response: ArrayResponse,
    resolution_level: float = _RESPONSE_LEVEL,
    ambiguity_level: float = _AMBIGUITY_LEVEL,
    azimuth_size: int = _AZIMUTH_SIZE,
    peak_neighborhood: int = _PEAK_NEIGHBORHOOD,
    minimum_peak_power: float = 0.0,
    expected: bool = False,
) -> ResponseAnalysis:
    """
    Return the main lobe, the ambiguities and the limits of a response.

    Parameters
    ----------
    response : ArrayResponse
        Response of a two-dimensional (rank-two) array.
    resolution_level : float, optional
        Relative power of the main lobe. Default is 0.5.
    ambiguity_level : float, optional
        Relative power of the remote regions. Default is 0.5.
    azimuth_size : int, optional
        Axis azimuths of the directional radii (at least four).
        Default is 180.
    peak_neighborhood : int, optional
        Odd size of the neighbourhood of the local maxima. Default 3.
    minimum_peak_power : float, optional
        Smallest power of the listed maxima. Default is 0.
    expected : bool, optional
        If True, analyse the response expected under the coordinate
        uncertainties.

    Returns
    -------
    ResponseAnalysis
        Main lobe, ambiguity structure and limits.

    Raises
    ------
    ValueError
        If the array is linear, a level is out of range or the main
        lobe reaches the grid boundary.
    """
    _two_dimensional(response)
    resolution_level = _level(resolution_level, "resolution_level")
    ambiguity_level = _level(ambiguity_level, "ambiguity_level")
    power = response.selected_power(expected)
    grid = response.grid
    azimuths = np.linspace(
        0.0, 180.0, _integer(azimuth_size, "azimuth_size", 4),
        endpoint=False,
    )
    mask = _main_lobe_mask(power, grid, resolution_level)
    main_lobe = MainLobe(
        level=resolution_level,
        mask=_readonly(mask, bool),
        area=float(np.count_nonzero(mask) * grid.step ** 2),
        azimuths=_readonly(azimuths),
        radii=_readonly(_radii(power, grid, resolution_level, azimuths)),
    )
    peaks = _ambiguity_peaks(power, grid, mask, peak_neighborhood,
                             minimum_peak_power)
    significant = tuple(
        peak for peak in peaks if peak.power >= ambiguity_level
    )
    sidelobe = float(peaks[0].power) if peaks else 0.0
    peak_distance = (
        min(peak.wavenumber for peak in significant) if significant
        else None
    )
    remote, regions = _remote_mask(power, grid, ambiguity_level)
    onset = _onset_distance(power, grid, remote, ambiguity_level)
    ambiguity = AmbiguityAnalysis(
        level=ambiguity_level,
        peaks=peaks,
        significant_peaks=significant,
        peak_sidelobe_level=sidelobe,
        peak_distance=peak_distance,
        onset_distance=onset,
        ambiguity_margin=1.0 - sidelobe,
        remote_region_count=regions,
    )
    geometry = response.geometry
    limits = ArrayLimits(
        resolution_level=resolution_level,
        ambiguity_level=ambiguity_level,
        azimuths=main_lobe.azimuths,
        main_lobe_radii=main_lobe.radii,
        resolution_separations=_readonly(main_lobe.resolution_separation),
        resolution_radius=main_lobe.worst_radius,
        resolution_separation=main_lobe.worst_resolution_separation,
        ambiguity_peak_distance=peak_distance,
        ambiguity_onset_distance=onset,
        aliasing_radius=ambiguity.aliasing_radius,
        rayleigh_proxy=geometry.rayleigh_proxy,
        alias_period_proxy=geometry.alias_period_proxy,
        nyquist_proxy=geometry.nyquist_proxy,
    )

    return ResponseAnalysis(response=response, main_lobe=main_lobe,
                            ambiguity=ambiguity, limits=limits,
                            expected=bool(expected))


def analyze_search_domain(
    response: ArrayResponse,
    maximum_wavenumber: float,
    ambiguity_level: float = _AMBIGUITY_LEVEL,
    peak_neighborhood: int = _PEAK_NEIGHBORHOOD,
    minimum_peak_power: float = 0.0,
    expected: bool = False,
) -> SearchDomainAnalysis:
    """
    Check the ambiguity of a full-azimuth circular search domain.

    Parameters
    ----------
    response : ArrayResponse
        Response of a two-dimensional array.
    maximum_wavenumber : float
        Radius ``K`` of the searched wavenumbers.
    ambiguity_level : float, optional
        Relative power of the remote regions. Default is 0.5.
    peak_neighborhood, minimum_peak_power : optional
        Options of the listed maxima (see :func:`analyze_response`).
    expected : bool, optional
        If True, use the expected response.

    Returns
    -------
    SearchDomainAnalysis
        Reachable ambiguity.

    Raises
    ------
    ValueError
        If the grid does not contain the disk ``|dk| <= 2 K``.
    """
    _two_dimensional(response)
    maximum_wavenumber = float(maximum_wavenumber)

    if not np.isfinite(maximum_wavenumber) or maximum_wavenumber <= 0.0:
        raise ValueError("maximum_wavenumber must be positive")

    ambiguity_level = _level(ambiguity_level, "ambiguity_level")
    radius = 2.0 * maximum_wavenumber
    grid = response.grid

    if radius > grid.maximum - grid.step:
        raise ValueError("the grid does not contain the differences of "
                         "the search domain: increase "
                         "ResponseGrid.maximum")

    power = response.selected_power(expected)
    remote, _ = _remote_mask(power, grid, ambiguity_level)
    kx, ky = grid.mesh()
    reachable = remote & (np.hypot(kx, ky) <= radius)
    largest = float(np.max(power[reachable])) if reachable.any() else 0.0
    peaks = _ambiguity_peaks(
        power, grid, _main_lobe_mask(power, grid, ambiguity_level),
        peak_neighborhood, minimum_peak_power,
    )

    return SearchDomainAnalysis(
        maximum_wavenumber=maximum_wavenumber,
        difference_radius=radius,
        ambiguity_level=ambiguity_level,
        maximum_ambiguity=largest,
        ambiguity_margin=1.0 - largest,
        is_safe=not bool(reachable.any()),
        ambiguity_peaks=tuple(
            peak for peak in peaks
            if peak.power >= ambiguity_level and peak.wavenumber <= radius
        ),
    )


def velocity_limits(frequencies: Any, limits: Any) -> VelocityLimits:
    """
    Convert the wavenumber limits to phase velocities.

    With ``c = 2 pi f / k``, the aliasing radius gives the minimum
    and the resolution radius the maximum velocity.

    Parameters
    ----------
    frequencies : array-like
        Increasing positive frequencies in Hz.
    limits : ArrayLimits or ResponseAnalysis
        Limits of :func:`analyze_response`.

    Returns
    -------
    VelocityLimits
        Velocity limits.

    Raises
    ------
    ValueError
        If no remote ambiguity was found in the grid (no upper
        wavenumber limit): extend the grid or change the level.
    """
    if isinstance(limits, ResponseAnalysis):
        limits = limits.limits

    if not isinstance(limits, ArrayLimits):
        raise TypeError("limits must be an ArrayLimits or a "
                        "ResponseAnalysis")

    frequency = np.asarray(frequencies, dtype=float)

    if (
        frequency.ndim != 1 or frequency.size == 0
        or not np.all(np.isfinite(frequency)) or np.any(frequency <= 0.0)
        or np.any(np.diff(frequency) <= 0.0)
    ):
        raise ValueError("frequencies must be increasing positive values")

    if limits.aliasing_radius is None:
        raise ValueError("no upper wavenumber limit: no remote ambiguity "
                         "was found in the response grid")

    return VelocityLimits(
        frequencies=_readonly(frequency),
        minimum=_readonly(2.0 * np.pi * frequency / limits.aliasing_radius),
        maximum=_readonly(2.0 * np.pi * frequency / limits.resolution_radius),
        resolution_level=limits.resolution_level,
        ambiguity_level=limits.ambiguity_level,
    )
