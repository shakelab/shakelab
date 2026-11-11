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
Figures of array processing (:mod:`shakelab.exploration.yarray`).

- :func:`plot_dispersion`: spectral maxima in frequency and phase
  velocity;
- :func:`plot_ellipticity`: signed Rayleigh ellipticity of the maxima
  versus frequency;
- :func:`plot_spectrum`: one beam-power spectrum (or its standard
  error, or the ellipticity of the joint processing) in wavenumber
  space or in azimuth and radial grid coordinates;
- :func:`plot_density`: density of the spectral maxima versus
  frequency (:mod:`shakelab.exploration.yarray.statistics`), with its
  peaks or summary;
- :func:`plot_array_geometry`, :func:`plot_array_response`,
  :func:`plot_resolution`, :func:`plot_search_domain` and
  :func:`plot_velocity_limits`: station layout, geometric response
  and limits (:mod:`shakelab.exploration.yarray.response`);
- :func:`plot_performance` and :func:`plot_performance_sweep`:
  Monte Carlo estimates and sweeps
  (:mod:`shakelab.exploration.yarray.performance`);
- :func:`plot_operational`: frequency-velocity envelope of an array
  (:mod:`shakelab.exploration.yarray.operational`).

The functions draw on a given matplotlib axis or on a new figure and
return ``(figure, axis)``; they never modify the results. matplotlib is
an optional dependency of ShakeLab, imported when a function is called.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Optional

import numpy as np

from shakelab.exploration.yarray.beamforming import (
    BeamformingResult,
    FrequencySpectrum,
)
from shakelab.exploration.yarray.config import GridQuantity
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.operational import OperationalAssessment
from shakelab.exploration.yarray.performance import (
    OutcomeClass,
    PerformanceResult,
    PerformanceSweepResult,
)
from shakelab.exploration.yarray.response import (
    ArrayGeometry,
    ArrayResponse,
    MainLobe,
    ResponseAnalysis,
    SearchDomainAnalysis,
    VelocityLimits,
)
from shakelab.exploration.yarray.statistics import (
    DensityCoordinate,
    DensityPeaks,
    DensitySummary,
    MaximumDensity,
)

__all__ = [
    "plot_array_geometry",
    "plot_array_response",
    "plot_density",
    "plot_dispersion",
    "plot_ellipticity",
    "plot_operational",
    "plot_performance",
    "plot_performance_sweep",
    "plot_resolution",
    "plot_search_domain",
    "plot_spectrum",
    "plot_velocity_limits",
]

# Quantities shown by :func:`plot_spectrum`.
_QUANTITIES = ("power", "standard_error", "ellipticity")

_RADIAL_LABELS = {
    GridQuantity.VELOCITY: "Phase velocity (m/s)",
    GridQuantity.SLOWNESS: "Slowness (s/m)",
    GridQuantity.WAVENUMBER: "Angular wavenumber (rad/m)",
}


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _pyplot() -> Any:
    """Return matplotlib.pyplot, with a hint if it is not installed."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "the figures need matplotlib: pip install matplotlib"
        ) from error

    return plt


def _axes(ax: Any, figsize: tuple) -> tuple:
    """Return the figure and axis to draw on."""
    if ax is not None:
        return ax.figure, ax

    return _pyplot().subplots(figsize=figsize)


def _check_result(result: Any) -> None:
    """Raise TypeError for anything but a BeamformingResult."""
    if not isinstance(result, BeamformingResult):
        raise TypeError("result must be a BeamformingResult")


def _component(result: BeamformingResult, component: Optional[str]) -> str:
    """Component to plot: given, or the only one processed."""
    if component is not None:
        if component not in result.components:
            raise ValueError(
                f"component {component!r} was not processed "
                f"({', '.join(result.components)})"
            )

        return component

    if len(result.components) > 1:
        raise ValueError(
            "choose the component: " + ", ".join(result.components)
        )

    return result.components[0]


def _maxima(
    result: BeamformingResult,
    component: str,
    ranks: Optional[Iterable],
) -> list:
    """Maxima of a component and ranks."""
    ranks = None if ranks is None else {int(rank) for rank in ranks}
    return [
        item for item in result.maxima
        if item.component == component
        and (ranks is None or item.rank in ranks)
    ]


def _plain_log_axis(axis: Any, which: str) -> None:
    """Logarithmic axis labelled 1, 2, 5 per decade as plain numbers."""
    from matplotlib.ticker import FuncFormatter, LogLocator

    def label(value: float, position: Any) -> str:
        mantissa = value / 10.0 ** np.floor(np.log10(value) + 1e-9)
        return f"{value:g}" if round(mantissa) in (1, 2, 5) else ""

    getattr(axis, f"set_{which}scale")("log")
    ticks = getattr(axis, f"{which}axis")
    ticks.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
    ticks.set_major_formatter(FuncFormatter(label))
    ticks.set_minor_formatter(FuncFormatter(lambda value, position: ""))


def _is_geometric(values: np.ndarray) -> bool:
    """Whether positive values are geometrically spaced."""
    if values.size < 3 or np.any(values <= 0.0):
        return False

    ratios = values[1:] / values[:-1]
    return bool(np.allclose(ratios, ratios[0], rtol=1e-6))


def _scatter_by_rank(axis: Any, maxima: list, value: Any) -> None:
    """Draw the maxima with one series per rank."""
    ranks = sorted({item.rank for item in maxima})

    for rank in ranks:
        selected = [item for item in maxima if item.rank == rank]
        axis.scatter(
            [item.frequency for item in selected],
            [value(item) for item in selected],
            s=12, alpha=0.7, label=f"rank {rank}",
        )

    if len(ranks) > 1:
        axis.legend()


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def plot_dispersion(
    result: BeamformingResult,
    component: Optional[str] = None,
    ranks: Optional[Iterable] = None,
    log_frequency: bool = True,
    log_velocity: bool = True,
    ax: Any = None,
) -> tuple:
    """
    Plot the spectral maxima in frequency and phase velocity.

    Parameters
    ----------
    result : BeamformingResult
        Result with maxima.
    component : str, optional
        Component; required if several were processed.
    ranks : iterable of int, optional
        Ranks to show. Default: all.
    log_frequency, log_velocity : bool, optional
        Logarithmic axes. Default True.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.

    Raises
    ------
    TypeError
        If ``result`` is not a BeamformingResult.
    ValueError
        If the component is ambiguous or not processed, or no maximum
        is selected.
    """
    _check_result(result)
    component = _component(result, component)
    maxima = _maxima(result, component, ranks)

    if not maxima:
        raise ValueError("no maximum matches the selection")

    figure, axis = _axes(ax, (8.0, 5.5))
    _scatter_by_rank(axis, maxima, lambda item: item.velocity)

    if log_frequency:
        _plain_log_axis(axis, "x")

    if log_velocity:
        _plain_log_axis(axis, "y")

    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel("Phase velocity (m/s)")
    axis.set_title(f"Dispersion, {component.replace('_', ' ')}")
    axis.grid(True, which="both", alpha=0.25)

    return figure, axis


def plot_ellipticity(
    result: BeamformingResult,
    ranks: Optional[Iterable] = None,
    angle: bool = False,
    log_frequency: bool = True,
    ax: Any = None,
) -> tuple:
    """
    Plot the signed Rayleigh ellipticity of the maxima.

    Parameters
    ----------
    result : BeamformingResult
        Result of the Rayleigh joint processing.
    ranks : iterable of int, optional
        Ranks to show. Default: all.
    angle : bool, optional
        Show the ellipticity angle ``atan(e)`` in degrees instead of
        ``e``. Default False.
    log_frequency : bool, optional
        Logarithmic frequency axis. Default True.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.

    Raises
    ------
    ValueError
        If the result has no ellipticity.
    """
    _check_result(result)
    maxima = [
        item for item in _maxima(result, "rayleigh_joint", ranks)
        if item.ellipticity is not None
    ]

    if not maxima:
        raise ValueError("the result has no ellipticity (rayleigh_joint)")

    figure, axis = _axes(ax, (8.0, 5.5))
    _scatter_by_rank(
        axis, maxima,
        (lambda item: item.ellipticity_angle) if angle
        else (lambda item: item.ellipticity),
    )
    axis.axhline(0.0, color="0.5", linewidth=0.8)

    if log_frequency:
        _plain_log_axis(axis, "x")

    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel(
        "Ellipticity angle (degrees)" if angle
        else "Signed ellipticity (radial / vertical)"
    )
    axis.set_title("Rayleigh ellipticity")
    axis.grid(True, which="both", alpha=0.25)

    return figure, axis


def plot_spectrum(
    result: BeamformingResult,
    frequency: float,
    block: int = 0,
    component: Optional[str] = None,
    view: str = "kxky",
    quantity: str = "power",
    normalize: bool = True,
    ranks: Optional[Iterable] = None,
    ax: Any = None,
) -> tuple:
    """
    Plot one spectrum of a block.

    Parameters
    ----------
    result : BeamformingResult
        Result with spectra (``output`` ``"spectra"`` or ``"both"``).
    frequency : float
        Frequency in Hz (the nearest analysis frequency).
    block : int, optional
        Covariance block. Default 0.
    component : str, optional
        Component; required if several were processed.
    view : {"kxky", "azimuth"}, optional
        Wavenumber plane (default), or azimuth and radial values of the
        grid.
    quantity : {"power", "standard_error", "ellipticity"}, optional
        Value shown. Default is the power.
    normalize : bool, optional
        Divide the power by its maximum. Default True.
    ranks : iterable of int, optional
        Maxima to mark. Default: all those of the block.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.

    Raises
    ------
    ValueError
        If the result has no spectra or the requested quantity, the
        options are unknown, or the block does not exist.
    """
    _check_result(result)

    if view not in ("kxky", "azimuth"):
        raise ValueError("view must be 'kxky' or 'azimuth'")

    if quantity not in _QUANTITIES:
        raise ValueError(f"quantity must be one of {_QUANTITIES}")

    component = _component(result, component)
    spectrum = result.spectrum(frequency, component)

    if not 0 <= block < spectrum.nblocks:
        raise ValueError(
            f"block {block} does not exist ({spectrum.nblocks} blocks at "
            f"{spectrum.frequency:g} Hz)"
        )

    values, label = _spectrum_values(
        spectrum, block, quantity, normalize, result.power_units
    )
    maxima = [
        item for item in _maxima(result, component, ranks)
        if item.frequency == spectrum.frequency and item.block == block
    ]
    figure, axis = _axes(ax, (7.0, 6.0) if view == "kxky" else (8.0, 5.5))
    azimuth = np.append(spectrum.azimuths, spectrum.azimuths[0] + 360.0)
    closed = np.concatenate([values, values[:1]])

    if view == "kxky":
        k = spectrum.wavenumbers()
        angle = np.radians(azimuth)[:, None]
        image = axis.pcolormesh(
            np.sin(angle) * k, np.cos(angle) * k, closed,
            shading="gouraud",
        )
        axis.plot([item.kx for item in maxima], [item.ky for item in maxima],
                  "wx", markersize=8)
        axis.set_aspect("equal")
        axis.set_xlabel("East wavenumber (rad/m)")
        axis.set_ylabel("North wavenumber (rad/m)")
    else:
        radial = spectrum.grid_values
        image = axis.pcolormesh(azimuth, radial, closed.T, shading="auto")
        axis.plot([item.azimuth for item in maxima],
                  [_radial_value(item, spectrum.quantity) for item in maxima],
                  "wx", markersize=8)
        axis.set_xlim(0.0, 360.0)

        if _is_geometric(radial):
            _plain_log_axis(axis, "y")

        axis.set_xlabel("Propagation azimuth (degrees)")
        axis.set_ylabel(_RADIAL_LABELS[spectrum.quantity])

    figure.colorbar(image, ax=axis, label=label)
    axis.set_title(
        f"{component.replace('_', ' ')}, {spectrum.frequency:.4g} Hz, "
        f"block {block} ({spectrum.times[block]})"
    )

    return figure, axis


def _spectrum_values(
    spectrum: FrequencySpectrum,
    block: int,
    quantity: str,
    normalize: bool,
    units: Optional[str],
) -> tuple:
    """Values and colour-bar label of a spectrum."""
    if quantity == "ellipticity":
        if spectrum.ellipticity is None:
            raise ValueError("the spectrum has no ellipticity")

        return spectrum.ellipticity[block], "Signed ellipticity"

    if quantity == "standard_error":
        if spectrum.power_standard_error is None:
            raise ValueError(
                "the spectrum has no standard error (empirical "
                "statistical uncertainty)"
            )

        values = spectrum.power_standard_error[block]
        name = "Standard error of the power"
    else:
        values = spectrum.power[block]
        name = "Power"

    if normalize:
        scale = np.nanmax(spectrum.power[block])
        return values / scale, f"{name} / maximum power"

    return values, name + ("" if units is None else f" ({units})")


def _radial_value(maximum: Any, quantity: GridQuantity) -> float:
    """Radial grid coordinate of a maximum."""
    if quantity is GridQuantity.VELOCITY:
        return maximum.velocity

    if quantity is GridQuantity.SLOWNESS:
        return maximum.slowness

    return maximum.wavenumber


def plot_density(
    density: MaximumDensity,
    peaks: Optional[DensityPeaks] = None,
    summary: Optional[DensitySummary] = None,
    scale: str = "linear",
    floor: float = 0.01,
    log_frequency: bool = True,
    ax: Any = None,
) -> tuple:
    """
    Plot the density of the spectral maxima versus frequency.

    The density of each frequency is divided by its maximum, so that
    weak modes stay visible. Logarithmic coordinates are shown as phase
    velocity on a logarithmic axis.

    Parameters
    ----------
    density : MaximumDensity
        Density (:func:`~shakelab.exploration.yarray.statistics.
        maximum_density`).
    peaks : DensityPeaks, optional
        Peaks to mark (circles sized by their relative height).
    summary : DensitySummary, optional
        Mode (line) and quantile interval (dashed lines) to draw.
    scale : {"linear", "log"}, optional
        Colour scale of the relative density. Default linear.
    floor : float, optional
        Lower limit of the logarithmic colour scale. Default 0.01.
    log_frequency : bool, optional
        Logarithmic frequency axis. Default True.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.

    Raises
    ------
    TypeError
        If ``density`` is not a MaximumDensity.
    ValueError
        If ``scale`` is unknown.
    """
    if not isinstance(density, MaximumDensity):
        raise TypeError("density must be a MaximumDensity")

    if scale not in ("linear", "log"):
        raise ValueError("scale must be 'linear' or 'log'")

    from matplotlib.colors import LogNorm, Normalize

    logarithmic = density.coordinate is not DensityCoordinate.WAVENUMBER

    def axis_values(values: np.ndarray) -> np.ndarray:
        """Plotted coordinate of density values."""
        if density.coordinate is DensityCoordinate.LOG_VELOCITY:
            return np.exp(values)

        if density.coordinate is DensityCoordinate.LOG_SLOWNESS:
            return np.exp(-values)

        return values

    figure, axis = _axes(ax, (9.0, 5.5))
    relative = density.relative_density(floor=floor if scale == "log"
                                        else 0.0)
    norm = (LogNorm(vmin=floor, vmax=1.0) if scale == "log"
            else Normalize(vmin=0.0, vmax=1.0))
    ordinate = axis_values(density.values)
    order = np.argsort(ordinate)
    image = axis.pcolormesh(
        density.frequencies, ordinate[order], relative[:, order].T,
        shading="nearest", norm=norm, cmap="magma_r",
    )
    figure.colorbar(image, ax=axis, label="Density / maximum")

    if peaks is not None:
        for column in range(peaks.values.shape[1]):
            valid = np.isfinite(peaks.values[:, column])
            axis.scatter(
                peaks.frequencies[valid],
                axis_values(peaks.values[valid, column]),
                s=15.0 + 45.0 * peaks.relative_density[valid, column],
                facecolors="none", edgecolors="tab:cyan",
                linewidths=1.2, label="density peaks" if column == 0
                else None,
            )

    if summary is not None:
        axis.plot(summary.frequencies, axis_values(summary.mode),
                  color="tab:green", label="mode")

        for values in (summary.lower, summary.upper):
            axis.plot(summary.frequencies, axis_values(values),
                      color="tab:green", linestyle="--", linewidth=0.8)

    if log_frequency:
        _plain_log_axis(axis, "x")

    if logarithmic:
        _plain_log_axis(axis, "y")
        axis.set_ylabel("Phase velocity (m/s)")
    else:
        axis.set_ylabel("Angular wavenumber (rad/m)")

    if peaks is not None or summary is not None:
        axis.legend(loc="best")

    axis.set_xlabel("Frequency (Hz)")
    axis.set_title(
        f"Density of maxima, {density.component.replace('_', ' ')} "
        f"({int(np.sum(density.counts))} maxima)"
    )

    return figure, axis


# ----------------------------------------------------------------------------
# Array response
# ----------------------------------------------------------------------------


def _response_values(power: np.ndarray, scale: str,
                     floor: float) -> tuple[np.ndarray, str]:
    """Return the displayed response and its label."""
    if scale == "linear":
        return power, "Relative power"

    if scale != "db":
        raise ValueError("scale must be 'db' or 'linear'")

    floor = float(floor)

    if not np.isfinite(floor) or floor >= 0.0:
        raise ValueError("floor must be negative (dB)")

    return (10.0 * np.log10(np.maximum(power, 10.0 ** (floor / 10.0))),
            "Relative power (dB)")


def _response_image(axis: Any, figure: Any, response: ArrayResponse,
                    power: np.ndarray, scale: str, floor: float) -> None:
    """Draw a response as an image with its colour bar."""
    display, label = _response_values(power, scale, floor)
    image = axis.pcolormesh(response.kx, response.ky, display,
                            shading="nearest", cmap="viridis")
    figure.colorbar(image, ax=axis, label=label)
    axis.set_aspect("equal", adjustable="box")


def plot_array_geometry(
    array: Any,
    principal_axes: bool = True,
    centered: bool = False,
    labels: Any = None,
    ax: Any = None,
) -> tuple:
    """
    Plot the horizontal layout of the stations.

    Parameters
    ----------
    array : ArrayData, ArrayGeometry or ArrayResponse
        Array, geometry or response.
    principal_axes : bool, optional
        Draw the principal directions of the geometry, two standard
        deviations of the coordinates long on each side. Default True.
    centered : bool, optional
        Coordinates relative to the centroid. Default False.
    labels : bool or sequence of str, optional
        Station labels; True for the identifiers of an ArrayData.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    names = None

    if isinstance(array, ArrayData):
        names = array.sids
        geometry = ArrayGeometry.from_array(array)
    elif isinstance(array, ArrayResponse):
        geometry = array.geometry
    elif isinstance(array, ArrayGeometry):
        geometry = array
    else:
        raise TypeError("array must be an ArrayData, ArrayGeometry or "
                        "ArrayResponse")

    if labels is True:
        if names is None:
            raise ValueError("labels=True needs an ArrayData")

        labels = names

    coordinates = (geometry.centered_coordinates if centered
                   else geometry.coordinates)
    centroid = np.zeros(2) if centered else geometry.centroid
    figure, axis = _axes(ax, (7.0, 7.0))
    axis.scatter(coordinates[:, 0], coordinates[:, 1], marker="^",
                 color="tab:blue", label="stations", zorder=3)
    axis.plot(centroid[0], centroid[1], "+", color="k", markersize=10,
              label="centroid")

    if labels:
        for (x, y), name in zip(coordinates, labels):
            axis.annotate(str(name), (x, y), textcoords="offset points",
                          xytext=(4, 4), fontsize=8)

    if principal_axes:
        for index in range(geometry.rank):
            half = (2.0 * np.sqrt(geometry.principal_values[index])
                    * geometry.principal_directions[:, index])
            axis.plot([centroid[0] - half[0], centroid[0] + half[0]],
                      [centroid[1] - half[1], centroid[1] + half[1]],
                      color="gray", linestyle="--", linewidth=1.0)

    axis.set_aspect("equal", adjustable="datalim")
    axis.set_xlabel("East")
    axis.set_ylabel("North")
    axis.set_title(f"Array geometry ({geometry.n_stations} stations)")
    axis.grid(True, alpha=0.25)
    axis.legend(loc="best")

    return figure, axis


def plot_array_response(
    response: Any,
    scale: str = "db",
    floor: float = -40.0,
    expected: bool = False,
    peaks: Optional[int] = 12,
    ax: Any = None,
) -> tuple:
    """
    Plot the geometric response in wavenumber space.

    With a ResponseAnalysis the main-lobe contour (solid), the
    ambiguity-level contour (dashed) and the strongest remote maxima
    (crosses) are drawn as well.

    Parameters
    ----------
    response : ArrayResponse or ResponseAnalysis
        Response, or its analysis.
    scale : {"db", "linear"}, optional
        Power scale. Default decibels.
    floor : float, optional
        Lowest level shown in decibels. Default -40.
    expected : bool, optional
        Show the response expected under the coordinate
        uncertainties. Default: as analysed.
    peaks : int, optional
        Remote maxima marked (None for all). Default 12.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    analysis = None

    if isinstance(response, ResponseAnalysis):
        analysis = response
        expected = expected or analysis.expected
        response = analysis.response

    if not isinstance(response, ArrayResponse):
        raise TypeError("response must be an ArrayResponse or a "
                        "ResponseAnalysis")

    power = response.selected_power(expected)
    figure, axis = _axes(ax, (7.5, 7.0))
    _response_image(axis, figure, response, power, scale, floor)

    if analysis is not None:
        axis.contour(response.kx, response.ky, power,
                     levels=[analysis.main_lobe.level], colors="w",
                     linewidths=1.2)
        axis.contour(response.kx, response.ky, power,
                     levels=[analysis.ambiguity.level], colors="w",
                     linewidths=0.8, linestyles="--")
        marked = analysis.ambiguity.peaks

        if peaks is not None:
            marked = marked[:max(int(peaks), 0)]

        for peak in marked:
            axis.plot(peak.kx, peak.ky, "x", color="tab:red",
                      markersize=7)

    axis.axhline(0.0, color="w", linewidth=0.5, alpha=0.5)
    axis.axvline(0.0, color="w", linewidth=0.5, alpha=0.5)
    axis.set_xlabel("$k_x$ (rad per unit)")
    axis.set_ylabel("$k_y$ (rad per unit)")
    axis.set_title("Expected array response" if expected
                   else "Array response")

    return figure, axis


def plot_resolution(
    lobe: Any,
    quantity: str = "radius",
    polar: bool = False,
    ax: Any = None,
) -> tuple:
    """
    Plot the directional size of the main lobe.

    Parameters
    ----------
    lobe : MainLobe or ResponseAnalysis
        Main lobe, or the analysis holding it.
    quantity : {"radius", "separation"}, optional
        One-sided radius (default) or full width.
    polar : bool, optional
        Polar diagram over 360 degrees instead of a curve over the
        axis azimuths ``[0, 180)``. Default False.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on (polar if ``polar``); a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if isinstance(lobe, ResponseAnalysis):
        lobe = lobe.main_lobe

    if not isinstance(lobe, MainLobe):
        raise TypeError("lobe must be a MainLobe or a ResponseAnalysis")

    if quantity == "radius":
        values, label = lobe.radii, "Main-lobe radius"
    elif quantity == "separation":
        values, label = lobe.resolution_separation, "Resolution separation"
    else:
        raise ValueError("quantity must be 'radius' or 'separation'")

    if polar:
        if ax is None:
            figure, axis = _pyplot().subplots(
                figsize=(7.0, 7.0), subplot_kw={"projection": "polar"}
            )
        elif getattr(ax, "name", None) != "polar":
            raise ValueError("polar=True needs a polar axis")
        else:
            figure, axis = ax.figure, ax

        angles = np.radians(np.concatenate(
            (lobe.azimuths, lobe.azimuths + 180.0, [360.0])
        ))
        axis.plot(angles, np.concatenate((values, values, values[:1])))
        axis.set_theta_zero_location("N")
        axis.set_theta_direction(-1)
    else:
        figure, axis = _axes(ax, (9.0, 5.5))
        axis.plot(np.append(lobe.azimuths, 180.0),
                  np.append(values, values[0]))
        axis.set_xlim(0.0, 180.0)
        axis.set_xlabel("Axis azimuth (degrees)")
        axis.set_ylabel(f"{label} (rad per unit)")
        axis.grid(True, alpha=0.25)

    axis.set_title(f"{label} at level {lobe.level:g}")

    return figure, axis


def plot_search_domain(
    analysis: ResponseAnalysis,
    domain: SearchDomainAnalysis,
    scale: str = "db",
    floor: float = -40.0,
    ax: Any = None,
) -> tuple:
    """
    Plot the difference disk of a search domain over the response.

    The circle has radius ``2 K``: the differences of the wavevectors
    of a full-azimuth search of radius ``K``.

    Parameters
    ----------
    analysis : ResponseAnalysis
        Response analysis.
    domain : SearchDomainAnalysis
        Search domain (:func:`~shakelab.exploration.yarray.response.
        analyze_search_domain`).
    scale : {"db", "linear"}, optional
        Power scale. Default decibels.
    floor : float, optional
        Lowest level shown in decibels. Default -40.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if not isinstance(analysis, ResponseAnalysis):
        raise TypeError("analysis must be a ResponseAnalysis")

    if not isinstance(domain, SearchDomainAnalysis):
        raise TypeError("domain must be a SearchDomainAnalysis")

    response = analysis.response
    power = response.selected_power(analysis.expected)
    figure, axis = _axes(ax, (7.5, 7.0))
    _response_image(axis, figure, response, power, scale, floor)
    axis.contour(response.kx, response.ky, power,
                 levels=[analysis.main_lobe.level], colors="w",
                 linewidths=1.0)
    axis.contour(response.kx, response.ky, power,
                 levels=[domain.ambiguity_level], colors="w",
                 linewidths=0.8, linestyles="--")
    angles = np.linspace(0.0, 2.0 * np.pi, 361)
    radius = domain.difference_radius
    axis.plot(radius * np.sin(angles), radius * np.cos(angles),
              color="tab:orange", linewidth=1.5,
              label=f"|dk| = 2K = {radius:.4g}")

    for peak in domain.ambiguity_peaks:
        axis.plot(peak.kx, peak.ky, "x", color="tab:red", markersize=8)

    axis.set_xlabel("$\\Delta k_x$ (rad per unit)")
    axis.set_ylabel("$\\Delta k_y$ (rad per unit)")
    axis.set_title("Search domain: "
                   + ("safe" if domain.is_safe else "ambiguous"))
    axis.legend(loc="upper right")

    return figure, axis


def plot_velocity_limits(
    limits: VelocityLimits,
    fill: bool = True,
    log_frequency: bool = True,
    log_velocity: bool = True,
    ax: Any = None,
) -> tuple:
    """
    Plot the geometric phase-velocity limits versus frequency.

    Parameters
    ----------
    limits : VelocityLimits
        Limits (:func:`~shakelab.exploration.yarray.response.
        velocity_limits`).
    fill : bool, optional
        Shade the usable band. Default True.
    log_frequency, log_velocity : bool, optional
        Logarithmic axes. Default True.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if not isinstance(limits, VelocityLimits):
        raise TypeError("limits must be a VelocityLimits")

    figure, axis = _axes(ax, (9.0, 5.5))
    axis.plot(limits.frequencies, limits.minimum, color="tab:red",
              label=f"aliasing (level {limits.ambiguity_level:g})")
    axis.plot(limits.frequencies, limits.maximum, color="tab:blue",
              label=f"resolution (level {limits.resolution_level:g})")

    if fill:
        axis.fill_between(limits.frequencies, limits.minimum,
                          limits.maximum, color="tab:green", alpha=0.15)

    if log_frequency:
        _plain_log_axis(axis, "x")

    if log_velocity:
        _plain_log_axis(axis, "y")

    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel("Phase velocity")
    axis.set_title("Geometric velocity limits")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(loc="best")

    return figure, axis


# ----------------------------------------------------------------------------
# Performance and operational assessment
# ----------------------------------------------------------------------------


_OUTCOME_COLORS = {
    OutcomeClass.LOCAL: "tab:green",
    OutcomeClass.ALIAS: "tab:red",
    OutcomeClass.OTHER: "tab:orange",
    OutcomeClass.BOUNDARY: "tab:gray",
}


def plot_performance(result: PerformanceResult, ax: Any = None) -> tuple:
    """
    Plot the Monte Carlo estimates in wavenumber space.

    The estimates are coloured by outcome; the true wavevector is a
    black star and the main lobe around it a dashed contour.

    Parameters
    ----------
    result : PerformanceResult
        Experiment (:func:`~shakelab.exploration.yarray.performance.
        evaluate_performance`).
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if not isinstance(result, PerformanceResult):
        raise TypeError("result must be a PerformanceResult")

    figure, axis = _axes(ax, (7.0, 7.0))
    samples = result.samples
    source = result.source

    for outcome, color in _OUTCOME_COLORS.items():
        mask = samples.mask(outcome)

        if np.any(mask):
            axis.scatter(
                samples.kx[mask], samples.ky[mask], s=12, color=color,
                alpha=0.6, label=f"{outcome.name.lower()} "
                f"({result.outcome.probability(outcome):.3f})",
            )

    response = result.response
    axis.contour(response.kx + source.kx, response.ky + source.ky,
                 response.power, levels=[result.analysis.main_lobe.level],
                 colors="k", linewidths=0.8, linestyles="--")
    axis.plot(source.kx, source.ky, "*", color="k", markersize=12,
              label="true")
    limit = result.config.grid.maximum
    angles = np.linspace(0.0, 2.0 * np.pi, 361)
    axis.plot(limit * np.sin(angles), limit * np.cos(angles), color="gray",
              linewidth=0.8)
    axis.set_xlim(-1.05 * limit, 1.05 * limit)
    axis.set_ylim(-1.05 * limit, 1.05 * limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("$k_x$ (rad per unit)")
    axis.set_ylabel("$k_y$ (rad per unit)")
    config = result.config
    axis.set_title(
        f"{config.estimator.method.value}, SNR {config.snr_db:g} dB, "
        f"{config.monte_carlo.realizations} realizations"
    )
    axis.legend(loc="upper right", fontsize=8)

    return figure, axis


_SWEEP_LABELS = {
    "local_probability": "LOCAL probability",
    "alias_probability": "ALIAS probability",
    "other_probability": "OTHER probability",
    "boundary_probability": "BOUNDARY probability",
    "local_rmse": "LOCAL wavevector RMSE",
    "radial_standard_deviation": "Radial standard deviation",
    "transverse_standard_deviation": "Transverse standard deviation",
}


def plot_performance_sweep(
    sweep: PerformanceSweepResult,
    quantity: str = "local_probability",
    ax: Any = None,
) -> tuple:
    """
    Plot a quantity of a sweep versus SNR and realizations.

    Parameters
    ----------
    sweep : PerformanceSweepResult
        Sweep (:func:`~shakelab.exploration.yarray.performance.
        performance_sweep`).
    quantity : str, optional
        ``"local_probability"`` (default), ``"alias_probability"``,
        ``"other_probability"``, ``"boundary_probability"``,
        ``"local_rmse"``, ``"radial_standard_deviation"`` or
        ``"transverse_standard_deviation"``.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if not isinstance(sweep, PerformanceSweepResult):
        raise TypeError("sweep must be a PerformanceSweepResult")

    if quantity not in _SWEEP_LABELS:
        raise ValueError(f"unknown quantity {quantity!r}")

    values = getattr(sweep, quantity)
    probability = quantity.endswith("probability")
    figure, axis = _axes(ax, (8.0, 5.5))
    image = axis.pcolormesh(
        np.arange(sweep.shape[1]), sweep.snr_db, values,
        shading="nearest", cmap="viridis",
        vmin=0.0 if probability else None, vmax=1.0 if probability else None,
    )
    figure.colorbar(image, ax=axis, label=_SWEEP_LABELS[quantity])

    for (i, j), value in np.ndenumerate(values):
        if np.isfinite(value):
            light = image.norm(value) > 0.6
            axis.text(j, sweep.snr_db[i], f"{value:.3g}", ha="center",
                      va="center", fontsize=8, color="k" if light else "w")

    axis.set_xticks(np.arange(sweep.shape[1]))
    axis.set_xticklabels([str(value) for value in sweep.realizations])
    axis.set_xlabel("Fourier realizations")
    axis.set_ylabel("SNR (dB)")
    axis.set_title(_SWEEP_LABELS[quantity])

    return figure, axis


_OPERATIONAL_LABELS = {
    "local_probability": "LOCAL probability (worst azimuth)",
    "alias_probability": "ALIAS probability (worst azimuth)",
    "relative_velocity_rmse": "Relative velocity RMSE (worst azimuth)",
}


def plot_operational(
    assessment: OperationalAssessment,
    quantity: str = "local_probability",
    log_frequency: bool = True,
    log_velocity: bool = True,
    ax: Any = None,
) -> tuple:
    """
    Plot the frequency-velocity envelope of an array.

    The colour is a worst-azimuth quantity; the usable points are
    circled.

    Parameters
    ----------
    assessment : OperationalAssessment
        Assessment (:func:`~shakelab.exploration.yarray.operational.
        assess_array`).
    quantity : str, optional
        ``"local_probability"`` (default), ``"alias_probability"`` or
        ``"relative_velocity_rmse"``.
    log_frequency, log_velocity : bool, optional
        Logarithmic axes. Default True.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on; a new figure otherwise.

    Returns
    -------
    tuple
        ``(figure, axis)``.
    """
    if not isinstance(assessment, OperationalAssessment):
        raise TypeError("assessment must be an OperationalAssessment")

    if quantity not in _OPERATIONAL_LABELS:
        raise ValueError(f"unknown quantity {quantity!r}")

    values = getattr(assessment, f"conservative_{quantity}")
    probability = quantity.endswith("probability")
    cmap = "viridis" if quantity == "local_probability" else "viridis_r"
    figure, axis = _axes(ax, (9.0, 5.5))
    image = axis.pcolormesh(
        assessment.frequencies, assessment.velocities, values.T,
        shading="nearest", cmap=cmap,
        vmin=0.0 if probability else None, vmax=1.0 if probability else None,
    )
    figure.colorbar(image, ax=axis, label=_OPERATIONAL_LABELS[quantity])
    usable = assessment.conservative_usable
    f, c = np.meshgrid(assessment.frequencies, assessment.velocities,
                       indexing="ij")

    if np.any(usable):
        axis.scatter(f[usable], c[usable], s=60, facecolors="none",
                     edgecolors="w", linewidths=1.5, label="usable")
        axis.legend(loc="best")

    if log_frequency:
        _plain_log_axis(axis, "x")

    if log_velocity:
        _plain_log_axis(axis, "y")

    template = assessment.performance_template
    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel("Phase velocity")
    axis.set_title(
        f"Operational envelope: {template.estimator.method.value}, SNR "
        f"{template.snr_db:g} dB, {template.monte_carlo.realizations} "
        "realizations"
    )

    return figure, axis
