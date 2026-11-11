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
  peaks or summary.

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
from shakelab.exploration.yarray.statistics import (
    DensityCoordinate,
    DensityPeaks,
    DensitySummary,
    MaximumDensity,
)

__all__ = [
    "plot_density",
    "plot_dispersion",
    "plot_ellipticity",
    "plot_spectrum",
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
