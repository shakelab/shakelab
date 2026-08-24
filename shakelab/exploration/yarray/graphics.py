# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
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
# You should have received a copy of the GNU Affero General Public
# License with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Graphics utilities for seismic-array processing results.

The functions in this module operate exclusively on
:class:`BeamformingResult` objects. They do not read waveform files,
perform beamforming, or modify processing results.

The plotting layer therefore remains independent from data ingestion
and numerical processing and can be reused by scripts, notebooks, and
future graphical interfaces.

Propagation azimuths are expressed clockwise from North.
"""

from __future__ import annotations

from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from .beamforming import (
    BeamformingResult,
    FrequencySpectrum,
    SpectralMaximum,
)
from .response import (
    AmbiguityAnalysis,
    ArrayGeometry,
    ArrayResponse,
    MainLobe,
    SearchDomainAnalysis,
    VelocityLimits,
)


FloatArray = NDArray[np.float64]


def _select_maxima(
    result: BeamformingResult,
    *,
    component: str | None = None,
    ranks: Iterable[int] | None = None,
) -> tuple[SpectralMaximum, ...]:
    """Return maxima matching component and rank filters."""
    if not isinstance(result, BeamformingResult):
        raise TypeError(
            "result must be a BeamformingResult instance."
        )

    selected = result.maxima

    if component is not None:
        selected = tuple(
            maximum
            for maximum in selected
            if maximum.component == component
        )

    if ranks is not None:
        rank_set = {
            int(rank)
            for rank in ranks
        }

        if not rank_set:
            raise ValueError(
                "ranks cannot be empty."
            )

        if any(rank < 1 for rank in rank_set):
            raise ValueError(
                "ranks must contain positive integers."
            )

        selected = tuple(
            maximum
            for maximum in selected
            if maximum.rank in rank_set
        )

    return selected


def _available_components(
    result: BeamformingResult,
) -> tuple[str, ...]:
    """Return components present in maxima or retained spectra."""
    components = {
        maximum.component
        for maximum in result.maxima
    }

    if result.spectra is not None:
        components.update(
            spectrum.component
            for spectrum in result.spectra
        )

    return tuple(
        sorted(components)
    )


def _resolve_component(
    result: BeamformingResult,
    component: str | None,
) -> str:
    """Resolve an omitted component only when it is unambiguous."""
    if component is not None:
        return component

    components = _available_components(
        result
    )

    if len(components) == 1:
        return components[0]

    if not components:
        raise ValueError(
            "The result contains no component information."
        )

    values = ", ".join(
        components
    )

    raise ValueError(
        "component must be specified when multiple components are "
        f"available: {values}."
    )


def _resolve_spectrum(
    result: BeamformingResult,
    *,
    frequency: float,
    block: int,
    component: str | None,
    nearest: bool,
) -> FrequencySpectrum:
    """Return one retained spectrum matching the requested selection."""
    if result.spectra is None:
        raise ValueError(
            "The BeamformingResult does not contain retained spectra."
        )

    component = _resolve_component(
        result,
        component,
    )

    candidates = tuple(
        spectrum
        for spectrum in result.spectra
        if spectrum.component == component
    )

    if not candidates:
        raise ValueError(
            f"No retained spectra are available for {component!r}."
        )

    frequency = float(
        frequency
    )

    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(
            "frequency must be finite and positive."
        )

    if nearest:
        spectrum = min(
            candidates,
            key=lambda item: abs(
                item.frequency - frequency
            ),
        )
    else:
        matches = tuple(
            spectrum
            for spectrum in candidates
            if np.isclose(
                spectrum.frequency,
                frequency,
                rtol=1e-12,
                atol=0.0,
            )
        )

        if len(matches) != 1:
            raise ValueError(
                f"No unique spectrum exists at {frequency:g} Hz "
                f"for component {component!r}."
            )

        spectrum = matches[0]

    if block < 0 or block >= spectrum.power.shape[0]:
        raise IndexError(
            f"Block {block} is outside the valid range "
            f"0..{spectrum.power.shape[0] - 1} for "
            f"{spectrum.frequency:g} Hz."
        )

    return spectrum


def _spectrum_maxima(
    result: BeamformingResult,
    *,
    spectrum: FrequencySpectrum,
    block: int,
    ranks: Iterable[int] | None = None,
) -> tuple[SpectralMaximum, ...]:
    """Return ranked maxima belonging to one retained spectrum."""
    maxima = tuple(
        maximum
        for maximum in result.maxima
        if (
            maximum.component == spectrum.component
            and maximum.block == block
            and np.isclose(
                maximum.frequency,
                spectrum.frequency,
                rtol=1e-12,
                atol=0.0,
            )
        )
    )

    if ranks is not None:
        rank_set = {
            int(rank)
            for rank in ranks
        }

        maxima = tuple(
            maximum
            for maximum in maxima
            if maximum.rank in rank_set
        )

    return tuple(
        sorted(
            maxima,
            key=lambda maximum: maximum.rank,
        )
    )


def _new_axes(
    ax: Axes | None,
    *,
    figsize: tuple[float, float],
) -> tuple[Figure, Axes]:
    """Return an existing axis or create a new figure and axis."""
    if ax is None:
        figure, axis = plt.subplots(
            figsize=figsize
        )
        return figure, axis

    return ax.figure, ax


def _normalised_power(
    power: FloatArray,
) -> FloatArray:
    """Return beam power normalized by its finite positive maximum."""
    power = np.asarray(
        power,
        dtype=float,
    )

    finite = np.isfinite(
        power
    )

    if not np.any(finite):
        raise ValueError(
            "Spectrum contains no finite beam-power values."
        )

    maximum = float(
        np.max(
            power[finite]
        )
    )

    if maximum <= 0.0:
        raise ValueError(
            "Spectrum maximum must be positive."
        )

    return power / maximum


def _closed_azimuth(
    azimuths: FloatArray,
    values: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Close a periodic azimuth axis at 360 degrees."""
    return (
        np.concatenate(
            (
                np.asarray(
                    azimuths,
                    dtype=float,
                ),
                np.array(
                    [360.0],
                    dtype=float,
                ),
            )
        ),
        np.concatenate(
            (
                np.asarray(values),
                np.asarray(values)[:1],
            ),
            axis=0,
        ),
    )


def plot_dispersion(
    result: BeamformingResult,
    *,
    component: str | None = None,
    ranks: Iterable[int] | None = None,
    ax: Axes | None = None,
    log_frequency: bool = True,
    log_velocity: bool = True,
) -> tuple[Figure, Axes]:
    """
    Plot extracted spectral maxima in frequency-velocity space.

    Parameters
    ----------
    result
        Beamforming result containing extracted maxima.
    component
        Component to display. It may be omitted when the result
        contains only one component.
    ranks
        Optional subset of one-based maximum ranks.
    ax
        Existing Matplotlib axis. A new figure is created when omitted.
    log_frequency
        Use a logarithmic frequency axis.
    log_velocity
        Use a logarithmic phase-velocity axis.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    component = _resolve_component(
        result,
        component,
    )

    maxima = _select_maxima(
        result,
        component=component,
        ranks=ranks,
    )

    if not maxima:
        raise ValueError(
            "No spectral maxima match the requested selection."
        )

    figure, axis = _new_axes(
        ax,
        figsize=(9.0, 6.0),
    )

    available_ranks = sorted(
        {
            maximum.rank
            for maximum in maxima
        }
    )

    for rank in available_ranks:
        selected = tuple(
            maximum
            for maximum in maxima
            if maximum.rank == rank
        )

        axis.scatter(
            [
                maximum.frequency
                for maximum in selected
            ],
            [
                maximum.velocity
                for maximum in selected
            ],
            label=f"Rank {rank}",
        )

    if log_frequency:
        axis.set_xscale(
            "log"
        )

    if log_velocity:
        axis.set_yscale(
            "log"
        )

    axis.set_xlabel(
        "Frequency (Hz)"
    )
    axis.set_ylabel(
        "Phase velocity"
    )
    axis.set_title(
        f"{component.replace('_', ' ').title()} — dispersion"
    )
    axis.grid(
        True,
        alpha=0.25,
    )

    if len(available_ranks) > 1:
        axis.legend()

    return figure, axis


def plot_ellipticity(
    result: BeamformingResult,
    *,
    component: str = "rayleigh_joint",
    ranks: Iterable[int] | None = None,
    ax: Axes | None = None,
    log_frequency: bool = True,
) -> tuple[Figure, Axes]:
    """
    Plot signed Rayleigh ellipticity versus frequency.

    The search may be sampled internally in ellipticity angle, but the
    displayed quantity is always the signed ellipticity ratio.
    """
    maxima = _select_maxima(
        result,
        component=component,
        ranks=ranks,
    )

    maxima = tuple(
        maximum
        for maximum in maxima
        if maximum.ellipticity is not None
    )

    if not maxima:
        raise ValueError(
            "No maxima with ellipticity information are available."
        )

    figure, axis = _new_axes(
        ax,
        figsize=(9.0, 6.0),
    )

    available_ranks = sorted(
        {
            maximum.rank
            for maximum in maxima
        }
    )

    for rank in available_ranks:
        selected = tuple(
            maximum
            for maximum in maxima
            if maximum.rank == rank
        )

        axis.scatter(
            [
                maximum.frequency
                for maximum in selected
            ],
            [
                maximum.ellipticity
                for maximum in selected
            ],
            label=f"Rank {rank}",
        )

    axis.axhline(
        0.0,
        linewidth=0.8,
    )

    if log_frequency:
        axis.set_xscale(
            "log"
        )

    axis.set_xlabel(
        "Frequency (Hz)"
    )
    axis.set_ylabel(
        "Signed Rayleigh ellipticity"
    )
    axis.set_title(
        "Rayleigh ellipticity"
    )
    axis.grid(
        True,
        alpha=0.25,
    )

    if len(available_ranks) > 1:
        axis.legend()

    return figure, axis


def plot_spectrum(
    result: BeamformingResult,
    *,
    frequency: float,
    block: int = 0,
    component: str | None = None,
    view: str = "kxky",
    nearest: bool = True,
    ranks: Iterable[int] | None = None,
    normalise: bool = True,
    annotate: bool = True,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """
    Plot one retained frequency-wavenumber beam-power spectrum.

    Parameters
    ----------
    result
        Beamforming result with retained spectra.
    frequency
        Requested analysis frequency.
    block
        Covariance-block index.
    component
        Component name. It may be omitted when unambiguous.
    view
        ``"kxky"`` for Cartesian wavenumber space or ``"azimuth"``
        for azimuth-radial coordinates.
    nearest
        Select the nearest retained frequency when True. Otherwise an
        exact processed frequency is required.
    ranks
        Optional subset of extracted maxima to mark.
    normalise
        Normalize power by its maximum before plotting.
    annotate
        Label marked maxima with their rank. Ellipticity is also shown
        when available.
    ax
        Existing Matplotlib axis.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    spectrum = _resolve_spectrum(
        result,
        frequency=frequency,
        block=block,
        component=component,
        nearest=nearest,
    )

    maxima = _spectrum_maxima(
        result,
        spectrum=spectrum,
        block=block,
        ranks=ranks,
    )

    power = np.asarray(
        spectrum.power[block],
        dtype=float,
    )

    if normalise:
        power = _normalised_power(
            power
        )

    if view == "kxky":
        return _plot_spectrum_kxky(
            spectrum,
            power,
            maxima=maxima,
            block=block,
            normalised=normalise,
            annotate=annotate,
            ax=ax,
        )

    if view == "azimuth":
        return _plot_spectrum_azimuth(
            spectrum,
            power,
            maxima=maxima,
            block=block,
            normalised=normalise,
            annotate=annotate,
            ax=ax,
        )

    raise ValueError(
        "view must be either 'kxky' or 'azimuth'."
    )


def _plot_spectrum_kxky(
    spectrum: FrequencySpectrum,
    power: FloatArray,
    *,
    maxima: tuple[SpectralMaximum, ...],
    block: int,
    normalised: bool,
    annotate: bool,
    ax: Axes | None,
) -> tuple[Figure, Axes]:
    """Plot one spectrum in Cartesian angular-wavenumber space."""
    azimuths, closed_power = _closed_azimuth(
        spectrum.azimuths,
        power,
    )

    wavenumbers = spectrum.wavenumbers()
    angles = np.deg2rad(
        azimuths
    )

    kx = (
        np.sin(angles)[:, None]
        * wavenumbers[None, :]
    )
    ky = (
        np.cos(angles)[:, None]
        * wavenumbers[None, :]
    )

    figure, axis = _new_axes(
        ax,
        figsize=(7.0, 7.0),
    )

    contour = axis.contourf(
        kx,
        ky,
        closed_power,
        levels=31,
    )

    figure.colorbar(
        contour,
        ax=axis,
        label=(
            "Normalized beam power"
            if normalised
            else "Beam power"
        ),
    )

    for maximum in maxima:
        axis.plot(
            maximum.kx,
            maximum.ky,
            marker="x",
            markersize=8,
            markeredgewidth=1.5,
        )

        if annotate:
            axis.annotate(
                _maximum_label(
                    maximum
                ),
                (
                    maximum.kx,
                    maximum.ky,
                ),
                xytext=(5, 5),
                textcoords="offset points",
            )

    axis.axhline(
        0.0,
        linewidth=0.8,
    )
    axis.axvline(
        0.0,
        linewidth=0.8,
    )
    axis.set_aspect(
        "equal",
        adjustable="box",
    )
    axis.set_xlabel(
        r"$k_x$"
    )
    axis.set_ylabel(
        r"$k_y$"
    )
    axis.set_title(
        f"{spectrum.component.replace('_', ' ').title()} — "
        f"{spectrum.frequency:g} Hz, block {block}"
    )

    return figure, axis


def _plot_spectrum_azimuth(
    spectrum: FrequencySpectrum,
    power: FloatArray,
    *,
    maxima: tuple[SpectralMaximum, ...],
    block: int,
    normalised: bool,
    annotate: bool,
    ax: Axes | None,
) -> tuple[Figure, Axes]:
    """Plot one spectrum in azimuth-radial coordinates."""
    azimuths, closed_power = _closed_azimuth(
        spectrum.azimuths,
        power,
    )

    radial = np.asarray(
        spectrum.grid_values,
        dtype=float,
    )

    figure, axis = _new_axes(
        ax,
        figsize=(9.0, 6.0),
    )

    mesh = axis.pcolormesh(
        azimuths,
        radial,
        closed_power.T,
        shading="auto",
    )

    figure.colorbar(
        mesh,
        ax=axis,
        label=(
            "Normalized beam power"
            if normalised
            else "Beam power"
        ),
    )

    radial_maxima = _maximum_radial_values(
        spectrum,
        maxima,
    )

    for maximum, radial_value in zip(
        maxima,
        radial_maxima,
    ):
        axis.plot(
            maximum.azimuth,
            radial_value,
            marker="x",
            markersize=8,
            markeredgewidth=1.5,
        )

        if annotate:
            axis.annotate(
                _maximum_label(
                    maximum
                ),
                (
                    maximum.azimuth,
                    radial_value,
                ),
                xytext=(5, 5),
                textcoords="offset points",
            )

    axis.set_xlim(
        0.0,
        360.0,
    )
    axis.set_xlabel(
        "Propagation azimuth (degrees)"
    )
    axis.set_ylabel(
        _radial_label(
            spectrum
        )
    )
    axis.set_title(
        f"{spectrum.component.replace('_', ' ').title()} — "
        f"{spectrum.frequency:g} Hz, block {block}"
    )

    if _is_logarithmic_grid(
        radial
    ):
        axis.set_yscale(
            "log"
        )

    return figure, axis


def _maximum_radial_values(
    spectrum: FrequencySpectrum,
    maxima: tuple[SpectralMaximum, ...],
) -> tuple[float, ...]:
    """Convert canonical maximum wavenumbers to spectrum radial units."""
    values: list[float] = []

    for maximum in maxima:
        if spectrum.quantity.value == "wavenumber":
            value = maximum.wavenumber

        elif spectrum.quantity.value == "velocity":
            value = maximum.velocity

        elif spectrum.quantity.value == "slowness":
            value = maximum.slowness

        else:
            raise ValueError(
                f"Unsupported grid quantity: {spectrum.quantity!r}."
            )

        values.append(
            float(value)
        )

    return tuple(values)


def _radial_label(
    spectrum: FrequencySpectrum,
) -> str:
    """Return a display label for the retained radial grid."""
    labels = {
        "velocity": "Phase velocity",
        "slowness": "Slowness",
        "wavenumber": "Angular wavenumber",
    }

    try:
        return labels[
            spectrum.quantity.value
        ]
    except KeyError as error:
        raise ValueError(
            f"Unsupported grid quantity: {spectrum.quantity!r}."
        ) from error


def _is_logarithmic_grid(
    values: FloatArray,
) -> bool:
    """Return whether positive radial values are geometrically spaced."""
    values = np.asarray(
        values,
        dtype=float,
    )

    if values.size < 3 or np.any(values <= 0.0):
        return False

    ratios = values[1:] / values[:-1]

    return bool(
        np.allclose(
            ratios,
            ratios[0],
            rtol=1e-6,
            atol=0.0,
        )
    )


def _maximum_label(
    maximum: SpectralMaximum,
) -> str:
    """Return concise annotation text for one spectral maximum."""
    label = str(
        maximum.rank
    )

    if maximum.ellipticity is not None:
        label += (
            f"\ne={maximum.ellipticity:.2f}"
        )

    return label
# ---------------------------------------------------------------------------
# Geometrical array-response graphics
# ---------------------------------------------------------------------------


def plot_array_geometry(
    geometry: ArrayGeometry,
    *,
    principal_axes: bool = True,
    centered: bool = False,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """
    Plot the horizontal geometry of a seismic array.

    Parameters
    ----------
    geometry
        Geometry descriptor returned by the array-response layer.
    principal_axes
        Draw the principal geometry directions when True.
    centered
        Plot centroid-centered coordinates instead of original
        coordinates.
    ax
        Existing Matplotlib axis. A new figure is created when omitted.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    if not isinstance(
        geometry,
        ArrayGeometry,
    ):
        raise TypeError(
            "geometry must be an ArrayGeometry instance."
        )

    coordinates = (
        geometry.centered_coordinates
        if centered
        else geometry.coordinates
    )

    centroid = (
        np.zeros(2, dtype=float)
        if centered
        else geometry.centroid
    )

    figure, axis = _new_axes(
        ax,
        figsize=(7.0, 7.0),
    )

    axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        marker="o",
        label="Stations",
    )

    axis.plot(
        centroid[0],
        centroid[1],
        marker="+",
        markersize=10,
        markeredgewidth=1.5,
        linestyle="none",
        label="Centroid",
    )

    if principal_axes:
        _plot_geometry_principal_axes(
            axis,
            geometry,
            centroid,
        )

    axis.set_aspect(
        "equal",
        adjustable="box",
    )
    axis.set_xlabel(
        "East coordinate"
    )
    axis.set_ylabel(
        "North coordinate"
    )
    axis.set_title(
        "Array geometry"
    )
    axis.grid(
        True,
        alpha=0.25,
    )
    axis.legend()

    return figure, axis


def plot_array_response(
    response: ArrayResponse,
    *,
    main_lobe: MainLobe | None = None,
    ambiguity: AmbiguityAnalysis | None = None,
    scale: str = "db",
    db_floor: float = -40.0,
    use_expected_response: bool = False,
    show_main_lobe: bool = True,
    show_ambiguity_regions: bool = True,
    show_peaks: bool = True,
    peak_number: int | None = 12,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """
    Plot the two-dimensional geometric ambiguity response.

    Parameters
    ----------
    response
        Geometry-only array response.
    main_lobe
        Optional central-lobe analysis. When supplied, its threshold
        contour can be displayed.
    ambiguity
        Optional ambiguity analysis. Its threshold contour and detected
        peaks can be displayed.
    scale
        ``"linear"`` or ``"db"``. Decibel values are relative power.
    db_floor
        Lower display limit used with ``scale="db"``.
    use_expected_response
        Display the expected response under coordinate uncertainty.
    show_main_lobe
        Draw the main-lobe threshold contour when available.
    show_ambiguity_regions
        Draw the ambiguity-level contour when available.
    show_peaks
        Mark detected ambiguity peaks when available.
    peak_number
        Maximum number of ambiguity peaks to mark. ``None`` displays all.
    ax
        Existing Matplotlib axis. A new figure is created when omitted.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    if not isinstance(
        response,
        ArrayResponse,
    ):
        raise TypeError(
            "response must be an ArrayResponse instance."
        )

    if main_lobe is not None and not isinstance(
        main_lobe,
        MainLobe,
    ):
        raise TypeError(
            "main_lobe must be a MainLobe instance or None."
        )

    if ambiguity is not None and not isinstance(
        ambiguity,
        AmbiguityAnalysis,
    ):
        raise TypeError(
            "ambiguity must be an AmbiguityAnalysis instance or None."
        )

    power = response.selected_power(
        use_expected_response=use_expected_response
    )

    display, label = _response_display_values(
        power,
        scale=scale,
        db_floor=db_floor,
    )

    figure, axis = _new_axes(
        ax,
        figsize=(7.5, 7.0),
    )

    mesh = axis.pcolormesh(
        response.kx,
        response.ky,
        display,
        shading="auto",
    )

    figure.colorbar(
        mesh,
        ax=axis,
        label=label,
    )

    if show_main_lobe and main_lobe is not None:
        axis.contour(
            response.kx,
            response.ky,
            power,
            levels=[main_lobe.level],
            linewidths=1.2,
        )

    if show_ambiguity_regions and ambiguity is not None:
        if (
            main_lobe is None
            or not np.isclose(
                ambiguity.level,
                main_lobe.level,
            )
        ):
            axis.contour(
                response.kx,
                response.ky,
                power,
                levels=[ambiguity.level],
                linewidths=1.0,
                linestyles="--",
            )

    if show_peaks and ambiguity is not None:
        peaks = ambiguity.peaks

        if peak_number is not None:
            if (
                isinstance(peak_number, bool)
                or not isinstance(peak_number, int)
                or peak_number < 1
            ):
                raise ValueError(
                    "peak_number must be a positive integer or None."
                )

            peaks = peaks[:peak_number]

        for peak in peaks:
            axis.plot(
                peak.kx,
                peak.ky,
                marker="x",
                markersize=7,
                markeredgewidth=1.2,
                linestyle="none",
            )

    axis.plot(
        0.0,
        0.0,
        marker="+",
        markersize=9,
        markeredgewidth=1.4,
        linestyle="none",
    )

    axis.axhline(
        0.0,
        linewidth=0.7,
    )
    axis.axvline(
        0.0,
        linewidth=0.7,
    )
    axis.set_aspect(
        "equal",
        adjustable="box",
    )
    axis.set_xlabel(
        r"$k_x$"
    )
    axis.set_ylabel(
        r"$k_y$"
    )

    response_name = (
        "Expected geometric response"
        if use_expected_response
        else "Geometric array response"
    )

    axis.set_title(
        response_name
    )

    return figure, axis


def plot_resolution(
    main_lobe: MainLobe,
    *,
    quantity: str = "radius",
    view: str = "cartesian",
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """
    Plot directional main-lobe resolution.

    Parameters
    ----------
    main_lobe
        Main-lobe analysis returned by :func:`analyze_response`.
    quantity
        ``"radius"`` for one-sided threshold radius or ``"separation"``
        for the corresponding full directional resolution separation.
    view
        ``"cartesian"`` for azimuth versus resolution or ``"polar"`` for
        a directional polar representation.
    ax
        Existing Matplotlib axis. For ``view="polar"`` it must be polar.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    if not isinstance(
        main_lobe,
        MainLobe,
    ):
        raise TypeError(
            "main_lobe must be a MainLobe instance."
        )

    if quantity == "radius":
        values = main_lobe.radii
        ylabel = "Resolution radius"
    elif quantity == "separation":
        values = main_lobe.resolution_separation
        ylabel = "Resolution separation"
    else:
        raise ValueError(
            "quantity must be either 'radius' or 'separation'."
        )

    if view == "cartesian":
        figure, axis = _new_axes(
            ax,
            figsize=(9.0, 5.5),
        )

        azimuths, closed_values = _closed_axis_response(
            main_lobe.azimuths,
            values,
        )

        axis.plot(
            azimuths,
            closed_values,
        )
        axis.set_xlim(
            0.0,
            180.0,
        )
        axis.set_xlabel(
            "Array-axis azimuth (degrees)"
        )
        axis.set_ylabel(
            ylabel
        )
        axis.grid(
            True,
            alpha=0.25,
        )

    elif view == "polar":
        figure, axis = _new_polar_axes(
            ax
        )

        angles, closed_values = _closed_polar_response(
            main_lobe.azimuths,
            values,
        )

        axis.plot(
            angles,
            closed_values,
        )
        axis.set_theta_zero_location(
            "N"
        )
        axis.set_theta_direction(
            -1
        )

    else:
        raise ValueError(
            "view must be either 'cartesian' or 'polar'."
        )

    axis.set_title(
        f"Directional {ylabel.lower()} at "
        f"power level {main_lobe.level:g}"
    )

    return figure, axis


def plot_search_domain(
    response: ArrayResponse,
    main_lobe: MainLobe,
    domain: SearchDomainAnalysis,
    *,
    scale: str = "db",
    db_floor: float = -40.0,
    use_expected_response: bool = False,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """
    Plot a circular search-domain difference disk over the response.

    The plotted circle has radius ``2*K`` because it represents the
    pairwise-difference domain corresponding to a full-azimuth search disk
    of radius ``K``.
    """
    if not isinstance(
        response,
        ArrayResponse,
    ):
        raise TypeError(
            "response must be an ArrayResponse instance."
        )

    if not isinstance(
        main_lobe,
        MainLobe,
    ):
        raise TypeError(
            "main_lobe must be a MainLobe instance."
        )

    if not isinstance(
        domain,
        SearchDomainAnalysis,
    ):
        raise TypeError(
            "domain must be a SearchDomainAnalysis instance."
        )

    power = response.selected_power(
        use_expected_response=use_expected_response
    )

    display, label = _response_display_values(
        power,
        scale=scale,
        db_floor=db_floor,
    )

    figure, axis = _new_axes(
        ax,
        figsize=(7.5, 7.0),
    )

    mesh = axis.pcolormesh(
        response.kx,
        response.ky,
        display,
        shading="auto",
    )

    figure.colorbar(
        mesh,
        ax=axis,
        label=label,
    )

    axis.contour(
        response.kx,
        response.ky,
        power,
        levels=[main_lobe.level],
        linewidths=1.1,
    )

    axis.contour(
        response.kx,
        response.ky,
        power,
        levels=[domain.ambiguity_level],
        linewidths=1.0,
        linestyles="--",
    )

    angles = np.linspace(
        0.0,
        2.0 * np.pi,
        361,
    )

    radius = domain.difference_radius

    axis.plot(
        radius * np.sin(angles),
        radius * np.cos(angles),
        linewidth=1.4,
        label="Difference-domain boundary",
    )

    for peak in domain.ambiguity_peaks:
        axis.plot(
            peak.kx,
            peak.ky,
            marker="x",
            markersize=8,
            markeredgewidth=1.3,
            linestyle="none",
        )

    axis.set_aspect(
        "equal",
        adjustable="box",
    )
    axis.set_xlabel(
        r"$\Delta k_x$"
    )
    axis.set_ylabel(
        r"$\Delta k_y$"
    )
    axis.set_title(
        "Search-domain ambiguity — "
        + (
            "safe"
            if domain.is_safe
            else "ambiguous"
        )
    )
    axis.legend()

    return figure, axis


def plot_velocity_limits(
    limits: VelocityLimits,
    *,
    ax: Axes | None = None,
    log_frequency: bool = True,
    log_velocity: bool = True,
    fill: bool = True,
) -> tuple[Figure, Axes]:
    """
    Plot geometry-only phase-velocity bounds versus frequency.

    Parameters
    ----------
    limits
        Velocity limits returned by :func:`velocity_limits`.
    ax
        Existing Matplotlib axis. A new figure is created when omitted.
    log_frequency
        Use a logarithmic frequency axis.
    log_velocity
        Use a logarithmic phase-velocity axis.
    fill
        Shade the geometry-only usable velocity region.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axis)``.
    """
    if not isinstance(
        limits,
        VelocityLimits,
    ):
        raise TypeError(
            "limits must be a VelocityLimits instance."
        )

    figure, axis = _new_axes(
        ax,
        figsize=(9.0, 6.0),
    )

    axis.plot(
        limits.frequencies,
        limits.minimum,
        label="Minimum velocity",
    )
    axis.plot(
        limits.frequencies,
        limits.maximum,
        label="Maximum velocity",
    )

    if fill:
        axis.fill_between(
            limits.frequencies,
            limits.minimum,
            limits.maximum,
            alpha=0.15,
        )

    if log_frequency:
        axis.set_xscale(
            "log"
        )

    if log_velocity:
        axis.set_yscale(
            "log"
        )

    axis.set_xlabel(
        "Frequency (Hz)"
    )
    axis.set_ylabel(
        "Phase velocity"
    )
    axis.set_title(
        "Geometry-only frequency-velocity limits"
    )
    axis.grid(
        True,
        alpha=0.25,
    )
    axis.legend()

    return figure, axis


def _plot_geometry_principal_axes(
    axis: Axes,
    geometry: ArrayGeometry,
    centroid: FloatArray,
) -> None:
    """Draw principal array-geometry directions."""
    coordinates = geometry.centered_coordinates

    if coordinates.size == 0:
        return

    scale = 0.6 * geometry.maximum_aperture

    for index in range(geometry.rank):
        direction = geometry.principal_directions[:, index]
        half = 0.5 * scale * direction

        axis.plot(
            [
                centroid[0] - half[0],
                centroid[0] + half[0],
            ],
            [
                centroid[1] - half[1],
                centroid[1] + half[1],
            ],
            linewidth=1.0,
            linestyle="--",
        )


def _response_display_values(
    power: FloatArray,
    *,
    scale: str,
    db_floor: float,
) -> tuple[FloatArray, str]:
    """Return response values transformed for display."""
    power = np.asarray(
        power,
        dtype=float,
    )

    if scale == "linear":
        return (
            power,
            "Normalized response power",
        )

    if scale != "db":
        raise ValueError(
            "scale must be either 'linear' or 'db'."
        )

    db_floor = float(
        db_floor
    )

    if (
        not np.isfinite(db_floor)
        or db_floor >= 0.0
    ):
        raise ValueError(
            "db_floor must be finite and negative."
        )

    floor_power = 10.0 ** (
        db_floor / 10.0
    )

    display = 10.0 * np.log10(
        np.maximum(
            power,
            floor_power,
        )
    )

    return (
        display,
        "Relative response power (dB)",
    )


def _closed_axis_response(
    azimuths: FloatArray,
    values: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Close a 180-degree axis-response representation."""
    azimuths = np.asarray(
        azimuths,
        dtype=float,
    )

    values = np.asarray(
        values,
        dtype=float,
    )

    return (
        np.concatenate(
            (
                azimuths,
                np.array([180.0]),
            )
        ),
        np.concatenate(
            (
                values,
                values[:1],
            )
        ),
    )


def _closed_polar_response(
    azimuths: FloatArray,
    values: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Expand an axis response to a closed 360-degree polar curve."""
    azimuths = np.asarray(
        azimuths,
        dtype=float,
    )

    values = np.asarray(
        values,
        dtype=float,
    )

    full_azimuth = np.concatenate(
        (
            azimuths,
            azimuths + 180.0,
            np.array([360.0]),
        )
    )

    full_values = np.concatenate(
        (
            values,
            values,
            values[:1],
        )
    )

    return (
        np.deg2rad(
            full_azimuth
        ),
        full_values,
    )


def _new_polar_axes(
    ax: Axes | None,
) -> tuple[Figure, Axes]:
    """Return an existing polar axis or create a new one."""
    if ax is None:
        figure, axis = plt.subplots(
            figsize=(7.0, 7.0),
            subplot_kw={
                "projection": "polar",
            },
        )

        return figure, axis

    if getattr(
        ax,
        "name",
        None,
    ) != "polar":
        raise ValueError(
            "A polar Matplotlib axis is required for view='polar'."
        )

    return ax.figure, ax
