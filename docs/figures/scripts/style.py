# ****************************************************************************
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
# ****************************************************************************
"""
Common style of the figures of the documentation.

The figures are PNG files on an opaque white background, 1600 pixels
wide for a full-width figure (8 inches at 200 dots per inch). Series
take the colors of a fixed categorical order, validated for colour
vision deficiencies; text, axes and grids are in neutral greys, and a
reference series (data of another source) is drawn in grey.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

__all__ = [
    "SERIES",
    "REFERENCE",
    "TEXT",
    "MUTED",
    "GRID",
    "WIDTH",
    "DPI",
    "REPO",
    "TEST_DATA",
    "FIGURES",
    "apply",
    "figure",
    "save",
    "label_panel",
]

# Categorical order: blue, orange, aqua (validated as a set of three)
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")
REFERENCE = "#8c8b86"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e4e3df"

WIDTH = 8.0
DPI = 200

REPO = Path(__file__).resolve().parents[3]
TEST_DATA = REPO / "tests" / "signals" / "data"
FIGURES = REPO / "docs" / "figures"


def apply() -> None:
    """Set the matplotlib parameters of the documentation figures."""
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "axes.labelcolor": MUTED,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 0.9,
        "text.color": TEXT,
    })


def figure(nrows: int = 1, ncols: int = 1, height: float = 3.0,
           **options) -> tuple:
    """
    Create a figure of the documentation width.

    Parameters
    ----------
    nrows : int, optional
        Number of rows of panels.
    ncols : int, optional
        Number of columns of panels.
    height : float, optional
        Height of the figure, in inches.
    **options
        Further arguments of ``matplotlib.pyplot.subplots``.

    Returns
    -------
    tuple
        The figure and its axes.
    """
    apply()
    return plt.subplots(nrows, ncols, figsize=(WIDTH, height),
                        layout="constrained", **options)


def save(fig: Figure, name: str) -> Path:
    """
    Save a figure under ``docs/figures`` and close it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    name : str
        Path relative to ``docs/figures``, e.g.
        ``"signals/waveforms_filters.png"``.

    Returns
    -------
    pathlib.Path
        Path of the saved file.
    """
    path = FIGURES / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def label_panel(ax: Axes, text: str) -> None:
    """
    Write a short label in the upper-left corner of a panel.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Panel.
    text : str
        Label.
    """
    ax.text(0.005, 0.97, text, transform=ax.transAxes, ha="left",
            va="top", fontsize=8, color=TEXT,
            bbox={"facecolor": "white", "edgecolor": "none",
                  "pad": 1.5})
