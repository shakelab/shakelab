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
Figures of the guide docs/signals/seedlink.md (schematic, no data).
"""

from __future__ import annotations

from pathlib import Path

from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import style

__all__ = ["FIGURES"]


def _box(ax, x: float, y: float, w: float, h: float, title: str,
         text: str, edge: str, face: str = "white") -> None:
    """A rounded box with a title and an optional second line."""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=0.12",
                                facecolor=face, edgecolor=edge, lw=1.4))
    if text:
        ax.text(x + w / 2, y + h * 0.62, title, ha="center",
                va="center", fontsize=9.5, fontweight="bold",
                color=style.TEXT)
        ax.text(x + w / 2, y + h * 0.3, text, ha="center", va="center",
                fontsize=8, color=style.MUTED)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color=style.TEXT)


def _arrow(ax, start: tuple, end: tuple, label: str = "",
           dashed: bool = False, rad: float = 0.0) -> None:
    """An arrow with an optional label above its middle."""
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=12,
        color=style.MUTED, lw=1.3, linestyle="--" if dashed else "-",
        connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2))
    if label:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.18,
                label, ha="center", va="bottom", fontsize=7.8,
                color=style.MUTED)


def callbacks() -> Path:
    """Flow of the records from the server to the destinations."""
    fig, ax = style.figure(1, 1, height=3.6)
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 7.2)
    ax.axis("off")
    blue, orange, aqua = style.SERIES

    # The receiving thread.
    ax.add_patch(FancyBboxPatch((3.05, 0.25), 8.05, 6.6,
                                boxstyle="round,pad=0,rounding_size=0.2",
                                facecolor="#f2f6fb", edgecolor="none"))
    ax.text(7.07, 6.62, "receiving thread  (run or start)", ha="center",
            va="center", fontsize=8.5, color=blue)

    w, h = 2.5, 1.3
    mid = 3.0
    _box(ax, 0.1, mid, 2.2, h, "SeedLink", "server", style.MUTED)
    _box(ax, 3.4, mid, w, h, "Client", "decode, select", blue)
    rows = {"archive": 5.0, "buffer": mid, "callbacks": 0.6}
    wide = 3.4
    _box(ax, 7.4, rows["archive"], wide, h, "Archive", "MiniSEED files",
         aqua)
    _box(ax, 7.4, rows["buffer"], wide, h, "Buffer", "StreamCollection",
         aqua)
    _box(ax, 7.4, rows["callbacks"], wide, h, "Your callbacks",
         "on_record, every, on_event", orange)
    _box(ax, 12.05, rows["archive"], 3.85, h, "Files", "read later with "
         "io.reader", style.MUTED)
    _box(ax, 12.05, rows["buffer"], 3.85, h, "Your other threads",
         "get_data(), status()", blue)
    _box(ax, 12.05, rows["callbacks"], 3.85, h, "Worker thread",
         "long processing", blue)

    centre = mid + h / 2
    _arrow(ax, (2.3, centre), (3.4, centre), "packets")
    for name, y in rows.items():
        _arrow(ax, (5.9, centre), (7.4, y + h / 2))
    _arrow(ax, (10.8, rows["archive"] + h / 2),
           (12.05, rows["archive"] + h / 2))
    _arrow(ax, (10.8, rows["buffer"] + h / 2),
           (12.05, rows["buffer"] + h / 2), "copy")
    _arrow(ax, (10.8, rows["callbacks"] + h / 2),
           (12.05, rows["callbacks"] + h / 2), "queue")
    _arrow(ax, (9.1, rows["callbacks"] + h), (9.1, rows["buffer"]),
           dashed=True)
    ax.text(9.25, 2.25, "client.buffer", ha="left", va="center",
            fontsize=7.5, color=style.MUTED)
    return style.save(fig, "signals/seedlink_callbacks.png")


FIGURES = {
    "callbacks": callbacks,
}
