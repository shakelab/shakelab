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
Figures of the guide docs/signals/io.md.

Data: the raw MiniSEED file of IV.MODE (INGV) from tests/signals/data;
the other figures are schematic.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

import style
from shakelab.signals.io import reader
from shakelab.signals.libio.mseed import msrawread

__all__ = ["FIGURES"]

# Lighter shade of each series colour, for alternate data records.
LIGHT = ("#9cc0ec", "#f5b49b", "#8fdcbf")


def _seconds(start_ns: int, origin_ns: int) -> float:
    """Seconds from an origin, both in Unix nanoseconds."""
    return (start_ns - origin_ns) / 1e9


# ---------------------------------------------------------------------------
# Figures

def mseed_records() -> Path:
    """MiniSEED data records and the ShakeLab records they form."""
    path = style.TEST_DATA / "MODE.mseed"
    raw = msrawread(path)
    data = reader(path)
    channels = data.sid
    origin = min(r.start_ns for r in raw)

    fig, axes = style.figure(3, 1, height=5.6, sharex=True,
                             height_ratios=(1.0, 1.1, 1.3))

    # Vertical component, to show where the signal changes fast.
    vertical = data["IV.MODE..HNZ"][0]
    ax = axes[0]
    ax.plot(vertical.time_axis(), vertical.data / 1e3, color=style.REFERENCE,
            lw=0.5)
    ax.set_ylabel("counts (×1000)")
    style.label_panel(ax, "IV.MODE..HNZ, raw counts")

    # Samples held by each 512-byte data record.
    ax = axes[1]
    for k, sid in enumerate(channels):
        own = [r for r in raw if r.sid == sid]
        start = [_seconds(r.start_ns, origin) for r in own]
        count = [len(r) for r in own]
        ax.step(start, count, where="post", color=style.SERIES[k],
                label=sid)
    ax.set_ylabel("samples per\ndata record")
    ax.set_ylim(0, 850)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    style.label_panel(ax, "STEIM2, 512 bytes per data record")

    # Data records (alternate shades) and the ShakeLab records.
    ax = axes[2]
    rows = []
    for k, sid in enumerate(channels):
        y = 2 * (len(channels) - k) - 1
        own = [r for r in raw if r.sid == sid]
        for j, record in enumerate(own):
            begin = _seconds(record.start_ns, origin)
            colour = style.SERIES[k] if j % 2 else LIGHT[k]
            ax.add_patch(Rectangle((begin, y), record.duration
                                   + record.delta, 0.7, color=colour,
                                   lw=0))
        rows.append((y + 0.35, f"{sid}: {len(own)} data records"))
        whole = data[sid][0]
        ax.add_patch(Rectangle((0.0, y - 0.95), whole.duration, 0.7,
                               color=style.SERIES[k], lw=0))
        rows.append((y - 0.6, "1 Record"))
    ax.set_yticks([y for y, _ in rows], [text for _, text in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.2, 2 * len(channels) + 0.9)
    ax.grid(False)
    ax.set_xlim(0, vertical.duration)
    ax.set_xlabel("seconds from 2012-05-20T02:03:23")
    style.label_panel(ax, "reader: data records joined into Records")
    return style.save(fig, "signals/io_mseed_records.png")


def sds_days() -> Path:
    """Day files of an SDS archive and a request after midnight."""
    fig, ax = style.figure(1, 1, height=2.6)
    minute = 60.0
    # Records of one channel, about 4 minutes each (schematic).
    starts = np.arange(-17.0, 19.0, 4.0) * minute + 50.0
    for start in starts:
        day = 140 if start < 0 else 141
        colour = style.SERIES[0] if day == 140 else style.SERIES[1]
        ax.add_patch(Rectangle((start, 1.0), 3.9 * minute, 0.8,
                               facecolor=colour, lw=0))
    # Request window.
    first, last = 2.0 * minute, 12.0 * minute
    ax.axvspan(first, last, color=style.GRID, zorder=0)
    ax.annotate("", xy=(first, 2.35), xytext=(last, 2.35),
                arrowprops={"arrowstyle": "<->", "color": style.MUTED,
                            "lw": 0.8})
    ax.text((first + last) / 2, 2.45, "starttime … endtime", ha="center",
            va="bottom", fontsize=8, color=style.MUTED)
    ax.axvline(0.0, color=style.TEXT, lw=0.8)
    ax.text(-0.3 * minute, 0.55, "IV.MODE..HNZ.D.2012.140", ha="right",
            va="center", fontsize=8, color=style.SERIES[0])
    ax.text(0.3 * minute, 0.55, "IV.MODE..HNZ.D.2012.141", ha="left",
            va="center", fontsize=8, color=style.SERIES[1])
    ax.text(0.0, 3.0, "midnight", ha="center", va="bottom", fontsize=8,
            zorder=3, bbox={"facecolor": "white", "edgecolor": "none",
                            "pad": 1.0})
    # The record across midnight.
    across = starts[starts < 0][-1]
    ax.add_patch(Rectangle((across, 1.0), 3.9 * minute, 0.8,
                           facecolor="none", edgecolor=style.TEXT,
                           lw=1.2))
    ax.annotate("begins on day 140,\nends on day 141",
                xy=(across + 2.0 * minute, 1.8),
                xytext=(-12.5 * minute, 2.55), fontsize=8,
                color=style.TEXT, ha="center",
                arrowprops={"arrowstyle": "->", "color": style.MUTED,
                            "lw": 0.8})
    ticks = np.arange(-16, 20, 4) * minute
    ax.set_xticks(ticks, [f"{int(t // minute):+d}" for t in ticks])
    ax.set_xlim(-17 * minute, 19 * minute)
    ax.set_ylim(0.2, 3.4)
    ax.set_yticks([1.4], ["records"])
    ax.tick_params(axis="y", length=0)
    ax.grid(False)
    ax.set_xlabel("minutes from midnight between day 140 and day 141 "
                  "(colour: the file holding the record)")
    return style.save(fig, "signals/io_sds_days.png")


def _box(ax, x: float, y: float, w: float, h: float, title: str,
         text: str, colour: str) -> None:
    """A rounded box with a title and a few lines of text."""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=0.08",
                                facecolor="white", edgecolor=colour,
                                lw=1.1))
    ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="top",
            fontsize=8.5, fontweight="bold", color=style.TEXT)
    ax.text(x + w / 2, y + h - 0.42, text, ha="center", va="top",
            fontsize=7.5, color=style.MUTED, linespacing=1.35)


def _arrow(ax, start: tuple, end: tuple) -> None:
    """An arrow between two points."""
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                 mutation_scale=9, color=style.MUTED,
                                 lw=0.9))


def reader_flow() -> Path:
    """The steps of io.reader."""
    fig, ax = style.figure(1, 1, height=3.3)
    ax.set_xlim(0, 16)
    ax.set_ylim(0.3, 6.6)
    ax.axis("off")
    blue, orange, aqua = style.SERIES
    w, h = 3.4, 2.0

    top = 4.4
    _box(ax, 0.1, top, w, h, "1. Sources",
         "file, pattern, folder,\nlist of them, list_file\n"
         "(~ and $VAR expanded;\nhidden files ignored)", blue)
    _box(ax, 0.1, 0.5, w, h + 0.3, "or an archive",
         "layout=\"sds\" + root,\nselect, starttime, endtime:\n"
         "only the files of the\nchannels and days\nrequested", blue)
    _box(ax, 4.2, top, w, h, "2. Format of each file",
         "format= if given,\nelse the extension,\n"
         "else the first 1024 bytes\n(is_mseed, is_sac, ...);\n"
         "MiniSEED for archives", orange)
    _box(ax, 4.2, 0.5, w, h + 0.3, "unknown format",
         "file named: ValueError\n\nfound in a folder or\n"
         "pattern: skipped,\none IoWarning", orange)
    _box(ax, 8.3, top, w, h, "3. Format module",
         "libio: msread, sacread,\ndynaread, ...\n"
         "options checked against\nthe reader signature", aqua)
    _box(ax, 8.3, 0.5, w, h + 0.3, "4. StreamCollection",
         "one Stream per channel;\ncontiguous records joined\n"
         "(Stream.append)", aqua)
    _box(ax, 12.4, top, w - 0.1, h, "5. Selection",
         "select: channel patterns\nstarttime, endtime:\n"
         "records dropped or cut", blue)
    _box(ax, 12.4, 0.5, w - 0.1, h + 0.3, "6. Metadata, result",
         "metadata= attached\n(set_metadata);\nadded to stream_collection\n"
         "if given; sorted", blue)

    mid = top + h / 2
    _arrow(ax, (0.1 + w, mid), (4.2, mid))
    _arrow(ax, (0.1 + w, 0.5 + (h + 0.3) / 2), (4.2, mid - 0.3))
    _arrow(ax, (4.2 + w, mid), (8.3, mid))
    _arrow(ax, (4.2 + w / 2, top), (4.2 + w / 2, 0.5 + h + 0.3))
    _arrow(ax, (8.3 + w / 2, top), (8.3 + w / 2, 0.5 + h + 0.3))
    _arrow(ax, (8.3 + w, 0.5 + (h + 0.3) / 2), (12.4, mid - 0.3))
    _arrow(ax, (12.4 + (w - 0.1) / 2, top),
           (12.4 + (w - 0.1) / 2, 0.5 + h + 0.3))
    return style.save(fig, "signals/io_reader_flow.png")


FIGURES = {
    "mseed_records": mseed_records,
    "reader_flow": reader_flow,
    "sds_days": sds_days,
}
