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
Figures of the guide docs/signals/waveforms.md.

Data: IV.MODE, 20 May 2012 Emilia earthquake (raw counts and
StationXML of INGV) and its processed records in the ESM database,
from tests/signals/data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import signal

import style
from shakelab.signals.base import Record, Stream, StreamCollection
from shakelab.signals.io import metadata_reader, reader

__all__ = ["FIGURES"]

ESM = style.TEST_DATA / "ESM"
ESM_NAME = "IV.MODE..{cha}.D.IT-2012-0008.{kind}.MP.ASC"


def _raw(metadata: bool = False) -> StreamCollection:
    """Raw counts of IV.MODE, optionally with the StationXML."""
    xml = style.TEST_DATA / "MODE.xml" if metadata else None
    return reader(style.TEST_DATA / "MODE.mseed", metadata=xml)


def _esm(channel: str, kind: str) -> Record:
    """Processed ESM record (ACC, VEL or DIS) of a channel."""
    name = ESM_NAME.format(cha=channel, kind=kind)
    return reader(ESM / "ASCII" / name)[0][0]


# ---------------------------------------------------------------------------
# Figures

def components() -> Path:
    """The three components of a strong-motion station, raw counts."""
    data = _raw()
    fig, axes = style.figure(3, 1, height=3.6, sharex=True)
    for ax, stream in zip(axes, data):
        record = stream[0]
        ax.plot(record.time_axis(), record.data, color=style.SERIES[0])
        ax.set_ylabel("counts")
        style.label_panel(ax, record.sid)
    axes[-1].set_xlabel(f"time from {data[0][0].time} (s)")
    axes[-1].set_xlim(0, data[0][0].duration)
    return style.save(fig, "signals/waveforms_components.png")


def tolerance() -> Path:
    """Joining rule: offsets of 0.3 and 0.6 samples."""
    fig, axes = style.figure(2, 1, height=2.6, sharex=True)
    cases = [(0.3, "δ = 0.3: m = 0, contiguous (joined)"),
             (0.6, "δ = 0.6: m = 1, one-sample gap")]
    for ax, (delta, text) in zip(axes, cases):
        ka = np.arange(-5, 0)
        kb = delta + np.arange(0, 5)
        ax.axvspan(-0.5, 0.5, color=style.GRID, zorder=0)
        ax.vlines(np.arange(-5, 5), -1, 1, color=style.GRID, lw=0.6)
        ax.plot(ka, np.zeros_like(ka), "o", color=style.SERIES[0],
                ms=6, label="record A")
        ax.plot(kb, np.zeros_like(kb), "o", color=style.SERIES[1],
                ms=6, label="record B")
        ax.plot([0], [0], "o", mfc="white", mec=style.MUTED, ms=7,
                label="next sample of A")
        ax.annotate("", xy=(delta, 0.45), xytext=(0, 0.45),
                    arrowprops={"arrowstyle": "->", "color": style.TEXT,
                                "lw": 0.8})
        ax.text(delta / 2, 0.6, "δ", ha="center", fontsize=9)
        ax.set_ylim(-1, 1)
        ax.set_yticks([])
        ax.grid(False)
        style.label_panel(ax, text)
    axes[0].fill_between([], [], color=style.GRID,
                         label="joined as contiguous")
    axes[0].legend(loc="lower right", ncols=4)
    axes[-1].set_xlabel("time from the next sample of A (samples)")
    axes[-1].set_xlim(-5.5, 5)
    return style.save(fig, "signals/waveforms_tolerance.png")


def filters() -> Path:
    """Zero-phase and causal Butterworth low-pass filters."""
    acc = _esm("HNZ", "ACC")
    corner, order = 1.0, 4
    sos = signal.butter(order, corner, "lowpass", fs=acc.rate,
                        output="sos")
    freq, h = signal.sosfreqz(sos, worN=4096, fs=acc.rate)
    fig, (left, right) = style.figure(1, 2, height=3.0,
                                      width_ratios=(1, 1.6))
    left.loglog(freq[1:], np.abs(h[1:]), color=style.SERIES[1],
                label="causal, |H|")
    left.loglog(freq[1:], np.abs(h[1:]) ** 2, color=style.SERIES[0],
                label="zero-phase, |H|²")
    left.axvline(corner, color=style.MUTED, lw=0.6, ls=":")
    for level, text, align in ((2 ** -0.5, "−3 dB", "bottom"),
                               (0.5, "−6 dB", "top")):
        left.axhline(level, color=style.MUTED, lw=0.6, ls=":")
        left.text(18, level, text, ha="right", va=align, fontsize=7,
                  color=style.MUTED)
    left.set_xlim(0.1, 20)
    left.set_ylim(1e-4, 2)
    left.set_xlabel("frequency (Hz)")
    left.set_ylabel("amplitude")
    left.legend(loc="lower left")
    left.set_title("low-pass, 1 Hz, order 4", loc="left")

    zero = acc.copy()
    zero.filter(lowpass=corner)
    causal = acc.copy()
    causal.filter(lowpass=corner, causal=True)
    t = acc.time_axis()
    right.plot(t, zero.data, color=style.SERIES[0], lw=1.4,
               label="zero-phase")
    right.plot(t, causal.data, color=style.SERIES[1], lw=1.4,
               label="causal")
    right.set_xlim(40, 56)
    right.set_xlabel(f"time from {acc.time} (s)")
    right.set_ylabel("acceleration (cm/s²)")
    right.legend(loc="lower right", ncols=2)
    right.set_title("IV.MODE..HNZ (ESM), low-pass 1 Hz", loc="left")
    return style.save(fig, "signals/waveforms_filters.png")


def gaps() -> Path:
    """A stream with a gap, and the record merged with NaN."""
    z = _raw()["IV.MODE..HNZ"][0]
    stream = Stream()
    for start, end in ((0.0, 59.99), (60.0, 119.99), (125.0, None)):
        stream.append(z.extract(start, end))
    merged = stream.extract(merge=True, enforce=True, fillvalue=np.nan)
    reference = stream.starttime
    fig, (top, bottom) = style.figure(2, 1, height=3.2, sharex=True)
    for k, record in enumerate(stream):
        top.plot(record.time_axis(reference), record.data,
                 color=style.SERIES[k], label=f"record {k}")
    gap_start = stream[0].endtime - reference
    gap_end = stream[1].starttime - reference
    for ax in (top, bottom):
        ax.axvspan(gap_start, gap_end, color=style.GRID, zorder=0)
        ax.set_ylabel("counts")
        ax.set_ylim(-30000, 0)
    top.legend(loc="lower right", ncols=2)
    style.label_panel(top, "Stream: two records")
    bottom.plot(merged.time_axis(reference), merged.data,
                color=style.SERIES[0])
    style.label_panel(bottom, "extract(merge=True, enforce=True, "
                              "fillvalue=np.nan)")
    bottom.set_xlabel(f"time from {reference} (s)")
    bottom.set_xlim(100, 145)
    return style.save(fig, "signals/waveforms_gaps.png")


def response() -> Path:
    """Response correction of the raw counts against ESM."""
    data = _raw(metadata=True)
    data.detrend()
    data.taper(5.0)
    data.deconvolve_response(prefilter=(0.05, 0.1, 30.0, 40.0))
    fig, axes = style.figure(3, 1, height=4.2, sharex=True)
    for ax, stream in zip(axes, data):
        record = stream[0]
        esm = _esm(record.sid.split(".")[-1], "ACC")
        reference = record.time
        ax.plot(esm.time_axis(reference), esm.data, color=style.REFERENCE,
                lw=1.6, label="ESM")
        ax.plot(record.time_axis(), record.data * 100,
                color=style.SERIES[0], lw=0.8, label="ShakeLab")
        ax.set_ylabel("cm/s²")
        style.label_panel(ax, record.sid)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper right", ncols=2)
    axes[-1].set_xlim(40, 60)
    axes[-1].set_xlabel(f"time from {data[0][0].time} (s)")
    return style.save(fig, "signals/waveforms_response.png")


def integration() -> Path:
    """Integration of the ESM acceleration against the ESM velocity."""
    acc = _esm("HNZ", "ACC")
    velocity = acc.copy()
    velocity.integrate(method="cum")
    displacement = velocity.copy()
    displacement.integrate(method="cum")
    panels = [(velocity, _esm("HNZ", "VEL"), "velocity (cm/s)"),
              (displacement, _esm("HNZ", "DIS"), "displacement (cm)")]
    fig, axes = style.figure(2, 1, height=3.4, sharex=True)
    for ax, (ours, esm, label) in zip(axes, panels):
        t = ours.time_axis()
        ax.plot(t, esm.data, color=style.REFERENCE, lw=1.6, label="ESM")
        ax.plot(t, ours.data, color=style.SERIES[0], lw=0.8,
                label="integrated acceleration")
        ax.set_ylabel(label)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper right", ncols=2)
    style.label_panel(axes[0], "IV.MODE..HNZ")
    axes[-1].set_xlim(30, 100)
    axes[-1].set_xlabel(f"time from {acc.time} (s)")
    return style.save(fig, "signals/waveforms_integration.png")


FIGURES = {
    "components": components,
    "tolerance": tolerance,
    "filters": filters,
    "gaps": gaps,
    "response": response,
    "integration": integration,
}
