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
Figures of the guide docs/signals/response.md.

Data: StationXML and raw counts of IV.MODE (INGV) and the RESP file of
IU.ANMO.00.BHZ, from tests/signals/data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import style
from shakelab.signals.fourier import cosine_taper, inverse_spectrum
from shakelab.signals.io import reader
from shakelab.signals.metadata import MetadataCollection
from shakelab.signals.response import (
    Sensitivity,
    StageCoefficients,
    StageGain,
    StagePoleZero,
    StageSet,
)

__all__ = ["FIGURES"]

RESP = style.TEST_DATA / "RESP"


def _mode() -> StageSet:
    """Response of IV.MODE..HNZ."""
    metadata = MetadataCollection(style.TEST_DATA / "MODE.xml")
    return metadata.get("IV.MODE..HNZ").response


def _anmo() -> StageSet:
    """Response of IU.ANMO.00.BHZ on 1 January 2015."""
    metadata = MetadataCollection(RESP / "RESP.IU.ANMO.00.BHZ")
    return metadata.get("IU.ANMO.00.BHZ", "2015-01-01").response


def _phase(value: np.ndarray, floor: float = 1e-3) -> np.ndarray:
    """Unwrapped phase in degrees, NaN where the amplitude is negligible."""
    phase = np.degrees(np.unwrap(np.angle(value)))
    amplitude = np.abs(value)
    return np.where(amplitude >= floor * amplitude.max(), phase, np.nan)


def _stage_response(response: StageSet, stage, frequency) -> np.ndarray:
    """Transfer function of one stage normalized to unit gain."""
    gain_frequency = next(
        (s.frequency for s in response.stage
         if isinstance(s, StageGain)
         and s.stage_number == stage.stage_number), None)
    value = stage.transfer_function(frequency)
    if gain_frequency is not None:
        value = value / abs(stage.transfer_function([gain_frequency])[0])
    return value


# ---------------------------------------------------------------------------
# Figures

def bode() -> Path:
    """Amplitude and phase of the complete and analog responses."""
    response = _mode()
    frequency = np.logspace(-2, np.log10(50.0), 600)
    full = response.transfer_function(frequency)
    analog = response.transfer_function(frequency, stages=("gain", "paz"))
    sens = response.sensitivity
    fig, (top, bottom) = style.figure(2, 1, height=4.0, sharex=True)
    top.loglog(frequency, np.abs(full), color=style.SERIES[0], lw=1.6,
               label="complete response")
    top.loglog(frequency, np.abs(analog), color=style.SERIES[1], ls="--",
               label="analog stages (gain, paz)")
    top.plot([sens.frequency], [sens.value], "o", ms=5,
             color=style.TEXT, label="declared sensitivity")
    top.set_ylim(1e4, 1e6)
    top.set_ylabel("amplitude (count per m/s²)")
    top.legend(loc="lower left")
    style.label_panel(top, "IV.MODE..HNZ")
    bottom.semilogx(frequency, _phase(full), color=style.SERIES[0],
                    lw=1.6)
    bottom.semilogx(frequency, _phase(analog), color=style.SERIES[1],
                    ls="--")
    bottom.set_ylabel("phase (degrees)")
    bottom.set_xlabel("frequency (Hz)")
    for ax in (top, bottom):
        ax.axvline(50.0, color=style.MUTED, lw=0.6, ls=":")
    bottom.text(48.0, bottom.get_ylim()[0], "Nyquist", ha="right",
                va="bottom", fontsize=7, color=style.MUTED)
    bottom.set_xlim(1e-2, 55.0)
    return style.save(fig, "signals/response_bode.png")


def stages() -> Path:
    """Contribution of each stage, normalized to unit gain."""
    response = _mode()
    frequency = np.logspace(-1, np.log10(100.0), 800)
    labels = {1: "stage 1, sensor (poles and zeros)",
              3: "stage 3, FIR, 800 → 200 Hz",
              4: "stage 4, FIR, 200 → 100 Hz"}
    fig, ax = style.figure(1, 1, height=3.0)
    k = 0
    for stage in response.stage:
        if stage.stage_number not in labels or isinstance(stage, StageGain):
            continue
        value = _stage_response(response, stage, frequency)
        if stage.stage_number == 4:
            value = np.where(frequency <= 100.0, value, np.nan)
        ax.loglog(frequency, np.abs(value), color=style.SERIES[k],
                  label=labels[stage.stage_number])
        k += 1
    total = response.transfer_function(frequency) / response.sensitivity.value
    ax.loglog(frequency, np.abs(total), color=style.TEXT, lw=1.4,
              ls="--", label="complete response / sensitivity")
    ax.axvline(50.0, color=style.MUTED, lw=0.6, ls=":")
    ax.text(48.0, 2e-4, "Nyquist of the data", ha="right", fontsize=7,
            color=style.MUTED)
    ax.set_ylim(1e-4, 3)
    ax.set_xlim(0.1, 100.0)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("normalized amplitude")
    ax.legend(loc="lower left")
    return style.save(fig, "signals/response_stages.png")


def fir_phase() -> Path:
    """Phase of a linear-phase FIR with and without the correction."""
    response = _mode()
    stage = next(s for s in response.stage
                 if isinstance(s, StageCoefficients)
                 and s.stage_number == 4)
    uncorrected = StageCoefficients(stage.to_dict())
    uncorrected.correction = 0.0
    frequency = np.linspace(0.0, 50.0, 1001)
    with_corr = stage.transfer_function(frequency)
    without = uncorrected.transfer_function(frequency)
    fig, ax = style.figure(1, 1, height=2.8)
    ax.plot(frequency, _phase(without), color=style.SERIES[1],
            label="without correction")
    ax.plot(frequency, _phase(with_corr),
            color=style.SERIES[0], lw=1.4,
            label=f"with the declared correction ({stage.correction} s)")
    ax.set_xlim(0, 50)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("phase (degrees)")
    ax.legend(loc="lower left")
    ax.set_title(f"IV.MODE..HNZ, stage 4: FIR of "
                 f"{np.size(stage.numerator)} coefficients at "
                 f"{stage.input_sample_rate:g} Hz, delay "
                 f"{stage.delay} s", loc="left")
    return style.save(fig, "signals/response_fir_phase.png")


def waterlevel() -> Path:
    """Regularized inverse of a broadband response."""
    response = _anmo()
    frequency = np.logspace(-4, 1, 800)
    h = response.transfer_function(frequency)
    s0 = response.sensitivity.value
    fig, ax = style.figure(1, 1, height=3.2)
    ax.loglog(frequency, s0 / np.abs(h), color=style.REFERENCE, lw=1.6,
              label="exact inverse")
    styles = (("water", style.SERIES[0], "-"),
              ("smooth", style.SERIES[1], "--"),
              ("sharp", style.SERIES[2], ":"))
    for method, color, line in styles:
        inverse = inverse_spectrum(h, waterlevel=40.0, method=method)
        ax.loglog(frequency, s0 * np.abs(inverse), color=color, ls=line,
                  lw=1.2, label=f'"{method}", 40 dB')
    window = cosine_taper(frequency, (0.005, 0.01, 8.0, 9.5))
    inverse = inverse_spectrum(h, waterlevel=40.0) * window
    ax.loglog(frequency, s0 * np.abs(inverse), color=style.TEXT, lw=1.4,
              label='"water" with prefilter (0.005, 0.01, 8, 9.5) Hz')
    ax.set_ylim(1e-3, 1e5)
    ax.set_xlim(1e-4, 10)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("|inverse| × sensitivity")
    ax.legend(loc="upper center")
    ax.set_title("IU.ANMO.00.BHZ (velocity sensor, 20 Hz)", loc="left")
    return style.save(fig, "signals/response_waterlevel.png")


def padding() -> Path:
    """Wrap-around at the ends of a short record."""
    data = reader(style.TEST_DATA / "MODE.mseed",
                  metadata=style.TEST_DATA / "MODE.xml")
    record = data["IV.MODE..HNZ"][0]
    prefilter = (0.05, 0.1, 30.0, 40.0)
    reference = record.copy()
    reference.detrend()
    reference.taper(5.0)
    reference.deconvolve_response(prefilter=prefilter)
    reference = reference.extract(40.0, 50.0)
    fig, ax = style.figure(1, 1, height=2.8)
    for k, padded in enumerate((False, True)):
        window = record.extract(40.0, 50.0)
        window.detrend()
        window.deconvolve_response(prefilter=prefilter, padding=padded)
        error = (window.data - reference.data) * 100.0
        ax.plot(window.time_axis(record.time), error,
                color=style.SERIES[1 - k],
                label="padding=True" if padded else "padding=False")
    ax.set_xlim(40, 50)
    ax.set_xlabel(f"time from {record.time} (s)")
    ax.set_ylabel("difference (cm/s²)")
    ax.legend(loc="lower left", ncols=2)
    ax.set_title("IV.MODE..HNZ: 10 s corrected on their own, minus the "
                 "same samples corrected in the full record (peak "
                 f"{reference.peak * 100:.1f} cm/s²)", loc="left")
    return style.save(fig, "signals/response_padding.png")


def _geophone() -> StageSet:
    """Response of a 1 Hz geophone, 100 V per m/s (as in the guide)."""
    w0, damping = 2 * np.pi * 1.0, 0.707
    poles = [complex(-damping * w0, w0 * np.sqrt(1 - damping ** 2)),
             complex(-damping * w0, -w0 * np.sqrt(1 - damping ** 2))]
    response = StageSet()
    response.append([
        StagePoleZero(poles=poles, zeros=[0, 0], normalization_factor=1.0,
                      normalization_frequency=10.0, input_units="m/s",
                      output_units="V", stage_number=1),
        StageGain(sensitivity=100.0, frequency=10.0, stage_number=1),
    ])
    response.sensitivity = Sensitivity(100.0, 10.0, "m/s", "V")
    return response


def sensitivity() -> Path:
    """Sensitivity-only and full correction of a simulated geophone."""
    data = reader(style.TEST_DATA / "MODE.mseed",
                  metadata=style.TEST_DATA / "MODE.xml")
    record = data["IV.MODE..HNZ"][0]
    prefilter = (0.05, 0.1, 30.0, 40.0)
    record.detrend()
    record.taper(5.0)
    record.deconvolve_response(prefilter=prefilter)
    truth = record.copy()
    truth.integrate()
    truth.head.units = "m/s"
    geophone = _geophone()
    recorded = truth.copy()
    recorded.convolve_response(geophone)
    quick = recorded.copy()
    quick.remove_sensitivity(geophone)
    full = recorded.copy()
    full.deconvolve_response(geophone, prefilter=prefilter)

    fig, (left, right) = style.figure(1, 2, height=3.0,
                                      width_ratios=(1, 2.2))
    frequency = np.logspace(-2, np.log10(50.0), 400)
    left.loglog(frequency, np.abs(geophone.transfer_function(frequency)),
                color=style.SERIES[0], lw=1.4, label="response |H|")
    left.axhline(geophone.sensitivity.value, color=style.SERIES[1],
                 ls="--", label="sensitivity")
    left.axvline(1.0, color=style.MUTED, lw=0.6, ls=":")
    left.set_ylim(1e-2, 300)
    left.set_xlim(1e-2, 50)
    left.set_xlabel("frequency (Hz)")
    left.set_ylabel("V per m/s")
    left.legend(loc="lower right")
    left.set_title("1 Hz geophone", loc="left")

    t = truth.time_axis()
    right.plot(t, truth.data * 100, color=style.REFERENCE, lw=2.2,
               label="true velocity")
    right.plot(t, quick.data * 100, color=style.SERIES[1], lw=1.0,
               label="remove_sensitivity()")
    right.plot(t, full.data * 100, color=style.SERIES[0], lw=1.0,
               label="deconvolve_response()")
    right.set_xlim(40, 75)
    right.set_xlabel(f"time from {truth.time} (s)")
    right.set_ylabel("velocity (cm/s)")
    handles, labels = right.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower right", ncols=3)
    right.set_title("IV.MODE..HNZ recorded by the geophone, then "
                    "corrected", loc="left")
    return style.save(fig, "signals/response_sensitivity.png")


FIGURES = {
    "bode": bode,
    "stages": stages,
    "fir_phase": fir_phase,
    "waterlevel": waterlevel,
    "padding": padding,
    "sensitivity": sensitivity,
}
