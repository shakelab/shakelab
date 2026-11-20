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
Figures of the guide docs/signals/spectral.md.

Data: the accelerogram of IV.MODE..HNZ processed by ESM and the raw
counts of IV.MODE (INGV), from tests/signals/data, and synthetic
signals.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import style
from shakelab.signals import fourier, windowing
from shakelab.signals.base import Record
from shakelab.signals.crossspectral import CrossSpectra
from shakelab.signals.io import reader

__all__ = ["FIGURES"]

ACC = ("ESM", "ASCII", "IV.MODE..HNZ.D.IT-2012-0008.ACC.MP.ASC")
CMAP = "Blues"


def _acc() -> Record:
    """Accelerogram of IV.MODE..HNZ processed by ESM (cm/s**2)."""
    return reader(style.TEST_DATA.joinpath(*ACC))[0][0]


def _amplitude(data: np.ndarray, delta: float, nfft: int = None,
               window: str = None) -> tuple:
    """One-sided amplitude calibrated so that a sinusoid A gives A."""
    size = len(data)
    taper = (np.ones(size) if window is None
             else windowing.taper(size, window))
    nfft = size if nfft is None else nfft
    spectrum = np.abs(np.fft.rfft(data * taper, nfft)) * 2 / taper.sum()
    return fourier.frequency_axis(delta, nfft), spectrum


def _image(ax, gram, floor: float = 60.0):
    """Draw a spectrogram in dB with the edge effects shaded."""
    db = gram.db(floor=floor)
    mesh = ax.pcolormesh(gram.time, gram.frequency, db, cmap=CMAP,
                         shading="auto", rasterized=True)
    shade = np.where(gram.coi_mask, 1.0, np.nan)
    ax.pcolormesh(gram.time, gram.frequency, shade, cmap="Greys",
                  vmin=0, vmax=3, alpha=0.55, shading="auto",
                  rasterized=True)
    return mesh


# ---------------------------------------------------------------------------
# Figures

def smoothing() -> Path:
    """Fourier amplitude spectrum, raw and smoothed."""
    acc = _acc()
    spectrum = acc.to_spectrum()
    ko = spectrum.copy()
    ko.logsmooth(method="konno-ohmachi", bandwidth=40.0, phase=False)
    gauss = spectrum.copy()
    gauss.logsmooth(sigma=0.2, phase=False)
    f = spectrum.frequency
    fig, ax = style.figure(1, 1, height=3.0)
    ax.loglog(f[1:], spectrum.amplitude[1:], color=style.REFERENCE,
              lw=0.5, label="Fourier amplitude")
    ax.loglog(f[1:], ko.amplitude[1:], color=style.SERIES[0], lw=1.6,
              label="Konno-Ohmachi, b = 40")
    ax.loglog(f[1:], gauss.amplitude[1:], color=style.SERIES[1], lw=1.2,
              ls="--", label="Gaussian, σ = 0.2")
    ax.set_xlim(0.05, 50)
    ax.set_ylim(1e-3, 30)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("amplitude (cm/s)")
    ax.legend(loc="lower center", ncols=3)
    ax.set_title("IV.MODE..HNZ (ESM), 200 s at 200 Hz", loc="left")
    return style.save(fig, "signals/spectral_smoothing.png")


def leakage() -> Path:
    """A sinusoid between two DFT frequencies, without and with taper."""
    delta, size = 0.01, 1000
    t = np.arange(size) * delta
    data = np.sin(2 * np.pi * 5.15 * t)
    fig, ax = style.figure(1, 1, height=2.8)
    for window, color, label in ((None, style.SERIES[1], "no taper"),
                                 ("hann", style.SERIES[0], "Hann taper")):
        f, amp = _amplitude(data, delta, window=window)
        ax.semilogy(f, amp, "o-", ms=2.5, lw=0.9, color=color, label=label)
    ax.axvline(5.15, color=style.MUTED, lw=0.6, ls=":")
    ax.set_xlim(0, 15)
    ax.set_ylim(1e-5, 2)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("amplitude")
    ax.legend(loc="upper right")
    ax.set_title("Sinusoid of amplitude 1 at 5.15 Hz, 10 s at 100 Hz "
                 "(DFT step 0.1 Hz)", loc="left")
    return style.save(fig, "signals/spectral_leakage.png")


def resolution() -> Path:
    """Resolution against zero padding: two close sinusoids."""
    delta = 0.01
    fig, ax = style.figure(1, 1, height=2.8)
    cases = (
        (2.0, None, style.REFERENCE, "2 s", "o"),
        (2.0, 2000, style.SERIES[1], "2 s padded to 20 s", None),
        (20.0, 20000, style.SERIES[0], "20 s", None),
    )
    for duration, nfft, color, label, marker in cases:
        t = np.arange(int(round(duration / delta))) * delta
        data = (np.sin(2 * np.pi * 5.0 * t)
                + np.sin(2 * np.pi * 5.3 * t))
        f, amp = _amplitude(data, delta, nfft=nfft, window="hann")
        ax.plot(f, amp, color=color, lw=1.4, marker=marker, ms=4,
                label=label)
    for value in (5.0, 5.3):
        ax.axvline(value, color=style.MUTED, lw=0.6, ls=":")
    ax.set_xlim(3.5, 7.0)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("amplitude")
    ax.legend(loc="upper right")
    ax.set_title("Two sinusoids of amplitude 1 at 5.0 and 5.3 Hz, "
                 "Hann taper", loc="left")
    return style.save(fig, "signals/spectral_resolution.png")


def transforms() -> Path:
    """The three time-frequency transforms of the same accelerogram."""
    acc = _acc().extract(25.0, 95.0)
    grams = (
        ("windowed Fourier transform, 2 s windows",
         acc.spectrogram("wft", window=2.0, overlap=0.9, fmin=0.2,
                         fmax=25.0)),
        ("Morlet wavelet transform, ω0 = 6",
         acc.spectrogram("cwt", fmin=0.2, fmax=25.0, nfreq=80)),
        ("Stockwell transform",
         acc.spectrogram("st", fmin=0.2, fmax=25.0, nfreq=80)),
    )
    fig, axes = style.figure(3, 1, height=6.0, sharex=True)
    for ax, (title, gram) in zip(axes, grams):
        mesh = _image(ax, gram)
        ax.set_yscale("log")
        ax.set_ylim(0.2, 25)
        ax.set_ylabel("frequency (Hz)")
        ax.grid(False)
        ax.set_title(title, loc="left")
    axes[-1].set_xlabel(f"time from {acc.time} (s)")
    bar = fig.colorbar(mesh, ax=axes, shrink=0.6)
    bar.set_label("dB (re 1 cm/s²)")
    return style.save(fig, "signals/spectral_transforms.png")


def window_length() -> Path:
    """Short and long windows of the windowed Fourier transform."""
    delta = 0.01
    t = np.arange(3000) * delta
    data = (np.sin(2 * np.pi * 4.0 * t) + np.sin(2 * np.pi * 4.6 * t))
    pulse = np.exp(-((t - 15.0) / 0.03) ** 2) * 20.0
    record = Record(data + pulse, delta=delta)
    fig, axes = style.figure(1, 2, height=3.0, sharey=True)
    for ax, window in zip(axes, (0.5, 5.0)):
        gram = record.spectrogram("wft", window=window, overlap=0.9,
                                  fmax=10.0)
        _image(ax, gram, floor=40.0)
        ax.set_ylim(0, 10)
        ax.grid(False)
        ax.set_xlabel("time (s)")
        ax.set_title(f"{window:g} s windows", loc="left")
    axes[0].set_ylabel("frequency (Hz)")
    fig.suptitle("Sinusoids at 4.0 and 4.6 Hz and a pulse at 15 s",
                 x=0.01, ha="left", fontsize=9)
    return style.save(fig, "signals/spectral_window_length.png")


def welch() -> Path:
    """Periodogram of one window against Welch's average."""
    data = reader(style.TEST_DATA / "MODE.mseed",
                  metadata=style.TEST_DATA / "MODE.xml")
    data.detrend()
    data.remove_sensitivity()
    noise = data["IV.MODE..HNZ"][0].extract(0.0, 35.0)
    grid = fourier.frequency_range(0.3, 40.0, 400)
    fig, ax = style.figure(1, 1, height=3.0)
    cases = ((34.0, None, style.SERIES[1], 0.9,
              "one 34 s window (periodogram)"),
             (5.0, 0.5, style.SERIES[0], 1.6,
              "5 s windows, 50% overlap: 13 averaged"))
    for duration, overlap, color, width, label in cases:
        spectra = CrossSpectra.from_data(
            noise.data[np.newaxis], noise.delta, grid, duration=duration,
            overlap=overlap, window="hann")
        psd = spectra.mean()[:, 0, 0].real
        ax.semilogx(grid, 10 * np.log10(psd), color=color, lw=width,
                    label=label)
    ax.set_xlim(0.3, 40)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (dB re 1 (m/s²)²/Hz)")
    ax.legend(loc="lower left")
    ax.set_title("IV.MODE..HNZ, 35 s of noise before the earthquake",
                 loc="left")
    return style.save(fig, "signals/spectral_welch.png")


FIGURES = {
    "smoothing": smoothing,
    "leakage": leakage,
    "resolution": resolution,
    "transforms": transforms,
    "window_length": window_length,
    "welch": welch,
}
