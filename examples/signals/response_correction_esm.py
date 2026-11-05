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
Instrument response correction compared with ESM.

The raw strong-motion records of station IV.MODE (Mw 6.1 Emilia
earthquake, 20 May 2012) are corrected for the complete instrument
response read from StationXML, and compared with the accelerations
processed by the Engineering Strong-Motion database (ESM).

Steps
-----
1. Read the StationXML metadata, check their consistency, and read the
   raw records (counts) assigning the metadata to them.
2. Read the ESM accelerations (DYNA 1.2 ASCII files, cm/s**2) and,
   from their header, the band-pass filter applied by ESM to each
   component.
3. Remove the mean, taper, filter in the same band and deconvolve the
   response.
4. Plot the corrected accelerations over the ESM ones on a common time
   axis, and their difference.

Run from the root of the repository::

    python examples/signals/response_correction_esm.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.signals.io import metadata_reader, reader
from shakelab.signals.libio.dyna import dynaread

DATA = Path(__file__).resolve().parents[2] / "tests" / "signals" / "data"
ESM = DATA / "ESM" / "ASCII"
COMPONENTS = ("HNE", "HNN", "HNZ")


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------

meta = metadata_reader(str(DATA / "MODE.xml"))
issues = meta.validate()
print("Metadata issues:", issues if issues else "none")

raw = reader(str(DATA / "MODE.mseed"), metadata=meta)

esm = {}
for cha in COMPONENTS:
    record = dynaread(ESM / f"IV.MODE..{cha}.D.IT-2012-0008.ACC.MP.ASC")
    header = record.head.extra["dyna"]
    esm[cha] = {
        "record": record,
        "band": (
            header["LOW_CUT_FREQUENCY_HZ"],
            header["HIGH_CUT_FREQUENCY_HZ"],
        ),
        "order": header["FILTER_ORDER"],
    }


# ----------------------------------------------------------------------------
# Response correction, with the ESM filter of each component
# ----------------------------------------------------------------------------

acc = {}
for cha in COMPONENTS:
    rec = raw[f"IV.MODE..{cha}"][0].copy()
    low, high = esm[cha]["band"]
    rec.remove_mean()
    rec.taper(5.0)
    rec.filter(highpass=low, lowpass=high, order=esm[cha]["order"])
    rec.deconvolve_response()
    acc[cha] = rec
    print(f"{cha}: band {low}-{high} Hz, PGA {rec.peak:.4f} m/s**2 "
          f"(ESM {esm[cha]['record'].peak / 100.0:.4f})")

# Common time reference for both data sets.
t0 = acc[COMPONENTS[0]].time


# ----------------------------------------------------------------------------
# Comparison
# ----------------------------------------------------------------------------

fig, axes = plt.subplots(len(COMPONENTS), 1, sharex=True, figsize=(11, 8))
for ax, cha in zip(axes, COMPONENTS):
    rec = acc[cha]
    ref = esm[cha]["record"]

    ax.plot(ref.time_axis(reference=t0), ref.data / 100.0, lw=1.2,
            color="0.6", label=f"ESM (PGA {ref.peak / 100.0:.4f})")
    ax.plot(rec.time_axis(reference=t0), rec.data, lw=0.6, color="C0",
            label=f"ShakeLab (PGA {rec.peak:.4f})")
    ax.set_ylabel(f"{cha} [m/s**2]")
    ax.legend(loc="upper right", fontsize=8)

axes[-1].set_xlabel(f"Time from {t0} [s]")
axes[-1].set_xlim(30, 90)
fig.suptitle("IV.MODE, 2012-05-20: response correction vs ESM")
fig.tight_layout()

# ESM is sampled at 200 Hz from the same start time: every second ESM
# sample coincides with a sample of the corrected records.
fig, axes = plt.subplots(len(COMPONENTS), 1, sharex=True, figsize=(11, 8))
for ax, cha in zip(axes, COMPONENTS):
    rec = acc[cha]
    ref = esm[cha]["record"].data[::2] / 100.0
    n = min(rec.nsamp, ref.size)

    diff = rec.data[:n] - ref[:n]
    ratio = np.sqrt(np.mean(diff ** 2) / np.mean(ref[:n] ** 2))
    ax.plot(rec.time_axis(reference=t0)[:n], diff, lw=0.6, color="C3")
    ax.set_ylabel(f"{cha} diff [m/s**2]")
    ax.set_title(f"RMS difference / RMS ESM = {ratio:.3f}", fontsize=8,
                 loc="right")

axes[-1].set_xlabel(f"Time from {t0} [s]")
fig.suptitle("ShakeLab minus ESM")
fig.tight_layout()

plt.show()
