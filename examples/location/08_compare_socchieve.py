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
Compare Socchieve residuals from two velocity models.

Use the same picks and starting geometry for both models. Evaluate
residuals at the catalog hypocenter and after relocation.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shakelab.libutils.geodeticN import WgsPoint, epicentral_distance
from shakelab.location import (
    HomogeneousEngine,
    Hypocenter,
    LayeredTravelTimeEngine,
    locate_event,
)
from shakelab.location.locator import prepare_picks_for_location
from shakelab.location.locutils import load_station_list
from shakelab.location.phases import (
    load_picks,
    normalize_phases,
    read_pun_origin,
)
from shakelab.location.profiles import read_model

data = Path("data")
event_file = data / "Picking" / "20240327211937"
stations = load_station_list(data / "Stations" / "Station_list_OGS.csv")
origin = read_pun_origin(event_file.with_suffix(".pun"))
raw_picks = normalize_phases(load_picks(event_file.with_suffix(".obs")))

# Keep the earliest P and S per station after sorting by arrival time.
earliest = {}
for pick in prepare_picks_for_location(raw_picks, stations):
    if pick.phase in ("P", "S"):
        earliest.setdefault((pick.station, pick.phase), pick)
picks = list(earliest.values())

# Perturb the catalog origin to provide a starting point for relocation.
initial = Hypocenter(
    WgsPoint(origin.longitude - 0.04, origin.latitude + 0.03,
             origin.location.elevation - 2000.0),
    origin.time + 0.5,
)
reference = WgsPoint(origin.longitude, origin.latitude, 0.0)
catalog = Hypocenter(origin.location, origin.time)

distances = np.array([
    epicentral_distance(catalog.location, stations[pick.station].location)
    / 1000.0
    for pick in picks
])
phases = np.array([pick.phase for pick in picks])

profile = data / "Model" / "friuli_velocity_profile_custom.csv"
models = (
    ("Homogeneous", HomogeneousEngine(5900.0, 3300.0)),
    ("Friuli, 20 layers", LayeredTravelTimeEngine(read_model(profile))),
)
fig, axes = plt.subplots(
    2, 2, figsize=(10, 6.6), sharex=True, sharey=True,
    layout="constrained",
)
fig.suptitle("Socchieve · Arrival-time residuals", fontsize=15)
fig.supxlabel("Catalog epicentral distance (km)")
fig.supylabel("Observed − predicted (s)")
station_count = len({pick.station for pick in picks})
print(f"{len(picks)} picks from {station_count} stations")
for column, (label, engine) in enumerate(models):
    # Zero iterations evaluates the catalog solution without moving it.
    _, fixed, _ = locate_event(
        picks, stations, engine, catalog, ref_geo=reference, max_iter=0,
    )
    estimated, free, _ = locate_event(
        picks, stations, engine, initial, ref_geo=reference,
        damping=0.05, max_iter=30, tol_space_m=10.0, tol_time_s=0.005,
    )
    for row, (name, result) in enumerate((
        ("Fixed catalog origin", fixed),
        ("After relocation", free),
    )):
        ax = axes[row, column]
        ax.set_title(f"{label} · {name}\nRMS {result.rms:.3f} s")
        ax.axhline(0, color="0.4", linewidth=0.8)
        ax.grid(axis="y", alpha=0.22)
        for phase, color, marker in (
            ("P", "#1565a5", "o"), ("S", "#ae480e", "^"),
        ):
            mask = phases == phase
            ax.scatter(distances[mask], result.residuals[mask],
                       color=color, marker=marker, s=32, label=phase)
    print(
        f"{label}: fixed RMS {fixed.rms:.3f} s, "
        f"relocated RMS {free.rms:.3f} s, "
        f"depth {estimated.depth_km:.2f} km, "
        f"origin-time shift {estimated.t0 - catalog.t0:+.3f} s, "
        f"converged: {free.converged} "
        f"({free.num_iterations} iterations)"
    )
axes[0, 0].legend(loc="lower left")
# To keep a copy of the figure, uncomment the next line.
# fig.savefig("socchieve.png", dpi=160)
plt.show()
