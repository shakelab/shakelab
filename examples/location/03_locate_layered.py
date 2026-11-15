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
Relocate Socchieve using a 1D layered velocity model.

Use direct arrivals and head waves from the supplied Friuli profile.
The input picks and initial geometry are those of the homogeneous run.
"""

from pathlib import Path

from shakelab.libutils.geodeticN import WgsPoint
from shakelab.location import (
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

profile = data / "Model" / "friuli_velocity_profile_custom.csv"
model = read_model(profile)
engine = LayeredTravelTimeEngine(model)
event, result, _ = locate_event(
    picks, stations, engine, initial, ref_geo=reference,
    damping=0.05, max_iter=30, tol_space_m=10, tol_time_s=0.005,
)
print(f"{model.num_layers} layers, converged: {result.converged}")
print(f"Longitude {event.lon:.4f}, latitude {event.lat:.4f}")
print(f"Depth {event.depth_km:.2f} km, origin {event.t0}")
print(f"RMS {result.rms:.3f} s for {result.num_picks} picks")
