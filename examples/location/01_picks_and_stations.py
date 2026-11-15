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
Read station coordinates, phase picks and a catalog hypocenter.

The sample files are from the Socchieve event of 27 March 2024.
The picks are normalized to P and S for the travel-time engines.
"""

from pathlib import Path

from shakelab.location.locutils import load_station_list
from shakelab.location.phases import (
    load_picks,
    normalize_phases,
    read_pun_origin,
)

data = Path("data")
event_file = data / "Picking" / "20240327211937"
stations = load_station_list(data / "Stations" / "Station_list_OGS.csv")
origin = read_pun_origin(event_file.with_suffix(".pun"))
picks = normalize_phases(load_picks(event_file.with_suffix(".obs")))
usable = [
    pick for pick in picks
    if pick.station in stations and pick.time is not None
    and pick.phase in ("P", "S")
]
print(f"Catalog origin {origin.time}, depth {origin.depth_km:.2f} km")
print(f"{len(usable)} timed P/S picks with station coordinates")
for pick in usable[:5]:
    station = stations[pick.station]
    print(f"{pick.station:5} {pick.phase} {pick.time}"
          f"  ({station.lon:.4f} E, {station.lat:.4f} N)")

# The same reader dispatches other pick formats by extension.
for extension in (".dat", ".xml", ".hpl"):
    alternate = load_picks(event_file.with_suffix(extension))
    print(f"{extension}: {len(alternate)} picks")
