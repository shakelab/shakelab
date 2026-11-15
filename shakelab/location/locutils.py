# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Utility helpers for loading station metadata for the location module.

Currently this module provides a simple reader for station lists
stored in a semicolon-separated CSV file, returning a dictionary of
:class:`shakelab.location.locator.Station` instances keyed by
station code.

The expected CSV structure is::

    code;network;latitude;longitude;elevation;...

where:
    * latitude, longitude are in degrees (WGS84);
    * elevation is in meters, positive upward;
    * the decimal separator may be either '.' or ','.

Extra columns are ignored. Lines that cannot be parsed are skipped.
"""

from __future__ import annotations

from typing import Dict

from shakelab.location.locator import Station  # type: ignore
from shakelab.libutils.geodeticN import WgsPoint  # type: ignore


def load_station_list(csv_path: str) -> Dict[str, Station]:
    """
    Load a station list from a semicolon-separated CSV file.

    The CSV file is expected to contain a header line, followed by
    one line per station with at least five columns:

        code;network;latitude;longitude;elevation;...

    Only the first, third, fourth and fifth columns are used:

        * code      -> station identifier (string),
        * latitude  -> latitude in degrees (float),
        * longitude -> longitude in degrees (float),
        * elevation -> elevation in meters, positive upward (float).

    Returns a dictionary mapping station codes to
    :class:`Station` objects, where each station stores its
    geographic position as a :class:`WgsPoint`.

    Parameters
    ----------
    csv_path : str
        Path to the station CSV file.

    Returns
    -------
    dict
        Dictionary ``{station_code: Station}``.

    Notes
    -----
    The reader is tolerant to European decimal commas: any ',' in
    the numeric fields is replaced by '.' before parsing. Lines
    with missing or non-numeric coordinates are silently skipped.
    """
    stations: Dict[str, Station] = {}

    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        # Skip header line
        next(f, None)

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split(";")
            if len(parts) < 5:
                # Not enough columns to extract coordinates
                continue

            code = parts[0].strip()
            if not code:
                continue

            lat_str = parts[2].strip().replace(",", ".")
            lon_str = parts[3].strip().replace(",", ".")
            elev_str = parts[4].strip().replace(",", ".")

            try:
                lat = float(lat_str)
                lon = float(lon_str)
                elev = float(elev_str)
            except ValueError:
                # Skip malformed numeric fields
                continue

            location = WgsPoint(latitude=lat, longitude=lon, elevation=elev)
            stations[code] = Station(station_id=code, location=location)

    return stations
