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
Text format of the FDSN station service.

The station service returns, with ``format=text``, one line per station
or channel, with fields separated by ``|`` and a header line starting
with ``#`` (fdsnws-station 1.1, section "Text output"):

- station level: ``Network | Station | Latitude | Longitude |
  Elevation | SiteName | StartTime | EndTime``;
- channel level: ``Network | Station | Location | Channel | Latitude |
  Longitude | Elevation | Depth | Azimuth | Dip | SensorDescription |
  Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime``.

It is much lighter than StationXML. The channel lines become channel
epochs with coordinates, burial depth (converted from m to km),
orientation, sampling rate, sensor description and the overall
sensitivity (``Scale`` per ``ScaleUnits`` at ``ScaleFreq``), without
response stages and without station epochs; the station lines become
station epochs.
"""

from __future__ import annotations

from typing import Any, Optional

from shakelab.libutils.ascii import read_text
from shakelab.signals import metadata as mdm
from shakelab.signals import response as rspm

__all__ = ["read_station_text"]

_STATION_FIELDS = 8
_CHANNEL_FIELDS = 17


def read_station_text(source: Any) -> mdm.MetadataCollection:
    """
    Read the text format of the FDSN station service.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        File path, content (text with line breaks, or bytes) or open
        file.

    Returns
    -------
    MetadataCollection
        Station epochs (station level) or channel epochs (channel
        level).

    Raises
    ------
    ValueError
        If a line has neither the station nor the channel fields.

    Examples
    --------
    >>> text = (
    ...     "#Network|Station|Latitude|Longitude|Elevation|SiteName|"
    ...     "StartTime|EndTime\\n"
    ...     "IV|MODE|44.8|11.1|10.0|Modena|2008-01-01T00:00:00|\\n"
    ... )
    >>> read_station_text(text).stations()
    [StationEpoch('IV.MODE', 2008-01-01T00:00:00.000000Z -> open)]
    """
    collection = mdm.MetadataCollection()

    for number, line in enumerate(read_text(source).splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [field.strip() for field in line.split("|")]
        if len(fields) == _CHANNEL_FIELDS:
            collection.append(_channel(fields))
        elif len(fields) == _STATION_FIELDS:
            collection.append(_station(fields))
        else:
            raise ValueError(
                f"line {number}: {len(fields)} fields; expected "
                f"{_STATION_FIELDS} (stations) or {_CHANNEL_FIELDS} "
                "(channels)"
            )

    return collection


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _station(fields: list) -> mdm.StationEpoch:
    """Station epoch of a station line."""
    network, station, latitude, longitude, elevation = fields[:5]
    return mdm.StationEpoch(
        f"{network}.{station}",
        starttime=fields[6] or None,
        endtime=fields[7] or None,
        location=_location(longitude, latitude, elevation),
        site=fields[5] or None,
    )


def _channel(fields: list) -> mdm.ChannelEpoch:
    """Channel epoch of a channel line."""
    network, station, location, channel = fields[:4]
    latitude, longitude, elevation, depth, azimuth, dip = fields[4:10]
    description, scale, frequency, units, rate = fields[10:15]

    sample_rate = _float(rate)
    stage_set = None
    if _float(scale) is not None:
        stage_set = rspm.StageSet(sample_rate=sample_rate)
        stage_set.sensitivity = rspm.Sensitivity(
            value=_float(scale),
            frequency=_float(frequency),
            input_units=units or None,
        )

    depth_m = _float(depth)
    azimuth_deg = _float(azimuth)
    dip_deg = _float(dip)

    return mdm.ChannelEpoch(
        f"{network}.{station}.{location}.{channel}",
        starttime=fields[15] or None,
        endtime=fields[16] or None,
        location=_location(longitude, latitude, elevation),
        depth=None if depth_m is None else depth_m / 1000.0,
        orientation=(
            None if azimuth_deg is None or dip_deg is None
            else (azimuth_deg, dip_deg)
        ),
        sample_rate=sample_rate,
        sensor=(
            mdm.Equipment(description=description) if description
            else None
        ),
        response=stage_set,
    )


def _location(longitude: str, latitude: str, elevation: str):
    """``(longitude, latitude, elevation)`` or None."""
    lon, lat = _float(longitude), _float(latitude)
    if lon is None or lat is None:
        return None
    elev = _float(elevation)
    return (lon, lat, 0.0 if elev is None else elev)


def _float(text: str) -> Optional[float]:
    """Value of a numeric field; None if empty."""
    return float(text) if text else None
