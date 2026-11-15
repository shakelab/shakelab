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
Phase and origin I/O utilities for ShakeLab location.

This module provides a unified in-memory representation for phase
picks (and event origins) used by the location subsystem, together
with reader functions for several file formats adopted by OGS and
Antelope.

Supported formats
-----------------
- QuakeML (.xml): full event description with picks.
- OBS (.obs): tabular phase picks (Antelope/OGS style).
- DAT (.dat): compact Antelope/OGS station-phase summary.
- HPL (.hpl): "hypo line" pick listing with event header.
- PUN (.pun): compact event header with origin information.

Time handling
-------------
All times are represented using :class:`shakelab.libutils.timeN.Date`.
The internal epoch and calendar logic are delegated to that class. In
practice, all times in these files are treated as UTC; any original
time-zone information is ignored or normalized to UTC on input.

Geographical representation
---------------------------
Event locations are represented using :class:`shakelab.libutils.
geodetic.WgsPoint`, which stores latitude, longitude and elevation
in meters.

In this module we adopt the following convention:

* ``latitude``: decimal degrees (WGS84).
* ``longitude``: decimal degrees (WGS84, positive east).
* ``elevation``: meters, **positive above** the reference surface
  (e.g. sea level or topography), and **negative at depth** for
  hypocenters.

The hypocentral depth in kilometers (positive downward) is provided
as a read-only property ``EventOrigin.depth_km`` for convenience and
backward compatibility.

Typical usage
-------------
    from shakelab.location.phases import (
        load_picks,
        read_obs_picks,
        read_pun_origin,
    )

    picks_obs = read_obs_picks("20240327211937.obs")
    origin = read_pun_origin("20240327211937.pun")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re
import xml.etree.ElementTree as ET

from shakelab.libutils.timeN import Date
from shakelab.libutils.geodeticN import WgsPoint  # type: ignore


# ---------------------------------------------------------------------------
# Common data structures
# ---------------------------------------------------------------------------


@dataclass
class PhasePick:
    """
    Container for a single phase pick.

    Attributes
    ----------
    station : str
        Station code.
    phase : str
        Phase name (e.g. "P", "S", "Pg", "Sg").
    time : Optional[Date]
        Absolute pick time, represented as a :class:`Date` instance,
        or ``None`` if not available.
    network : Optional[str]
        Network code, if known.
    location : Optional[str]
        Location code, if known.
    channel : Optional[str]
        Channel code, if known.
    onset : Optional[str]
        Onset quality / type (e.g. "e", "i").
    polarity : Optional[str]
        First-motion polarity (e.g. "+", "-", "undecidable").
    weight : Optional[float]
        Numeric weight, if available. A common choice is ``1 / sigma``,
        where ``sigma`` is the expected picking uncertainty in seconds.
    picker : Optional[str]
        Picker code or evaluation mode (e.g. "GAU", "manual").
    extra : Dict[str, object]
        Additional, format-specific information.

    Notes
    -----
    This data structure is intentionally minimal and independent from
    any particular location algorithm. The :mod:`locator` module uses
    only ``station``, ``phase``, ``time`` and ``weight`` directly; all
    other fields are treated as metadata.
    """

    station: str
    phase: str
    time: Optional[Date]

    network: Optional[str] = None
    location: Optional[str] = None
    channel: Optional[str] = None

    onset: Optional[str] = None
    polarity: Optional[str] = None
    weight: Optional[float] = None
    picker: Optional[str] = None

    extra: Dict[str, object] = field(default_factory=dict)


@dataclass
class EventOrigin:
    """
    Basic information about an event origin.

    Attributes
    ----------
    time : Date
        Origin time, represented as a :class:`Date` instance.
    location : WgsPoint
        Epicentral / hypocentral position in WGS84 coordinates.
        The latitude and longitude are in decimal degrees, while
        ``location.elevation`` is in meters. Elevation is **negative**
        when the hypocenter lies below the reference surface.
    magnitude : Optional[float]
        Scalar magnitude, if present in the file.
    extra : Dict[str, object]
        Additional, format-specific information (gaps, errors, etc.).

    Notes
    -----
    Older versions of ShakeLab used separate attributes
    ``latitude``, ``longitude`` and ``depth_km``. For backward
    compatibility, these are still exposed as read-only properties
    derived from :attr:`location`:

    * ``latitude``  → ``location.latitude``
    * ``longitude`` → ``location.longitude``
    * ``depth_km``  → ``-location.elevation / 1000.0``
    """

    time: Date
    location: WgsPoint
    magnitude: Optional[float] = None
    extra: Dict[str, object] = field(default_factory=dict)

    @property
    def latitude(self) -> float:
        """Epicenter latitude (degrees, read-only)."""
        return self.location.latitude

    @property
    def longitude(self) -> float:
        """Epicenter longitude (degrees, read-only)."""
        return self.location.longitude

    @property
    def depth_km(self) -> float:
        """
        Hypocentral depth in kilometers (positive downward, read-only).

        The value is derived from ``-location.elevation / 1000.0``,
        assuming that ``location.elevation`` is negative at depth.
        """
        return -self.location.elevation / 1000.0


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _build_date(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    seconds: float,
) -> Date:
    """
    Build a :class:`Date` from integer fields and seconds.

    Parameters
    ----------
    year, month, day : int
        Calendar date.
    hour, minute : int
        Time of day.
    seconds : float
        Seconds (can include a fractional part).

    Returns
    -------
    Date
        A new :class:`Date` instance with the requested components.
    """
    return Date([year, month, day, hour, minute, float(seconds)])


def _parse_pun_datetime_code(code: str) -> Tuple[int, int, int, int, int]:
    """
    Parse the 11-digit Antelope/OGS date-time code used in .pun/.dat files.

    The code has the form::

        f y y M M d d h h m m

    where the first digit ``f`` is a flag and the remaining 10 digits
    encode year (2 digits), month, day, hour, minute.

    Example
    -------
    ``"12403272119" -> 2024-03-27 21:19``
    """
    if len(code) != 11 or not code.isdigit():
        raise ValueError(f"Invalid datetime code: {code!r}")

    yy = int(code[1:3])
    mm = int(code[3:5])
    dd = int(code[5:7])
    hh = int(code[7:9])
    mi = int(code[9:11])
    year = 2000 + yy

    return year, mm, dd, hh, mi


def _parse_hpl_date_token(token: str) -> Tuple[int, int, int]:
    """
    Parse the date token used in .hpl header lines.

    The token has the form::

        f f y y M M d d

    where the first two digits are flags and the remaining digits
    encode year (2 digits), month, and day.

    Example
    -------
    ``"11240327" -> 2024-03-27``
    """
    if len(token) != 8 or not token.isdigit():
        raise ValueError(f"Invalid HPL date token: {token!r}")

    yy = int(token[2:4])
    mm = int(token[4:6])
    dd = int(token[6:8])
    year = 2000 + yy

    return year, mm, dd


def _parse_deg_min(value: str) -> float:
    """
    Convert a ``deg-min`` string (e.g. ``"46-21.53"``) to decimal degrees.

    The format used in OGS .pun/.hpl files is::

        ``DD-MM.mm``

    where ``MM.mm`` are arc minutes. The sign is taken from the degree
    part, so negative degrees produce negative coordinates.
    """
    if "-" not in value:
        # Already in decimal degrees
        return float(value)

    deg_str, min_str = value.split("-", 1)
    deg = float(deg_str)
    minutes = float(min_str)
    sign = -1.0 if deg < 0.0 else 1.0
    return sign * (abs(deg) + minutes / 60.0)


# ---------------------------------------------------------------------------
# QuakeML (.xml)
# ---------------------------------------------------------------------------


def read_quakeml_picks(path: str | Path) -> List[PhasePick]:
    """
    Read all picks from a QuakeML file.

    The time values parsed from QuakeML are converted to
    :class:`Date` instances. Any explicit time-zone information is
    interpreted and normalized to UTC by the :class:`Date`
    constructor, which internally handles ISO8601 strings.

    Parameters
    ----------
    path : str or Path
        Path to the QuakeML (.xml) file.

    Returns
    -------
    List[PhasePick]
        List of picks found in the file.
    """
    path = Path(path)
    tree = ET.parse(path)
    root = tree.getroot()

    picks: List[PhasePick] = []

    # Use wildcard namespaces so we do not depend on a specific prefix.
    for pick_elem in root.findall(".//{*}pick"):
        # Time
        time_elem = pick_elem.find(".//{*}time/{*}value")
        if time_elem is None or time_elem.text is None:
            pick_time: Optional[Date] = None
        else:
            t_str = time_elem.text.strip()
            # Date understands ISO 8601, including trailing 'Z'
            pick_time = Date(t_str)

        # Station information
        wf = pick_elem.find(".//{*}waveformID")
        network = wf.get("networkCode") if wf is not None else None
        station = wf.get("stationCode") if wf is not None else None
        location = wf.get("locationCode") if wf is not None else None
        channel = wf.get("channelCode") if wf is not None else None

        if not station:
            # Station code is essential; skip broken entries
            continue

        # Phase and metadata
        phase_elem = pick_elem.find(".//{*}phaseHint")
        onset_elem = pick_elem.find(".//{*}onset")
        pol_elem = pick_elem.find(".//{*}polarity")
        eval_mode = pick_elem.find(".//{*}evaluationMode")
        eval_status = pick_elem.find(".//{*}evaluationStatus")

        phase = phase_elem.text.strip() if phase_elem is not None else "?"
        onset = onset_elem.text.strip() if onset_elem is not None else None
        polarity = pol_elem.text.strip() if pol_elem is not None else None
        picker = eval_mode.text.strip() if eval_mode is not None else None

        extra: Dict[str, object] = {}
        if eval_status is not None and eval_status.text:
            extra["evaluation_status"] = eval_status.text.strip()

        picks.append(
            PhasePick(
                station=station,
                phase=phase,
                time=pick_time,
                network=network,
                location=location,
                channel=channel,
                onset=onset,
                polarity=polarity,
                weight=None,
                picker=picker,
                extra=extra,
            )
        )

    return picks


# ---------------------------------------------------------------------------
# OBS (.obs)
# ---------------------------------------------------------------------------


def read_obs_picks(path: str | Path) -> List[PhasePick]:
    """
    Read phase picks from an Antelope/OGS-style .obs file.

    The format in use is a whitespace-delimited table of the form::

        STATION ? ? ? PHASE ? YYYYMMDD HHMM SEC PICKER ERR R1 R2 R3

    Only the fields relevant to phase picking (station, phase, time,
    picker, error) are interpreted. All times are interpreted as UTC.

    Parameters
    ----------
    path : str or Path
        Path to the .obs file.

    Returns
    -------
    List[PhasePick]
        List of picks found in the file.
    """
    path = Path(path)
    picks: List[PhasePick] = []

    with path.open("r", encoding="ascii", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 9:
                # Unexpected line; skip
                continue

            station = parts[0]
            phase = parts[4]
            date_str = parts[6]
            hm_str = parts[7]
            sec_str = parts[8]

            try:
                year = int(date_str[0:4])
                month = int(date_str[4:6])
                day = int(date_str[6:8])
                hour = int(hm_str[0:2])
                minute = int(hm_str[2:4])
                seconds = float(sec_str)
                pick_time = _build_date(
                    year, month, day, hour, minute, seconds
                )
            except Exception:
                # If the time cannot be parsed, skip this row
                continue

            picker = parts[9] if len(parts) > 9 else None
            try:
                err = float(parts[10]) if len(parts) > 10 else None
            except ValueError:
                err = None

            extra: Dict[str, object] = {}
            if err is not None:
                extra["std_error"] = err

            picks.append(
                PhasePick(
                    station=station,
                    phase=phase,
                    time=pick_time,
                    picker=picker,
                    extra=extra,
                )
            )

    return picks


# ---------------------------------------------------------------------------
# DAT (.dat)
# ---------------------------------------------------------------------------


_S_PHASE_RE = re.compile(r"(\d{4})([A-Za-z\?])S\s*([0-4])")


def read_dat_picks(path: str | Path) -> List[PhasePick]:
    """
    Read phase picks from an OGS .dat file.

    Each line encodes station information and one or two phases in a
    compact fixed-width format. The relevant layout is::

        0-3   : station code
        4     : onset code (e, i, ...)
        5     : phase letter ('P' for primary; 'S' only used in S field)
        6     : quality / polarity code
        7     : weight digit
        8-18  : date-time code (see _parse_pun_datetime_code)
        19    : space
        20-23 : P seconds * 100 (e.g. 3970 -> 39.70 s), or spaces
        24-...: optional S block "SSSSeS w"

    Both P and S picks are returned when present. Absolute times are
    reconstructed from the embedded date-time code and the seconds
    fields.

    Parameters
    ----------
    path : str or Path
        Path to the .dat file.

    Returns
    -------
    List[PhasePick]
        List of picks (P and S) extracted from the file.
    """
    path = Path(path)
    picks: List[PhasePick] = []

    with path.open("r", encoding="ascii", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue

            # Ensure the line is long enough for the fixed fields
            padded = line.rstrip("\n")
            if len(padded) < 24:
                continue

            station = padded[0:4].strip()
            if not station:
                continue

            onset = padded[4].strip() or None
            # char 5 is always 'P' for the primary entry
            quality_code = padded[6].strip() or None
            weight_char = padded[7].strip() or None
            weight = (
                float(weight_char)
                if weight_char and weight_char.isdigit()
                else None
            )

            dt_code = padded[8:19]
            try:
                year, month, day, hour, minute = _parse_pun_datetime_code(
                    dt_code
                )
            except ValueError:
                continue

            # P pick
            p_sec_field = padded[20:24].strip()
            if p_sec_field:
                try:
                    seconds = int(p_sec_field) / 100.0
                    p_time = _build_date(
                        year, month, day, hour, minute, seconds
                    )
                    picks.append(
                        PhasePick(
                            station=station,
                            phase="P",
                            time=p_time,
                            onset=onset,
                            polarity=None,
                            weight=weight,
                            picker=None,
                            extra={
                                "quality_code": quality_code,
                                "source": "DAT",
                            },
                        )
                    )
                except Exception:
                    # Ignore malformed numeric fields
                    pass

            # Optional S pick encoded in the tail, pattern "SSSSeS w"
            tail = padded[24:]
            m = _S_PHASE_RE.search(tail)
            if m:
                s_code, s_onset, s_weight = m.groups()
                try:
                    s_seconds = int(s_code) / 100.0
                    s_time = _build_date(
                        year, month, day, hour, minute, s_seconds
                    )
                    picks.append(
                        PhasePick(
                            station=station,
                            phase="S",
                            time=s_time,
                            onset=s_onset,
                            polarity=None,
                            weight=float(s_weight),
                            picker=None,
                            extra={
                                "primary_onset": onset,
                                "quality_code": quality_code,
                                "source": "DAT",
                            },
                        )
                    )
                except Exception:
                    # Ignore malformed S entries
                    pass

    return picks


# ---------------------------------------------------------------------------
# HPL (.hpl)
# ---------------------------------------------------------------------------


def read_hpl_picks(path: str | Path) -> Tuple[EventOrigin, List[PhasePick]]:
    """
    Read event origin and phase picks from an OGS .hpl file.

    The HPL file contains:
      * a header line with the event location name,
      * a second (often empty) line,
      * a third line with origin date/time and hypocenter,
      * subsequent lines with per-station picks.

    This function parses the origin and the primary pick on each row.
    Some rows also contain a secondary S pick in later columns.

    Parameters
    ----------
    path : str or Path
        Path to the .hpl file.

    Returns
    -------
    (EventOrigin, List[PhasePick])
        The origin information and list of picks.

    Notes
    -----
    The depth value stored in the file is given in kilometers. In the
    returned :class:`EventOrigin` the same depth is represented by a
    negative elevation in meters:

    ``location.elevation = - depth_km * 1000.0``

    The original depth in km is preserved in ``origin.extra`` under
    the key ``"depth_km_original"`` and is also accessible via
    ``origin.depth_km``.
    """
    path = Path(path)
    with path.open("r", encoding="ascii", errors="replace") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError("HPL file is too short to contain a header.")

    # Header line 2: origin info
    header = lines[2].strip()
    header_parts = header.split()
    if len(header_parts) < 7:
        raise ValueError("Cannot parse HPL header line.")

    date_token = header_parts[0]
    time_hm = header_parts[1]
    time_sec = header_parts[2]
    lat_str = header_parts[3]
    lon_str = header_parts[4]
    depth_km = float(header_parts[5])
    magnitude = float(header_parts[6])

    year, month, day = _parse_hpl_date_token(date_token)
    hour = int(time_hm[0:2])
    minute = int(time_hm[2:4])
    seconds = float(time_sec)

    origin_time = _build_date(year, month, day, hour, minute, seconds)
    latitude = _parse_deg_min(lat_str)
    longitude = _parse_deg_min(lon_str)

    # Epicentral location with elevation representing depth (negative)
    elevation = -depth_km * 1000.0
    location = WgsPoint(
        longitude=longitude,
        latitude=latitude,
        elevation=elevation,
    )

    origin_extra: Dict[str, object] = {
        "depth_km_original": depth_km,
        "source": "HPL",
    }

    origin = EventOrigin(
        time=origin_time,
        location=location,
        magnitude=magnitude,
        extra=origin_extra,
    )

    picks: List[PhasePick] = []

    # Station picks start at line index 3
    for line in lines[3:]:
        line = line.rstrip("\n")
        if not line.strip():
            continue

        parts = line.split()
        # Expected layout (simplified):
        #  0: branch index
        #  1: station
        #  2: distance (km)
        #  3: azimuth (deg)
        #  4: takeoff (deg)
        #  5: phase code (e.g. "eP+1" or "eP 2" with blank polarity)
        #  6: HHMM
        #  7: SS.ss
        if len(parts) < 8:
            continue

        station = parts[1]
        dist_km = float(parts[2])
        azimuth = float(parts[3])
        takeoff = float(parts[4])
        phase_code = parts[5]
        time_index = 6
        if len(parts[6]) == 1 and parts[6].isdigit():
            if len(parts) < 9:
                continue
            phase_code += " " + parts[6]
            time_index += 1
        hhmm = parts[time_index]
        sec_str = parts[time_index + 1]

        if phase_code[0] in "PS":
            onset = None
        else:
            onset, phase_code = phase_code[0], phase_code[1:]
        phase = phase_code[:1] or "?"
        polarity = phase_code[1:2].strip() or None
        weight_code = phase_code[2:3]
        weight = float(weight_code) if weight_code.isdigit() else None

        hour = int(hhmm[0:2])
        minute = int(hhmm[2:4])
        seconds = float(sec_str)

        pick_time = _build_date(year, month, day, hour, minute, seconds)

        extra = {
            "distance_km": dist_km,
            "azimuth_deg": azimuth,
            "takeoff_deg": takeoff,
            "source": "HPL",
        }

        picks.append(
            PhasePick(
                station=station,
                phase=phase,
                time=pick_time,
                onset=onset,
                polarity=polarity,
                weight=weight,
                picker=None,
                extra=extra,
            )
        )

    return origin, picks


# ---------------------------------------------------------------------------
# PUN (.pun)
# ---------------------------------------------------------------------------


def read_pun_origin(path: str | Path) -> EventOrigin:
    """
    Read event origin information from an OGS .pun file.

    The .pun file contains two main lines:

      1. A header line with column labels.
      2. A line with encoded date/time, coordinates, depth, magnitude,
         and additional quality parameters.

    Only the core origin information is extracted here.

    Parameters
    ----------
    path : str or Path
        Path to the .pun file.

    Returns
    -------
    EventOrigin
        Parsed origin information.

    Notes
    -----
    The depth value stored in the file is given in kilometers. In the
    returned :class:`EventOrigin` the same depth is represented by a
    negative elevation in meters:

    ``location.elevation = - depth_km * 1000.0``

    The original depth in km is preserved in ``origin.extra`` under
    the key ``"depth_km_original"`` and is also accessible via
    ``origin.depth_km``.
    """
    path = Path(path)
    with path.open("r", encoding="ascii", errors="replace") as f:
        lines = [ln.strip() for ln in f.readlines()]

    if len(lines) < 2:
        raise ValueError("PUN file is too short to contain origin data.")

    origin_line = lines[1]
    parts = origin_line.split()
    if len(parts) < 7:
        raise ValueError("Cannot parse PUN origin line.")

    datetime_code = parts[0]
    sec_str = parts[1]
    lat_str = parts[2]
    lon_str = parts[3]
    depth_km = float(parts[4])
    magnitude = float(parts[5])

    year, month, day, hour, minute = _parse_pun_datetime_code(
        datetime_code
    )
    seconds = float(sec_str)
    origin_time = _build_date(year, month, day, hour, minute, seconds)
    latitude = _parse_deg_min(lat_str)
    longitude = _parse_deg_min(lon_str)

    elevation = -depth_km * 1000.0
    location = WgsPoint(
        longitude=longitude,
        latitude=latitude,
        elevation=elevation,
    )

    # Collect a few additional parameters when present
    extra: Dict[str, object] = {"depth_km_original": depth_km, "source": "PUN"}
    if len(parts) >= 7:
        try:
            extra["gap_deg"] = float(parts[6])
        except ValueError:
            pass
    if len(parts) >= 8:
        try:
            extra["num_sta"] = int(parts[7])
        except ValueError:
            pass

    return EventOrigin(
        time=origin_time,
        location=location,
        magnitude=magnitude,
        extra=extra,
    )


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def load_picks(path: str | Path) -> List[PhasePick]:
    """
    Dispatch to the appropriate reader based on file extension.

    Parameters
    ----------
    path : str or Path
        File path. The extension determines which reader is used.

    Returns
    -------
    List[PhasePick]
        List of picks read from the file.

    Raises
    ------
    ValueError
        If the extension is unknown or unsupported.
    """
    path = Path(path)
    ext = path.suffix.lower()

    if ext == ".xml":
        return read_quakeml_picks(path)
    if ext == ".obs":
        return read_obs_picks(path)
    if ext == ".dat":
        return read_dat_picks(path)
    if ext == ".hpl":
        # For convenience, discard the origin here
        _, picks = read_hpl_picks(path)
        return picks

    raise ValueError(f"Unsupported pick file extension: {ext!r}")


def normalize_phases(picks: List[PhasePick]) -> List[PhasePick]:
    """
    Normalize detailed phase names (Pg, Pn, Sg, Sn, etc.) to "P" or "S".

    This helper is intentionally simple:
        * any phase starting with "P" or "p" becomes "P"
        * any phase starting with "S" or "s" becomes "S"
        * all other phases are left unchanged (or could optionally be removed)

    Parameters
    ----------
    picks : List[PhasePick]
        List of picks to normalize (in place).

    Returns
    -------
    List[PhasePick]
        The same list, for convenience.

    Notes
    -----
    - The function modifies the PhasePick instances in place.
    - This is intended as a lightweight compatibility layer for engines
      that only support fundamental phase families ("P" and "S").
    """
    for p in picks:
        if not p.phase:
            continue
        ph = p.phase.strip().upper()
        if ph.startswith("P"):
            p.phase = "P"
        elif ph.startswith("S"):
            p.phase = "S"
        # else: leave unchanged (e.g. Lg, Rg, IAML, ...)

    return picks
