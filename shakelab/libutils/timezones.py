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
Time zone table and UTC offset helpers.

:data:`TIMEZONE` lists common time zones, keyed by short codes derived
from the Windows time zone names (e.g. ``"WET"`` for *W. Europe
Standard Time*, ``"JST"`` for *Jordan Standard Time*). Each entry has

- ``label``: display label with the standard offset and main cities;
- ``delay``: standard (non-daylight-saving) offset from UTC in hours;
- ``iana``: corresponding IANA time zone name (e.g. ``"Europe/Berlin"``).

The helpers convert a zone into an offset, optionally at a given date
(daylight saving time included), and into the ``"+HH:MM"`` strings
accepted by :meth:`shakelab.libutils.timeN.Date.to_iso8601`.

Conventions
-----------
- The codes are not the usual civil abbreviations: several of them
  (``JST``, ``CEST``, ``IST``, ``CST``, ``AST``, ``BST``, ``GST``)
  denote different zones elsewhere. Prefer the IANA names when
  exchanging data.
- Daylight saving time is computed with the standard-library
  :mod:`zoneinfo` module, which uses the system time zone database. If
  the database is not available (e.g. on Windows without the optional
  ``tzdata`` package), table zones fall back to their standard offset
  with a warning; no extra dependency is required.
- Without a date, a table zone returns its standard offset and an IANA
  name returns its current offset.
"""

from __future__ import annotations

import warnings
from datetime import datetime, timedelta, timezone
from typing import Any, Optional
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from shakelab.libutils.timeN import Date

__all__ = [
    "TIMEZONE",
    "utc_offset",
    "offset_string",
    "format_offset",
    "find_timezone",
]


# ----------------------------------------------------------------------------
# Offset formatting
# ----------------------------------------------------------------------------


def format_offset(hours: float) -> str:
    """
    Format an offset in hours as an ISO 8601 ``"+HH:MM"`` string.

    Parameters
    ----------
    hours : float
        Offset from UTC in hours (e.g. ``5.75`` or ``-3.5``); it is
        rounded to the minute.

    Returns
    -------
    str
        Offset such as ``"+05:45"`` or ``"-03:30"``.

    Raises
    ------
    ValueError
        If the offset is not strictly within one day.
    """
    minutes = round(float(hours) * 60.0)
    if abs(minutes) >= 24 * 60:
        raise ValueError(f"offset out of range: {hours!r} hours")

    sign = "-" if minutes < 0 else "+"
    hh, mm = divmod(abs(minutes), 60)
    return f"{sign}{hh:02d}:{mm:02d}"


def _label(delay: float, cities: str) -> str:
    if delay == 0.0:
        return f"(UTC) {cities}"
    return f"(UTC{format_offset(delay)}) {cities}"


# ----------------------------------------------------------------------------
# Time zone table
# ----------------------------------------------------------------------------

# (code, standard offset in hours, IANA name, cities)
_ZONES = (
    ("ABST", 3.0, "Asia/Riyadh", "Kuwait, Riyadh"),
    ("ACST", 9.5, "Australia/Darwin", "Darwin"),
    ("AEST", 10.0, "Australia/Sydney", "Canberra, Melbourne, Sydney"),
    ("AFT", 4.5, "Asia/Kabul", "Kabul"),
    ("AK", -9.0, "America/Anchorage", "Alaska"),
    ("AMT", 4.0, "Asia/Yerevan", "Yerevan"),
    ("ARBST", 4.0, "Asia/Dubai", "Abu Dhabi, Muscat"),
    ("ARST", 3.0, "Asia/Baghdad", "Baghdad"),
    ("ART", -3.0, "America/Argentina/Buenos_Aires", "Buenos Aires"),
    ("AST", -4.0, "America/Halifax", "Atlantic Time (Canada)"),
    ("AWST", 8.0, "Australia/Perth", "Perth"),
    ("AZOST", -1.0, "Atlantic/Azores", "Azores"),
    ("AZT", 4.0, "Asia/Baku", "Baku"),
    ("BIT", -12.0, "Etc/GMT+12", "International Date Line West"),
    ("BST", 6.0, "Asia/Dhaka", "Dhaka"),
    ("BTT", 5.0, "Asia/Almaty", "Astana"),
    ("CAST", -6.0, "America/Guatemala", "Central America"),
    ("CAUST", 9.5, "Australia/Adelaide", "Adelaide"),
    ("CBST", -4.0, "America/Cuiaba", "Cuiaba"),
    ("CCST", -6.0, "America/Regina", "Saskatchewan"),
    ("CEST", 1.0, "Europe/Warsaw", "Sarajevo, Skopje, Warsaw, Zagreb"),
    ("CET", 1.0, "Europe/Budapest",
     "Belgrade, Bratislava, Budapest, Ljubljana, Prague"),
    ("CHN", 8.0, "Asia/Shanghai",
     "Beijing, Chongqing, Hong Kong, Urumqi"),
    ("CST", -6.0, "America/Mexico_City",
     "Guadalajara, Mexico City, Monterrey"),
    ("CT", -6.0, "America/Chicago", "Central Time (US & Canada)"),
    ("CVT", -1.0, "Atlantic/Cape_Verde", "Cape Verde Is."),
    ("EAST", 10.0, "Australia/Brisbane", "Brisbane"),
    ("EAT", 3.0, "Africa/Nairobi", "Nairobi"),
    ("ECT", 1.0, "Africa/Lagos", "West Central Africa"),
    ("EEST", 3.0, "Europe/Minsk", "Minsk"),
    ("EET", 2.0, "Europe/Helsinki",
     "Helsinki, Kyiv, Riga, Sofia, Tallinn, Vilnius"),
    ("EGST", 2.0, "Africa/Cairo", "Cairo"),
    ("ESAST", -3.0, "America/Sao_Paulo", "Brasilia"),
    ("EST", -5.0, "America/Indiana/Indianapolis", "Indiana (East)"),
    ("ET", -5.0, "America/New_York", "Eastern Time (US & Canada)"),
    ("FJT", 12.0, "Pacific/Fiji", "Fiji"),
    ("GET", 4.0, "Asia/Tbilisi", "Tbilisi"),
    ("GMT", 0.0, "Europe/London", "Dublin, Edinburgh, Lisbon, London"),
    ("GNST", -2.0, "America/Nuuk", "Greenland"),
    ("GST", 0.0, "Atlantic/Reykjavik", "Monrovia, Reykjavik"),
    ("GTBST", 2.0, "Europe/Bucharest", "Athens, Bucharest"),
    ("HAST", -10.0, "Pacific/Honolulu", "Hawaii"),
    ("IRKT", 8.0, "Asia/Irkutsk", "Irkutsk"),
    ("IRST", 3.5, "Asia/Tehran", "Tehran"),
    ("ISST", 2.0, "Asia/Jerusalem", "Jerusalem"),
    ("IST", 5.5, "Asia/Kolkata",
     "Chennai, Kolkata, Mumbai, New Delhi"),
    ("JST", 3.0, "Asia/Amman", "Amman"),
    ("KRAT", 7.0, "Asia/Krasnoyarsk", "Krasnoyarsk"),
    ("KST", 9.0, "Asia/Seoul", "Seoul"),
    ("MEST", 2.0, "Asia/Beirut", "Beirut"),
    ("MOST", 1.0, "Africa/Casablanca", "Casablanca"),
    ("MSK", 3.0, "Europe/Moscow", "Moscow, St. Petersburg, Volgograd"),
    ("MST", -7.0, "America/Phoenix", "Arizona"),
    ("MSTM", -7.0, "America/Mazatlan", "La Paz, Mazatlan"),
    ("MT", -7.0, "America/Denver", "Mountain Time (US & Canada)"),
    ("MUT", 4.0, "Indian/Mauritius", "Port Louis"),
    ("MVST", -3.0, "America/Montevideo", "Montevideo"),
    ("MYST", 6.5, "Asia/Yangon", "Yangon (Rangoon)"),
    ("NCAST", 7.0, "Asia/Novosibirsk", "Novosibirsk"),
    ("NMST", 2.0, "Africa/Windhoek", "Windhoek"),
    ("NPT", 5.75, "Asia/Kathmandu", "Kathmandu"),
    ("NST", -3.5, "America/St_Johns", "Newfoundland"),
    ("NZST", 12.0, "Pacific/Auckland", "Auckland, Wellington"),
    ("PETT", 12.0, "Asia/Kamchatka", "Petropavlovsk-Kamchatsky"),
    ("PHOT", 13.0, "Pacific/Tongatapu", "Nuku'alofa"),
    ("PKT", 5.0, "Asia/Karachi", "Islamabad, Karachi"),
    ("PRST", -3.0, "America/Asuncion", "Asuncion"),
    ("PSAST", -4.0, "America/Santiago", "Santiago"),
    ("PST", -8.0, "America/Tijuana", "Baja California"),
    ("PT", -8.0, "America/Los_Angeles", "Pacific Time (US & Canada)"),
    ("RST", 1.0, "Europe/Paris", "Brussels, Copenhagen, Madrid, Paris"),
    ("SAEST", -3.0, "America/Cayenne", "Cayenne, Fortaleza"),
    ("SAPST", -5.0, "America/Bogota", "Bogota, Lima, Quito"),
    ("SAST", 2.0, "Africa/Johannesburg", "Harare, Pretoria"),
    ("SAWST", -4.0, "America/La_Paz",
     "Georgetown, La Paz, Manaus, San Juan"),
    ("SBT", 11.0, "Pacific/Guadalcanal",
     "Magadan, Solomon Is., New Caledonia"),
    ("SLT", 5.5, "Asia/Colombo", "Sri Jayawardenepura"),
    ("SMST", 13.0, "Pacific/Apia", "Samoa"),
    ("SNST", 8.0, "Asia/Singapore", "Kuala Lumpur, Singapore"),
    ("SST", 3.0, "Asia/Damascus", "Damascus"),
    ("TAST", 10.0, "Australia/Hobart", "Hobart"),
    ("THA", 7.0, "Asia/Bangkok", "Bangkok, Hanoi, Jakarta"),
    ("TIST", 8.0, "Asia/Taipei", "Taipei"),
    ("TST", 9.0, "Asia/Tokyo", "Osaka, Sapporo, Tokyo"),
    ("UST", 8.0, "Asia/Ulaanbaatar", "Ulaanbaatar"),
    ("UTC", 0.0, "UTC", "Coordinated Universal Time"),
    ("UTC12", 12.0, "Etc/GMT-12", "Coordinated Universal Time+12"),
    ("VLAT", 10.0, "Asia/Vladivostok", "Vladivostok"),
    ("VST", -4.0, "America/Caracas", "Caracas"),
    ("WAST", 5.0, "Asia/Tashkent", "Tashkent"),
    ("WET", 1.0, "Europe/Berlin",
     "Amsterdam, Berlin, Bern, Rome, Stockholm, Vienna"),
    ("WPST", 10.0, "Pacific/Port_Moresby", "Guam, Port Moresby"),
    ("YAKT", 9.0, "Asia/Yakutsk", "Yakutsk"),
    ("YEKT", 5.0, "Asia/Yekaterinburg", "Ekaterinburg"),
)

TIMEZONE: dict[str, dict[str, Any]] = {
    code: {"label": _label(delay, cities), "delay": delay, "iana": iana}
    for code, delay, iana, cities in _ZONES
}
"""Time zones by code: ``label``, standard ``delay`` (h) and ``iana``."""


# ----------------------------------------------------------------------------
# Offsets
# ----------------------------------------------------------------------------


def utc_offset(zone: str, date: Optional[Any] = None) -> float:
    """
    Return the offset of a time zone from UTC.

    Parameters
    ----------
    zone : str
        Code of :data:`TIMEZONE` (e.g. ``"WET"``) or IANA name (e.g.
        ``"Europe/Rome"``).
    date : Date or datetime.datetime, optional
        Instant at which the offset is evaluated, as a
        :class:`~shakelab.libutils.timeN.Date` or a
        :class:`datetime.datetime` (naive values are taken as UTC).
        If None, a table zone returns its standard offset and an IANA
        name its current offset.

    Returns
    -------
    float
        Offset in hours, daylight saving time included when a date is
        given.

    Raises
    ------
    ValueError
        If the zone is unknown.
    TypeError
        If ``zone`` is not a string or ``date`` has an unsupported type.

    Warns
    -----
    RuntimeWarning
        If the time zone database is not available and the standard
        offset of a table zone is returned instead.
    """
    if not isinstance(zone, str):
        raise TypeError("zone must be a string")

    entry = TIMEZONE.get(zone)
    if entry is not None and date is None:
        return entry["delay"]

    when = _to_utc_datetime(date)
    name = entry["iana"] if entry is not None else zone

    try:
        tzinfo = ZoneInfo(name)
    except (ZoneInfoNotFoundError, ValueError):
        if entry is None:
            raise ValueError(f"unknown time zone: {zone!r}") from None
        warnings.warn(
            f"time zone data for {name!r} not available; "
            f"using the standard offset of {zone!r}",
            RuntimeWarning,
            stacklevel=2,
        )
        return entry["delay"]

    offset = when.astimezone(tzinfo).utcoffset()
    return offset.total_seconds() / 3600.0


def offset_string(zone: str, date: Optional[Any] = None) -> str:
    """
    Return the offset of a time zone as a ``"+HH:MM"`` string.

    The result can be passed as ``timezone`` to
    :meth:`shakelab.libutils.timeN.Date.to_iso8601`.

    Parameters
    ----------
    zone : str
        Code of :data:`TIMEZONE` or IANA name.
    date : Date or datetime.datetime, optional
        Instant at which the offset is evaluated (see
        :func:`utc_offset`).

    Returns
    -------
    str
        Offset such as ``"+02:00"``.

    Examples
    --------
    >>> date = Date("2026-07-01T10:00:00Z")
    >>> offset = offset_string("Europe/Rome", date)
    >>> date.to_iso8601(precision=0, timezone=offset)
    '2026-07-01T12:00:00+02:00'
    """
    return format_offset(utc_offset(zone, date))


def _to_utc_datetime(date: Optional[Any]) -> datetime:
    """Convert a Date, a datetime or None (now) to an aware UTC datetime."""
    if date is None:
        return datetime.now(timezone.utc)

    if isinstance(date, Date):
        epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
        return epoch + timedelta(microseconds=date.to_unix_ns() // 1000)

    if isinstance(date, datetime):
        if date.tzinfo is None:
            return date.replace(tzinfo=timezone.utc)
        return date.astimezone(timezone.utc)

    raise TypeError("date must be a Date, a datetime or None")


# ----------------------------------------------------------------------------
# Search
# ----------------------------------------------------------------------------


def find_timezone(text: str) -> list[str]:
    """
    Find table zones by code, city, label or IANA name.

    Parameters
    ----------
    text : str
        Case-insensitive text to search for (e.g. ``"rome"``).

    Returns
    -------
    list of str
        Matching codes, in table order.
    """
    needle = text.strip().lower()
    if not needle:
        return []

    return [
        code for code, entry in TIMEZONE.items()
        if needle in code.lower()
        or needle in entry["label"].lower()
        or needle in entry["iana"].lower()
    ]
