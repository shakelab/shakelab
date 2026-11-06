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
Source identifier codes and patterns.

A source identifier (sid) is ``NET.STA.LOC.CHA``: network, station,
location and channel codes; the location may be empty. A station code
is ``NET.STA``.

Patterns select identifiers code by code, so that a joker never spans
a dot:

- ``*`` matches any sequence of characters, ``?`` one character and
  ``[...]`` one of a set of characters (``HH[ZNE]``);
- ``{A,B,...}`` gives alternative words (``OX.{CAE,BAD}.*.HHZ``), also
  inside a code (``HH{Z,N}``);
- codes missing at the end match anything: ``"IV.MODE"`` is every
  channel of MODE, ``"IV"`` every channel of the network;
- an empty code matches the empty location only (``"IV.MODE..HNZ"``).

The same syntax is used by :meth:`MetadataCollection.select
<shakelab.signals.metadata.MetadataCollection.select>`, by the SeedLink
client and by other selections of channels.
"""

from __future__ import annotations

import fnmatch
from typing import Any, NamedTuple, Optional

__all__ = [
    "SidPattern",
    "expand_braces",
    "match_code",
    "match_sid",
    "parse_sid_pattern",
]


# ----------------------------------------------------------------------------
# Patterns
# ----------------------------------------------------------------------------


class SidPattern(NamedTuple):
    """
    Pattern of a source identifier, one joker pattern per code.

    Attributes
    ----------
    network, station, location, channel : str
        Patterns of the codes (``*``, ``?``, ``[...]``).
    """

    network: str
    station: str
    location: str
    channel: str

    def __str__(self) -> str:
        return ".".join(self)

    def matches(self, sid: Optional[str]) -> bool:
        """
        Return whether a source identifier matches the pattern.

        Parameters
        ----------
        sid : str or None
            ``NET.STA.LOC.CHA``; None and malformed identifiers never
            match.

        Returns
        -------
        bool
            True if the identifier matches.
        """
        if not isinstance(sid, str):
            return False
        codes = sid.split(".")
        if len(codes) != 4:
            return False
        return all(
            fnmatch.fnmatchcase(code, pattern)
            for code, pattern in zip(codes, self)
        )

    @property
    def has_station_joker(self) -> bool:
        """Whether the network or station pattern has jokers."""
        return _has_joker(self.network + self.station)


def expand_braces(pattern: str) -> list:
    """
    Expand the ``{A,B,...}`` alternatives of a pattern.

    Parameters
    ----------
    pattern : str
        Pattern, possibly with nested alternatives.

    Returns
    -------
    list of str
        Patterns without braces, in order, without repetitions.

    Raises
    ------
    ValueError
        If the braces are not balanced.

    Examples
    --------
    >>> expand_braces("OX.{CAE,BAD}.*.HH{Z,N}")
    ['OX.CAE.*.HHZ', 'OX.CAE.*.HHN', 'OX.BAD.*.HHZ', 'OX.BAD.*.HHN']
    """
    if not isinstance(pattern, str):
        raise TypeError("pattern must be a string")

    start = pattern.find("{")
    if start < 0:
        if "}" in pattern:
            raise ValueError(f"unbalanced braces in {pattern!r}")
        return [pattern]

    depth = 0
    parts = []
    begin = start + 1
    for index in range(start, len(pattern)):
        char = pattern[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                parts.append(pattern[begin:index])
                end = index
                break
        elif char == "," and depth == 1:
            parts.append(pattern[begin:index])
            begin = index + 1
    else:
        raise ValueError(f"unbalanced braces in {pattern!r}")

    prefix, suffix = pattern[:start], pattern[end + 1:]
    if "}" in prefix:
        raise ValueError(f"unbalanced braces in {pattern!r}")

    result = []
    for part in parts:
        for expanded in expand_braces(prefix + part + suffix):
            if expanded not in result:
                result.append(expanded)
    return result


def parse_sid_pattern(pattern: str, parts: int = 4) -> list:
    """
    Parse a source identifier pattern.

    Parameters
    ----------
    pattern : str
        Pattern (see the module notes).
    parts : int, optional
        Number of codes: 4 for channels (default), 2 for stations.

    Returns
    -------
    list of SidPattern or list of tuple
        One pattern per brace alternative: :class:`SidPattern` for
        channels, ``(network, station)`` tuples for stations.

    Raises
    ------
    TypeError
        If the pattern is not a string.
    ValueError
        If the pattern is empty, has too many codes, an empty network,
        station or channel code, spaces or unbalanced braces.
    """
    if not isinstance(pattern, str):
        raise TypeError("patterns must be strings")
    if parts not in (2, 4):
        raise ValueError("parts must be 2 or 4")

    text = pattern.strip()
    if not text:
        raise ValueError("empty pattern")

    result = []
    for expanded in expand_braces(text):
        codes = expanded.split(".")
        if len(codes) > parts:
            raise ValueError(
                f"invalid pattern {pattern!r}: more than {parts} codes"
            )
        codes += ["*"] * (parts - len(codes))
        for index, code in enumerate(codes):
            if not code and index != 2:
                raise ValueError(f"invalid pattern {pattern!r}: empty code")
            if any(char.isspace() for char in code):
                raise ValueError(f"invalid pattern {pattern!r}: spaces")
        parsed = SidPattern(*codes) if parts == 4 else tuple(codes)
        if parsed not in result:
            result.append(parsed)
    return result


def match_sid(sid: Optional[str], pattern: Any) -> bool:
    """
    Return whether a source identifier matches a pattern.

    Parameters
    ----------
    sid : str or None
        ``NET.STA.LOC.CHA``.
    pattern : str or SidPattern
        Pattern (see the module notes).

    Returns
    -------
    bool
        True if the identifier matches.

    Examples
    --------
    >>> match_sid("OX.CAE..HHZ", "OX.{CAE,BAD}")
    True
    >>> match_sid("OX.CAE..HHZ", "*Z")
    False
    """
    if isinstance(pattern, SidPattern):
        return pattern.matches(sid)
    return any(item.matches(sid) for item in parse_sid_pattern(pattern))


def match_code(code: Optional[str], pattern: str) -> bool:
    """
    Return whether a station code ``NET.STA`` matches a pattern.

    Parameters
    ----------
    code : str or None
        Station code.
    pattern : str
        Pattern of at most two codes (``"IV"``, ``"IV.MO*"``,
        ``"{IV,OX}.*"``).

    Returns
    -------
    bool
        True if the code matches.
    """
    if not isinstance(code, str):
        return False
    codes = code.split(".")
    if len(codes) != 2:
        return False
    return any(
        all(fnmatch.fnmatchcase(value, item) for value, item in zip(codes, p))
        for p in parse_sid_pattern(pattern, parts=2)
    )


def _has_joker(text: str) -> bool:
    return any(char in text for char in "*?[")
