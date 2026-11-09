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
Channel selections of FDSN web service requests.

Channels are selected with the patterns of :mod:`shakelab.signals.codes`
(``NET.STA.LOC.CHA`` with ``*``, ``?``, ``[...]`` and ``{A,B}``, codes
missing at the end matching anything). The FDSN services know only
``*``, ``?`` and lists separated by commas, so each pattern is
translated into codes the server understands:

- an empty location becomes ``--``;
- a set ``[...]`` becomes the list of its words (``HH[ZNE]`` →
  ``HHZ,HHN,HHE``), or ``?`` if it is negated or too large;
- the alternatives ``{A,B}`` become separate requests.

The server may return more than asked (a ``?`` in place of a set):
the client keeps only the channels that match the original patterns.

One request is sent with GET; several, or requests with different time
windows, are sent with POST, one line per request.
"""

from __future__ import annotations

import itertools
from typing import Any, NamedTuple, Optional

from shakelab.libutils.timeN import Date
from shakelab.signals.codes import SidPattern, parse_sid_pattern
from shakelab.signals.fdsnws.transport import format_value

__all__ = [
    "Request",
    "get_params",
    "parse_selection",
    "post_body",
    "server_codes",
    "to_date",
]

# Largest number of words a set ``[...]`` is expanded into.
MAX_SET_WORDS = 100


class Request(NamedTuple):
    """
    One channel pattern and its time window.

    Attributes
    ----------
    pattern : SidPattern
        Channel pattern, without braces.
    starttime, endtime : Date or None
        Time window; None if not given.
    """

    pattern: SidPattern
    starttime: Optional[Date]
    endtime: Optional[Date]

    def matches(self, sid: Optional[str]) -> bool:
        """
        Return whether a source identifier matches the pattern.

        Parameters
        ----------
        sid : str or None
            ``NET.STA.LOC.CHA``.

        Returns
        -------
        bool
            True if the identifier matches.
        """
        return self.pattern.matches(sid)


# ----------------------------------------------------------------------------
# Selections
# ----------------------------------------------------------------------------


def parse_selection(
    select: Any,
    starttime: Any = None,
    endtime: Any = None,
) -> list:
    """
    Return the requests of a channel selection.

    Parameters
    ----------
    select : str or sequence
        A pattern, a sequence of patterns, or a list whose items are
        patterns or ``(pattern, starttime, endtime)`` tuples with their
        own time windows (None takes the window of the arguments).
    starttime, endtime : Date or str, optional
        Time window of the patterns without one (ISO 8601 strings or
        :class:`Date`).

    Returns
    -------
    list of Request
        One request per pattern and brace alternative, in order,
        without repetitions.

    Raises
    ------
    TypeError
        If the selection or a time is of a wrong type.
    ValueError
        If a pattern is invalid, a time cannot be parsed, or a time
        window ends before it starts.

    Examples
    --------
    >>> [str(r.pattern) for r in parse_selection("IV.{MODE,CAVE}..HN?")]
    ['IV.MODE..HN?', 'IV.CAVE..HN?']
    """
    start = to_date(starttime, "starttime")
    end = to_date(endtime, "endtime")

    if isinstance(select, str):
        items = [select]
    elif isinstance(select, (list, tuple)):
        items = list(select)
    else:
        raise TypeError(
            "select must be a pattern or a sequence of patterns"
        )
    if not items:
        raise ValueError("empty selection")

    requests = []
    for item in items:
        if isinstance(item, str):
            pattern, item_start, item_end = item, start, end
        elif isinstance(item, (list, tuple)) and len(item) == 3:
            pattern = item[0]
            item_start = (
                start if item[1] is None else to_date(item[1], "starttime")
            )
            item_end = end if item[2] is None else to_date(item[2], "endtime")
        else:
            raise TypeError(
                "selection items must be patterns or "
                "(pattern, starttime, endtime) tuples"
            )
        if (
            item_start is not None and item_end is not None
            and item_end < item_start
        ):
            raise ValueError(f"{pattern}: endtime before starttime")
        for parsed in parse_sid_pattern(pattern):
            request = Request(parsed, item_start, item_end)
            if request not in requests:
                requests.append(request)
    return requests


def server_codes(pattern: SidPattern) -> tuple:
    """
    Return the codes of a pattern as the FDSN services expect them.

    Parameters
    ----------
    pattern : SidPattern
        Channel pattern without braces.

    Returns
    -------
    tuple of str
        Network, station, location and channel codes, with lists
        separated by commas and ``--`` for the empty location.

    Examples
    --------
    >>> server_codes(parse_sid_pattern("IV.MODE..HN[ZNE]")[0])
    ('IV', 'MODE', '--', 'HNZ,HNN,HNE')
    """
    codes = []
    for index, code in enumerate(pattern):
        if index == 2 and code == "":
            codes.append("--")
            continue
        words = _expand_sets(code)
        codes.append(",".join(words))
    return tuple(codes)


def get_params(requests: list) -> Optional[dict]:
    """
    Return the GET parameters of the requests, if one line suffices.

    Parameters
    ----------
    requests : list of Request
        Requests of a selection.

    Returns
    -------
    dict or None
        ``network``, ``station``, ``location``, ``channel``,
        ``starttime`` and ``endtime`` (None values for missing times);
        None if the requests need POST (more than one request).
    """
    if len(requests) != 1:
        return None
    request = requests[0]
    network, station, location, channel = server_codes(request.pattern)
    return {
        "network": network,
        "station": station,
        "location": location,
        "channel": channel,
        "starttime": request.starttime,
        "endtime": request.endtime,
    }


def post_body(requests: list, params: Optional[dict] = None) -> str:
    """
    Return the body of a POST request.

    Parameters
    ----------
    requests : list of Request
        Requests of a selection.
    params : dict, optional
        Other parameters, written first as ``key=value`` lines (None
        values are left out).

    Returns
    -------
    str
        The parameter lines, then one line ``NET STA LOC CHA START
        END`` for each combination of the listed codes; a missing time
        is written as ``*``.
    """
    lines = []
    for key, value in (params or {}).items():
        if value is not None:
            lines.append(f"{key}={format_value(value)}")

    for request in requests:
        start = (
            "*" if request.starttime is None
            else format_value(request.starttime)
        )
        end = (
            "*" if request.endtime is None
            else format_value(request.endtime)
        )
        lists = [code.split(",") for code in server_codes(request.pattern)]
        for codes in itertools.product(*lists):
            line = f"{' '.join(codes)} {start} {end}"
            if line not in lines:
                lines.append(line)
    return "\n".join(lines) + "\n"


def to_date(value: Any, name: str = "time") -> Optional[Date]:
    """
    Return a time as a Date.

    Parameters
    ----------
    value : Date, str or None
        :class:`Date` or ISO 8601 string; None is kept.
    name : str, optional
        Name used in the error messages.

    Returns
    -------
    Date or None
        The time.

    Raises
    ------
    TypeError
        If the value is neither a Date nor a string.
    ValueError
        If the string cannot be parsed.
    """
    if value is None or isinstance(value, Date):
        return value
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a Date or an ISO 8601 string")
    return Date(value)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _expand_sets(code: str) -> list:
    """Words of a code with its sets ``[...]`` expanded, or ``?``."""
    parts = []
    index = 0
    while index < len(code):
        char = code[index]
        if char != "[":
            parts.append([char])
            index += 1
            continue
        close = code.find("]", index + 2)
        if close < 0:
            parts.append([char])
            index += 1
            continue
        parts.append(_set_chars(code[index + 1:close]))
        index = close + 1

    count = 1
    for options in parts:
        count *= len(options)
    if count > MAX_SET_WORDS:
        parts = [options if len(options) == 1 else ["?"]
                 for options in parts]
    return ["".join(word) for word in itertools.product(*parts)]


def _set_chars(content: str) -> list:
    """Characters of a set; ``["?"]`` for a negated set."""
    if content.startswith(("!", "^")):
        return ["?"]
    chars: list = []
    index = 0
    while index < len(content):
        if (
            index + 2 < len(content) and content[index + 1] == "-"
        ):
            first, last = ord(content[index]), ord(content[index + 2])
            for code in range(first, last + 1):
                if chr(code) not in chars:
                    chars.append(chr(code))
            index += 3
        else:
            if content[index] not in chars:
                chars.append(content[index])
            index += 1
    return chars or ["?"]
