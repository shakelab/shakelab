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
SeedLink client for real-time and recent waveform data.

The package has two layers:

- :mod:`~shakelab.signals.seedlink.protocol`: connection, negotiation
  of the protocol version (3.1 or 4), commands and packets;
- :mod:`~shakelab.signals.seedlink.archiver`: MiniSEED archive files;
- :mod:`~shakelab.signals.seedlink.client`: the :class:`Client`, which
  requests channels and distributes the records to a buffer, an archive
  and callbacks, and the :func:`read`, :func:`archive` and
  :func:`records` functions.

Examples
--------
>>> from shakelab.signals import seedlink
>>> sc = seedlink.read("rtserve.example.org", "IV.MODE.*.HH?",
...                    duration=60)
"""

from __future__ import annotations

from shakelab.signals.seedlink.archiver import Archive
from shakelab.signals.seedlink.client import (
    Client,
    SeedLinkEvent,
    archive,
    read,
    records,
)
from shakelab.signals.seedlink.protocol import DEFAULT_PORT, SeedLinkError

__all__ = [
    "DEFAULT_PORT",
    "Archive",
    "Client",
    "SeedLinkError",
    "SeedLinkEvent",
    "archive",
    "read",
    "records",
]
