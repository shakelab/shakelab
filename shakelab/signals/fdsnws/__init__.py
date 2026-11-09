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
Client of the FDSN web services.

The package has four layers:

- :mod:`~shakelab.signals.fdsnws.transport`: HTTP requests (standard
  library only), retries, errors and parameter formatting;
- :mod:`~shakelab.signals.fdsnws.registry`: the data centers and the
  addresses of their services;
- :mod:`~shakelab.signals.fdsnws.selection`: channel patterns
  translated into requests;
- :mod:`~shakelab.signals.fdsnws.client`: the :class:`Client`, with
  waveforms, metadata, events and data availability; the text format
  of the station service is read by
  :mod:`~shakelab.signals.fdsnws.stationtext`.

Examples
--------
>>> from shakelab.signals.fdsnws import Client
>>> client = Client("INGV")
>>> sc = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:03:23",
...                           "2012-05-20T02:06:23", metadata=True)
>>> sc.deconvolve_response()
"""

from __future__ import annotations

from shakelab.signals.fdsnws.client import Client, TimeSpan
from shakelab.signals.fdsnws.registry import DataCenter, datacenters
from shakelab.signals.fdsnws.stationtext import read_station_text
from shakelab.signals.fdsnws.transport import FdsnError, FdsnWarning

__all__ = [
    "Client",
    "DataCenter",
    "FdsnError",
    "FdsnWarning",
    "TimeSpan",
    "datacenters",
    "read_station_text",
]
