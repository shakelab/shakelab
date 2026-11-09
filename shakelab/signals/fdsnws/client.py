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

A :class:`Client` talks to one data center, named in the registry
(:mod:`~shakelab.signals.fdsnws.registry`) or given by its address:

- :meth:`Client.get_waveforms` (dataselect): waveforms as a
  :class:`~shakelab.signals.base.StreamCollection`, optionally with
  their channel metadata, corrected for the instrument response or
  filtered, or the MiniSEED records written unchanged to a file;
- :meth:`Client.get_metadata` (station): channel and station epochs as
  a :class:`~shakelab.signals.metadata.MetadataCollection`;
- :meth:`Client.get_events` (event): an earthquake catalogue;
- :meth:`Client.get_availability` (availability): the time spans of
  the data held by the data center;
- :meth:`Client.version`: the version of a service.

Channels are selected with the patterns of :mod:`shakelab.signals.codes`
(see :mod:`~shakelab.signals.fdsnws.selection`); geographic constraints
follow the conventions of ShakeLab: ``region=(lonmin, lonmax, latmin,
latmax)``, ``center=(lon, lat)``, radii in degrees and depths in km.

Corrected downloads (``correct=True``) request a margin ``pad`` before
and after the window, remove a linear trend, taper the margins, correct
the response in the frequency domain with zero padding (and the
optional ``prefilter``), filter if asked, and cut the records to the
window: the edge effects fall in the discarded margins.

Examples
--------
>>> from shakelab.signals.fdsnws import Client
>>> client = Client("INGV")
>>> sc = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:03:23",
...                           "2012-05-20T02:06:23", correct=True)
>>> md = client.get_metadata("IV.MO*", level="channel")
"""

from __future__ import annotations

import json
import os
import warnings
from numbers import Real
from typing import Any, NamedTuple, Optional

from shakelab.libutils.timeN import Date
from shakelab.signals.base import StreamCollection
from shakelab.signals.codes import match_code
from shakelab.signals.fdsnws.registry import (
    SERVICES,
    find_datacenter,
    service_root,
    service_urls,
)
from shakelab.signals.fdsnws.selection import (
    Request,
    get_params,
    parse_selection,
    post_body,
    to_date,
)
from shakelab.signals.fdsnws.transport import (
    DEFAULT_RETRIES,
    DEFAULT_TIMEOUT,
    FdsnError,
    FdsnWarning,
    request,
)
from shakelab.signals.metadata import MetadataCollection

__all__ = ["Client", "TimeSpan"]

# Margin of corrected or filtered downloads: a fraction of the window,
# at least PAD_MINIMUM seconds.
PAD_FRACTION = 0.1
PAD_MINIMUM = 30.0

_LEVELS = ("station", "channel", "response")
_MERGE = ("samplerate", "quality", "overlap")


class TimeSpan(NamedTuple):
    """
    Interval of available data of a channel.

    Attributes
    ----------
    sid : str
        Source identifier ``NET.STA.LOC.CHA``.
    starttime, endtime : Date
        First and last sample of the interval.
    quality : str or None
        Quality code of the data, if reported.
    sample_rate : float or None
        Sampling rate in Hz, if reported.
    """

    sid: str
    starttime: Date
    endtime: Date
    quality: Optional[str] = None
    sample_rate: Optional[float] = None


class Client:
    """
    Client of the FDSN web services of a data center.

    Parameters
    ----------
    datacenter : str, optional
        Name of a data center of the registry (see
        :func:`~shakelab.signals.fdsnws.registry.datacenters`; case and
        former names such as ``IRIS`` or ``ORFEUS`` are accepted), or
        the ``http(s)://`` address of a host with the standard paths
        ``/fdsnws/<service>/1/``.
    station, dataselect, event, availability : str, optional
        Address of a service, replacing that of ``datacenter`` (also
        for services with non-standard paths). At least one service
        must be known.
    timeout : float, optional
        Seconds without data before a request is abandoned. Default
        is 120.
    retries : int, optional
        Attempts after a failure worth retrying (timeouts, connection
        failures, busy servers). Default is 3.

    Attributes
    ----------
    name : str or None
        Name of the data center, if from the registry.
    services : dict
        Address of each service.

    Raises
    ------
    TypeError
        If an argument is of a wrong type.
    ValueError
        If the data center is unknown, an address is invalid, or no
        service is given.

    Examples
    --------
    >>> Client("ingv").services["station"]
    'https://webservices.ingv.it/fdsnws/station/1/'
    """

    def __init__(
        self,
        datacenter: Optional[str] = None,
        *,
        station: Optional[str] = None,
        dataselect: Optional[str] = None,
        event: Optional[str] = None,
        availability: Optional[str] = None,
        timeout: Real = DEFAULT_TIMEOUT,
        retries: int = DEFAULT_RETRIES,
    ) -> None:
        if isinstance(timeout, bool) or not isinstance(timeout, Real):
            raise TypeError("timeout must be a real number")
        if timeout <= 0:
            raise ValueError("timeout must be positive")
        if isinstance(retries, bool) or not isinstance(retries, int):
            raise TypeError("retries must be an integer")
        if retries < 0:
            raise ValueError("retries must be non-negative")

        self.name: Optional[str] = None
        self.services: dict = {}

        if datacenter is not None:
            if not isinstance(datacenter, str):
                raise TypeError("datacenter must be a string")
            if datacenter.strip().lower().startswith(("http://", "https://")):
                self.services = service_urls(datacenter)
            else:
                center = find_datacenter(datacenter)
                self.name = center.name
                self.services = dict(center.services)

        explicit = {
            "station": station,
            "dataselect": dataselect,
            "event": event,
            "availability": availability,
        }
        for service, url in explicit.items():
            if url is not None:
                self.services[service] = service_root(url)

        if not self.services:
            raise ValueError("give a data center or a service address")

        self.timeout = float(timeout)
        self.retries = retries

    def __repr__(self) -> str:
        if self.name is not None:
            return f"Client({self.name!r})"
        return f"Client(services={sorted(self.services)})"

    # -- services -------------------------------------------------------

    def version(self, service: str = "dataselect") -> str:
        """
        Return the version of a service.

        Parameters
        ----------
        service : {"dataselect", "station", "event", "availability"}
            Service to ask. Default is ``"dataselect"``.

        Returns
        -------
        str
            Version reported by the service (e.g. ``"1.1.64"``).

        Raises
        ------
        ValueError
            If the data center does not offer the service.
        FdsnError
            If the request fails.
        """
        content = request(
            self._url(service, "version"),
            timeout=self.timeout,
            retries=self.retries,
        )
        return "" if content is None else content.decode("utf-8").strip()

    def get_waveforms(
        self,
        select: Any,
        starttime: Any = None,
        endtime: Any = None,
        *,
        metadata: bool = False,
        correct: bool = False,
        pad: Optional[Real] = None,
        highpass: Optional[Real] = None,
        lowpass: Optional[Real] = None,
        prefilter: Optional[tuple] = None,
        waterlevel: Real = 60.0,
        chunk: Optional[Real] = None,
        quality: Optional[str] = None,
        minimumlength: Optional[Real] = None,
        longestonly: bool = False,
        mseed_backend: str = "auto",
        stream_collection: Optional[StreamCollection] = None,
        target: Any = None,
        owrite: bool = False,
    ) -> Any:
        """
        Download waveforms.

        Parameters
        ----------
        select : str or sequence
            Channel patterns, or ``(pattern, starttime, endtime)``
            tuples (see
            :func:`~shakelab.signals.fdsnws.selection.parse_selection`).
        starttime, endtime : Date or str, optional
            Time window of the patterns without one; every pattern
            needs a window.
        metadata : bool, optional
            If True, also download the channel metadata (with the
            responses) and assign them to the records.
        correct : bool, optional
            If True, correct the records for the instrument response
            (implies ``metadata``), as described in the module notes.
            Records without a response are left out with a
            :class:`FdsnWarning`.
        pad : float, optional
            Margin in seconds requested before and after each window,
            and tapered, when the records are corrected or filtered.
            Default: 10 % of the window, at least 30 s.
        highpass, lowpass : float, optional
            Corner frequencies in Hz of a Butterworth filter applied
            to the padded records (after the correction, if any), as
            :meth:`~shakelab.signals.base.Record.filter`.
        prefilter : tuple of float, optional
            Corner frequencies ``(f1, f2, f3, f4)`` in Hz of the
            frequency window of the correction (see
            :meth:`~shakelab.signals.base.Record.deconvolve_response`);
            only with ``correct``.
        waterlevel : float, optional
            Water level of the correction in dB. Default is 60.
        chunk : float, optional
            Longest window, in seconds, of a single request: longer
            windows are split, so that the servers do not refuse them
            as too large.
        quality : str, optional
            Data quality code (``D``, ``R``, ``Q``, ``M`` or ``B``).
        minimumlength : float, optional
            Shortest continuous segment, in seconds, that the server
            returns.
        longestonly : bool, optional
            If True, the server returns only the longest continuous
            segment of each channel.
        mseed_backend : {"auto", "libmseed", "python"}, optional
            MiniSEED reader (see :func:`shakelab.signals.io.reader`).
        stream_collection : StreamCollection, optional
            Collection to which the records are added; a new one if
            omitted. Metadata and corrections are applied only to the
            new records.
        target : str or path-like, optional
            File where the MiniSEED records are written as received,
            in place of returning them; incompatible with
            ``metadata``, ``correct`` and the filters.
        owrite : bool, optional
            If True, ``target`` may be replaced. Default is False.

        Returns
        -------
        StreamCollection or int
            The collection (empty if there are no data), or the number
            of bytes written to ``target``.

        Raises
        ------
        TypeError
            If an argument is of a wrong type.
        ValueError
            If a window is missing, the options are incompatible, or
            the data center has no dataselect service.
        FileExistsError
            If ``target`` exists and ``owrite`` is False.
        FdsnError
            If a request fails.
        """
        _check_bool(metadata, "metadata")
        _check_bool(correct, "correct")
        _check_bool(longestonly, "longestonly")
        _check_bool(owrite, "owrite")

        requests = parse_selection(select, starttime, endtime)
        for item in requests:
            if item.starttime is None or item.endtime is None:
                raise ValueError(
                    f"{item.pattern}: waveforms need starttime and endtime"
                )

        filtering = highpass is not None or lowpass is not None
        processing = correct or filtering
        if prefilter is not None and not correct:
            raise ValueError("prefilter is used only with correct=True")
        if target is not None and (metadata or processing):
            raise ValueError(
                "target writes the records as received: it cannot be "
                "combined with metadata, correct or filters"
            )

        if processing:
            margin = _pad(pad, requests)
            requests = [
                Request(item.pattern, item.starttime - margin,
                        item.endtime + margin)
                for item in requests
            ]
        elif pad is not None:
            raise ValueError("pad is used only to correct or filter")
        else:
            margin = 0.0

        params = {
            "quality": quality,
            "minimumlength": minimumlength,
            "longestonly": True if longestonly else None,
        }
        groups = _split(requests, _check_positive(chunk, "chunk"))

        if target is not None:
            return self._download(
                "dataselect", groups, params, target, owrite
            )

        if stream_collection is None:
            stream_collection = StreamCollection()
        elif not isinstance(stream_collection, StreamCollection):
            raise TypeError("stream_collection must be a StreamCollection")

        from shakelab.signals.io import read_mseed_bytes

        received = StreamCollection()
        for group in groups:
            content = self._query("dataselect", "query", group, params)
            if not content:
                continue
            if len(groups) == 1:
                read_mseed_bytes(content, mseed_backend=mseed_backend,
                                 stream_collection=received)
            else:
                piece = read_mseed_bytes(content,
                                         mseed_backend=mseed_backend)
                _add_piece(piece, group, requests, received)

        pairs = _select_records(received, requests)
        if not pairs:
            return stream_collection

        if metadata or correct:
            channels = self.get_metadata(
                [
                    (str(item.pattern), item.starttime, item.endtime)
                    for item in requests
                ],
                level="response",
            )
            missing = []
            for record, _ in pairs:
                if record.sid in channels:
                    record.set_metadata(channels)
                elif record.sid not in missing:
                    missing.append(record.sid)
            if missing and not correct:
                warnings.warn(
                    "no metadata for " + ", ".join(missing),
                    FdsnWarning,
                    stacklevel=2,
                )

        if processing:
            pairs = _process(
                pairs, margin, correct, highpass, lowpass,
                prefilter, waterlevel,
            )

        for record, item in pairs:
            if margin:
                start = item.starttime + margin
                end = item.endtime - margin
                if not record.overlaps(start, end):
                    continue
                record = record.cut(start, end, inplace=False)
            stream_collection.append(record)

        return stream_collection

    def get_metadata(
        self,
        select: Any = "*",
        starttime: Any = None,
        endtime: Any = None,
        *,
        level: str = "channel",
        text: bool = False,
        region: Optional[tuple] = None,
        center: Optional[tuple] = None,
        radius: Optional[tuple] = None,
        startbefore: Any = None,
        startafter: Any = None,
        endbefore: Any = None,
        endafter: Any = None,
        includerestricted: Optional[bool] = None,
        matchtimeseries: Optional[bool] = None,
        updatedafter: Any = None,
        target: Any = None,
        owrite: bool = False,
    ) -> Any:
        """
        Download station and channel metadata.

        Parameters
        ----------
        select : str or sequence, optional
            Channel patterns, or ``(pattern, starttime, endtime)``
            tuples; at station level only the network and station
            codes count. Default: every channel.
        starttime, endtime : Date or str, optional
            Epochs active in this window.
        level : {"station", "channel", "response"}, optional
            Detail of the metadata: station epochs, channel epochs, or
            channel epochs with their responses. Default is
            ``"channel"``.
        text : bool, optional
            If True, request the text format: lighter, with
            coordinates, orientation, sampling rate and overall
            sensitivity only (see
            :mod:`~shakelab.signals.stationtext`); not at
            response level.
        region : tuple of float, optional
            ``(lonmin, lonmax, latmin, latmax)`` in degrees.
        center : tuple of float, optional
            ``(lon, lat)`` of a circular region, with ``radius``.
        radius : tuple of float, optional
            ``(rmin, rmax)`` distance from ``center`` in degrees;
            either may be None.
        startbefore, startafter, endbefore, endafter : Date or str, optional
            Constraints on the start and end of the epochs.
        includerestricted : bool, optional
            Whether restricted channels are listed (server default:
            yes).
        matchtimeseries : bool, optional
            If True, only channels with data in the archive.
        updatedafter : Date or str, optional
            Only metadata updated after this time.
        target : str or path-like, optional
            File where the reply (StationXML or text) is written, in
            place of returning it.
        owrite : bool, optional
            If True, ``target`` may be replaced. Default is False.

        Returns
        -------
        MetadataCollection or int
            The metadata (empty if there are none), or the number of
            bytes written to ``target``.

        Raises
        ------
        TypeError
            If an argument is of a wrong type.
        ValueError
            If an argument is invalid, or the data center has no
            station service.
        FileExistsError
            If ``target`` exists and ``owrite`` is False.
        FdsnError
            If the request fails.
        """
        _check_bool(text, "text")
        _check_bool(owrite, "owrite")
        if level not in _LEVELS:
            raise ValueError(f"level must be one of {', '.join(_LEVELS)}")
        if text and level == "response":
            raise ValueError("the text format has no response level")

        requests = parse_selection(select, starttime, endtime)
        params = {
            "level": level,
            "format": "text" if text else None,
            "startbefore": to_date(startbefore, "startbefore"),
            "startafter": to_date(startafter, "startafter"),
            "endbefore": to_date(endbefore, "endbefore"),
            "endafter": to_date(endafter, "endafter"),
            "includerestricted": _optional_bool(
                includerestricted, "includerestricted"
            ),
            "matchtimeseries": _optional_bool(
                matchtimeseries, "matchtimeseries"
            ),
            "updatedafter": to_date(updatedafter, "updatedafter"),
            **_geographic(region, center, radius),
        }

        if target is not None:
            return self._download("station", [requests], params, target,
                                  owrite)

        content = self._query("station", "query", requests, params)
        collection = MetadataCollection()
        if not content:
            return collection

        if text:
            from shakelab.signals.stationtext import read_station_text

            received = read_station_text(content)
        else:
            from shakelab.signals.stationxml import read_stationxml

            received = read_stationxml(content)

        if level == "station":
            for station in received.stations():
                if any(
                    match_code(
                        station.code,
                        f"{item.pattern.network}.{item.pattern.station}",
                    )
                    for item in requests
                ):
                    collection.append(station)
        else:
            for epoch in received:
                if any(item.matches(epoch.sid) for item in requests):
                    collection.append(epoch)
        return collection

    def get_events(
        self,
        starttime: Any = None,
        endtime: Any = None,
        *,
        region: Optional[tuple] = None,
        center: Optional[tuple] = None,
        radius: Optional[tuple] = None,
        depth: Optional[tuple] = None,
        magnitude: Optional[tuple] = None,
        magnitudetype: Optional[str] = None,
        eventtype: Optional[str] = None,
        catalog: Optional[str] = None,
        contributor: Optional[str] = None,
        eventid: Optional[str] = None,
        limit: Optional[int] = None,
        offset: Optional[int] = None,
        orderby: Optional[str] = None,
        includeallorigins: Optional[bool] = None,
        includeallmagnitudes: Optional[bool] = None,
        includearrivals: Optional[bool] = None,
        updatedafter: Any = None,
        target: Any = None,
        owrite: bool = False,
    ) -> Any:
        """
        Download an earthquake catalogue (QuakeML).

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Origin time window.
        region : tuple of float, optional
            ``(lonmin, lonmax, latmin, latmax)`` in degrees.
        center : tuple of float, optional
            ``(lon, lat)`` of a circular region, with ``radius``.
        radius : tuple of float, optional
            ``(rmin, rmax)`` distance from ``center`` in degrees;
            either may be None.
        depth : tuple of float, optional
            ``(min, max)`` depth in km, positive downward; either may
            be None.
        magnitude : tuple of float, optional
            ``(min, max)`` magnitude; either may be None.
        magnitudetype : str, optional
            Magnitude type used by ``magnitude`` (e.g. ``"ML"``).
        eventtype : str, optional
            Event types, separated by commas.
        catalog, contributor : str, optional
            Catalogue and contributor of the solutions.
        eventid : str, optional
            Identifier of one event.
        limit, offset : int, optional
            Maximum number of events, and index of the first one
            (starting from 1).
        orderby : str, optional
            ``"time"`` (default of the services), ``"time-asc"``,
            ``"magnitude"`` or ``"magnitude-asc"``.
        includeallorigins, includeallmagnitudes, includearrivals : bool
            Whether to include every origin, every magnitude and the
            phase arrivals (server defaults: no).
        updatedafter : Date or str, optional
            Only events updated after this time.
        target : str or path-like, optional
            File where the QuakeML reply is written, in place of
            returning it.
        owrite : bool, optional
            If True, ``target`` may be replaced. Default is False.

        Returns
        -------
        EqDatabase or int
            The catalogue (see :mod:`shakelab.seismicity.catalogue`;
            empty if there are no events), or the number of bytes
            written to ``target``.

        Raises
        ------
        TypeError
            If an argument is of a wrong type.
        ValueError
            If an argument is invalid, or the data center has no event
            service.
        FileExistsError
            If ``target`` exists and ``owrite`` is False.
        FdsnError
            If the request fails.
        """
        _check_bool(owrite, "owrite")
        depth_min, depth_max = _range(depth, "depth")
        mag_min, mag_max = _range(magnitude, "magnitude")
        for name, value in (("limit", limit), ("offset", offset)):
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int)
                or value < 1
            ):
                raise ValueError(f"{name} must be a positive integer")

        params = {
            "starttime": to_date(starttime, "starttime"),
            "endtime": to_date(endtime, "endtime"),
            **_geographic(region, center, radius),
            "mindepth": depth_min,
            "maxdepth": depth_max,
            "minmagnitude": mag_min,
            "maxmagnitude": mag_max,
            "magnitudetype": _optional_str(magnitudetype, "magnitudetype"),
            "eventtype": _optional_str(eventtype, "eventtype"),
            "catalog": _optional_str(catalog, "catalog"),
            "contributor": _optional_str(contributor, "contributor"),
            "eventid": _optional_str(eventid, "eventid"),
            "limit": limit,
            "offset": offset,
            "orderby": _optional_str(orderby, "orderby"),
            "includeallorigins": _optional_bool(
                includeallorigins, "includeallorigins"
            ),
            "includeallmagnitudes": _optional_bool(
                includeallmagnitudes, "includeallmagnitudes"
            ),
            "includearrivals": _optional_bool(
                includearrivals, "includearrivals"
            ),
            "updatedafter": to_date(updatedafter, "updatedafter"),
        }

        if target is not None:
            return self._download("event", [None], params, target, owrite)

        # The catalogue classes belong to another subpackage.
        from shakelab.seismicity.catalogue import EqDatabase
        from shakelab.seismicity.quakeml import parse_quakeml

        content = self._query("event", "query", None, params)
        if not content:
            return EqDatabase()
        return parse_quakeml(content.decode("utf-8"))

    def get_availability(
        self,
        select: Any,
        starttime: Any = None,
        endtime: Any = None,
        *,
        extent: bool = False,
        merge: Any = None,
        mergegaps: Optional[Real] = None,
        quality: Optional[str] = None,
        includerestricted: Optional[bool] = None,
    ) -> list:
        """
        Return the time spans of the available data.

        Parameters
        ----------
        select : str or sequence
            Channel patterns, or ``(pattern, starttime, endtime)``
            tuples.
        starttime, endtime : Date or str, optional
            Time window of the patterns without one.
        extent : bool, optional
            If True, only the first and last sample of each channel
            (method ``extent``); otherwise every continuous span
            (method ``query``).
        merge : str or sequence of str, optional
            Spans merged across ``"samplerate"``, ``"quality"`` and,
            not with ``extent``, ``"overlap"``.
        mergegaps : float, optional
            Gaps up to this length in seconds are ignored (not with
            ``extent``).
        quality : str, optional
            Data quality code.
        includerestricted : bool, optional
            Whether restricted data are included (server default: no).

        Returns
        -------
        list of TimeSpan
            The spans, in the order of the server.

        Raises
        ------
        TypeError
            If an argument is of a wrong type.
        ValueError
            If an argument is invalid, or the data center has no
            availability service.
        FdsnError
            If the request fails.
        """
        _check_bool(extent, "extent")
        if merge is not None:
            values = [merge] if isinstance(merge, str) else list(merge)
            for value in values:
                if value not in _MERGE or (extent and value == "overlap"):
                    raise ValueError(f"invalid merge value: {value!r}")
            merge = ",".join(values)
        if extent and mergegaps is not None:
            raise ValueError("mergegaps is not used with extent")

        requests = parse_selection(select, starttime, endtime)
        params = {
            "merge": merge,
            "mergegaps": mergegaps,
            "quality": _optional_str(quality, "quality"),
            "includerestricted": _optional_bool(
                includerestricted, "includerestricted"
            ),
            "format": "json",
        }
        method = "extent" if extent else "query"
        content = self._query("availability", method, requests, params)
        if not content:
            return []

        try:
            data = json.loads(content.decode("utf-8"))
            spans = _parse_availability(data, extent)
        except (ValueError, KeyError, TypeError) as error:
            raise FdsnError(
                f"invalid availability reply: {error}"
            ) from None
        return [
            span for span in spans
            if any(item.matches(span.sid) for item in requests)
        ]

    # -- helpers --------------------------------------------------------

    def _url(self, service: str, method: str) -> str:
        """Address of a method of a service."""
        if service not in SERVICES:
            raise ValueError(
                f"service must be one of {', '.join(SERVICES)}"
            )
        if service not in self.services:
            name = self.name or "this client"
            raise ValueError(f"{name} has no {service} service")
        return self.services[service] + method

    def _query(
        self,
        service: str,
        method: str,
        requests: Optional[list],
        params: dict,
        target: Any = None,
    ) -> Any:
        """Send one request: GET for one selection line, else POST."""
        url = self._url(service, method)
        if requests is None:
            return request(url, params=params, target=target,
                           timeout=self.timeout, retries=self.retries)

        line = get_params(requests)
        if line is not None:
            return request(url, params={**line, **params}, target=target,
                           timeout=self.timeout, retries=self.retries)

        body = post_body(requests, params)
        return request(url, body=body, target=target,
                       timeout=self.timeout, retries=self.retries)

    def _download(
        self,
        service: str,
        groups: list,
        params: dict,
        target: Any,
        owrite: bool,
    ) -> int:
        """Write the replies of the request groups to one file."""
        if not isinstance(target, (str, os.PathLike)):
            raise TypeError("target must be a path")
        path = os.fspath(target)
        if os.path.exists(path) and not owrite:
            raise FileExistsError(f"{path} exists (use owrite=True)")

        if len(groups) == 1:
            size = self._query(service, "query", groups[0], params, path)
            return size or 0

        size = 0
        piece = f"{path}.piece"
        try:
            with open(f"{path}.part", "wb") as output:
                for group in groups:
                    written = self._query(
                        service, "query", group, params, piece
                    )
                    if not written:
                        continue
                    with open(piece, "rb") as fobj:
                        while True:
                            block = fobj.read(1 << 16)
                            if not block:
                                break
                            output.write(block)
                    size += written
            if size:
                os.replace(f"{path}.part", path)
        finally:
            for name in (piece, f"{path}.part"):
                if os.path.exists(name):
                    os.remove(name)
        return size


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _pad(pad: Optional[Real], requests: list) -> float:
    """Margin of processed downloads, in seconds."""
    if pad is None:
        longest = max(item.endtime - item.starttime for item in requests)
        return max(PAD_FRACTION * longest, PAD_MINIMUM)
    if isinstance(pad, bool) or not isinstance(pad, Real):
        raise TypeError("pad must be a real number")
    if not pad >= 0:
        raise ValueError("pad must be non-negative")
    return float(pad)


def _split(requests: list, chunk: Optional[float]) -> list:
    """Groups of requests, with windows of at most ``chunk`` seconds."""
    if chunk is None:
        return [requests]
    pieces: list = []
    for item in requests:
        start = item.starttime
        index = 0
        while True:
            end = start + chunk
            last = end >= item.endtime
            if index == len(pieces):
                pieces.append([])
            pieces[index].append(
                Request(item.pattern, start, item.endtime if last else end)
            )
            if last:
                break
            start, index = end, index + 1
    return pieces


def _add_piece(
    piece: StreamCollection,
    group: list,
    requests: list,
    received: StreamCollection,
) -> None:
    """
    Add the records of a chunk, cut to its window.

    The servers return whole MiniSEED records, so that consecutive
    chunks overlap: each chunk is cut to its window, without the sample
    at its end unless it is the end of the request, and the pieces join
    again in ``received``.
    """
    ends = {(item.pattern, item.endtime) for item in requests}
    for stream in piece:
        for record in stream:
            for item in group:
                if not item.matches(record.sid):
                    continue
                end = item.endtime
                if (item.pattern, end) not in ends:
                    end = end - 0.5 * record.delta
                if not record.overlaps(item.starttime, end):
                    continue
                received.append(
                    record.cut(item.starttime, end, inplace=False)
                )


def _select_records(received: StreamCollection, requests: list) -> list:
    """Records matching each request, cut to its window, with it."""
    pairs = []
    for stream in received:
        for record in stream:
            for item in requests:
                if not item.matches(record.sid):
                    continue
                if not record.overlaps(item.starttime, item.endtime):
                    continue
                pairs.append(
                    (
                        record.cut(
                            item.starttime, item.endtime, inplace=False
                        ),
                        item,
                    )
                )
    return pairs


def _process(
    pairs: list,
    margin: float,
    correct: bool,
    highpass: Optional[Real],
    lowpass: Optional[Real],
    prefilter: Optional[tuple],
    waterlevel: Real,
) -> list:
    """Detrend, taper, correct and filter the padded records."""
    kept = []
    missing = []
    for record, item in pairs:
        record.detrend()
        if margin:
            record.taper(margin)
        if correct:
            try:
                record.deconvolve_response(
                    waterlevel=waterlevel,
                    prefilter=prefilter,
                    padding=True,
                )
            except ValueError as error:
                missing.append(f"{record.sid} ({error})")
                continue
        if highpass is not None or lowpass is not None:
            record.filter(highpass=highpass, lowpass=lowpass)
        kept.append((record, item))

    if missing:
        warnings.warn(
            "records left out, not corrected: " + "; ".join(missing),
            FdsnWarning,
            stacklevel=3,
        )
    return kept


def _parse_availability(data: dict, extent: bool) -> list:
    """Time spans of the JSON reply of the availability service."""
    spans = []
    for source in data.get("datasources", []):
        sid = ".".join(
            source.get(key, "")
            for key in ("network", "station", "location", "channel")
        )
        quality = source.get("quality")
        rate = source.get("samplerate")
        rate = None if rate is None else float(rate)
        if extent:
            intervals = [(source["earliest"], source["latest"])]
        else:
            intervals = source.get("timespans", [])
        for start, end in intervals:
            spans.append(
                TimeSpan(sid, Date(start), Date(end), quality, rate)
            )
    return spans


def _geographic(
    region: Optional[tuple],
    center: Optional[tuple],
    radius: Optional[tuple],
) -> dict:
    """Parameters of a rectangular or circular region."""
    params: dict = {}
    if region is not None:
        values = _reals(region, 4, "region")
        params.update(
            minlongitude=values[0],
            maxlongitude=values[1],
            minlatitude=values[2],
            maxlatitude=values[3],
        )
    if center is not None:
        lon, lat = _reals(center, 2, "center")
        params.update(longitude=lon, latitude=lat)
        rmin, rmax = _range(radius, "radius")
        params.update(minradius=rmin, maxradius=rmax)
    elif radius is not None:
        raise ValueError("radius needs center")
    return params


def _reals(value: Any, size: int, name: str) -> tuple:
    """A tuple of ``size`` real numbers."""
    if (
        not isinstance(value, (tuple, list)) or len(value) != size
        or any(
            isinstance(item, bool) or not isinstance(item, Real)
            for item in value
        )
    ):
        raise TypeError(f"{name} must be a sequence of {size} numbers")
    return tuple(float(item) for item in value)


def _range(value: Any, name: str) -> tuple:
    """``(min, max)`` with None allowed for either bound."""
    if value is None:
        return None, None
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise TypeError(f"{name} must be a (min, max) pair")
    bounds = []
    for item in value:
        if item is None:
            bounds.append(None)
        elif isinstance(item, bool) or not isinstance(item, Real):
            raise TypeError(f"{name} bounds must be numbers or None")
        else:
            bounds.append(float(item))
    if None not in bounds and bounds[1] < bounds[0]:
        raise ValueError(f"{name}: max is smaller than min")
    return tuple(bounds)


def _check_bool(value: Any, name: str) -> None:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a boolean")


def _optional_bool(value: Any, name: str) -> Optional[bool]:
    if value is not None:
        _check_bool(value, name)
    return value


def _optional_str(value: Any, name: str) -> Optional[str]:
    if value is not None and not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    return value


def _check_positive(value: Any, name: str) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    if not value > 0:
        raise ValueError(f"{name} must be positive")
    return float(value)
