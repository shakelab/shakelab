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
SeedLink client.

A :class:`Client` connects to a SeedLink server, requests the streams
selected with :meth:`Client.select` and distributes the records it
receives to its destinations:

- a buffer (:meth:`Client.set_buffer`): a
  :class:`~shakelab.signals.base.StreamCollection`, unlimited or
  circular (records older than a given duration are discarded);
- callbacks: at every record (:meth:`Client.on_record`), at regular
  intervals of wall-clock or data time (:meth:`Client.every`), and at
  connection events and data gaps (:meth:`Client.on_event`).

The client runs in the calling thread (:meth:`Client.run`, blocking)
or in a background thread (:meth:`Client.start`, :meth:`Client.stop`).
It reconnects after network failures and resumes each station from the
packet after the last one received, so that no data are lost while the
server keeps them. Callbacks run in the thread that receives the data;
other threads read the buffer with :meth:`Client.get_data`.

:func:`read` collects a time window, or a given duration, in one call.

Examples
--------
>>> from shakelab.signals import seedlink
>>> client = seedlink.Client("rtserve.example.org")
>>> client.select(["IV.MODE.*.HH?", "OX.*.*.HHZ"])
>>> client.set_buffer(max_duration=300)
>>> client.every(10.0, lambda client: print(client.buffer))
>>> client.run(duration=60)
"""

from __future__ import annotations

import fnmatch
import logging
import threading
import time as _time
import warnings
from dataclasses import dataclass
from numbers import Real
from typing import Any, Callable, Iterator, NamedTuple, Optional

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.signals import metadata as _metadata
from shakelab.signals.base import (
    MERGE_TOLERANCE,
    Record,
    StreamCollection,
    _compatible_meta,
    _set_metadata,
)
from shakelab.signals.seedlink.protocol import (
    DEFAULT_PORT,
    Connection,
    Packet,
    SeedLinkError,
    data_command,
    is_newer,
    next_sequence,
    normalize_protocol,
    select_commands,
    start_command,
    station_command,
)

__all__ = ["Client", "SeedLinkEvent", "read"]

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

# Seconds between checks of the stop request and of the timers while
# no data arrive.
POLL_INTERVAL = 0.5

# Seconds to wait for the connection and for the replies to commands.
CONNECT_TIMEOUT = 30.0

# Delay before reconnecting: doubled after every failed attempt, from
# the first to the second value (seconds).
RECONNECT_DELAY = (1.0, 60.0)

# Seconds without packets after which a time window that ended in the
# past is considered complete.
WINDOW_IDLE = 10.0

# Event kinds reported to the on_event callbacks.
EVENT_KINDS = (
    "connected", "disconnected", "gap", "overlap", "end",
)

_NS = 10 ** 9


# ----------------------------------------------------------------------------
# Events and selection
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class SeedLinkEvent:
    """
    Event reported to the callbacks of :meth:`Client.on_event`.

    Attributes
    ----------
    kind : str
        ``"connected"``, ``"disconnected"`` (the client reconnects if
        allowed), ``"gap"`` or ``"overlap"`` (between two records of a
        channel), ``"end"`` (the server has sent all the requested
        data).
    time : Date
        Wall-clock time of the event (UTC).
    message : str
        Description.
    sid : str or None
        Channel of a gap or overlap.
    start, end : Date or None
        Interval of a gap or overlap: from the expected start of the
        record to its actual start.
    """

    kind: str
    time: Date
    message: str = ""
    sid: Optional[str] = None
    start: Optional[Date] = None
    end: Optional[Date] = None


class _Pattern(NamedTuple):
    """Channel pattern ``NET.STA.LOC.CHA`` with ``*`` and ``?``."""

    network: str
    station: str
    location: str
    channel: str

    def matches(self, sid: Optional[str]) -> bool:
        """Whether a source identifier matches the pattern."""
        if sid is None:
            return False
        parts = sid.split(".")
        if len(parts) != 4:
            return False
        return all(
            fnmatch.fnmatchcase(code, pattern)
            for code, pattern in zip(parts, self)
        )


@dataclass
class _Timer:
    """Callback of :meth:`Client.every`."""

    interval: float
    func: Callable
    clock: str
    due: Optional[float] = None


@dataclass
class _RunSettings:
    """Stop conditions and time window of a run."""

    duration: Optional[float]
    start_ns: Optional[int]
    end_ns: Optional[int]
    fetch: bool


class _IdleTimeout(ConnectionError):
    """No data received for longer than the timeout."""


# ----------------------------------------------------------------------------
# Client
# ----------------------------------------------------------------------------


class Client:
    """
    SeedLink client.

    Parameters
    ----------
    host : str
        Server host name or address.
    port : int, optional
        Server port. Default is 18000.
    protocol : {"auto", "3.1", "4.0"}, optional
        Protocol version. ``"auto"`` (default) uses version 4 when the
        server supports it, version 3.1 otherwise.
    metadata : MetadataCollection, path or sequence of paths, optional
        Channel metadata assigned to the records as they arrive (see
        :func:`shakelab.signals.io.read_metadata`).
    timeout : float, optional
        Seconds without data after which the connection is considered
        lost and reopened. Default is 600 s.
    keepalive : float or None, optional
        Seconds without data after which the client sends a keep-alive
        request, so that idle connections are not closed by firewalls.
        Default is 60 s; None disables it.
    reconnect : bool, optional
        Whether to reconnect after a network failure. Default is True.
        A failure of the first connection is always raised.

    Attributes
    ----------
    host : str
    port : int
    server : str or None
        Server software and version, from the last connection.
    organization : str or None
        Station or data center description, from the last connection.
    protocol : str or None
        Protocol version of the last connection.

    Notes
    -----
    Records reach the destinations in the order they arrive. They are
    filtered against the patterns of :meth:`select`, whatever the
    server sends.

    The client keeps the sequence number of the last packet of each
    station: after a reconnection, and in later runs of the same
    client, every station resumes from the next packet, as long as the
    server still holds it. Repeated packets are discarded.
    """

    def __init__(
        self,
        host: str,
        port: int = DEFAULT_PORT,
        protocol: str = "auto",
        metadata: Any = None,
        timeout: Real = 600.0,
        keepalive: Optional[Real] = 60.0,
        reconnect: bool = True,
    ) -> None:
        if not isinstance(host, str):
            raise TypeError("host must be a string")
        if isinstance(port, bool) or not isinstance(port, int):
            raise TypeError("port must be an integer")
        if not isinstance(reconnect, bool):
            raise TypeError("reconnect must be a boolean")

        self.host = host
        self.port = port
        self.timeout = _positive(timeout, "timeout")
        self.keepalive = (
            None if keepalive is None
            else _positive(keepalive, "keepalive")
        )
        self.reconnect = reconnect
        self._requested_protocol = normalize_protocol(protocol)
        self._metadata = _as_metadata(metadata)

        self.server: Optional[str] = None
        self.organization: Optional[str] = None
        self.protocol: Optional[str] = None

        self._patterns: list = []
        self._buffer: Optional[StreamCollection] = None
        self._pending: dict = {}
        self._record_callbacks: list = []
        self._timers: list = []
        self._event_callbacks: list = []

        self._sequence: dict = {}
        self._next_start: dict = {}
        self._run_end: dict = {}
        self._metadata_cache: dict = {}
        self._missing_metadata: set = set()

        self._lock = threading.RLock()
        self._stop_request = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._active = False
        self._error: Optional[BaseException] = None

    def __repr__(self) -> str:
        state = "running" if self.running else "idle"
        return (
            f"Client({self.host}:{self.port}, {len(self._patterns)} "
            f"pattern(s), {state})"
        )

    def __enter__(self) -> Client:
        return self

    def __exit__(self, *exc_info: Any) -> None:
        if self.running:
            self.stop()

    # -- configuration ---------------------------------------------------

    def select(self, patterns: Any) -> None:
        """
        Add channels to request.

        Parameters
        ----------
        patterns : str or sequence of str
            Source identifier patterns ``NET.STA.LOC.CHA`` with the
            ``*`` (any sequence) and ``?`` (one character) jokers, and
            ``[...]`` for one of a set of characters (as in
            ``HH[ZNE]``). Missing trailing codes match anything:
            ``"IV.MODE"`` requests every channel of MODE. An empty
            location (``"IV.MODE..HHZ"``) is the blank location code.

        Raises
        ------
        TypeError
            If a pattern is not a string.
        ValueError
            If a pattern is malformed.
        RuntimeError
            If the client is running.

        Examples
        --------
        >>> client.select("IV.MODE.*.HH?")
        >>> client.select(["IV.CAE*.*.HHZ", "OX.*.*.HH[ZN]"])
        """
        self._require_idle()
        if isinstance(patterns, str):
            patterns = [patterns]
        if not isinstance(patterns, (list, tuple)):
            raise TypeError("patterns must be a string or a sequence")
        parsed = [_parse_pattern(pattern) for pattern in patterns]
        for pattern in parsed:
            if pattern not in self._patterns:
                self._patterns.append(pattern)

    @property
    def selection(self) -> list:
        """Selected patterns, as ``NET.STA.LOC.CHA`` strings."""
        return [".".join(pattern) for pattern in self._patterns]

    def set_buffer(
        self,
        collection: Optional[StreamCollection] = None,
        max_duration: Optional[Real] = None,
    ) -> None:
        """
        Set the buffer that receives the records.

        Parameters
        ----------
        collection : StreamCollection, optional
            Existing collection to which the records are appended. If
            omitted, a new one is created.
        max_duration : float, optional
            For a new collection, time span in seconds kept for each
            channel (circular buffer). If None (default), the buffer
            keeps all the data.

        Raises
        ------
        TypeError
            If ``collection`` is not a StreamCollection.
        ValueError
            If both arguments are given.
        RuntimeError
            If the client is running.

        Notes
        -----
        An unlimited buffer grows with the data received: for long
        runs, set ``max_duration``. A circular buffer of an existing
        collection follows its own ``max_duration``.
        """
        self._require_idle()
        if collection is not None:
            if not isinstance(collection, StreamCollection):
                raise TypeError("collection must be a StreamCollection")
            if max_duration is not None:
                raise ValueError(
                    "max_duration applies to a new collection only"
                )
        else:
            collection = StreamCollection(max_duration=max_duration)
        self._buffer = collection

    @property
    def buffer(self) -> Optional[StreamCollection]:
        """
        Buffer of the records received (None if not set).

        Read it from the callbacks or when the client is not running.
        While the client runs in the background (:meth:`start`), use
        :meth:`get_data` from other threads.
        """
        with self._lock:
            self._flush()
            return self._buffer

    def on_record(self, func: Callable) -> Callable:
        """
        Call a function at every record received.

        Parameters
        ----------
        func : callable
            Called as ``func(client, record)`` after the record is added
            to the buffer.

        Returns
        -------
        callable
            ``func`` itself, so that the method can be used as a
            decorator.
        """
        self._require_idle()
        _require_callable(func)
        self._record_callbacks.append(func)
        return func

    def every(
        self,
        interval: Real,
        func: Callable,
        clock: str = "wall",
    ) -> Callable:
        """
        Call a function at regular intervals.

        Parameters
        ----------
        interval : float
            Interval in seconds.
        func : callable
            Called as ``func(client)``.
        clock : {"wall", "data"}, optional
            ``"wall"`` (default) counts wall-clock time from the start
            of the run. ``"data"`` counts the time of the data: the
            function is called when the end of the latest record passes
            the next multiple of ``interval`` from the first record,
            also when a backlog arrives faster than real time.

        Returns
        -------
        callable
            ``func`` itself.
        """
        self._require_idle()
        _require_callable(func)
        if clock not in ("wall", "data"):
            raise ValueError("clock must be 'wall' or 'data'")
        self._timers.append(
            _Timer(_positive(interval, "interval"), func, clock)
        )
        return func

    def on_event(self, func: Callable) -> Callable:
        """
        Call a function at connection events and data gaps.

        Parameters
        ----------
        func : callable
            Called as ``func(client, event)`` with a
            :class:`SeedLinkEvent`.

        Returns
        -------
        callable
            ``func`` itself.
        """
        self._require_idle()
        _require_callable(func)
        self._event_callbacks.append(func)
        return func

    # -- running ---------------------------------------------------------

    @property
    def running(self) -> bool:
        """Whether the client is receiving data."""
        return self._active or (
            self._thread is not None and self._thread.is_alive()
        )

    def run(
        self,
        duration: Optional[Real] = None,
        start: Any = None,
        end: Any = None,
        fetch: bool = False,
    ) -> Optional[StreamCollection]:
        """
        Receive data in the calling thread until a stop condition.

        Parameters
        ----------
        duration : float, optional
            Wall-clock seconds after which the run stops.
        start : Date or str, optional
            Start of the data requested: the server sends the data it
            still holds from this time on. If omitted, the transfer
            starts from the next packet (or resumes from the last one
            received in an earlier run).
        end : Date or str, optional
            End of the data requested. A past end time implies
            ``fetch=True``. The run stops when the wall clock and every
            channel received have passed it, or when no packet arrives
            for 10 s after it.
        fetch : bool, optional
            Dial-up mode: the server sends the data it holds and then
            ends the session. Default is False (real time).

        Returns
        -------
        StreamCollection or None
            The buffer (see :attr:`buffer`).

        Raises
        ------
        SeedLinkError
            If the server rejects the request.
        OSError
            If the server cannot be reached, or the connection is lost
            and ``reconnect`` is False.

        Notes
        -----
        Without any stop condition the run lasts until :meth:`stop` is
        called (e.g. by a callback) or the user interrupts it. If no
        buffer and no callback are set, an unlimited buffer is created.
        Packets that overlap the time window are kept whole.
        """
        if self.running:
            raise RuntimeError("the client is already running")
        settings = self._prepare(duration, start, end, fetch)
        self._stop_request.clear()
        self._loop(settings)
        return self.buffer

    def start(
        self,
        duration: Optional[Real] = None,
        start: Any = None,
        end: Any = None,
        fetch: bool = False,
    ) -> None:
        """
        Receive data in a background thread.

        The parameters are those of :meth:`run`. The thread ends at a
        stop condition or at :meth:`stop`; an error that ends it is
        raised again by :meth:`stop` or :meth:`wait`.
        """
        if self.running:
            raise RuntimeError("the client is already running")
        settings = self._prepare(duration, start, end, fetch)
        self._stop_request.clear()
        self._error = None
        self._thread = threading.Thread(
            target=self._background,
            args=(settings,),
            name=f"seedlink-{self.host}",
            daemon=True,
        )
        self._thread.start()

    def stop(self, timeout: Optional[Real] = None) -> None:
        """
        Stop receiving data.

        Called from a callback, the run stops after the callback
        returns. Called from another thread, it waits for the
        background thread to end.

        Parameters
        ----------
        timeout : float, optional
            Seconds to wait for the background thread.

        Raises
        ------
        Exception
            The error that ended the background thread, if any.
        """
        self._stop_request.set()
        self.wait(timeout)

    def wait(self, timeout: Optional[Real] = None) -> bool:
        """
        Wait for the background thread to end.

        Parameters
        ----------
        timeout : float, optional
            Seconds to wait; None waits indefinitely.

        Returns
        -------
        bool
            True if the client is no longer running.

        Raises
        ------
        Exception
            The error that ended the background thread, if any.
        """
        thread = self._thread
        if thread is not None and thread is not threading.current_thread():
            thread.join(timeout)
        if self._error is not None:
            error, self._error = self._error, None
            raise error
        return not (thread is not None and thread.is_alive())

    def get_data(
        self,
        sid: str = "*",
        start: Any = None,
        end: Any = None,
    ) -> StreamCollection:
        """
        Return a copy of the buffer, safe from any thread.

        Parameters
        ----------
        sid : str, optional
            Pattern of the source identifiers to copy, with the ``*``
            and ``?`` jokers. Default is every channel.
        start, end : Date or str, optional
            Time window to extract.

        Returns
        -------
        StreamCollection
            Independent copy of the selected data.

        Raises
        ------
        RuntimeError
            If no buffer is set.
        """
        if not isinstance(sid, str):
            raise TypeError("sid must be a string")
        with self._lock:
            self._flush()
            if self._buffer is None:
                raise RuntimeError("no buffer is set")
            collection = StreamCollection(
                max_duration=self._buffer.max_duration
            )
            collection.stream = [
                stream.copy()
                for stream in self._buffer
                if stream.sid is not None
                and fnmatch.fnmatchcase(stream.sid, sid)
            ]
        if start is not None or end is not None:
            collection = collection.extract(
                starttime=_as_date(start, "start"),
                endtime=_as_date(end, "end"),
            )
        return collection

    # -- run internals ---------------------------------------------------

    def _require_idle(self) -> None:
        if self.running:
            raise RuntimeError(
                "the configuration cannot change while the client runs"
            )

    def _prepare(self, duration, start, end, fetch) -> _RunSettings:
        """Validate the run parameters."""
        if not self._patterns:
            raise ValueError("no channel selected: call select() first")
        if duration is not None:
            duration = _positive(duration, "duration")
        if not isinstance(fetch, bool):
            raise TypeError("fetch must be a boolean")

        start = _as_date(start, "start")
        end = _as_date(end, "end")
        start_ns = start.to_unix_ns() if start is not None else None
        end_ns = end.to_unix_ns() if end is not None else None
        if start_ns is not None and end_ns is not None:
            if end_ns <= start_ns:
                raise ValueError("end must follow start")
        if end_ns is not None and end_ns <= Date.now().to_unix_ns():
            fetch = True

        if (
            self._buffer is None
            and not self._record_callbacks
            and not self._timers
        ):
            self._buffer = StreamCollection()
            logger.info("no destination set: using an unlimited buffer")

        return _RunSettings(duration, start_ns, end_ns, fetch)

    def _background(self, settings: _RunSettings) -> None:
        """Body of the background thread."""
        try:
            self._loop(settings)
        except BaseException as error:  # reported by stop() and wait()
            logger.error("SeedLink client stopped by an error: %s", error)
            self._error = error

    def _loop(self, settings: _RunSettings) -> None:
        """Connect, receive and reconnect until a stop condition."""
        self._active = True
        deadline = (
            _time.monotonic() + settings.duration
            if settings.duration is not None else None
        )
        for timer in self._timers:
            timer.due = None
        with self._lock:
            self._run_end.clear()
        delay = RECONNECT_DELAY[0]
        connected = False

        try:
            while not self._stop_request.is_set():
                connection = Connection(
                    self.host, self.port, timeout=CONNECT_TIMEOUT
                )
                try:
                    self._open(connection, settings)
                    connected = True
                    delay = RECONNECT_DELAY[0]
                    if self._receive(connection, settings, deadline):
                        return
                except OSError as error:
                    if not connected or not self.reconnect:
                        raise
                    self._emit("disconnected", message=str(error))
                    logger.warning(
                        "SeedLink connection lost (%s): reconnecting in "
                        "%.0f s", error, delay,
                    )
                finally:
                    connection.close()

                wait = delay
                if deadline is not None:
                    wait = min(wait, deadline - _time.monotonic())
                    if wait <= 0:
                        return
                if self._stop_request.wait(wait):
                    return
                delay = min(2 * delay, RECONNECT_DELAY[1])
        finally:
            self._active = False
            with self._lock:
                self._flush()

    def _open(self, connection: Connection, settings: _RunSettings) -> None:
        """Connect, negotiate and request the selected channels."""
        connection.open(self._requested_protocol)
        self.server = connection.server
        self.organization = connection.organization
        self.protocol = connection.version
        logger.info(
            "connected to %s:%d (%s), protocol %s",
            self.host, self.port, connection.server, connection.version,
        )

        for command in self._request_commands(connection, settings):
            connection.command(command)
        connection.send(
            start_command(connection.version, settings.fetch)
        )
        connection.set_poll(POLL_INTERVAL)
        self._emit("connected", message=connection.server)

    def _request_commands(
        self,
        connection: Connection,
        settings: _RunSettings,
    ) -> list:
        """
        Commands that request the selected channels.

        Stations already received resume from their next packet and are
        requested first; the station patterns follow (in version 4 they
        skip the stations already requested).
        """
        version = connection.version
        groups = self._station_groups(connection)
        commands = []
        explicit = set()

        def request(network, station, selectors, sequence):
            commands.append(
                station_command(version, network, station)
            )
            commands.extend(select_commands(version, selectors))
            commands.append(data_command(
                version,
                sequence=sequence,
                start_ns=settings.start_ns,
                end_ns=settings.end_ns,
                fetch=settings.fetch,
            ))

        for key in sorted(self._sequence):
            network, station = key.split("_", 1)
            selectors = []
            for (net, sta), group in groups.items():
                if (
                    fnmatch.fnmatchcase(network, net)
                    and fnmatch.fnmatchcase(station, sta)
                ):
                    selectors.extend(group)
            if selectors:
                sequence = next_sequence(
                    version, self._sequence[key]
                )
                request(network, station, selectors, sequence)
                explicit.add(key)

        for (network, station), selectors in groups.items():
            if f"{network}_{station}" in explicit:
                continue
            request(network, station, selectors, None)

        if not commands:
            raise SeedLinkError(
                "no station on the server matches the selection"
            )
        return commands

    def _station_groups(self, connection: Connection) -> dict:
        """
        Selected ``(location, channel)`` patterns by station pattern.

        Version 3 servers without the NSWILDCARD capability do not
        accept jokers in the station: they are expanded with the list of
        the stations of the server.
        """
        groups: dict = {}
        for pattern in self._patterns:
            key = (pattern.network, pattern.station)
            groups.setdefault(key, []).append(
                (pattern.location, pattern.channel)
            )

        if (
            connection.version == "4.0"
            or "NSWILDCARD" in connection.capabilities
            or not any(_has_joker("".join(key)) for key in groups)
        ):
            return groups

        available = _v3_stations(connection.info("STATIONS"))
        expanded: dict = {}
        for (network, station), selectors in groups.items():
            if not _has_joker(network + station):
                expanded.setdefault((network, station), []).extend(
                    selectors
                )
                continue
            matches = [
                key for key in available
                if fnmatch.fnmatchcase(key[0], network)
                and fnmatch.fnmatchcase(key[1], station)
            ]
            if not matches:
                logger.warning(
                    "no station of the server matches %s.%s",
                    network, station,
                )
            for key in matches:
                expanded.setdefault(key, []).extend(selectors)
        return expanded

    def _receive(
        self,
        connection: Connection,
        settings: _RunSettings,
        deadline: Optional[float],
    ) -> bool:
        """
        Receive packets until a stop condition (True) or a failure.
        """
        last_packet = _time.monotonic()
        last_keepalive = last_packet

        while True:
            now = _time.monotonic()
            if self._stop_request.is_set():
                return True
            if deadline is not None and now >= deadline:
                return True
            if settings.end_ns is not None and self._window_complete(
                settings.end_ns, now - last_packet
            ):
                return True
            self._run_wall_timers(now)

            try:
                packet = connection.read_packet()
            except SeedLinkError as error:
                # The stream is out of step: start a new connection.
                raise ConnectionError(f"protocol error: {error}") from None
            except TimeoutError:
                idle = _time.monotonic() - last_packet
                if idle > self.timeout:
                    raise _IdleTimeout(
                        f"no data for {idle:.0f} s"
                    ) from None
                if (
                    self.keepalive is not None
                    and idle > self.keepalive
                    and _time.monotonic() - last_keepalive > self.keepalive
                ):
                    connection.keepalive()
                    last_keepalive = _time.monotonic()
                continue

            last_packet = _time.monotonic()
            if packet is None:
                self._emit(
                    "end", message="the server has sent all the data"
                )
                return True
            self._handle_packet(packet, connection.version)

    def _window_complete(self, end_ns: int, idle: float) -> bool:
        """
        Whether the requested data up to ``end_ns`` have arrived.

        True when the wall clock has passed the end and every channel
        received in this run has reached it, or no packet has arrived
        for :data:`WINDOW_IDLE` seconds (the server has no more data).
        """
        if Date.now().to_unix_ns() < end_ns:
            return False
        if idle >= min(WINDOW_IDLE, self.timeout):
            return True
        with self._lock:
            reached = list(self._run_end.values())
        return bool(reached) and min(reached) >= end_ns

    # -- packets and records ---------------------------------------------

    def _handle_packet(self, packet: Packet, version: str) -> None:
        """Decode a packet and deliver its records."""
        if packet.kind != "data":
            logger.debug("%s packet skipped", packet.kind)
            return

        records, station = _decode(packet)
        if packet.sequence is not None and station is not None:
            last = self._sequence.get(station)
            if last is not None and not is_newer(
                version, packet.sequence, last
            ):
                logger.debug(
                    "repeated packet %d of %s skipped",
                    packet.sequence, station,
                )
                return
            self._sequence[station] = packet.sequence

        for record in records:
            self._deliver(record)

    def _deliver(self, record: Record) -> None:
        """Send a record to the buffer and to the callbacks."""
        sid = record.sid
        if not any(pattern.matches(sid) for pattern in self._patterns):
            logger.debug("record of %s not selected, skipped", sid)
            return

        if self._metadata is not None:
            self._assign_metadata(record)

        start_ns = record.time.to_unix_ns()
        end_ns = start_ns + round(record.nsamp * record.delta * _NS)
        self._check_continuity(record, start_ns)

        with self._lock:
            self._next_start[sid] = max(
                end_ns, self._next_start.get(sid, end_ns)
            )
            self._run_end[sid] = self._next_start[sid]
            if self._buffer is not None:
                self._store(record)

        for func in self._record_callbacks:
            func(self, record)
        self._run_data_timers(end_ns)

    def _store(self, record: Record) -> None:
        """
        Add a record to the buffer (lock held).

        Records for an unlimited buffer wait in a list, joined in one
        array when the buffer is read: appending each packet to the
        stream would copy its whole data every time.
        """
        if self._buffer.max_duration is None:
            self._pending.setdefault(record.sid, []).append(record)
            return
        try:
            self._buffer.append(record)
        except ValueError as error:
            logger.warning("record of %s not buffered: %s", record.sid, error)

    def _flush(self) -> None:
        """Move the waiting records into the buffer (lock held)."""
        if not self._pending or self._buffer is None:
            return
        for records in self._pending.values():
            for record in _join_contiguous(records):
                self._buffer.append(record)
        self._pending.clear()

    def _assign_metadata(self, record: Record) -> None:
        """Assign the channel metadata (one warning per channel)."""
        try:
            _set_metadata(record, self._metadata, self._metadata_cache)
        except KeyError:
            if record.sid not in self._missing_metadata:
                self._missing_metadata.add(record.sid)
                warnings.warn(
                    f"no metadata for {record.sid}: its records keep "
                    "empty metadata",
                    _metadata.MetadataWarning,
                    stacklevel=2,
                )

    def _check_continuity(self, record: Record, start_ns: int) -> None:
        """Report a gap or an overlap with the previous record."""
        expected = self._next_start.get(record.sid)
        if expected is None:
            return
        offset = (start_ns - expected) / (record.delta * _NS)
        if abs(offset) <= MERGE_TOLERANCE:
            return
        kind = "gap" if offset > 0 else "overlap"
        self._emit(
            kind,
            message=f"{record.sid}: {kind} of {abs(offset):.1f} samples",
            sid=record.sid,
            start=Date.from_unix_ns(expected),
            end=record.time,
        )

    # -- callbacks -------------------------------------------------------

    def _emit(self, kind: str, **fields: Any) -> None:
        """Report an event to the logger and to the callbacks."""
        event = SeedLinkEvent(kind=kind, time=Date.now(), **fields)
        logger.debug("event %s: %s", kind, event.message)
        for func in self._event_callbacks:
            func(self, event)

    def _run_wall_timers(self, now: float) -> None:
        """Call the wall-clock timers that are due."""
        for timer in self._timers:
            if timer.clock != "wall":
                continue
            if timer.due is None:
                timer.due = now + timer.interval
            elif now >= timer.due:
                while timer.due <= now:
                    timer.due += timer.interval
                timer.func(self)

    def _run_data_timers(self, end_ns: int) -> None:
        """Call the data-time timers that are due."""
        data_time = end_ns / _NS
        for timer in self._timers:
            if timer.clock != "data":
                continue
            if timer.due is None:
                timer.due = data_time + timer.interval
            elif data_time >= timer.due:
                while timer.due <= data_time:
                    timer.due += timer.interval
                timer.func(self)


# ----------------------------------------------------------------------------
# Convenience function
# ----------------------------------------------------------------------------


def read(
    host: str,
    select: Any,
    duration: Optional[Real] = None,
    start: Any = None,
    end: Any = None,
    port: int = DEFAULT_PORT,
    protocol: str = "auto",
    metadata: Any = None,
    timeout: Real = 600.0,
) -> StreamCollection:
    """
    Collect data from a SeedLink server.

    Parameters
    ----------
    host : str
        Server host name or address.
    select : str or sequence of str
        Channel patterns (see :meth:`Client.select`).
    duration : float, optional
        Wall-clock seconds of data to collect from now (or from
        ``start``).
    start, end : Date or str, optional
        Time window. The server must still hold the past data.
    port, protocol, metadata, timeout
        See :class:`Client`.

    Returns
    -------
    StreamCollection
        The data received, trimmed to the time window if given.

    Raises
    ------
    ValueError
        If neither ``duration`` nor ``end`` is given.

    Examples
    --------
    >>> sc = read("rtserve.example.org", "IV.MODE.*.HH?", duration=60)
    >>> sc = read("rtserve.example.org", ["IV.MODE.*.HHZ"],
    ...           start="2026-11-05T08:00:00", end="2026-11-05T08:10:00")
    """
    if duration is None and end is None:
        raise ValueError("duration or end is required")
    client = Client(
        host, port=port, protocol=protocol, metadata=metadata,
        timeout=timeout,
    )
    client.select(select)
    client.set_buffer()
    collection = client.run(duration=duration, start=start, end=end)
    if start is not None or end is not None:
        collection = collection.extract(
            starttime=_as_date(start, "start"),
            endtime=_as_date(end, "end"),
        )
    return collection


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _parse_pattern(pattern: Any) -> _Pattern:
    """Parse a channel pattern ``NET.STA.LOC.CHA``."""
    if not isinstance(pattern, str):
        raise TypeError("channel patterns must be strings")
    text = pattern.strip()
    parts = text.split(".")
    if not text or len(parts) > 4:
        raise ValueError(
            f"invalid channel pattern {pattern!r}: expected NET.STA.LOC.CHA"
        )
    parts += ["*"] * (4 - len(parts))

    for index, part in enumerate(parts):
        if not part and index != 2:
            raise ValueError(
                f"invalid channel pattern {pattern!r}: empty code"
            )
        if any(char.isspace() for char in part):
            raise ValueError(
                f"invalid channel pattern {pattern!r}: spaces"
            )
    return _Pattern(*parts)


def _has_joker(text: str) -> bool:
    return any(char in text for char in "*?[")


def _v3_stations(xml_text: str) -> list:
    """``(network, station)`` pairs of a version 3 INFO STATIONS."""
    from shakelab.libutils.xml import read_xml

    root = read_xml(xml_text.strip())
    return [
        (element.get("network", ""), element.get("name", ""))
        for element in root.iter("station")
    ]


def _decode(packet: Packet) -> tuple:
    """
    Records of a data packet, and the station ``NET_STA``.

    Invalid packets are logged and skipped.
    """
    if packet.format == "2":
        from shakelab.signals.libio.mseed import MSRecord

        try:
            msrecord = MSRecord.from_bytes(packet.payload)
        except ValueError as error:
            logger.warning("invalid MiniSEED packet skipped: %s", error)
            return [], packet.station
        station = packet.station or "_".join((
            msrecord.header["NETWORK_CODE"].strip(),
            msrecord.header["STATION_CODE"].strip(),
        ))
        if not msrecord.is_waveform:
            return [], station
        return [msrecord.to_shakelab()], station

    if packet.format == "3":
        try:
            from shakelab.signals.libio import cymseed
        except ImportError:
            raise NotImplementedError(
                "MiniSEED 3 packets need the cymseed3 package"
            ) from None
        collection = cymseed.msread(packet.payload)
        records = [record for stream in collection for record in stream]
        return records, packet.station

    logger.debug("packet of format %s skipped", packet.format)
    return [], packet.station


def _join_contiguous(records: list) -> Iterator[Record]:
    """Join runs of contiguous records into single records."""
    group = [records[0]]
    for record in records[1:]:
        if _contiguous(group[-1], record):
            group.append(record)
        else:
            yield _concatenate(group)
            group = [record]
    yield _concatenate(group)


def _contiguous(first: Record, second: Record) -> bool:
    """Whether ``second`` continues ``first`` without gap or overlap."""
    delta = first.delta
    if (
        not np.isclose(delta, second.delta, rtol=1e-9, atol=0.0)
        or first.data.dtype != second.data.dtype
        or not _compatible_meta(first.head.meta, second.head.meta)
    ):
        return False
    expected = first.time.to_unix_ns() + first.nsamp * delta * _NS
    offset = (second.time.to_unix_ns() - expected) / (delta * _NS)
    return abs(offset) <= MERGE_TOLERANCE


def _concatenate(group: list) -> Record:
    """Single record with the samples of contiguous records."""
    if len(group) == 1:
        return group[0]
    record = group[0].copy()
    record.data = np.concatenate([item.data for item in group])
    return record


def _as_metadata(metadata: Any) -> Any:
    """Metadata collection from a collection, path or paths."""
    if metadata is None:
        return None
    if isinstance(metadata, _metadata.MetadataCollection):
        return metadata
    from shakelab.signals.io import read_metadata

    return read_metadata(metadata)


def _as_date(value: Any, name: str) -> Optional[Date]:
    """Date from a Date, an ISO string or None."""
    if value is None or isinstance(value, Date):
        return value
    try:
        return Date(value)
    except (TypeError, ValueError) as error:
        raise type(error)(f"invalid {name}: {error}") from None


def _positive(value: Any, name: str) -> float:
    """Validate a positive finite number of seconds."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a number")
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _require_callable(func: Any) -> None:
    if not callable(func):
        raise TypeError("the callback must be callable")
