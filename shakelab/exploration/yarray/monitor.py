# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Real-time array processing.

An :class:`ArrayMonitor` processes, at regular intervals, the latest
data of an array arriving through the SeedLink client
(:mod:`shakelab.signals.seedlink`). It is a function of the client,
called by :meth:`Client.every <shakelab.signals.seedlink.Client.every>`:

.. code-block:: python

    client = seedlink.Client("server", metadata=meta)
    client.select("XX.A*..HH?")
    client.set_buffer(max_duration=900)
    monitor = ArrayMonitor(Beamformer(config, cache=True), window=600)
    client.every(600, monitor)
    client.start()
    ...
    monitor.status(), monitor.latest(), monitor.result(start, end)

Design choices
--------------
The client calls its functions in the thread that receives the data,
so a long computation there would stop the reception. A call of the
monitor only copies the data of the buffer
(:meth:`Client.get_data`) and queues them; the processing runs in a
worker thread of the monitor. The queue is short (one pending analysis
by default): when the processing is slower than the calls, the oldest
pending analysis is dropped and reported (:meth:`on_skip`), so that
the monitor never accumulates delay. With ``background=False`` the
processing runs in the call itself (light configurations, tests,
replays of archives with ``clock="data"``).

All the query methods (:meth:`status`, :meth:`latest`, :meth:`results`,
:meth:`result`, :meth:`array`) can be called from any thread while the
monitor runs: they return immutable objects or copies.

Windows whose data cannot be processed (no data, too few stations,
too short, singular matrices) are skipped and reported; other errors
stop the worker, are logged and raised again by :meth:`close`.
"""

from __future__ import annotations

import logging
import queue
import threading
import time
from collections import deque
from typing import Any, Callable, Optional

import numpy as np

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
)
from shakelab.exploration.yarray.config import BeamformingConfig
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.io import write_result
from shakelab.libutils.timeN import Date
from shakelab.signals.base import StreamCollection

__all__ = ["ArrayMonitor"]

_LOGGER = logging.getLogger(__name__)

# Errors of a window that are reported as skipped analyses.
_SKIPPED_ERRORS = (ValueError, np.linalg.LinAlgError)


class ArrayMonitor:
    """
    Periodic array processing of the data of a SeedLink client.

    Parameters
    ----------
    beamformer : Beamformer or BeamformingConfig
        Processing applied to each window (a configuration is wrapped
        in a Beamformer with cache).
    window : float
        Duration of the analysed data, in seconds, ending at the latest
        time covered by all the selected channels.
    select : str or list of str, optional
        Channels of the buffer (patterns ``NET.STA.LOC.CHA``). Default:
        every channel.
    gaps : {"drop", "raise"}, optional
        Stations with gaps in the window: excluded (default) or the
        window is skipped.
    min_stations : int, optional
        Fewest stations for an analysis. Default is 3.
    output : str or path-like, optional
        File to which each result is appended (pickle, see
        :func:`~.io.write_result`).
    keep : int, optional
        Results kept in memory (the latest). Default: all.
    background : bool, optional
        Process in a worker thread (default) or in the call.
    queue_size : int, optional
        Analyses waiting for the worker. Default is 1.
    coordinates, reference, orientations, coordinate_covariance : optional
        Passed to :meth:`ArrayData.from_collection` (positions and
        orientations otherwise come from the metadata of the records).

    Raises
    ------
    TypeError
        If ``beamformer`` is neither a Beamformer nor a configuration.
    ValueError
        If ``window``, ``min_stations``, ``keep`` or ``queue_size`` are
        not positive.
    """

    def __init__(
        self,
        beamformer: Any,
        window: float,
        select: Any = None,
        gaps: str = "drop",
        min_stations: int = 3,
        output: Any = None,
        keep: Optional[int] = None,
        background: bool = True,
        queue_size: int = 1,
        coordinates: Any = None,
        reference: Any = None,
        orientations: Any = None,
        coordinate_covariance: Any = None,
    ) -> None:
        if isinstance(beamformer, BeamformingConfig):
            beamformer = Beamformer(beamformer, cache=True)

        if not isinstance(beamformer, Beamformer):
            raise TypeError("beamformer must be a Beamformer or a config")

        if not float(window) > 0.0:
            raise ValueError("window must be positive")

        for name, value in (("min_stations", min_stations),
                            ("queue_size", queue_size),
                            ("keep", 1 if keep is None else keep)):
            if int(value) < 1:
                raise ValueError(f"{name} must be at least 1")

        self.beamformer = beamformer
        self.window = float(window)
        self.select = select
        self.gaps = gaps
        self.min_stations = int(min_stations)
        self.output = output
        self.background = bool(background)
        self._array_options = {
            "coordinates": coordinates,
            "reference": reference,
            "orientations": orientations,
            "coordinate_covariance": coordinate_covariance,
        }
        self._lock = threading.Lock()
        self._results: deque = deque(maxlen=keep)
        self._queue: queue.Queue = queue.Queue(maxsize=int(queue_size))
        self._thread: Optional[threading.Thread] = None
        self._on_result: list = []
        self._on_skip: list = []
        self._array: Optional[ArrayData] = None
        self._error: Optional[BaseException] = None
        self._counts = {"calls": 0, "analyses": 0, "skipped": 0}
        self._last: dict = {}

    def __repr__(self) -> str:
        return (
            f"ArrayMonitor(window={self.window:g} s, "
            f"{self._counts['analyses']} analyses)"
        )

    def __enter__(self) -> ArrayMonitor:
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    # -- callbacks ---------------------------------------------------------

    def on_result(self, func: Callable) -> Callable:
        """
        Call a function after each analysis.

        Parameters
        ----------
        func : callable
            Called as ``func(monitor, result)`` in the processing thread.

        Returns
        -------
        callable
            ``func`` itself (usable as a decorator).
        """
        self._on_result.append(func)
        return func

    def on_skip(self, func: Callable) -> Callable:
        """
        Call a function when an analysis is skipped.

        Parameters
        ----------
        func : callable
            Called as ``func(monitor, reason)`` with a text.

        Returns
        -------
        callable
            ``func`` itself (usable as a decorator).
        """
        self._on_skip.append(func)
        return func

    def on_event(self, client: Any, event: Any) -> None:
        """
        Log the connection events of the client.

        Pass it to :meth:`Client.on_event`.

        Parameters
        ----------
        client : Client
            SeedLink client.
        event : SeedLinkEvent
            Event.
        """
        _LOGGER.info("seedlink event: %s", event)

    # -- processing --------------------------------------------------------

    def __call__(self, client: Any) -> None:
        """
        Copy the latest window of the buffer and process it.

        Parameters
        ----------
        client : Client
            SeedLink client (or any object with ``get_data(select)``
            returning a StreamCollection).
        """
        self._raise_error()

        with self._lock:
            self._counts["calls"] += 1

        pattern = self.select if isinstance(self.select, str) else "*"
        collection = client.get_data(pattern)
        window = self._window(collection)

        if window is None:
            self._skip("no data in the buffer")
            return

        if not self.background:
            self.process(collection, *window)
            return

        self._start_worker()

        while True:
            try:
                self._queue.put_nowait((collection, window))
                return
            except queue.Full:
                try:
                    _, dropped = self._queue.get_nowait()
                except queue.Empty:
                    continue

                self._skip(
                    f"window {dropped[0]} - {dropped[1]} dropped: the "
                    "processing is slower than the calls"
                )

    def process(
        self,
        collection: StreamCollection,
        starttime: Any = None,
        endtime: Any = None,
    ) -> Optional[BeamformingResult]:
        """
        Process a window of data now, in the calling thread.

        Parameters
        ----------
        collection : StreamCollection
            Data of the array.
        starttime, endtime : Date or str, optional
            Window; default: the last ``window`` seconds covered by all
            channels.

        Returns
        -------
        BeamformingResult or None
            Result, or None if the window was skipped.
        """
        if starttime is None or endtime is None:
            window = self._window(collection)

            if window is None:
                self._skip("no data")
                return None

            starttime, endtime = window

        begin = time.perf_counter()

        try:
            array = ArrayData.from_collection(
                collection, select=self.select, starttime=starttime,
                endtime=endtime, gaps=self.gaps, **self._array_options,
            )

            if array.nsta < self.min_stations:
                self._skip(
                    f"window {starttime} - {endtime}: {array.nsta} "
                    f"stations (at least {self.min_stations} required)"
                )
                return None

            result = self.beamformer.process(array)

        except _SKIPPED_ERRORS as error:
            self._skip(f"window {starttime} - {endtime}: {error}")
            return None

        if not any(result.block_times.values()):
            self._skip(
                f"window {starttime} - {endtime}: too short for one "
                "covariance block"
            )
            return None

        if self.output is not None:
            write_result(result, self.output, append=True)

        with self._lock:
            self._results.append((Date(starttime), Date(endtime), result))
            self._array = array
            self._counts["analyses"] += 1
            self._last = {
                "window": (Date(starttime), Date(endtime)),
                "duration": time.perf_counter() - begin,
                "stations": array.nsta,
                "excluded": array.excluded,
            }

        for func in self._on_result:
            func(self, result)

        return result

    def close(self, wait: bool = True) -> None:
        """
        Stop the worker thread.

        Parameters
        ----------
        wait : bool, optional
            Wait for the pending analyses (default) or drop them.

        Raises
        ------
        RuntimeError
            If the worker stopped on an unexpected error.
        """
        thread = self._thread

        if thread is not None and thread.is_alive():
            if not wait:
                while not self._queue.empty():
                    try:
                        self._queue.get_nowait()
                    except queue.Empty:
                        break

            self._queue.put(None)
            thread.join()

        self._thread = None
        self._raise_error()

    # -- queries -----------------------------------------------------------

    @property
    def running(self) -> bool:
        """Whether the worker thread is running."""
        return self._thread is not None and self._thread.is_alive()

    def status(self) -> dict:
        """
        Return the state of the monitor.

        Returns
        -------
        dict
            ``running``, ``calls``, ``analyses``, ``skipped``,
            ``pending`` (queued analyses), ``last_window``,
            ``last_duration`` (seconds), ``stations``, ``excluded``
            (stations left out of the last analysis) and ``error``.
        """
        with self._lock:
            return {
                "running": self.running,
                **self._counts,
                "pending": self._queue.qsize(),
                "last_window": self._last.get("window"),
                "last_duration": self._last.get("duration"),
                "stations": self._last.get("stations"),
                "excluded": self._last.get("excluded", ()),
                "error": None if self._error is None else str(self._error),
            }

    def latest(self) -> Optional[BeamformingResult]:
        """
        Return the latest result.

        Returns
        -------
        BeamformingResult or None
            Latest result, None before the first analysis.
        """
        with self._lock:
            return self._results[-1][2] if self._results else None

    def results(self, starttime: Any = None, endtime: Any = None) -> list:
        """
        Return the results kept, optionally within a time interval.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Interval; results whose window overlaps it are returned.

        Returns
        -------
        list of BeamformingResult
            Results in time order of the analyses.
        """
        start = None if starttime is None else Date(starttime)
        end = None if endtime is None else Date(endtime)

        with self._lock:
            items = list(self._results)

        return [
            result for first, last, result in items
            if (start is None or last >= start)
            and (end is None or first <= end)
        ]

    def result(
        self,
        starttime: Any = None,
        endtime: Any = None,
    ) -> Optional[BeamformingResult]:
        """
        Return the results kept, joined into one.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Interval (see :meth:`results`).

        Returns
        -------
        BeamformingResult or None
            Joined result (:meth:`BeamformingResult.concatenate`), None
            if there is none.
        """
        items = self.results(starttime, endtime)
        return BeamformingResult.concatenate(items) if items else None

    def array(self) -> Optional[ArrayData]:
        """
        Return the data of the latest analysis.

        Returns
        -------
        ArrayData or None
            Immutable snapshot, None before the first analysis.
        """
        with self._lock:
            return self._array

    # -- helpers -----------------------------------------------------------

    def _window(self, collection: StreamCollection) -> Optional[tuple]:
        """Last ``window`` seconds covered by all selected channels."""
        ends = [
            stream.endtime for stream in collection
            if stream.records()
        ]

        if not ends:
            return None

        end = min(ends)
        return end - self.window, end

    def _skip(self, reason: str) -> None:
        """Count, log and report a skipped analysis."""
        with self._lock:
            self._counts["skipped"] += 1

        _LOGGER.warning("array analysis skipped: %s", reason)

        for func in self._on_skip:
            func(self, reason)

    def _start_worker(self) -> None:
        """Start the worker thread if needed."""
        if self._thread is None or not self._thread.is_alive():
            self._thread = threading.Thread(
                target=self._work, name="ArrayMonitor", daemon=True
            )
            self._thread.start()

    def _work(self) -> None:
        """Process the queued windows until the stop marker."""
        while True:
            item = self._queue.get()

            if item is None:
                return

            collection, (start, end) = item

            try:
                self.process(collection, start, end)
            except Exception as error:  # stop the worker, see close()
                _LOGGER.exception("array monitor stopped")

                with self._lock:
                    self._error = error

                return

    def _raise_error(self) -> None:
        """Raise the error that stopped the worker, if any."""
        if self._error is not None:
            raise RuntimeError(
                f"the array monitor stopped: {self._error}"
            ) from self._error
