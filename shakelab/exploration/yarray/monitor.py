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

A monitor processes, at regular intervals, the latest data of an array
arriving through the SeedLink client
(:mod:`shakelab.signals.seedlink`). It is a function of the client,
called by :meth:`Client.every <shakelab.signals.seedlink.Client.every>`:

.. code-block:: python

    client = seedlink.Client("server", metadata=meta)
    client.select("XX.A*..HH?")
    client.set_buffer(max_duration=900)
    monitor = BeamformingMonitor(config, window=600)
    client.every(600, monitor)
    client.start()
    ...
    monitor.status(), monitor.latest(), monitor.result(start, end)

There is one monitor for each processing: :class:`BeamformingMonitor`
(:class:`~.beamforming.Beamformer`, results
:class:`~.beamforming.BeamformingResult`) and :class:`FDDMonitor`
(:class:`~.fdd.FDD`, results :class:`~.fdd.FDDResult`). Both derive
from :class:`ArrayMonitor`, which holds the buffer window, the worker
thread, the callbacks and the queries, and is not used directly.

Design choices
--------------
The client calls its functions in the thread that receives the data,
so a long computation there would stop the reception. A call of the
monitor only copies the data of the buffer
(:meth:`Client.get_data`) and queues them; the processing runs in a
worker thread of the monitor. The queue is short (one pending analysis
by default): when the processing is slower than the calls, the oldest
pending analysis is dropped and reported (:meth:`ArrayMonitor.on_skip`),
so that the monitor never accumulates delay. With ``background=False``
the processing runs in the call itself (light configurations, tests,
replays of archives with ``clock="data"``).

All the query methods (:meth:`~ArrayMonitor.status`,
:meth:`~ArrayMonitor.latest`, :meth:`~ArrayMonitor.results`,
:meth:`~ArrayMonitor.result`, :meth:`~ArrayMonitor.array`) can be
called from any thread while the monitor runs: they return copies.

Windows whose data cannot be processed (no data, too few stations,
too short, singular matrices, unknown positions for the beamforming)
are skipped and reported; other errors stop the worker, are logged and
raised again by :meth:`~ArrayMonitor.close`.

The results of a frequency domain decomposition can be joined only
when they have the same channels: an :class:`FDDMonitor` keeps the
stations of its first result (or those given with ``stations``) and
skips the windows in which one of them is missing; by default it also
skips the windows in which a station has gaps (``gaps="raise"``).
"""

from __future__ import annotations

import json
import logging
import queue
import threading
import time
from collections import deque
from typing import Any, Callable, Optional, Union

import numpy as np

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
)
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.fdd import FDD, FDDResult
from shakelab.exploration.yarray.io import result_path
from shakelab.libutils.timeN import Date
from shakelab.signals.base import StreamCollection

__all__ = ["ArrayMonitor", "BeamformingMonitor", "FDDMonitor"]

_LOGGER = logging.getLogger(__name__)

# Errors of a window that are reported as skipped analyses.
_SKIPPED_ERRORS = (ValueError, np.linalg.LinAlgError)


# ----------------------------------------------------------------------------
# Base monitor
# ----------------------------------------------------------------------------


class ArrayMonitor:
    """
    Periodic array processing of the data of a SeedLink client.

    Base of :class:`BeamformingMonitor` and :class:`FDDMonitor`, which
    are the classes to use; the parameters below are common to both.

    Parameters
    ----------
    config : configuration, str or path-like
        Configuration of the processing of each window, or a JSON file
        (the monitor keeps a copy).
    window : float
        Duration of the analysed data, in seconds, ending at the latest
        time covered by all the selected channels.
    select : str or list of str, optional
        Channels of the buffer (patterns ``NET.STA.LOC.CHA``). Default:
        every channel.
    gaps : {"drop", "raise"}, optional
        Stations with gaps in the window: excluded, or the window is
        skipped.
    min_stations : int, optional
        Fewest stations for an analysis. Default is 3.
    output : str or path-like, optional
        Result file (JSON Lines) to which each result is appended (see
        the ``write`` method of the results; ``.json`` is refused).
    arrays : bool, optional
        Write also the arrays of each result to ``output`` (spectra of
        the beamforming, all singular values and kept vectors of the
        FDD). Default True.
    log : str or path-like, optional
        File (JSON Lines) to which a line is appended for each
        analysis, skipped window or error: ``time`` (wall clock),
        ``event`` (``"result"``, ``"skip"``, ``"error"``), ``window``,
        and ``duration``, ``stations``, ``excluded`` or ``reason``.
    keep : int, optional
        Results kept in memory (the latest). Default: all.
    background : bool, optional
        Process in a worker thread (default) or in the call.
    queue_size : int, optional
        Analyses waiting for the worker. Default is 1.
    coordinates, reference, orientations, coordinate_covariance : optional
        Passed to :class:`ArrayData` (positions and orientations
        otherwise come from the metadata of the records).

    Raises
    ------
    TypeError
        If the class is used directly, or ``config`` is neither a
        configuration nor a file.
    ValueError
        If the configuration is incomplete, ``window``,
        ``min_stations``, ``keep`` or ``queue_size`` are not positive,
        or ``output`` or ``log`` end in ``.json``.
    """

    # Class of the results (set by the derived classes).
    _result_type: type = type(None)

    def __init__(
        self,
        config: Any,
        window: float,
        select: Any = None,
        gaps: str = "drop",
        min_stations: int = 3,
        output: Any = None,
        arrays: bool = True,
        log: Any = None,
        keep: Optional[int] = None,
        background: bool = True,
        queue_size: int = 1,
        coordinates: Any = None,
        reference: Any = None,
        orientations: Any = None,
        coordinate_covariance: Any = None,
    ) -> None:
        if type(self) is ArrayMonitor:
            raise TypeError(
                "ArrayMonitor is the base of the monitors: use "
                "BeamformingMonitor or FDDMonitor"
            )

        processor = self._processor_for(config)
        processor.config.check_complete()

        if not float(window) > 0.0:
            raise ValueError("window must be positive")

        for name, value in (("min_stations", min_stations),
                            ("queue_size", queue_size),
                            ("keep", 1 if keep is None else keep)):
            if int(value) < 1:
                raise ValueError(f"{name} must be at least 1")

        if gaps not in ("drop", "raise"):
            raise ValueError("gaps must be 'drop' or 'raise'")

        output = None if output is None else result_path(output)
        log = None if log is None else result_path(log)

        self._processor = processor
        self.window = float(window)
        self.select = select
        self.gaps = gaps
        self.min_stations = int(min_stations)
        self.output = output
        self.arrays = bool(arrays)
        self.log = log
        self.background = bool(background)
        self._array_options = {
            "coordinates": coordinates,
            "reference": reference,
            "orientations": orientations,
            "coordinate_covariance": coordinate_covariance,
        }
        self._lock = threading.Lock()
        self._log_lock = threading.Lock()
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
            f"{type(self).__name__}(window={self.window:g} s, "
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

        first = min(
            stream.starttime for stream in collection if stream.records()
        )

        if first > window[0]:
            self._skip(
                f"the buffer covers {window[1] - first:.2f} s, the window "
                f"is {self.window:g} s", window,
            )
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
                    "processing is slower than the calls", dropped,
                )

    def process(
        self,
        collection: StreamCollection,
        starttime: Any = None,
        endtime: Any = None,
    ) -> Any:
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
        BeamformingResult, FDDResult or None
            Result, or None if the window was skipped.
        """
        if starttime is None or endtime is None:
            window = self._window(collection)

            if window is None:
                self._skip("no data")
                return None

            starttime, endtime = window

        begin = time.perf_counter()
        label = f"window {starttime} - {endtime}"
        window = (Date(starttime), Date(endtime))

        try:
            array = ArrayData(
                self._waveforms(collection), select=self.select,
                starttime=starttime, endtime=endtime, gaps=self.gaps,
                **self._array_options,
            )

            if array.nsta < self.min_stations:
                self._skip(
                    f"{label}: {array.nsta} stations (at least "
                    f"{self.min_stations} required)", window,
                )
                return None

            array, reason = self._prepare(array)

            if reason is not None:
                self._skip(f"{label}: {reason}", window)
                return None

            self._processor.array = array
            result = self._processor.process(gaps=self.gaps)

        except _SKIPPED_ERRORS as error:
            self._skip(f"{label}: {error}", window)
            return None

        reason = self._check(result)

        if reason is not None:
            self._skip(f"{label}: {reason}", window)
            return None

        if self.output is not None:
            result.write(self.output, append=True, arrays=self.arrays)

        last = {
            "window": window,
            "duration": time.perf_counter() - begin,
            "stations": array.nsta,
            "excluded": array.excluded,
        }

        with self._lock:
            self._results.append((*window, result))
            self._array = array
            self._counts["analyses"] += 1
            self._last = last

        self._log("result", **last)

        self._accepted(result)

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
    def config(self) -> Any:
        """Copy of the configuration of the processing."""
        return self._processor.config.copy()

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

    def latest(self) -> Any:
        """
        Return the latest result.

        Returns
        -------
        BeamformingResult, FDDResult or None
            Latest result, None before the first analysis.
        """
        with self._lock:
            latest = self._results[-1][2] if self._results else None

        return None if latest is None else latest.copy()

    def results(self, starttime: Any = None, endtime: Any = None) -> list:
        """
        Return the results kept, optionally within a time interval.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Interval; results whose window overlaps it are returned.

        Returns
        -------
        list
            Results in time order of the analyses.
        """
        start = None if starttime is None else Date(starttime)
        end = None if endtime is None else Date(endtime)

        with self._lock:
            items = list(self._results)

        return [
            result.copy() for first, last, result in items
            if (start is None or last >= start)
            and (end is None or first <= end)
        ]

    def result(self, starttime: Any = None, endtime: Any = None) -> Any:
        """
        Return the results kept, joined into one.

        Parameters
        ----------
        starttime, endtime : Date or str, optional
            Interval (see :meth:`results`).

        Returns
        -------
        BeamformingResult, FDDResult or None
            Joined result (``append`` method of the results), None if
            there is none.
        """
        items = self.results(starttime, endtime)

        if not items:
            return None

        joined = self._result_type()
        joined.append(items)
        return joined

    def array(self) -> Optional[ArrayData]:
        """
        Return the data of the latest analysis.

        Returns
        -------
        ArrayData or None
            Copy of the array, None before the first analysis.
        """
        with self._lock:
            array = self._array

        return None if array is None else array.copy()

    # -- processing steps (derived classes) --------------------------------

    def _processor_for(self, config: Any) -> Any:
        """Processor of the monitor, with a copy of the configuration."""
        raise NotImplementedError

    def _waveforms(self, collection: StreamCollection) -> Any:
        """Waveforms of the buffer to process."""
        return collection

    def _prepare(self, array: ArrayData) -> tuple:
        """Array to process and the reason to skip it (or None)."""
        return array, None

    def _check(self, result: Any) -> Optional[str]:
        """Reason to skip a result (None to accept it)."""
        return None

    def _accepted(self, result: Any) -> None:
        """Record an accepted result."""

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

    def _skip(self, reason: str, window: Optional[tuple] = None) -> None:
        """Count, log and report a skipped analysis."""
        with self._lock:
            self._counts["skipped"] += 1

        _LOGGER.warning("array analysis skipped: %s", reason)
        self._log("skip", window=window, reason=reason)

        for func in self._on_skip:
            func(self, reason)

    def _start_worker(self) -> None:
        """Start the worker thread if needed."""
        if self._thread is None or not self._thread.is_alive():
            self._thread = threading.Thread(
                target=self._work, name=type(self).__name__, daemon=True
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

                self._log("error", window=(start, end), reason=repr(error))
                return

    def _log(self, event: str, window: Any = None, **fields: Any) -> None:
        """Append an event to the log file, if any."""
        if self.log is None:
            return

        record = {
            "time": Date.now().to_iso8601(precision=6),
            "event": event,
            "window": None if window is None else [
                Date(value).to_iso8601(precision=9) for value in window
            ],
            **fields,
        }

        if "excluded" in record:
            record["excluded"] = [list(item) for item in record["excluded"]]

        with self._log_lock:
            self.log.parent.mkdir(parents=True, exist_ok=True)

            with self.log.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(record) + "\n")

    def _raise_error(self) -> None:
        """Raise the error that stopped the worker, if any."""
        if self._error is not None:
            raise RuntimeError(
                f"the array monitor stopped: {self._error}"
            ) from self._error


# ----------------------------------------------------------------------------
# Beamforming monitor
# ----------------------------------------------------------------------------


class BeamformingMonitor(ArrayMonitor):
    """
    Periodic beamforming of the data of a SeedLink client.

    Each window is processed by a :class:`~.beamforming.Beamformer`
    (with the cache of the steering vectors); the results are
    :class:`~.beamforming.BeamformingResult`. The positions of the
    stations are required: windows with stations of unknown position
    are skipped.

    Parameters
    ----------
    config : BeamformingConfig, str or path-like
        Configuration of the beamforming, or a JSON file (the monitor
        keeps a copy).
    window : float
        Duration of the analysed data, in seconds.
    **options
        Options of :class:`ArrayMonitor` (``select``, ``gaps``, default
        ``"drop"``, ``min_stations``, ``output``, ``arrays``, ``log``,
        ``keep``, ``background``, ``queue_size``, ``coordinates``,
        ``reference``, ``orientations``, ``coordinate_covariance``).

    Raises
    ------
    TypeError
        If ``config`` is neither a BeamformingConfig nor a file.
    ValueError
        If the configuration is incomplete or an option is invalid.

    Examples
    --------
    >>> monitor = BeamformingMonitor(config, window=600, output="bf.jsonl")
    >>> client.every(600, monitor)
    """

    _result_type = BeamformingResult

    def __init__(self, config: Any, window: float, **options: Any) -> None:
        super().__init__(config, window, **options)

    def _processor_for(self, config: Any) -> Beamformer:
        """Beamformer with a copy of the configuration."""
        return Beamformer(config=config, cache=True)

    def _check(self, result: BeamformingResult) -> Optional[str]:
        """Skip the windows too short for one covariance block."""
        if not any(result.block_times.values()):
            return "too short for one covariance block"

        return None


# ----------------------------------------------------------------------------
# FDD monitor
# ----------------------------------------------------------------------------


class FDDMonitor(ArrayMonitor):
    """
    Periodic frequency domain decomposition of the data of a client.

    Each window is processed by an :class:`~.fdd.FDD`; the results are
    :class:`~.fdd.FDDResult`. Results can be joined (:meth:`result`, the
    ``output`` file) only with the same channels: the monitor processes
    always the same stations, those given with ``stations`` or else
    those of its first result. Other stations of the buffer are
    ignored; a window in which one of them is missing (or, with the
    default ``gaps="raise"``, has gaps) is skipped. Positions are not
    required.

    Parameters
    ----------
    config : FDDConfig, str or path-like
        Configuration of the decomposition, or a JSON file (the monitor
        keeps a copy).
    window : float
        Duration of the analysed data, in seconds.
    stations : str or sequence of str, optional
        Stations to process (``NET.STA.LOC``), in this order. Default:
        those of the first result.
    gaps : {"raise", "drop"}, optional
        Stations with gaps in the window: the window is skipped
        (default) or they are excluded (then the window is skipped
        anyway if they are among the stations processed).
    **options
        Other options of :class:`ArrayMonitor` (``select``,
        ``min_stations``, ``output``, ``arrays``, ``log``, ``keep``,
        ``background``, ``queue_size``, ``coordinates``,
        ``reference``, ``orientations``, ``coordinate_covariance``).

    Raises
    ------
    TypeError
        If ``config`` is neither an FDDConfig nor a file, or
        ``stations`` is not a string or a sequence of strings.
    ValueError
        If the configuration is incomplete, an option is invalid or
        ``stations`` is empty, has an empty name or repeats a station.

    Examples
    --------
    >>> monitor = FDDMonitor(config, window=1800, output="fdd.jsonl")
    >>> client.every(300, monitor)
    >>> monitor.result().mode(1.2).info()
    """

    _result_type = FDDResult

    def __init__(
        self,
        config: Any,
        window: float,
        stations: Union[str, Any, None] = None,
        gaps: str = "raise",
        **options: Any,
    ) -> None:
        super().__init__(config, window, gaps=gaps, **options)
        self._stations: Optional[tuple] = None

        if stations is not None:
            if isinstance(stations, str):
                stations = [stations]

            try:
                stations = tuple(stations)
            except TypeError:
                raise TypeError(
                    "stations must be a string or a sequence of strings"
                ) from None

            if not all(isinstance(sid, str) for sid in stations):
                raise TypeError(
                    "stations must be a string or a sequence of strings"
                )

            stations = tuple(sid.strip() for sid in stations)

            if (
                not stations or not all(stations)
                or len(set(stations)) != len(stations)
            ):
                raise ValueError(
                    "stations must be non-empty and without repetitions"
                )

            self._stations = stations

    @property
    def stations(self) -> Optional[tuple[str, ...]]:
        """Stations processed (None before the first result)."""
        with self._lock:
            return self._stations

    def _processor_for(self, config: Any) -> FDD:
        """FDD processor with a copy of the configuration."""
        return FDD(config=config)

    def _waveforms(self, collection: StreamCollection) -> Any:
        """Records of the stations processed (all before the first)."""
        stations = self.stations

        if stations is None:
            return collection

        wanted = set(stations)
        return [
            record for stream in collection
            if stream.sid is not None
            and stream.sid.rsplit(".", 1)[0] in wanted
            for record in stream.records()
        ]

    def _prepare(self, array: ArrayData) -> tuple:
        """Array of the stations processed, in their order."""
        stations = self.stations

        if stations is None:
            return array, None

        missing = [sid for sid in stations if sid not in array]

        if missing:
            return array, (
                "stations " + ", ".join(missing) + " missing (the "
                "results keep the same stations)"
            )

        if array.sids == stations:
            return array, None

        return array._derived([array[sid] for sid in stations], []), None

    def _check(self, result: FDDResult) -> Optional[str]:
        """Skip the windows too short for one block."""
        if not result.nblocks:
            return "too short for one block"

        return None

    def _accepted(self, result: FDDResult) -> None:
        """Keep the stations of the first result."""
        with self._lock:
            if self._stations is None:
                self._stations = tuple(dict.fromkeys(result.stations))
