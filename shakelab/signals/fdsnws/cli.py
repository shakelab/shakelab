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
Command-line client of the FDSN web services.

Subcommands:

- ``waveforms``: MiniSEED records of a time window, written as
  received, or corrected for the instrument response (``--correct``)
  or filtered and written in any format of :mod:`shakelab.signals.io`;
  optionally with the StationXML of the channels (``--metadata``);
- ``metadata``: StationXML or text written to a file, or a table of the
  channels or stations;
- ``events``: QuakeML written to a file, or a table of the events;
- ``availability``: time spans of the available data;
- ``datacenters``: the known data centers and their services.

The data center is a name of the FDSN registry (``INGV``,
``EarthScope``, ``GEOFON``...) or the ``http(s)://`` address of a host;
channels are patterns ``NET.STA.LOC.CHA`` (see
:mod:`shakelab.signals.codes`); times are ISO 8601 in UTC.

Examples
--------
.. code-block:: shell

    slab-fdsn waveforms INGV "IV.MODE..HN?" \\
        --start 2012-05-20T02:03:23 --end 2012-05-20T02:06:23 \\
        -o mode.mseed --metadata mode.xml
    slab-fdsn waveforms INGV "IV.MODE..HN?" \\
        --start 2012-05-20T02:04:00 --end 2012-05-20T02:05:30 \\
        --correct --lowpass 40 -o "{net}.{sta}.{loc}.{cha}.sac"
    slab-fdsn metadata INGV "IV.MO*" --level station
    slab-fdsn events USGS --start 2012-05-01T00:00:00 \\
        --end 2012-06-30T00:00:00 --minmag 5.5 \\
        --region 9 13 43.5 46
    slab-fdsn datacenters

The command ``slab-fdsn`` is installed with ShakeLab; the same client
runs as ``python -m shakelab.signals.fdsnws``.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from typing import Optional

from shakelab.signals.fdsnws.client import Client
from shakelab.signals.fdsnws.registry import datacenters
from shakelab.signals.fdsnws.transport import DEFAULT_TIMEOUT, FdsnError

__all__ = ["main"]


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------


def main(argv: Optional[list] = None, prog: Optional[str] = None) -> int:
    """
    Run the command-line client.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; ``sys.argv[1:]`` if omitted.
    prog : str, optional
        Program name shown in the help.

    Returns
    -------
    int
        Exit status: 0 on success, 1 on errors or without data, 130
        when interrupted.
    """
    parser = _parser(prog)
    args = parser.parse_args(argv)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            warnings.showwarning = _show_warning
            return args.command(args)
    except KeyboardInterrupt:
        return 130
    except (OSError, FdsnError, ValueError, TypeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


# ----------------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------------


def _waveforms(args) -> int:
    """Download waveforms."""
    client = _client(args)
    processing = (
        args.correct or args.highpass is not None
        or args.lowpass is not None
    )

    if args.metadata is not None:
        size = client.get_metadata(
            args.patterns, args.start, args.end, level="response",
            target=args.metadata, owrite=args.owrite,
        )
        if not size:
            print("no metadata", file=sys.stderr)

    if not processing:
        if args.format not in (None, "mseed"):
            raise ValueError(
                "--format needs --correct or a filter: without them the "
                "records are written as received (MiniSEED)"
            )
        size = client.get_waveforms(
            args.patterns, args.start, args.end, chunk=args.chunk,
            quality=args.quality, target=args.output, owrite=args.owrite,
        )
        return _report(size, f"{size} bytes written to {args.output}")

    collection = client.get_waveforms(
        args.patterns, args.start, args.end,
        correct=args.correct, pad=args.pad, highpass=args.highpass,
        lowpass=args.lowpass,
        prefilter=None if args.prefilter is None else tuple(args.prefilter),
        waterlevel=args.waterlevel, chunk=args.chunk, quality=args.quality,
    )
    records = collection.records()
    if not records:
        return _report(0, "")

    from shakelab.signals.io import writer

    writer(collection, args.output, format=args.format, owrite=args.owrite)
    return _report(len(records), f"{len(records)} record(s) written")


def _metadata(args) -> int:
    """Download metadata, to a file or as a table."""
    client = _client(args)
    patterns = args.patterns or "*"
    if args.output is not None:
        size = client.get_metadata(
            patterns, args.start, args.end, level=args.level,
            text=args.text, target=args.output, owrite=args.owrite,
        )
        return _report(size, f"{size} bytes written to {args.output}")

    collection = client.get_metadata(
        patterns, args.start, args.end, level=args.level, text=args.text,
    )
    if args.level == "station":
        stations = collection.stations()
        for station in stations:
            print(
                f"{station.code:<12} {_point(station.location)}  "
                f"{_time(station.starttime):<27} "
                f"{_time(station.endtime):<27} {station.site or ''}"
            )
        return _report(len(stations), f"{len(stations)} station(s)")

    epochs = list(collection)
    for epoch in epochs:
        rate = "" if epoch.sample_rate is None else f"{epoch.sample_rate:g}"
        print(
            f"{epoch.sid:<16} {_point(epoch.location)}  {rate:>7} Hz  "
            f"{_time(epoch.starttime):<27} {_time(epoch.endtime)}"
        )
    return _report(len(epochs), f"{len(epochs)} channel epoch(s)")


def _events(args) -> int:
    """Download a catalogue, to a file or as a table."""
    client = _client(args)
    options = dict(
        region=None if args.region is None else tuple(args.region),
        center=None if args.center is None else tuple(args.center),
        radius=None if args.maxradius is None else (None, args.maxradius),
        depth=_range(args.mindepth, args.maxdepth),
        magnitude=_range(args.minmag, args.maxmag),
        magnitudetype=args.magtype,
        limit=args.limit,
        orderby=args.orderby,
    )
    if args.output is not None:
        size = client.get_events(
            args.start, args.end, target=args.output, owrite=args.owrite,
            **options,
        )
        return _report(size, f"{size} bytes written to {args.output}")

    catalogue = client.get_events(args.start, args.end, **options)
    for event in catalogue.event:
        location = event.location.prime
        magnitude = event.magnitude.prime
        size = "" if magnitude is None else (
            f"{magnitude['MagSize']:.1f} {magnitude['MagType'] or ''}"
        )
        print(
            f"{str(location.date):<27} {location['Longitude']:>9.4f} "
            f"{location['Latitude']:>8.4f} {location['Depth']:>6.1f} km  "
            f"{size:<8} {event.id}"
        )
    count = len(catalogue.event)
    return _report(count, f"{count} event(s)")


def _availability(args) -> int:
    """List the spans of available data."""
    client = _client(args)
    spans = client.get_availability(
        args.patterns, args.start, args.end, extent=args.extent,
        merge=args.merge, mergegaps=args.mergegaps,
    )
    for span in spans:
        print(
            f"{span.sid:<16} {_time(span.starttime):<27} "
            f"{_time(span.endtime):<27} {span.quality or ''}"
        )
    return _report(len(spans), f"{len(spans)} span(s)")


def _datacenters(args) -> int:
    """List the data centers and their services."""
    centers = datacenters(refresh=args.refresh)
    for name, center in centers.items():
        if args.service and args.service not in center.services:
            continue
        print(f"{name:<12} {', '.join(sorted(center.services))}")
        if args.urls:
            for service, url in sorted(center.services.items()):
                print(f"{'':<12} {service:<12} {url}")
    return 0


# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------


def _parser(prog: Optional[str] = None) -> argparse.ArgumentParser:
    """Parser of the command line."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="FDSN web service client of ShakeLab.",
    )
    commands = parser.add_subparsers(dest="name", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "datacenter",
        help="data center of the FDSN registry (e.g. INGV) or "
        "http(s):// address of a host",
    )
    common.add_argument(
        "--timeout", type=float, default=DEFAULT_TIMEOUT,
        help=f"seconds without data (default: {DEFAULT_TIMEOUT:g})",
    )

    window = argparse.ArgumentParser(add_help=False)
    window.add_argument("--start", help="start time (ISO 8601, UTC)")
    window.add_argument("--end", help="end time (ISO 8601, UTC)")

    output = argparse.ArgumentParser(add_help=False)
    output.add_argument(
        "--owrite", action="store_true", help="replace existing files",
    )

    waveforms = commands.add_parser(
        "waveforms", parents=[common, window, output],
        help="waveforms of a time window",
    )
    waveforms.add_argument(
        "patterns", nargs="+", metavar="PATTERN",
        help="channel patterns NET.STA.LOC.CHA, e.g. 'IV.MODE..HN?'",
    )
    waveforms.add_argument(
        "-o", "--output", required=True,
        help="file, directory or file name template (see "
        "shakelab.signals.io.writer)",
    )
    waveforms.add_argument(
        "--format", help="output format with --correct or filters "
        "(default: from the extension)",
    )
    waveforms.add_argument(
        "--metadata", metavar="FILE",
        help="also write the StationXML of the channels",
    )
    waveforms.add_argument(
        "--correct", action="store_true",
        help="correct the instrument response",
    )
    waveforms.add_argument(
        "--pad", type=float,
        help="margin in seconds (default: 10%% of the window, at least "
        "30 s)",
    )
    waveforms.add_argument("--highpass", type=float, help="Hz")
    waveforms.add_argument("--lowpass", type=float, help="Hz")
    waveforms.add_argument(
        "--prefilter", type=float, nargs=4, metavar=("F1", "F2", "F3", "F4"),
        help="frequency window of the correction, Hz",
    )
    waveforms.add_argument(
        "--waterlevel", type=float, default=60.0,
        help="water level of the correction in dB (default: 60)",
    )
    waveforms.add_argument(
        "--chunk", type=float, metavar="SECONDS",
        help="longest window of a single request",
    )
    waveforms.add_argument("--quality", help="data quality code")
    waveforms.set_defaults(command=_waveforms)

    metadata = commands.add_parser(
        "metadata", parents=[common, window, output],
        help="station and channel metadata",
    )
    metadata.add_argument(
        "patterns", nargs="*", metavar="PATTERN",
        help="channel patterns (default: all)",
    )
    metadata.add_argument(
        "--level", default="channel",
        choices=("station", "channel", "response"),
        help="detail of the metadata (default: channel)",
    )
    metadata.add_argument(
        "--text", action="store_true", help="text format (not response)",
    )
    metadata.add_argument(
        "-o", "--output", help="file where the reply is written",
    )
    metadata.set_defaults(command=_metadata)

    events = commands.add_parser(
        "events", parents=[common, window, output],
        help="earthquake catalogue",
    )
    events.add_argument("--minmag", type=float)
    events.add_argument("--maxmag", type=float)
    events.add_argument("--magtype", help="magnitude type, e.g. Mw")
    events.add_argument("--mindepth", type=float, help="km")
    events.add_argument("--maxdepth", type=float, help="km")
    events.add_argument(
        "--region", type=float, nargs=4,
        metavar=("LONMIN", "LONMAX", "LATMIN", "LATMAX"),
    )
    events.add_argument(
        "--center", type=float, nargs=2, metavar=("LON", "LAT"),
    )
    events.add_argument(
        "--maxradius", type=float, help="degrees from --center",
    )
    events.add_argument("--limit", type=int)
    events.add_argument(
        "--orderby",
        choices=("time", "time-asc", "magnitude", "magnitude-asc"),
    )
    events.add_argument("-o", "--output", help="QuakeML file to write")
    events.set_defaults(command=_events)

    availability = commands.add_parser(
        "availability", parents=[common, window],
        help="time spans of the available data",
    )
    availability.add_argument(
        "patterns", nargs="+", metavar="PATTERN", help="channel patterns",
    )
    availability.add_argument(
        "--extent", action="store_true",
        help="first and last sample of each channel only",
    )
    availability.add_argument(
        "--merge", nargs="+", choices=("samplerate", "quality", "overlap"),
    )
    availability.add_argument(
        "--mergegaps", type=float, metavar="SECONDS",
        help="ignore shorter gaps",
    )
    availability.set_defaults(command=_availability)

    centers = commands.add_parser(
        "datacenters", help="known data centers and their services",
    )
    centers.add_argument(
        "--refresh", action="store_true",
        help="read the FDSN registry instead of the built-in table",
    )
    centers.add_argument(
        "--service",
        choices=("dataselect", "station", "event", "availability"),
        help="only the data centers offering this service",
    )
    centers.add_argument(
        "--urls", action="store_true", help="show the service addresses",
    )
    centers.set_defaults(command=_datacenters)

    return parser


def _client(args) -> Client:
    """Client of the data center of the command line."""
    return Client(args.datacenter, timeout=args.timeout)


def _range(low: Optional[float], high: Optional[float]):
    return None if low is None and high is None else (low, high)


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------


def _report(count: int, message: str) -> int:
    """Summary on standard error; exit status 1 without data."""
    if not count:
        print("no data", file=sys.stderr)
        return 1
    print(message, file=sys.stderr)
    return 0


def _point(point) -> str:
    if point is None:
        return f"{'-':>9} {'-':>8} {'-':>7}"
    return (
        f"{point.longitude:>9.4f} {point.latitude:>8.4f} "
        f"{point.elevation:>7.1f}"
    )


def _time(value) -> str:
    return str(value) if value is not None else "-"


def _show_warning(message, category, filename, lineno, file=None,
                  line=None) -> None:
    print(f"warning: {message}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(prog="python -m shakelab.signals.fdsnws.cli"))
