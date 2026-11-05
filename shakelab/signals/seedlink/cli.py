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
Command-line SeedLink client.

Subcommands:

- ``info``: stations, or data streams with their time span, offered by
  a server;
- ``monitor``: one line per record received, or a table of the channel
  latencies at regular intervals;
- ``get``: the records of a time window (or of a duration from now) in
  a MiniSEED file;
- ``archive``: continuous archive of MiniSEED files (template or SDS),
  with a state file to resume after a restart.

The server is given as ``HOST[:PORT]`` (default port 18000); channels
as patterns ``NET.STA.LOC.CHA`` (see :mod:`shakelab.signals.codes`) or
a file of patterns.

Examples
--------
.. code-block:: shell

    slab-seedlink info HOST:PORT --streams OX.CAE
    slab-seedlink monitor HOST:PORT "OX.*.*.HHZ" \\
        --status 10
    slab-seedlink get HOST:PORT "OX.CAE..HH?" \\
        --start 2026-11-05T12:00:00 --end 2026-11-05T12:05:00 \\
        -o cae.mseed
    slab-seedlink archive HOST:PORT "OX.*.*.HH?" \\
        --sds archive --state ogs.state

The command ``slab-seedlink`` is installed with ShakeLab; the same
client runs as ``python -m shakelab.signals.seedlink``.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Optional

from shakelab.signals.seedlink.client import Client
from shakelab.signals.seedlink.protocol import DEFAULT_PORT, SeedLinkError

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
        Exit status: 0 on success, 1 on errors, 130 when interrupted.
    """
    parser = _parser(prog)
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    try:
        host, port = _server(args.server)
        return args.command(args, host, port)
    except KeyboardInterrupt:
        return 130
    except (OSError, SeedLinkError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


# ----------------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------------


def _info(args, host: str, port: int) -> int:
    """List the stations or the streams of the server."""
    client = _client(args, host, port)
    pattern = args.pattern or "*"
    if args.streams:
        streams = client.streams(pattern)
        for stream in streams:
            print(f"{stream.sid:<16} {_time(stream.start):<28} "
                  f"{_time(stream.end)}")
        _count(streams, "stream")
    else:
        stations = client.stations(pattern)
        for station in stations:
            print(f"{station.code:<10} {station.description}")
        _count(stations, "station")
    return 0


def _monitor(args, host: str, port: int) -> int:
    """Print the records received, or the channel latencies."""
    client = _client(args, host, port)
    _select(client, args)
    if args.status:
        client.every(args.status, _print_status)
    else:
        client.on_record(_print_record)
    client.run(duration=args.duration)
    return 0


def _get(args, host: str, port: int) -> int:
    """Write the records of a time window in a MiniSEED file."""
    if args.duration is None and args.end is None:
        raise ValueError("give --end or --duration")
    if os.path.exists(args.output):
        raise ValueError(f"{args.output} exists")
    client = _client(args, host, port)
    _select(client, args)
    client.set_archive(args.output.replace("{", "{{").replace("}", "}}"))
    client.run(duration=args.duration, start=args.start, end=args.end)
    if not os.path.exists(args.output):
        print("no data received", file=sys.stderr)
        return 1
    return 0


def _archive(args, host: str, port: int) -> int:
    """Archive the records continuously."""
    client = _client(args, host, port)
    _select(client, args)
    if args.sds is not None:
        client.set_archive(args.sds, layout="sds")
    else:
        client.set_archive(args.path, rotate=_rotate(args.rotate))
    client.run(duration=args.duration)
    return 0


# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------


def _parser(prog: Optional[str] = None) -> argparse.ArgumentParser:
    """Parser of the command line."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="SeedLink client of ShakeLab.",
    )
    commands = parser.add_subparsers(dest="name", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "server", help=f"HOST[:PORT] of the server (port {DEFAULT_PORT} "
        "if omitted)",
    )
    common.add_argument(
        "--protocol", default="auto", choices=("auto", "3.1", "4.0"),
        help="protocol version (default: auto)",
    )
    common.add_argument(
        "--timeout", type=float, default=600.0,
        help="seconds without data before reconnecting (default: 600)",
    )
    common.add_argument(
        "-v", "--verbose", action="store_true",
        help="report connections and events",
    )

    channels = argparse.ArgumentParser(add_help=False)
    channels.add_argument(
        "patterns", nargs="*", metavar="PATTERN",
        help="channel patterns NET.STA.LOC.CHA, e.g. 'OX.{CAE,BAD}.*.HHZ'",
    )
    channels.add_argument(
        "--file", help="file of channel patterns (one or more per line)",
    )
    channels.add_argument(
        "--duration", type=float,
        help="stop after this many seconds",
    )

    info = commands.add_parser(
        "info", parents=[common], help="stations or streams of the server",
    )
    info.add_argument("pattern", nargs="?", help="station or channel pattern")
    info.add_argument(
        "--streams", action="store_true",
        help="list the data streams and their time span",
    )
    info.set_defaults(command=_info)

    monitor = commands.add_parser(
        "monitor", parents=[common, channels],
        help="print the records received",
    )
    monitor.add_argument(
        "--status", type=float, metavar="SECONDS",
        help="print a table of the channel latencies at this interval "
        "instead of one line per record",
    )
    monitor.set_defaults(command=_monitor)

    get = commands.add_parser(
        "get", parents=[common, channels],
        help="records of a time window in a MiniSEED file",
    )
    get.add_argument("--start", help="start time (ISO 8601, UTC)")
    get.add_argument("--end", help="end time (ISO 8601, UTC)")
    get.add_argument(
        "-o", "--output", required=True, help="MiniSEED file to write",
    )
    get.set_defaults(command=_get)

    archive = commands.add_parser(
        "archive", parents=[common, channels],
        help="continuous MiniSEED archive",
    )
    target = archive.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--sds", metavar="DIR", help="SeisComP Data Structure root",
    )
    target.add_argument(
        "--path", metavar="TEMPLATE",
        help="path template, e.g. 'data/{net}.{sta}.{loc}.{cha}.mseed'",
    )
    archive.add_argument(
        "--rotate",
        help="file duration with --path: seconds or e.g. 10min, 1h, 1d",
    )
    archive.add_argument(
        "--state", metavar="FILE",
        help="state file to resume after a restart",
    )
    archive.set_defaults(command=_archive)

    return parser


def _server(text: str) -> tuple:
    """Host and port of ``HOST[:PORT]``."""
    host, separator, port = text.rpartition(":")
    if not separator:
        return text, DEFAULT_PORT
    if not port.isdigit() or not host:
        raise ValueError(f"invalid server {text!r}: use HOST[:PORT]")
    return host, int(port)


def _client(args, host: str, port: int) -> Client:
    """Client configured by the common options."""
    return Client(
        host,
        port=port,
        protocol=args.protocol,
        timeout=args.timeout,
        state_file=getattr(args, "state", None),
    )


def _select(client: Client, args) -> None:
    """Select the channels of the command line."""
    if not args.patterns and args.file is None:
        raise ValueError("give channel patterns or --file")
    client.select(args.patterns or None, file=args.file)


def _rotate(text: Optional[str]):
    """Rotation of the archive: seconds or a string with a unit."""
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return text


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------


def _print_record(client: Client, record) -> None:
    """One line per record."""
    print(
        f"{record.sid:<16} {record.time}  {record.nsamp:>5d} samples  "
        f"{record.rate:g} Hz",
        flush=True,
    )


def _print_status(client: Client) -> None:
    """Table of the channel latencies."""
    client.info()
    print(flush=True)


def _time(value) -> str:
    return str(value) if value is not None else "-"


def _count(items: list, name: str) -> None:
    print(f"{len(items)} {name}(s)", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(prog="python -m shakelab.signals.seedlink.cli"))
