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
SeedLink protocol: connection, handshake, commands and packets.

This module speaks the SeedLink protocol over a TCP connection, in the
two versions in use:

- version 4, specified by the FDSN (``docs.fdsn.org/projects/
  seedlink``): packets with an ``SE`` header (format, subformat,
  payload length, 64-bit sequence number, station identifier), MiniSEED
  2 or 3 payloads, INFO responses in JSON, time windows in ``DATA`` and
  dial-up mode with ``ENDFETCH``;
- version 3.1, the de facto standard of SeisComP and RingServer, which
  has no formal specification: 8-byte ``SL`` headers with a 24-bit
  hexadecimal sequence number followed by a 512-byte MiniSEED 2 record,
  INFO responses as XML in MiniSEED log records, ``TIME`` and ``FETCH``
  commands.

The version is negotiated from the ``HELLO`` response. Earlier
versions and the single-station mode are not supported.

The module is the low layer of :mod:`shakelab.signals.seedlink`: it
knows commands and packets, not waveforms. Packets are decoded and
distributed by :class:`~shakelab.signals.seedlink.client.Client`.
"""

from __future__ import annotations

import datetime as _dt
import re
import socket
import struct
from dataclasses import dataclass
from importlib import metadata as _importlib_metadata
from typing import Any, Optional

__all__ = [
    "DEFAULT_PORT",
    "PROTOCOLS",
    "Connection",
    "Packet",
    "SeedLinkError",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

DEFAULT_PORT = 18000

# Protocol versions spoken by the client.
PROTOCOLS = ("3.1", "4.0")

# Version 3: "SL" and six hexadecimal digits, or "SLINFO" and a flag,
# followed by a 512-byte MiniSEED 2 record.
V3_HEADER_SIZE = 8
V3_RECORD_LENGTH = 512
V3_SEQUENCE_MODULUS = 2 ** 24

# Version 4: signature, format, subformat, payload length, sequence
# number and station identifier length (little-endian), then the
# station identifier and the payload.
V4_HEADER = struct.Struct("<2sccIQB")

# Sent by the server, without line terminator, at the end of the data
# of a dial-up session.
END_SIGNAL = b"END"

# Longest command line accepted by servers, terminator included.
MAX_COMMAND_LENGTH = 255

# Longest reply line accepted from servers.
MAX_LINE_LENGTH = 4096

_RECV_SIZE = 65536
_HELLO_VERSION = re.compile(r"SeedLink\s+v(\d+)\.(\d+)", re.IGNORECASE)


# ----------------------------------------------------------------------------
# Errors and packets
# ----------------------------------------------------------------------------


class SeedLinkError(RuntimeError):
    """
    Error reported by a SeedLink server, or protocol violation.

    Parameters
    ----------
    message : str
        Description of the error.
    code : str, optional
        Error code of the server reply (version 4, or version 3 with
        extended replies), such as ``"ARGUMENTS"`` or ``"UNSUPPORTED"``.
    """

    def __init__(self, message: str, code: Optional[str] = None) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class Packet:
    """
    SeedLink packet.

    Attributes
    ----------
    kind : str
        ``"data"`` (MiniSEED payload of a station), ``"info"`` (response
        to an INFO request) or ``"other"`` (any other version 4 payload).
    sequence : int or None
        Sequence number of the packet in the station (None for INFO).
    station : str or None
        Station identifier ``NET_STA`` declared in the header (version
        4 only; in version 3 it is read from the MiniSEED record).
    format : str
        Payload format code: ``"2"`` and ``"3"`` for MiniSEED 2 and 3,
        ``"J"`` for JSON, ``"X"`` for XML.
    subformat : str
        Payload subformat code: ``"D"`` for data, ``"I"`` for INFO,
        ``"E"`` for an INFO error (version 4).
    payload : bytes
        Packet content.
    last : bool
        False for version 3 INFO packets followed by others of the same
        response.
    """

    kind: str
    sequence: Optional[int]
    station: Optional[str]
    format: str
    subformat: str
    payload: bytes
    last: bool = True


# ----------------------------------------------------------------------------
# Connection
# ----------------------------------------------------------------------------


class Connection:
    """
    TCP connection to a SeedLink server.

    Parameters
    ----------
    host : str
        Server host name or address.
    port : int, optional
        Server port. Default is :data:`DEFAULT_PORT`.
    timeout : float, optional
        Timeout in seconds for connecting and for the replies to the
        commands. Default is 30 s.

    Attributes
    ----------
    version : str or None
        Negotiated protocol version (``"3.1"`` or ``"4.0"``), None
        before :meth:`open`.
    capabilities : frozenset of str
        Capability flags announced in the ``HELLO`` response.
    server : str
        First line of the ``HELLO`` response (software and version).
    organization : str
        Second line of the ``HELLO`` response (station or data center).

    Notes
    -----
    The connection is opened with :meth:`open`, which negotiates the
    protocol version, configured with :meth:`command` and switched to
    data transfer by sending ``END`` or ``ENDFETCH`` with :meth:`send`.
    Packets are then read with :meth:`read_packet`. Partial packets
    survive a read timeout, so that reads can be polled with a short
    timeout (:meth:`set_poll`).
    """

    def __init__(
        self,
        host: str,
        port: int = DEFAULT_PORT,
        timeout: float = 30.0,
    ) -> None:
        if not isinstance(host, str):
            raise TypeError("host must be a string")
        if isinstance(port, bool) or not isinstance(port, int):
            raise TypeError("port must be an integer")

        self.host = host
        self.port = port
        self.timeout = float(timeout)
        self.version: Optional[str] = None
        self.capabilities: frozenset = frozenset()
        self.server = ""
        self.organization = ""
        self._socket: Optional[socket.socket] = None
        self._buffer = bytearray()

    def __repr__(self) -> str:
        return (
            f"Connection({self.host}:{self.port}, "
            f"version={self.version})"
        )

    # -- handshake -------------------------------------------------------

    def open(self, protocol: str = "auto") -> None:
        """
        Connect and negotiate the protocol version.

        Parameters
        ----------
        protocol : {"auto", "3.1", "4.0"}, optional
            Version to use. ``"auto"`` (default) selects version 4 when
            the server supports it, version 3.1 otherwise.

        Raises
        ------
        OSError
            If the server cannot be reached.
        SeedLinkError
            If the server does not support the requested version.
        """
        requested = normalize_protocol(protocol)

        try:
            self._socket = socket.create_connection(
                (self.host, self.port), timeout=self.timeout
            )
        except OSError as error:
            message = (
                f"cannot connect to the SeedLink server "
                f"{self.host}:{self.port}: {error.strerror or error}"
            )
            if error.errno is None:
                raise type(error)(message) from error
            raise type(error)(error.errno, message) from error
        self._socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
        self._buffer.clear()

        self.send("HELLO")
        self.server = self._read_line()
        self.organization = self._read_line()

        offered, self.capabilities = parse_hello(self.server)
        self.version = choose_version(offered, requested)

        if self.version == "4.0":
            self.command("SLPROTO 4.0")
            self.command(f"USERAGENT shakelab/{_shakelab_version()}")

    def close(self) -> None:
        """
        Send ``BYE`` and close the connection (if open).
        """
        if self._socket is None:
            return
        try:
            self._socket.sendall(b"BYE\r\n")
        except OSError:
            pass
        try:
            self._socket.close()
        finally:
            self._socket = None

    def send(self, line: str) -> None:
        """
        Send a command line without waiting for a reply.

        Parameters
        ----------
        line : str
            Command, without terminator.

        Raises
        ------
        ValueError
            If the command is too long or not ASCII.
        """
        data = line.encode("ascii") + b"\r\n"
        if len(data) > MAX_COMMAND_LENGTH:
            raise ValueError(f"SeedLink command too long: {line!r}")
        self._require_socket().sendall(data)

    def command(self, line: str) -> None:
        """
        Send a command and check the ``OK`` reply.

        Parameters
        ----------
        line : str
            Command, without terminator.

        Raises
        ------
        SeedLinkError
            If the server replies with an error, or with an unexpected
            line.
        """
        self.send(line)
        reply = self._read_line()

        if reply == "OK":
            return
        if reply.startswith("ERROR"):
            fields = reply.split(None, 2)
            code = fields[1] if len(fields) > 1 else None
            text = fields[2] if len(fields) > 2 else ""
            message = f"{line!r} rejected by the server"
            if code:
                message += f": {code}"
            if text:
                message += f" {text}"
            raise SeedLinkError(message, code=code)

        raise SeedLinkError(f"unexpected reply to {line!r}: {reply!r}")

    def info(self, item: str) -> str:
        """
        Request INFO during the handshake.

        Parameters
        ----------
        item : str
            INFO item, such as ``"ID"``, ``"STATIONS"`` or
            ``"STREAMS"``.

        Returns
        -------
        str
            Response: XML (version 3) or JSON (version 4) text.

        Raises
        ------
        SeedLinkError
            If the server reports an error.
        """
        self.send(f"INFO {item}")
        parts = []
        while True:
            packet = self.read_packet()
            if packet is None:
                raise SeedLinkError("END received instead of INFO")
            if packet.kind != "info":
                continue
            parts.append(_info_text(packet))
            if packet.subformat == "E":
                raise SeedLinkError(f"INFO {item} failed: {parts[-1]}")
            if packet.last:
                return "".join(parts)

    def keepalive(self) -> None:
        """
        Send an ``INFO ID`` request during data transfer.

        The response arrives as an INFO packet among the data packets.
        """
        self.send("INFO ID")

    # -- data transfer ---------------------------------------------------

    def set_poll(self, timeout: Optional[float]) -> None:
        """
        Set the timeout of the reads in :meth:`read_packet`.

        Parameters
        ----------
        timeout : float or None
            Seconds to wait for data before ``TimeoutError``; None
            blocks.
        """
        self._require_socket().settimeout(timeout)

    def read_packet(self) -> Optional[Packet]:
        """
        Read the next packet.

        Returns
        -------
        Packet or None
            Packet, or None when the server signals the end of a
            dial-up session.

        Raises
        ------
        TimeoutError
            If no complete packet arrives within the poll timeout (the
            bytes received so far are kept for the next call).
        ConnectionError
            If the server closes the connection.
        SeedLinkError
            If the bytes received are not a valid packet.
        """
        self._fill(len(END_SIGNAL))
        if self._buffer[:len(END_SIGNAL)] == END_SIGNAL:
            del self._buffer[:len(END_SIGNAL)]
            return None

        if self.version == "4.0":
            return self._read_v4()
        return self._read_v3()

    def _read_v3(self) -> Packet:
        """Read a version 3 packet."""
        self._fill(V3_HEADER_SIZE)
        header = bytes(self._buffer[:V3_HEADER_SIZE])

        if header.startswith(b"SLINFO"):
            payload = self._take(V3_HEADER_SIZE, V3_RECORD_LENGTH)
            return Packet(
                kind="info",
                sequence=None,
                station=None,
                format="2",
                subformat="I",
                payload=payload,
                last=header[7:8] != b"*",
            )

        if not header.startswith(b"SL"):
            raise SeedLinkError(f"invalid SeedLink packet header {header!r}")
        try:
            sequence = int(header[2:], 16)
        except ValueError:
            raise SeedLinkError(
                f"invalid SeedLink sequence number {header!r}"
            ) from None

        payload = self._take(V3_HEADER_SIZE, V3_RECORD_LENGTH)
        return Packet(
            kind="data",
            sequence=sequence,
            station=None,
            format="2",
            subformat="D",
            payload=payload,
        )

    def _read_v4(self) -> Packet:
        """Read a version 4 packet."""
        self._fill(V4_HEADER.size)
        signature, fmt, sub, length, sequence, id_length = (
            V4_HEADER.unpack_from(self._buffer)
        )
        if signature != b"SE":
            raise SeedLinkError(
                "invalid SeedLink packet signature "
                f"{bytes(self._buffer[:2])!r}"
            )

        header_size = V4_HEADER.size + id_length
        self._fill(header_size)
        station = bytes(self._buffer[V4_HEADER.size:header_size]).decode(
            "ascii", errors="replace"
        )
        payload = self._take(header_size, length)

        fmt = fmt.decode("ascii", errors="replace")
        sub = sub.decode("ascii", errors="replace")
        if fmt in ("J", "X"):
            kind = "info"
        elif sub == "D" and fmt in ("2", "3"):
            kind = "data"
        else:
            kind = "other"

        return Packet(
            kind=kind,
            sequence=sequence if kind != "info" else None,
            station=station or None,
            format=fmt,
            subformat=sub,
            payload=payload,
        )

    # -- low-level reading -----------------------------------------------

    def _require_socket(self) -> socket.socket:
        if self._socket is None:
            raise ConnectionError("the SeedLink connection is not open")
        return self._socket

    def _fill(self, size: int) -> None:
        """Receive until the buffer holds at least ``size`` bytes."""
        sock = self._require_socket()
        while len(self._buffer) < size:
            chunk = sock.recv(_RECV_SIZE)
            if not chunk:
                raise ConnectionError(
                    "connection closed by the SeedLink server"
                )
            self._buffer += chunk

    def _take(self, header_size: int, payload_size: int) -> bytes:
        """Remove a packet from the buffer and return its payload."""
        total = header_size + payload_size
        self._fill(total)
        payload = bytes(self._buffer[header_size:total])
        del self._buffer[:total]
        return payload

    def _read_line(self) -> str:
        """Read a reply line (terminated by LF, optionally after CR)."""
        sock = self._require_socket()
        sock.settimeout(self.timeout)
        while True:
            end = self._buffer.find(b"\n")
            if end >= 0:
                line = bytes(self._buffer[:end])
                del self._buffer[:end + 1]
                return line.rstrip(b"\r").decode("ascii", errors="replace")
            if len(self._buffer) > MAX_LINE_LENGTH:
                raise SeedLinkError("SeedLink reply line too long")
            chunk = sock.recv(_RECV_SIZE)
            if not chunk:
                raise ConnectionError(
                    "connection closed by the SeedLink server"
                )
            self._buffer += chunk


# ----------------------------------------------------------------------------
# Version negotiation
# ----------------------------------------------------------------------------


def normalize_protocol(protocol: Any) -> str:
    """
    Return the canonical name of a protocol selection.

    Parameters
    ----------
    protocol : str or int
        ``"auto"``, a version (``"3.1"``, ``"4.0"``) or its major
        number (``"3"``, ``"4"``, ``3``, ``4``).

    Returns
    -------
    str
        ``"auto"``, ``"3.1"`` or ``"4.0"``.
    """
    if isinstance(protocol, bool) or not isinstance(protocol, (str, int)):
        raise TypeError("protocol must be a string")
    text = str(protocol).strip().lower()
    if text == "auto":
        return "auto"
    for version in PROTOCOLS:
        if text in (version, version.split(".")[0]):
            return version
    raise ValueError(
        f"unsupported SeedLink protocol {protocol!r}: use 'auto', "
        + ", ".join(repr(version) for version in PROTOCOLS)
    )


def parse_hello(line: str) -> tuple:
    """
    Parse the first line of a ``HELLO`` response.

    Parameters
    ----------
    line : str
        For example ``"SeedLink v3.1 (2017.052 RingServer) ::
        SLPROTO:3.1 CAP EXTREPLY NSWILDCARD BATCH WS:13"``.

    Returns
    -------
    versions : list of str
        Protocol versions offered (``SLPROTO`` entries, or the version
        in the banner if there are none).
    capabilities : frozenset of str
        Other flags after ``::``.
    """
    banner, _, flags = line.partition("::")
    tokens = flags.split()

    versions = [
        token.split(":", 1)[1]
        for token in tokens
        if token.upper().startswith("SLPROTO:")
    ]
    capabilities = frozenset(
        token for token in tokens
        if not token.upper().startswith("SLPROTO:")
    )

    if not versions:
        match = _HELLO_VERSION.search(banner)
        if match:
            versions = [f"{match.group(1)}.{match.group(2)}"]

    return versions, capabilities


def choose_version(offered: list, requested: str = "auto") -> str:
    """
    Choose the protocol version for a session.

    Parameters
    ----------
    offered : list of str
        Versions offered by the server (see :func:`parse_hello`).
    requested : {"auto", "3.1", "4.0"}
        Version requested by the user.

    Returns
    -------
    str
        ``"3.1"`` or ``"4.0"``.

    Raises
    ------
    SeedLinkError
        If no supported version is offered, or the requested one is not.
    """
    majors = set()
    for version in offered:
        major = version.split(".")[0]
        if major.isdigit():
            majors.add(int(major))

    if requested == "auto":
        if 4 in majors:
            return "4.0"
        if 3 in majors:
            return "3.1"
    elif int(requested.split(".")[0]) in majors:
        return requested

    wanted = "version 3 or 4" if requested == "auto" else requested
    raise SeedLinkError(
        f"the server offers SeedLink {', '.join(offered) or 'unknown'}, "
        f"not {wanted}",
        code="UNSUPPORTED",
    )


# ----------------------------------------------------------------------------
# Commands
# ----------------------------------------------------------------------------


def station_command(version: str, network: str, station: str) -> str:
    """
    Return the ``STATION`` command of a station or station pattern.
    """
    if version == "4.0":
        return f"STATION {network}_{station}"
    return f"STATION {station} {network}"


def select_commands(version: str, selectors: list) -> list:
    """
    Return the ``SELECT`` commands of a station.

    Parameters
    ----------
    version : str
        Protocol version.
    selectors : list of tuple
        ``(location, channel)`` patterns with the ``*`` and ``?``
        jokers; an empty location is the blank location code.

    Returns
    -------
    list of str
        Commands; empty if every stream is requested. Patterns that the
        protocol cannot express exactly are widened: the client filters
        the records it receives in any case.
    """
    patterns = []
    for location, channel in selectors:
        pattern = (
            _v4_selector(location, channel)
            if version == "4.0"
            else _v3_selector(location, channel)
        )
        if pattern is None:
            return []
        if pattern not in patterns:
            patterns.append(pattern)
    return [f"SELECT {pattern}" for pattern in patterns]


def _v3_selector(location: str, channel: str) -> Optional[str]:
    """
    Version 3 selector ``LLCCC.T``, or None for every stream.

    An omitted location matches any location (the blank one cannot be
    written), and ``*`` cannot be written in a channel.
    """
    explicit_location = (
        len(location) == 2 and "*" not in location
    )
    explicit_channel = len(channel) == 3 and "*" not in channel

    if explicit_channel:
        pattern = channel
    elif explicit_location:
        pattern = "???"
    else:
        return None

    if explicit_location:
        pattern = location + pattern
    return pattern + ".D"


def _v4_selector(location: str, channel: str) -> Optional[str]:
    """
    Version 4 stream pattern ``LOC_B_S_SS``, or None for every stream.
    """
    if len(channel) == 3 and "*" not in channel:
        stream = "_".join(channel)
    else:
        stream = "*"
    if location == "*" and stream == "*":
        return None
    return f"{location}_{stream}"


def data_command(
    version: str,
    sequence: Optional[int] = None,
    start_ns: Optional[int] = None,
    end_ns: Optional[int] = None,
    fetch: bool = False,
) -> str:
    """
    Return the command that starts the data of the current station.

    Parameters
    ----------
    version : str
        Protocol version.
    sequence : int, optional
        Sequence number of the first packet wanted.
    start_ns, end_ns : int, optional
        Time window, in nanoseconds from the Unix epoch.
    fetch : bool, optional
        Dial-up mode (version 3 ``FETCH``; in version 4 the mode is set
        by ``ENDFETCH``).

    Notes
    -----
    Version 3 cannot combine a sequence number and an end time: the
    end time is dropped and the client stops the transfer.
    """
    if version == "4.0":
        parts = ["DATA"]
        if sequence is not None:
            parts.append(str(sequence))
        elif start_ns is not None:
            parts.append("ALL")
        if start_ns is not None:
            parts.append(_v4_time(start_ns))
            if end_ns is not None:
                parts.append(_v4_time(end_ns))
        return " ".join(parts)

    if sequence is None and start_ns is not None:
        command = f"TIME {_v3_time(start_ns)}"
        if end_ns is not None:
            command += f" {_v3_time(end_ns, round_up=True)}"
        return command

    command = "FETCH" if fetch else "DATA"
    if sequence is not None:
        command += f" {sequence % V3_SEQUENCE_MODULUS:06X}"
        if start_ns is not None:
            command += f" {_v3_time(start_ns)}"
    return command


def start_command(version: str, fetch: bool = False) -> str:
    """
    Return the command that ends the handshake.
    """
    if fetch and version == "4.0":
        return "ENDFETCH"
    return "END"


def next_sequence(version: str, sequence: int) -> int:
    """
    Return the sequence number following ``sequence``.
    """
    if version == "4.0":
        return sequence + 1
    return (sequence + 1) % V3_SEQUENCE_MODULUS


def is_newer(version: str, sequence: int, last: int) -> bool:
    """
    Return whether ``sequence`` follows ``last`` in a station.

    Version 3 numbers wrap around after 24 bits: a number is newer if
    it is less than half the range ahead.
    """
    if version == "4.0":
        return sequence > last
    ahead = (sequence - last) % V3_SEQUENCE_MODULUS
    return 0 < ahead < V3_SEQUENCE_MODULUS // 2


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _utc(ns: int) -> _dt.datetime:
    """UTC datetime (microseconds) of a Unix time in nanoseconds."""
    epoch = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
    return epoch + _dt.timedelta(microseconds=ns // 1000)


def _v3_time(ns: int, round_up: bool = False) -> str:
    """Version 3 time ``YYYY,MM,DD,hh,mm,ss`` (whole seconds)."""
    seconds, rest = divmod(ns, 10 ** 9)
    if round_up and rest:
        seconds += 1
    return _utc(seconds * 10 ** 9).strftime("%Y,%m,%d,%H,%M,%S")


def _v4_time(ns: int) -> str:
    """Version 4 time ``YYYY-MM-DDThh:mm:ss.ffffffZ``."""
    return _utc(ns).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _info_text(packet: Packet) -> str:
    """Text of an INFO packet."""
    if packet.format == "2":
        # Version 3: XML in the ASCII data of a MiniSEED log record.
        from shakelab.signals.libio.mseed import MSRecord

        record = MSRecord.from_bytes(packet.payload)
        return record.data if isinstance(record.data, str) else ""
    return packet.payload.decode("utf-8", errors="replace")


def _shakelab_version() -> str:
    """Installed ShakeLab version, for the ``USERAGENT`` command."""
    try:
        return _importlib_metadata.version("shakelab")
    except _importlib_metadata.PackageNotFoundError:
        return "unknown"
