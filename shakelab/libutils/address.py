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
Server addresses of the network clients.

Every client of ShakeLab that connects to a server by host name and
port reads the address with :func:`parse_address`, so that the same
text means the same server everywhere (in Python and on the command
line):

- ``host``, ``host:port``, an IPv4 address with an optional port, an
  IPv6 address alone (``::1``) or in brackets with a port
  (``[::1]:18000``);
- for clients of web services, also ``scheme://host[:port][/path]``;
  without a scheme the client's default applies (``https`` for the
  FDSN web services);
- a ``port`` given as a separate argument always takes precedence over
  the one written in the address.

Examples
--------
>>> from shakelab.libutils.address import parse_address
>>> parse_address("rtserve.example.org:18001", default_port=18000)
Address(host='rtserve.example.org', port=18001, scheme=None, path='')
>>> str(parse_address("[::1]", default_port=18000))
'[::1]:18000'
>>> address = parse_address("example.org:8080/fdsnws", port=8081,
...                         schemes=("http", "https"),
...                         default_scheme="https")
>>> str(address)
'https://example.org:8081/fdsnws'
"""

from __future__ import annotations

import ipaddress
import re
from typing import Any, NamedTuple, Optional

__all__ = [
    "Address",
    "parse_address",
]

# Letters, digits, hyphens, underscores and dots: host names and IPv4
# addresses (a stricter check is left to the name resolution).
_HOST_NAME = re.compile(r"^[A-Za-z0-9_](?:[A-Za-z0-9_.-]*[A-Za-z0-9_])?$")


# ---------------------------------------------------------------------------
class Address(NamedTuple):
    """
    Address of a server.

    Attributes
    ----------
    host : str
        Host name or IP address (IPv6 without brackets).
    port : int or None
        Port, or None for the default port of the scheme.
    scheme : str or None
        Scheme in lower case (``"https"``), or None for protocols that
        have none (SeedLink).
    path : str
        Path prefix starting with ``/``, or ``""``.
    """

    host: str
    port: Optional[int]
    scheme: Optional[str] = None
    path: str = ""

    @property
    def netloc(self) -> str:
        """Host and port as written in an address (``host:port``)."""
        host = f"[{self.host}]" if ":" in self.host else self.host
        return host if self.port is None else f"{host}:{self.port}"

    def __str__(self) -> str:
        if self.scheme is None:
            return self.netloc + self.path
        return f"{self.scheme}://{self.netloc}{self.path}"


def parse_address(
    text: str,
    *,
    port: Optional[int] = None,
    default_port: Optional[int] = None,
    schemes: tuple = (),
    default_scheme: Optional[str] = None,
) -> Address:
    """
    Read the address of a server.

    Parameters
    ----------
    text : str
        ``host``, ``host:port``, ``[ipv6]:port`` or, if ``schemes`` is
        not empty, ``scheme://host[:port][/path]`` (the scheme and the
        path are optional).
    port : int, optional
        Port that replaces the one in ``text``.
    default_port : int, optional
        Port when neither ``text`` nor ``port`` gives one. Default is
        None (the default port of the scheme).
    schemes : tuple of str, optional
        Accepted schemes, in lower case. Empty (default) for protocols
        without a scheme: ``text`` must not have one, nor a path.
    default_scheme : str, optional
        Scheme of an address written without one. Required when
        ``schemes`` is not empty.

    Returns
    -------
    Address
        The address.

    Raises
    ------
    TypeError
        If ``text`` is not a string or ``port`` not an integer.
    ValueError
        If the address is empty or malformed, has a scheme that is not
        accepted, or a port outside 1-65535.
    """
    if not isinstance(text, str):
        raise TypeError("server addresses must be strings")
    if schemes and default_scheme not in schemes:
        raise ValueError("default_scheme must be one of schemes")
    if port is not None:
        port = _check_port(port, "port")
    address = text.strip()
    if not address:
        raise ValueError("empty server address")
    if any(char.isspace() for char in address):
        raise ValueError(f"blank in the server address {text!r}")

    scheme = default_scheme if schemes else None
    name, separator, rest = address.partition("://")
    if separator:
        if not schemes:
            raise ValueError(
                f"{text!r}: write the server as HOST[:PORT], without "
                f"a scheme"
            )
        scheme = name.lower()
        if scheme not in schemes:
            raise ValueError(
                f"{text!r}: the scheme must be "
                f"{' or '.join(s + '://' for s in schemes)}"
            )
        address = rest

    if any(char in address for char in "?#@"):
        raise ValueError(
            f"{text!r}: queries, fragments and user names are not "
            f"accepted in server addresses"
        )
    netloc, slash, path = address.partition("/")
    path = slash + path
    if path and not schemes:
        raise ValueError(f"{text!r}: no path is accepted after the port")
    path = path.rstrip("/")

    host, written = _split_netloc(netloc, text)
    if port is None and written is not None:
        if not written.isdigit() or not 1 <= int(written) <= 65535:
            raise ValueError(f"invalid port in {text!r}")
        port = int(written)
    if port is None and default_port is not None:
        port = _check_port(default_port, "default_port")
    return Address(host=host, port=port, scheme=scheme, path=path)


# ---------------------------------------------------------------------------

def _split_netloc(netloc: str, text: str) -> tuple:
    """Host and port (as text, or None) of ``host[:port]``."""
    if netloc.startswith("["):
        host, bracket, rest = netloc[1:].partition("]")
        if not bracket or (rest and not rest.startswith(":")):
            raise ValueError(f"invalid IPv6 address in {text!r}")
        _check_ipv6(host, text)
        return host, (rest[1:] if rest else None)
    if netloc.count(":") > 1:
        # An IPv6 address without brackets cannot carry a port.
        _check_ipv6(netloc, text)
        return netloc, None
    host, colon, port = netloc.partition(":")
    if not host:
        raise ValueError(f"server address without host: {text!r}")
    if not _HOST_NAME.match(host):
        raise ValueError(f"invalid host name in {text!r}")
    return host, (port if colon else None)


def _check_ipv6(host: str, text: str) -> None:
    try:
        ipaddress.IPv6Address(host)
    except ValueError:
        raise ValueError(f"invalid IPv6 address in {text!r}") from None


def _check_port(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if not 1 <= value <= 65535:
        raise ValueError(f"{name} must be in 1-65535: {value}")
    return value
