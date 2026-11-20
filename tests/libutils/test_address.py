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
Unit tests for :mod:`shakelab.libutils.address`.

The tests check the forms of address accepted by the network clients,
the precedence of a separate port, and the errors.
"""

from __future__ import annotations

import pytest

from shakelab.libutils.address import Address, parse_address

WEB = {"schemes": ("http", "https"), "default_scheme": "https"}


@pytest.mark.parametrize("text, host, port", [
    ("rtserve.example.org", "rtserve.example.org", 18000),
    ("rtserve.example.org:18001", "rtserve.example.org", 18001),
    (" localhost:5 ", "localhost", 5),
    ("192.168.1.10:8080", "192.168.1.10", 8080),
    ("::1", "::1", 18000),
    ("[::1]", "::1", 18000),
    ("[2001:db8::1]:18001", "2001:db8::1", 18001),
])
def test_host_and_port(text, host, port):
    address = parse_address(text, default_port=18000)
    assert address == Address(host, port)


def test_separate_port_takes_precedence():
    address = parse_address("example.org:18001", port=18002,
                            default_port=18000)
    assert address.port == 18002
    assert parse_address("example.org", port=18002).port == 18002
    assert parse_address("example.org").port is None


@pytest.mark.parametrize("text, expected", [
    ("example.org", "https://example.org"),
    ("example.org:8080", "https://example.org:8080"),
    ("example.org:8080/prefix/", "https://example.org:8080/prefix"),
    ("http://158.110.30.85:8080", "http://158.110.30.85:8080"),
    ("HTTPS://Example.org/fdsnws", "https://Example.org/fdsnws"),
    ("[::1]:8080", "https://[::1]:8080"),
])
def test_web_addresses(text, expected):
    assert str(parse_address(text, **WEB)) == expected


def test_netloc():
    assert parse_address("::1", default_port=1).netloc == "[::1]:1"
    assert parse_address("example.org").netloc == "example.org"


@pytest.mark.parametrize("text, message", [
    ("", "empty"),
    ("example.org:", "invalid port"),
    ("example.org:x", "invalid port"),
    ("example.org:0", "invalid port"),
    ("example.org:65536", "invalid port"),
    ("[::1", "IPv6"),
    ("[::1]x", "IPv6"),
    ("1::2::3", "IPv6"),
    (":18000", "without host"),
    ("-example.org", "host name"),
    ("example org", "blank"),
    ("user@example.org", "user names"),
    ("http://example.org", "without a scheme"),
    ("example.org/path", "no path"),
])
def test_invalid_addresses(text, message):
    with pytest.raises(ValueError, match=message):
        parse_address(text, default_port=18000)


def test_invalid_web_addresses():
    with pytest.raises(ValueError, match="http:// or https://"):
        parse_address("ftp://example.org", **WEB)
    with pytest.raises(ValueError, match="queries"):
        parse_address("example.org/x?y=1", **WEB)


def test_types():
    with pytest.raises(TypeError):
        parse_address(18000)
    with pytest.raises(TypeError):
        parse_address("example.org", port="18000")
    with pytest.raises(TypeError):
        parse_address("example.org", port=True)
    with pytest.raises(ValueError, match="1-65535"):
        parse_address("example.org", port=0)
    with pytest.raises(ValueError, match="default_scheme"):
        parse_address("example.org", schemes=("http",))
