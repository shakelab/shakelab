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
Tests of :mod:`shakelab.signals.fdsnws.registry`.

The built-in table, the names and aliases, the addresses of hosts and
services, and the refresh from the FDSN registry (served by the fake
server).
"""

from __future__ import annotations

import json

import pytest

from shakelab.signals.fdsnws import registry
from shakelab.signals.fdsnws.transport import FdsnError

REGISTRY = {
    "datacenters": [
        {
            "name": "ALPHA",
            "website": "https://alpha.example.org",
            "repositories": [
                {
                    "name": "archive",
                    "services": [
                        {"name": "fdsnws-station-1",
                         "url": "https://alpha.example.org/fdsnws/station/1"},
                        {"name": "fdsnws-dataselect-1",
                         "url": "https://alpha.example.org/fdsnws/"
                                "dataselect/1/query"},
                        {"name": "eidaws-wfcatalog",
                         "url": "https://alpha.example.org/eidaws/"},
                    ],
                },
                {
                    "name": "near-real-time",
                    "services": [
                        {"name": "fdsnws-dataselect-1",
                         "url": "https://nrt.example.org/fdsnws/"
                                "dataselect/1/"},
                    ],
                },
            ],
        },
        {"name": "EMPTY", "repositories": [{"name": "x", "services": []}]},
        {
            "name": "INGV",
            "repositories": [
                {"services": [
                    {"name": "fdsnws-station-1",
                     "url": "https://webservices.ingv.it/fdsnws/station/1/"},
                ]},
            ],
        },
    ]
}


@pytest.fixture
def restore_table(monkeypatch):
    """Keep the table of the session unchanged by the test."""
    monkeypatch.setattr(registry, "_CURRENT", dict(registry._CURRENT))


# ----------------------------------------------------------------------------
# Built-in table
# ----------------------------------------------------------------------------


def test_table_addresses():
    """Every service address is http(s) and ends with a slash."""
    centers = registry.datacenters()
    assert list(centers) == sorted(centers, key=str.lower)
    for center in centers.values():
        assert center.services
        assert set(center.services) <= set(registry.SERVICES)
        for url in center.services.values():
            assert url.startswith(("http://", "https://"))
            assert url.endswith("/")


@pytest.mark.parametrize(
    "name, expected",
    [
        ("INGV", "INGV"),
        ("ingv", "INGV"),
        (" Geofon ", "GEOFON"),
        ("IRIS", "EarthScope"),
        ("irisdmc", "EarthScope"),
        ("ORFEUS", "ODC"),
        ("resif", "EPOSFR"),
    ],
)
def test_find_datacenter(name, expected):
    """Names without regard to case, and former names."""
    assert registry.find_datacenter(name).name == expected


def test_find_datacenter_errors():
    """Unknown names and wrong types are rejected."""
    with pytest.raises(ValueError, match="unknown data center"):
        registry.find_datacenter("NOWHERE")
    with pytest.raises(TypeError):
        registry.find_datacenter(1)


def test_non_standard_services():
    """Services on other hosts or paths are kept as registered."""
    usp = registry.find_datacenter("USP").services
    assert usp["event"].startswith("http://www.moho.iag.usp.br/")
    assert usp["station"].startswith("http://seisrequest.iag.usp.br/")
    ign = registry.find_datacenter("IGN").services
    assert ign["station"].endswith("/fdsnws/station/")
    assert "dataselect" not in registry.find_datacenter("USGS").services


# ----------------------------------------------------------------------------
# Addresses
# ----------------------------------------------------------------------------


def test_service_urls():
    """Standard paths of a host."""
    urls = registry.service_urls("https://example.org:8080/")
    assert urls["station"] == "https://example.org:8080/fdsnws/station/1/"
    assert set(urls) == set(registry.SERVICES)


@pytest.mark.parametrize(
    "url, expected",
    [
        ("https://a.org/fdsnws/station/1", "https://a.org/fdsnws/station/1/"),
        ("https://a.org/fdsnws/station/1/query",
         "https://a.org/fdsnws/station/1/"),
        ("http://a.org/fdsnws/dataselect/1/queryauth",
         "http://a.org/fdsnws/dataselect/1/"),
    ],
)
def test_service_root(url, expected):
    """Final slash added, method names removed."""
    assert registry.service_root(url) == expected


def test_address_errors():
    """Only http and https addresses."""
    with pytest.raises(ValueError):
        registry.service_root("ftp://a.org/fdsnws/station/1/")
    with pytest.raises(TypeError):
        registry.service_urls(None)


# ----------------------------------------------------------------------------
# FDSN registry
# ----------------------------------------------------------------------------


def test_parse_registry():
    """First address of each service, additions, empty centers out."""
    table = registry.parse_registry(REGISTRY)
    alpha = table["ALPHA"]
    assert alpha.website == "https://alpha.example.org"
    assert alpha.services == {
        "station": "https://alpha.example.org/fdsnws/station/1/",
        "dataselect": "https://alpha.example.org/fdsnws/dataselect/1/",
    }
    assert "EMPTY" not in table
    assert table["INGV"].services["event"].startswith(
        "https://webservices.ingv.it/"
    )
    assert "OGS-SCP" in table


def test_refresh(server, monkeypatch, restore_table):
    """The table is replaced by the registry for the session."""
    server.replies["/registry"] = lambda req: (
        200, json.dumps(REGISTRY).encode(), {}
    )
    monkeypatch.setattr(registry, "REGISTRY_URL", server.url + "/registry")

    centers = registry.datacenters(refresh=True)
    assert "ALPHA" in centers
    assert "GEOFON" not in centers
    assert registry.find_datacenter("alpha").name == "ALPHA"


def test_refresh_errors(server, monkeypatch, restore_table):
    """An invalid registry raises FdsnError and keeps the table."""
    server.replies["/registry"] = lambda req: (200, b"[]", {})
    monkeypatch.setattr(registry, "REGISTRY_URL", server.url + "/registry")
    with pytest.raises(FdsnError, match="invalid FDSN registry"):
        registry.datacenters(refresh=True)
    assert "GEOFON" in registry.datacenters()
