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
Data centers offering FDSN web services.

The built-in table follows the FDSN registry of data centers
(https://www.fdsn.org/datacenters/, read on 9 November 2026), with the
address of each service: the services of a data center may live on
different hosts or follow non-standard paths. A few services missing
from the registry are added (the event service of INGV, the servers of
OGS). :func:`datacenters` can refresh the table from the registry.

Names are those of the registry; they are matched without regard to
case, and some former names are accepted as aliases (``IRIS`` for
EarthScope, ``ORFEUS`` for ODC, ``RESIF`` for EPOSFR).
"""

from __future__ import annotations

import json
from typing import Any, NamedTuple, Optional

__all__ = [
    "DataCenter",
    "REGISTRY_URL",
    "SERVICES",
    "datacenters",
    "find_datacenter",
    "parse_registry",
    "service_root",
    "service_urls",
]

REGISTRY_URL = "https://www.fdsn.org/ws/datacenters/1/query"

# Services used by the client, with their names in the registry.
SERVICES = ("dataselect", "station", "event", "availability")


class DataCenter(NamedTuple):
    """
    Data center and the addresses of its FDSN web services.

    Attributes
    ----------
    name : str
        Name, as in the FDSN registry.
    services : dict
        Address of each service (``"dataselect"``, ``"station"``,
        ``"event"``, ``"availability"``), ending with ``/``; the
        methods (``query``, ``version``...) are appended to it.
    website : str or None
        Web site of the data center.
    """

    name: str
    services: dict
    website: Optional[str] = None


def _standard(host: str, services: tuple) -> dict:
    """Standard addresses ``<host>/fdsnws/<service>/1/``."""
    return {name: f"{host}/fdsnws/{name}/1/" for name in services}


_DSA = ("dataselect", "station", "availability")
_DS = ("dataselect", "station")

_TABLE = [
    DataCenter(
        "AusPass",
        _standard("https://auspass.edu.au", _DSA + ("event",)),
        "https://auspass.edu.au",
    ),
    DataCenter(
        "BGR", _standard("https://eida.bgr.de", _DSA),
        "https://www.bgr.bund.de",
    ),
    DataCenter(
        "BGS", _standard("https://eida.bgs.ac.uk", _DSA),
        "https://eida.bgs.ac.uk",
    ),
    DataCenter(
        "EarthScope", _standard("https://service.earthscope.org", _DSA),
        "https://www.earthscope.org/data/",
    ),
    DataCenter(
        "EMSC", _standard("https://www.seismicportal.eu", ("event",)),
        "https://www.emsc-csem.org/",
    ),
    DataCenter(
        "EPOSFR",
        {
            **_standard("https://seisdata.epos-france.fr", _DSA),
            "event": "https://api.franceseisme.fr/fdsnws/event/1/",
        },
        "https://seismology.epos-france.fr",
    ),
    DataCenter(
        "ETH", _standard("https://eida.ethz.ch", _DSA + ("event",)),
        "https://www.ethz.ch",
    ),
    DataCenter(
        "GEOFON", _standard("https://geofon.gfz.de", _DSA + ("event",)),
        "https://geofon.gfz.de/",
    ),
    DataCenter(
        "GeoNet",
        _standard("https://service.geonet.org.nz", _DS + ("event",)),
        "https://www.geonet.org.nz/",
    ),
    DataCenter(
        "GEOSHAKE",
        _standard("https://api.geoshake.org", ("station", "event")),
        "https://geoshake.org/open-data",
    ),
    DataCenter(
        "ICGC", _standard("https://ws.icgc.cat", _DSA),
        "https://www.icgc.cat/en/terratremols",
    ),
    DataCenter(
        "IESDMC", _standard("http://batsws.earth.sinica.edu.tw", _DS),
        "http://batsws.earth.sinica.edu.tw/fdsnws/",
    ),
    DataCenter(
        "IGN",
        {
            "dataselect": (
                "http://fdsnws.sismologia.ign.es/fdsnws/dataselect/"
            ),
            "station": "http://fdsnws.sismologia.ign.es/fdsnws/station/",
        },
        "http://www.ign.es",
    ),
    DataCenter(
        "INGV",
        _standard("https://webservices.ingv.it", _DSA + ("event",)),
        "https://www.eida.ingv.it",
    ),
    DataCenter(
        "IPGP", _standard("http://ws.ipgp.fr", _DS),
        "http://datacenter.ipgp.fr",
    ),
    DataCenter(
        "ISC", _standard("http://www.isc.ac.uk", ("event",)),
        "http://www.isc.ac.uk",
    ),
    DataCenter(
        "KAGSR", _standard("http://sdis.emsd.ru", ("station", "event")),
        "http://sdis.emsd.ru",
    ),
    DataCenter(
        "KIT",
        {
            "dataselect": "http://ws.gpi.kit.edu/fdsnws/dataselect/",
            "station": "http://ws.gpi.kit.edu/fdsnws/station/1/",
        },
        "http://ws.gpi.kit.edu",
    ),
    DataCenter(
        "KOERI", _standard("https://eida.koeri.boun.edu.tr", _DSA),
        "https://www.koeri.boun.edu.tr",
    ),
    DataCenter(
        "LMU",
        _standard(
            "https://erde.geophysik.uni-muenchen.de", _DSA + ("event",)
        ),
        "https://www.geophysik.uni-muenchen.de",
    ),
    DataCenter(
        "NCEDC", _standard("https://service.ncedc.org", _DSA + ("event",)),
        "https://ncedc.org",
    ),
    DataCenter(
        "NIEP", _standard("https://eida-sc3.infp.ro", _DSA),
        "https://www.infp.ro",
    ),
    DataCenter(
        "NOA", _standard("https://eida.gein.noa.gr", _DSA + ("event",)),
        "https://www.noa.gr",
    ),
    DataCenter(
        "NRCAN", _standard("https://earthquakescanada.nrcan.gc.ca", _DS),
        "https://earthquakescanada.nrcan.gc.ca/stndon/CNDC/index-en.php",
    ),
    DataCenter(
        "ODC", _standard("https://www.orfeus-eu.org", _DSA),
        "https://www.orfeus-eu.org/",
    ),
    DataCenter(
        "OGS-ANT", _standard("http://158.110.30.202:5600", _DS), None,
    ),
    DataCenter(
        "OGS-SCP", _standard("http://158.110.30.85:8080", _DS), None,
    ),
    DataCenter(
        "OSDC", _standard("http://osdc.obslab.cn", _DS),
        "https://obslab.sustech.edu.cn/mda/",
    ),
    DataCenter(
        "RASPISHAKE", _standard("https://data.raspberryshake.org", _DS),
        "http://raspberryshake.org/",
    ),
    DataCenter(
        "RSBR", _standard("http://rsbr.iag.usp.br", _DS),
        "http://rsbr.on.br/",
    ),
    DataCenter(
        "SCEDC",
        _standard("https://service.scedc.caltech.edu", _DSA + ("event",)),
        "https://scedc.caltech.edu",
    ),
    DataCenter(
        "UIB-NORSAR", _standard("https://eida.geo.uib.no", _DSA),
        "https://eida.geo.uib.no",
    ),
    DataCenter(
        "USGS", _standard("https://earthquake.usgs.gov", ("event",)),
        "https://earthquake.usgs.gov/",
    ),
    DataCenter(
        "USP",
        {
            **_standard("http://seisrequest.iag.usp.br", _DS),
            "event": "http://www.moho.iag.usp.br/fdsnws/event/1/",
        },
        "http://www.moho.iag.usp.br/",
    ),
]

# Services added to those of the registry when it is refreshed.
_EXTRA = {
    "INGV": {"event": "https://webservices.ingv.it/fdsnws/event/1/"},
    "OGS-ANT": _standard("http://158.110.30.202:5600", _DS),
    "OGS-SCP": _standard("http://158.110.30.85:8080", _DS),
}

_ALIASES = {
    "IRIS": "EarthScope",
    "IRISDMC": "EarthScope",
    "ORFEUS": "ODC",
    "RESIF": "EPOSFR",
}

_CURRENT = {center.name: center for center in _TABLE}


# ----------------------------------------------------------------------------
# Registry
# ----------------------------------------------------------------------------


def datacenters(refresh: bool = False, timeout: float = 30.0) -> dict:
    """
    Return the known data centers.

    Parameters
    ----------
    refresh : bool, optional
        If True, read the FDSN registry (:data:`REGISTRY_URL`) and use
        it in place of the built-in table for the rest of the session,
        with the additions listed in the module notes. Default is
        False.
    timeout : float, optional
        Timeout of the registry request, in seconds.

    Returns
    -------
    dict
        :class:`DataCenter` by name, in alphabetical order (without
        regard to case).

    Raises
    ------
    FdsnError
        If the registry cannot be read.
    """
    global _CURRENT
    if refresh:
        from shakelab.signals.fdsnws.transport import FdsnError, request

        content = request(REGISTRY_URL, timeout=timeout, retries=1)
        if content is None:
            raise FdsnError("the FDSN registry returned no data")
        try:
            table = parse_registry(json.loads(content.decode("utf-8")))
        except (ValueError, KeyError, TypeError) as error:
            raise FdsnError(f"invalid FDSN registry: {error}") from None
        _CURRENT = table
    return dict(sorted(_CURRENT.items(), key=lambda item: item[0].lower()))


def find_datacenter(name: str) -> DataCenter:
    """
    Return a data center by name.

    Parameters
    ----------
    name : str
        Name of the registry or alias, without regard to case.

    Returns
    -------
    DataCenter
        The data center.

    Raises
    ------
    TypeError
        If ``name`` is not a string.
    ValueError
        If the name is unknown.
    """
    if not isinstance(name, str):
        raise TypeError("data center names must be strings")
    key = name.strip().upper()
    for alias, target in _ALIASES.items():
        if key == alias:
            key = target.upper()
    for center in _CURRENT.values():
        if center.name.upper() == key:
            return center
    raise ValueError(
        f"unknown data center {name!r}; known: "
        f"{', '.join(datacenters())}"
    )


def service_urls(base: str) -> dict:
    """
    Return the standard service addresses of a host.

    Parameters
    ----------
    base : str
        ``http://`` or ``https://`` address of the host, with an
        optional port and path prefix (e.g. ``https://example.org``).

    Returns
    -------
    dict
        ``<base>/fdsnws/<service>/1/`` for each of :data:`SERVICES`.

    Raises
    ------
    ValueError
        If the address does not start with ``http://`` or
        ``https://``.
    """
    return _standard(_check_url(base).rstrip("/"), SERVICES)


def service_root(url: str) -> str:
    """
    Return the address of a service, ready to append a method to.

    Parameters
    ----------
    url : str
        Address of a service (``.../fdsnws/station/1/``), with or
        without the final ``/`` or a method name (``query``,
        ``queryauth``).

    Returns
    -------
    str
        The address ending with ``/``.

    Raises
    ------
    TypeError
        If ``url`` is not a string.
    ValueError
        If it does not start with ``http://`` or ``https://``.
    """
    url = _check_url(url)
    for method in ("queryauth", "query"):
        if url.endswith("/" + method):
            url = url[: -len(method)]
    return url if url.endswith("/") else url + "/"


def parse_registry(data: dict) -> dict:
    """
    Return the data centers of the JSON content of the FDSN registry.

    Parameters
    ----------
    data : dict
        Decoded reply of :data:`REGISTRY_URL`.

    Returns
    -------
    dict
        :class:`DataCenter` by name. The first address of each service
        is kept (a data center may list one per repository); data
        centers without any of :data:`SERVICES` are left out.
    """
    table = {}
    for center in data["datacenters"]:
        name = center["name"]
        services: dict = {}
        for repository in center.get("repositories", []):
            for service in repository.get("services", []):
                kind = _service_kind(service.get("name", ""))
                url = service.get("url")
                if kind is None or not url or kind in services:
                    continue
                try:
                    services[kind] = service_root(url)
                except ValueError:
                    continue
        services.update(_EXTRA.get(name, {}))
        if services:
            table[name] = DataCenter(name, services, center.get("website"))
    for name, services in _EXTRA.items():
        if name not in table:
            table[name] = DataCenter(name, dict(services), None)
    return table


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _service_kind(name: str) -> Optional[str]:
    """Service of a registry name such as ``fdsnws-station-1``."""
    for kind in SERVICES:
        if name == f"fdsnws-{kind}-1":
            return kind
    return None


def _check_url(url: Any) -> str:
    if not isinstance(url, str):
        raise TypeError("service addresses must be strings")
    url = url.strip()
    if not url.lower().startswith(("http://", "https://")):
        raise ValueError(
            f"service addresses must start with http:// or https://: "
            f"{url!r}"
        )
    return url
