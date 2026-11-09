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
Tests of :mod:`shakelab.signals.fdsnws.client`.

The client talks to the fake server of ``fdsn_server.py``, which
returns ``MODE.mseed`` (IV.MODE..HN?, 2012-05-20 02:03:23, 180 s at
100 Hz) and ``MODE.xml`` whatever the selection: the tests check the
requests sent, the filtering and cutting done by the client, the
metadata, the response correction (against the PGA of the ESM
processing of the same record), the raw downloads, the catalogue and
the availability.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals import io
from shakelab.signals.base import Record, StreamCollection
from shakelab.signals.fdsnws import Client, FdsnError, FdsnWarning, TimeSpan
from shakelab.signals.metadata import MetadataCollection

from fdsn_server import SIGNALS_DATA, read_file

T0 = "2012-05-20T02:04:00"
T1 = "2012-05-20T02:05:30"

# PGA of IV.MODE..HNZ from the ESM processing (m/s**2).
ESM_PGA_HNZ = 0.3273


@pytest.fixture
def client(server):
    return Client(server.url, retries=1, timeout=10.0)


@pytest.fixture(autouse=True)
def mseed_backend(monkeypatch):
    """Records read with the default backend of signals.io."""
    try:
        import cymseed3  # noqa: F401
    except ImportError:
        monkeypatch.setattr(io, "USE_LIBMSEED", False)


def _records(collection: StreamCollection) -> list:
    return [record for stream in collection for record in stream]


# ----------------------------------------------------------------------------
# Client
# ----------------------------------------------------------------------------


def test_client_from_registry():
    """Registry names, aliases and explicit services."""
    client = Client("ingv")
    assert client.name == "INGV"
    assert repr(client) == "Client('INGV')"
    assert client.services["dataselect"] == (
        "https://webservices.ingv.it/fdsnws/dataselect/1/"
    )

    client = Client("IRIS", event="https://example.org/fdsnws/event/1")
    assert client.name == "EarthScope"
    assert client.services["event"] == "https://example.org/fdsnws/event/1/"


def test_client_from_services():
    """Services only, without a data center."""
    client = Client(station="https://example.org/fdsnws/station/1/query")
    assert client.services == {
        "station": "https://example.org/fdsnws/station/1/"
    }
    with pytest.raises(ValueError, match="has no dataselect"):
        client.get_waveforms("IV.MODE", T0, T1)


@pytest.mark.parametrize(
    "args, kwargs, error",
    [
        ((), {}, ValueError),
        (("NOWHERE",), {}, ValueError),
        ((1,), {}, TypeError),
        (("INGV",), {"timeout": 0}, ValueError),
        (("INGV",), {"timeout": "10"}, TypeError),
        (("INGV",), {"retries": -1}, ValueError),
        (("INGV",), {"retries": 1.0}, TypeError),
        ((), {"station": "ftp://example.org"}, ValueError),
    ],
)
def test_client_errors(args, kwargs, error):
    """Invalid data centers and options are rejected."""
    with pytest.raises(error):
        Client(*args, **kwargs)


def test_version(client, server):
    """The version method of a service."""
    assert client.version("station") == "1.1.0"
    assert server.requests[-1].path == "/fdsnws/station/1/version"
    with pytest.raises(ValueError):
        client.version("seedlink")


# ----------------------------------------------------------------------------
# Waveforms
# ----------------------------------------------------------------------------


def test_waveforms_get_request(client, server):
    """One pattern: GET with the codes of the server."""
    collection = client.get_waveforms("IV.MODE..HN?", T0, T1)

    query = server.requests[-1].query
    assert server.requests[-1].method == "GET"
    assert query == {
        "network": "IV",
        "station": "MODE",
        "location": "--",
        "channel": "HN?",
        "starttime": "2012-05-20T02:04:00.000000",
        "endtime": "2012-05-20T02:05:30.000000",
    }
    records = _records(collection)
    assert [record.sid for record in records] == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ"
    ]
    for record in records:
        assert record.time == Date(T0)
        assert record.nsamp == 9001


def test_waveforms_are_filtered_and_cut(client):
    """The client keeps only what matches the pattern and the window."""
    collection = client.get_waveforms(
        "IV.MODE..HN[ZE]", "2012-05-20T02:04:00.004", T1
    )
    records = _records(collection)
    assert [record.sid for record in records] == [
        "IV.MODE..HNE", "IV.MODE..HNZ"
    ]
    assert records[0].time == Date("2012-05-20T02:04:00.01")


def test_waveforms_match_raw_data(client):
    """The samples are those of the file."""
    collection = client.get_waveforms("IV.MODE..HNZ", T0, T1)
    expected = io.reader(
        os.path.join(SIGNALS_DATA, "MODE.mseed"),
        select="IV.MODE..HNZ", starttime=T0, endtime=T1,
    )
    np.testing.assert_array_equal(
        _records(collection)[0].data, _records(expected)[0].data
    )


def test_waveforms_bulk_post(client, server):
    """Several patterns or windows: POST, one line each."""
    collection = client.get_waveforms(
        ["IV.MODE..HNZ", ("IV.MODE..HNE", "2012-05-20T02:05:00", None)],
        T0, T1, quality="D", longestonly=True,
    )
    request = server.requests[-1]
    assert request.method == "POST"
    assert request.body.splitlines() == [
        "quality=D",
        "longestonly=true",
        "IV MODE -- HNZ 2012-05-20T02:04:00.000000 "
        "2012-05-20T02:05:30.000000",
        "IV MODE -- HNE 2012-05-20T02:05:00.000000 "
        "2012-05-20T02:05:30.000000",
    ]
    records = _records(collection)
    assert [(r.sid, r.nsamp) for r in records] == [
        ("IV.MODE..HNE", 3001), ("IV.MODE..HNZ", 9001)
    ]


def _window_reply(tmp_path):
    """Reply with the records of MODE.mseed cut to the window asked."""
    def reply(request):
        query = request.query
        collection = io.reader(
            os.path.join(SIGNALS_DATA, "MODE.mseed"),
            select=f"IV.MODE..{query['channel']}",
            starttime=query["starttime"], endtime=query["endtime"],
        )
        path = tmp_path / f"reply{len(os.listdir(tmp_path))}.mseed"
        io.writer(collection, path, format="mseed")
        return 200, path.read_bytes(), {}
    return reply


def test_waveforms_chunks(client, server, tmp_path):
    """Long windows are split; the pieces are joined again."""
    server.replies["/fdsnws/dataselect/1/query"] = _window_reply(tmp_path)
    collection = client.get_waveforms("IV.MODE..HNZ", T0, T1, chunk=40.0)
    windows = [
        (r.query["starttime"], r.query["endtime"]) for r in server.requests
    ]
    assert windows == [
        ("2012-05-20T02:04:00.000000", "2012-05-20T02:04:40.000000"),
        ("2012-05-20T02:04:40.000000", "2012-05-20T02:05:20.000000"),
        ("2012-05-20T02:05:20.000000", "2012-05-20T02:05:30.000000"),
    ]
    (record,) = _records(collection)
    assert record.nsamp == 9001


def test_chunks_with_overlapping_records(client):
    """
    Whole records returned by consecutive chunks are not duplicated.

    The fake server returns the whole file to every chunk, as a server
    returns the MiniSEED records that overlap each chunk.
    """
    whole = client.get_waveforms("IV.MODE..HNZ", T0, T1)
    chunked = client.get_waveforms("IV.MODE..HNZ", T0, T1, chunk=25.0)
    (expected,) = _records(whole)
    (record,) = _records(chunked)
    assert record.time == expected.time
    np.testing.assert_array_equal(record.data, expected.data)


def test_waveforms_into_collection(client):
    """New records are added to a given collection."""
    collection = StreamCollection()
    collection.append(Record(data=np.zeros(5), delta=1.0, sid="XX.A..BHZ"))
    result = client.get_waveforms(
        "IV.MODE..HNZ", T0, T1, stream_collection=collection
    )
    assert result is collection
    assert sorted(r.sid for r in _records(collection)) == [
        "IV.MODE..HNZ", "XX.A..BHZ"
    ]


def test_waveforms_without_data(client, server):
    """No data: an empty collection."""
    server.replies["/fdsnws/dataselect/1/query"] = lambda req: (204, b"", {})
    collection = client.get_waveforms("IV.MODE..HNZ", T0, T1)
    assert len(_records(collection)) == 0


def test_waveforms_with_metadata(client, server):
    """Metadata of the channels are requested and assigned."""
    collection = client.get_waveforms("IV.MODE..HNZ", T0, T1, metadata=True)
    request = server.requests[-1]
    assert request.path == "/fdsnws/station/1/query"
    assert request.query["level"] == "response"
    assert request.query["location"] == "--"

    (record,) = _records(collection)
    assert record.head.meta.sid == "IV.MODE..HNZ"
    assert record.response is not None
    assert record.location.latitude == pytest.approx(44.629722)


def test_waveforms_corrected(client, server):
    """
    Corrected records: margins requested, physical units, PGA of ESM.

    The window keeps the strong motion of IV.MODE..HNZ; the PGA must
    agree with the ESM processing within 0.5 %.
    """
    collection = client.get_waveforms(
        "IV.MODE..HNZ", T0, T1, correct=True, pad=30.0
    )
    data_request = server.requests[-2]
    assert data_request.query["starttime"] == "2012-05-20T02:03:30.000000"
    assert data_request.query["endtime"] == "2012-05-20T02:06:00.000000"

    (record,) = _records(collection)
    assert record.time == Date(T0)
    assert record.nsamp == 9001
    assert record.head.units == "m/s**2"
    pga = np.abs(record.data).max()
    assert pga == pytest.approx(ESM_PGA_HNZ, rel=5e-3)


def test_waveforms_corrected_and_filtered(client):
    """Prefilter and band-pass filter keep the PGA close to ESM."""
    collection = client.get_waveforms(
        "IV.MODE..HNZ", T0, T1, correct=True, pad=30.0,
        prefilter=(0.02, 0.05, 40.0, 45.0), highpass=0.05, lowpass=40.0,
    )
    (record,) = _records(collection)
    pga = np.abs(record.data).max()
    assert pga == pytest.approx(ESM_PGA_HNZ, rel=1e-2)


def test_default_pad(client, server):
    """Default margin: 10 % of the window, at least 30 s."""
    client.get_waveforms("IV.MODE..HNZ", T0, T1, highpass=0.1)
    query = server.requests[-1].query
    assert query["starttime"] == "2012-05-20T02:03:30.000000"


def test_uncorrectable_records_are_left_out(client, server):
    """Records without a response are dropped with a warning."""
    with open(os.path.join(SIGNALS_DATA, "MODE.xml"), "rb") as fobj:
        xml = fobj.read().replace(b'code="HNE"', b'code="HXE"')
    server.replies["/fdsnws/station/1/query"] = lambda req: (200, xml, {})

    with pytest.warns(FdsnWarning, match="IV.MODE..HNE"):
        collection = client.get_waveforms(
            "IV.MODE..HN?", T0, T1, correct=True, pad=30.0
        )
    assert [r.sid for r in _records(collection)] == [
        "IV.MODE..HNN", "IV.MODE..HNZ"
    ]


@pytest.mark.parametrize(
    "kwargs, error",
    [
        ({"prefilter": (0.1, 0.2, 10.0, 20.0)}, ValueError),
        ({"pad": 10.0}, ValueError),
        ({"correct": True, "pad": -1.0}, ValueError),
        ({"correct": True, "pad": "1"}, TypeError),
        ({"chunk": 0}, ValueError),
        ({"correct": 1}, TypeError),
        ({"stream_collection": []}, TypeError),
    ],
)
def test_waveform_option_errors(client, kwargs, error):
    """Inconsistent or wrong options are rejected."""
    with pytest.raises(error):
        client.get_waveforms("IV.MODE..HNZ", T0, T1, **kwargs)


def test_waveforms_need_a_window(client):
    """Every pattern needs starttime and endtime."""
    with pytest.raises(ValueError, match="need starttime"):
        client.get_waveforms("IV.MODE..HNZ", T0)


# ----------------------------------------------------------------------------
# Raw downloads
# ----------------------------------------------------------------------------


def test_waveforms_to_file(client, server, tmp_path):
    """The records are written as received."""
    target = tmp_path / "mode.mseed"
    server.compress = True
    size = client.get_waveforms("IV.MODE", T0, T1, target=target)

    raw = read_file(SIGNALS_DATA, "MODE.mseed")
    assert size == len(raw)
    assert target.read_bytes() == raw

    with pytest.raises(FileExistsError):
        client.get_waveforms("IV.MODE", T0, T1, target=target)
    assert client.get_waveforms(
        "IV.MODE", T0, T1, target=target, owrite=True
    ) == len(raw)


def test_chunks_to_file(client, tmp_path):
    """The replies of the chunks are concatenated in one file."""
    target = tmp_path / "mode.mseed"
    raw = read_file(SIGNALS_DATA, "MODE.mseed")
    size = client.get_waveforms("IV.MODE", T0, T1, chunk=60.0,
                                target=target)
    assert size == 2 * len(raw)
    assert target.read_bytes() == raw + raw
    assert os.listdir(tmp_path) == ["mode.mseed"]


def test_raw_download_options(client, tmp_path):
    """Raw downloads cannot be processed."""
    with pytest.raises(ValueError, match="as received"):
        client.get_waveforms("IV.MODE", T0, T1, target=tmp_path / "x",
                             metadata=True)
    with pytest.raises(TypeError):
        client.get_waveforms("IV.MODE", T0, T1, target=3)


def test_metadata_to_file(client, tmp_path):
    """StationXML written as received."""
    target = tmp_path / "mode.xml"
    size = client.get_metadata("IV.MODE", level="response", target=target)
    raw = read_file(SIGNALS_DATA, "MODE.xml")
    assert size == len(raw)
    assert target.read_bytes() == raw


# ----------------------------------------------------------------------------
# Metadata
# ----------------------------------------------------------------------------


def test_metadata(client, server):
    """Channel level by default, filtered by the client."""
    collection = client.get_metadata("IV.MODE..HN[ZN]", T0)
    assert isinstance(collection, MetadataCollection)
    assert sorted(collection.sid) == ["IV.MODE..HNN", "IV.MODE..HNZ"]
    query = server.requests[-1].query
    assert query["level"] == "channel"
    assert query["channel"] == "HNZ,HNN"
    assert query["starttime"] == "2012-05-20T02:04:00.000000"
    assert "endtime" not in query


def test_metadata_text(client, server):
    """Text format at channel and station level."""
    channels = client.get_metadata("IV.MODE..HNZ", text=True)
    assert server.requests[-1].query["format"] == "text"
    assert channels.sid == ["IV.MODE..HNZ"]

    stations = client.get_metadata("IV.MOX*", level="station", text=True)
    assert [s.code for s in stations.stations()] == ["IV.MOXX"]


def test_metadata_station_level(client):
    """Station epochs only."""
    collection = client.get_metadata("IV.MODE", level="station")
    assert len(collection) == 0
    assert [s.code for s in collection.stations()] == ["IV.MODE"]


def test_metadata_geographic(client, server):
    """Regions in the order of ShakeLab, sent with FDSN names."""
    client.get_metadata(region=(10.0, 12.0, 44.0, 45.5))
    query = server.requests[-1].query
    assert (query["minlongitude"], query["maxlongitude"]) == ("10", "12")
    assert (query["minlatitude"], query["maxlatitude"]) == ("44", "45.5")

    client.get_metadata(center=(10.95, 44.63), radius=(None, 0.5),
                        includerestricted=False)
    query = server.requests[-1].query
    assert query["longitude"] == "10.95"
    assert query["latitude"] == "44.63"
    assert query["maxradius"] == "0.5"
    assert "minradius" not in query
    assert query["includerestricted"] == "false"


@pytest.mark.parametrize(
    "kwargs, error",
    [
        ({"level": "network"}, ValueError),
        ({"level": "response", "text": True}, ValueError),
        ({"region": (1.0, 2.0, 3.0)}, TypeError),
        ({"radius": (0.0, 1.0)}, ValueError),
        ({"center": (1.0, 2.0), "radius": (2.0, 1.0)}, ValueError),
        ({"matchtimeseries": "yes"}, TypeError),
    ],
)
def test_metadata_errors(client, kwargs, error):
    """Invalid options are rejected."""
    with pytest.raises(error):
        client.get_metadata(**kwargs)


def test_metadata_without_data(client, server):
    """No metadata: an empty collection."""
    server.replies["/fdsnws/station/1/query"] = lambda req: (204, b"", {})
    assert len(client.get_metadata("XX")) == 0


# ----------------------------------------------------------------------------
# Events
# ----------------------------------------------------------------------------


def test_events(client, server):
    """QuakeML catalogue, depths in km."""
    catalogue = client.get_events(
        "2012-05-01T00:00:00", "2012-06-30T00:00:00",
        magnitude=(5.5, None), depth=(0, 30), orderby="time-asc",
        limit=10,
    )
    query = server.requests[-1].query
    assert query == {
        "starttime": "2012-05-01T00:00:00.000000",
        "endtime": "2012-06-30T00:00:00.000000",
        "mindepth": "0",
        "maxdepth": "30",
        "minmagnitude": "5.5",
        "orderby": "time-asc",
        "limit": "10",
    }
    assert len(catalogue.event) == 2


def test_events_without_data(client, server):
    """No events: an empty catalogue."""
    server.replies["/fdsnws/event/1/query"] = lambda req: (204, b"", {})
    assert len(client.get_events().event) == 0


@pytest.mark.parametrize(
    "kwargs, error",
    [
        ({"magnitude": 5.0}, TypeError),
        ({"limit": 0}, ValueError),
        ({"catalog": 1}, TypeError),
        ({"includearrivals": 1}, TypeError),
    ],
)
def test_event_errors(client, kwargs, error):
    """Invalid options are rejected."""
    with pytest.raises(error):
        client.get_events(**kwargs)


# ----------------------------------------------------------------------------
# Availability
# ----------------------------------------------------------------------------


def test_availability(client, server):
    """Spans of the matching channels, JSON format."""
    spans = client.get_availability("IV.MODE..HNZ", T0, T1,
                                    merge=["quality", "samplerate"])
    query = server.requests[-1].query
    assert query["format"] == "json"
    assert query["merge"] == "quality,samplerate"
    assert server.requests[-1].path == "/fdsnws/availability/1/query"
    assert spans == [
        TimeSpan("IV.MODE..HNZ", Date("2012-05-20T02:03:23"),
                 Date("2012-05-20T02:06:23"), "D", 100.0),
        TimeSpan("IV.MODE..HNZ", Date("2012-05-20T02:10:00"),
                 Date("2012-05-20T02:12:00"), "D", 100.0),
    ]


def test_availability_extent(client, server):
    """First and last sample of each channel."""
    spans = client.get_availability("IV.MODE", T0, T1, extent=True)
    assert server.requests[-1].path == "/fdsnws/availability/1/extent"
    assert spans == [
        TimeSpan("IV.MODE..HNZ", Date("2012-05-20T02:03:23"),
                 Date("2012-05-20T02:12:00"), "D", 100.0),
    ]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"merge": "gaps"},
        {"extent": True, "merge": "overlap"},
        {"extent": True, "mergegaps": 1.0},
    ],
)
def test_availability_errors(client, kwargs):
    """Invalid merge options are rejected."""
    with pytest.raises(ValueError):
        client.get_availability("IV.MODE", T0, T1, **kwargs)


def test_invalid_availability_reply(client, server):
    """A malformed reply raises FdsnError."""
    server.replies["/fdsnws/availability/1/query"] = lambda req: (
        200, b"not json", {}
    )
    with pytest.raises(FdsnError, match="invalid availability"):
        client.get_availability("IV.MODE", T0, T1)


def test_server_errors_reach_the_user(client, server):
    """An error of the service is raised with its message."""
    server.replies["/fdsnws/dataselect/1/query"] = lambda req: (
        400, b"Error 400: Bad Request\nstarttime is invalid\n", {}
    )
    with pytest.raises(FdsnError, match="Error 400") as info:
        client.get_waveforms("IV.MODE", T0, T1)
    assert "starttime is invalid" in info.value.detail
