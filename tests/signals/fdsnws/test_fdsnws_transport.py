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
Tests of :mod:`shakelab.signals.fdsnws.transport`.

The requests go to the fake server of ``fdsn_server.py``: GET and POST,
compressed replies, retries of busy servers, errors with the message
of the server, redirects of POST requests and replies written to files.
"""

from __future__ import annotations

import os

import pytest

from shakelab.libutils.timeN import Date
from shakelab.signals.fdsnws import transport
from shakelab.signals.fdsnws.transport import FdsnError, format_value, request

from fdsn_server import SIGNALS_DATA


def _mseed() -> bytes:
    with open(os.path.join(SIGNALS_DATA, "MODE.mseed"), "rb") as fobj:
        return fobj.read()


# ----------------------------------------------------------------------------
# Parameter values
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected",
    [
        (Date("2012-05-20T02:03:23.5Z"), "2012-05-20T02:03:23.500000"),
        (True, "true"),
        (False, "false"),
        (7, "7"),
        (2.0, "2"),
        (1e-05, "0.00001"),
        (-0.5, "-0.5"),
        (123456789.125, "123456789.125"),
        (1e20, "100000000000000000000"),
        ("HN?", "HN?"),
    ],
)
def test_format_value(value, expected):
    """Times without zone, booleans in lower case, decimal reals."""
    assert format_value(value) == expected


def test_format_value_errors():
    """Non-finite reals and other types are rejected."""
    with pytest.raises(ValueError):
        format_value(float("nan"))
    with pytest.raises(TypeError):
        format_value([1, 2])


# ----------------------------------------------------------------------------
# Requests
# ----------------------------------------------------------------------------


def test_get_request(server):
    """Parameters go in the URL; None values are left out."""
    url = server.url + "/fdsnws/dataselect/1/query"
    content = request(url, params={"network": "IV", "quality": None,
                                   "channel": "HN?,HH?"})
    assert content == _mseed()

    received = server.requests[-1]
    assert received.method == "GET"
    assert received.query == {"network": "IV", "channel": "HN?,HH?"}
    assert received.headers["User-Agent"].startswith("ShakeLab/")
    assert received.headers["Accept-Encoding"] == "gzip"


def test_post_request(server):
    """The body is sent as text."""
    url = server.url + "/fdsnws/dataselect/1/query"
    body = "quality=B\nIV MODE -- HNZ * *\n"
    assert request(url, body=body) == _mseed()

    received = server.requests[-1]
    assert received.method == "POST"
    assert received.body == body
    assert received.headers["Content-Type"] == "text/plain"


def test_params_and_body_are_exclusive(server):
    """A request has either parameters or a body."""
    with pytest.raises(ValueError):
        request(server.url, params={"a": "b"}, body="x")


def test_compressed_reply(server):
    """Gzip replies are decompressed."""
    server.compress = True
    url = server.url + "/fdsnws/dataselect/1/query"
    assert request(url) == _mseed()


def test_no_data(server):
    """Status 204 means no data."""
    server.replies["/empty"] = lambda req: (204, b"", {})
    assert request(server.url + "/empty") is None


def test_redirected_post(server):
    """A POST redirected with 307 is sent again with its body."""
    url = server.url + "/moved/fdsnws/dataselect/1/query"
    body = "IV MODE -- HNZ * *\n"
    assert request(url, body=body) == _mseed()
    assert [r.method for r in server.requests] == ["POST", "POST"]
    assert server.requests[-1].body == body


# ----------------------------------------------------------------------------
# Errors and retries
# ----------------------------------------------------------------------------


def test_busy_server_is_retried(server):
    """Statuses 503 are retried, honoring Retry-After."""
    server.failures = [503, 503]
    server.retry_after = "0"
    url = server.url + "/fdsnws/dataselect/1/query"
    assert request(url, retries=2) == _mseed()
    assert len(server.requests) == 3


def test_retries_exhausted(server):
    """After the retries, the last error is raised."""
    server.failures = [503, 503, 503]
    url = server.url + "/fdsnws/dataselect/1/query"
    with pytest.raises(FdsnError) as info:
        request(url, retries=1)
    assert info.value.status == 503
    assert len(server.requests) == 2


def test_client_errors_are_not_retried(server):
    """A 400 is raised at once, with the message of the server."""
    server.replies["/bad"] = lambda req: (
        400, b"Error 400: Bad Request\n\nUnknown parameter: foo\n", {}
    )
    with pytest.raises(FdsnError) as info:
        request(server.url + "/bad", retries=3)
    assert info.value.status == 400
    assert "Unknown parameter: foo" in info.value.detail
    assert "Error 400" in str(info.value)
    assert len(server.requests) == 1


def test_too_large_request_hint(server):
    """A 413 suggests splitting the request."""
    server.replies["/large"] = lambda req: (413, b"Too much data\n", {})
    with pytest.raises(FdsnError, match="chunk"):
        request(server.url + "/large")


def test_unreachable_server():
    """Connection failures become FdsnError after the retries."""
    with pytest.raises(FdsnError, match="cannot reach"):
        request("http://127.0.0.1:9/fdsnws/station/1/query", retries=1,
                timeout=2)


# ----------------------------------------------------------------------------
# Replies written to files
# ----------------------------------------------------------------------------


def test_reply_to_file(server, tmp_path):
    """The reply is written as received; no partial file remains."""
    server.compress = True
    target = tmp_path / "data.mseed"
    url = server.url + "/fdsnws/dataselect/1/query"
    size = request(url, target=target)

    assert size == len(_mseed())
    assert target.read_bytes() == _mseed()
    assert os.listdir(tmp_path) == ["data.mseed"]


def test_no_data_writes_no_file(server, tmp_path):
    """Without data, no file is created."""
    server.replies["/empty"] = lambda req: (204, b"", {})
    target = tmp_path / "data.mseed"
    assert request(server.url + "/empty", target=target) is None
    assert os.listdir(tmp_path) == []


def test_failed_download_leaves_no_file(server, tmp_path):
    """An error leaves neither the file nor its temporary copy."""
    server.failures = [500]
    target = tmp_path / "data.mseed"
    with pytest.raises(FdsnError):
        request(server.url + "/fdsnws/dataselect/1/query", target=target,
                retries=0)
    assert os.listdir(tmp_path) == []


def test_user_agent():
    """The User-Agent names ShakeLab and its version."""
    assert transport.user_agent().startswith("ShakeLab/")
