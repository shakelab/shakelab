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
Fake FDSN web services on a local socket, for the tests of
:mod:`shakelab.signals.fdsnws`.

The server answers the ``query``, ``extent`` and ``version`` methods of
the dataselect, station, event and availability services under the
standard paths, with the files of ``tests/signals/data`` and
``tests/signals/fdsnws/data``; it records every request. Its behavior
can be changed per test:

- ``replies``: a path mapped to a function ``f(request) -> (status,
  body, headers)``, which replaces the default reply;
- ``failures``: statuses returned, in order, before the normal replies
  (e.g. ``[503, 503]``), with an optional ``Retry-After``;
- ``compress``: gzip the replies;
- ``/moved/...``: answered with a 307 redirect to ``/...``.

The data replies ignore the selection (the client filters what it
receives), except that ``format=text`` selects the text files and
``level=station`` the station list.
"""

from __future__ import annotations

import gzip
import os
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Callable, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
SIGNALS_DATA = os.path.join(HERE, "..", "data")
DATA = os.path.join(HERE, "data")


def read_file(*parts: str) -> bytes:
    """Content of a file."""
    with open(os.path.join(*parts), "rb") as fobj:
        return fobj.read()


class FakeRequest:
    """One request received by the server."""

    def __init__(self, method: str, path: str, query: dict, body: str,
                 headers: dict) -> None:
        self.method = method
        self.path = path
        self.query = query
        self.body = body
        self.headers = headers

    def __repr__(self) -> str:
        return f"FakeRequest({self.method} {self.path} {self.query})"


class FakeFdsnServer:
    """Fake FDSN server on ``127.0.0.1`` and a free port."""

    def __init__(self) -> None:
        self.requests: list = []
        self.replies: dict = {}
        self.failures: list = []
        self.retry_after: Optional[str] = None
        self.compress = False
        self._lock = threading.Lock()

        server = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args: Any) -> None:
                pass

            def do_GET(self) -> None:
                server._handle(self, "GET")

            def do_POST(self) -> None:
                server._handle(self, "POST")

        self._httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._thread = threading.Thread(
            target=self._httpd.serve_forever,
            kwargs={"poll_interval": 0.02},
            daemon=True,
        )
        self._thread.start()

    @property
    def url(self) -> str:
        """Address of the server, ``http://127.0.0.1:<port>``."""
        return f"http://127.0.0.1:{self._httpd.server_address[1]}"

    def close(self) -> None:
        """Stop the server."""
        self._httpd.shutdown()
        self._httpd.server_close()

    # -- replies --------------------------------------------------------

    def _handle(self, handler: BaseHTTPRequestHandler, method: str) -> None:
        parsed = urllib.parse.urlsplit(handler.path)
        query = {
            key: values[0]
            for key, values in urllib.parse.parse_qs(
                parsed.query, keep_blank_values=True
            ).items()
        }
        length = int(handler.headers.get("Content-Length") or 0)
        body = handler.rfile.read(length).decode("utf-8") if length else ""
        request = FakeRequest(method, parsed.path, query, body,
                              dict(handler.headers))
        with self._lock:
            self.requests.append(request)
            failure = self.failures.pop(0) if self.failures else None

        if failure is not None:
            headers = {}
            if self.retry_after is not None:
                headers["Retry-After"] = self.retry_after
            self._send(handler, failure, b"Service busy\n", headers)
            return

        if parsed.path.startswith("/moved/"):
            location = parsed.path[len("/moved"):]
            if parsed.query:
                location += "?" + parsed.query
            self._send(handler, 307, b"", {"Location": location})
            return

        reply: Callable = self.replies.get(parsed.path, self._default)
        status, content, headers = reply(request)
        self._send(handler, status, content, headers)

    def _send(self, handler: BaseHTTPRequestHandler, status: int,
              content: bytes, headers: dict) -> None:
        if self.compress and content:
            content = gzip.compress(content)
            headers = {**headers, "Content-Encoding": "gzip"}
        handler.send_response(status)
        for key, value in headers.items():
            handler.send_header(key, value)
        handler.send_header("Content-Length", str(len(content)))
        handler.end_headers()
        if status != 204:
            handler.wfile.write(content)

    def _default(self, request: FakeRequest) -> tuple:
        path = request.path
        params = dict(request.query)
        for line in request.body.splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                params[key.strip()] = value.strip()

        if path.endswith("/version"):
            return 200, b"1.1.0\n", {"Content-Type": "text/plain"}

        if path == "/fdsnws/dataselect/1/query":
            return (200, read_file(SIGNALS_DATA, "MODE.mseed"),
                    {"Content-Type": "application/vnd.fdsn.mseed"})

        if path == "/fdsnws/station/1/query":
            if params.get("format") == "text":
                name = (
                    "stations.txt" if params.get("level") == "station"
                    else "channels.txt"
                )
                return (200, read_file(DATA, name),
                        {"Content-Type": "text/plain"})
            return (200, read_file(SIGNALS_DATA, "MODE.xml"),
                    {"Content-Type": "application/xml"})

        if path == "/fdsnws/event/1/query":
            return (200, read_file(DATA, "events.xml"),
                    {"Content-Type": "application/xml"})

        if path == "/fdsnws/availability/1/query":
            return (200, read_file(DATA, "availability.json"),
                    {"Content-Type": "application/json"})

        if path == "/fdsnws/availability/1/extent":
            return (200, read_file(DATA, "extent.json"),
                    {"Content-Type": "application/json"})

        return 404, b"Error 404: Not Found\n", {"Content-Type": "text/plain"}
