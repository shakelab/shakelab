# protocol.py
# -*- coding: utf-8 -*-
# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
ShakeLab - Engineering - ShakeScenario TCP protocol helpers.

This module defines a simple framed protocol:

- Each message is a JSON document encoded as UTF-8 bytes.
- It is preceded by an 8-byte unsigned big-endian length header.

This avoids ambiguity with TCP segmentation (recv() may return partial data).

"""

from __future__ import annotations

import json
import socket
import struct
from dataclasses import dataclass
from typing import Any


_LEN_STRUCT = struct.Struct(">Q")  # 8-byte unsigned length (big-endian)


class ProtocolError(RuntimeError):
    """Raised on protocol framing/decoding errors."""


def _recvall(sock: socket.socket, n: int) -> bytes:
    """
    Receive exactly n bytes from a socket.

    Raises
    ------
    EOFError
        If the peer closes the connection before n bytes are read.
    """
    chunks: list[bytes] = []
    remaining = n
    while remaining > 0:
        data = sock.recv(remaining)
        if not data:
            raise EOFError("Connection closed by peer.")
        chunks.append(data)
        remaining -= len(data)
    return b"".join(chunks)


def send_message(sock: socket.socket, obj: dict[str, Any]) -> None:
    """
    Send a framed JSON message.

    Parameters
    ----------
    sock
        Connected TCP socket.
    obj
        JSON-serializable dictionary.
    """
    payload = json.dumps(obj, separators=(",", ":"), ensure_ascii=False).encode(
        "utf-8"
    )
    header = _LEN_STRUCT.pack(len(payload))
    sock.sendall(header + payload)


def recv_message(sock: socket.socket) -> dict[str, Any]:
    """
    Receive a framed JSON message.

    Returns
    -------
    dict
        Decoded JSON object.

    Raises
    ------
    EOFError
        If the connection closes.
    ProtocolError
        If the message cannot be decoded or is not a JSON object.
    """
    header = _recvall(sock, _LEN_STRUCT.size)
    (size,) = _LEN_STRUCT.unpack(header)

    if size == 0:
        raise ProtocolError("Empty message.")
    if size > 128 * 1024 * 1024:
        raise ProtocolError("Message too large.")

    body = _recvall(sock, int(size))
    try:
        obj = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProtocolError(str(exc)) from exc

    if not isinstance(obj, dict):
        raise ProtocolError("Message must be a JSON object.")
    return obj


@dataclass(frozen=True)
class ProtocolClient:
    """
    Small helper for request/response interaction.

    Parameters
    ----------
    host, port
        Server address.
    timeout
        Socket timeout (seconds).

    """

    host: str
    port: int
    timeout: float = 10.0
    api_version: int = 1

    def request(
        self,
        op: str,
        payload: dict[str, Any],
        req_id: str,
    ) -> dict[str, Any]:
        """
        Send a request and wait for its response.

        Returns the "result" dict if ok=True, otherwise raises RuntimeError.
        """
        req = {"v": self.api_version, "op": op, "req_id": req_id, "payload": payload}

        with socket.create_connection((self.host, self.port), self.timeout) as sock:
            sock.settimeout(self.timeout)
            send_message(sock, req)
            resp = recv_message(sock)

        if not resp.get("ok", False):
            err = resp.get("error", {})
            code = err.get("code", "ERROR")
            msg = err.get("message", "Unknown error.")
            raise RuntimeError(f"{code}: {msg}")

        result = resp.get("result", {})
        if isinstance(result, dict):
            return result
        return {"result": result}

