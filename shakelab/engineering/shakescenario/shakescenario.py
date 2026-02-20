# shakescenario.py
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
ShakeLab - ShakeScenario client (CLI).

This module implements a CLI client for the shakescenario-server.

Operations
----------
- ping
- submit
- list
- get
- wait
- cancel
- delete
- reset

Notes
-----
This client speaks the length-prefixed JSON TCP protocol defined in
protocol.py.

"""

from __future__ import annotations

import argparse
import json
import time
import uuid
from pathlib import Path
from typing import Any

from protocol import ProtocolClient


def _json_load_file(path: str) -> Any:
    p = Path(path).expanduser().resolve()
    return json.loads(p.read_text(encoding="utf-8"))


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="shakescenario",
        description="ShakeScenario client CLI (TCP).",
    )
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=6000)

    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("ping", help="Ping server and print capabilities.")

    p_submit = sub.add_parser("submit", help="Submit a new scenario job.")
    p_submit.add_argument("--params-json", help="Path to a JSON params file.")
    p_submit.add_argument("--tag", default=None)
    # Convenience event args (optional; merged into params)
    p_submit.add_argument("--mag", type=float, default=None)
    p_submit.add_argument("--lon", type=float, default=None)
    p_submit.add_argument("--lat", type=float, default=None)
    p_submit.add_argument("--depth", type=float, default=None)
    p_submit.add_argument("--strike", type=float, default=None)
    p_submit.add_argument("--dip", type=float, default=None)
    p_submit.add_argument("--rake", type=float, default=None)

    p_list = sub.add_parser("list", help="List jobs.")
    p_list.add_argument("--status", default=None)
    p_list.add_argument("--limit", type=int, default=50)
    p_list.add_argument("--offset", type=int, default=0)

    p_get = sub.add_parser("get", help="Get job details.")
    p_get.add_argument("job_id", type=int)

    p_wait = sub.add_parser("wait", help="Wait for job completion.")
    p_wait.add_argument("job_id", type=int)
    p_wait.add_argument("--timeout", type=float, default=600.0)
    p_wait.add_argument("--interval", type=float, default=2.0)

    p_cancel = sub.add_parser("cancel", help="Cancel a job (best effort).")
    p_cancel.add_argument("job_id", type=int)

    p_delete = sub.add_parser("delete", help="Delete job from database.")
    p_delete.add_argument("job_id", type=int)
    p_delete.add_argument("--purge", action="store_true")

    p_reset = sub.add_parser("reset", help="Reset database (dangerous).")
    p_reset.add_argument("--yes-i-know", action="store_true")

    return p


def _req_id() -> str:
    return str(uuid.uuid4())


def _print_json(obj: Any) -> None:
    print(json.dumps(obj, indent=2, sort_keys=True))


def main() -> None:
    args = _build_parser().parse_args()
    client = ProtocolClient(host=args.host, port=args.port, timeout=10.0)

    if args.cmd == "ping":
        res = client.request("ping", payload={}, req_id=_req_id())
        _print_json(res)
        return

    if args.cmd == "submit":
        params: dict[str, Any] = {}

        if args.params_json:
            loaded = _json_load_file(args.params_json)
            if not isinstance(loaded, dict):
                raise ValueError("params-json must contain a JSON object.")
            params.update(loaded)

        event = params.get("event", {})
        if not isinstance(event, dict):
            event = {}
        for k, v in (
            ("magnitude", args.mag),
            ("longitude", args.lon),
            ("latitude", args.lat),
            ("depth_km", args.depth),
            ("strike", args.strike),
            ("dip", args.dip),
            ("rake", args.rake),
        ):
            if v is not None:
                event[k] = v

        if event:
            params["event"] = event

        payload = {"params": params}
        if args.tag is not None:
            payload["tag"] = args.tag

        res = client.request("submit", payload=payload, req_id=_req_id())
        _print_json(res)
        return

    if args.cmd == "list":
        payload = {
            "status": args.status,
            "limit": args.limit,
            "offset": args.offset,
        }
        res = client.request("list", payload=payload, req_id=_req_id())
        _print_json(res)
        return

    if args.cmd == "get":
        res = client.request(
            "get",
            payload={"job_id": args.job_id},
            req_id=_req_id(),
        )
        _print_json(res)
        return

    if args.cmd == "wait":
        t0 = time.time()
        while True:
            res = client.request(
                "get",
                payload={"job_id": args.job_id},
                req_id=_req_id(),
            )
            status = res.get("status")
            if status in ("completed", "failed", "canceled"):
                _print_json(res)
                return

            if (time.time() - t0) > args.timeout:
                raise TimeoutError("Timeout waiting for job completion.")

            time.sleep(max(0.1, args.interval))

    if args.cmd == "cancel":
        res = client.request(
            "cancel",
            payload={"job_id": args.job_id},
            req_id=_req_id(),
        )
        _print_json(res)
        return

    if args.cmd == "delete":
        res = client.request(
            "delete",
            payload={"job_id": args.job_id, "purge": args.purge},
            req_id=_req_id(),
        )
        _print_json(res)
        return

    if args.cmd == "reset":
        if not args.yes_i_know:
            raise ValueError("Use --yes-i-know to confirm reset.")
        res = client.request(
            "reset",
            payload={"confirm": True},
            req_id=_req_id(),
        )
        _print_json(res)
        return


if __name__ == "__main__":
    main()