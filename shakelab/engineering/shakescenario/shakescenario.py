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
"""ShakeScenario client (CLI + reusable API).

This module provides:
- :class:`ShakeScenarioClient` for programmatic use (returns Python objects).
- a thin CLI wrapper around it (prints human-friendly output by default).

The client speaks the length-prefixed JSON TCP protocol implemented in
:mod:`protocol`.

"""

from __future__ import annotations

import argparse
import json
import sys
import time
import uuid
from pathlib import Path
from typing import Any, Iterable, Sequence

from protocol import ProtocolClient


class ShakeScenarioClient:
    """High-level client for the ShakeScenario TCP API.

    Parameters
    ----------
    host, port
        Server address.
    timeout
        Socket timeout (seconds).
    """

    def __init__(self, host: str, port: int, timeout: float = 10.0) -> None:
        self._client = ProtocolClient(host=host, port=port, timeout=timeout)

    @staticmethod
    def _req_id() -> str:
        return str(uuid.uuid4())

    def ping(self) -> dict[str, Any]:
        return self._client.request("ping", payload={}, req_id=self._req_id())

    def submit(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self._client.request(
            "submit",
            payload=payload,
            req_id=self._req_id(),
        )

    def list_jobs(
        self,
        status: str | None = None,
        limit: int = 50,
        offset: int = 0,
    ) -> dict[str, Any]:
        return self._client.request(
            "list",
            payload={"status": status, "limit": int(limit), "offset": int(offset)},
            req_id=self._req_id(),
        )

    def models_list(self) -> dict[str, Any]:
        return self._client.request(
            "models.list",
            payload={},
            req_id=self._req_id(),
        )

    def get_job(self, job_id: int) -> dict[str, Any]:
        return self._client.request(
            "get",
            payload={"job_id": int(job_id)},
            req_id=self._req_id(),
        )

    def cancel_job(self, job_id: int) -> dict[str, Any]:
        return self._client.request(
            "cancel",
            payload={"job_id": int(job_id)},
            req_id=self._req_id(),
        )

    def delete_job(self, job_id: int, purge: bool = False) -> dict[str, Any]:
        return self._client.request(
            "delete",
            payload={"job_id": int(job_id), "purge": bool(purge)},
            req_id=self._req_id(),
        )

    def reset(self) -> dict[str, Any]:
        return self._client.request(
            "reset",
            payload={"confirm": True},
            req_id=self._req_id(),
        )

    def wait_job(
        self,
        job_id: int,
        timeout: float = 600.0,
        interval: float = 2.0,
    ) -> dict[str, Any]:
        """Poll job status until completion/failure/cancel or timeout."""
        t0 = time.time()
        last_status: str | None = None
        while True:
            job = self.get_job(job_id)
            status = str(job.get("status", ""))

            if status and status != last_status:
                last_status = status

            if status in ("completed", "failed", "canceled"):
                return job

            if (time.time() - t0) > float(timeout):
                raise TimeoutError("Timeout waiting for job completion.")

            time.sleep(max(0.1, float(interval)))


def _safe_getattr(ns: argparse.Namespace, name: str, default: Any = None) -> Any:
    return getattr(ns, name, default)


def _json_load_file(path: str) -> Any:
    p = Path(path).expanduser().resolve()
    return json.loads(p.read_text(encoding="utf-8"))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="shakescenario",
        description="ShakeScenario client CLI (TCP).",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=6000)
    parser.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format (default: text).",
    )

    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("ping", help="Ping server and print capabilities.")

    p_submit = sub.add_parser("submit", help="Submit a new scenario job.")
    p_submit.add_argument("--tag", default=None)
    p_submit.add_argument("--request-json", dest="request_json", default=None)
    p_submit.add_argument("--model-id", dest="model_id", default=None)
    p_submit.add_argument("--gmpe", default=None)
    p_submit.add_argument("--distance-approx", dest="distance_approx", default=None)
    p_submit.add_argument("--mag", type=float, default=None)
    p_submit.add_argument("--lon", type=float, default=None)
    p_submit.add_argument("--lat", type=float, default=None)
    p_submit.add_argument("--depth", type=float, default=None)

    p_list = sub.add_parser("list", help="List jobs.")
    p_list.add_argument("--status", default=None)
    p_list.add_argument("--limit", type=int, default=50)
    p_list.add_argument("--offset", type=int, default=0)

    p_models = sub.add_parser("models", help="Model registry operations.")
    models_sub = p_models.add_subparsers(dest="models_cmd", required=True)
    models_sub.add_parser("list", help="List model_id available on server.")

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
    p_reset.add_argument("--yes-i-know", dest="yes_i_know", action="store_true")

    return parser


def _print_json(obj: Any) -> None:
    print(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False))


def _get_nested(obj: Any, path: Sequence[str], default: Any = None) -> Any:
    cur: Any = obj
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _format_float(value: Any, ndigits: int) -> str:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "-"
    return f"{f:.{ndigits}f}"


def _elevation_to_depth_km(elev_m: Any) -> str:
    try:
        elev = float(elev_m)
    except (TypeError, ValueError):
        return "-"
    return _format_float((-elev) / 1000.0, 2)


def _print_table(rows: list[list[str]], headers: list[str]) -> None:
    if not rows:
        print("No results.")
        return

    widths = [len(h) for h in headers]
    for r in rows:
        for i, cell in enumerate(r):
            widths[i] = max(widths[i], len(cell))

    def fmt_row(r: Iterable[str]) -> str:
        return "  ".join(cell.ljust(widths[i]) for i, cell in enumerate(r))

    print(fmt_row(headers))
    print(fmt_row(["-" * w for w in widths]))
    for r in rows:
        print(fmt_row(r))


def _job_key_fields(job: dict[str, Any]) -> dict[str, str]:
    params = job.get("params", {})
    model_id = _get_nested(params, ["models", "model_id"], "-")
    mag = _get_nested(params, ["scenario", "event", "magnitude"], None)
    lon = _get_nested(params, ["scenario", "event", "hypocentre", "longitude"], None)
    lat = _get_nested(params, ["scenario", "event", "hypocentre", "latitude"], None)
    elev = _get_nested(params, ["scenario", "event", "hypocentre", "elevation"], None)
    gmpe = _get_nested(params, ["scenario", "ground_motion", "gmpe_name"], "-")

    return {
        "model_id": str(model_id) if model_id is not None else "-",
        "mag": _format_float(mag, 2) if mag is not None else "-",
        "lon": _format_float(lon, 4) if lon is not None else "-",
        "lat": _format_float(lat, 4) if lat is not None else "-",
        "depth": _elevation_to_depth_km(elev),
        "gmpe": str(gmpe) if gmpe is not None else "-",
    }


def _print_ping_text(res: dict[str, Any]) -> None:
    print(f"Server: {res.get('name', '-')}")
    print(f"API version: {res.get('api_version', '-')}")
    print(f"Protocol: {res.get('protocol', '-')}")


def _print_models_list_text(res: dict[str, Any]) -> None:
    mids = res.get("model_ids", [])
    if not isinstance(mids, list) or not mids:
        print("No models.")
        return
    for m in mids:
        print(str(m))


def _print_list_text(res: dict[str, Any]) -> None:
    jobs = res.get("jobs", [])
    if not isinstance(jobs, list) or not jobs:
        print("No jobs.")
        return

    headers = [
        "id",
        "status",
        "created_at",
        "tag",
        "model_id",
        "mag",
        "lon",
        "lat",
        "depth_km",
    ]
    rows: list[list[str]] = []
    for j in jobs:
        if not isinstance(j, dict):
            continue
        k = _job_key_fields(j)
        rows.append(
            [
                str(j.get("id", "-")),
                str(j.get("status", "-")),
                str(j.get("created_at") or "-"),
                str(j.get("tag") or "-"),
                k["model_id"],
                k["mag"],
                k["lon"],
                k["lat"],
                k["depth"],
            ]
        )
    _print_table(rows, headers)


def _print_get_text(job: dict[str, Any]) -> None:
    k = _job_key_fields(job)

    print(f"Job {job.get('id', '-')}: {job.get('status', '-')}")
    if job.get("tag"):
        print(f"Tag: {job.get('tag')}")
    print(f"Created: {job.get('created_at') or '-'}")
    print(f"Started: {job.get('started_at') or '-'}")
    print(f"Ended: {job.get('ended_at') or '-'}")

    print("")
    print("Scenario")
    print(f"  model_id: {k['model_id']}")
    print(f"  magnitude: {k['mag']}")
    print(f"  lon/lat: {k['lon']} / {k['lat']}")
    print(f"  depth_km: {k['depth']}")
    if k["gmpe"] != "-":
        print(f"  gmpe: {k['gmpe']}")

    meta = job.get("result_meta", {})
    if isinstance(meta, dict) and meta:
        workdir = meta.get("workdir")
        manifest = meta.get("manifest")
        if workdir or manifest:
            print("")
            print("Artifacts")
        if workdir:
            print(f"  workdir: {workdir}")
        if isinstance(manifest, dict) and manifest:
            keys = ", ".join(sorted(manifest.keys()))
            print(f"  manifest keys: {keys}")

    err = job.get("error")
    if err:
        print("")
        print("Error")
        print(str(err))


def _print_submit_text(res: dict[str, Any]) -> None:
    job_id = res.get("job_id")
    status = res.get("status")
    print(f"Submitted job {job_id} ({status}).")
    print(f"Next: python3 shakescenario.py wait {job_id}")


def _print_simple_kv_text(res: dict[str, Any]) -> None:
    for k, v in res.items():
        print(f"{k}: {v}")


def _print_result(obj: dict[str, Any], fmt: str, kind: str) -> None:
    if fmt == "json":
        _print_json(obj)
        return

    if kind == "ping":
        _print_ping_text(obj)
        return
    if kind == "submit":
        _print_submit_text(obj)
        return
    if kind == "list":
        _print_list_text(obj)
        return
    if kind == "models.list":
        _print_models_list_text(obj)
        return
    if kind == "get":
        _print_get_text(obj)
        return

    _print_simple_kv_text(obj)


def _merge_submit_payload(
    args: argparse.Namespace,
) -> dict[str, Any]:
    payload: dict[str, Any] = {}

    req_path = _safe_getattr(args, "request_json", None)
    if req_path:
        loaded = _json_load_file(req_path)
        if not isinstance(loaded, dict):
            raise ValueError("--request-json must contain a JSON object.")
        payload = loaded

    tag = _safe_getattr(args, "tag", None)
    if tag is not None:
        payload["tag"] = tag

    model_id = _safe_getattr(args, "model_id", None)
    if model_id is not None:
        models = payload.get("models", {})
        if not isinstance(models, dict):
            models = {}
        models["model_id"] = model_id
        payload["models"] = models

    scenario = payload.get("scenario", {})
    if not isinstance(scenario, dict):
        scenario = {}

    event = scenario.get("event", {})
    if not isinstance(event, dict):
        event = {}

    hypoc = event.get("hypocentre", {})
    if not isinstance(hypoc, dict):
        hypoc = {}

    mag = _safe_getattr(args, "mag", None)
    lon = _safe_getattr(args, "lon", None)
    lat = _safe_getattr(args, "lat", None)
    depth = _safe_getattr(args, "depth", None)

    if mag is not None:
        event["magnitude"] = float(mag)
    if lon is not None:
        hypoc["longitude"] = float(lon)
    if lat is not None:
        hypoc["latitude"] = float(lat)
    if depth is not None:
        hypoc["elevation"] = -1000.0 * float(depth)

    if hypoc:
        event["hypocentre"] = hypoc
    if event:
        scenario["event"] = event

    gm = scenario.get("ground_motion", {})
    if not isinstance(gm, dict):
        gm = {}

    gmpe = _safe_getattr(args, "gmpe", None)
    dist_approx = _safe_getattr(args, "distance_approx", None)

    if gmpe is not None:
        gm["provider"] = "gmpe"
        gm["gmpe_name"] = str(gmpe)
    if dist_approx is not None:
        gm["distance_approx"] = str(dist_approx)
    if gm:
        scenario["ground_motion"] = gm

    if scenario:
        payload["scenario"] = scenario

    return payload


def main() -> None:
    args = _build_parser().parse_args()
    client = ShakeScenarioClient(host=args.host, port=args.port, timeout=10.0)

    try:
        if args.cmd == "ping":
            res = client.ping()
            _print_result(res, args.format, "ping")
            return

        if args.cmd == "submit":
            payload = _merge_submit_payload(args)
            res = client.submit(payload)
            _print_result(res, args.format, "submit")
            return

        if args.cmd == "list":
            res = client.list_jobs(
                status=args.status,
                limit=args.limit,
                offset=args.offset,
            )
            _print_result(res, args.format, "list")
            return

        if args.cmd == "models":
            if args.models_cmd == "list":
                res = client.models_list()
                _print_result(res, args.format, "models.list")
                return

        if args.cmd == "get":
            res = client.get_job(args.job_id)
            _print_result(res, args.format, "get")
            return

        if args.cmd == "wait":
            job = client.wait_job(
                args.job_id,
                timeout=args.timeout,
                interval=args.interval,
            )
            _print_result(job, args.format, "get")
            return

        if args.cmd == "cancel":
            res = client.cancel_job(args.job_id)
            _print_result(res, args.format, "other")
            return

        if args.cmd == "delete":
            res = client.delete_job(args.job_id, purge=args.purge)
            _print_result(res, args.format, "other")
            return

        if args.cmd == "reset":
            if not args.yes_i_know:
                raise ValueError("Use --yes-i-know to confirm reset.")
            res = client.reset()
            _print_result(res, args.format, "other")
            return

        raise RuntimeError("Unknown command.")

    except Exception as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
