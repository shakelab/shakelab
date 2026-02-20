# shakescenario-server.py
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
ShakeLab - ShakeScenario service (TCP server).

This module provides a small TCP server that exposes a job-based API for
running rapid damage/impact scenarios.

Protocol
--------
- Length-prefixed JSON messages over TCP.
- Each request/response is a JSON object with:
  - v: int (API version)
  - op: str (operation)
  - req_id: str (client request id)
  - payload: dict (operation parameters)

Persistence
-----------
Jobs are persisted in a sqlite3 database. Job artifacts are stored on the
filesystem in a per-job directory under the configured workdir.

Notes
-----
- This is intentionally "simple layout" (flat files) to fit inside
  shakelab.engineering, while keeping the service self-contained.

"""

from __future__ import annotations

import argparse
import json
import socket
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path
from typing import Any

from database import JobDatabase
from models import JobStatus, ServerInfo
from protocol import ProtocolError, recv_message, send_message


def _safe_json(obj: Any) -> Any:
    """
    Return an object suitable for JSON serialization.

    This is a defensive helper to avoid crashes when writing "result" fields.
    """
    try:
        json.dumps(obj)
        return obj
    except TypeError:
        return str(obj)


class ShakeScenarioServer:
    """
    TCP server for ShakeScenario jobs.

    Parameters
    ----------
    host, port
        Bind address and port.
    db
        Job database wrapper (sqlite3).
    workdir
        Root directory for job artifacts.
    workers
        Number of scenario workers executing jobs in parallel.

    """

    api_version = 1
    server_name = "shakescenario-server"

    def __init__(
        self,
        host: str,
        port: int,
        db: JobDatabase,
        workdir: Path,
        workers: int = 2,
    ) -> None:
        self._host = host
        self._port = port
        self._db = db
        self._workdir = workdir
        self._workdir.mkdir(parents=True, exist_ok=True)

        self._sock: socket.socket | None = None
        self._stop = threading.Event()

        self._executor = ThreadPoolExecutor(max_workers=max(1, workers))

    def serve_forever(self) -> None:
        """Start accepting connections."""
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((self._host, self._port))
        sock.listen(64)
        self._sock = sock

        try:
            while not self._stop.is_set():
                conn, addr = sock.accept()
                t = threading.Thread(
                    target=self._handle_connection,
                    args=(conn, addr),
                    daemon=True,
                )
                t.start()
        finally:
            try:
                sock.close()
            except OSError:
                pass
            self._executor.shutdown(wait=False, cancel_futures=True)

    def shutdown(self) -> None:
        """Request server shutdown."""
        self._stop.set()
        if self._sock is not None:
            try:
                self._sock.close()
            except OSError:
                pass

    def _handle_connection(self, conn: socket.socket, addr: Any) -> None:
        """Handle a single TCP client connection."""
        try:
            with conn:
                while True:
                    try:
                        req = recv_message(conn)
                    except EOFError:
                        return
                    except ProtocolError as exc:
                        self._send_error(
                            conn,
                            req_id=None,
                            code="BAD_PROTOCOL",
                            message=str(exc),
                        )
                        return

                    req_id = req.get("req_id")
                    op = req.get("op")
                    v = req.get("v")
                    payload = req.get("payload", {})

                    if v != self.api_version:
                        self._send_error(
                            conn,
                            req_id=req_id,
                            code="BAD_VERSION",
                            message=f"Unsupported API version: {v}",
                        )
                        continue

                    if not isinstance(op, str):
                        self._send_error(
                            conn,
                            req_id=req_id,
                            code="INVALID_REQUEST",
                            message="Missing or invalid 'op'.",
                        )
                        continue

                    try:
                        res = self._dispatch(op, payload)
                    except Exception as exc:  # noqa: BLE001
                        self._send_error(
                            conn,
                            req_id=req_id,
                            code="SERVER_ERROR",
                            message=str(exc),
                        )
                        continue

                    send_message(
                        conn,
                        {
                            "v": self.api_version,
                            "req_id": req_id,
                            "ok": True,
                            "result": _safe_json(res),
                        },
                    )
        except Exception:  # noqa: BLE001
            # Intentionally silent: avoid crashing the server thread.
            return

    def _send_error(
        self,
        conn: socket.socket,
        req_id: str | None,
        code: str,
        message: str,
    ) -> None:
        send_message(
            conn,
            {
                "v": self.api_version,
                "req_id": req_id,
                "ok": False,
                "error": {"code": code, "message": message},
            },
        )

    def _dispatch(self, op: str, payload: dict[str, Any]) -> dict[str, Any]:
        """Dispatch an operation."""
        if op == "ping":
            return self._op_ping()

        if op == "submit":
            return self._op_submit(payload)

        if op == "list":
            return self._op_list(payload)

        if op == "get":
            return self._op_get(payload)

        if op == "cancel":
            return self._op_cancel(payload)

        if op == "delete":
            return self._op_delete(payload)

        if op == "reset":
            return self._op_reset(payload)

        raise ValueError(f"Unknown operation: {op}")

    def _op_ping(self) -> dict[str, Any]:
        info = ServerInfo(
            name=self.server_name,
            api_version=self.api_version,
            protocol="json-tcp-length-prefixed",
        )
        return asdict(info)

    def _op_submit(self, payload: dict[str, Any]) -> dict[str, Any]:
        """
        Submit a new scenario job.

        Expected payload keys (v1 baseline)
        ----------------------------------
        - params: dict
            Scenario parameters (event, config, options, etc.).
        - tag: str, optional
            User-defined label for quick filtering.
        """
        params = payload.get("params", {})
        tag = payload.get("tag")

        if not isinstance(params, dict):
            raise ValueError("'params' must be a dict.")
        if tag is not None and not isinstance(tag, str):
            raise ValueError("'tag' must be a string or omitted.")

        job_id = self._db.create_job(params=params, tag=tag)
        job_dir = self._job_dir(job_id)
        job_dir.mkdir(parents=True, exist_ok=True)

        self._executor.submit(self._run_job, job_id, params, job_dir)

        return {"job_id": job_id, "status": JobStatus.QUEUED.value}

    def _op_list(self, payload: dict[str, Any]) -> dict[str, Any]:
        status = payload.get("status")
        limit = payload.get("limit", 50)
        offset = payload.get("offset", 0)

        if status is not None and not isinstance(status, str):
            raise ValueError("'status' must be a string or omitted.")

        rows = self._db.list_jobs(
            status=status,
            limit=int(limit),
            offset=int(offset),
        )
        return {"jobs": rows}

    def _op_get(self, payload: dict[str, Any]) -> dict[str, Any]:
        job_id = payload.get("job_id")
        if job_id is None:
            raise ValueError("Missing 'job_id'.")
        job = self._db.get_job(int(job_id))
        if job is None:
            raise ValueError(f"Job not found: {job_id}")
        return job

    def _op_cancel(self, payload: dict[str, Any]) -> dict[str, Any]:
        job_id = payload.get("job_id")
        if job_id is None:
            raise ValueError("Missing 'job_id'.")

        job = self._db.get_job(int(job_id))
        if job is None:
            raise ValueError(f"Job not found: {job_id}")

        status = job["status"]
        if status in (JobStatus.COMPLETED.value, JobStatus.FAILED.value):
            return {"job_id": int(job_id), "status": status}

        # v1: best-effort cancellation. If queued -> mark canceled.
        if status == JobStatus.QUEUED.value:
            self._db.update_status(int(job_id), JobStatus.CANCELED)
            return {"job_id": int(job_id), "status": JobStatus.CANCELED.value}

        # If running, we do not force-kill in v1.
        return {"job_id": int(job_id), "status": status}

    def _op_delete(self, payload: dict[str, Any]) -> dict[str, Any]:
        job_id = payload.get("job_id")
        purge = bool(payload.get("purge", False))

        if job_id is None:
            raise ValueError("Missing 'job_id'.")

        ok = self._db.delete_job(int(job_id))
        if not ok:
            raise ValueError(f"Job not found: {job_id}")

        if purge:
            job_dir = self._job_dir(int(job_id))
            if job_dir.exists():
                for p in sorted(job_dir.rglob("*"), reverse=True):
                    if p.is_file():
                        p.unlink(missing_ok=True)
                    else:
                        p.rmdir()
                job_dir.rmdir()

        return {"job_id": int(job_id), "deleted": True, "purge": purge}

    def _op_reset(self, payload: dict[str, Any]) -> dict[str, Any]:
        confirm = bool(payload.get("confirm", False))
        if not confirm:
            raise ValueError("Reset requires confirm=true.")
        self._db.reset()
        return {"reset": True}

    def _job_dir(self, job_id: int) -> Path:
        return self._workdir / f"job_{job_id:06d}"

    def _run_job(
        self,
        job_id: int,
        params: dict[str, Any],
        job_dir: Path,
    ) -> None:
        """
        Execute a job.

        This is the integration point where ShakeScenario will call the
        real pipeline (ground motion + fragility + impact).

        For v1 bootstrap, we write metadata and a minimal manifest.
        """
        # If someone canceled while queued, do not start.
        job = self._db.get_job(job_id)
        if job is None:
            return
        if job["status"] == JobStatus.CANCELED.value:
            return

        self._db.update_status(job_id, JobStatus.RUNNING)

        try:
            # --- Placeholder computation ---
            # Replace with actual ShakeScenario compute call.
            meta = {
                "job_id": job_id,
                "params": params,
                "note": "Placeholder job runner; replace with real pipeline.",
            }
            (job_dir / "meta.json").write_text(
                json.dumps(meta, indent=2, sort_keys=True),
                encoding="utf-8",
            )

            manifest = {
                "job_id": job_id,
                "files": ["meta.json"],
            }
            (job_dir / "result_manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True),
                encoding="utf-8",
            )
            # --- End placeholder computation ---

            self._db.update_status(job_id, JobStatus.COMPLETED)
            self._db.update_result_meta(
                job_id,
                {"workdir": str(job_dir), "manifest": "result_manifest.json"},
            )

        except Exception as exc:  # noqa: BLE001
            self._db.update_status(job_id, JobStatus.FAILED)
            self._db.update_error(job_id, str(exc))


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="shakescenario-server",
        description="ShakeScenario TCP server (job-based).",
    )
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=6000)
    p.add_argument("--db", default="shakescenario.db")
    p.add_argument("--workdir", default="./runs")
    p.add_argument("--workers", type=int, default=2)
    return p


def main() -> None:
    args = _build_arg_parser().parse_args()

    db_path = Path(args.db).expanduser().resolve()
    workdir = Path(args.workdir).expanduser().resolve()

    db = JobDatabase(db_path)
    db.initialize()

    srv = ShakeScenarioServer(
        host=args.host,
        port=args.port,
        db=db,
        workdir=workdir,
        workers=args.workers,
    )
    srv.serve_forever()


if __name__ == "__main__":
    main()
