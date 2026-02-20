# database.py
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
ShakeLab - ShakeScenario sqlite3 persistence.

This module provides a thin sqlite wrapper for job persistence.

Schema (v1)
-----------
jobs(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  created_at REAL NOT NULL,
  started_at REAL,
  ended_at REAL,
  status TEXT NOT NULL,
  tag TEXT,
  params_json TEXT NOT NULL,
  result_meta_json TEXT,
  error TEXT
)

Notes
-----
- sqlite3 is used from multiple threads; we keep a single connection
  (check_same_thread=False) and serialize access with a lock.

"""

from __future__ import annotations

import json
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any

from models import JobStatus


class JobDatabase:
    """sqlite-backed job database."""

    def __init__(self, path: Path) -> None:
        self._path = Path(path)
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(
            str(self._path),
            check_same_thread=False,
        )
        self._conn.row_factory = sqlite3.Row

    def initialize(self) -> None:
        """Create tables if needed."""
        with self._lock, self._conn:
            self._conn.execute(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  created_at REAL NOT NULL,
                  started_at REAL,
                  ended_at REAL,
                  status TEXT NOT NULL,
                  tag TEXT,
                  params_json TEXT NOT NULL,
                  result_meta_json TEXT,
                  error TEXT
                )
                """
            )
            self._conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status)")
            self._conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_created ON jobs(created_at)")

    def reset(self) -> None:
        """Delete all jobs."""
        with self._lock, self._conn:
            self._conn.execute("DELETE FROM jobs")

    def create_job(self, params: dict[str, Any], tag: str | None = None) -> int:
        """Insert a new job and return its id."""
        now = time.time()
        params_json = json.dumps(params, ensure_ascii=False)

        with self._lock, self._conn:
            cur = self._conn.execute(
                """
                INSERT INTO jobs (
                  created_at, status, tag, params_json
                ) VALUES (?, ?, ?, ?)
                """,
                (now, JobStatus.QUEUED.value, tag, params_json),
            )
            return int(cur.lastrowid)

    def update_status(self, job_id: int, status: JobStatus) -> None:
        """Update job status and timestamps."""
        now = time.time()

        started_at = None
        ended_at = None
        if status == JobStatus.RUNNING:
            started_at = now
        if status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELED):
            ended_at = now

        with self._lock, self._conn:
            row = self._conn.execute(
                "SELECT started_at FROM jobs WHERE id=?",
                (job_id,),
            ).fetchone()
            if row is None:
                return

            old_started = row["started_at"]
            if old_started is not None:
                started_at = old_started

            self._conn.execute(
                """
                UPDATE jobs
                SET status=?,
                    started_at=COALESCE(?, started_at),
                    ended_at=COALESCE(?, ended_at)
                WHERE id=?
                """,
                (status.value, started_at, ended_at, job_id),
            )

    def update_result_meta(self, job_id: int, meta: dict[str, Any]) -> None:
        """Attach result metadata (e.g., workdir, manifest)."""
        meta_json = json.dumps(meta, ensure_ascii=False)
        with self._lock, self._conn:
            self._conn.execute(
                "UPDATE jobs SET result_meta_json=? WHERE id=?",
                (meta_json, job_id),
            )

    def update_error(self, job_id: int, error: str) -> None:
        """Attach an error message."""
        with self._lock, self._conn:
            self._conn.execute(
                "UPDATE jobs SET error=? WHERE id=?",
                (error, job_id),
            )

    def get_job(self, job_id: int) -> dict[str, Any] | None:
        """Fetch a single job as a dict."""
        with self._lock, self._conn:
            row = self._conn.execute(
                "SELECT * FROM jobs WHERE id=?",
                (job_id,),
            ).fetchone()

        if row is None:
            return None
        return self._row_to_dict(row)

    def list_jobs(
        self,
        status: str | None = None,
        limit: int = 50,
        offset: int = 0,
    ) -> list[dict[str, Any]]:
        """List jobs with optional status filter."""
        limit = max(1, min(int(limit), 500))
        offset = max(0, int(offset))

        if status is None:
            q = "SELECT * FROM jobs ORDER BY id DESC LIMIT ? OFFSET ?"
            args = (limit, offset)
        else:
            q = """
            SELECT * FROM jobs
            WHERE status=?
            ORDER BY id DESC
            LIMIT ? OFFSET ?
            """
            args = (status, limit, offset)

        with self._lock, self._conn:
            rows = self._conn.execute(q, args).fetchall()

        return [self._row_to_dict(r) for r in rows]

    def delete_job(self, job_id: int) -> bool:
        """Delete a job and return True if it existed."""
        with self._lock, self._conn:
            cur = self._conn.execute("DELETE FROM jobs WHERE id=?", (job_id,))
            return cur.rowcount > 0

    @staticmethod
    def _row_to_dict(row: sqlite3.Row) -> dict[str, Any]:
        out: dict[str, Any] = dict(row)

        try:
            out["params"] = json.loads(out.pop("params_json") or "{}")
        except json.JSONDecodeError:
            out["params"] = {}

        meta_json = out.pop("result_meta_json", None)
        if meta_json:
            try:
                out["result_meta"] = json.loads(meta_json)
            except json.JSONDecodeError:
                out["result_meta"] = {}
        else:
            out["result_meta"] = {}

        return out

