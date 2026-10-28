"""Regression tests for WebUI scenario error propagation."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from shakelab.engineering.shakescenario.webui.backend.config import (
    WebUIConfig,
)
from shakelab.engineering.shakescenario.webui.backend.repository import (
    ScenarioRepository,
)


class ScenarioErrorTests(unittest.TestCase):
    """Check error propagation and compatibility with historical runs."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        home = Path(self.directory.name)
        config = WebUIConfig(
            home=home,
            server_host="127.0.0.1",
            server_port=5001,
            server_timeout=10.0,
            web_host="127.0.0.1",
            web_port=8000,
            debug=False,
        )
        self.repository = ScenarioRepository(config)
        config.runs_dir.mkdir()

    def make_run(self, name, manifest, with_log=False):
        run_dir = self.repository.config.runs_dir / name
        run_dir.mkdir()
        (run_dir / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8",
        )
        if with_log:
            (run_dir / "logs").mkdir()
            (run_dir / "logs" / "error.log").write_text(
                "Traceback details", encoding="utf-8",
            )

    def test_failed_run_exposes_manifest_error(self):
        error = {
            "type": "ValueError",
            "message": "Missing required GMPE parameter: vs30",
        }
        self.make_run(
            "job_000001",
            {"status": "failed", "error": error},
            with_log=True,
        )
        run = self.repository.read_run("job_000001").to_dict()
        self.assertEqual(run["status"], "failed")
        self.assertEqual(run["error"], error)
        self.assertTrue(run["has_error"])
        self.assertEqual(
            json.loads(json.dumps(run))["error"]["message"],
            error["message"],
        )

    def test_legacy_run_without_error_field(self):
        self.make_run("job_000002", {"status": "failed"}, with_log=True)
        run = self.repository.read_run("job_000002").to_dict()
        self.assertEqual(run["status"], "failed")
        self.assertIsNone(run["error"])
        self.assertTrue(run["has_error"])

    def test_completed_run_has_no_error(self):
        self.make_run("job_000003", {"status": "completed"})
        run = self.repository.read_run("job_000003").to_dict()
        self.assertEqual(run["status"], "completed")
        self.assertIsNone(run["error"])
        self.assertFalse(run["has_error"])

    def test_invalid_error_payload_is_ignored(self):
        self.make_run(
            "job_000004",
            {"status": "failed", "error": "invalid legacy format"},
        )
        run = self.repository.read_run("job_000004").to_dict()
        self.assertIsNone(run["error"])

    def test_list_runs_preserves_error_details(self):
        self.make_run(
            "job_000005",
            {
                "status": "failed",
                "error": {"type": "ValueError", "message": "Missing rake"},
            },
        )
        runs = [run.to_dict() for run in self.repository.list_runs()]
        self.assertEqual(len(runs), 1)
        self.assertEqual(runs[0]["error"]["message"], "Missing rake")


if __name__ == "__main__":
    unittest.main()
