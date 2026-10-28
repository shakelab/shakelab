"""HTTP regression tests for the ShakeScenario WebUI API."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from shakelab.engineering.shakescenario.webui.app import create_app
from shakelab.engineering.shakescenario.webui.backend.config import (
    WebUIConfig,
)


@pytest.fixture
def client(tmp_path: Path):
    """Create a Flask client backed by an isolated temporary archive."""
    (tmp_path / "runs").mkdir()
    config = WebUIConfig(
        home=tmp_path,
        server_host="127.0.0.1",
        server_port=6000,
        server_timeout=1.0,
        web_host="127.0.0.1",
        web_port=8000,
        debug=False,
    )
    app = create_app(config)
    app.config["TESTING"] = True
    return app.test_client()


def _make_run(home: Path, job_id: str, manifest: dict) -> None:
    run_dir = home / "runs" / job_id
    run_dir.mkdir()
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8",
    )


def test_failed_run_error_in_detail(client, tmp_path):
    """The detail endpoint exposes the structured failure message."""
    error = {"type": "ValueError", "message": "Missing rake"}
    _make_run(
        tmp_path, "job_000001",
        {"status": "failed", "error": error},
    )
    response = client.get("/api/v1/runs/job_000001")
    assert response.status_code == 200
    data = response.get_json()
    assert data["status"] == "failed"
    assert data["error"] == error


def test_failed_run_error_in_archive(client, tmp_path):
    """The archive endpoint preserves structured diagnostics."""
    _make_run(
        tmp_path, "job_000002",
        {
            "status": "failed",
            "error": {"type": "RuntimeError", "message": "GMPE failed"},
        },
    )
    response = client.get("/api/v1/runs")
    assert response.status_code == 200
    runs = response.get_json()["runs"]
    assert len(runs) == 1
    assert runs[0]["error"]["message"] == "GMPE failed"


def test_legacy_run_detail_has_null_error(client, tmp_path):
    """Historical runs without structured diagnostics remain readable."""
    _make_run(tmp_path, "job_000003", {"status": "failed"})
    response = client.get("/api/v1/runs/job_000003")
    assert response.status_code == 200
    assert response.get_json()["error"] is None


def test_import_and_factory_do_not_require_runtime_config(tmp_path):
    """Explicit config allows isolated application construction."""
    config = WebUIConfig(
        home=tmp_path,
        server_host="127.0.0.1",
        server_port=6000,
        server_timeout=1.0,
        web_host="127.0.0.1",
        web_port=8000,
        debug=False,
    )
    (tmp_path / "runs").mkdir()
    app = create_app(config)
    assert app.config["WEBUI_CONFIG"] == config
