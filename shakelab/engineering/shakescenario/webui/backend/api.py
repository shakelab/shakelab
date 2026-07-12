"""Flask API routes for the ShakeScenario WebUI."""

from __future__ import annotations

import json

from flask import Blueprint, jsonify, request, Response

from .client import Client
from .config import WebUIConfig
from .repository import ScenarioRepository


def create_api_blueprint(config: WebUIConfig) -> Blueprint:
    """Create and configure the WebUI API blueprint."""
    api = Blueprint("api", __name__, url_prefix="/api/v1")
    repository = ScenarioRepository(config)
    client = Client(config)

    # ------------------------------------------------------------------
    # Backend and server status
    # ------------------------------------------------------------------

    @api.get("/status")
    def status():
        """Return WebUI backend status."""
        return jsonify({
            "status": "online",
            "home": str(config.home),
            "runs_dir": str(config.runs_dir),
            "models_dir": str(config.models_dir),
            "database": str(config.database_path),
            "runs_dir_exists": config.runs_dir.exists(),
            "models_dir_exists": config.models_dir.exists(),
        })

    @api.get("/server")
    def server():
        """Return ShakeScenario server status."""
        server_info = {
            "host": config.server_host,
            "port": config.server_port,
            "timeout": config.server_timeout,
        }

        try:
            capabilities = client.ping()

            return jsonify({
                "status": "online",
                "server": server_info,
                "capabilities": capabilities,
            })

        except Exception as exc:
            return jsonify({
                "status": "offline",
                "server": server_info,
                "error": str(exc),
            }), 503

    @api.get("/models")
    def models():
        """Return registered ShakeScenario models."""
        try:
            return jsonify(client.list_models())

        except Exception as exc:
            return jsonify({
                "error": "server_unavailable",
                "message": str(exc),
            }), 503

    # ------------------------------------------------------------------
    # Scenario submission
    # ------------------------------------------------------------------

    @api.post("/scenarios")
    def create_scenario():
        """Submit a new scenario to the ShakeScenario server."""
        data = request.get_json(silent=True) or {}

        required = [
            "tag",
            "model_id",
            "origin_time",
            "magnitude",
            "longitude",
            "latitude",
            "depth",
        ]

        missing = [
            name for name in required
            if data.get(name) in (None, "")
        ]

        if missing:
            return jsonify({
                "error": "invalid_request",
                "message": "Missing required scenario parameters.",
                "missing": missing,
            }), 400

        try:
            response = client.submit(
                tag=data.get("tag"),
                model_id=data["model_id"],
                origin_time=data["origin_time"],
                magnitude=data["magnitude"],
                longitude=data["longitude"],
                latitude=data["latitude"],
                depth=data["depth"],
                gmpe=data.get("gmpe"),
                distance_approx=data.get("distance_approx"),
            )

            return jsonify(response), 202

        except Exception as exc:
            return jsonify({
                "error": "submission_failed",
                "message": str(exc),
            }), 502

    # ------------------------------------------------------------------
    # Scenario archive
    # ------------------------------------------------------------------

    @api.get("/runs")
    def runs():
        """Return available scenario runs."""
        items = [
            run.to_dict()
            for run in repository.list_runs()
        ]

        return jsonify({
            "runs": items,
        })

    @api.get("/runs/<job_id>")
    def run_detail(job_id: str):
        """Return details for one scenario run."""
        run = repository.read_run(job_id)

        return jsonify(run.to_dict())

    @api.delete("/runs/<job_id>")
    def delete_run(job_id: str):
        """Delete one scenario run."""
        try:
            deleted = repository.delete_run(job_id)

        except ValueError as exc:
            return jsonify({
                "error": "invalid_run",
                "message": str(exc),
            }), 400

        if not deleted:
            return jsonify({
                "error": "not_found",
                "message": f"Scenario run not found: {job_id}",
            }), 404

        return jsonify({
            "job_id": job_id,
            "deleted": True,
        })

    @api.get("/runs/<job_id>/error")
    def run_error(job_id: str):
        """Return error text for a failed run, if available."""
        error_text = repository.read_error(job_id)

        return jsonify({
            "job_id": job_id,
            "error": error_text,
        })

    # ------------------------------------------------------------------
    # Layers and impact data
    # ------------------------------------------------------------------

    @api.get("/runs/<job_id>/layers")
    def run_layers(job_id: str):
        """Return available layers for one run."""
        layers = [
            layer.to_dict()
            for layer in repository.list_layers(job_id)
        ]

        return jsonify({
            "layers": layers,
        })

    @api.get("/runs/<job_id>/layer/<layer_name>")
    def run_layer(job_id: str, layer_name: str):
        """Return one GeoJSON layer for one run."""
        return jsonify(repository.read_layer(job_id, layer_name))

    @api.get("/runs/<job_id>/impact-assets")
    def run_impact_assets(job_id: str):
        """Return asset-level impact summary for one run."""
        items = [
            item.to_dict()
            for item in repository.read_impact_assets_summary(job_id)
        ]

        return jsonify({
            "assets": items,
        })

    # ------------------------------------------------------------------
    # Downloads
    # ------------------------------------------------------------------

    @api.get("/runs/<job_id>/download/json")
    def download_run_json(job_id: str):
        """Download canonical impact assets JSON for one run."""
        data = repository.read_impact_assets_object(job_id)

        return Response(
            json.dumps(data, indent=2, ensure_ascii=False),
            mimetype="application/json",
            headers={
                "Content-Disposition": (
                    f"attachment; filename={job_id}_impact_assets.json"
                )
            },
        )

    @api.get("/runs/<job_id>/download/geojson")
    def download_run_geojson(job_id: str):
        """Download generated GeoJSON damage layer for one run."""
        data = repository.read_layer(job_id, "damage")

        return Response(
            json.dumps(data, indent=2, ensure_ascii=False),
            mimetype="application/geo+json",
            headers={
                "Content-Disposition": (
                    f"attachment; filename={job_id}_damage.geojson"
                )
            },
        )

    return api