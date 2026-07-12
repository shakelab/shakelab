"""Filesystem repository for ShakeScenario WebUI data."""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any

from .config import WebUIConfig
from .models import (
    EventInfo,
    JsonDict,
    AssetImpact,
    ResultLayer,
    ScenarioRun,
)
from .layers import build_impact_geojson

LOGGER = logging.getLogger(__name__)


class ScenarioRepository:
    """Read ShakeScenario runs, models, and outputs from the filesystem."""

    def __init__(self, config: WebUIConfig):
        self.config = config

    def list_runs(self) -> list[ScenarioRun]:
        """Return summaries of available calculated runs."""
        runs_dir = self.config.runs_dir

        if not runs_dir.exists():
            return []

        runs = []

        for path in sorted(runs_dir.iterdir(), reverse=True):
            if not path.is_dir():
                continue

            if not path.name.startswith("job_"):
                continue

            runs.append(self.read_run(path.name))

        return runs

    def read_run(self, job_id: str) -> ScenarioRun:
        """Read one scenario run summary."""
        run_dir = self._run_dir(job_id)
        manifest = self._read_json(run_dir / "manifest.json", default={})
        request = self._read_json(run_dir / "request.json", default={})

        error_path = run_dir / "logs" / "error.log"
        legacy_error_path = run_dir / "error.txt"
        has_error = self._has_non_empty_file(error_path)
        has_error = has_error or self._has_non_empty_file(legacy_error_path)

        event_data = self._extract_event_data(request, manifest)
        hypocentre = self._first_dict(
            event_data,
            ["hypocentre", "hypocenter", "source", "location"],
        )

        event = EventInfo(
            origin_time=self._first_value(
                event_data,
                ["origin_time", "time", "datetime"],
            ),
            magnitude=self._as_float(self._first_value(
                event_data,
                ["magnitude", "mag", "mw", "Mw", "ML", "ml"],
            )),
            depth_km=self._extract_depth_km(event_data, hypocentre),
            latitude=self._as_float(self._first_value(
                hypocentre or event_data,
                ["latitude", "lat"],
            )),
            longitude=self._as_float(self._first_value(
                hypocentre or event_data,
                ["longitude", "lon", "lng"],
            )),
            place=self._first_value(
                event_data,
                ["place", "location", "description", "name"],
            ),
        )

        status = self._first_value(
            manifest,
            ["status", "state"],
            default="failed" if has_error else "completed",
        )

        model_id = self._first_value(
            manifest,
            ["model_id", "model", "exposure_model"],
        )

        if model_id is None:
            model_id = self._first_value(
                self._first_dict(request, ["models", "model"]),
                ["model_id", "id", "name"],
            )

        return ScenarioRun(
            job_id=job_id,
            status=str(status),
            model_id=model_id,
            event=event,
            created_at=self._first_value(
                manifest,
                ["created_at", "creation_time", "started_at"],
            ),
            completed_at=self._first_value(
                manifest,
                ["completed_at", "end_time", "finished_at"],
            ),
            has_error=has_error,
        )

    def list_models(self) -> list[JsonDict]:
        """Return available input models."""
        models_dir = self.config.models_dir

        if not models_dir.exists():
            return []

        models = []

        for path in sorted(models_dir.iterdir()):
            if path.is_dir():
                models.append({
                    "model_id": path.name,
                    "path": str(path),
                })

        return models

    def list_layers(self, job_id: str) -> list[ResultLayer]:
        """Return available layers for a run."""
        run_dir = self._run_dir(job_id)
        manifest = self._read_json(run_dir / "manifest.json", default={})

        layers = self._layers_from_manifest(manifest)

        if layers:
            return layers

        layers = self._guess_layers(run_dir)

        if layers:
            return layers

        if self._impact_assets_path(job_id).exists():
            return [ResultLayer(
                name="damage",
                layer_type="geojson",
                path=None,
                title="Damage",
            )]

        return []

    def read_layer(self, job_id: str, layer_name: str) -> JsonDict:
        """Read or build one GeoJSON result layer."""
        run_dir = self._run_dir(job_id)

        for layer in self.list_layers(job_id):
            if layer.name != layer_name:
                continue

            if layer.path:
                return self._read_json(run_dir / layer.path, default={})

        if layer_name in ("damage", "impact", "assets"):
            return self._build_damage_layer(job_id)

        return {
            "type": "FeatureCollection",
            "features": [],
        }

    def read_impact_assets_summary(self, job_id: str) -> list[AssetImpact]:
        """Read asset-level impact summary for the impact table."""
        assets = self._read_impact_assets(job_id)
        items = []

        for asset in assets:
            if not isinstance(asset, dict):
                continue

            asset_id = str(asset.get("id", ""))
            asset_name = asset.get("name")
    
            properties = self._damage_properties(asset_id, asset)
    
            items.append(AssetImpact(
                name=(
                    str(asset_name)
                    if asset_name is not None
                    else None
                ),
                asset_id=asset_id,
                damaged_buildings=properties.get("damage"),
                d4_d5=properties.get("damage_d4_d5"),
                pga=properties.get("pga"),
                intensity=None,
                raw=properties,
            ))

        items.sort(
            key=lambda item: self._as_float(item.d4_d5) or 0.0,
            reverse=True,
        )

        return items[:10]

    def read_error(self, job_id: str) -> str | None:
        """Read the error text associated with a failed run."""
        for path in [
            self._run_dir(job_id) / "logs" / "error.log",
            self._run_dir(job_id) / "error.txt",
        ]:
            if path.exists():
                return path.read_text(
                    encoding="utf-8",
                    errors="replace",
                )

        return None

    def delete_run(self, job_id: str) -> bool:
        """Delete one scenario run directory.

        The deletion is restricted to directories inside ``runs_dir`` whose
        name starts with ``job_``. This prevents accidental removal of files
        outside the scenario archive.
        """
        run_dir = self._validated_run_dir(job_id)

        if not run_dir.exists():
            return False

        if not run_dir.is_dir():
            raise ValueError(f"Run path is not a directory: {job_id}")

        shutil.rmtree(run_dir)
        return True

    def _validated_run_dir(self, job_id: str) -> Path:
        """Return a validated run directory path."""
        if not job_id.startswith("job_"):
            raise ValueError(f"Invalid run identifier: {job_id}")

        runs_dir = self.config.runs_dir.resolve()
        run_dir = (runs_dir / job_id).resolve()

        try:
            run_dir.relative_to(runs_dir)
        except ValueError as exc:
            raise ValueError(f"Invalid run path: {job_id}") from exc

        return run_dir

    def _run_dir(self, job_id: str) -> Path:
        """Return the directory associated with a run."""
        return self.config.runs_dir / job_id

    def _impact_assets_path(self, job_id: str) -> Path:
        """Return the canonical impact-assets path for a run."""
        run_dir = self._run_dir(job_id)
        manifest = self._read_json(run_dir / "manifest.json", default={})
        artifacts = manifest.get("artifacts", {})

        if isinstance(artifacts, dict):
            value = artifacts.get("impact_assets")
            if value:
                return run_dir / str(value)

        return run_dir / "results" / "impact_assets.json"

    def _impact_summary_path(self, job_id: str) -> Path:
        """Return the canonical impact-summary path for a run."""
        run_dir = self._run_dir(job_id)
        manifest = self._read_json(run_dir / "manifest.json", default={})
        artifacts = manifest.get("artifacts", {})

        if isinstance(artifacts, dict):
            value = artifacts.get("impact_summary")
            if value:
                return run_dir / str(value)

        return run_dir / "results" / "impact_summary.json"

    def _read_impact_assets(self, job_id: str) -> list[JsonDict]:
        """Read asset-level impact results."""
        impact = self._read_json(
            self._impact_assets_path(job_id),
            default={},
        )

        if isinstance(impact, list):
            return [item for item in impact if isinstance(item, dict)]

        if not isinstance(impact, dict):
            return []

        assets = impact.get("assets", [])

        if not isinstance(assets, list):
            return []

        return [item for item in assets if isinstance(item, dict)]

    def _read_impact_summary(self, job_id: str) -> JsonDict:
        """Read aggregate impact summary."""
        summary = self._read_json(
            self._impact_summary_path(job_id),
            default={},
        )

        if isinstance(summary, dict):
            return summary

        return {}

    def _model_geometry_path(self, job_id: str) -> Path | None:
        """Return the geometry file associated with a run model."""
        run = self.read_run(job_id)

        if not run.model_id:
            return None

        return self.config.models_dir / run.model_id / "geometry.geojson"

    def _read_model_geometry(self, job_id: str) -> JsonDict:
        """Read the geometry GeoJSON associated with a run."""
        path = self._model_geometry_path(job_id)

        if path is None:
            return {
                "type": "FeatureCollection",
                "features": [],
            }

        return self._read_json(path, default={
            "type": "FeatureCollection",
            "features": [],
        })

    def read_impact_assets_object(self, job_id: str) -> JsonDict:
        """Read the canonical impact-assets JSON object."""
        return self._read_impact_assets_object(job_id)

    def _read_impact_assets_object(self, job_id: str) -> JsonDict:
        """Read the canonical impact-assets JSON object."""
        impact = self._read_json(
            self._impact_assets_path(job_id),
            default={},
        )

        if isinstance(impact, dict):
            return impact

        if isinstance(impact, list):
            return {
                "assets": [
                    item for item in impact
                    if isinstance(item, dict)
                ]
            }

        return {
            "assets": [],
        }

    def _build_damage_layer(self, job_id: str) -> JsonDict:
        """Build the damage GeoJSON layer for one run."""
        impact_assets = self._read_impact_assets_object(job_id)
        geometry = self._read_model_geometry(job_id)

        try:
            return build_impact_geojson(
                impact_assets,
                geometry,
                missing_geometry="skip",
            )

        except ValueError as error:
            LOGGER.warning(
                "Unable to build polygon damage layer for %s: %s. "
                "Falling back to point damage layer.",
                job_id,
                error,
            )

        return self._build_point_damage_layer(job_id)

    def _build_point_damage_layer(self, job_id: str) -> JsonDict:
        """Build a point GeoJSON damage layer from asset impacts."""
        features = []

        for asset in self._read_impact_assets(job_id):
            asset_id = asset.get("id")
            location = asset.get("reference_location") or {}

            lon = self._as_float(location.get("longitude"))
            lat = self._as_float(location.get("latitude"))

            if asset_id is None or lon is None or lat is None:
                continue

            features.append({
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [lon, lat],
                },
                "properties": self._damage_properties(
                    str(asset_id),
                    asset,
                ),
            })

        return {
            "type": "FeatureCollection",
            "features": features,
        }

    def _damage_properties(
        self,
        asset_id: str,
        impact: JsonDict,
    ) -> JsonDict:
        """Return flattened damage properties for one asset."""
        damage = impact.get("damage") or {}
        expected = damage.get("expected_counts") or {}
        probabilities = damage.get("probabilities") or {}
        ground_motion = impact.get("ground_motion") or {}
        pga = ground_motion.get("PGA") or {}
        location = impact.get("reference_location") or {}

        d0 = self._as_number(expected.get("D0")) or 0
        d1 = self._as_number(expected.get("D1")) or 0
        d2 = self._as_number(expected.get("D2")) or 0
        d3 = self._as_number(expected.get("D3")) or 0
        d4 = self._as_number(expected.get("D4")) or 0
        d5 = self._as_number(expected.get("D5")) or 0
        gt_last = self._as_number(expected.get("GT_LAST")) or 0

        pga_median = self._as_float(pga.get("median"))
        severe = d4 + d5
        damage_total = severe + gt_last

        return {
            "id": asset_id,
            "asset_id": asset_id,
            "name": impact.get("name"),
            "n_units": self._as_number(impact.get("n_units")),
            "longitude": self._as_float(location.get("longitude")),
            "latitude": self._as_float(location.get("latitude")),
            "pga": pga_median,
            "PGA": pga_median,
            "damage": damage_total,
            "damage_d4_d5": severe,
            "damage_gt_last": gt_last,
            "D0": d0,
            "D1": d1,
            "D2": d2,
            "D3": d3,
            "D4": d4,
            "D5": d5,
            "GT_LAST": gt_last,
            "prob_D0": self._as_float(probabilities.get("D0")),
            "prob_D1": self._as_float(probabilities.get("D1")),
            "prob_D2": self._as_float(probabilities.get("D2")),
            "prob_D3": self._as_float(probabilities.get("D3")),
            "prob_D4": self._as_float(probabilities.get("D4")),
            "prob_D5": self._as_float(probabilities.get("D5")),
            "prob_GT_LAST": self._as_float(
                probabilities.get("GT_LAST"),
            ),
        }

    @staticmethod
    def _extract_event_data(
        request: JsonDict,
        manifest: JsonDict,
    ) -> JsonDict:
        """Extract event information from request or manifest."""
        scenario = ScenarioRepository._first_dict(request, ["scenario"])
        event = ScenarioRepository._first_dict(scenario, ["event"])

        if event:
            return event

        return ScenarioRepository._first_dict(
            manifest,
            ["event", "earthquake", "source"],
        )

    @staticmethod
    def _extract_depth_km(
        event_data: JsonDict,
        hypocentre: JsonDict,
    ) -> float | None:
        """Extract hypocentral depth in kilometres."""
        depth = ScenarioRepository._first_value(
            event_data,
            ["depth_km", "depth", "depth_in_km"],
        )

        if depth is not None:
            return ScenarioRepository._as_float(depth)

        elevation = ScenarioRepository._first_value(
            hypocentre,
            ["elevation", "elev", "z"],
        )
        elevation = ScenarioRepository._as_float(elevation)

        if elevation is None:
            return None

        if elevation < 0:
            return abs(elevation) / 1000.0

        return None

    @staticmethod
    def _layers_from_manifest(manifest: JsonDict) -> list[ResultLayer]:
        """Extract declared GeoJSON layers from manifest artifacts."""
        artifacts = manifest.get("artifacts", {})

        if not isinstance(artifacts, dict):
            return []

        layers_data = artifacts.get("layers")

        if not isinstance(layers_data, list):
            return []

        layers = []

        for item in layers_data:
            if not isinstance(item, dict):
                continue

            name = item.get("name")
            path = item.get("path")

            if not name or not path:
                continue

            layers.append(ResultLayer(
                name=str(name),
                layer_type=str(item.get("type", "geojson")),
                path=str(path),
                title=item.get("title"),
            ))

        return layers

    @staticmethod
    def _guess_layers(run_dir: Path) -> list[ResultLayer]:
        """Guess available layers if no manifest declaration exists."""
        layers = []
        candidate_dirs = [
            run_dir / "results" / "layers",
            run_dir / "layers",
            run_dir,
        ]

        for directory in candidate_dirs:
            if not directory.exists():
                continue

            for path in sorted(directory.glob("*.geojson")):
                layers.append(ResultLayer(
                    name=path.stem,
                    layer_type="geojson",
                    path=str(path.relative_to(run_dir)),
                    title=path.stem.replace("_", " ").title(),
                ))

        return layers

    @staticmethod
    def _read_json(path: Path, default: Any = None) -> Any:
        """Read a JSON file, returning a default value if absent."""
        if not path.exists():
            return default

        with path.open("r", encoding="utf-8") as file:
            return json.load(file)

    @staticmethod
    def _has_non_empty_file(path: Path) -> bool:
        """Return true if path exists and contains data."""
        return path.exists() and path.stat().st_size > 0

    @staticmethod
    def _first_value(
        data: JsonDict,
        keys: list[str],
        default: Any = None,
    ) -> Any:
        """Return the first available value from a dictionary."""
        if not isinstance(data, dict):
            return default

        for key in keys:
            if key in data and data[key] is not None:
                return data[key]

        return default

    @staticmethod
    def _first_dict(data: JsonDict, keys: list[str]) -> JsonDict:
        """Return the first dictionary value from a dictionary."""
        value = ScenarioRepository._first_value(data, keys, default={})

        if isinstance(value, dict):
            return value

        return {}

    @staticmethod
    def _as_float(value: Any) -> float | None:
        """Convert a value to float when possible."""
        if value is None or value == "":
            return None

        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _as_number(value: Any) -> float | int | None:
        """Convert a value to int or float when possible."""
        if value is None or value == "":
            return None

        try:
            number = float(value)
        except (TypeError, ValueError):
            return None

        if number.is_integer():
            return int(number)

        return number
