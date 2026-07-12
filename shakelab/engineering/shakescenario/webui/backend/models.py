"""Data models for the ShakeScenario WebUI backend."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


JsonDict = dict[str, Any]


@dataclass(frozen=True)
class EventInfo:
    """Earthquake source information associated with a scenario run."""

    origin_time: str | None = None
    magnitude: float | None = None
    depth_km: float | None = None
    latitude: float | None = None
    longitude: float | None = None
    place: str | None = None

    def to_dict(self) -> JsonDict:
        """Return a JSON-serializable representation."""
        return {
            "origin_time": self.origin_time,
            "magnitude": self.magnitude,
            "depth_km": self.depth_km,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "place": self.place,
        }


@dataclass(frozen=True)
class ScenarioRun:
    """Summary information for a calculated ShakeScenario run."""

    job_id: str
    status: str = "unknown"
    model_id: str | None = None
    event: EventInfo = field(default_factory=EventInfo)
    created_at: str | None = None
    completed_at: str | None = None
    has_error: bool = False

    def to_dict(self) -> JsonDict:
        """Return a JSON-serializable representation."""
        return {
            "job_id": self.job_id,
            "status": self.status,
            "model_id": self.model_id,
            "event": self.event.to_dict(),
            "created_at": self.created_at,
            "completed_at": self.completed_at,
            "has_error": self.has_error,
        }


@dataclass(frozen=True)
class ResultLayer:
    """Description of an available result layer."""

    name: str
    layer_type: str = "geojson"
    path: str | None = None
    title: str | None = None

    def to_dict(self) -> JsonDict:
        """Return a JSON-serializable representation."""
        return {
            "name": self.name,
            "type": self.layer_type,
            "path": self.path,
            "title": self.title,
        }


@dataclass(frozen=True)
class AssetImpact:
    """Asset-level impact summary."""

    name: str | None = None
    asset_id: str | None = None
    damaged_buildings: float | int | None = None
    d4_d5: float | int | None = None
    pga: float | None = None
    intensity: str | float | None = None
    raw: JsonDict = field(default_factory=dict)

    def to_dict(self) -> JsonDict:
        """Return a JSON-serializable representation."""
        data = dict(self.raw)

        data.update({
            "name": self.name,
            "asset_id": self.asset_id,
            "damaged_buildings": self.damaged_buildings,
            "d4_d5": self.d4_d5,
            "pga": self.pga,
            "intensity": self.intensity,
        })

        return data