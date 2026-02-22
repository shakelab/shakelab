# models.py
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
ShakeLab - ShakeScenario data models.

This module defines small, stable datatypes shared by server/client/database.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any
import json
import re


class JobStatus(str, Enum):
    """Job status values (persisted)."""

    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELED = "canceled"


@dataclass(frozen=True)
class ServerInfo:
    """Server info returned by ping."""

    name: str
    api_version: int
    protocol: str


# ---------------------------------------------------------------------------
# Server configuration + model resolution
# ---------------------------------------------------------------------------

_MODEL_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")


@dataclass(frozen=True)
class ServerPaths:
    """Server path configuration."""

    db: Path
    workdir: Path
    model_root: Path


@dataclass(frozen=True)
class GroundMotionDefaults:
    """Default ground-motion provider configuration."""

    provider: str
    gmpe_name: str | None = None
    distance_approx: str | None = None


@dataclass(frozen=True)
class ServerDefaults:
    """Default values applied when submit payload omits them."""

    model_id: str
    ground_motion: GroundMotionDefaults
    impact_config: dict[str, Any]


@dataclass(frozen=True)
class ServerConfig:
    """Resolved server configuration."""

    schema_version: str
    paths: ServerPaths
    server: dict[str, Any]
    defaults: ServerDefaults


def load_server_config(path: str | Path) -> ServerConfig:
    """
    Load and validate the server configuration from a JSON file.

    Parameters
    ----------
    path
        Path to the JSON configuration file.

    Returns
    -------
    ServerConfig
        Parsed and validated configuration.

    Raises
    ------
    ValueError
        If the configuration is invalid.

    Notes
    -----
    The JSON file is intentionally strict (no comments). Validation happens
    at server startup to fail fast with explicit error messages.
    """
    cfg_path = Path(path).expanduser().resolve()
    try:
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(f"Cannot read config: {cfg_path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON config: {cfg_path}") from exc

    if not isinstance(raw, dict):
        raise ValueError("Config must be a JSON object.")

    schema_version = str(raw.get("schema_version", "")).strip()
    if not schema_version:
        raise ValueError("Missing 'schema_version' in config.")

    paths = _parse_paths(raw.get("paths"))
    server = _parse_server(raw.get("server"))
    defaults = _parse_defaults(raw.get("defaults"))

    _validate_paths(paths)
    _validate_defaults(paths, defaults)

    return ServerConfig(
        schema_version=schema_version,
        paths=paths,
        server=server,
        defaults=defaults,
    )


def _parse_paths(obj: Any) -> ServerPaths:
    if not isinstance(obj, dict):
        raise ValueError("Missing or invalid 'paths' section in config.")

    db = obj.get("db")
    workdir = obj.get("workdir")
    model_root = obj.get("model_root")

    if not isinstance(db, str) or not db:
        raise ValueError("paths.db must be a non-empty string.")
    if not isinstance(workdir, str) or not workdir:
        raise ValueError("paths.workdir must be a non-empty string.")
    if not isinstance(model_root, str) or not model_root:
        raise ValueError("paths.model_root must be a non-empty string.")

    return ServerPaths(
        db=Path(db).expanduser().resolve(),
        workdir=Path(workdir).expanduser().resolve(),
        model_root=Path(model_root).expanduser().resolve(),
    )


def _parse_server(obj: Any) -> dict[str, Any]:
    if obj is None:
        return {}
    if not isinstance(obj, dict):
        raise ValueError("server must be a JSON object or omitted.")
    return dict(obj)


def _parse_defaults(obj: Any) -> ServerDefaults:
    if not isinstance(obj, dict):
        raise ValueError("Missing or invalid 'defaults' section in config.")

    model_id = obj.get("model_id")
    if not isinstance(model_id, str) or not model_id:
        raise ValueError("defaults.model_id must be a non-empty string.")
    _validate_model_id(model_id)

    gm = obj.get("ground_motion")
    if not isinstance(gm, dict):
        raise ValueError("defaults.ground_motion must be a JSON object.")

    provider = gm.get("provider")
    if not isinstance(provider, str) or not provider:
        raise ValueError("defaults.ground_motion.provider must be a string.")

    gmpe_name = gm.get("gmpe_name")
    if gmpe_name is not None and not isinstance(gmpe_name, str):
        raise ValueError("defaults.ground_motion.gmpe_name must be a string.")

    distance_approx = gm.get("distance_approx")
    if distance_approx is not None and not isinstance(distance_approx, str):
        raise ValueError(
            "defaults.ground_motion.distance_approx must be a string."
        )

    impact_cfg = obj.get("impact_config", {})
    if not isinstance(impact_cfg, dict):
        raise ValueError("defaults.impact_config must be a JSON object.")

    return ServerDefaults(
        model_id=model_id,
        ground_motion=GroundMotionDefaults(
            provider=provider,
            gmpe_name=gmpe_name,
            distance_approx=distance_approx,
        ),
        impact_config=dict(impact_cfg),
    )


def _validate_paths(paths: ServerPaths) -> None:
    if paths.workdir.exists() and not paths.workdir.is_dir():
        raise ValueError(f"workdir is not a directory: {paths.workdir}")
    if paths.model_root.exists() and not paths.model_root.is_dir():
        raise ValueError(f"model_root is not a directory: {paths.model_root}")


def _validate_defaults(paths: ServerPaths, defaults: ServerDefaults) -> None:
    model_dir = resolve_model_dir(paths.model_root, defaults.model_id)
    _validate_model_dir(model_dir)


def _validate_model_id(model_id: str) -> None:
    if not _MODEL_ID_RE.match(model_id):
        raise ValueError(
            "Invalid model_id. Allowed: [A-Za-z0-9._-], "
            "max length 64, no slashes."
        )


def resolve_model_dir(model_root: Path, model_id: str) -> Path:
    """
    Resolve a model_id to a directory under model_root.

    The resolution is hardened against path traversal.
    """
    _validate_model_id(model_id)

    root = model_root.expanduser().resolve()
    model_dir = (root / model_id).resolve()

    try:
        model_dir.relative_to(root)
    except ValueError as exc:
        raise ValueError("model_id escapes model_root.") from exc

    return model_dir


def _validate_model_dir(model_dir: Path) -> None:
    """
    Validate that a model directory exists and contains required files.

    Required files (v1)
    -------------------
    - exposure.json
    - fragility.json
    - taxonomy_tree.json
    """
    if not model_dir.exists():
        raise ValueError(f"MODEL_NOT_FOUND: {model_dir}")
    if not model_dir.is_dir():
        raise ValueError(f"MODEL_INVALID: not a directory: {model_dir}")

    required = ("exposure.json", "fragility.json", "taxonomy_tree.json")
    missing: list[str] = []
    for name in required:
        p = model_dir / name
        if not p.exists() or not p.is_file():
            missing.append(name)

    if missing:
        miss = ", ".join(missing)
        raise ValueError(f"MODEL_INVALID: missing files: {miss}")


def model_paths(model_root: Path, model_id: str) -> dict[str, Path]:
    """
    Return required model file paths for a given model_id.
    """
    model_dir = resolve_model_dir(model_root, model_id)
    _validate_model_dir(model_dir)
    return {
        "model_dir": model_dir,
        "exposure_path": model_dir / "exposure.json",
        "fragility_path": model_dir / "fragility.json",
        "taxonomy_tree_path": model_dir / "taxonomy_tree.json",
    }


def list_models(model_root: Path) -> list[str]:
    """
    List valid model_id directories under model_root.

    A model_id is considered valid if it matches the model_id pattern and
    contains the required files.
    """
    root = model_root.expanduser().resolve()
    if not root.exists():
        return []
    if not root.is_dir():
        return []

    out: list[str] = []
    for p in sorted(root.iterdir()):
        if not p.is_dir():
            continue
        model_id = p.name
        try:
            _validate_model_id(model_id)
            _validate_model_dir(p)
        except ValueError:
            continue
        out.append(model_id)

    return out

