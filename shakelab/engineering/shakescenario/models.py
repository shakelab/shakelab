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
_MODEL_MANIFEST_NAME = "manifest.json"


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


@dataclass(frozen=True)
class ScenarioModel:
    """Loaded ShakeScenario model components."""

    model_id: str
    model_dir: Path
    exposure: Any
    taxonomy_tree: Any
    fragility: Any


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


def _load_model_manifest(model_dir: Path) -> dict[str, Any]:
    """
    Load and validate the model manifest.

    Parameters
    ----------
    model_dir
        Path to the model directory.

    Returns
    -------
    dict
        Parsed manifest content.

    Raises
    ------
    ValueError
        If the manifest is missing, invalid, or malformed.
    """
    manifest_path = model_dir / _MODEL_MANIFEST_NAME

    if not manifest_path.exists() or not manifest_path.is_file():
        raise ValueError(
            f"MODEL_INVALID: missing file: {_MODEL_MANIFEST_NAME}"
        )

    try:
        manifest = json.loads(
            manifest_path.read_text(encoding="utf-8")
        )
    except OSError as exc:
        raise ValueError(
            f"Cannot read model manifest: {manifest_path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON model manifest: {manifest_path}"
        ) from exc

    if not isinstance(manifest, dict):
        raise ValueError("Model manifest must be a JSON object.")

    schema_version = manifest.get("schema_version")
    if schema_version != "1.0.0":
        raise ValueError("Model manifest schema_version must be '1.0.0'.")

    files = manifest.get("files")
    if not isinstance(files, dict):
        raise ValueError("Model manifest must contain a 'files' object.")

    return manifest


def _resolve_model_file(model_dir: Path, relpath: Any, key: str) -> Path:
    """
    Resolve a manifest file path relative to the model directory.

    Parameters
    ----------
    model_dir
        Path to the model directory.
    relpath
        Relative path read from the manifest.
    key
        Manifest key used for error messages.

    Returns
    -------
    Path
        Resolved file path.

    Raises
    ------
    ValueError
        If the path is invalid, absolute, escapes model_dir, or is missing.
    """
    if not isinstance(relpath, str) or not relpath.strip():
        raise ValueError(f"files.{key} must be a non-empty string.")

    rel = Path(relpath)

    if rel.is_absolute():
        raise ValueError(f"files.{key} must be relative, not absolute.")

    path = (model_dir / rel).resolve()

    try:
        path.relative_to(model_dir.resolve())
    except ValueError as exc:
        raise ValueError(f"files.{key} escapes model directory.") from exc

    if not path.exists() or not path.is_file():
        raise ValueError(f"MODEL_INVALID: missing file: {relpath}")

    return path


def _manifest_model_paths(model_dir: Path) -> dict[str, Any]:
    """
    Return model file paths declared by the model manifest.

    Parameters
    ----------
    model_dir
        Path to the model directory.

    Returns
    -------
    dict
        Dictionary with model_dir, manifest_path, exposure_path,
        taxonomy_tree_path, and fragility_paths. The fragility_paths entry is
        a list of Path objects.
    """
    manifest = _load_model_manifest(model_dir)
    files = manifest["files"]

    exposure_path = _resolve_model_file(
        model_dir,
        files.get("exposure"),
        "exposure",
    )

    taxonomy_tree_path = _resolve_model_file(
        model_dir,
        files.get("taxonomy_tree"),
        "taxonomy_tree",
    )

    fragility_raw = files.get("fragility")
    if isinstance(fragility_raw, str):
        fragility_raw = [fragility_raw]

    if not isinstance(fragility_raw, list) or not fragility_raw:
        raise ValueError(
            "files.fragility must be a non-empty string or array."
        )

    fragility_paths = []
    for i, relpath in enumerate(fragility_raw):
        fragility_paths.append(
            _resolve_model_file(
                model_dir,
                relpath,
                f"fragility[{i}]",
            )
        )

    return {
        "model_dir": model_dir,
        "manifest_path": model_dir / _MODEL_MANIFEST_NAME,
        "exposure_path": exposure_path,
        "taxonomy_tree_path": taxonomy_tree_path,
        "fragility_paths": fragility_paths,
    }


def _validate_model_dir(model_dir: Path) -> None:
    """
    Validate that a model directory exists and contains a valid manifest.

    Required files (v1)
    -------------------
    - manifest.json

    The manifest must declare:
    - files.exposure
    - files.taxonomy_tree
    - files.fragility

    The declared file paths must be relative to the model directory.
    """
    if not model_dir.exists():
        raise ValueError(f"MODEL_NOT_FOUND: {model_dir}")
    if not model_dir.is_dir():
        raise ValueError(f"MODEL_INVALID: not a directory: {model_dir}")

    _manifest_model_paths(model_dir)


def model_paths(model_root: Path, model_id: str) -> dict[str, Any]:
    """
    Return model file paths for a given model_id.

    Paths are resolved from the model manifest.
    """
    model_dir = resolve_model_dir(model_root, model_id)

    if not model_dir.exists():
        raise ValueError(f"MODEL_NOT_FOUND: {model_dir}")
    if not model_dir.is_dir():
        raise ValueError(f"MODEL_INVALID: not a directory: {model_dir}")

    return _manifest_model_paths(model_dir)


def load_model(model_root: Path, model_id: str) -> ScenarioModel:
    """
    Load all model components for a given model_id.

    Parameters
    ----------
    model_root
        Root directory containing model directories.
    model_id
        Model identifier.

    Returns
    -------
    ScenarioModel
        Loaded exposure model, taxonomy tree, and fragility collection.

    Raises
    ------
    ValueError
        If the model directory or declared model files are invalid.
    """
    paths = model_paths(model_root, model_id)

    # Import here to keep server startup light and avoid import cycles.
    from shakelab.engineering.exposure.exposure import ExposureModel
    from shakelab.engineering.fragility.fragility import (
        FragilityCollection,
    )
    from shakelab.engineering.taxonomy.taxonomy_tree import (
        TaxonomyTree,
    )

    exposure = ExposureModel.from_json(
        str(paths["exposure_path"]),
        validate=True,
    )
    taxonomy_tree = TaxonomyTree.from_json(
        str(paths["taxonomy_tree_path"])
    )
    fragility = FragilityCollection.from_json(
        paths["fragility_paths"]
    )

    return ScenarioModel(
        model_id=model_id,
        model_dir=paths["model_dir"],
        exposure=exposure,
        taxonomy_tree=taxonomy_tree,
        fragility=fragility,
    )


def list_models(model_root: Path) -> list[str]:
    """
    List valid model_id directories under model_root.

    A model_id is considered valid if it matches the model_id pattern and
    contains a valid model manifest.
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

