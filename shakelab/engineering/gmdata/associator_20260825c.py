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
ShakeLab Engineering ground-motion data association and resolution.

This module defines the engineering-level ground-motion configuration used
by impact calculations. It intentionally builds on the generic ground-motion
facilities provided by :mod:`shakelab.groundmotion.providers` without moving or
duplicating them.

Conceptual layers
-----------------
1) GroundMotionModel
   Pure configuration data loaded from a ShakeLabGroundMotion JSON file.

2) GroundMotionRuntime
   Runtime layer that instantiates each configured provider exactly once,
   binds providers to one ScenarioEvent through GroundMotionContext objects,
   and resolves the provider associated with each exposure asset.

3) GroundMotionResolution
   Immutable result of resolving one asset. It contains:
   - provider configuration id,
   - shared GroundMotionContext,
   - assignment-specific parameters.

Schema version 1.0.0
--------------------
A configuration contains:

- type = "ShakeLabGroundMotion"
- schema_version = "1.0.0"
- optional metadata
- a non-empty providers list

Each provider contains:
- id: unique configured-provider identifier
- provider: backend id registered in GroundMotionProvider
- optional default: boolean
- optional config: provider initialization parameters
- optional assignments

Each assignment contains:
- assets: non-empty list of exposure asset ids
- optional provider-specific fields at the same JSON level

For the ``imt`` provider, each assignment must contain exactly one
``station_code``.

Resolution rules
----------------
- An asset may appear in at most one explicit assignment.
- At most one provider may have default=true.
- Explicit assignment takes precedence over the default provider.
- If no explicit assignment exists:
    - the default provider is used if defined;
    - otherwise resolution fails.
- If an ExposureModel is supplied during validation, all referenced asset ids
  must exist in it.
- If no default provider exists and an ExposureModel is supplied, all exposure
  assets must be explicitly assigned.
- Schema v1 resolves exactly one provider per asset.

Notes
-----
- ScenarioEvent is deliberately not part of the JSON configuration.
- IMT is deliberately not part of the configuration: it is requested by the
  downstream fragility model.
- Multi-provider resolution for one asset is intentionally not implemented in
  schema v1. Its scientific semantics must be defined before extending the
  schema.
"""

from __future__ import annotations

import json

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Union

from shakelab.groundmotion.providers import (
    GroundMotionContext,
    GroundMotionProvider,
    ScenarioEvent,
)


EXPECTED_TYPE = "ShakeLabGroundMotion"
EXPECTED_SCHEMA_VERSION = "1.0.0"

__all__ = [
    "GroundMotionAssignment",
    "GroundMotionProviderConfig",
    "GroundMotionResolution",
    "GroundMotionModel",
    "GroundMotionRuntime",
    "load_ground_motion_model",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _is_non_empty_string(value: Any) -> bool:
    """Return True if *value* is a non-empty string."""
    return isinstance(value, str) and bool(value.strip())


def _as_dict(value: Any, field_name: str) -> Dict[str, Any]:
    """Return a plain dictionary or raise a descriptive ValueError."""
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"{field_name} must be a JSON object.")
    return dict(value)


def _as_string_list(value: Any, field_name: str) -> List[str]:
    """Validate and normalize a non-empty list of non-empty strings."""
    if not isinstance(value, list) or not value:
        raise ValueError(f"{field_name} must be a non-empty array.")

    out: List[str] = []
    for item in value:
        if not _is_non_empty_string(item):
            raise ValueError(
                f"{field_name} must contain only non-empty strings."
            )
        out.append(str(item).strip())

    return out


def _assignment_parameters(
    assignment: Mapping[str, Any],
    field_name: str,
) -> Dict[str, Any]:
    """
    Normalize provider-specific assignment fields.

    ``assets`` is the only structural assignment key. Every other direct
    field is treated as provider-specific runtime data.

    The legacy nested ``parameters`` object remains accepted for backward
    compatibility, but duplicate definitions are rejected.
    """
    if not isinstance(assignment, Mapping):
        raise ValueError(
            f"{field_name} must be a JSON object."
        )

    legacy = _as_dict(
        assignment.get("parameters", {}),
        f"{field_name}.parameters",
    )

    direct = {
        str(key): value
        for key, value in assignment.items()
        if key not in {"assets", "parameters"}
    }

    duplicate = set(legacy) & set(direct)

    if duplicate:
        raise ValueError(
            f"{field_name} defines fields both directly and inside "
            f"'parameters': {sorted(duplicate)}."
        )

    parameters = dict(legacy)
    parameters.update(direct)

    return parameters


def _resolve_provider_config(
    provider_id: str,
    config: Mapping[str, Any],
    *,
    base_directory: Path,
) -> Dict[str, Any]:
    """
    Normalize provider configuration requiring file-path resolution.

    Relative paths for file-backed providers are interpreted relative to the
    ShakeLabGroundMotion configuration file.
    """
    out = dict(config)

    if provider_id == "imt" and "root_path" in out:
        value = out["root_path"]

        if not _is_non_empty_string(value):
            raise ValueError(
                "IMT provider config.root_path must be a "
                "non-empty string."
            )

        resolved = Path(str(value)).expanduser()

        if not resolved.is_absolute():
            resolved = base_directory / resolved

        out["root_path"] = str(resolved.resolve())

    return out


# ---------------------------------------------------------------------------
# Pure configuration model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GroundMotionAssignment:
    """
    Assignment of one configured provider to one or more exposure assets.

    Parameters
    ----------
    assets
        Exposure asset identifiers assigned to the provider.
    parameters
        Provider-specific runtime parameters normalized by the loader.
        The public JSON representation is flat. For example,
        ``{"assets": ["BLDG_001"], "station_code": "UD01"}`` becomes
        ``parameters={"station_code": "UD01"}`` internally.
    """

    assets: tuple[str, ...]
    parameters: Mapping[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        """Validate assignment fields."""
        if not self.assets:
            raise ValueError(
                "GroundMotionAssignment.assets must not be empty."
            )

        for asset_id in self.assets:
            if not _is_non_empty_string(asset_id):
                raise ValueError(
                    "GroundMotionAssignment.assets must contain "
                    "non-empty strings."
                )

        if not isinstance(self.parameters, Mapping):
            raise ValueError(
                "GroundMotionAssignment.parameters must be a mapping."
            )


@dataclass(frozen=True)
class GroundMotionProviderConfig:
    """
    Configuration for one shared ground-motion provider.

    Parameters
    ----------
    id
        Unique identifier of this configured provider.
    provider
        Backend identifier registered in GroundMotionProvider.
    default
        If True, this provider is used for assets without explicit assignment.
    config
        Shared initialization parameters passed to GroundMotionProvider.from_id.
    assignments
        Optional explicit asset assignments for this provider.
    """

    id: str
    provider: str
    default: bool = False
    config: Mapping[str, Any] = field(default_factory=dict)
    assignments: tuple[GroundMotionAssignment, ...] = field(
        default_factory=tuple
    )

    def validate(self) -> None:
        """Validate provider configuration fields."""
        if not _is_non_empty_string(self.id):
            raise ValueError(
                "GroundMotionProviderConfig.id must be a non-empty string."
            )

        if not _is_non_empty_string(self.provider):
            raise ValueError(
                "GroundMotionProviderConfig.provider must be a "
                "non-empty string."
            )

        if not isinstance(self.default, bool):
            raise ValueError(
                "GroundMotionProviderConfig.default must be a boolean."
            )

        if not isinstance(self.config, Mapping):
            raise ValueError(
                "GroundMotionProviderConfig.config must be a mapping."
            )

        if not isinstance(self.assignments, tuple):
            raise ValueError(
                "GroundMotionProviderConfig.assignments must be a tuple."
            )

        for assignment in self.assignments:
            if not isinstance(assignment, GroundMotionAssignment):
                raise TypeError(
                    "GroundMotionProviderConfig.assignments must contain "
                    "GroundMotionAssignment objects."
                )
            assignment.validate()

        if self.provider == "imt":
            if self.default:
                raise ValueError(
                    "The 'imt' provider cannot be a default provider. "
                    "Each asset must be explicitly associated with one "
                    "station_code."
                )

            for assignment in self.assignments:
                station_code = assignment.parameters.get(
                    "station_code"
                )

                if not _is_non_empty_string(station_code):
                    raise ValueError(
                        "Each 'imt' assignment must define exactly one "
                        "non-empty station_code."
                    )

                if "station_codes" in assignment.parameters:
                    raise ValueError(
                        "Multiple-station association is not supported "
                        "in schema v1."
                    )


@dataclass(frozen=True)
class GroundMotionResolution:
    """
    Runtime resolution for one exposure asset.

    Parameters
    ----------
    provider_id
        Configured provider identifier.
    context
        Shared GroundMotionContext associated with the configured provider.
    parameters
        Assignment-specific runtime parameters.
    explicit
        True if the asset was explicitly assigned; False if the default
        provider was used.
    """

    provider_id: str
    context: GroundMotionContext
    parameters: Mapping[str, Any]
    explicit: bool


@dataclass
class GroundMotionModel:
    """
    Pure engineering-level ground-motion configuration.

    The model stores only configuration data and precomputed lookup indexes.
    It does not instantiate providers and does not bind them to an event.
    """

    type: str
    schema_version: str
    providers: List[GroundMotionProviderConfig]
    metadata: Dict[str, Any] = field(default_factory=dict)

    _provider_by_id: Dict[str, GroundMotionProviderConfig] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    _assignment_by_asset: Dict[
        str,
        tuple[str, Mapping[str, Any]],
    ] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    _default_provider_id: Optional[str] = field(
        default=None,
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        self._build_indexes()

    def _build_indexes(self) -> None:
        """Build provider and asset-resolution indexes."""
        self._provider_by_id = {}
        self._assignment_by_asset = {}
        self._default_provider_id = None

        for provider_cfg in self.providers:
            if provider_cfg.id not in self._provider_by_id:
                self._provider_by_id[provider_cfg.id] = provider_cfg

            if provider_cfg.default and self._default_provider_id is None:
                self._default_provider_id = provider_cfg.id

            for assignment in provider_cfg.assignments:
                for asset_id in assignment.assets:
                    if asset_id not in self._assignment_by_asset:
                        self._assignment_by_asset[asset_id] = (
                            provider_cfg.id,
                            dict(assignment.parameters),
                        )

    def validate(
        self,
        exposure_model: Optional[Any] = None,
        *,
        validate_provider_ids: bool = True,
    ) -> None:
        """
        Validate the complete configuration.

        Parameters
        ----------
        exposure_model
            Optional ExposureModel used for cross-validation of assignment
            asset ids and full coverage when no default provider is defined.
        validate_provider_ids
            If True, verify that each provider backend is currently registered
            in GroundMotionProvider.
        """
        if self.type != EXPECTED_TYPE:
            raise ValueError(
                f"GroundMotionModel.type must be {EXPECTED_TYPE!r}, "
                f"got {self.type!r}."
            )

        if self.schema_version != EXPECTED_SCHEMA_VERSION:
            raise ValueError(
                "GroundMotionModel.schema_version must be "
                f"{EXPECTED_SCHEMA_VERSION!r}, "
                f"got {self.schema_version!r}."
            )

        if not isinstance(self.metadata, dict):
            raise ValueError(
                "GroundMotionModel.metadata must be a dictionary."
            )

        if not isinstance(self.providers, list) or not self.providers:
            raise ValueError(
                "GroundMotionModel.providers must be a non-empty list."
            )

        provider_ids: set[str] = set()
        assigned_assets: set[str] = set()
        default_ids: List[str] = []

        available_provider_ids = set(GroundMotionProvider.available_ids())

        for provider_cfg in self.providers:
            if not isinstance(provider_cfg, GroundMotionProviderConfig):
                raise TypeError(
                    "GroundMotionModel.providers must contain "
                    "GroundMotionProviderConfig objects."
                )

            provider_cfg.validate()

            if provider_cfg.id in provider_ids:
                raise ValueError(
                    "Duplicate ground-motion provider id: "
                    f"{provider_cfg.id!r}."
                )
            provider_ids.add(provider_cfg.id)

            if validate_provider_ids:
                if provider_cfg.provider not in available_provider_ids:
                    raise ValueError(
                        "Unknown ground-motion provider backend "
                        f"{provider_cfg.provider!r} for configured provider "
                        f"{provider_cfg.id!r}. Available backends: "
                        f"{sorted(available_provider_ids)}."
                    )

            if provider_cfg.default:
                default_ids.append(provider_cfg.id)

            for assignment in provider_cfg.assignments:
                for asset_id in assignment.assets:
                    if asset_id in assigned_assets:
                        raise ValueError(
                            "Asset assigned more than once in ground-motion "
                            f"configuration: {asset_id!r}."
                        )
                    assigned_assets.add(asset_id)

        if len(default_ids) > 1:
            raise ValueError(
                "Only one ground-motion provider may have default=true. "
                f"Found: {default_ids}."
            )

        if exposure_model is not None:
            exposure_ids = self._exposure_asset_ids(exposure_model)

            unknown = assigned_assets - exposure_ids
            if unknown:
                raise ValueError(
                    "Ground-motion assignments reference unknown exposure "
                    f"assets: {sorted(unknown)}."
                )

            if not default_ids:
                missing = exposure_ids - assigned_assets
                if missing:
                    raise ValueError(
                        "Exposure assets without a ground-motion assignment "
                        "and no default provider is defined: "
                        f"{sorted(missing)}."
                    )

        self._build_indexes()

    @staticmethod
    def _exposure_asset_ids(exposure_model: Any) -> set[str]:
        """Extract exposure asset identifiers for cross-validation."""
        assets = getattr(exposure_model, "assets", None)
        if not isinstance(assets, Iterable):
            raise TypeError(
                "exposure_model must expose an iterable 'assets' attribute."
            )

        out: set[str] = set()

        for asset in assets:
            asset_id = getattr(asset, "id", None)
            if not _is_non_empty_string(asset_id):
                raise ValueError(
                    "Exposure contains an asset without a valid string id."
                )

            asset_key = str(asset_id).strip()

            if asset_key in out:
                raise ValueError(
                    f"Duplicate exposure asset id: {asset_key!r}."
                )

            out.add(asset_key)

        return out

    @property
    def default_provider_id(self) -> Optional[str]:
        """Return the configured default provider id, if any."""
        return self._default_provider_id

    def provider_config(
        self,
        provider_id: str,
    ) -> GroundMotionProviderConfig:
        """Return a provider configuration by id."""
        key = str(provider_id)

        try:
            return self._provider_by_id[key]
        except KeyError as exc:
            raise KeyError(
                f"Unknown ground-motion provider id: {key!r}."
            ) from exc

    def resolve_config(
        self,
        asset_id: str,
    ) -> tuple[GroundMotionProviderConfig, Dict[str, Any], bool]:
        """
        Resolve provider configuration for an asset.

        Returns
        -------
        provider_config, parameters, explicit
            The resolved provider configuration, assignment-specific
            parameters, and whether the assignment was explicit.
        """
        asset_key = str(asset_id)

        assigned = self._assignment_by_asset.get(asset_key)
        if assigned is not None:
            provider_id, parameters = assigned
            return (
                self.provider_config(provider_id),
                dict(parameters),
                True,
            )

        if self._default_provider_id is None:
            raise KeyError(
                "No ground-motion provider assigned to asset "
                f"{asset_key!r} and no default provider is defined."
            )

        return (
            self.provider_config(self._default_provider_id),
            {},
            False,
        )

    def runtime(
        self,
        event: ScenarioEvent,
    ) -> "GroundMotionRuntime":
        """
        Build a runtime bound to one ScenarioEvent.

        Each configured provider is instantiated exactly once.
        """
        return GroundMotionRuntime(
            model=self,
            event=event,
        )

    @classmethod
    def from_json(
        cls,
        json_path: Union[str, Path],
        *,
        exposure_model: Optional[Any] = None,
        validate: bool = True,
        validate_provider_ids: bool = True,
    ) -> "GroundMotionModel":
        """Load a ShakeLabGroundMotion model from JSON."""
        return load_ground_motion_model(
            json_path,
            exposure_model=exposure_model,
            validate=validate,
            validate_provider_ids=validate_provider_ids,
        )


# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------


class GroundMotionRuntime:
    """
    Runtime ground-motion resolver bound to one ScenarioEvent.

    One GroundMotionProvider object and one GroundMotionContext object are
    created for each configured provider. They are shared by every asset that
    resolves to that provider.
    """

    def __init__(
        self,
        model: GroundMotionModel,
        event: ScenarioEvent,
    ) -> None:
        if not isinstance(model, GroundMotionModel):
            raise TypeError("model must be a GroundMotionModel.")

        if not isinstance(event, ScenarioEvent):
            raise TypeError("event must be a ScenarioEvent.")

        self._model = model
        self._event = event

        self._providers: Dict[str, Any] = {}
        self._contexts: Dict[str, GroundMotionContext] = {}

        self._build_runtime()

    @property
    def model(self) -> GroundMotionModel:
        """Return the underlying pure configuration model."""
        return self._model

    @property
    def event(self) -> ScenarioEvent:
        """Return the ScenarioEvent associated with this runtime."""
        return self._event

    def _build_runtime(self) -> None:
        """Instantiate every configured provider exactly once."""
        for provider_cfg in self._model.providers:
            provider = GroundMotionProvider.from_id(
                provider_cfg.provider,
                config=dict(provider_cfg.config),
            )

            context = GroundMotionContext(
                event=self._event,
                provider=provider,
            )

            self._providers[provider_cfg.id] = provider
            self._contexts[provider_cfg.id] = context

    def provider(self, provider_id: str) -> Any:
        """Return the shared instantiated provider for a configured id."""
        key = str(provider_id)

        try:
            return self._providers[key]
        except KeyError as exc:
            raise KeyError(
                f"Unknown runtime ground-motion provider id: {key!r}."
            ) from exc

    def context(self, provider_id: str) -> GroundMotionContext:
        """Return the shared GroundMotionContext for a configured id."""
        key = str(provider_id)

        try:
            return self._contexts[key]
        except KeyError as exc:
            raise KeyError(
                f"Unknown runtime ground-motion provider id: {key!r}."
            ) from exc

    def resolve(self, asset_id: str) -> GroundMotionResolution:
        """Resolve the runtime ground-motion context for an asset."""
        provider_cfg, parameters, explicit = self._model.resolve_config(
            asset_id
        )

        return GroundMotionResolution(
            provider_id=provider_cfg.id,
            context=self.context(provider_cfg.id),
            parameters=parameters,
            explicit=explicit,
        )

    def evaluate_at_site(
        self,
        asset_id: str,
        imt: str,
        lon: float,
        lat: float,
        elevation_m: float = 0.0,
        **kwargs: Any,
    ) -> tuple[float, float]:
        """
        Resolve an asset and evaluate ground motion at one site.

        Assignment-specific provider data are normalized by the loader and
        merged with runtime keyword arguments. Explicit runtime keyword
        arguments take precedence.
        """
        resolved = self.resolve(asset_id)

        runtime_kwargs = dict(resolved.parameters)
        runtime_kwargs.update(kwargs)

        return resolved.context.evaluate_at_site(
            imt=imt,
            lon=lon,
            lat=lat,
            elevation_m=elevation_m,
            **runtime_kwargs,
        )

    def evaluate_at_asset(
        self,
        asset: Any,
        imt: str,
        **kwargs: Any,
    ) -> tuple[float, float]:
        """
        Resolve and evaluate ground motion at an exposure asset.

        The method uses asset.id for provider resolution and
        asset.reference_location for coordinates and elevation.
        """
        asset_id = getattr(asset, "id", None)
        if not _is_non_empty_string(asset_id):
            raise ValueError("Asset has no valid string id.")

        location = getattr(asset, "reference_location", None)
        if not isinstance(location, Mapping):
            raise ValueError(
                f"Asset {asset_id!r} has no valid reference_location."
            )

        try:
            lon = float(location["longitude"])
            lat = float(location["latitude"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Asset {asset_id!r} has invalid reference_location."
            ) from exc

        try:
            elevation_m = float(location.get("elevation", 0.0) or 0.0)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Asset {asset_id!r} has invalid elevation."
            ) from exc

        return self.evaluate_at_site(
            asset_id=str(asset_id),
            imt=imt,
            lon=lon,
            lat=lat,
            elevation_m=elevation_m,
            **kwargs,
        )


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_ground_motion_model(
    json_path: Union[str, Path],
    *,
    exposure_model: Optional[Any] = None,
    validate: bool = True,
    validate_provider_ids: bool = True,
) -> GroundMotionModel:
    """
    Load a ShakeLabGroundMotion JSON configuration.
    """
    path = Path(json_path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    if not path.is_file():
        raise ValueError(f"Not a file: {path}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(f"Cannot read ground-motion file: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON ground-motion file: {path}"
        ) from exc

    if not isinstance(data, Mapping):
        raise ValueError(
            "Ground-motion configuration root must be a JSON object."
        )

    for key in ("type", "schema_version", "providers"):
        if key not in data:
            raise ValueError(
                f"Missing required field {key!r} in ground-motion file."
            )

    metadata = _as_dict(
        data.get("metadata", {}),
        "metadata",
    )

    providers_raw = data.get("providers")
    if not isinstance(providers_raw, list) or not providers_raw:
        raise ValueError(
            "Field 'providers' must be a non-empty array."
        )

    providers: List[GroundMotionProviderConfig] = []

    for provider_index, provider_raw in enumerate(providers_raw):
        prefix = f"providers[{provider_index}]"

        if not isinstance(provider_raw, Mapping):
            raise ValueError(f"{prefix} must be a JSON object.")

        provider_id = provider_raw.get("id")
        provider_backend = provider_raw.get("provider")
        default = provider_raw.get("default", False)

        if not _is_non_empty_string(provider_id):
            raise ValueError(
                f"{prefix}.id must be a non-empty string."
            )

        if not _is_non_empty_string(provider_backend):
            raise ValueError(
                f"{prefix}.provider must be a non-empty string."
            )

        if not isinstance(default, bool):
            raise ValueError(
                f"{prefix}.default must be a boolean."
            )

        config = _as_dict(
            provider_raw.get("config", {}),
            f"{prefix}.config",
        )

        config = _resolve_provider_config(
            str(provider_backend).strip(),
            config,
            base_directory=path.parent,
        )

        assignments_raw = provider_raw.get("assignments", [])
        if assignments_raw is None:
            assignments_raw = []

        if not isinstance(assignments_raw, list):
            raise ValueError(
                f"{prefix}.assignments must be an array."
            )

        assignments: List[GroundMotionAssignment] = []

        for assignment_index, assignment_raw in enumerate(assignments_raw):
            aprefix = f"{prefix}.assignments[{assignment_index}]"

            if not isinstance(assignment_raw, Mapping):
                raise ValueError(
                    f"{aprefix} must be a JSON object."
                )

            if "assets" not in assignment_raw:
                raise ValueError(
                    f"Missing required field {aprefix}.assets."
                )

            assets = _as_string_list(
                assignment_raw.get("assets"),
                f"{aprefix}.assets",
            )

            parameters = _assignment_parameters(
                assignment_raw,
                aprefix,
            )

            assignments.append(
                GroundMotionAssignment(
                    assets=tuple(assets),
                    parameters=parameters,
                )
            )

        providers.append(
            GroundMotionProviderConfig(
                id=str(provider_id).strip(),
                provider=str(provider_backend).strip(),
                default=default,
                config=config,
                assignments=tuple(assignments),
            )
        )

    model = GroundMotionModel(
        type=str(data.get("type")),
        schema_version=str(data.get("schema_version")),
        metadata=metadata,
        providers=providers,
    )

    if validate:
        model.validate(
            exposure_model=exposure_model,
            validate_provider_ids=validate_provider_ids,
        )

    return model
