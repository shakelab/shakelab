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
General ground-motion interface and providers.

This module defines a provider-independent interface for computing intensity
measures at arbitrary sites for a scenario event.

The scientific providers implemented here are intentionally independent of
the engineering exposure/fragility/impact layers. Engineering-level provider configuration and asset routing are
implemented in ``shakelab.engineering.gmdata.associator``.

Providers
---------
``gmpe``
    Analytical/statistical GMPEs implemented in ShakeLab.

``shakemap``
    Precomputed ShakeMap products read through ``shakelab.gmmodel.usgs``.

``imt``
    Precomputed ShakeScenarioIMT datasets addressed by station code.

``plugin``
    Placeholder for external or numerical ground-motion backends.

Output convention
-----------------
Providers return ``(im, sigma_ln)``:

- ``im`` is the median intensity measure in linear space;
- ``sigma_ln`` is the standard deviation of ln(IM).

A fully general unit system is not yet enforced by this module. Individual
providers must document their output convention explicitly.

For the current ShakeMap provider, the default output convention is chosen to
match the existing ShakeLab fragility models:

- PGA and spectral acceleration: ``g``;
- PGV: ``cm/s``;
- MMI: dimensionless.

ShakeMap XML acceleration values expressed in percent-g are therefore divided
by 100. This policy is explicit and configurable pending a general ShakeLab
unit framework.
"""

from __future__ import annotations

import json
import re

from dataclasses import dataclass
from math import exp, isfinite
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

from shakelab.libutils.constants import GRAVITY
from shakelab.libutils.geodeticN.primitives import WgsPoint


# ---------------------------------------------------------------------------
# Scenario primitives
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ScenarioEvent:
    """
    Seismic event container used by ground-motion providers.

    Parameters
    ----------
    magnitude
        Event magnitude.
    hypocentre
        Hypocentral location as WgsPoint.
    origin_time
        Optional origin time.
    mechanism
        Optional focal-mechanism mapping.
    event_id
        Optional event identifier. Providers backed by event-specific
        products, such as ShakeMap, may require it.
    magnitude_type
        Optional magnitude type.
    """

    magnitude: float
    hypocentre: WgsPoint
    origin_time: Optional[Any] = None
    mechanism: Optional[Mapping[str, float]] = None
    event_id: Optional[str] = None
    magnitude_type: Optional[str] = None

    def __post_init__(self) -> None:
        self.validate()

    @property
    def longitude(self) -> float:
        """Return hypocentre longitude in degrees."""
        return float(self.hypocentre.longitude)

    @property
    def latitude(self) -> float:
        """Return hypocentre latitude in degrees."""
        return float(self.hypocentre.latitude)

    @property
    def depth_km(self) -> float:
        """Return hypocentral depth in km, positive downward."""
        return -float(self.hypocentre.elevation) / 1000.0

    @property
    def epicentre(self) -> WgsPoint:
        """Return surface projection of the hypocentre."""
        return WgsPoint(
            longitude=self.longitude,
            latitude=self.latitude,
            elevation=0.0,
        )

    @classmethod
    def from_lonlat_depth(
        cls,
        magnitude: float,
        longitude: float,
        latitude: float,
        depth_km: float,
        **kwargs: Any,
    ) -> "ScenarioEvent":
        """Build an event from longitude, latitude and depth."""
        return cls(
            magnitude=float(magnitude),
            hypocentre=WgsPoint(
                longitude=float(longitude),
                latitude=float(latitude),
                elevation=-float(depth_km) * 1000.0,
            ),
            **kwargs,
        )

    def validate(self) -> None:
        """Validate basic event geometry and metadata."""
        if (
            not isinstance(self.magnitude, (int, float))
            or isinstance(self.magnitude, bool)
        ):
            raise TypeError("magnitude must be a number.")

        if not isfinite(float(self.magnitude)):
            raise ValueError("magnitude must be finite.")

        if not isinstance(self.hypocentre, WgsPoint):
            raise TypeError("hypocentre must be a WgsPoint.")

        lon = float(self.hypocentre.longitude)
        lat = float(self.hypocentre.latitude)
        elev = float(self.hypocentre.elevation)

        if not (isfinite(lon) and isfinite(lat) and isfinite(elev)):
            raise ValueError(
                "hypocentre longitude, latitude and elevation "
                "must be finite."
            )

        if not -180.0 <= lon <= 180.0:
            raise ValueError("longitude out of range [-180, 180].")

        if not -90.0 <= lat <= 90.0:
            raise ValueError("latitude out of range [-90, 90].")

        if self.mechanism is not None and not isinstance(
            self.mechanism,
            Mapping,
        ):
            raise TypeError("mechanism must be a mapping or None.")

        if self.event_id is not None:
            if not isinstance(self.event_id, str):
                raise TypeError("event_id must be a string or None.")
            if not self.event_id.strip():
                raise ValueError("event_id must not be empty.")


# ---------------------------------------------------------------------------
# Provider registry and interface
# ---------------------------------------------------------------------------


class GroundMotionProvider:
    """Factory/registry for ground-motion evaluation backends."""

    _REGISTRY: Dict[str, Callable[..., "_BaseProvider"]] = {}

    @classmethod
    def register(
        cls,
        provider_id: str,
    ) -> Callable[..., Any]:
        """Register a provider constructor."""
        if not isinstance(provider_id, str) or not provider_id.strip():
            raise ValueError("provider_id must be a non-empty string.")

        key = provider_id.strip()

        def _decorator(
            provider: Callable[..., "_BaseProvider"],
        ) -> Callable[..., Any]:
            cls._REGISTRY[key] = provider
            return provider

        return _decorator

    @classmethod
    def available_ids(cls) -> List[str]:
        """Return registered provider identifiers."""
        return sorted(cls._REGISTRY.keys())

    @classmethod
    def from_id(
        cls,
        provider_id: str,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """Instantiate a registered provider."""
        if provider_id not in cls._REGISTRY:
            raise KeyError(
                f"Unknown provider_id: {provider_id!r}. Available: "
                f"{', '.join(cls.available_ids())}"
            )

        if config is None:
            cfg: Dict[str, Any] = {}
        else:
            if not isinstance(config, Mapping):
                raise TypeError("config must be a mapping.")
            cfg = dict(config)

        return cls._REGISTRY[provider_id](**cfg)

    @classmethod
    def gmpe(
        cls,
        gmpe_name: str,
        *,
        distance_approx: str = "ellipsoid",
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """Convenience constructor for the GMPE provider."""
        cfg = dict(config or {})
        cfg["gmpe_name"] = gmpe_name
        cfg["distance_approx"] = distance_approx
        return cls.from_id("gmpe", config=cfg)

    @classmethod
    def shakemap(
        cls,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """Convenience constructor for the ShakeMap provider."""
        return cls.from_id("shakemap", config=config)

    @classmethod
    def imt(
        cls,
        root_path: str,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """Convenience constructor for the precomputed IMT provider."""
        cfg = dict(config or {})
        cfg["root_path"] = root_path
        return cls.from_id("imt", config=cfg)

    @classmethod
    def plugin(
        cls,
        plugin_id: str,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """Convenience constructor for plugin providers."""
        cfg = dict(config or {})
        cfg["plugin_id"] = plugin_id
        return cls.from_id("plugin", config=cfg)


class _BaseProvider:
    """Internal provider interface."""

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """Evaluate median IM and sigma_ln."""
        raise NotImplementedError


# ---------------------------------------------------------------------------
# GMPE provider
# ---------------------------------------------------------------------------


@GroundMotionProvider.register("gmpe")
class _GmpeProvider(_BaseProvider):
    """Ground-motion provider backed by the ShakeLab GMPE registry."""

    def __init__(
        self,
        gmpe_name: str,
        distance_approx: str = "ellipsoid",
        **gmpe_kwargs: Any,
    ) -> None:
        from shakelab.gmmodel.gmpe.registry import create_gmpe

        self._gmpe = create_gmpe(
            gmpe_name,
            **gmpe_kwargs,
        )

        approx = str(distance_approx).strip().lower()
        self._distance_approx = approx or "ellipsoid"

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """Evaluate a GMPE at all requested sites."""
        _ = kwargs

        metric = self._get_distance_metric()

        values: List[float] = []
        sigma_ln: List[float] = []

        for site in sites:
            if not isinstance(site, WgsPoint):
                raise TypeError("sites must contain WgsPoint objects.")

            distance_km = self._compute_distance_km(
                event,
                site,
                metric,
                approx=self._distance_approx,
            )

            mean_ln, sig_ln = self._gmpe.ground_motion(
                imt,
                event.magnitude,
                distance_km,
            )

            mean_ln_f = float(mean_ln)
            sig_ln_f = float(sig_ln)

            if isfinite(mean_ln_f):
                values.append(exp(mean_ln_f))
            else:
                values.append(float("nan"))

            if not isfinite(sig_ln_f) or sig_ln_f < 0.0:
                sig_ln_f = float("nan")

            sigma_ln.append(sig_ln_f)

        return values, sigma_ln

    def _get_distance_metric(self) -> str:
        metric = str(
            getattr(self._gmpe, "DISTANCE_METRIC")
        ).strip().lower()

        if metric in {"repi", "epicentral"}:
            return "epicentral"

        if metric in {"rhypo", "hypocentral"}:
            return "hypocentral"

        raise ValueError(
            f"Unsupported GMPE distance metric: {metric!r}. "
            "Expected 'epicentral/repi' or 'hypocentral/rhypo'."
        )

    @staticmethod
    def _compute_distance_km(
        event: ScenarioEvent,
        site: WgsPoint,
        metric: str,
        approx: str = "ellipsoid",
    ) -> float:
        if metric == "epicentral":
            return (
                float(
                    event.epicentre.epicentral_distance_to(site)
                )
                / 1000.0
            )

        if metric == "hypocentral":
            return (
                float(
                    event.hypocentre.hypocentral_distance_to(
                        site,
                        approx=approx,
                    )
                )
                / 1000.0
            )

        raise ValueError(
            f"Unsupported distance metric: {metric!r}."
        )


# ---------------------------------------------------------------------------
# ShakeMap provider
# ---------------------------------------------------------------------------


@GroundMotionProvider.register("shakemap")
class _ShakeMapProvider(_BaseProvider):
    """
    Provider backed by event-specific ShakeMap XML products.

    Parameters
    ----------
    root_path
        Root directory containing one subdirectory per ShakeMap event id.
        For ``event.event_id == "183907"``, the provider reads
        ``root_path/183907``.
    interp
        Spatial interpolation method: ``linear`` or ``nearest``.
    outside
        Behavior outside the ShakeMap domain: ``raise`` or ``nan``.
    load_uncertainty
        Load ``uncertainty.xml`` when present.
    require_uncertainty
        If True, fail when the uncertainty field for the requested IMT is
        unavailable. If False, missing uncertainty is returned as 0.0.
    require_info
        Require ``info.json`` in every event directory.
    check_event_id
        Compare event id from ``info.json`` with the requested event id when
        metadata are available.
    acceleration_unit
        Output acceleration unit. Supported: ``g`` and ``percent_g``.
        The default ``g`` matches current ShakeLab fragility examples.
    velocity_unit
        Output velocity unit. Supported: ``cm/s`` and ``m/s``.
    """

    _IMT_FIELDS = {
        "PGA": ("PGA", "STDPGA", "acceleration"),
        "PGV": ("PGV", "STDPGV", "velocity"),
        "MMI": ("MMI", "STDMMI", "intensity"),
        "SA(0.3)": ("PSA03", "STDPSA03", "acceleration"),
        "SA(1.0)": ("PSA10", "STDPSA10", "acceleration"),
        "SA(3.0)": ("PSA30", "STDPSA30", "acceleration"),
    }

    def __init__(
        self,
        root_path: str,
        interp: str = "linear",
        outside: str = "raise",
        load_uncertainty: bool = True,
        require_uncertainty: bool = False,
        require_info: bool = False,
        acceleration_unit: str = "g",
        velocity_unit: str = "cm/s",
        **kwargs: Any,
    ) -> None:
        self.root_path = Path(root_path).expanduser().resolve()

        if not self.root_path.exists():
            raise FileNotFoundError(
                f"ShakeMap root directory not found: "
                f"{self.root_path}"
            )

        if not self.root_path.is_dir():
            raise ValueError(
                f"ShakeMap root path is not a directory: "
                f"{self.root_path}"
            )

        interp_key = str(interp).strip().lower()
        if interp_key == "bilinear":
            interp_key = "linear"

        if interp_key not in {"linear", "nearest"}:
            raise ValueError(
                "ShakeMap interp must be 'linear' or 'nearest'."
            )

        outside_key = str(outside).strip().lower()
        if outside_key not in {"raise", "nan"}:
            raise ValueError(
                "ShakeMap outside must be 'raise' or 'nan'."
            )

        acceleration_key = _normalize_acceleration_unit(
            acceleration_unit
        )
        velocity_key = _normalize_velocity_unit(
            velocity_unit
        )

        self.interp = interp_key
        self.outside = outside_key
        self.load_uncertainty = bool(load_uncertainty)
        self.require_uncertainty = bool(require_uncertainty)
        self.require_info = bool(require_info)
        self.check_event_id = bool(check_event_id)
        self.acceleration_unit = acceleration_key
        self.velocity_unit = velocity_key
        self.kwargs = dict(kwargs)

        self._cache: Dict[str, Any] = {}

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Interpolate ShakeMap median IM and sigma_ln at requested sites.

        Runtime keyword arguments may override ``interp`` and ``outside``.
        """
        from shakelab.gmmodel.usgs import ShakeMapEvent

        _ = ShakeMapEvent

        event_id = self._event_id(event)
        shakemap = self._load_event(event_id)

        imt_key = _normalize_imt(imt)
        field, sigma_field, quantity = self._field_mapping(imt_key)

        if shakemap.gm_data is None:
            raise RuntimeError(
                f"ShakeMap event {event_id!r} has no grid data."
            )

        if not shakemap.gm_data.has_param(field):
            available = ", ".join(shakemap.gm_data.field_names)
            raise KeyError(
                f"ShakeMap field {field!r} required for IMT "
                f"{imt_key!r} is unavailable. Available: {available}."
            )

        has_sigma = shakemap.gm_data.has_param(sigma_field)

        if self.require_uncertainty and not has_sigma:
            raise KeyError(
                f"ShakeMap uncertainty field {sigma_field!r} "
                f"required for IMT {imt_key!r} is unavailable."
            )

        parameters = [field]
        if has_sigma:
            parameters.append(sigma_field)

        interp = str(
            kwargs.get("interp", self.interp)
        ).strip().lower()
        if interp == "bilinear":
            interp = "linear"

        outside = str(
            kwargs.get("outside", self.outside)
        ).strip().lower()

        coordinates: List[Tuple[float, float]] = []

        for site in sites:
            if not isinstance(site, WgsPoint):
                raise TypeError(
                    "sites must contain WgsPoint objects."
                )
            coordinates.append(
                (
                    float(site.longitude),
                    float(site.latitude),
                )
            )

        interpolated = shakemap.get_ground_motion(
            coordinates,
            parameters=parameters,
            method=interp,
            outside=outside,
        )

        native_unit = shakemap.units.get(field, "")

        values = [
            _convert_shakemap_value(
                row[field],
                quantity=quantity,
                native_unit=native_unit,
                acceleration_unit=self.acceleration_unit,
                velocity_unit=self.velocity_unit,
            )
            for row in interpolated
        ]

        if has_sigma:
            sigma_ln = [
                _validate_sigma_ln(row[sigma_field])
                for row in interpolated
            ]
        else:
            sigma_ln = [
                0.0
                for _ in interpolated
            ]

        return values, sigma_ln

    def event_directory(
        self,
        event_id: str,
    ) -> Path:
        """Resolve an event id safely under the configured root."""
        key = _validate_event_directory_id(event_id)

        event_dir = (self.root_path / key).resolve()

        try:
            event_dir.relative_to(self.root_path)
        except ValueError as exc:
            raise ValueError(
                "ShakeMap event directory escapes root_path."
            ) from exc

        return event_dir

    def clear_cache(
        self,
        event_id: Optional[str] = None,
    ) -> None:
        """Clear all cached events or one selected event."""
        if event_id is None:
            self._cache.clear()
            return

        self._cache.pop(str(event_id), None)

    def _event_id(
        self,
        event: ScenarioEvent,
    ) -> str:
        if event.event_id is None:
            raise ValueError(
                "ShakeMap provider requires ScenarioEvent.event_id."
            )

        return _validate_event_directory_id(
            event.event_id
        )

    def _load_event(
        self,
        event_id: str,
    ) -> Any:
        from shakelab.gmmodel.usgs import ShakeMapEvent

        cached = self._cache.get(event_id)
        if cached is not None:
            return cached

        event_dir = self.event_directory(event_id)

        if not event_dir.exists():
            raise FileNotFoundError(
                f"ShakeMap event directory not found: {event_dir}"
            )

        if not event_dir.is_dir():
            raise ValueError(
                f"ShakeMap event path is not a directory: "
                f"{event_dir}"
            )

        shakemap = ShakeMapEvent(
            folder=event_dir,
            load_stations=False,
            load_uncertainty=self.load_uncertainty,
            require_info=self.require_info,
        )

        if (
            self.check_event_id
            and shakemap.info is not None
            and shakemap.event_id not in {"", "N/A"}
            and str(shakemap.event_id) != event_id
        ):
            raise ValueError(
                "ShakeMap event-id mismatch: requested "
                f"{event_id!r}, info.json contains "
                f"{shakemap.event_id!r}."
            )

        self._cache[event_id] = shakemap
        return shakemap

    @classmethod
    def _field_mapping(
        cls,
        imt: str,
    ) -> Tuple[str, str, str]:
        try:
            return cls._IMT_FIELDS[imt]
        except KeyError as exc:
            supported = ", ".join(cls._IMT_FIELDS.keys())
            raise ValueError(
                f"Unsupported ShakeMap IMT: {imt!r}. "
                f"Supported: {supported}."
            ) from exc


# ---------------------------------------------------------------------------
# Precomputed IMT provider
# ---------------------------------------------------------------------------


IMT_EXPECTED_TYPE = "ShakeScenarioIMT"
IMT_EXPECTED_SCHEMA_VERSION = "1.0.0"


@GroundMotionProvider.register("imt")
class _IMTProvider(_BaseProvider):
    """
    Provider backed by event-specific ShakeScenarioIMT JSON products.

    Parameters
    ----------
    root_path
        Root directory containing one JSON file per event id. For
        ``event.event_id == "friuli_1976"``, the provider reads
        ``root_path/friuli_1976.json``.
    acceleration_unit
        Output unit for PGA and SA. Supported values are ``g`` and ``m/s2``.
        The default ``g`` matches the current ShakeLab fragility convention.
    velocity_unit
        Output unit for PGV. Supported values are ``cm/s`` and ``m/s``.
        The default ``cm/s`` matches the current ShakeLab fragility
        convention.
    period_tolerance
        Absolute tolerance in seconds used when matching requested SA periods
        against the period vector stored in the selected dataset.

    Notes
    -----
    The station associated with an exposure asset is supplied at runtime
    through the ``station_code`` keyword argument. Site coordinates are not
    used to select or interpolate values.

    Event datasets are loaded lazily and cached by event id. The provider is
    deterministic and therefore returns ``sigma_ln = 0.0``.
    """

    def __init__(
        self,
        root_path: str,
        acceleration_unit: str = "g",
        velocity_unit: str = "cm/s",
        period_tolerance: float = 1.0e-9,
        **kwargs: Any,
    ) -> None:
        _ = kwargs

        self.root_path = Path(root_path).expanduser().resolve()

        if not self.root_path.exists():
            raise FileNotFoundError(
                "ShakeScenarioIMT root directory not found: "
                f"{self.root_path}"
            )

        if not self.root_path.is_dir():
            raise ValueError(
                "ShakeScenarioIMT root path is not a directory: "
                f"{self.root_path}"
            )

        self.acceleration_unit = _normalize_imt_acceleration_unit(
            acceleration_unit
        )
        self.velocity_unit = _normalize_velocity_unit(
            velocity_unit
        )

        self.period_tolerance = float(period_tolerance)

        if (
            not isfinite(self.period_tolerance)
            or self.period_tolerance < 0.0
        ):
            raise ValueError(
                "period_tolerance must be finite and non-negative."
            )

        self._cache: Dict[str, Dict[str, Any]] = {}

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Return one precomputed IMT value for the assigned station.

        The common provider interface still receives ``sites`` even though
        this provider selects values exclusively through ``station_code``.
        """
        event_id = self._event_id(event)
        dataset = self._load_event(event_id)

        station_code = kwargs.get("station_code")
        station = self._station(
            dataset,
            station_code,
        )

        value = self._station_imt_value(
            dataset,
            station,
            imt,
        )

        count = len(sites)

        return (
            [value for _ in range(count)],
            [0.0 for _ in range(count)],
        )

    def event_file(
        self,
        event_id: str,
    ) -> Path:
        """Resolve one event-specific IMT file safely under root_path."""
        key = _validate_imt_event_id(event_id)

        event_file = (self.root_path / f"{key}.json").resolve()

        try:
            event_file.relative_to(self.root_path)
        except ValueError as exc:
            raise ValueError(
                "ShakeScenarioIMT event file escapes root_path."
            ) from exc

        return event_file

    def clear_cache(
        self,
        event_id: Optional[str] = None,
    ) -> None:
        """Clear all cached IMT datasets or one selected event."""
        if event_id is None:
            self._cache.clear()
            return

        self._cache.pop(str(event_id), None)

    @staticmethod
    def _event_id(
        event: ScenarioEvent,
    ) -> str:
        """Return the event id required to select an IMT dataset."""
        if event.event_id is None:
            raise ValueError(
                "IMT provider requires ScenarioEvent.event_id."
            )

        return _validate_imt_event_id(
            event.event_id
        )

    def _load_event(
        self,
        event_id: str,
    ) -> Dict[str, Any]:
        """Load, validate and cache one event-specific IMT dataset."""
        cached = self._cache.get(event_id)

        if cached is not None:
            return cached

        event_file = self.event_file(event_id)

        if not event_file.exists():
            raise FileNotFoundError(
                f"ShakeScenarioIMT event file not found: {event_file}"
            )

        if not event_file.is_file():
            raise ValueError(
                "ShakeScenarioIMT event path is not a file: "
                f"{event_file}"
            )

        try:
            data = json.loads(
                event_file.read_text(encoding="utf-8")
            )
        except OSError as exc:
            raise ValueError(
                f"Cannot read ShakeScenarioIMT file: {event_file}"
            ) from exc
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid ShakeScenarioIMT JSON file: {event_file}"
            ) from exc

        if not isinstance(data, Mapping):
            raise ValueError(
                "ShakeScenarioIMT root must be a JSON object."
            )

        if data.get("type") != IMT_EXPECTED_TYPE:
            raise ValueError(
                "Invalid IMT dataset type: "
                f"{data.get('type')!r}; expected "
                f"{IMT_EXPECTED_TYPE!r}."
            )

        if data.get("schema_version") != IMT_EXPECTED_SCHEMA_VERSION:
            raise ValueError(
                "Unsupported ShakeScenarioIMT schema_version: "
                f"{data.get('schema_version')!r}; expected "
                f"{IMT_EXPECTED_SCHEMA_VERSION!r}."
            )

        if not isinstance(data.get("imts"), Mapping):
            raise ValueError(
                "ShakeScenarioIMT must contain an 'imts' object."
            )

        if not isinstance(data.get("stations"), list):
            raise ValueError(
                "ShakeScenarioIMT must contain a 'stations' array."
            )

        dataset = {
            "path": event_file,
            "imts": self._build_imt_definitions(data),
            "stations": self._build_station_index(data),
        }

        self._cache[event_id] = dataset

        return dataset

    @staticmethod
    def _build_imt_definitions(
        data: Mapping[str, Any],
    ) -> Dict[str, Dict[str, Any]]:
        """Validate and normalize available IMT definitions."""
        definitions: Dict[str, Dict[str, Any]] = {}

        for key, raw_definition in data["imts"].items():
            name = str(key).strip().upper()

            if not name:
                raise ValueError(
                    "ShakeScenarioIMT contains an empty IMT name."
                )

            if not isinstance(raw_definition, Mapping):
                raise ValueError(
                    f"IMT definition {key!r} must be a JSON object."
                )

            definition = dict(raw_definition)

            units = definition.get("units")

            if not isinstance(units, str) or not units.strip():
                raise ValueError(
                    f"IMT {key!r} must define non-empty 'units'."
                )

            if name == "SA":
                periods = definition.get("periods")

                if not isinstance(periods, list) or not periods:
                    raise ValueError(
                        "SA definition must contain a non-empty "
                        "'periods' array."
                    )

                normalized_periods: List[float] = []

                for period in periods:
                    value = float(period)

                    if not isfinite(value) or value <= 0.0:
                        raise ValueError(
                            "SA periods must be finite and positive."
                        )

                    normalized_periods.append(value)

                definition["periods"] = normalized_periods

            definitions[name] = definition

        return definitions

    @staticmethod
    def _build_station_index(
        data: Mapping[str, Any],
    ) -> Dict[str, Dict[str, Any]]:
        """Build a station-code lookup index for one dataset."""
        index: Dict[str, Dict[str, Any]] = {}

        for position, raw_station in enumerate(data["stations"]):
            if not isinstance(raw_station, Mapping):
                raise ValueError(
                    f"stations[{position}] must be a JSON object."
                )

            station_code = raw_station.get("station_code")

            if (
                not isinstance(station_code, str)
                or not station_code.strip()
            ):
                raise ValueError(
                    f"stations[{position}].station_code must be a "
                    "non-empty string."
                )

            key = station_code.strip()

            if key in index:
                raise ValueError(
                    f"Duplicate IMT station_code: {key!r}."
                )

            values = raw_station.get("values")

            if not isinstance(values, Mapping):
                raise ValueError(
                    f"Station {key!r} must contain a 'values' object."
                )

            station = dict(raw_station)
            station["values"] = dict(values)

            index[key] = station

        return index

    @staticmethod
    def _station(
        dataset: Mapping[str, Any],
        station_code: Any,
    ) -> Dict[str, Any]:
        """Return one station from the selected precomputed dataset."""
        if not isinstance(station_code, str):
            raise TypeError(
                "IMT provider requires assignment field "
                "'station_code' as a string."
            )

        key = station_code.strip()

        if not key:
            raise ValueError(
                "IMT provider requires a non-empty 'station_code'."
            )

        try:
            return dataset["stations"][key]
        except KeyError as exc:
            path = Path(dataset["path"])

            raise KeyError(
                f"Station {key!r} not found in ShakeScenarioIMT "
                f"dataset {path.name!r}."
            ) from exc

    def _station_imt_value(
        self,
        dataset: Mapping[str, Any],
        station: Mapping[str, Any],
        imt: str,
    ) -> float:
        """Extract and convert one requested station IMT."""
        normalized, period = _parse_precomputed_imt(imt)
        values = station["values"]

        if normalized in {"PGA", "PGV"}:
            if normalized not in values:
                raise KeyError(
                    f"Station {station['station_code']!r} has no "
                    f"{normalized} value."
                )

            value = _finite_float(
                values[normalized],
                f"{station['station_code']}.{normalized}",
            )

            definition = self._imt_definition(
                dataset,
                normalized,
            )

            if normalized == "PGA":
                return self._convert_acceleration(
                    value,
                    definition["units"],
                )

            return self._convert_velocity(
                value,
                definition["units"],
            )

        definition = self._imt_definition(
            dataset,
            "SA",
        )
        raw_values = values.get("SA")

        if not isinstance(raw_values, list):
            raise ValueError(
                f"Station {station['station_code']!r} has no valid "
                "SA array."
            )

        periods = definition["periods"]

        if len(raw_values) != len(periods):
            raise ValueError(
                f"Station {station['station_code']!r} SA length "
                "does not match the dataset period vector."
            )

        index = self._period_index(
            float(period),
            periods,
        )

        value = _finite_float(
            raw_values[index],
            f"{station['station_code']}.SA({float(period):g})",
        )

        return self._convert_acceleration(
            value,
            definition["units"],
        )

    @staticmethod
    def _imt_definition(
        dataset: Mapping[str, Any],
        imt: str,
    ) -> Dict[str, Any]:
        """Return one validated IMT definition."""
        try:
            return dataset["imts"][imt]
        except KeyError as exc:
            path = Path(dataset["path"])

            raise ValueError(
                f"IMT {imt!r} is not available in "
                f"{path.name!r}."
            ) from exc

    def _period_index(
        self,
        requested: float,
        periods: Sequence[float],
    ) -> int:
        """Return the index of the requested spectral period."""
        best_index: Optional[int] = None
        best_delta: Optional[float] = None

        for index, period in enumerate(periods):
            delta = abs(float(period) - requested)

            if best_delta is None or delta < best_delta:
                best_index = index
                best_delta = delta

        if (
            best_index is None
            or best_delta is None
            or best_delta > self.period_tolerance
        ):
            available = ", ".join(
                f"{float(period):g}"
                for period in periods
            )

            raise ValueError(
                f"SA period {requested:g} s is not available. "
                f"Available periods: {available} s."
            )

        return best_index

    def _convert_acceleration(
        self,
        value: float,
        input_unit: str,
    ) -> float:
        """Convert acceleration to the configured output unit."""
        source = _normalize_imt_acceleration_unit(input_unit)

        if source == self.acceleration_unit:
            return value

        if source == "m/s2" and self.acceleration_unit == "g":
            return value / GRAVITY

        if source == "g" and self.acceleration_unit == "m/s2":
            return value * GRAVITY

        raise ValueError(
            "Unsupported acceleration conversion: "
            f"{source!r} -> {self.acceleration_unit!r}."
        )

    def _convert_velocity(
        self,
        value: float,
        input_unit: str,
    ) -> float:
        """Convert velocity to the configured output unit."""
        source = _normalize_velocity_unit(input_unit)

        if source == self.velocity_unit:
            return value

        if source == "m/s" and self.velocity_unit == "cm/s":
            return value * 100.0

        if source == "cm/s" and self.velocity_unit == "m/s":
            return value / 100.0

        raise ValueError(
            "Unsupported velocity conversion: "
            f"{source!r} -> {self.velocity_unit!r}."
        )


# ---------------------------------------------------------------------------
# Plugin provider
# ---------------------------------------------------------------------------


@GroundMotionProvider.register("plugin")
class _PluginProvider(_BaseProvider):
    """Placeholder provider for external or numerical backends."""

    def __init__(
        self,
        plugin_id: str,
        **kwargs: Any,
    ) -> None:
        self.plugin_id = plugin_id
        self.kwargs = dict(kwargs)

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """Return placeholder values until plugin dispatch is implemented."""
        _ = (imt, event, kwargs)

        values = [0.0 for _ in sites]
        sigma_ln = [0.0 for _ in sites]

        return values, sigma_ln


# ---------------------------------------------------------------------------
# Ground-motion context
# ---------------------------------------------------------------------------


@dataclass
class GroundMotionContext:
    """Bind one ScenarioEvent to one ground-motion provider."""

    event: ScenarioEvent
    provider: _BaseProvider

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """Evaluate one IMT at multiple sites."""
        return self.provider.evaluate(
            imt,
            sites,
            self.event,
            **kwargs,
        )

    def evaluate_at_site(
        self,
        imt: str,
        lon: float,
        lat: float,
        elevation_m: float = 0.0,
        **kwargs: Any,
    ) -> Tuple[float, float]:
        """Evaluate one IMT at one WGS84 location."""
        site = WgsPoint(
            longitude=float(lon),
            latitude=float(lat),
            elevation=float(elevation_m),
        )

        values, sigma_ln = self.evaluate(
            imt,
            [site],
            **kwargs,
        )

        return (
            float(values[0]),
            float(sigma_ln[0]),
        )


# ---------------------------------------------------------------------------
# Precomputed IMT helpers
# ---------------------------------------------------------------------------


def _parse_precomputed_imt(imt: str) -> Tuple[str, float | None]:
    """Normalize PGA/PGV and parse arbitrary SA(period) requests."""
    key = str(imt).strip().upper().replace(" ", "")

    if key in {"PGA", "PGV"}:
        return key, None

    match = re.fullmatch(
        r"SA\((\d+(?:\.\d+)?)\)",
        key,
    )

    if match is None:
        raise ValueError(
            f"Unsupported IMT syntax: {imt!r}. "
            "Expected PGA, PGV or SA(period)."
        )

    period = float(match.group(1))

    if not isfinite(period) or period <= 0.0:
        raise ValueError(
            "SA period must be finite and positive."
        )

    return "SA", period


def _finite_float(
    value: Any,
    field_name: str,
) -> float:
    """Convert a JSON numeric value to a finite float."""
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{field_name} must be numeric."
        ) from exc

    if not isfinite(number):
        raise ValueError(
            f"{field_name} must be finite."
        )

    return number

def _validate_imt_event_id(
    event_id: str,
) -> str:
    """Validate an event id used as one IMT JSON filename."""
    if not isinstance(event_id, str):
        raise TypeError(
            "IMT event_id must be a string."
        )

    key = event_id.strip()

    if not key:
        raise ValueError(
            "IMT event_id must not be empty."
        )

    if key in {".", ".."}:
        raise ValueError(
            "Invalid IMT event_id."
        )

    if "/" in key or "\\" in key:
        raise ValueError(
            "IMT event_id must be a file stem, not a path."
        )

    return key


# ---------------------------------------------------------------------------
# ShakeMap helpers
# ---------------------------------------------------------------------------


def _normalize_imt(imt: str) -> str:
    """Normalize supported IMT spelling."""
    key = str(imt).strip().upper().replace(" ", "")

    if key in {"PGA", "PGV", "MMI"}:
        return key

    match = re.fullmatch(
        r"SA\((\d+(?:\.\d+)?)\)",
        key,
    )

    if match is None:
        return key

    period = float(match.group(1))

    if abs(period - 0.3) <= 1e-9:
        return "SA(0.3)"

    if abs(period - 1.0) <= 1e-9:
        return "SA(1.0)"

    if abs(period - 3.0) <= 1e-9:
        return "SA(3.0)"

    return f"SA({period:g})"


def _validate_event_directory_id(
    event_id: str,
) -> str:
    """Validate an event id used as one directory name."""
    if not isinstance(event_id, str):
        raise TypeError(
            "ShakeMap event_id must be a string."
        )

    key = event_id.strip()

    if not key:
        raise ValueError(
            "ShakeMap event_id must not be empty."
        )

    if key in {".", ".."}:
        raise ValueError(
            "Invalid ShakeMap event_id."
        )

    if "/" in key or "\\" in key:
        raise ValueError(
            "ShakeMap event_id must be a directory name, "
            "not a path."
        )

    return key


def _normalize_imt_acceleration_unit(
    unit: str,
) -> str:
    """
    Normalize acceleration units supported by the IMT provider.

    Supported canonical units are ``g`` and ``m/s2``.
    """
    key = str(unit).strip().lower().replace(" ", "")

    aliases = {
        "g": "g",
        "m/s2": "m/s2",
        "m/s^2": "m/s2",
        "m/s²": "m/s2",
        "mps2": "m/s2",
    }

    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            "IMT acceleration_unit must be 'g' or 'm/s2'."
        ) from exc


def _normalize_acceleration_unit(
    unit: str,
) -> str:
    key = str(unit).strip().lower().replace(" ", "")

    aliases = {
        "g": "g",
        "percent_g": "percent_g",
        "percentg": "percent_g",
        "%g": "percent_g",
        "pctg": "percent_g",
    }

    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            "acceleration_unit must be 'g' or 'percent_g'."
        ) from exc


def _normalize_velocity_unit(
    unit: str,
) -> str:
    key = str(unit).strip().lower().replace(" ", "")

    aliases = {
        "cm/s": "cm/s",
        "cm/sec": "cm/s",
        "cms": "cm/s",
        "m/s": "m/s",
        "m/sec": "m/s",
        "ms": "m/s",
    }

    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            "velocity_unit must be 'cm/s' or 'm/s'."
        ) from exc


def _normalize_native_unit(
    unit: str,
) -> str:
    key = str(unit).strip().lower().replace(" ", "")

    aliases = {
        "": "",
        "g": "g",
        "%g": "percent_g",
        "pctg": "percent_g",
        "percentg": "percent_g",
        "percent_g": "percent_g",
        "cm/s": "cm/s",
        "cm/sec": "cm/s",
        "cms": "cm/s",
        "m/s": "m/s",
        "m/sec": "m/s",
        "ms": "m/s",
        "intensity": "dimensionless",
        "mmi": "dimensionless",
        "unitless": "dimensionless",
    }

    return aliases.get(key, key)


def _convert_shakemap_value(
    value: float,
    *,
    quantity: str,
    native_unit: str,
    acceleration_unit: str,
    velocity_unit: str,
) -> float:
    """Convert one native ShakeMap value to configured output units."""
    numeric = float(value)

    if not isfinite(numeric):
        return numeric

    native = _normalize_native_unit(
        native_unit
    )

    if quantity == "acceleration":
        # ShakeMap XML commonly stores PGA/PSA in percent-g.
        if native == "":
            native = "percent_g"

        if native == acceleration_unit:
            return numeric

        if native == "percent_g" and acceleration_unit == "g":
            return numeric / 100.0

        if native == "g" and acceleration_unit == "percent_g":
            return numeric * 100.0

        raise ValueError(
            "Unsupported ShakeMap acceleration-unit conversion: "
            f"{native_unit!r} -> {acceleration_unit!r}."
        )

    if quantity == "velocity":
        # ShakeMap XML commonly stores PGV in cm/s.
        if native == "":
            native = "cm/s"

        if native == velocity_unit:
            return numeric

        if native == "cm/s" and velocity_unit == "m/s":
            return numeric / 100.0

        if native == "m/s" and velocity_unit == "cm/s":
            return numeric * 100.0

        raise ValueError(
            "Unsupported ShakeMap velocity-unit conversion: "
            f"{native_unit!r} -> {velocity_unit!r}."
        )

    if quantity == "intensity":
        return numeric

    raise ValueError(
        f"Unsupported ShakeMap quantity: {quantity!r}."
    )


def _validate_sigma_ln(
    value: float,
) -> float:
    sigma = float(value)

    if not isfinite(sigma) or sigma < 0.0:
        raise ValueError(
            f"Invalid ShakeMap sigma_ln value: {value!r}."
        )

    return sigma


__all__ = [
    "ScenarioEvent",
    "GroundMotionProvider",
    "GroundMotionContext",
]
