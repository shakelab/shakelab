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
ShakeLab Exposure data model (ShakeLabExposure schema_version 1.0.0).

This module defines a compact representation of an exposure model used by
ShakeLab for seismic impact scenarios. The model is stored as JSON and
organized in three layers:

1) Exposure (root)
   - type: fixed string "ShakeLabExposure"
   - schema_version: fixed string "1.0.0" for this implementation
   - metadata: free-form dictionary with descriptive information
   - assets: list of Asset objects

2) Asset (georeferenced exposure unit)
   Each asset represents either a single building or an aggregated unit
   (e.g., census zone, grid cell, custom polygon). Assets may include an
   optional geometry and optional site attributes.

3) Typology (construction class within an asset)
   Typologies represent one or more buildings of the same construction
   class within an asset. Only (taxonomy, count) are mandatory.

Validation philosophy
---------------------
Validation is pragmatic and schema-aligned: enforce required keys, basic
type constraints, and key cross-field rules (e.g., geometry consistency
with aggregated flag when geometry is provided). Optional fields are
validated only when present.

Public API
----------
- load_exposure(path): read JSON -> Exposure and validate
- save_exposure(exposure, path): validate (optional) -> write JSON
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


ALLOWED_GEOM_TYPES: set[str] = {"Point", "Polygon"}
EXPECTED_TYPE: str = "ShakeLabExposure"
EXPECTED_SCHEMA_VERSION: str = "1.0.0"

CODE_LEVEL_ENUM: set[str] = {"none", "low", "moderate", "high", "retrofit"}

__all__ = ["Typology", "Asset", "Exposure", "load_exposure", "save_exposure"]


def _drop_none(d: Dict[str, Any]) -> Dict[str, Any]:
    """Return a shallow copy of *d* without keys whose value is None."""
    return {k: v for k, v in d.items() if v is not None}


def _is_number(x: Any) -> bool:
    """Return True if x is int/float (but not bool)."""
    return isinstance(x, (int, float)) and not isinstance(x, bool)


@dataclass
class Typology:
    """
    Construction typology for one asset.

    Mandatory:
    - taxonomy: non-empty string
    - count: integer >= 1

    Other fields are optional and validated only if present (not None).
    """

    taxonomy: str
    count: int

    usage: Optional[str] = None
    building_type: Optional[str] = None
    code_level: Optional[str] = None
    occupants: Optional[Dict[str, float]] = None
    period: Optional[Dict[str, Optional[int]]] = None
    replacement_cost: Optional[float] = None
    stories: Optional[int] = None
    damage_state: Optional[str] = None

    def validate(self, aggregated: bool) -> None:
        """
        Validate a typology with permissive rules.

        Parameters
        ----------
        aggregated
            Whether the parent asset is aggregated. Reserved for future
            cross-field rules.
        """
        if not isinstance(self.taxonomy, str) or not self.taxonomy.strip():
            raise ValueError("Typology.taxonomy must be a non-empty string.")

        if not isinstance(self.count, int) or self.count < 1:
            raise ValueError("Typology.count must be an integer >= 1.")

        if self.usage is not None:
            if not isinstance(self.usage, str):
                raise ValueError("Typology.usage must be a string or None.")

        if self.building_type is not None:
            if not isinstance(self.building_type, str):
                raise ValueError(
                    "Typology.building_type must be a string or None."
                )

        if self.code_level is not None:
            if self.code_level not in CODE_LEVEL_ENUM:
                allowed = ", ".join(sorted(CODE_LEVEL_ENUM))
                raise ValueError(
                    f"Typology.code_level must be one of: {allowed}, or None."
                )

        if self.occupants is not None:
            if not isinstance(self.occupants, dict):
                raise ValueError("Typology.occupants must be a dict or None.")
            for key in ("day", "night"):
                if key not in self.occupants:
                    raise ValueError(
                        f"Typology.occupants must include '{key}'."
                    )
                val = self.occupants[key]
                if not _is_number(val) or val < 0:
                    raise ValueError(
                        f"Typology.occupants.{key} must be a number >= 0."
                    )

        if self.period is not None:
            if not isinstance(self.period, dict):
                raise ValueError("Typology.period must be a dict or None.")

            unknown = set(self.period.keys()) - {"start", "end"}
            if unknown:
                raise ValueError(
                    f"Typology.period contains unknown keys: {sorted(unknown)}."
                )

            start = self.period.get("start")
            end = self.period.get("end")

            if start is not None and not isinstance(start, int):
                raise ValueError("Typology.period.start must be int or None.")
            if end is not None and not isinstance(end, int):
                raise ValueError("Typology.period.end must be int or None.")

        if self.replacement_cost is not None:
            if not _is_number(self.replacement_cost) or self.replacement_cost < 0:
                raise ValueError(
                    "Typology.replacement_cost must be a number >= 0 or None."
                )

        if self.stories is not None:
            if not isinstance(self.stories, int) or self.stories < 1:
                raise ValueError("Typology.stories must be int >= 1 or None.")

        if self.damage_state is not None:
            if not isinstance(self.damage_state, str):
                raise ValueError("Typology.damage_state must be a string or None.")

        _ = aggregated

    def to_dict(self, include_nulls: bool = False) -> Dict[str, Any]:
        """Serialize the typology to a JSON-compatible dict."""
        d = asdict(self)
        return d if include_nulls else _drop_none(d)


@dataclass
class Asset:
    """
    Georeferenced exposure unit (single building or aggregated area).

    Mandatory:
    - id: non-empty string
    - aggregated: boolean
    - reference_location: dict with longitude/latitude numbers
    - typologies: non-empty list[Typology]

    Optional:
    - name: string or None
    - critical: bool or None
    - geometry: GeoJSON-like dict (validated if present)
    - reference_geology: dict or None (validated if present)
    - aggregation_area: number >= 0 or None
    """

    id: str
    aggregated: bool
    reference_location: Dict[str, Any]
    typologies: List[Typology] = field(default_factory=list)

    name: Optional[str] = None
    aggregation_area: Optional[float] = None
    critical: Optional[bool] = None
    geometry: Optional[Dict[str, Any]] = None
    reference_geology: Optional[Dict[str, Any]] = None

    def validate(self) -> None:
        """Validate an asset with permissive, schema-aligned rules."""
        if not isinstance(self.id, str) or not self.id.strip():
            raise ValueError("Asset.id must be a non-empty string.")

        if not isinstance(self.aggregated, bool):
            raise ValueError("Asset.aggregated must be a boolean.")

        if self.name is not None and not isinstance(self.name, str):
            raise ValueError("Asset.name must be a string or None.")

        if self.critical is not None and not isinstance(self.critical, bool):
            raise ValueError("Asset.critical must be a boolean or None.")

        if self.aggregation_area is not None:
            if not _is_number(self.aggregation_area) or self.aggregation_area < 0:
                raise ValueError(
                    "Asset.aggregation_area must be number >= 0 or None."
                )

        if self.geometry is not None:
            if not isinstance(self.geometry, dict):
                raise ValueError("Asset.geometry must be an object or None.")
            gtype = self.geometry.get("type")
            if gtype not in ALLOWED_GEOM_TYPES:
                allowed = ", ".join(sorted(ALLOWED_GEOM_TYPES))
                raise ValueError(
                    f"Asset {self.id}: geometry.type '{gtype}' not allowed. "
                    f"Allowed: {allowed}"
                )
            if "coordinates" not in self.geometry:
                raise ValueError(
                    f"Asset {self.id}: geometry.coordinates is required."
                )

            if self.aggregated and gtype != "Polygon":
                raise ValueError(
                    f"Asset {self.id}: aggregated assets require Polygon."
                )
            if (not self.aggregated) and gtype != "Point":
                raise ValueError(
                    f"Asset {self.id}: single buildings require Point."
                )

        if not isinstance(self.reference_location, dict):
            raise ValueError(
                f"Asset {self.id}: reference_location must be an object."
            )
        for k in ("longitude", "latitude"):
            if k not in self.reference_location:
                raise ValueError(
                    f"Asset {self.id}: reference_location.{k} required."
                )
            if not _is_number(self.reference_location[k]):
                raise ValueError(
                    f"Asset {self.id}: reference_location.{k} must be a number."
                )

        if "elevation" in self.reference_location:
            elev = self.reference_location["elevation"]
            if elev is not None and not _is_number(elev):
                raise ValueError(f"Asset {self.id}: elevation must be a number.")

        if self.reference_geology is not None:
            if not isinstance(self.reference_geology, dict):
                raise ValueError(
                    f"Asset {self.id}: reference_geology must be an object or "
                    f"None."
                )
            if "vs30" in self.reference_geology:
                vs30 = self.reference_geology.get("vs30")
                if vs30 is not None and (not _is_number(vs30) or vs30 < 0):
                    raise ValueError(
                        f"Asset {self.id}: vs30 must be >= 0 or None."
                    )

        if not self.typologies:
            raise ValueError(f"Asset {self.id}: at least one typology is required.")

        for t in self.typologies:
            t.validate(aggregated=self.aggregated)

    def to_dict(self, include_nulls: bool = False) -> Dict[str, Any]:
        """Serialize the asset to a JSON-compatible dict."""
        d: Dict[str, Any] = {
            "id": self.id,
            "name": self.name,
            "aggregated": self.aggregated,
            "aggregation_area": self.aggregation_area,
            "critical": self.critical,
            "geometry": self.geometry,
            "reference_location": self.reference_location,
            "reference_geology": self.reference_geology,
            "typologies": [t.to_dict(include_nulls=include_nulls)
                           for t in self.typologies],
        }
        return d if include_nulls else _drop_none(d)


@dataclass
class ExposureModel:
    """
    Root container for a ShakeLabExposure model (schema_version 1.0.0).

    Mandatory:
    - type
    - schema_version
    - metadata (non-empty dict)
    - assets (non-empty list[Asset])
    """

    type: str
    schema_version: str
    metadata: Dict[str, Any]
    assets: List[Asset]

    def validate(self) -> None:
        """Validate the root exposure container."""
        if self.type != EXPECTED_TYPE:
            raise ValueError(
                f"Exposure.type must be '{EXPECTED_TYPE}', got '{self.type}'."
            )
        if self.schema_version != EXPECTED_SCHEMA_VERSION:
            raise ValueError(
                f"Exposure.schema_version must be '{EXPECTED_SCHEMA_VERSION}'."
            )
        if not isinstance(self.metadata, dict) or not self.metadata:
            raise ValueError("Exposure.metadata must be a non-empty dict.")
        if not self.assets:
            raise ValueError("Exposure.assets must not be empty.")
        for asset in self.assets:
            asset.validate()

    def to_dict(self, include_nulls: bool = False) -> Dict[str, Any]:
        """Serialize the exposure container to a JSON-compatible dict."""
        d = {
            "type": self.type,
            "schema_version": self.schema_version,
            "metadata": self.metadata,
            "assets": [a.to_dict(include_nulls=include_nulls)
                       for a in self.assets],
        }
        return d if include_nulls else _drop_none(d)

    def list_taxonomies(self) -> List[str]:
        """
        Return the list of unique taxonomies present in the exposure.
    
        Taxonomies are derived from the typologies attached to each asset.
        The list is returned in sorted order.
    
        Returns
        -------
        list of str
            Sorted list of unique typology taxonomies.
        """
        tax = {
            typ.taxonomy
            for asset in self.assets
            for typ in asset.typologies
            if isinstance(typ.taxonomy, str) and typ.taxonomy.strip()
        }
        return sorted(tax)

    @classmethod
    def from_json(
        cls,
        json_path: Union[str, Path],
        validate: bool = True,
    ) -> "Exposure":
        """
        Load an exposure dataset from a JSON file.
    
        This is a convenience wrapper around `load_exposure`.
    
        Parameters
        ----------
        json_path
            Path to the exposure JSON file.
        validate
            If True, run full model validation.
    
        Returns
        -------
        Exposure
            Parsed Exposure instance.
        """
        return load_exposure(json_path, validate=validate)

    def to_json(
        self,
        json_path: Union[str, Path],
        validate: bool = True,
        *,
        include_nulls: bool = False,
        indent: int = 2,
        sort_keys: bool = False,
        ensure_ascii: bool = False,
    ) -> None:
        """
        Write the exposure dataset to a JSON file.
    
        This is a convenience wrapper around `save_exposure`.
    
        Parameters
        ----------
        json_path
            Output path for the JSON file.
        validate
            If True, run validation before writing.
        include_nulls
            If True, include keys with null values.
        indent
            Indentation level for JSON formatting.
        sort_keys
            Sort JSON keys alphabetically.
        ensure_ascii
            Escape non-ASCII characters if True.
        """
        save_exposure(
            self,
            json_path,
            validate=validate,
            include_nulls=include_nulls,
            indent=indent,
            sort_keys=sort_keys,
            ensure_ascii=ensure_ascii,
        )


def load_exposure(
    json_path: Union[str, Path],
    validate: bool = True,
) -> ExposureModel:
    """
    Load a ShakeLabExposure JSON file (schema_version 1.0.0) into an
    Exposure object.

    Parameters
    ----------
    json_path
        Path to the JSON file.
    validate
        If True, run the datamodel validation (ExposureModel.validate()).

    Returns
    -------
    Exposure
        Parsed Exposure instance.
    """
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    for key in ("type", "schema_version", "metadata", "assets"):
        if key not in data:
            raise ValueError(f"Missing required field '{key}' in file.")

    assets_in = data.get("assets")
    if not isinstance(assets_in, list) or not assets_in:
        raise ValueError("Field 'assets' must be a non-empty list.")

    assets: List[Asset] = []
    for a in assets_in:
        if not isinstance(a, dict):
            raise ValueError("Each asset must be an object (dict).")

        for key in ("id", "aggregated", "reference_location", "typologies"):
            if key not in a:
                raise ValueError(
                    f"Missing required field 'assets[].{key}' in file."
                )

        typs_in = a.get("typologies")
        if not isinstance(typs_in, list) or not typs_in:
            raise ValueError("Field 'assets[].typologies' must be non-empty.")

        typologies: List[Typology] = []
        for t in typs_in:
            if not isinstance(t, dict):
                raise ValueError("Each typology must be an object (dict).")

            for key in ("taxonomy", "count"):
                if key not in t:
                    raise ValueError(
                        f"Missing required field 'typologies[].{key}'."
                    )

            typology = Typology(
                taxonomy=t.get("taxonomy"),
                count=int(t.get("count")),
                usage=t.get("usage"),
                building_type=t.get("building_type"),
                code_level=t.get("code_level"),
                occupants=t.get("occupants"),
                period=t.get("period"),
                replacement_cost=t.get("replacement_cost"),
                stories=t.get("stories"),
                damage_state=t.get("damage_state"),
            )
            typologies.append(typology)

        asset = Asset(
            id=a.get("id"),
            name=a.get("name"),
            aggregated=bool(a.get("aggregated")),
            aggregation_area=a.get("aggregation_area"),
            critical=a.get("critical"),
            geometry=a.get("geometry"),
            reference_location=a.get("reference_location"),
            reference_geology=a.get("reference_geology"),
            typologies=typologies,
        )
        assets.append(asset)

    exposure = ExposureModel(
        type=data.get("type"),
        schema_version=data.get("schema_version"),
        metadata=data.get("metadata"),
        assets=assets,
    )

    if validate:
        exposure.validate()

    return exposure


def save_exposure(
    exposure: ExposureModel,
    json_path: Union[str, Path],
    validate: bool = True,
    *,
    include_nulls: bool = False,
    indent: int = 2,
    sort_keys: bool = False,
    ensure_ascii: bool = False,
) -> None:
    """
    Save a ShakeLabExposure object to a JSON file.

    Parameters
    ----------
    exposure
        Exposure instance to serialize.
    json_path
        Output path for the JSON file.
    validate
        If True, call `exposure.validate()` before writing.
    include_nulls
        If True, preserve keys with null values in the output JSON.
    indent
        Indentation level passed to `json.dump`.
    sort_keys
        Sort JSON keys alphabetically.
    ensure_ascii
        If True, escape non-ASCII characters.
    """
    if validate:
        exposure.validate()

    path = Path(json_path)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)

    data = exposure.to_dict(include_nulls=include_nulls)

    with path.open("w", encoding="utf-8") as f:
        json.dump(
            data,
            f,
            ensure_ascii=ensure_ascii,
            indent=indent,
            sort_keys=sort_keys,
        )
        f.write("\n")
