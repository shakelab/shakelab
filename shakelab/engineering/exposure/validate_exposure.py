#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ShakeLabExposure validator (schema_version 1.0.0).

This validator is intentionally pragmatic: it enforces the minimum
structural constraints needed to ensure downstream components can rely
on a stable shape, while allowing partially-populated exposure models.

Key rules (asset level)
-----------------------
Required:
- id (non-empty string)
- aggregated (boolean)
- reference_location (object with numeric longitude/latitude)
- typologies (non-empty array)

Note: geometry is optional; when provided it must be a GeoJSON-like Point or
Polygon (and MultiPolygon in the future)

Key rules (typology level)
--------------------------
Required:
- taxonomy (non-empty string)
- count (int >= 1)

Other fields are optional; when present they are type-checked and validated.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from typing import Any, Dict, List


CODE_LEVEL = {"none", "low", "moderate", "high", "retrofit"}
GEOM_TYPES = {"Point", "Polygon"}


class ValidationError(Exception):
    """Raised when validation fails (not used by default CLI flow)."""


@dataclass(frozen=True)
class ValidationIssue:
    """
    Single validation issue.

    Attributes
    ----------
    path
        JSON-like path to the failing field.
    message
        Human-readable description of the issue.
    """

    path: str
    message: str

    def __str__(self) -> str:
        return f"{self.path}: {self.message}"


def load_json(path: str) -> Dict[str, Any]:
    """Load a JSON file as a dictionary."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Root JSON value must be an object.")
    return data


def _is_number(x: Any) -> bool:
    """Return True if x is int/float (but not bool)."""
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _add(issues: List[ValidationIssue], path: str, msg: str) -> None:
    """Append a validation issue."""
    issues.append(ValidationIssue(path=path, message=msg))


def _is_coord2(x: Any) -> bool:
    """Check a [lon, lat] coordinate pair."""
    return (
        isinstance(x, list)
        and len(x) == 2
        and _is_number(x[0])
        and _is_number(x[1])
    )


def _is_ring(x: Any) -> bool:
    """
    Check a polygon ring: list of coord2, length >= 4, closed.

    Minimal check: first coordinate equals last coordinate.
    """
    if not isinstance(x, list) or len(x) < 4:
        return False
    if not all(_is_coord2(c) for c in x):
        return False
    return x[0] == x[-1]


def validate_root(data: Dict[str, Any], issues: List[ValidationIssue]) -> None:
    """Validate top-level keys and minimal metadata."""
    if data.get("type") != "ShakeLabExposure":
        _add(issues, "root.type", "must be 'ShakeLabExposure'")

    if data.get("schema_version") != "1.0.0":
        _add(issues, "root.schema_version", "must be '1.0.0'")

    meta = data.get("metadata")
    if not isinstance(meta, dict):
        _add(issues, "root.metadata", "must be an object")
        return

    name = meta.get("name")
    if not isinstance(name, str) or not name.strip():
        _add(issues, "metadata.name", "is required and must be a string")

    date = meta.get("date")
    if not isinstance(date, str) or not date.strip():
        _add(issues, "metadata.date", "is required and must be a string")


def validate_geometry(
    geom: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """
    Validate optional geometry.

    When present, geometry must be a GeoJSON-like Point or Polygon.
    This function does not enforce CRS.
    """
    if geom is None:
        return

    if not isinstance(geom, dict):
        _add(issues, path, "must be an object or null")
        return

    gtype = geom.get("type")
    coords = geom.get("coordinates")

    if gtype not in GEOM_TYPES:
        _add(issues, path + ".type", "must be Point or Polygon")
        return

    if "coordinates" not in geom:
        _add(issues, path + ".coordinates", "is required")
        return

    if gtype == "Point":
        if not _is_coord2(coords):
            _add(issues, path + ".coordinates", "must be [lon, lat]")
        return

    if not isinstance(coords, list) or len(coords) < 1:
        _add(issues, path + ".coordinates", "must be a list of rings")
        return

    if not all(_is_ring(r) for r in coords):
        _add(issues, path + ".coordinates", "rings must be closed and valid")


def validate_reference_location(
    loc: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """Validate reference_location (required)."""
    if not isinstance(loc, dict):
        _add(issues, path, "must be an object")
        return

    for k in ("longitude", "latitude"):
        if k not in loc or not _is_number(loc[k]):
            _add(issues, f"{path}.{k}", "must be a number")

    if "elevation" in loc and loc["elevation"] is not None:
        if not _is_number(loc["elevation"]):
            _add(issues, path + ".elevation", "must be a number or null")


def validate_reference_geology(
    geo: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """
    Validate optional reference_geology.

    If missing or null: OK.
    If present: must be an object. Only vs30 is checked (>= 0) when given.
    """
    if geo is None:
        return

    if not isinstance(geo, dict):
        _add(issues, path, "must be an object or null")
        return

    if "vs30" in geo and geo["vs30"] is not None:
        if not _is_number(geo["vs30"]) or geo["vs30"] < 0:
            _add(issues, path + ".vs30", "must be number >= 0 or null")


def validate_occupants(
    occ: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """Validate occupants object if present."""
    if occ is None:
        return

    if not isinstance(occ, dict):
        _add(issues, path, "must be an object or null")
        return

    for k in ("day", "night"):
        if k not in occ or not _is_number(occ[k]) or occ[k] < 0:
            _add(issues, f"{path}.{k}", "must be number >= 0")


def validate_period(
    per: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """Validate period object if present."""
    if per is None:
        return

    if not isinstance(per, dict):
        _add(issues, path, "must be an object or null")
        return

    if "start" in per and per["start"] is not None:
        if not isinstance(per["start"], int):
            _add(issues, path + ".start", "must be integer or null")

    if "end" in per and per["end"] is not None:
        if not isinstance(per["end"], int):
            _add(issues, path + ".end", "must be integer or null")

    unknown = set(per.keys()) - {"start", "end"}
    if unknown:
        _add(
            issues,
            path,
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_typology(
    t: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """
    Validate a typology.

    Required:
    - taxonomy (non-empty string)
    - count (int >= 1)

    Optional fields are validated when present.
    Unknown keys are not allowed (schema-aligned).
    """
    if not isinstance(t, dict):
        _add(issues, path, "must be an object")
        return

    tax = t.get("taxonomy")
    if not isinstance(tax, str) or not tax.strip():
        _add(issues, path + ".taxonomy", "is required and must be a string")

    count = t.get("count")
    if not isinstance(count, int) or count < 1:
        _add(issues, path + ".count", "is required and must be int >= 1")

    for k in ("usage", "building_type"):
        if k in t:
            val = t.get(k)
            if val is not None and (not isinstance(val, str) or not val.strip()):
                _add(
                    issues,
                    path + f".{k}",
                    "must be a non-empty string or null",
                )

    if "code_level" in t:
        cl = t.get("code_level")
        if cl is not None and cl not in CODE_LEVEL:
            _add(
                issues,
                path + ".code_level",
                f"must be one of {sorted(CODE_LEVEL)} or null",
            )

    if "occupants" in t:
        validate_occupants(t.get("occupants"), path + ".occupants", issues)

    if "period" in t:
        validate_period(t.get("period"), path + ".period", issues)

    if "replacement_cost" in t:
        rc = t.get("replacement_cost")
        if rc is not None and (not _is_number(rc) or rc < 0):
            _add(
                issues,
                path + ".replacement_cost",
                "must be number >= 0 or null",
            )

    if "stories" in t:
        st = t.get("stories")
        if st is not None and (not isinstance(st, int) or st < 1):
            _add(issues, path + ".stories", "must be int >= 1 or null")

    if "damage_state" in t:
        ds = t.get("damage_state")
        if ds is not None and (not isinstance(ds, str) or not ds.strip()):
            _add(
                issues,
                path + ".damage_state",
                "must be a non-empty string or null",
            )

    allowed = {
        "taxonomy",
        "count",
        "usage",
        "building_type",
        "code_level",
        "occupants",
        "period",
        "replacement_cost",
        "stories",
        "damage_state",
    }
    unknown = set(t.keys()) - allowed
    if unknown:
        _add(
            issues,
            path,
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_asset(a: Any, path: str, issues: List[ValidationIssue]) -> None:
    """
    Validate an asset.

    Required:
    - id (non-empty string)
    - aggregated (boolean)
    - reference_location (lon/lat numbers)
    - typologies (non-empty array)

    Optional:
    - name (string or null; empty string allowed)
    - critical (boolean or null)
    - reference_geology (object or null; {} allowed)
    - geometry (object or null; validated if present)
    - aggregation_area (number >= 0 or null)
    """
    if not isinstance(a, dict):
        _add(issues, path, "must be an object")
        return

    asset_id = a.get("id")
    if not isinstance(asset_id, str) or not asset_id.strip():
        _add(issues, f"{path}.id", "is required and must be a string")

    aggregated = a.get("aggregated")
    if not isinstance(aggregated, bool):
        _add(issues, path + ".aggregated", "is required and must be boolean")

    if "name" in a:
        nm = a.get("name")
        if nm is not None and not isinstance(nm, str):
            _add(issues, path + ".name", "must be a string or null")

    if "critical" in a:
        cr = a.get("critical")
        if cr is not None and not isinstance(cr, bool):
            _add(issues, path + ".critical", "must be boolean or null")

    if "aggregation_area" in a:
        aa = a.get("aggregation_area")
        if aa is not None and (not _is_number(aa) or aa < 0):
            _add(
                issues,
                path + ".aggregation_area",
                "must be number >= 0 or null",
            )

    if "geometry" in a:
        validate_geometry(a.get("geometry"), path + ".geometry", issues)
        geom = a.get("geometry")
        if isinstance(geom, dict) and isinstance(aggregated, bool):
            gtype = geom.get("type")
            if aggregated and gtype != "Polygon":
                _add(
                    issues,
                    path + ".geometry.type",
                    "must be Polygon when aggregated=true",
                )
            if (not aggregated) and gtype != "Point":
                _add(
                    issues,
                    path + ".geometry.type",
                    "must be Point when aggregated=false",
                )

    validate_reference_location(
        a.get("reference_location"),
        path + ".reference_location",
        issues,
    )

    if "reference_geology" in a:
        validate_reference_geology(
            a.get("reference_geology"),
            path + ".reference_geology",
            issues,
        )

    typs = a.get("typologies")
    if not isinstance(typs, list) or not typs:
        _add(issues, path + ".typologies", "must be a non-empty array")
        return

    for i, t in enumerate(typs):
        validate_typology(t, f"{path}.typologies[{i}]", issues)

    allowed = {
        "id",
        "name",
        "aggregated",
        "aggregation_area",
        "critical",
        "geometry",
        "reference_location",
        "reference_geology",
        "typologies",
    }
    unknown = set(a.keys()) - allowed
    if unknown:
        _add(
            issues,
            path,
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_exposure(data: Dict[str, Any]) -> List[ValidationIssue]:
    """Validate a ShakeLabExposure object and return a list of issues."""
    issues: List[ValidationIssue] = []

    validate_root(data, issues)

    assets = data.get("assets")
    if not isinstance(assets, list) or not assets:
        _add(issues, "root.assets", "must be a non-empty array")
        return issues

    for i, a in enumerate(assets):
        validate_asset(a, f"assets[{i}]", issues)

    return issues


def main() -> None:
    """CLI entry point."""
    if len(sys.argv) != 2:
        print("Usage: validate_exposure.py <exposure.json>")
        sys.exit(2)

    try:
        data = load_json(sys.argv[1])
    except Exception as exc:
        print(f"ERROR: cannot read JSON: {exc}")
        sys.exit(1)

    issues = validate_exposure(data)
    if issues:
        print("Validation FAILED:")
        for issue in issues:
            print(f" - {issue}")
        sys.exit(1)

    print("Validation OK.")


if __name__ == "__main__":
    main()
