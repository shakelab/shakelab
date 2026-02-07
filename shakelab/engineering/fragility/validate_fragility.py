#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ShakeLabFragility validator (schema_version 1.0.0).

This validator is intentionally pragmatic: it enforces the minimum
structural constraints needed to ensure downstream components can rely
on a stable shape, while allowing additional metadata fields.

Key rules (root)
----------------
Required:
- type == "ShakeLabFragility"
- schema_version == "1.0.0"
- metadata.name (non-empty string)
- metadata.date (non-empty string; ISO date recommended)
- models (non-empty array)

Key rules (model)
-----------------
Required:
- id, taxonomy, imt, model_type, damage_scale, im_bounds
- damage_scale.levels must be a non-empty list of unique strings
- model_type:
  - "lognormal_continuous" requires parameters[level] = {theta, beta}
  - "discrete" requires tables[level] = {im, poe, [log_im]}

Additional pragmatic checks
--------------------------
- model ids must be unique within the collection
- parameters/tables must cover exactly damage_scale.levels
- lognormal theta > 0, beta > 0
- discrete tables:
  - len(im) == len(poe) >= 2
  - im strictly increasing
  - poe in [0, 1]
  - if log_im=True, im must be > 0
- im_bounds: min >= 0, max > min
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set


MODEL_TYPES = {"lognormal_continuous", "discrete"}


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


def _add(issues: List[ValidationIssue], path: str, msg: str) -> None:
    """Append a validation issue."""
    issues.append(ValidationIssue(path=path, message=msg))


def _is_number(x: Any) -> bool:
    """Return True if x is int/float (but not bool)."""
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def validate_root(data: Dict[str, Any], issues: List[ValidationIssue]) -> None:
    """Validate top-level keys and minimal metadata."""
    if data.get("type") != "ShakeLabFragility":
        _add(issues, "root.type", "must be 'ShakeLabFragility'")

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


def _read_levels(model: Dict[str, Any], issues: List[ValidationIssue],
                 path: str) -> Optional[List[str]]:
    """Read and validate damage_scale.levels."""
    ds = model.get("damage_scale")
    if not isinstance(ds, dict):
        _add(issues, path + ".damage_scale", "must be an object")
        return None

    ds_id = ds.get("id")
    if not isinstance(ds_id, str) or not ds_id.strip():
        _add(issues, path + ".damage_scale.id", "must be a non-empty string")

    levels = ds.get("levels")
    if not isinstance(levels, list) or not levels:
        _add(issues, path + ".damage_scale.levels", "must be a non-empty array")
        return None

    out: List[str] = []
    for i, lv in enumerate(levels):
        if not isinstance(lv, str) or not lv.strip():
            _add(
                issues,
                f"{path}.damage_scale.levels[{i}]",
                "must be a non-empty string",
            )
        else:
            out.append(lv)

    if len(set(out)) != len(out):
        _add(issues, path + ".damage_scale.levels", "must not contain duplicates")

    return out if out else None


def _validate_im_bounds(model: Dict[str, Any], issues: List[ValidationIssue],
                        path: str) -> None:
    """Validate im_bounds {min,max}."""
    bnd = model.get("im_bounds")
    if not isinstance(bnd, dict):
        _add(issues, path + ".im_bounds", "must be an object")
        return

    vmin = bnd.get("min")
    vmax = bnd.get("max")

    if not _is_number(vmin) or vmin < 0:
        _add(issues, path + ".im_bounds.min", "must be number >= 0")

    if not _is_number(vmax) or vmax <= 0:
        _add(issues, path + ".im_bounds.max", "must be number > 0")

    if _is_number(vmin) and _is_number(vmax) and vmax <= vmin:
        _add(issues, path + ".im_bounds", "max must be > min")


def _validate_required_str(model: Dict[str, Any], key: str,
                           issues: List[ValidationIssue], path: str) -> None:
    """Validate that a required field is a non-empty string."""
    val = model.get(key)
    if not isinstance(val, str) or not val.strip():
        _add(issues, f"{path}.{key}", "is required and must be a string")


def _validate_lognormal_params(
    params: Any,
    levels: List[str],
    issues: List[ValidationIssue],
    path: str,
) -> None:
    """Validate parameters for lognormal_continuous model."""
    if not isinstance(params, dict):
        _add(issues, path + ".parameters", "must be an object")
        return

    keys = set(params.keys())
    expected = set(levels)
    if keys != expected:
        missing = sorted(expected - keys)
        extra = sorted(keys - expected)
        if missing:
            _add(issues, path + ".parameters", f"missing levels: {missing}")
        if extra:
            _add(issues, path + ".parameters", f"unknown levels: {extra}")

    for lv in levels:
        p = params.get(lv)
        lv_path = path + f".parameters.{lv}"
        if not isinstance(p, dict):
            _add(issues, lv_path, "must be an object with theta/beta")
            continue

        theta = p.get("theta")
        beta = p.get("beta")

        if not _is_number(theta) or theta <= 0:
            _add(issues, lv_path + ".theta", "must be number > 0")

        if not _is_number(beta) or beta <= 0:
            _add(issues, lv_path + ".beta", "must be number > 0")


def _is_strictly_increasing(x: List[float]) -> bool:
    """Check strict monotonic increase."""
    for i in range(1, len(x)):
        if not (x[i] > x[i - 1]):
            return False
    return True


def _validate_discrete_tables(
    tables: Any,
    levels: List[str],
    issues: List[ValidationIssue],
    path: str,
) -> None:
    """Validate tables for discrete model."""
    if not isinstance(tables, dict):
        _add(issues, path + ".tables", "must be an object")
        return

    keys = set(tables.keys())
    expected = set(levels)
    if keys != expected:
        missing = sorted(expected - keys)
        extra = sorted(keys - expected)
        if missing:
            _add(issues, path + ".tables", f"missing levels: {missing}")
        if extra:
            _add(issues, path + ".tables", f"unknown levels: {extra}")

    for lv in levels:
        t = tables.get(lv)
        lv_path = path + f".tables.{lv}"
        if not isinstance(t, dict):
            _add(issues, lv_path, "must be an object with im/poe")
            continue

        im = t.get("im")
        poe = t.get("poe")
        log_im = t.get("log_im", False)

        if not isinstance(log_im, bool):
            _add(issues, lv_path + ".log_im", "must be boolean")
            log_im = False

        if not isinstance(im, list) or len(im) < 2:
            _add(issues, lv_path + ".im", "must be an array (len >= 2)")
            continue

        if not isinstance(poe, list) or len(poe) < 2:
            _add(issues, lv_path + ".poe", "must be an array (len >= 2)")
            continue

        if len(im) != len(poe):
            _add(issues, lv_path, "im and poe must have same length")
            continue

        im_vals: List[float] = []
        for i, v in enumerate(im):
            if not _is_number(v):
                _add(issues, f"{lv_path}.im[{i}]", "must be a number")
            else:
                im_vals.append(float(v))

        for i, v in enumerate(poe):
            if not _is_number(v) or v < 0 or v > 1:
                _add(
                    issues,
                    f"{lv_path}.poe[{i}]",
                    "must be a number in [0, 1]",
                )

        if len(im_vals) == len(im) and not _is_strictly_increasing(im_vals):
            _add(issues, lv_path + ".im", "must be strictly increasing")

        if log_im and any(v <= 0 for v in im_vals):
            _add(issues, lv_path + ".im", "must be > 0 when log_im=true")


def validate_model(model: Any, path: str, issues: List[ValidationIssue]) -> None:
    """Validate a single model entry."""
    if not isinstance(model, dict):
        _add(issues, path, "must be an object")
        return

    _validate_required_str(model, "id", issues, path)
    _validate_required_str(model, "taxonomy", issues, path)
    _validate_required_str(model, "imt", issues, path)

    mtype = model.get("model_type")
    if mtype not in MODEL_TYPES:
        _add(
            issues,
            path + ".model_type",
            f"must be one of {sorted(MODEL_TYPES)}",
        )

    levels = _read_levels(model, issues, path)
    _validate_im_bounds(model, issues, path)

    if mtype == "lognormal_continuous":
        _validate_lognormal_params(model.get("parameters"), levels or [],
                                   issues, path)

    if mtype == "discrete":
        _validate_discrete_tables(model.get("tables"), levels or [],
                                  issues, path)


def validate_fragility(data: Dict[str, Any]) -> List[ValidationIssue]:
    """Validate a ShakeLabFragility object and return a list of issues."""
    issues: List[ValidationIssue] = []

    validate_root(data, issues)

    models = data.get("models")
    if not isinstance(models, list) or not models:
        _add(issues, "root.models", "must be a non-empty array")
        return issues

    seen: Set[str] = set()
    for i, m in enumerate(models):
        path = f"models[{i}]"
        validate_model(m, path, issues)
        if isinstance(m, dict):
            mid = m.get("id")
            if isinstance(mid, str) and mid.strip():
                if mid in seen:
                    _add(issues, path + ".id", "duplicate model id")
                seen.add(mid)

    return issues


def main() -> None:
    """CLI entry point."""
    if len(sys.argv) != 2:
        print("Usage: validate_fragility.py <fragility.json>")
        sys.exit(2)

    try:
        data = load_json(sys.argv[1])
    except Exception as exc:
        print(f"ERROR: cannot read JSON: {exc}")
        sys.exit(1)

    issues = validate_fragility(data)
    if issues:
        print("Validation FAILED:")
        for issue in issues:
            print(f" - {issue}")
        sys.exit(1)

    print("Validation OK.")


if __name__ == "__main__":
    main()
