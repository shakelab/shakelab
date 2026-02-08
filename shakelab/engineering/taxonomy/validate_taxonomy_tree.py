#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ShakeLabTaxonomyTree validator (schema_version 1.0.0).

This validator mirrors the pragmatic approach used for ShakeLabExposure:
it enforces a stable minimum structure and rejects unknown keys, while
keeping runtime policies (e.g., weight normalization, fallback behavior)
outside the JSON file.

Key rules (root level)
----------------------
Required:
- type == "ShakeLabTaxonomyTree"
- schema_version == "1.0.0"
- metadata (object with non-empty name and date)
- mappings (non-empty object)

Key rules (mapping entry)
-------------------------
Required:
- fragility_models (non-empty array)

Each fragility model reference:
- id (non-empty string)
- weight (optional, number > 0)

Unknown keys are not allowed (schema-aligned).
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from typing import Any, Dict, List


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
    if data.get("type") != "ShakeLabTaxonomyTree":
        _add(issues, "root.type", "must be 'ShakeLabTaxonomyTree'")

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

    allowed = {"name", "date", "description", "version", "source", "license"}
    for k in meta.keys():
        if k not in allowed:
            # metadata in exposure schema allows additionalProperties=true.
            # Here we keep it permissive as well: do NOT error.
            pass

    allowed_root = {"type", "schema_version", "metadata", "mappings"}
    unknown = set(data.keys()) - allowed_root
    if unknown:
        _add(
            issues,
            "root",
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_fragility_ref(
    ref: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """Validate a single fragility model reference."""
    if not isinstance(ref, dict):
        _add(issues, path, "must be an object")
        return

    model_id = ref.get("id")
    if not isinstance(model_id, str) or not model_id.strip():
        _add(issues, path + ".id", "is required and must be a string")

    if "weight" in ref:
        w = ref.get("weight")
        if w is None or (not _is_number(w)) or w <= 0:
            _add(issues, path + ".weight", "must be a number > 0")

    allowed = {"id", "weight"}
    unknown = set(ref.keys()) - allowed
    if unknown:
        _add(
            issues,
            path,
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_mapping_entry(
    entry: Any,
    path: str,
    issues: List[ValidationIssue],
) -> None:
    """Validate a mapping entry for one taxonomy."""
    if not isinstance(entry, dict):
        _add(issues, path, "must be an object")
        return

    models = entry.get("fragility_models")
    if not isinstance(models, list) or not models:
        _add(issues, path + ".fragility_models", "must be a non-empty array")
    else:
        for i, ref in enumerate(models):
            validate_fragility_ref(ref, f"{path}.fragility_models[{i}]", issues)

    if "notes" in entry:
        notes = entry.get("notes")
        if not isinstance(notes, str) or not notes.strip():
            _add(issues, path + ".notes", "must be a non-empty string")

    allowed = {"fragility_models", "notes"}
    unknown = set(entry.keys()) - allowed
    if unknown:
        _add(
            issues,
            path,
            f"unknown keys not allowed: {sorted(unknown)}",
        )


def validate_taxonomy_tree(data: Dict[str, Any]) -> List[ValidationIssue]:
    """Validate a ShakeLabTaxonomyTree object and return a list of issues."""
    issues: List[ValidationIssue] = []

    validate_root(data, issues)

    mappings = data.get("mappings")
    if not isinstance(mappings, dict) or not mappings:
        _add(issues, "root.mappings", "must be a non-empty object")
        return issues

    for taxonomy, entry in mappings.items():
        if not isinstance(taxonomy, str) or not taxonomy.strip():
            _add(issues, "mappings", "taxonomy keys must be non-empty strings")
            continue
        validate_mapping_entry(entry, f"mappings.{taxonomy}", issues)

    return issues


def main() -> None:
    """CLI entry point."""
    if len(sys.argv) != 2:
        print("Usage: validate_taxonomy_tree.py <taxonomy_tree.json>")
        sys.exit(2)

    try:
        data = load_json(sys.argv[1])
    except Exception as exc:
        print(f"ERROR: cannot read JSON: {exc}")
        sys.exit(1)

    issues = validate_taxonomy_tree(data)
    if issues:
        print("Validation FAILED:")
        for issue in issues:
            print(f" - {issue}")
        sys.exit(1)

    print("Validation OK.")


if __name__ == "__main__":
    main()
