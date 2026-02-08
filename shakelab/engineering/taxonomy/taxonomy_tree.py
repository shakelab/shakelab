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
Taxonomy-to-fragility mapping for seismic impact scenarios.

This module defines an explicit taxonomy tree used to map exposure taxonomies
to fragility models. The mapping is fully explicit: every taxonomy string
must be listed in the JSON file (no wildcard, no match rules, no priority).

Design overview
---------------
The implementation is organized in two layers.

1) Fragility references
   - FragilityRef: (model_id, weight) pair. Weight is optional and defaults
     to 1.0.

2) Mapping container
   - TaxonomyTree: mapping taxonomy -> list[FragilityRef]. It supports JSON
     I/O and can "resolve" references into actual FragilityModel objects when
     a FragilityCollection is provided.

JSON database
-------------
The expected JSON database uses a stable top-level structure:
- type: "ShakeLabTaxonomyTree"
- schema_version: "1.0.0"
- metadata: object (at least name, date)
- mappings: object mapping taxonomy -> entry
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, TYPE_CHECKING


if TYPE_CHECKING:
    from .fragility import FragilityCollection, FragilityModel


__all__ = [
    "FragilityRef",
    "TaxonomyTree",
]


_SCHEMA_TYPE = "ShakeLabTaxonomyTree"
_SCHEMA_VERSION = "1.0.0"


def _is_number(x: Any) -> bool:
    """Return True for int/float-like values (excluding bool)."""
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _require_str(dct: Mapping[str, Any], key: str, path: str) -> str:
    """Read a required non-empty string."""
    val = dct.get(key)
    if not isinstance(val, str) or not val.strip():
        raise ValueError(f"{path}.{key} must be a non-empty string.")
    return val.strip()


def _require_dict(dct: Mapping[str, Any], key: str, path: str) -> Dict[str, Any]:
    """Read a required object."""
    val = dct.get(key)
    if not isinstance(val, dict):
        raise ValueError(f"{path}.{key} must be an object.")
    return val


def _require_list(dct: Mapping[str, Any], key: str, path: str) -> List[Any]:
    """Read a required list."""
    val = dct.get(key)
    if not isinstance(val, list):
        raise ValueError(f"{path}.{key} must be an array.")
    return val


@dataclass(frozen=True)
class FragilityRef:
    """
    Reference to a fragility model.

    Parameters
    ----------
    id
        Fragility model identifier.
    weight
        Optional weight for logic-tree mixtures. Must be > 0.
        If omitted, it defaults to 1.0.
    """

    id: str
    weight: float = 1.0

    def __post_init__(self) -> None:
        if not isinstance(self.id, str) or not self.id.strip():
            raise ValueError("id must be a non-empty string.")
        if not _is_number(self.weight) or float(self.weight) <= 0.0:
            raise ValueError("weight must be a number > 0.")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "FragilityRef":
        """Build a FragilityRef from a dictionary."""
        mid = _require_str(data, "id", path)
        w = data.get("weight", 1.0)
        if not _is_number(w) or float(w) <= 0.0:
            raise ValueError(f"{path}.weight must be a number > 0.")
        extra = set(dict(data).keys()) - {"id", "weight"}
        if extra:
            raise ValueError(f"{path} has unknown keys: {sorted(extra)}")
        return cls(id=mid, weight=float(w))

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the reference to a JSON-ready dictionary."""
        out: Dict[str, Any] = {"id": self.id}
        if float(self.weight) != 1.0:
            out["weight"] = float(self.weight)
        return out


@dataclass
class TaxonomyTree:
    """
    Explicit mapping from exposure taxonomies to fragility models.

    Parameters
    ----------
    metadata
        Free-form metadata (name, date, description, ...).
    mappings
        Mapping taxonomy -> list of FragilityRef. Taxonomies must be
        non-empty strings, and each list must be non-empty.
    """

    metadata: Dict[str, Any] = field(default_factory=dict)
    mappings: Dict[str, List[FragilityRef]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.metadata, dict):
            raise ValueError("metadata must be a dict.")
        if not isinstance(self.mappings, dict) or not self.mappings:
            raise ValueError("mappings must be a non-empty dict.")

        for tax, refs in self.mappings.items():
            if not isinstance(tax, str) or not tax.strip():
                raise ValueError("Taxonomy keys must be non-empty strings.")
            if not isinstance(refs, list) or not refs:
                raise ValueError(f"mappings[{tax}] must be a non-empty list.")
            for ref in refs:
                if not isinstance(ref, FragilityRef):
                    raise ValueError("All mapping entries must be FragilityRef.")

    def __getitem__(self, taxonomy: str) -> List[FragilityRef]:
        """Return references for a taxonomy (dict-style access)."""
        return self.mappings[taxonomy]

    def __contains__(self, taxonomy: object) -> bool:
        """Return True if taxonomy exists in the tree."""
        return taxonomy in self.mappings

    def __len__(self) -> int:
        """Return the number of mapped taxonomies."""
        return len(self.mappings)

    def list_taxonomies(self) -> List[str]:
        """Return all mapped taxonomies (sorted)."""
        return sorted(self.mappings.keys())

    def get(self, taxonomy: str) -> List[FragilityRef]:
        """Get references for a taxonomy."""
        return self.mappings[taxonomy]

    def resolve(
        self,
        taxonomy: str,
        fragilities: "FragilityCollection",
    ) -> List[Tuple["FragilityModel", float]]:
        """
        Resolve references into fragility models.

        Parameters
        ----------
        taxonomy
            Exposure taxonomy string.
        fragilities
            Fragility model collection used to resolve model ids.

        Returns
        -------
        list of (FragilityModel, float)
            Resolved models with associated weights.
        """
        refs = self.mappings[taxonomy]
        out: List[Tuple["FragilityModel", float]] = []
        for ref in refs:
            out.append((fragilities.get(ref.id), float(ref.weight)))
        return out

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "TaxonomyTree":
        """Build a taxonomy tree from a JSON-ready dictionary."""
        if not isinstance(data, dict):
            raise ValueError("Root JSON value must be an object.")

        tval = data.get("type")
        if tval != _SCHEMA_TYPE:
            raise ValueError(f"root.type must be '{_SCHEMA_TYPE}'.")

        sval = data.get("schema_version")
        if sval != _SCHEMA_VERSION:
            msg = f"root.schema_version must be '{_SCHEMA_VERSION}'."
            raise ValueError(msg)

        meta = data.get("metadata")
        if not isinstance(meta, dict):
            raise ValueError("root.metadata must be an object.")

        mappings_raw = data.get("mappings")
        if not isinstance(mappings_raw, dict) or not mappings_raw:
            raise ValueError("root.mappings must be a non-empty object.")

        mappings: Dict[str, List[FragilityRef]] = {}
        for tax, entry in mappings_raw.items():
            if not isinstance(tax, str) or not tax.strip():
                raise ValueError("mappings keys must be non-empty strings.")
            if not isinstance(entry, dict):
                raise ValueError(f"mappings.{tax} must be an object.")
            refs_raw = _require_list(entry, "fragility_models", f"mappings.{tax}")
            if not refs_raw:
                raise ValueError(f"mappings.{tax}.fragility_models is empty.")

            refs: List[FragilityRef] = []
            for i, rdata in enumerate(refs_raw):
                if not isinstance(rdata, dict):
                    msg = f"mappings.{tax}.fragility_models[{i}] must be object."
                    raise ValueError(msg)
                ref = FragilityRef.from_dict(
                    rdata,
                    path=f"mappings.{tax}.fragility_models[{i}]",
                )
                refs.append(ref)

            extra = set(entry.keys()) - {"fragility_models", "notes"}
            if extra:
                raise ValueError(f"mappings.{tax} has unknown keys: {sorted(extra)}")

            mappings[tax.strip()] = refs

        return cls(metadata=dict(meta), mappings=mappings)

    @classmethod
    def from_json(cls, path: str) -> "TaxonomyTree":
        """Load a taxonomy tree from a JSON file."""
        with open(path, "r", encoding="utf-8") as fobj:
            data = json.load(fobj)
        return cls.from_dict(data)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the taxonomy tree to a JSON-ready dictionary."""
        mappings: Dict[str, Any] = {}
        for tax in self.list_taxonomies():
            refs = [ref.to_dict() for ref in self.mappings[tax]]
            mappings[tax] = {"fragility_models": refs}
        return {
            "type": _SCHEMA_TYPE,
            "schema_version": _SCHEMA_VERSION,
            "metadata": dict(self.metadata),
            "mappings": mappings,
        }

    def to_json(self, path: str, indent: int = 2) -> None:
        """Write the taxonomy tree to a JSON file."""
        with open(path, "w", encoding="utf-8") as fobj:
            json.dump(self.to_dict(), fobj, indent=indent, ensure_ascii=False)
            fobj.write("\n")
