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
GMPE registry utilities for ShakeLab.

This module loads a JSON registry (registry.json) shipped with the package and
provides helpers to resolve a GMPE by name or alias.

Registry format
---------------
The JSON file must follow this schema:

{
  "gmpemap": {
    "<CanonicalName>": {
      "pointer": "<module.path:ClassName>",
      "alias": ["<Alias1>", "<Alias2>", ...]
    },
    ...
  }
}

Design notes
------------
- The registry is explicit (no module scanning).
- GMPE classes are imported lazily (only when requested).
- Aliases are resolved to a canonical GMPE name with collision checks.
- The JSON is loaded via importlib.resources to work for installed packages.

Optional override
-----------------
If the environment variable SHAKELAB_GMPE_REGISTRY is set to an absolute
filesystem path, that JSON file is loaded instead of the packaged resource.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from importlib import import_module
from importlib import resources as importlib_resources
from typing import Any, Dict, List, Mapping, Optional, Tuple, Type

from .base import GMPE


_ENV_REGISTRY_PATH = "SHAKELAB_GMPE_REGISTRY"
_DEFAULT_RESOURCE = "registry.json"


@dataclass(frozen=True)
class _RegistryData:
    """
    Normalized registry data.
    """

    gmpemap: Dict[str, str]
    aliases: Dict[str, str]


# Caches
_REGISTRY_CACHE: Optional[_RegistryData] = None
_CLASS_CACHE: Dict[str, Type[GMPE]] = {}


def list_gmpes() -> List[str]:
    """
    Return the list of canonical GMPE names available in the registry.
    """
    reg = _load_registry()
    return sorted(reg.gmpemap.keys())


def resolve_name(name: str) -> str:
    """
    Resolve a name or alias to its canonical GMPE name.

    Parameters
    ----------
    name
        Canonical name or alias.

    Returns
    -------
    str
        Canonical GMPE name.

    Raises
    ------
    KeyError
        If `name` is unknown.
    """
    if not isinstance(name, str) or not name.strip():
        raise KeyError("GMPE name must be a non-empty string.")

    reg = _load_registry()
    if name in reg.gmpemap:
        return name

    canonical = reg.aliases.get(name)
    if canonical is None:
        available = ", ".join(list_gmpes())
        raise KeyError(
            f"Unknown GMPE {name!r}. Available: {available}"
        )
    return canonical


def get_gmpe_class(name: str) -> Type[GMPE]:
    """
    Resolve a GMPE name or alias to its class (lazy import).

    Parameters
    ----------
    name
        Canonical GMPE name or alias.

    Returns
    -------
    Type[GMPE]
        The resolved GMPE class.

    Raises
    ------
    KeyError
        If `name` is unknown.
    ImportError
        If the module/class cannot be imported.
    TypeError
        If the resolved object is not a GMPE subclass.
    """
    canonical = resolve_name(name)

    if canonical in _CLASS_CACHE:
        return _CLASS_CACHE[canonical]

    reg = _load_registry()
    pointer = reg.gmpemap[canonical]
    module_name, class_name = _split_pointer(pointer)

    module = import_module(module_name)
    cls = getattr(module, class_name, None)

    if cls is None:
        raise ImportError(
            f"GMPE class {class_name!r} not found in module "
            f"{module_name!r} (pointer={pointer!r})."
        )

    if not isinstance(cls, type) or not issubclass(cls, GMPE):
        raise TypeError(
            f"Resolved object {module_name}:{class_name} is not a "
            f"GMPE subclass."
        )

    _CLASS_CACHE[canonical] = cls
    return cls


def create_gmpe(name: str, **kwargs: Any) -> GMPE:
    """
    Instantiate a GMPE by canonical name or alias.

    Parameters
    ----------
    name
        Canonical GMPE name or alias.
    **kwargs
        Forwarded to the GMPE constructor.

    Returns
    -------
    GMPE
        GMPE instance.
    """
    cls = get_gmpe_class(name)
    return cls(**kwargs)


def clear_caches() -> None:
    """
    Clear registry and class caches (mainly for tests).
    """
    global _REGISTRY_CACHE
    _REGISTRY_CACHE = None
    _CLASS_CACHE.clear()


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

def _load_registry() -> _RegistryData:
    """
    Load and normalize the registry from JSON (with caching).
    """
    global _REGISTRY_CACHE
    if _REGISTRY_CACHE is not None:
        return _REGISTRY_CACHE

    raw = _read_registry_json()
    reg = _normalize_registry(raw)

    _REGISTRY_CACHE = reg
    return reg


def _read_registry_json() -> Mapping[str, Any]:
    """
    Read the registry JSON from an override path or package resources.
    """
    override = os.environ.get(_ENV_REGISTRY_PATH, "").strip()
    if override:
        with open(override, "r", encoding="utf-8") as f:
            return json.load(f)

    # Load packaged registry.json (same package as this module)
    data = importlib_resources.files(__package__).joinpath(_DEFAULT_RESOURCE)
    with data.open("r", encoding="utf-8") as f:
        return json.load(f)


def _normalize_registry(raw: Mapping[str, Any]) -> _RegistryData:
    """
    Validate and normalize the raw JSON registry into fast lookup maps.

    Returns
    -------
    _RegistryData
        - gmpemap: canonical -> pointer
        - aliases: alias -> canonical
    """
    if not isinstance(raw, Mapping):
        raise TypeError("registry.json must contain a JSON object at top level.")

    gmpemap_raw = raw.get("gmpemap")
    if not isinstance(gmpemap_raw, Mapping):
        raise TypeError("registry.json must contain a 'gmpemap' object.")

    gmpemap: Dict[str, str] = {}
    aliases: Dict[str, str] = {}

    # Collect canonical pointers first
    for canonical, entry in gmpemap_raw.items():
        if not isinstance(canonical, str) or not canonical.strip():
            raise ValueError("Canonical GMPE names must be non-empty strings.")

        if not isinstance(entry, Mapping):
            raise TypeError(
                f"Registry entry for {canonical!r} must be an object."
            )

        pointer = entry.get("pointer")
        if not isinstance(pointer, str) or not pointer.strip():
            raise ValueError(
                f"Registry entry for {canonical!r} must contain a non-empty "
                f"'pointer' string."
            )

        # Validate pointer format early
        _split_pointer(pointer)

        if canonical in gmpemap:
            raise ValueError(f"Duplicate canonical GMPE name: {canonical!r}")

        gmpemap[canonical] = pointer

    # Build alias map
    for canonical, entry in gmpemap_raw.items():
        alias_list = entry.get("alias", [])
        if alias_list is None:
            alias_list = []

        if not isinstance(alias_list, list):
            raise TypeError(
                f"'alias' for {canonical!r} must be a list of strings."
            )

        for alias in alias_list:
            if not isinstance(alias, str) or not alias.strip():
                raise ValueError(
                    f"Invalid alias for {canonical!r}: {alias!r}"
                )

            if alias in gmpemap and alias != canonical:
                raise ValueError(
                    f"Alias {alias!r} collides with a canonical GMPE name."
                )

            prev = aliases.get(alias)
            if prev is not None and prev != canonical:
                raise ValueError(
                    f"Alias {alias!r} is defined for multiple GMPEs: "
                    f"{prev!r}, {canonical!r}"
                )

            aliases[alias] = canonical

    return _RegistryData(gmpemap=gmpemap, aliases=aliases)


def _split_pointer(pointer: str) -> Tuple[str, str]:
    """
    Split and validate a 'module.path:ClassName' pointer.
    """
    if ":" not in pointer:
        raise ValueError(
            f"Invalid GMPE pointer {pointer!r}. Expected 'module:ClassName'."
        )

    module_name, class_name = pointer.split(":", 1)
    module_name = module_name.strip()
    class_name = class_name.strip()

    if not module_name or not class_name:
        raise ValueError(
            f"Invalid GMPE pointer {pointer!r}. Expected 'module:ClassName'."
        )

    return module_name, class_name