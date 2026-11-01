# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
Architecture tests for :mod:`shakelab.signals`.

The modules of the package must form a hierarchy without import cycles
(``docs/development.md``, Section 3): every module imports only lower
modules when it is loaded. Imports placed inside functions are allowed
for convenience delegation to higher modules (e.g.
``StreamCollection.read`` calling :mod:`shakelab.signals.io`) and are
not considered here.
"""

from __future__ import annotations

import ast
import importlib
import warnings
from pathlib import Path

import pytest

PACKAGE = "shakelab.signals"
ROOT = Path(importlib.import_module(PACKAGE).__file__).parent

# Modules that cannot be imported at present (they import modules
# removed from the repository); they will be revised separately.
BROKEN = {
    f"{PACKAGE}.metadata",
    f"{PACKAGE}.seedlink",
}


def _module_name(path: Path) -> str:
    parts = list(path.relative_to(ROOT.parent.parent).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _modules() -> dict[str, Path]:
    """Python modules of the package, excluding data folders."""
    return {
        _module_name(path): path
        for path in ROOT.rglob("*.py")
        if "data" not in path.relative_to(ROOT).parts
    }


def _load_time_imports(name: str, path: Path, modules) -> set[str]:
    """Package modules imported at module level by ``name``."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    targets = set()

    for node in tree.body:
        if isinstance(node, ast.Import):
            candidates = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            base = node.module
            if node.level:
                anchor = name.split(".")[: -node.level]
                base = ".".join(anchor + [node.module])
            candidates = [base] + [
                f"{base}.{alias.name}" for alias in node.names
            ]
        else:
            continue

        targets.update(
            candidate for candidate in candidates
            if candidate in modules and candidate != name
        )

    return targets


def _cycles(graph: dict[str, set[str]]) -> list[list[str]]:
    """Return the import cycles of a dependency graph."""
    cycles = []
    state: dict[str, int] = {}
    stack: list[str] = []

    def visit(node):
        state[node] = 1
        stack.append(node)
        for target in sorted(graph.get(node, ())):
            if state.get(target) == 1:
                cycles.append(stack[stack.index(target):] + [target])
            elif target not in state:
                visit(target)
        stack.pop()
        state[node] = 2

    for node in sorted(graph):
        if node not in state:
            visit(node)

    return cycles


# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    reason="import cycles between base, fourier, response and "
    "stationxml, being removed in the revision of signals",
)
def test_no_import_cycles():
    """Module-level imports within the package form no cycle."""
    modules = _modules()
    graph = {
        name: _load_time_imports(name, path, modules)
        for name, path in modules.items()
    }

    cycles = _cycles(graph)

    assert not cycles, "import cycles: " + "; ".join(
        " -> ".join(cycle) for cycle in cycles
    )


def test_cycle_detection():
    """The cycle finder reports cycles and accepts acyclic graphs."""
    assert _cycles({"a": {"b"}, "b": {"c"}, "c": set()}) == []
    assert _cycles({"a": {"b"}, "b": {"a"}}) == [["a", "b", "a"]]


def _import_params():
    """Module names, with the known broken modules marked."""
    params = []
    for name in sorted(_modules()):
        if name in BROKEN:
            # Skipped rather than expected to fail: importing them would
            # still load heavy dependencies (seedlink imports matplotlib).
            params.append(pytest.param(name, marks=pytest.mark.skip(
                reason="imports modules removed from the repository",
            )))
        else:
            params.append(name)
    return params


@pytest.mark.parametrize("name", _import_params())
def test_module_imports(name):
    """Every module of the package can be imported."""
    with warnings.catch_warnings():
        # Experimental modules announce themselves with a FutureWarning.
        warnings.simplefilter("ignore", FutureWarning)
        importlib.import_module(name)
