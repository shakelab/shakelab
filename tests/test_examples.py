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
Smoke test of the examples.

Every script in ``examples/`` is run with a non-interactive matplotlib
backend, so that the examples stay in line with the library. The test
is skipped if matplotlib is not installed.
"""

from __future__ import annotations

import runpy
import warnings
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
SCRIPTS = sorted(EXAMPLES.rglob("*.py"))


@pytest.mark.parametrize(
    "script", SCRIPTS, ids=[str(p.relative_to(EXAMPLES)) for p in SCRIPTS]
)
def test_example_runs(script, monkeypatch):
    """The example runs to the end without errors."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    monkeypatch.chdir(EXAMPLES.parent)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        runpy.run_path(str(script), run_name="__main__")
    plt.close("all")
