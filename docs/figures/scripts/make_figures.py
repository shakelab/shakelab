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
Regenerate the figures of the documentation.

Usage, from any folder, with ShakeLab and matplotlib installed:

    python docs/figures/scripts/make_figures.py              # all
    python docs/figures/scripts/make_figures.py signals_waveforms
    python docs/figures/scripts/make_figures.py signals_waveforms filters

The images are written in docs/figures, organized as the guides.
"""

from __future__ import annotations

import importlib
import sys

__all__ = ["GUIDES", "main"]

# One module per guide, each with a dictionary FIGURES of functions
GUIDES = ("signals_waveforms", "signals_response", "signals_spectral")


def main(argv: list[str]) -> int:
    """Make the figures named on the command line (all by default)."""
    guides = argv[:1] or list(GUIDES)
    names = argv[1:]
    for guide in guides:
        if guide not in GUIDES:
            print(f"unknown guide {guide!r}; known: {', '.join(GUIDES)}")
            return 1
        figures = importlib.import_module(guide).FIGURES
        for name in names or figures:
            if name not in figures:
                print(f"unknown figure {name!r} of {guide}; "
                      f"known: {', '.join(figures)}")
                return 1
            print(figures[name]())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
