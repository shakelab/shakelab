# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
# ****************************************************************************
"""
Location and travel-time utilities for ShakeLab.

This subpackage groups together:

* Absolute Geiger hypocenter location (:mod:`locator`)
* Relative-time location from Δt between station pairs
  (:mod:`relative_locator`) – currently being refactored
* Travel-time engines (:mod:`ttengine`)
* 1D layered traveltime solver (:mod:`ttsolver`)
* Published 1D velocity profiles (:mod:`profiles`)

Typical high-level entry points re-exported here include:

* :class:`Station`, :class:`Hypocenter`
* :class:`CartStation`, :class:`CartHypocenter`
* :class:`PhasePick`
* :class:`TravelTimeEngine`
* :func:`locate_event`
* :class:`LayeredModel1D`, :func:`traveltime_solver`
"""

from .locator import (  # noqa: F401
    Station,
    Hypocenter,
    CartStation,
    CartHypocenter,
    LocationResult,
    build_design_matrices,
    build_normal_equations,
    compute_covariance,
    geiger_locate,
    locate_event,
)

from .phases import PhasePick  # noqa: F401

# ---------------------------------------------------------------------------
# Travel-time engines (new module)
# ---------------------------------------------------------------------------
from .ttengine import (  # noqa: F401
    TravelTimeEngine,
    HomogeneousEngine,
    LayeredTravelTimeEngine,
)

# ---------------------------------------------------------------------------
# 1D layered traveltime solver (new module)
# ---------------------------------------------------------------------------
try:
    from .ttsolver import (  # noqa: F401
        LayeredModel1D,
        solve_traveltime as traveltime_solver,
    )
except Exception:  # pragma: no cover
    LayeredModel1D = None  # type: ignore
    traveltime_solver = None  # type: ignore

__all__ = [
    # Geographic / Cartesian data structures
    "Station",
    "Hypocenter",
    "CartStation",
    "CartHypocenter",
    # Picks
    "PhasePick",
    # Core location machinery (absolute times)
    "LocationResult",
    "TravelTimeEngine",
    "HomogeneousEngine",
    "LayeredTravelTimeEngine",
    "build_design_matrices",
    "build_normal_equations",
    "compute_covariance",
    "geiger_locate",
    "locate_event",
    # Traveltime tools (if available)
    "LayeredModel1D",
    "traveltime_solver",
]
