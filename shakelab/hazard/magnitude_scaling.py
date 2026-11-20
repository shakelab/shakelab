# ****************************************************************************
#
# Copyright (C) 2019-2020, ShakeLab Developers.
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
#
# ****************************************************************************
"""Legacy access to magnitude scaling.

New code should import shakelab.groundmotion.scaling. This module keeps
the historical class name and scalar tuple-returning methods without
duplicating regression coefficients.
"""

from shakelab.groundmotion.scaling import WellsCoppersmith1994


_STYLE_ALIASES = {
    "g": "all", "generic": "all", "all": "all",
    "n": "normal", "normal": "normal",
    "r": "reverse", "reverse": "reverse",
    "s": "strike-slip", "strike-slip": "strike-slip",
}


def _model(style):
    try:
        canonical = _STYLE_ALIASES[style]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Unknown fault style: {style!r}.") from exc
    return WellsCoppersmith1994(canonical)


class wells_coppersmith_1994:
    """Compatibility interface; prefer WellsCoppersmith1994."""

    @staticmethod
    def area_to_magnitude(area, type="generic"):
        """Return predicted Mw and residual sigma in magnitude units."""
        estimate = _model(type).area_to_magnitude(area)
        return estimate.median, estimate.sigma

    @staticmethod
    def magnitude_to_area(mag, type="generic"):
        """Return median area (km**2) and residual log10-area sigma."""
        estimate = _model(type).magnitude_to_area(mag)
        return estimate.median, estimate.sigma
