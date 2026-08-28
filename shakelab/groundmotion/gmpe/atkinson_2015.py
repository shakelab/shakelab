# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Atkinson (2015) induced-seismicity ground-motion prediction equations.
"""

from __future__ import annotations

import numpy as _np
from scipy.constants import g as _g

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class Atkinson2015(GMPE):
    """
    Atkinson (2015) GMPE for induced seismicity.

    Reference
    ---------
    Atkinson, G. M. (2015).
    Ground-Motion Prediction Equation for Small-to-Moderate Events at
    Short Hypocentral Distances, with Application to Induced-Seismicity
    Hazards. Bulletin of the Seismological Society of America, 105(2).

    Notes
    -----
    The model requires moment magnitude and hypocentral distance only.
    No explicit site parameter is required; the model was derived using
    B/C site amplification factors.

    PGA and SA are returned in g. PGV is returned in m/s. Mean ground
    motion and aleatory standard deviations are expressed in natural-log
    space.
    """

    REFERENCE_VELOCITY = None
    MAGNITUDE_TYPE = "Mw"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
    }

    REQUIRES_SITE_PARAMETERS = set()

    REQUIRES_DISTANCES = {
        "rhypo",
    }

    _COEFF_FILE = "atkinson_2015.json"
    _COEFF_SET = "default"

    _DISTANCE_SATURATION = "default"

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the Atkinson (2015) model.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        C = self.get_coefficients(imt)

        mag = self._finite_scalar(
            rupture.mag,
            "mag",
        )
        rhypo = self._finite_site_array(
            sites.rhypo,
            "rhypo",
        )

        if _np.any(rhypo < 0.0):
            raise ValueError(
                "rhypo must be non-negative."
            )

        effective_depth = self._effective_depth(
            mag,
        )
        effective_distance = _np.sqrt(
            rhypo**2 + effective_depth**2
        )

        mean_log10 = (
            C["c0"]
            + C["c1"] * mag
            + C["c2"] * mag**2
            + C["c3"]
            * _np.log10(
                effective_distance
            )
            + C["c4"]
            * effective_distance
        )

        mean = self._convert_mean(
            imt,
            mean_log10,
        )

        ln10 = _np.log(10.0)

        sigma = C["sigma"] * ln10
        tau = C["tau"] * ln10
        phi = C["phi"] * ln10

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    def _effective_depth(
        self,
        mag: float,
    ) -> float:
        """
        Return the magnitude-dependent effective depth in km.
        """
        if self._DISTANCE_SATURATION == "default":
            value = 10.0 ** (
                -1.72 + 0.43 * mag
            )
        elif self._DISTANCE_SATURATION == "alternative":
            value = 10.0 ** (
                -0.28 + 0.19 * mag
            )
        else:
            raise ValueError(
                "Unknown distance-saturation model."
            )

        return max(
            1.0,
            value,
        )

    @staticmethod
    def _convert_mean(
        imt: str,
        mean_log10,
    ):
        """
        Convert model-native median to ShakeLab natural-log units.
        """
        if (
            imt == "PGA"
            or imt.startswith("SA-")
        ):
            return _np.log(
                (
                    10.0 ** (
                        mean_log10 - 2.0
                    )
                )
                / _g
            )

        if imt == "PGV":
            # OpenQuake retains PGV in cm/s; ShakeLab uses m/s.
            return (
                mean_log10 * _np.log(10.0)
                - _np.log(100.0)
            )

        raise KeyError(
            f"Unsupported intensity measure type: {imt!r}."
        )

    @staticmethod
    def _finite_scalar(
        value,
        name: str,
    ) -> float:
        """
        Return a finite scalar float.
        """
        array = _np.asarray(
            value,
            dtype=float,
        )

        if array.ndim != 0:
            raise ValueError(
                f"{name} must be scalar."
            )

        number = float(array)

        if not _np.isfinite(number):
            raise ValueError(
                f"{name} must be finite."
            )

        return number

    @staticmethod
    def _finite_site_array(
        value,
        name: str,
    ) -> _np.ndarray:
        """
        Return a finite one-dimensional site array.
        """
        array = _np.asarray(
            value,
            dtype=float,
        )

        if array.ndim != 1:
            raise ValueError(
                f"{name} must be one-dimensional."
            )

        if not _np.all(
            _np.isfinite(array)
        ):
            raise ValueError(
                f"{name} must contain finite values."
            )

        return array


class Atkinson2015AltDistSat(Atkinson2015):
    """
    Atkinson (2015) with alternative stronger distance saturation.

    The subclass uses the alternative effective-depth relation described
    in Atkinson (2015), while retaining the same regression coefficients.
    """

    _DISTANCE_SATURATION = "alternative"


__all__ = [
    "Atkinson2015",
    "Atkinson2015AltDistSat",
]
