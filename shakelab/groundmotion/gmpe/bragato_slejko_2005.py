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
Bragato and Slejko (2005) ground-motion prediction equations.
"""

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class BragatoSlejko2005(GMPE):
    """
    Ground-motion model by Bragato and Slejko (2005).

    Reference
    ---------
    Bragato, P. L. and Slejko, D. (2005).
    Empirical Ground-Motion Attenuation Relations for the Eastern Alps
    in the Magnitude Range 2.5-6.3.
    Bulletin of the Seismological Society of America,
    95(1), 252-276.

    Notes
    -----
    The reference model uses local magnitude (ML) and epicentral
    distance. PGA is returned in g and PGV in m/s.

    The median and standard deviation are returned in natural-log space.
    """

    REFERENCE_VELOCITY = 800.0
    MAGNITUDE_TYPE = "ML"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
    }

    REQUIRES_SITE_PARAMETERS = set()

    REQUIRES_DISTANCES = {
        "repi",
    }

    _COEFF_FILE = "bragato_slejko_2005.json"
    _COEFF_SET = "epicentral"

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the epicentral-distance model.

        Parameters
        ----------
        imt
            Intensity-measure type.
        rupture
            Rupture context containing ``mag``.
        sites
            Sites context containing ``repi`` in km.

        Returns
        -------
        GroundMotionResult
            Median ground motion and total logarithmic standard deviation.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        return self._evaluate(
            imt,
            rupture.mag,
            sites.repi,
        )

    def _evaluate(
        self,
        imt,
        mag,
        dist,
    ) -> GroundMotionResult:
        """
        Evaluate the Bragato-Slejko functional form.
        """
        C = self.get_coefficients(imt)

        mag = _np.asarray(
            mag,
            dtype=float,
        )
        dist = _np.asarray(
            dist,
            dtype=float,
        )

        if _np.any(dist < 0.0):
            raise ValueError(
                "Distance must be non-negative."
            )

        # Effective distance.
        r = _np.sqrt(
            dist**2 + C["h"]**2
        )

        # Mean in base-10 logarithmic space.
        mean = (
            C["a"]
            + (C["b"] + C["c"] * mag) * mag
        )
        mean += (
            C["d"] + C["e"] * mag**3
        ) * _np.log10(r)

        # Standard deviation in base-10 logarithmic space.
        sigma = C["s"]

        # Convert mean and standard deviation to natural-log space.
        mean = mean * _np.log(10.0)
        sigma = sigma * _np.log(10.0)

        if imt == "PGV":
            # Convert from cm/s to m/s.
            mean = mean - _np.log(100.0)

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
        )


class BragatoSlejko2005JB(BragatoSlejko2005):
    """
    Joyner-Boore distance variant of Bragato and Slejko (2005).
    """

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_SET = "joyner-boore"

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the Joyner-Boore distance model.

        Parameters
        ----------
        imt
            Intensity-measure type.
        rupture
            Rupture context containing ``mag``.
        sites
            Sites context containing ``rjb`` in km.

        Returns
        -------
        GroundMotionResult
            Median ground motion and total logarithmic standard deviation.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        return self._evaluate(
            imt,
            rupture.mag,
            sites.rjb,
        )


__all__ = [
    "BragatoSlejko2005",
    "BragatoSlejko2005JB",
]
