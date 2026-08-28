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
Bindi et al. (2011) ground-motion prediction equation.
"""

import numpy as _np
from scipy.constants import g

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class BindiEtAl2011(GMPE):
    """
    Ground-motion model by Bindi et al. (2011).

    Reference
    ---------
    Bindi, D., Pacor, F., Luzi, L., Puglia, R., Massa, M.,
    Ameri, G., and Paolucci, R. (2011).
    Ground motion prediction equations derived from the Italian
    strong motion database.
    Bulletin of Earthquake Engineering.
    DOI: 10.1007/s10518-011-9313-z.

    Notes
    -----
    The model uses moment magnitude (Mw), Joyner-Boore distance (Rjb),
    Vs30 and fault rake.

    PGA and SA are returned in g. PGV is returned in m/s. Mean ground
    motion and aleatory standard deviations are expressed in natural-log
    space.

    EC8 site class E follows the OpenQuake implementation convention and
    is selected by setting Vs30 to 0 m/s.
    """

    REFERENCE_VELOCITY = 800.0
    MAGNITUDE_TYPE = "Mw"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
        "rake",
    }

    REQUIRES_SITE_PARAMETERS = {
        "vs30",
    }

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_FILE = "bindi_et_al_2011.json"
    _COEFF_SET = "default"

    _MREF = 5.0
    _RREF = 1.0
    _MH = 6.75
    _B3 = 0.0

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the ground-motion model.

        Parameters
        ----------
        imt
            Intensity-measure type.
        rupture
            Rupture context containing ``mag`` and ``rake``.
        sites
            Sites context containing ``rjb`` in km and ``vs30`` in m/s.

        Returns
        -------
        GroundMotionResult
            Median ground motion and total/inter-/intra-event logarithmic
            standard deviations.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        C = self.get_coefficients(imt)

        mag = _np.asarray(
            rupture.mag,
            dtype=float,
        )
        rake = _np.asarray(
            rupture.rake,
            dtype=float,
        )
        rjb = _np.asarray(
            sites.rjb,
            dtype=float,
        )
        vs30 = _np.asarray(
            sites.vs30,
            dtype=float,
        )

        if _np.any(rjb < 0.0):
            raise ValueError(
                "Joyner-Boore distance must be non-negative."
            )

        if _np.any(
            (rake < -180.0)
            | (rake > 180.0)
        ):
            raise ValueError(
                "Rake must be within [-180, 180] degrees."
            )

        mean_log10 = (
            self._magnitude_term(
                mag,
                C,
            )
            + self._distance_term(
                mag,
                rjb,
                C,
            )
            + self._site_term(
                vs30,
                C,
            )
            + self._mechanism_term(
                rake,
                C,
            )
        )

        mean = self._convert_mean(
            imt,
            mean_log10,
        )

        ln10 = _np.log(10.0)

        sigma = C["SigmaTot"] * ln10
        tau = C["SigmaB"] * ln10
        phi = C["SigmaW"] * ln10

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    @classmethod
    def _magnitude_term(
        cls,
        mag,
        C,
    ):
        """
        Return the magnitude scaling term.
        """
        delta = mag - cls._MH

        return _np.where(
            mag <= cls._MH,
            (
                C["e1"]
                + C["b1"] * delta
                + C["b2"] * delta**2
            ),
            (
                C["e1"]
                + cls._B3 * delta
            ),
        )

    @classmethod
    def _distance_term(
        cls,
        mag,
        rjb,
        C,
    ):
        """
        Return the distance scaling term.
        """
        r = _np.sqrt(
            rjb**2 + C["h"]**2
        )

        return (
            (
                C["c1"]
                + C["c2"] * (mag - cls._MREF)
            )
            * _np.log10(r / cls._RREF)
            - C["c3"] * (r - cls._RREF)
        )

    @staticmethod
    def _site_term(
        vs30,
        C,
    ):
        """
        Return the EC8 site amplification term.
        """
        ssa = _np.zeros_like(
            vs30,
            dtype=float,
        )
        ssb = _np.zeros_like(
            vs30,
            dtype=float,
        )
        ssc = _np.zeros_like(
            vs30,
            dtype=float,
        )
        ssd = _np.zeros_like(
            vs30,
            dtype=float,
        )
        sse = _np.zeros_like(
            vs30,
            dtype=float,
        )

        # Class E is encoded as Vs30 = 0, matching OpenQuake.
        sse[_np.fabs(vs30) < 1.0e-10] = 1.0

        # Class D: Vs30 < 180 m/s.
        ssd[
            (vs30 >= 1.0e-10)
            & (vs30 < 180.0)
        ] = 1.0

        # Class C: 180 <= Vs30 < 360 m/s.
        ssc[
            (vs30 >= 180.0)
            & (vs30 < 360.0)
        ] = 1.0

        # Class B: 360 <= Vs30 < 800 m/s.
        ssb[
            (vs30 >= 360.0)
            & (vs30 < 800.0)
        ] = 1.0

        # Class A: Vs30 >= 800 m/s.
        ssa[
            vs30 >= 800.0
        ] = 1.0

        return (
            C["sA"] * ssa
            + C["sB"] * ssb
            + C["sC"] * ssc
            + C["sD"] * ssd
            + C["sE"] * sse
        )

    @staticmethod
    def _mechanism_term(
        rake,
        C,
    ):
        """
        Return the style-of-faulting term from the rake angle.
        """
        ss = _np.zeros_like(
            rake,
            dtype=float,
        )
        ns = _np.zeros_like(
            rake,
            dtype=float,
        )
        rs = _np.zeros_like(
            rake,
            dtype=float,
        )

        ss[
            (rake > -30.0)
            & (rake <= 30.0)
        ] = 1.0

        ss[
            (rake > 150.0)
            | (rake <= -150.0)
        ] = 1.0

        rs[
            (rake > 30.0)
            & (rake <= 150.0)
        ] = 1.0

        ns[
            (rake > -150.0)
            & (rake <= -30.0)
        ] = 1.0

        return (
            C["f1"] * ns
            + C["f2"] * rs
            + C["f3"] * ss
        )

    @staticmethod
    def _convert_mean(
        imt: str,
        mean_log10,
    ):
        """
        Convert model-native mean to ShakeLab natural-log units.
        """
        if imt == "PGA" or imt.startswith("SA-"):
            # Native acceleration is in cm/s^2.
            return _np.log(
                (10.0 ** (mean_log10 - 2.0))
                / g
            )

        if imt == "PGV":
            # Native PGV is in cm/s; ShakeLab uses m/s.
            return (
                mean_log10 * _np.log(10.0)
                - _np.log(100.0)
            )

        raise KeyError(
            f"Unsupported intensity measure type: {imt!r}"
        )


__all__ = [
    "BindiEtAl2011",
]
