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
Boore et al. (2014) BSSA14 ground-motion prediction equation.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class BooreEtAl2014(GMPE):
    """
    Boore et al. (2014) NGA-West2 ground-motion model (BSSA14).

    Reference
    ---------
    Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M.
    (2014). NGA-West2 equations for predicting PGA, PGV, and 5%-damped
    PSA for shallow crustal earthquakes.
    Earthquake Spectra, 30(3), 1057-1085.

    Notes
    -----
    This class implements the nominal BSSA14 model without a basin-depth
    term. It corresponds to the OpenQuake ``BooreEtAl2014`` model with
    ``region='nobasin'``.

    The model requires moment magnitude, rake, Vs30 and Joyner-Boore
    distance. The style-of-faulting term can be disabled with ``sof=False``.

    PGA and SA are returned in g. PGV is returned in m/s. Mean ground
    motion and aleatory standard deviations are expressed in natural-log
    space.
    """

    REFERENCE_VELOCITY = 760.0
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

    _COEFF_FILE = "boore_et_al_2014.json"
    _COEFF_SET = "default"

    _MREF = 4.5
    _RREF = 1.0
    _VREF = 760.0
    _F1 = 0.0
    _F3 = 0.1
    _V1 = 225.0
    _V2 = 300.0

    def __init__(
        self,
        sof: bool = True,
    ):
        """
        Initialize the nominal BSSA14 model.

        Parameters
        ----------
        sof
            If True, include the style-of-faulting term. If False, use
            the unspecified-style coefficient ``e0``.
        """
        if not isinstance(sof, bool):
            raise TypeError("sof must be a boolean.")

        self.sof = sof

        super().__init__()

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate BSSA14 for one rupture and many sites.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        mag = self._finite_scalar(
            rupture.mag,
            "mag",
        )
        rake = self._finite_scalar(
            rupture.rake,
            "rake",
        )
        rjb = self._finite_site_array(
            sites.rjb,
            "rjb",
        )
        vs30 = self._finite_site_array(
            sites.vs30,
            "vs30",
        )

        if rake < -180.0 or rake > 180.0:
            raise ValueError(
                "rake must be within [-180, 180] degrees."
            )

        if _np.any(rjb < 0.0):
            raise ValueError(
                "rjb must be non-negative."
            )

        if _np.any(vs30 <= 0.0):
            raise ValueError(
                "vs30 must be strictly positive."
            )

        C = self.get_coefficients(imt)
        C_pga = self.get_coefficients("PGA")

        pga_rock = _np.exp(
            self._magnitude_term(
                C_pga,
                mag,
                rake,
            )
            + self._path_term(
                C_pga,
                mag,
                rjb,
            )
        )

        mean = (
            self._magnitude_term(
                C,
                mag,
                rake,
            )
            + self._path_term(
                C,
                mag,
                rjb,
            )
            + self._linear_site_term(
                C,
                vs30,
            )
            + self._nonlinear_site_term(
                C,
                vs30,
                pga_rock,
            )
        )

        if imt == "PGV":
            # OpenQuake returns PGV in cm/s; ShakeLab uses m/s.
            mean = mean - _np.log(100.0)

        tau = self._inter_event_tau(
            C,
            mag,
        )
        phi = self._intra_event_phi(
            C,
            mag,
            rjb,
            vs30,
        )
        sigma = _np.sqrt(
            tau**2 + phi**2
        )

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    def _magnitude_term(
        self,
        C: dict,
        mag: float,
        rake: float,
    ) -> float:
        """
        Return magnitude and style-of-faulting scaling.
        """
        delta = mag - C["Mh"]

        if mag <= C["Mh"]:
            magnitude = (
                C["e4"] * delta
                + C["e5"] * delta**2
            )
        else:
            magnitude = (
                C["e6"] * delta
            )

        if not self.sof:
            return (
                C["e0"] + magnitude
            )

        return (
            self._style_of_faulting_term(
                C,
                rake,
            )
            + magnitude
        )

    @staticmethod
    def _style_of_faulting_term(
        C: dict,
        rake: float,
    ) -> float:
        """
        Return BSSA14 style-of-faulting scaling.

        The classification follows OpenQuake: strike-slip for rake
        within 30 degrees of horizontal (including +/-180), reverse
        for 30 < rake < 150 degrees, and normal otherwise.
        """
        strike_slip = (
            abs(rake) <= 30.0
            or (
                180.0 - abs(rake)
            ) <= 30.0
        )

        reverse = (
            rake > 30.0
            and rake < 150.0
        )

        if strike_slip:
            return C["e1"]

        if reverse:
            return C["e3"]

        return C["e2"]

    @classmethod
    def _path_term(
        cls,
        C: dict,
        mag: float,
        rjb: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return BSSA14 path scaling.
        """
        rval = _np.sqrt(
            rjb**2 + C["h"]**2
        )

        geometric = (
            C["c1"]
            + C["c2"]
            * (
                mag - cls._MREF
            )
        ) * _np.log(
            rval / cls._RREF
        )

        anelastic = (
            C["c3"] + C["Dc3"]
        ) * (
            rval - cls._RREF
        )

        return geometric + anelastic

    @classmethod
    def _linear_site_term(
        cls,
        C: dict,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return linear site scaling.
        """
        velocity_ratio = (
            vs30 / cls._VREF
        )

        velocity_ratio = velocity_ratio.copy()

        idx = vs30 > C["Vc"]

        velocity_ratio[idx] = (
            C["Vc"] / cls._VREF
        )

        return (
            C["c"]
            * _np.log(
                velocity_ratio
            )
        )

    @classmethod
    def _nonlinear_site_term(
        cls,
        C: dict,
        vs30: _np.ndarray,
        pga_rock: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return nonlinear site scaling.
        """
        velocity = _np.minimum(
            vs30,
            cls._VREF,
        )

        f2 = (
            C["f4"]
            * (
                _np.exp(
                    C["f5"]
                    * (
                        velocity - 360.0
                    )
                )
                - _np.exp(
                    C["f5"] * 400.0
                )
            )
        )

        return (
            cls._F1
            + f2
            * _np.log(
                (
                    pga_rock
                    + cls._F3
                )
                / cls._F3
            )
        )

    @staticmethod
    def _inter_event_tau(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return magnitude-dependent inter-event tau.
        """
        if mag <= 4.5:
            return C["tau1"]

        if mag >= 5.5:
            return C["tau2"]

        return (
            C["tau1"]
            + (
                C["tau2"]
                - C["tau1"]
            )
            * (
                mag - 4.5
            )
        )

    @classmethod
    def _intra_event_phi(
        cls,
        C: dict,
        mag: float,
        rjb: _np.ndarray,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return magnitude-, distance- and site-dependent intra-event phi.
        """
        if mag <= 4.5:
            base = C["f1"]

        elif mag >= 5.5:
            base = C["f2"]

        else:
            base = (
                C["f1"]
                + (
                    C["f2"]
                    - C["f1"]
                )
                * (
                    mag - 4.5
                )
            )

        phi = _np.full_like(
            rjb,
            base,
            dtype=float,
        )

        idx = rjb > C["R2"]

        phi[idx] += C["DfR"]

        idx = (
            (rjb > C["R1"])
            & (rjb <= C["R2"])
        )

        phi[idx] += (
            C["DfR"]
            * _np.log(
                rjb[idx] / C["R1"]
            )
            / _np.log(
                C["R2"] / C["R1"]
            )
        )

        idx = vs30 <= cls._V1

        phi[idx] -= C["DfV"]

        idx = (
            (vs30 >= cls._V1)
            & (vs30 <= cls._V2)
        )

        phi[idx] -= (
            C["DfV"]
            * _np.log(
                cls._V2 / vs30[idx]
            )
            / _np.log(
                cls._V2 / cls._V1
            )
        )

        return phi

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


__all__ = [
    "BooreEtAl2014",
]
