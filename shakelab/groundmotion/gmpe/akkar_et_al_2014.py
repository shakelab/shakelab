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
Akkar, Sandikkaya and Bommer (2014) ground-motion models.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class _AkkarEtAl2014Base(GMPE):
    """
    Shared implementation of Akkar et al. (2014).

    Concrete subclasses select the distance metric and the corresponding
    regression coefficient set.
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

    _COEFF_FILE = "akkar_et_al_2014.json"

    _C1 = 6.75
    _DISTANCE_PARAMETER = None

    def __init__(
        self,
        adjustment_factor: float = 1.0,
    ):
        """
        Initialize the model.

        Parameters
        ----------
        adjustment_factor
            Multiplicative adjustment applied to median ground motion.
            The default value 1.0 leaves the model unchanged.
        """
        factor = float(adjustment_factor)

        if not _np.isfinite(factor):
            raise ValueError(
                "adjustment_factor must be finite."
            )

        if factor <= 0.0:
            raise ValueError(
                "adjustment_factor must be strictly positive."
            )

        if self._DISTANCE_PARAMETER is None:
            raise TypeError(
                "A concrete Akkar et al. (2014) subclass is required."
            )

        self.adjustment_factor = factor
        self._ln_adjustment = _np.log(factor)

        super().__init__()

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the selected Akkar et al. (2014) model.
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

        if rake < -180.0 or rake > 180.0:
            raise ValueError(
                "rake must be within [-180, 180] degrees."
            )

        distance = self._finite_site_array(
            getattr(
                sites,
                self._DISTANCE_PARAMETER,
            ),
            self._DISTANCE_PARAMETER,
        )
        vs30 = self._finite_site_array(
            sites.vs30,
            "vs30",
        )

        if _np.any(distance < 0.0):
            raise ValueError(
                f"{self._DISTANCE_PARAMETER} must be non-negative."
            )

        if _np.any(vs30 <= 0.0):
            raise ValueError(
                "vs30 must be strictly positive."
            )

        C = self.get_coefficients(imt)
        C_pga = self.get_coefficients("PGA")

        rock_pga = _np.exp(
            self._mean_without_site(
                C_pga,
                mag,
                rake,
                distance,
            )
        )

        mean = (
            self._mean_without_site(
                C,
                mag,
                rake,
                distance,
            )
            + self._nonlinear_site_term(
                C,
                rock_pga,
                vs30,
            )
            + self._ln_adjustment
        )

        if imt == "PGV":
            # OpenQuake retains PGV in cm/s; ShakeLab uses m/s.
            mean = mean - _np.log(100.0)

        phi = C["sigma"]
        tau = C["tau"]
        sigma = _np.sqrt(
            phi**2 + tau**2
        )

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    @classmethod
    def _mean_without_site(
        cls,
        C: dict,
        mag: float,
        rake: float,
        distance: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return median ground motion excluding nonlinear site effects.
        """
        return (
            C["a1"]
            + cls._linear_magnitude_term(
                C,
                mag,
            )
            + cls._quadratic_magnitude_term(
                C,
                mag,
            )
            + cls._distance_term(
                C,
                mag,
                distance,
            )
            + cls._faulting_style_term(
                C,
                rake,
            )
        )

    @classmethod
    def _linear_magnitude_term(
        cls,
        C: dict,
        mag: float,
    ) -> float:
        """
        Return piecewise linear magnitude scaling.
        """
        if mag <= cls._C1:
            return (
                C["a2"]
                * (
                    mag - cls._C1
                )
            )

        return (
            C["a7"]
            * (
                mag - cls._C1
            )
        )

    @staticmethod
    def _quadratic_magnitude_term(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return quadratic magnitude scaling.
        """
        return (
            C["a3"]
            * (
                8.5 - mag
            )**2
        )

    @classmethod
    def _distance_term(
        cls,
        C: dict,
        mag: float,
        distance: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return logarithmic distance attenuation.
        """
        return (
            (
                C["a4"]
                + C["a5"]
                * (
                    mag - cls._C1
                )
            )
            * _np.log(
                _np.sqrt(
                    distance**2
                    + C["a6"]**2
                )
            )
        )

    @staticmethod
    def _faulting_style_term(
        C: dict,
        rake: float,
    ) -> float:
        """
        Return normal/reverse style-of-faulting scaling.

        The classification follows OpenQuake exactly:
        normal for -135 < rake < -45 degrees and reverse for
        45 < rake < 135 degrees. All other rake values contribute zero.
        """
        normal = (
            rake > -135.0
            and rake < -45.0
        )
        reverse = (
            rake > 45.0
            and rake < 135.0
        )

        return (
            C["a8"] * float(normal)
            + C["a9"] * float(reverse)
        )

    @staticmethod
    def _nonlinear_site_term(
        C: dict,
        rock_pga: _np.ndarray,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return nonlinear site amplification.
        """
        vref = C["Vref"]
        vcon = C["Vcon"]

        site = _np.zeros_like(
            vs30,
            dtype=float,
        )

        idx = vs30 < vref
        ratio = (
            vs30[idx] / vref
        )

        site[idx] = (
            C["b1"]
            * _np.log(ratio)
            + C["b2"]
            * _np.log(
                (
                    rock_pga[idx]
                    + C["c"]
                    * ratio**C["n"]
                )
                / (
                    (
                        rock_pga[idx]
                        + C["c"]
                    )
                    * ratio**C["n"]
                )
            )
        )

        idx = (
            (vs30 >= vref)
            & (vs30 <= vcon)
        )

        site[idx] = (
            C["b1"]
            * _np.log(
                vs30[idx] / vref
            )
        )

        idx = vs30 > vcon

        site[idx] = (
            C["b1"]
            * _np.log(
                vcon / vref
            )
        )

        return site

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


class AkkarEtAlRjb2014(_AkkarEtAl2014Base):
    """
    Akkar, Sandikkaya and Bommer (2014), Joyner-Boore variant.
    """

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_SET = "rjb"
    _DISTANCE_PARAMETER = "rjb"


class AkkarEtAlRepi2014(_AkkarEtAl2014Base):
    """
    Akkar, Sandikkaya and Bommer (2014), epicentral-distance variant.
    """

    REQUIRES_DISTANCES = {
        "repi",
    }

    _COEFF_SET = "repi"
    _DISTANCE_PARAMETER = "repi"


class AkkarEtAlRhypo2014(_AkkarEtAl2014Base):
    """
    Akkar, Sandikkaya and Bommer (2014), hypocentral-distance variant.
    """

    REQUIRES_DISTANCES = {
        "rhypo",
    }

    _COEFF_SET = "rhypo"
    _DISTANCE_PARAMETER = "rhypo"


__all__ = [
    "AkkarEtAlRjb2014",
    "AkkarEtAlRepi2014",
    "AkkarEtAlRhypo2014",
]
