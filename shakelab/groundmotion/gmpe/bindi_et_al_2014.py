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
Bindi et al. (2014) Pan-European ground-motion prediction equations.

The module provides the continuous-Vs30 variants based on Joyner-Boore
and hypocentral distance.
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


class _BindiEtAl2014Base(GMPE):
    """
    Shared implementation of the continuous-Vs30 Bindi et al. (2014) GMPE.

    Subclasses select the distance metric and the corresponding regression
    coefficients.
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

    _COEFF_FILE = "bindi_et_al_2014.json"

    _MREF = 5.5
    _MH = 6.75
    _RREF = 1.0
    _VREF = 800.0

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
                "A concrete Bindi et al. (2014) subclass is required."
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
        Evaluate the selected Bindi et al. (2014) model.
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

        mean_log10 = (
            self._magnitude_term(
                mag,
                C,
            )
            + self._distance_term(
                mag,
                distance,
                C,
            )
            + self._site_term(
                vs30,
                C,
            )
            + self._style_of_faulting_term(
                rake,
                C,
            )
        )

        mean = (
            self._convert_mean(
                imt,
                mean_log10,
            )
            + self._ln_adjustment
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

    @classmethod
    def _magnitude_term(
        cls,
        mag: float,
        C: dict,
    ) -> float:
        """
        Return magnitude scaling.
        """
        delta = mag - cls._MH

        if mag < cls._MH:
            return (
                C["e1"]
                + C["b1"] * delta
                + C["b2"] * delta**2
            )

        return (
            C["e1"]
            + C["b3"] * delta
        )

    @classmethod
    def _distance_term(
        cls,
        mag: float,
        distance: _np.ndarray,
        C: dict,
    ) -> _np.ndarray:
        """
        Return distance scaling.
        """
        adjusted = _np.sqrt(
            distance**2 + C["h"]**2
        )

        return (
            (
                C["c1"]
                + C["c2"]
                * (
                    mag - cls._MREF
                )
            )
            * _np.log10(
                adjusted / cls._RREF
            )
            - C["c3"]
            * (
                adjusted - cls._RREF
            )
        )

    @classmethod
    def _site_term(
        cls,
        vs30: _np.ndarray,
        C: dict,
    ) -> _np.ndarray:
        """
        Return continuous-Vs30 site amplification.
        """
        return (
            C["gamma"]
            * _np.log10(
                vs30 / cls._VREF
            )
        )

    @staticmethod
    def _style_of_faulting_term(
        rake: float,
        C: dict,
    ) -> float:
        """
        Return the style-of-faulting term.

        The rake classification follows OpenQuake exactly:
        strike-slip for |rake| <= 30 degrees or within 30 degrees
        of +/-180 degrees, reverse for 30 < rake < 150 degrees,
        and normal otherwise.
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

        if reverse:
            return C["sofR"]

        if strike_slip:
            return C["sofS"]

        return C["sofN"]

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


class BindiEtAl2014Rjb(_BindiEtAl2014Base):
    """
    Bindi et al. (2014) using Joyner-Boore distance.

    Reference
    ---------
    Bindi, D., Massa, M., Luzi, L., Ameri, G., Pacor, F.,
    Puglia, R., and Augliera, P. (2014).
    Pan-European ground motion prediction equations for the average
    horizontal component of PGA, PGV and 5%-damped PSA at spectral
    periods of up to 3.0 s using the RESORCE dataset.
    Bulletin of Earthquake Engineering, 12, 391-430.
    """

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_SET = "rjb"
    _DISTANCE_PARAMETER = "rjb"


class BindiEtAl2014Rhyp(_BindiEtAl2014Base):
    """
    Bindi et al. (2014) using hypocentral distance.

    This variant uses the dedicated hypocentral-distance regression
    coefficients from the electronic supplementary material.
    """

    REQUIRES_DISTANCES = {
        "rhypo",
    }

    _COEFF_SET = "rhypo"
    _DISTANCE_PARAMETER = "rhypo"


__all__ = [
    "BindiEtAl2014Rjb",
    "BindiEtAl2014Rhyp",
]
