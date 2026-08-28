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
Lanzano et al. (2019) ground-motion prediction equations for Italy.
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


class _LanzanoEtAl2019Base(GMPE):
    """
    Shared implementation of the nominal Lanzano et al. (2019) GMPE.

    Concrete subclasses select the distance metric and corresponding
    coefficient set.
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

    _COEFF_FILE = "lanzano_et_al_2019.json"
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
                "A concrete Lanzano et al. (2019) subclass is required."
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
        Evaluate the selected Lanzano et al. (2019) model.
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
                C,
                mag,
            )
            + self._distance_term(
                C,
                mag,
                distance,
            )
            + self._site_term(
                C,
                vs30,
            )
            + self._mechanism_term(
                C,
                rake,
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

        tau = C["tau"] * ln10
        phi = (
            _np.sqrt(
                C["phi_S2S"]**2
                + C["phi_0"]**2
            )
            * ln10
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

    @staticmethod
    def _magnitude_term(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return piecewise magnitude scaling.
        """
        delta = mag - C["Mh"]

        if mag <= C["Mh"]:
            return (
                C["a"]
                + C["b1"] * delta
            )

        return (
            C["a"]
            + C["b2"] * delta
        )

    @staticmethod
    def _distance_term(
        C: dict,
        mag: float,
        distance: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return distance attenuation.
        """
        effective = _np.sqrt(
            distance**2 + C["h"]**2
        )

        return (
            (
                C["c1"]
                * (
                    mag - C["Mref"]
                )
                + C["c2"]
            )
            * _np.log10(effective)
            + C["c3"] * effective
        )

    @staticmethod
    def _site_term(
        C: dict,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return continuous Vs30 site amplification.
        """
        velocity = _np.minimum(
            vs30,
            1500.0,
        )

        return (
            C["k"]
            * _np.log10(
                velocity / 800.0
            )
        )

    @staticmethod
    def _mechanism_term(
        C: dict,
        rake: float,
    ) -> float:
        """
        Return style-of-faulting scaling.

        OpenQuake uses strike-slip and reverse dummy variables; normal
        faulting is the reference category.
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

        return (
            C["f1"] * float(strike_slip)
            + C["f2"] * float(reverse)
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


class LanzanoEtAl2019Rjb(_LanzanoEtAl2019Base):
    """
    Lanzano et al. (2019), Joyner-Boore distance variant.

    The model predicts RotD50 ground motion for shallow crustal
    earthquakes in Italy.
    """

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_SET = "rjb"
    _DISTANCE_PARAMETER = "rjb"


class LanzanoEtAl2019Rrup(_LanzanoEtAl2019Base):
    """
    Lanzano et al. (2019), rupture-distance variant.

    The model predicts RotD50 ground motion for shallow crustal
    earthquakes in Italy.
    """

    REQUIRES_DISTANCES = {
        "rrup",
    }

    _COEFF_SET = "rrup"
    _DISTANCE_PARAMETER = "rrup"


__all__ = [
    "LanzanoEtAl2019Rjb",
    "LanzanoEtAl2019Rrup",
]
