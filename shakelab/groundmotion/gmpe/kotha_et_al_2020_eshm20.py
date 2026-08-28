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
Kotha et al. (2020) ESHM20 ground-motion prediction equation.

This implementation follows the ESHM20 adaptation implemented in
OpenQuake hazardlib while using the ShakeLab GMPE interface and units.
"""

from __future__ import annotations

import json as _json
import os as _os
from collections.abc import Mapping as _Mapping

import numpy as _np
from scipy.constants import g as _g

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class KothaEtAl2020ESHM20(GMPE):
    """
    Kotha et al. (2020), ESHM20 adaptation.

    References
    ----------
    Kotha, S. R., Weatherill, G., Bindi, D., and Cotton, F. (2020).
    A regionally-adaptable ground-motion model for shallow crustal
    earthquakes in Europe. Bulletin of Earthquake Engineering,
    18, 4091-4125.

    Weatherill, G., Kotha, S. R., and Cotton, F. (2020).
    A regionally-adaptable scaled-backbone ground motion logic tree
    for shallow seismicity in Europe: application to the 2020
    European Seismic Hazard Model. Bulletin of Earthquake Engineering,
    18, 5087-5117.

    Notes
    -----
    Required rupture parameters are moment magnitude and hypocentral
    depth. Required site parameters are Vs30, a flag specifying whether
    Vs30 is measured, and the ESHM20 residual attenuation region.
    Joyner-Boore distance is required.

    ``eshm20_region`` must be an integer from 0 to 5. Region 0 applies
    the generic Kotha et al. (2020) attenuation coefficient; regions
    1-5 apply the ESHM20 regional attenuation coefficients.

    The model returns PGA and SA in g, and PGV in m/s. Mean ground
    motion and aleatory standard deviations are expressed in natural-log
    space.
    """

    MAGNITUDE_TYPE = "Mw"
    REFERENCE_VELOCITY = None

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
        "hypo_depth",
    }

    REQUIRES_SITE_PARAMETERS = {
        "vs30",
        "vs30measured",
        "eshm20_region",
    }

    REQUIRES_DISTANCES = {
        "rjb",
    }

    _COEFF_FILE = "kotha_et_al_2020_eshm20.json"
    _COEFF_SET = "default"

    _MREF = 4.5
    _RREF = 30.0
    _MH = 5.7

    _H_SHALLOW = 4.0
    _H_INTERMEDIATE = 8.0
    _H_DEEP = 12.0

    _GLOBAL_TAU = {
        "PGV": {
            "tau1": 0.3733,
            "tau2": 0.3639,
            "tau3": 0.3434,
            "tau4": 0.3236,
        },
        "SA": {
            "tau1": 0.4518,
            "tau2": 0.4270,
            "tau3": 0.3863,
            "tau4": 0.3508,
        },
    }

    def __init__(
        self,
        sigma_mu_epsilon: float = 0.0,
        c3_epsilon: float = 0.0,
        ergodic: bool = True,
        dl2l: _Mapping | None = None,
        c3: _Mapping | None = None,
    ):
        """
        Initialize the ESHM20 model.

        Parameters
        ----------
        sigma_mu_epsilon
            Number of epistemic standard deviations applied to the
            source-region scaling term.
        c3_epsilon
            Number of standard deviations applied to the residual
            attenuation coefficient.
        ergodic
            If True, include the site-to-site component in phi.
        dl2l
            Optional mapping of IMT to direct source-region adjustment
            in natural-log units. If supplied, it takes precedence over
            ``sigma_mu_epsilon``.
        c3
            Optional mapping of IMT to direct c3 values. If supplied,
            it overrides ESHM20 attenuation regionalization.
        """
        self.sigma_mu_epsilon = self._finite_scalar(
            sigma_mu_epsilon,
            "sigma_mu_epsilon",
        )
        self.c3_epsilon = self._finite_scalar(
            c3_epsilon,
            "c3_epsilon",
        )

        if not isinstance(ergodic, bool):
            raise TypeError("ergodic must be a boolean.")

        self.ergodic = ergodic
        self.dl2l = self._normalize_adjustment_mapping(
            dl2l,
            "dl2l",
        )
        self.c3 = self._normalize_adjustment_mapping(
            c3,
            "c3",
        )

        super().__init__()

        self._c3_regions = self._load_coefficient_set(
            "c3_regions"
        )
        self._hetero_phi0 = self._load_coefficient_set(
            "hetero_phi0"
        )
        self._sigma_mu = self._load_coefficient_set(
            "sigma_mu"
        )

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate the ESHM20 model.
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
        hypo_depth = self._finite_scalar(
            rupture.hypo_depth,
            "hypo_depth",
        )

        if hypo_depth < 0.0:
            raise ValueError(
                "hypo_depth must be non-negative."
            )

        rjb = self._finite_site_array(
            sites.rjb,
            "rjb",
        )
        vs30 = self._finite_site_array(
            sites.vs30,
            "vs30",
        )
        measured = self._boolean_site_array(
            sites.vs30measured,
            "vs30measured",
        )
        region = self._region_array(
            sites.eshm20_region,
        )

        if _np.any(rjb < 0.0):
            raise ValueError(
                "rjb must be non-negative."
            )

        if _np.any(vs30 <= 0.0):
            raise ValueError(
                "vs30 must be strictly positive."
            )

        mean = (
            self._magnitude_term(
                mag,
                C,
            )
            + self._distance_term(
                imt,
                mag,
                hypo_depth,
                rjb,
                region,
                C,
            )
            + self._site_term(
                vs30,
                measured,
                C,
            )
        )

        mean = self._convert_mean(
            imt,
            mean,
        )

        if self.dl2l is not None:
            mean = mean + self._mapping_value(
                self.dl2l,
                imt,
                "dl2l",
            )

        elif self.sigma_mu_epsilon != 0.0:
            mean = mean + (
                self.sigma_mu_epsilon
                * self._sigma_mu_adjustment(
                    imt,
                    mag,
                    hypo_depth,
                    C,
                )
            )

        tau = self._tau(
            imt,
            mag,
        )
        phi0 = self._phi0(
            imt,
            mag,
        )

        if self.ergodic:
            phi_s2s = _np.where(
                measured,
                C["phi_s2s_obs"],
                C["phi_s2s_inf"],
            )
            phi = _np.sqrt(
                phi0**2 + phi_s2s**2
            )
        else:
            phi = phi0

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
        mag: float,
        C: dict,
    ) -> float:
        """
        Return magnitude scaling.
        """
        delta = mag - self._MH

        if mag <= self._MH:
            return (
                C["e1"]
                + C["b1"] * delta
                + C["b2"] * delta**2
            )

        return (
            C["e1"]
            + C["b3"] * delta
        )

    def _distance_term(
        self,
        imt: str,
        mag: float,
        hypo_depth: float,
        rjb: _np.ndarray,
        region: _np.ndarray,
        C: dict,
    ) -> _np.ndarray:
        """
        Return distance attenuation term.
        """
        h = self._pseudo_depth(
            hypo_depth
        )

        rval = _np.sqrt(
            rjb**2 + h**2
        )
        rref = _np.sqrt(
            self._RREF**2 + h**2
        )

        c3 = self._c3_values(
            imt,
            region,
            C,
        )

        return (
            (
                C["c1"]
                + C["c2"] * (
                    mag - self._MREF
                )
            )
            * _np.log(
                rval / rref
            )
            + c3
            * (
                rval - rref
            )
            / 100.0
        )

    def _c3_values(
        self,
        imt: str,
        region: _np.ndarray,
        C: dict,
    ) -> _np.ndarray:
        """
        Return site-specific residual attenuation coefficients.
        """
        if self.c3 is not None:
            value = self._mapping_value(
                self.c3,
                imt,
                "c3",
            )
            return _np.full(
                region.shape,
                value,
                dtype=float,
            )

        c3 = _np.full(
            region.shape,
            C["c3"],
            dtype=float,
        )
        tau_c3 = _np.full(
            region.shape,
            C["tau_c3"],
            dtype=float,
        )

        if imt != "PGV":
            try:
                CR = self._c3_regions[imt]
            except KeyError as exc:
                raise KeyError(
                    f"No ESHM20 c3 regional coefficients "
                    f"available for {imt!r}."
                ) from exc

            for region_id in range(1, 6):
                idx = region == region_id

                if not _np.any(idx):
                    continue

                c3[idx] = CR[
                    f"region_{region_id}"
                ]
                tau_c3[idx] = CR[
                    f"tau_region_{region_id}"
                ]

        return (
            c3
            + self.c3_epsilon * tau_c3
        )

    @staticmethod
    def _site_term(
        vs30: _np.ndarray,
        measured: _np.ndarray,
        C: dict,
    ) -> _np.ndarray:
        """
        Return ESHM20 Vs30 site amplification.
        """
        velocity = _np.minimum(
            vs30,
            1100.0,
        )

        amplification = _np.empty_like(
            velocity,
            dtype=float,
        )

        amplification[measured] = (
            C["d0_obs"]
            + C["d1_obs"]
            * _np.log(
                velocity[measured]
            )
        )

        inferred = ~measured

        amplification[inferred] = (
            C["d0_inf"]
            + C["d1_inf"]
            * _np.log(
                velocity[inferred]
            )
        )

        return amplification

    @staticmethod
    def _convert_mean(
        imt: str,
        mean: _np.ndarray,
    ) -> _np.ndarray:
        """
        Convert model-native median to ShakeLab units.
        """
        if (
            imt == "PGA"
            or imt.startswith("SA-")
        ):
            return (
                mean
                - _np.log(
                    100.0 * _g
                )
            )

        if imt == "PGV":
            # OpenQuake retains cm/s for PGV; ShakeLab uses m/s.
            return (
                mean
                - _np.log(100.0)
            )

        raise KeyError(
            f"Unsupported intensity measure type: {imt!r}."
        )

    def _tau(
        self,
        imt: str,
        mag: float,
    ) -> float:
        """
        Return Al Atik (2015) global heteroskedastic tau.
        """
        if imt == "PGV":
            C = self._GLOBAL_TAU["PGV"]
        else:
            C = self._GLOBAL_TAU["SA"]

        if mag <= 4.5:
            return C["tau1"]

        if mag <= 5.0:
            return self._linear_interpolation(
                mag,
                4.5,
                5.0,
                C["tau1"],
                C["tau2"],
            )

        if mag <= 5.5:
            return self._linear_interpolation(
                mag,
                5.0,
                5.5,
                C["tau2"],
                C["tau3"],
            )

        if mag <= 6.5:
            return self._linear_interpolation(
                mag,
                5.5,
                6.5,
                C["tau3"],
                C["tau4"],
            )

        return C["tau4"]

    def _phi0(
        self,
        imt: str,
        mag: float,
    ) -> float:
        """
        Return ESHM20 heteroskedastic single-station phi.
        """
        try:
            C = self._hetero_phi0[imt]
        except KeyError as exc:
            raise KeyError(
                f"No heteroskedastic phi0 coefficients "
                f"available for {imt!r}."
            ) from exc

        if mag <= 5.0:
            return C["a"]

        if mag > 6.5:
            return C["b"]

        return (
            C["a"]
            + (
                mag - 5.0
            )
            * (
                C["b"] - C["a"]
            )
            / 1.5
        )

    def _sigma_mu_adjustment(
        self,
        imt: str,
        mag: float,
        hypo_depth: float,
        C: dict,
    ) -> float:
        """
        Return epistemic source-region adjustment sigma_mu.
        """
        try:
            CS = self._sigma_mu[imt]
        except KeyError as exc:
            raise KeyError(
                f"No sigma_mu coefficients available "
                f"for {imt!r}."
            ) from exc

        if hypo_depth < 10.0:
            upper = CS[
                "sigma_mu_m8_shallow"
            ]
            lower = CS[
                "sigma_mu_m7p4_shallow"
            ]

        elif hypo_depth >= 20.0:
            upper = CS[
                "sigma_mu_m8_deep"
            ]
            lower = CS[
                "sigma_mu_m7p4_deep"
            ]

        else:
            upper = CS[
                "sigma_mu_m8_intermediate"
            ]
            lower = CS[
                "sigma_mu_m7p4_intermediate"
            ]

        if mag < 7.4:
            return C["tau_l2l"]

        if mag >= 8.0:
            return max(
                C["tau_l2l"],
                upper,
            )

        sigma_mu = self._linear_interpolation(
            mag,
            7.4,
            8.0,
            lower,
            upper,
        )

        return max(
            C["tau_l2l"],
            sigma_mu,
        )

    @classmethod
    def _pseudo_depth(
        cls,
        hypo_depth: float,
    ) -> float:
        """
        Return depth-dependent pseudo-depth h in km.
        """
        if hypo_depth <= 10.0:
            return cls._H_SHALLOW

        if hypo_depth > 20.0:
            return cls._H_DEEP

        return cls._H_INTERMEDIATE

    def _load_coefficient_set(
        self,
        coeff_set: str,
    ) -> dict[str, dict[str, float]]:
        """
        Load an auxiliary coefficient set from the model JSON file.
        """
        path = _os.path.join(
            _os.path.dirname(__file__),
            "data",
            self._COEFF_FILE,
        )

        with open(
            path,
            encoding="utf-8",
        ) as file:
            raw = _json.load(file)[
                "coefficients"
            ][coeff_set]

        keys = list(raw["keys"])

        return {
            imt: dict(
                zip(
                    keys,
                    [
                        float(value)
                        for value in values
                    ],
                )
            )
            for imt, values in raw["type"].items()
        }

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

    @staticmethod
    def _boolean_site_array(
        value,
        name: str,
    ) -> _np.ndarray:
        """
        Return a boolean site array, accepting bool or 0/1 values.
        """
        array = _np.asarray(value)

        if array.ndim != 1:
            raise ValueError(
                f"{name} must be one-dimensional."
            )

        if array.dtype == bool:
            return array.astype(
                bool,
                copy=False,
            )

        numeric = _np.asarray(
            array,
            dtype=float,
        )

        if not _np.all(
            _np.isfinite(numeric)
        ):
            raise ValueError(
                f"{name} must contain finite values."
            )

        valid = (
            (numeric == 0.0)
            | (numeric == 1.0)
        )

        if not _np.all(valid):
            raise ValueError(
                f"{name} must contain only boolean "
                "or 0/1 values."
            )

        return numeric.astype(bool)

    @staticmethod
    def _region_array(
        value,
    ) -> _np.ndarray:
        """
        Return validated ESHM20 attenuation-region identifiers.
        """
        numeric = _np.asarray(
            value,
            dtype=float,
        )

        if numeric.ndim != 1:
            raise ValueError(
                "eshm20_region must be one-dimensional."
            )

        if not _np.all(
            _np.isfinite(numeric)
        ):
            raise ValueError(
                "eshm20_region must contain finite values."
            )

        rounded = _np.rint(numeric)

        if not _np.all(
            numeric == rounded
        ):
            raise ValueError(
                "eshm20_region must contain integers."
            )

        region = rounded.astype(int)

        if _np.any(
            (region < 0)
            | (region > 5)
        ):
            raise ValueError(
                "eshm20_region must be in the range 0-5."
            )

        return region

    @staticmethod
    def _normalize_adjustment_mapping(
        value,
        name: str,
    ) -> dict[str, float] | None:
        """
        Normalize an optional IMT-to-float mapping.
        """
        if value is None:
            return None

        if not isinstance(value, _Mapping):
            raise TypeError(
                f"{name} must be a mapping or None."
            )

        output = {}

        for imt, adjustment in value.items():
            key = str(imt).strip()

            if not key:
                raise ValueError(
                    f"{name} contains an empty IMT key."
                )

            number = float(adjustment)

            if not _np.isfinite(number):
                raise ValueError(
                    f"{name}[{key!r}] must be finite."
                )

            output[key] = number

        return output

    @staticmethod
    def _mapping_value(
        mapping: dict[str, float],
        imt: str,
        name: str,
    ) -> float:
        """
        Return one IMT-specific model adjustment.
        """
        try:
            return mapping[imt]
        except KeyError as exc:
            available = ", ".join(
                sorted(mapping)
            )
            raise KeyError(
                f"{name} has no value for {imt!r}. "
                f"Available: {available}."
            ) from exc

    @staticmethod
    def _linear_interpolation(
        x: float,
        x0: float,
        x1: float,
        y0: float,
        y1: float,
    ) -> float:
        """
        Return linear interpolation between two points.
        """
        return (
            y0
            + (
                y1 - y0
            )
            * (
                x - x0
            )
            / (
                x1 - x0
            )
        )


__all__ = [
    "KothaEtAl2020ESHM20",
]
