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
Campbell and Bozorgnia (2014) NGA-West2 ground-motion models.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class CampbellBozorgnia2014(GMPE):
    """
    Campbell and Bozorgnia (2014), CB14.

    Reference
    ---------
    Campbell, K. W., and Bozorgnia, Y. (2014).
    NGA-West2 Ground Motion Model for the Average Horizontal Components
    of PGA, PGV, and 5%-Damped Linear Acceleration Response Spectra.
    Earthquake Spectra, 30(3), 1087-1115.

    Notes
    -----
    This implementation includes the nominal CB14 formulation, nonlinear
    shallow-site response, hanging-wall effects, basin-depth effects,
    hypocentral-depth scaling, fault-dip scaling and the complete CB14
    aleatory variability model.

    The original CB14 California basin/site formulation is used unless a
    JapanSite subclass is selected. HighQ and LowQ subclasses use the
    corresponding OpenQuake ``Dc20`` attenuation adjustments.

    Optional OpenQuake estimation of ztor, width and hypocentral depth is
    deliberately not included: ShakeLab requires those rupture parameters
    explicitly rather than deriving them inside the GMPE. USGS basin
    scaling and CyberShake basin adjustments are also excluded because
    they are external extensions to the original CB14 formulation.

    PGA and SA are returned in g. PGV is returned in m/s. Medians and
    standard deviations are expressed in natural-log space.
    """

    REFERENCE_VELOCITY = 1100.0
    MAGNITUDE_TYPE = "Mw"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
        "rake",
        "dip",
        "ztor",
        "width",
        "hypo_depth",
    }

    REQUIRES_SITE_PARAMETERS = {
        "vs30",
        "z2pt5",
    }

    REQUIRES_DISTANCES = {
        "rrup",
        "rjb",
        "rx",
    }

    _COEFF_FILE = "campbell_bozorgnia_2014.json"
    _COEFF_SET = "default"

    _JAPAN_SITE = False

    _H4 = 1.0
    _C = 1.88
    _N = 1.18

    def __init__(
        self,
        sigma_mu_epsilon: float = 0.0,
    ):
        """
        Initialize CB14.

        Parameters
        ----------
        sigma_mu_epsilon
            Number of epistemic standard deviations added to the median
            using the same USGS-2014 adjustment adopted by OpenQuake.
        """
        epsilon = float(sigma_mu_epsilon)

        if not _np.isfinite(epsilon):
            raise ValueError(
                "sigma_mu_epsilon must be finite."
            )

        self.sigma_mu_epsilon = epsilon

        super().__init__()

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate CB14 for one rupture and many sites.
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
        dip = self._finite_scalar(
            rupture.dip,
            "dip",
        )
        ztor = self._finite_scalar(
            rupture.ztor,
            "ztor",
        )
        width = self._finite_scalar(
            rupture.width,
            "width",
        )
        hypo_depth = self._finite_scalar(
            rupture.hypo_depth,
            "hypo_depth",
        )

        rrup = self._finite_site_array(
            sites.rrup,
            "rrup",
        )
        rjb = self._finite_site_array(
            sites.rjb,
            "rjb",
        )
        rx = self._finite_site_array(
            sites.rx,
            "rx",
        )
        vs30 = self._finite_site_array(
            sites.vs30,
            "vs30",
        )
        z2pt5 = self._finite_site_array(
            sites.z2pt5,
            "z2pt5",
        )

        self._validate_physical_inputs(
            rake,
            dip,
            ztor,
            width,
            hypo_depth,
            rrup,
            rjb,
            vs30,
            z2pt5,
        )

        C = self.get_coefficients(imt)
        C_pga = self.get_coefficients("PGA")

        pga1100 = _np.exp(
            self._mean_values(
                C_pga,
                mag,
                rake,
                dip,
                ztor,
                width,
                hypo_depth,
                rrup,
                rjb,
                rx,
                None,
                None,
                reference=True,
            )
        )

        mean = self._mean_values(
            C,
            mag,
            rake,
            dip,
            ztor,
            width,
            hypo_depth,
            rrup,
            rjb,
            rx,
            vs30,
            z2pt5,
            pga_rock=pga1100,
        )

        epistemic = (
            self.sigma_mu_epsilon
            * self._epistemic_sigma(
                mag,
                rrup,
            )
        )
        mean = mean + epistemic

        if (
            imt.startswith("SA-")
            and self._imt_period(imt) < 0.25
        ):
            pga = self._mean_values(
                C_pga,
                mag,
                rake,
                dip,
                ztor,
                width,
                hypo_depth,
                rrup,
                rjb,
                rx,
                vs30,
                z2pt5,
                pga_rock=pga1100,
            )

            idx = mean <= pga
            mean[idx] = pga[idx]

            # Reproduce the current OpenQuake CB14 ordering exactly:
            # epistemic adjustment is added again after the short-period
            # PGA floor operation.
            mean = mean + epistemic

        if imt == "PGV":
            # OpenQuake retains PGV in cm/s; ShakeLab uses m/s.
            mean = mean - _np.log(100.0)

        tau_pga_b = self._tau_ln_y(
            C_pga,
            mag,
        )
        phi_pga_b = _np.sqrt(
            self._phi_ln_y(
                C_pga,
                mag,
            )**2
            - C_pga["philnAF"]**2
        )

        sigma, tau, phi = self._stddevs(
            C,
            mag,
            vs30,
            pga1100,
            tau_pga_b,
            phi_pga_b,
        )

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    def _mean_values(
        self,
        C: dict,
        mag: float,
        rake: float,
        dip: float,
        ztor: float,
        width: float,
        hypo_depth: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        rx: _np.ndarray,
        vs30,
        z2pt5,
        pga_rock=None,
        reference: bool = False,
    ) -> _np.ndarray:
        """
        Return CB14 median before optional epistemic adjustment.
        """
        if reference:
            velocity = _np.full_like(
                rrup,
                1100.0,
                dtype=float,
            )
            pga_input = None
        else:
            velocity = vs30
            pga_input = pga_rock

        return (
            self._magnitude_term(
                C,
                mag,
            )
            + self._geometric_attenuation_term(
                C,
                mag,
                rrup,
            )
            + self._style_of_faulting_term(
                C,
                mag,
                rake,
            )
            + self._hanging_wall_term(
                C,
                mag,
                dip,
                ztor,
                width,
                rrup,
                rjb,
                rx,
            )
            + self._shallow_site_response_term(
                C,
                velocity,
                pga_input,
            )
            + self._basin_term(
                C,
                velocity,
                z2pt5,
                reference=reference,
            )
            + self._hypocentral_depth_term(
                C,
                mag,
                hypo_depth,
            )
            + self._fault_dip_term(
                C,
                mag,
                dip,
            )
            + self._anelastic_attenuation_term(
                C,
                rrup,
            )
        )

    @staticmethod
    def _magnitude_term(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return CB14 magnitude scaling.
        """
        result = (
            C["c0"]
            + C["c1"] * mag
        )

        if mag > 6.5:
            result += (
                C["c2"] * (mag - 4.5)
                + C["c3"] * (mag - 5.5)
                + C["c4"] * (mag - 6.5)
            )
        elif mag > 5.5:
            result += (
                C["c2"] * (mag - 4.5)
                + C["c3"] * (mag - 5.5)
            )
        elif mag > 4.5:
            result += (
                C["c2"] * (mag - 4.5)
            )

        return result

    @staticmethod
    def _geometric_attenuation_term(
        C: dict,
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return geometric attenuation.
        """
        return (
            C["c5"]
            + C["c6"] * mag
        ) * _np.log(
            _np.sqrt(
                rrup**2
                + C["c7"]**2
            )
        )

    @staticmethod
    def _style_of_faulting_term(
        C: dict,
        mag: float,
        rake: float,
    ) -> float:
        """
        Return style-of-faulting scaling.
        """
        reverse = (
            rake > 30.0
            and rake < 150.0
        )
        normal = (
            rake > -150.0
            and rake < -30.0
        )

        factor = mag - 4.5

        if mag <= 4.5:
            factor = 0.0
        elif mag > 5.5:
            factor = 1.0

        return (
            C["c8"] * float(reverse)
            + C["c9"] * float(normal)
        ) * factor

    @classmethod
    def _hanging_wall_term(
        cls,
        C: dict,
        mag: float,
        dip: float,
        ztor: float,
        width: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        rx: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return CB14 hanging-wall scaling.
        """
        r1 = (
            width
            * _np.cos(
                _np.radians(dip)
            )
        )
        r2 = (
            62.0 * mag
            - 350.0
        )

        f_rx = _np.zeros_like(
            rx,
            dtype=float,
        )

        idx = (
            (rx >= 0.0)
            & (rx < r1)
        )

        if _np.any(idx):
            ratio = rx[idx] / r1

            f_rx[idx] = (
                C["h1"]
                + C["h2"] * ratio
                + C["h3"] * ratio**2
            )

        idx = rx >= r1

        if _np.any(idx):
            drx = (
                rx[idx] - r1
            ) / (
                r2 - r1
            )

            values = (
                cls._H4
                + C["h5"] * drx
                + C["h6"] * drx**2
            )

            values[values < 0.0] = 0.0
            f_rx[idx] = values

        f_rrup = _np.ones_like(
            rrup,
            dtype=float,
        )

        idx = rrup > 0.0
        f_rrup[idx] = (
            rrup[idx] - rjb[idx]
        ) / rrup[idx]

        f_mag = (
            (mag - 5.5)
            * (
                1.0
                + C["a2"]
                * (
                    mag - 6.5
                )
            )
        )

        if mag < 5.5:
            f_mag = 0.0
        elif mag > 6.5:
            f_mag = (
                1.0
                + C["a2"]
                * (
                    mag - 6.5
                )
            )

        f_ztor = (
            1.0 - 0.06 * ztor
        )

        if ztor > 16.66:
            f_ztor = 0.0

        f_dip = (
            90.0 - dip
        ) / 45.0

        return (
            C["c10"]
            * f_rx
            * f_rrup
            * f_mag
            * f_ztor
            * f_dip
        )

    def _shallow_site_response_term(
        self,
        C: dict,
        vs30: _np.ndarray,
        pga_rock,
    ) -> _np.ndarray:
        """
        Return shallow-site response.
        """
        ratio = (
            vs30 / C["k1"]
        )

        site_global = (
            C["c11"]
            * _np.log(ratio)
        )

        idx = vs30 > C["k1"]

        site_global[idx] += (
            C["k2"]
            * self._N
            * _np.log(
                ratio[idx]
            )
        )

        idx = ~idx

        if _np.any(idx):
            if pga_rock is None:
                rock = _np.zeros_like(
                    vs30,
                    dtype=float,
                )
            else:
                rock = pga_rock

            site_global[idx] += (
                C["k2"]
                * (
                    _np.log(
                        rock[idx]
                        + self._C
                        * ratio[idx]**self._N
                    )
                    - _np.log(
                        rock[idx]
                        + self._C
                    )
                )
            )

        if not self._JAPAN_SITE:
            return site_global

        site_japan = (
            C["c13"]
            + C["k2"] * self._N
        ) * _np.log(ratio)

        idx = vs30 <= 200.0

        site_japan[idx] += (
            C["c12"]
            + C["k2"] * self._N
        ) * (
            _np.log(
                ratio[idx]
            )
            - _np.log(
                200.0 / C["k1"]
            )
        )

        return (
            site_global
            + site_japan
        )

    def _basin_term(
        self,
        C: dict,
        vs30: _np.ndarray,
        z2pt5,
        reference: bool = False,
    ) -> _np.ndarray:
        """
        Return original CB14 basin-depth term.
        """
        sj = self._JAPAN_SITE

        z_ref = (
            self._z2pt5_reference(
                1100.0,
            )
            * _np.ones_like(
                vs30,
                dtype=float,
            )
        )

        if reference:
            depth = z_ref
        else:
            depth = z2pt5.copy()

            missing = depth == -999.0

            if _np.any(missing):
                depth[missing] = (
                    self._z2pt5_reference(
                        vs30[missing],
                    )
                )

        out = _np.zeros_like(
            vs30,
            dtype=float,
        )

        idx = depth < 1.0

        out[idx] = (
            C["c14"]
            + C["c15"] * float(sj)
        ) * (
            depth[idx] - 1.0
        )

        idx = depth > 3.0

        out[idx] = (
            C["c16"]
            * C["k3"]
            * _np.exp(-0.75)
            * (
                1.0
                - _np.exp(
                    -0.25
                    * (
                        depth[idx] - 3.0
                    )
                )
            )
        )

        return out

    def _z2pt5_reference(
        self,
        vs30,
    ):
        """
        Return CB14 reference Z2.5 in km.
        """
        if self._JAPAN_SITE:
            return _np.exp(
                5.359
                - 1.102
                * _np.log(vs30)
            )

        return _np.exp(
            7.089
            - 1.144
            * _np.log(vs30)
        )

    @staticmethod
    def _hypocentral_depth_term(
        C: dict,
        mag: float,
        hypo_depth: float,
    ) -> float:
        """
        Return hypocentral-depth scaling.
        """
        depth_factor = _np.clip(
            hypo_depth - 7.0,
            0.0,
            13.0,
        )

        mag_factor = (
            C["c17"]
            + (
                C["c18"] - C["c17"]
            )
            * (
                mag - 5.5
            )
        )

        if mag <= 5.5:
            mag_factor = C["c17"]
        elif mag > 6.5:
            mag_factor = C["c18"]

        return (
            depth_factor
            * mag_factor
        )

    @staticmethod
    def _fault_dip_term(
        C: dict,
        mag: float,
        dip: float,
    ) -> float:
        """
        Return fault-dip scaling.
        """
        result = (
            C["c19"]
            * (
                5.5 - mag
            )
            * dip
        )

        if mag < 4.5:
            result = (
                C["c19"] * dip
            )
        elif mag > 5.5:
            result = 0.0

        return result

    @staticmethod
    def _anelastic_attenuation_term(
        C: dict,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return anelastic attenuation beyond 80 km.
        """
        out = _np.zeros_like(
            rrup,
            dtype=float,
        )

        idx = rrup >= 80.0

        out[idx] = (
            C["c20"] + C["Dc20"]
        ) * (
            rrup[idx] - 80.0
        )

        return out

    @classmethod
    def _alpha(
        cls,
        C: dict,
        vs30: _np.ndarray,
        pga_rock: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return linearized site-amplification derivative.
        """
        alpha = _np.zeros_like(
            pga_rock,
            dtype=float,
        )

        idx = vs30 < C["k1"]

        if _np.any(idx):
            ratio = (
                vs30[idx] / C["k1"]
            )

            af1 = (
                pga_rock[idx]
                + cls._C
                * ratio**cls._N
            )
            af2 = (
                pga_rock[idx]
                + cls._C
            )

            alpha[idx] = (
                C["k2"]
                * pga_rock[idx]
                * (
                    1.0 / af1
                    - 1.0 / af2
                )
            )

        return alpha

    @staticmethod
    def _tau_ln_y(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return basement-rock inter-event standard deviation.
        """
        result = (
            C["tau2"]
            + (
                C["tau1"] - C["tau2"]
            )
            * (
                5.5 - mag
            )
        )

        if mag <= 4.5:
            result = C["tau1"]
        elif mag >= 5.5:
            result = C["tau2"]

        return result

    @staticmethod
    def _phi_ln_y(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return basement-rock intra-event standard deviation.
        """
        result = (
            C["phi2"]
            + (
                C["phi1"] - C["phi2"]
            )
            * (
                5.5 - mag
            )
        )

        if mag <= 4.5:
            result = C["phi1"]
        elif mag >= 5.5:
            result = C["phi2"]

        return result

    @classmethod
    def _stddevs(
        cls,
        C: dict,
        mag: float,
        vs30: _np.ndarray,
        pga1100: _np.ndarray,
        tau_pga_b: float,
        phi_pga_b: float,
    ):
        """
        Return total, inter-event and intra-event variability.
        """
        tau_y_b = cls._tau_ln_y(
            C,
            mag,
        )
        phi_y_b = _np.sqrt(
            cls._phi_ln_y(
                C,
                mag,
            )**2
            - C["philnAF"]**2
        )

        alpha = cls._alpha(
            C,
            vs30,
            pga1100,
        )

        tau = _np.sqrt(
            tau_y_b**2
            + alpha**2 * tau_pga_b**2
            + 2.0
            * alpha
            * C["rholny"]
            * tau_y_b
            * tau_pga_b
        )

        phi = _np.sqrt(
            phi_y_b**2
            + C["philnAF"]**2
            + alpha**2 * phi_pga_b**2
            + 2.0
            * alpha
            * C["rholny"]
            * phi_y_b
            * phi_pga_b
        )

        sigma = _np.sqrt(
            tau**2 + phi**2
        )

        return sigma, tau, phi

    @staticmethod
    def _epistemic_sigma(
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return OpenQuake's USGS-2014 epistemic adjustment.
        """
        n = 2.0

        if 5.0 <= mag < 6.0:
            return _np.where(
                rrup <= 10.0,
                0.4 * _np.sqrt(n / 11.0),
                _np.where(
                    rrup < 30.0,
                    0.4 * _np.sqrt(n / 38.0),
                    0.4 * _np.sqrt(n / 94.0),
                ),
            )

        if 6.0 <= mag < 7.0:
            return _np.where(
                rrup <= 10.0,
                0.4 * _np.sqrt(n / 2.0),
                _np.where(
                    rrup < 30.0,
                    0.4 * _np.sqrt(n / 7.0),
                    0.4 * _np.sqrt(n / 13.0),
                ),
            )

        return _np.where(
            rrup <= 10.0,
            0.4 * _np.sqrt(n / 2.0),
            _np.where(
                rrup < 30.0,
                0.4 * _np.sqrt(n / 2.0),
                0.4 * _np.sqrt(n / 4.0),
            ),
        )

    @staticmethod
    def _imt_period(
        imt: str,
    ):
        """
        Return SA period or None for PGA/PGV.
        """
        if imt.startswith("SA-"):
            return float(
                imt.split(
                    "-",
                    1,
                )[1]
            )

        return None

    @staticmethod
    def _validate_physical_inputs(
        rake: float,
        dip: float,
        ztor: float,
        width: float,
        hypo_depth: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        vs30: _np.ndarray,
        z2pt5: _np.ndarray,
    ) -> None:
        """
        Validate physical CB14 inputs.
        """
        if rake < -180.0 or rake > 180.0:
            raise ValueError(
                "rake must be within [-180, 180] degrees."
            )

        if dip <= 0.0 or dip > 90.0:
            raise ValueError(
                "dip must be within (0, 90] degrees."
            )

        if ztor < 0.0:
            raise ValueError(
                "ztor must be non-negative."
            )

        if width <= 0.0:
            raise ValueError(
                "width must be strictly positive."
            )

        if hypo_depth < 0.0:
            raise ValueError(
                "hypo_depth must be non-negative."
            )

        if _np.any(rrup < 0.0):
            raise ValueError(
                "rrup must be non-negative."
            )

        if _np.any(rjb < 0.0):
            raise ValueError(
                "rjb must be non-negative."
            )

        if _np.any(vs30 <= 0.0):
            raise ValueError(
                "vs30 must be strictly positive."
            )

        valid_depth = (
            (z2pt5 >= 0.0)
            | (z2pt5 == -999.0)
        )

        if not _np.all(valid_depth):
            raise ValueError(
                "z2pt5 must be non-negative or -999 if unavailable."
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


class CampbellBozorgnia2014HighQ(CampbellBozorgnia2014):
    """
    CB14 high-Q attenuation variant.
    """

    _COEFF_SET = "high_q"


class CampbellBozorgnia2014LowQ(CampbellBozorgnia2014):
    """
    CB14 low-Q attenuation variant.
    """

    _COEFF_SET = "low_q"


class CampbellBozorgnia2014JapanSite(CampbellBozorgnia2014):
    """
    CB14 with the Japan shallow-site and basin model.
    """

    _JAPAN_SITE = True


class CampbellBozorgnia2014HighQJapanSite(
    CampbellBozorgnia2014JapanSite,
):
    """
    CB14 high-Q attenuation with the Japan site model.
    """

    _COEFF_SET = "high_q"


class CampbellBozorgnia2014LowQJapanSite(
    CampbellBozorgnia2014JapanSite,
):
    """
    CB14 low-Q attenuation with the Japan site model.
    """

    _COEFF_SET = "low_q"


__all__ = [
    "CampbellBozorgnia2014",
    "CampbellBozorgnia2014HighQ",
    "CampbellBozorgnia2014LowQ",
    "CampbellBozorgnia2014JapanSite",
    "CampbellBozorgnia2014HighQJapanSite",
    "CampbellBozorgnia2014LowQJapanSite",
]
