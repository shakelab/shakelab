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
Chiou and Youngs (2014) NGA-West2 ground-motion models.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class ChiouYoungs2014(GMPE):
    """
    Chiou and Youngs (2014), CY14.

    Reference
    ---------
    Chiou, B. S.-J., and Youngs, R. R. (2014).
    Update of the Chiou and Youngs NGA Model for the Average Horizontal
    Component of Peak Ground Motion and Response Spectra.
    Earthquake Spectra, 30(3), 1117-1153.

    Notes
    -----
    This implementation reproduces the original CY14 formulation used by
    OpenQuake, including source scaling, hanging-wall effects, nonlinear
    site response, basin-depth effects, regional attenuation adjustments
    and the nonlinear aleatory variability model.

    Regional subclasses provide the Japan, Italy and Wenchuan variants.
    The near-fault subclass explicitly requires ``rcdpp`` and activates
    the CY14 directivity term.

    Later OpenQuake extensions such as Boore et al. (2022) host-to-target
    source/path corrections, long-period delta-c1 adjustments, USGS basin
    scaling, CyberShake basin corrections and ACME adaptations are not
    included here.

    PGA and SA are returned in g. PGV is returned in m/s. Medians and
    standard deviations are expressed in natural-log space.
    """

    REFERENCE_VELOCITY = 1130.0
    MAGNITUDE_TYPE = "Mw"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
        "rake",
        "dip",
        "ztor",
    }

    REQUIRES_SITE_PARAMETERS = {
        "vs30",
        "vs30measured",
        "z1pt0",
    }

    REQUIRES_DISTANCES = {
        "rrup",
        "rjb",
        "rx",
    }

    _COEFF_FILE = "chiou_youngs_2014.json"
    _COEFF_SET = "default"

    _REGION = "CAL"
    _USE_DIRECTIVITY = False

    _C2 = 1.06
    _C4 = -2.1
    _C4A = -0.5
    _CRB = 50.0
    _C11 = 0.0
    _PHI6 = 300.0
    _PHI6_JPN = 800.0

    def __init__(
        self,
        sigma_mu_epsilon: float = 0.0,
        use_hanging_wall: bool = True,
    ):
        """
        Initialize CY14.

        Parameters
        ----------
        sigma_mu_epsilon
            Number of epistemic standard deviations added to the median,
            following the OpenQuake USGS-2014 epistemic adjustment.
        use_hanging_wall
            If True, include the CY14 hanging-wall term.
        """
        epsilon = float(sigma_mu_epsilon)

        if not _np.isfinite(epsilon):
            raise ValueError(
                "sigma_mu_epsilon must be finite."
            )

        if not isinstance(
            use_hanging_wall,
            bool,
        ):
            raise TypeError(
                "use_hanging_wall must be a boolean."
            )

        self.sigma_mu_epsilon = epsilon
        self.use_hanging_wall = use_hanging_wall

        super().__init__()

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate CY14 for one rupture and many sites.
        """
        self.validate_contexts(
            rupture,
            sites,
        )

        C = self.get_coefficients(imt)
        C_pga = self.get_coefficients("PGA")

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
        z1pt0 = self._finite_site_array(
            sites.z1pt0,
            "z1pt0",
        )
        vs30measured = self._boolean_site_array(
            sites.vs30measured,
            "vs30measured",
        )

        rcdpp = None

        if self._USE_DIRECTIVITY:
            rcdpp = self._finite_site_array(
                sites.rcdpp,
                "rcdpp",
            )

        self._validate_physical_inputs(
            rake,
            dip,
            ztor,
            rrup,
            rjb,
            vs30,
            z1pt0,
        )

        pga_mean, _, _, _ = self._mean_stddevs(
            C_pga,
            "PGA",
            mag,
            rake,
            dip,
            ztor,
            rrup,
            rjb,
            rx,
            vs30,
            z1pt0,
            vs30measured,
            rcdpp,
        )

        if imt == "PGA":
            mean = pga_mean
            _, sigma, tau, phi = self._mean_stddevs(
                C,
                imt,
                mag,
                rake,
                dip,
                ztor,
                rrup,
                rjb,
                rx,
                vs30,
                z1pt0,
                vs30measured,
                rcdpp,
            )
        else:
            mean, sigma, tau, phi = self._mean_stddevs(
                C,
                imt,
                mag,
                rake,
                dip,
                ztor,
                rrup,
                rjb,
                rx,
                vs30,
                z1pt0,
                vs30measured,
                rcdpp,
            )

            if (
                imt.startswith("SA-")
                and self._imt_period(imt) <= 0.3
            ):
                mean = _np.maximum(
                    mean,
                    pga_mean,
                )

        mean = (
            mean
            + self.sigma_mu_epsilon
            * self._epistemic_sigma(
                mag,
                rrup,
            )
        )

        if imt == "PGV":
            # OpenQuake retains PGV in cm/s; ShakeLab uses m/s.
            mean = mean - _np.log(100.0)

        return GroundMotionResult(
            mean=mean,
            sigma=sigma,
            tau=tau,
            phi=phi,
        )

    def _mean_stddevs(
        self,
        C: dict,
        imt: str,
        mag: float,
        rake: float,
        dip: float,
        ztor: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        rx: _np.ndarray,
        vs30: _np.ndarray,
        z1pt0: _np.ndarray,
        vs30measured: _np.ndarray,
        rcdpp,
    ):
        """
        Return median and standard deviations before epistemic adjustment.
        """
        ln_y_ref = self._reference_rock_mean(
            C,
            mag,
            rake,
            dip,
            ztor,
            rrup,
            rjb,
            rx,
            rcdpp,
        )
        y_ref = _np.exp(
            ln_y_ref
        )

        linear = self._linear_site_term(
            C,
            vs30,
        )
        nonlinear, scaling = self._nonlinear_site_term(
            C,
            vs30,
            y_ref,
        )
        basin = self._basin_term(
            C,
            vs30,
            z1pt0,
        )

        mean = (
            ln_y_ref
            + linear
            + nonlinear
            + basin
        )

        sigma, tau, phi = self._stddevs(
            C,
            mag,
            vs30,
            vs30measured,
            y_ref,
            scaling,
        )

        return mean, sigma, tau, phi

    def _reference_rock_mean(
        self,
        C: dict,
        mag: float,
        rake: float,
        dip: float,
        ztor: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        rx: _np.ndarray,
        rcdpp,
    ) -> _np.ndarray:
        """
        Return CY14 median on reference rock.
        """
        centered_ztor = self._centered_ztor(
            mag,
            rake,
            ztor,
        )

        out = (
            C["c1"]
            + self._magnitude_scaling(
                C,
                mag,
            )
            + self._source_scaling(
                C,
                mag,
                rake,
                dip,
                centered_ztor,
            )
            + self._geometric_spreading(
                C,
                mag,
                rrup,
            )
            + self._far_field_scaling(
                C,
                mag,
                rrup,
            )
        )

        if self.use_hanging_wall:
            out = (
                out
                + self._hanging_wall_term(
                    C,
                    dip,
                    ztor,
                    rrup,
                    rjb,
                    rx,
                )
            )

        if self._USE_DIRECTIVITY:
            out = (
                out
                + self._directivity_term(
                    C,
                    mag,
                    rrup,
                    rcdpp,
                )
            )

        return out

    @classmethod
    def _magnitude_scaling(
        cls,
        C: dict,
        mag: float,
    ) -> float:
        """
        Return CY14 magnitude scaling.
        """
        exponent = _np.log(
            1.0
            + _np.exp(
                C["cn"]
                * (
                    C["cm"] - mag
                )
            )
        )

        return (
            cls._C2
            * (
                mag - 6.0
            )
            + (
                cls._C2 - C["c3"]
            )
            / C["cn"]
            * exponent
        )

    @classmethod
    def _source_scaling(
        cls,
        C: dict,
        mag: float,
        rake: float,
        dip: float,
        centered_ztor: float,
    ) -> float:
        """
        Return style-of-faulting, Ztor and dip source scaling.
        """
        coshm = _np.cosh(
            2.0
            * max(
                mag - 4.5,
                0.0,
            )
        )

        reverse = (
            rake >= 30.0
            and rake <= 150.0
        )
        normal = (
            rake >= -120.0
            and rake <= -60.0
        )

        source = 0.0

        if reverse:
            source += (
                C["c1a"]
                + C["c1c"] / coshm
            )

        if normal:
            source += (
                C["c1b"]
                + C["c1d"] / coshm
            )

        source += (
            C["c7"]
            + C["c7b"] / coshm
        ) * centered_ztor

        source += (
            cls._C11
            + C["c11b"] / coshm
        ) * (
            _np.cos(
                _np.radians(dip)
            )**2
        )

        return source

    @classmethod
    def _geometric_spreading(
        cls,
        C: dict,
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return near-field geometric spreading.
        """
        return (
            cls._C4
            * _np.log(
                rrup
                + C["c5"]
                * _np.cosh(
                    C["c6"]
                    * max(
                        mag - C["chm"],
                        0.0,
                    )
                )
            )
        )

    def _far_field_scaling(
        self,
        C: dict,
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return regional far-field distance scaling.
        """
        geometric = (
            self._C4A - self._C4
        ) * _np.log(
            _np.sqrt(
                rrup**2
                + self._CRB**2
            )
        )

        gamma = (
            C["cg1"]
            + C["cg2"]
            / _np.cosh(
                max(
                    mag - C["cg3"],
                    0.0,
                )
            )
        )

        path = (
            gamma * rrup
        )

        if (
            self._REGION in {"JPN", "ITA"}
            and mag > 6.0
            and mag < 6.9
        ):
            path = (
                path * C["gjpit"]
            )

        if self._REGION == "WEN":
            path = (
                path * C["gwn"]
            )

        return geometric + path

    @staticmethod
    def _hanging_wall_term(
        C: dict,
        dip: float,
        ztor: float,
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        rx: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return CY14 hanging-wall scaling.
        """
        out = _np.zeros_like(
            rrup,
            dtype=float,
        )

        idx = rx >= 0.0

        if _np.any(idx):
            distance = (
                1.0
                - _np.sqrt(
                    rjb[idx]**2
                    + ztor**2
                )
                / (
                    rrup[idx] + 1.0
                )
            )

            distance *= (
                C["c9a"]
                + (
                    1.0 - C["c9a"]
                )
                * _np.tanh(
                    rx[idx] / C["c9b"]
                )
            )

            out[idx] = (
                C["c9"]
                * _np.cos(
                    _np.radians(dip)
                )
                * distance
            )

        return out

    @staticmethod
    def _directivity_term(
        C: dict,
        mag: float,
        rrup: _np.ndarray,
        rcdpp: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return CY14 near-fault directivity term.
        """
        magnitude = _np.clip(
            (
                mag - 5.5
            ) / 0.8,
            0.0,
            1.0,
        )

        distance = (
            1.0
            - _np.maximum(
                rrup - 40.0,
                0.0,
            ) / 30.0
        )
        distance = _np.maximum(
            distance,
            0.0,
        )

        return (
            C["c8"]
            * distance
            * _np.exp(
                -C["c8a"]
                * (
                    mag - C["c8b"]
                )**2
            )
            * magnitude
            * rcdpp
        )

    def _linear_site_term(
        self,
        C: dict,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return linear site amplification.
        """
        coefficient = (
            C["phi1jp"]
            if self._REGION == "JPN"
            else C["phi1"]
        )

        return (
            coefficient
            * _np.minimum(
                _np.log(
                    vs30
                    / self.REFERENCE_VELOCITY
                ),
                0.0,
            )
        )

    def _nonlinear_site_term(
        self,
        C: dict,
        vs30: _np.ndarray,
        y_ref: _np.ndarray,
    ):
        """
        Return nonlinear site term and nonlinear scaling factor.
        """
        velocity = _np.minimum(
            vs30,
            self.REFERENCE_VELOCITY,
        )

        scaling = (
            C["phi2"]
            * (
                _np.exp(
                    C["phi3"]
                    * (
                        velocity - 360.0
                    )
                )
                - _np.exp(
                    C["phi3"]
                    * (
                        self.REFERENCE_VELOCITY
                        - 360.0
                    )
                )
            )
        )

        nonlinear = (
            _np.log(
                (
                    y_ref
                    + C["phi4"]
                )
                / C["phi4"]
            )
            * scaling
        )

        return nonlinear, scaling

    def _basin_term(
        self,
        C: dict,
        vs30: _np.ndarray,
        z1pt0: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return original CY14 basin-depth term.
        """
        depth = z1pt0.copy()

        missing = (
            depth == -999.0
        )

        if _np.any(missing):
            depth[missing] = (
                self._z1pt0_reference(
                    vs30[missing],
                )
            )

        centered = (
            depth
            - self._z1pt0_reference(
                vs30,
            )
        )

        centered[depth <= 0.0] = 0.0

        if self._REGION == "JPN":
            return (
                C["phi5jp"]
                * (
                    1.0
                    - _np.exp(
                        -centered
                        / self._PHI6_JPN
                    )
                )
            )

        return (
            C["phi5"]
            * (
                1.0
                - _np.exp(
                    -centered
                    / self._PHI6
                )
            )
        )

    def _z1pt0_reference(
        self,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return regional reference Z1.0 in metres.
        """
        if self._REGION == "JPN":
            mean = (
                -5.23
                / 2.0
                * _np.log(
                    (
                        vs30**2
                        + 412.39**2
                    )
                    / (
                        1360.0**2
                        + 412.39**2
                    )
                )
            )
        else:
            mean = (
                -7.15
                / 4.0
                * _np.log(
                    (
                        vs30**4
                        + 570.94**4
                    )
                    / (
                        1360.0**4
                        + 570.94**4
                    )
                )
            )

        return _np.exp(mean)

    @staticmethod
    def _centered_ztor(
        mag: float,
        rake: float,
        ztor: float,
    ) -> float:
        """
        Return Ztor centered on the CY14 magnitude-dependent mean.
        """
        mean_ztor = max(
            2.673
            - 1.136
            * max(
                mag - 4.970,
                0.0,
            ),
            0.0,
        )**2

        if (
            rake >= 30.0
            and rake <= 150.0
        ):
            mean_ztor = max(
                2.704
                - 1.226
                * max(
                    mag - 5.849,
                    0.0,
                ),
                0.0,
            )**2

        return ztor - mean_ztor

    @staticmethod
    def _tau(
        C: dict,
        mag: float,
    ) -> float:
        """
        Return between-event variability before site propagation.
        """
        magnitude = _np.clip(
            mag - 5.0,
            0.0,
            1.5,
        )

        return (
            C["tau1"]
            + (
                C["tau2"]
                - C["tau1"]
            )
            / 1.5
            * magnitude
        )

    @staticmethod
    def _phi(
        C: dict,
        mag: float,
        vs30measured: _np.ndarray,
        nonlinear: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return within-event variability after nonlinear propagation.
        """
        base = _np.full_like(
            nonlinear,
            C["sig3"],
            dtype=float,
        )
        base[vs30measured] = 0.7

        site_factor = _np.sqrt(
            base
            + (
                1.0 + nonlinear
            )**2
        )

        magnitude = (
            C["sig1"]
            + (
                C["sig2"]
                - C["sig1"]
            )
            * _np.clip(
                mag - 5.0,
                0.0,
                1.5,
            )
            / 1.5
        )

        return (
            magnitude * site_factor
        )

    @classmethod
    def _stddevs(
        cls,
        C: dict,
        mag: float,
        vs30: _np.ndarray,
        vs30measured: _np.ndarray,
        y_ref: _np.ndarray,
        nonlinear_scaling: _np.ndarray,
    ):
        """
        Return total, inter-event and intra-event variability.
        """
        del vs30

        nonlinear = (
            nonlinear_scaling
            * (
                y_ref
                / (
                    y_ref + C["phi4"]
                )
            )
        )

        tau_base = cls._tau(
            C,
            mag,
        )
        phi = cls._phi(
            C,
            mag,
            vs30measured,
            nonlinear,
        )

        tau = _np.abs(
            (
                1.0 + nonlinear
            )
            * tau_base
        )

        sigma = _np.sqrt(
            (
                1.0 + nonlinear
            )**2
            * tau_base**2
            + phi**2
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
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        vs30: _np.ndarray,
        z1pt0: _np.ndarray,
    ) -> None:
        """
        Validate physical CY14 inputs.
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
            (z1pt0 >= 0.0)
            | (z1pt0 == -999.0)
        )

        if not _np.all(valid_depth):
            raise ValueError(
                "z1pt0 must be non-negative or -999 if unavailable."
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

    @staticmethod
    def _boolean_site_array(
        value,
        name: str,
    ) -> _np.ndarray:
        """
        Return a boolean site array, accepting bool or numeric 0/1.
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
            value,
            dtype=float,
        )

        if not _np.all(
            _np.isfinite(numeric)
        ):
            raise ValueError(
                f"{name} must contain finite values."
            )

        if not _np.all(
            (numeric == 0.0)
            | (numeric == 1.0)
        ):
            raise ValueError(
                f"{name} must contain only boolean or 0/1 values."
            )

        return numeric.astype(bool)


class ChiouYoungs2014Japan(ChiouYoungs2014):
    """
    CY14 Japan regional attenuation and site/basin model.
    """

    _REGION = "JPN"


class ChiouYoungs2014Italy(ChiouYoungs2014):
    """
    CY14 Italy far-field attenuation with California site model.
    """

    _REGION = "ITA"


class ChiouYoungs2014Wenchuan(ChiouYoungs2014):
    """
    CY14 Wenchuan far-field attenuation variant.

    This regional adjustment was calibrated for the Mw 7.9 Wenchuan
    earthquake and should be used outside that case with caution.
    """

    _REGION = "WEN"


class ChiouYoungs2014NearFaultEffect(ChiouYoungs2014):
    """
    CY14 with the near-fault directivity prediction term.
    """

    REQUIRES_DISTANCES = {
        "rrup",
        "rjb",
        "rx",
        "rcdpp",
    }

    _USE_DIRECTIVITY = True


__all__ = [
    "ChiouYoungs2014",
    "ChiouYoungs2014Japan",
    "ChiouYoungs2014Italy",
    "ChiouYoungs2014Wenchuan",
    "ChiouYoungs2014NearFaultEffect",
]
