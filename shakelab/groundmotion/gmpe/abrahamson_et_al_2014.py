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
Abrahamson, Silva and Kamai (2014) NGA-West2 ground-motion model.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class AbrahamsonEtAl2014(GMPE):
    """
    Abrahamson, Silva and Kamai (2014), ASK14.

    Reference
    ---------
    Abrahamson, N. A., Silva, W. J., and Kamai, R. (2014).
    Summary of the ASK14 Ground Motion Relation for Active Crustal
    Regions. Earthquake Spectra, 30(3), 1025-1055.

    Notes
    -----
    This implementation reproduces the nominal ASK14 model, including
    hanging-wall, nonlinear site, depth-to-top, basin-depth and aleatory
    variability terms. Regional subclasses implement the CHN, JPN and
    TWN adjustments used by OpenQuake.

    The optional USGS basin scaling and CyberShake basin adjustments are
    deliberately not included here because they are external extensions
    to the original ASK14 formulation.

    PGA and SA are returned in g. PGV is returned in m/s. All medians and
    standard deviations are expressed in natural-log space.
    """

    REFERENCE_VELOCITY = 1180.0
    MAGNITUDE_TYPE = "Mw"

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
        "rake",
        "dip",
        "ztor",
        "width",
    }

    REQUIRES_SITE_PARAMETERS = {
        "vs30",
        "z1pt0",
        "vs30measured",
    }

    REQUIRES_DISTANCES = {
        "rrup",
        "rjb",
        "rx",
        "ry0",
    }

    _COEFF_FILE = "abrahamson_et_al_2014.json"
    _COEFF_SET = "default"

    _REGION = None

    _N = 1.5
    _M2 = 5.0
    _H1 = 0.25
    _H2 = 1.50
    _H3 = -0.75
    _METRES_PER_KM = 1000.0

    def __init__(
        self,
        sigma_mu_epsilon: float = 0.0,
    ):
        """
        Initialize ASK14.

        Parameters
        ----------
        sigma_mu_epsilon
            Number of epistemic standard deviations added to the median,
            following the USGS-2014 epistemic adjustment implemented by
            OpenQuake. The default value 0.0 leaves the median unchanged.
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
        Evaluate ASK14 for one rupture and many sites.
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
        ry0 = self._finite_site_array(
            sites.ry0,
            "ry0",
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

        self._validate_physical_inputs(
            rake,
            dip,
            ztor,
            width,
            rrup,
            rjb,
            ry0,
            vs30,
        )

        # Rjb is required by ASK14/OpenQuake even though the current
        # analytical terms are expressed through Rrup, Rx and Ry0.
        del rjb

        sa1180 = _np.exp(
            self._reference_rock_mean(
                C,
                imt,
                mag,
                rake,
                dip,
                ztor,
                width,
                rrup,
                rx,
                ry0,
            )
        )

        mean = (
            self._basic_term(
                C,
                mag,
                rrup,
            )
            + self._hanging_wall_term(
                C,
                mag,
                dip,
                ztor,
                width,
                rx,
                ry0,
            )
            + self._site_response_term(
                C,
                imt,
                vs30,
                sa1180,
            )
            + self._top_of_rupture_term(
                C,
                ztor,
            )
            + self._faulting_style_term(
                C,
                mag,
                rake,
            )
            + self._basin_term(
                C,
                vs30,
                z1pt0,
            )
            + self._regional_term(
                C,
                imt,
                vs30,
                rrup,
            )
            + self.sigma_mu_epsilon
            * self._epistemic_sigma(
                mag,
                rrup,
            )
        )

        if imt == "PGV":
            # OpenQuake retains PGV in cm/s; ShakeLab uses m/s.
            mean = mean - _np.log(100.0)

        tau = self._inter_event_std(
            C,
            mag,
            sa1180,
            vs30,
        )
        phi = self._intra_event_std(
            C,
            mag,
            sa1180,
            vs30,
            vs30measured,
            rrup,
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

    @classmethod
    def _basic_term(
        cls,
        C: dict,
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return the ASK14 basic magnitude-distance term.
        """
        if mag > 5.0:
            c4m = C["c4"]
        elif mag < 4.0:
            c4m = 1.0
        else:
            c4m = (
                C["c4"]
                - (
                    C["c4"] - 1.0
                )
                * (
                    5.0 - mag
                )
            )

        rval = _np.sqrt(
            rrup**2 + c4m**2
        )

        base = (
            C["a1"]
            + C["a17"] * rrup
        )

        if mag >= C["m1"]:
            base = (
                base
                + C["a5"]
                * (
                    mag - C["m1"]
                )
                + C["a8"]
                * (
                    8.5 - mag
                )**2
                + (
                    C["a2"]
                    + C["a3"]
                    * (
                        mag - C["m1"]
                    )
                )
                * _np.log(rval)
            )

        elif mag >= cls._M2:
            base = (
                base
                + C["a4"]
                * (
                    mag - C["m1"]
                )
                + C["a8"]
                * (
                    8.5 - mag
                )**2
                + (
                    C["a2"]
                    + C["a3"]
                    * (
                        mag - C["m1"]
                    )
                )
                * _np.log(rval)
            )

        else:
            base = (
                base
                + C["a4"]
                * (
                    cls._M2 - C["m1"]
                )
                + C["a8"]
                * (
                    8.5 - cls._M2
                )**2
                + C["a6"]
                * (
                    mag - cls._M2
                )
                + C["a7"]
                * (
                    mag - cls._M2
                )**2
                + (
                    C["a2"]
                    + C["a3"]
                    * (
                        cls._M2 - C["m1"]
                    )
                )
                * _np.log(rval)
            )

        return base

    @staticmethod
    def _faulting_style_term(
        C: dict,
        mag: float,
        rake: float,
    ) -> float:
        """
        Return reverse and normal style-of-faulting scaling.
        """
        scale = _np.clip(
            mag - 4.0,
            0.0,
            1.0,
        )

        reverse = (
            rake > 30.0
            and rake < 150.0
        )
        normal = (
            rake > -150.0
            and rake < -30.0
        )

        return (
            C["a11"]
            * scale
            * float(reverse)
            + C["a12"]
            * scale
            * float(normal)
        )

    @classmethod
    def _hanging_wall_term(
        cls,
        C: dict,
        mag: float,
        dip: float,
        ztor: float,
        width: float,
        rx: _np.ndarray,
        ry0: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return the ASK14 hanging-wall term.
        """
        if dip == 90.0:
            return _np.zeros_like(
                rx,
                dtype=float,
            )

        fhw = _np.zeros_like(
            rx,
            dtype=float,
        )
        fhw[rx > 0.0] = 1.0

        t1 = (
            60.0 / 45.0
            if dip <= 30.0
            else (
                90.0 - dip
            ) / 45.0
        )

        a2hw = 0.2

        if mag > 6.5:
            t2 = (
                1.0
                + a2hw
                * (
                    mag - 6.5
                )
            )
        elif mag <= 5.5:
            t2 = 0.0
        else:
            t2 = (
                1.0
                + a2hw
                * (
                    mag - 6.5
                )
                - (
                    1.0 - a2hw
                )
                * (
                    mag - 6.5
                )**2
            )

        r1 = (
            width
            * _np.cos(
                _np.radians(dip)
            )
        )
        r2 = 3.0 * r1

        t3 = _np.zeros_like(
            rx,
            dtype=float,
        )

        idx = rx < r1
        ratio = rx[idx] / r1

        t3[idx] = (
            cls._H1
            + cls._H2 * ratio
            + cls._H3 * ratio**2
        )

        idx = (
            (rx >= r1)
            & (rx <= r2)
        )

        t3[idx] = (
            1.0
            - (
                rx[idx] - r1
            )
            / (
                r2 - r1
            )
        )

        t4 = (
            1.0 - ztor**2 / 100.0
            if ztor <= 10.0
            else 0.0
        )

        t5 = _np.zeros_like(
            rx,
            dtype=float,
        )

        ry1 = (
            rx
            * _np.tan(
                _np.radians(20.0)
            )
        )

        delta = ry0 - ry1

        idx = delta <= 0.0
        t5[idx] = 1.0

        idx = (
            (delta > 0.0)
            & (delta < 5.0)
        )

        t5[idx] = (
            1.0
            - delta[idx] / 5.0
        )

        return (
            fhw
            * C["a13"]
            * t1
            * t2
            * t3
            * t4
            * t5
        )

    @staticmethod
    def _top_of_rupture_term(
        C: dict,
        ztor: float,
    ) -> float:
        """
        Return depth-to-top-of-rupture scaling.
        """
        return (
            C["a15"]
            * min(
                ztor / 20.0,
                1.0,
            )
        )

    @classmethod
    def _site_response_term(
        cls,
        C: dict,
        imt: str,
        vs30: _np.ndarray,
        sa1180: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return ASK14 nonlinear site response.
        """
        vs30star = cls._vs30_star(
            vs30,
            imt,
        )
        ratio = (
            vs30star / C["vlin"]
        )

        site = _np.zeros_like(
            vs30,
            dtype=float,
        )

        idx = vs30 >= C["vlin"]

        site[idx] = (
            C["a10"]
            + C["b"] * cls._N
        ) * _np.log(
            ratio[idx]
        )

        idx = vs30 < C["vlin"]

        site[idx] = (
            C["a10"]
            * _np.log(
                ratio[idx]
            )
            - C["b"]
            * _np.log(
                sa1180[idx] + C["c"]
            )
            + C["b"]
            * _np.log(
                sa1180[idx]
                + C["c"]
                * ratio[idx]**cls._N
            )
        )

        return site

    @classmethod
    def _reference_rock_mean(
        cls,
        C: dict,
        imt: str,
        mag: float,
        rake: float,
        dip: float,
        ztor: float,
        width: float,
        rrup: _np.ndarray,
        rx: _np.ndarray,
        ry0: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return median motion for Vs30 = 1180 m/s.
        """
        vs30 = _np.full_like(
            rrup,
            1180.0,
            dtype=float,
        )
        ref_iml = _np.zeros_like(
            rrup,
            dtype=float,
        )

        return (
            cls._basic_term(
                C,
                mag,
                rrup,
            )
            + cls._faulting_style_term(
                C,
                mag,
                rake,
            )
            + cls._site_response_term(
                C,
                imt,
                vs30,
                ref_iml,
            )
            + cls._hanging_wall_term(
                C,
                mag,
                dip,
                ztor,
                width,
                rx,
                ry0,
            )
            + cls._top_of_rupture_term(
                C,
                ztor,
            )
            + cls._regional_term_static(
                cls._REGION,
                C,
                imt,
                vs30,
                rrup,
            )
        )

    def _basin_term(
        self,
        C: dict,
        vs30: _np.ndarray,
        z1pt0: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return the original ASK14 basin-depth term.
        """
        z1ref = self._z1pt0_reference(
            vs30,
        )

        z10 = (
            z1pt0.copy()
            / self._METRES_PER_KM
        )

        missing = z10 < 0.0
        z10[missing] = z1ref[missing]

        factor = _np.log(
            (
                z10 + 0.01
            )
            / (
                z1ref + 0.01
            )
        )

        f2 = _np.interp(
            vs30,
            [
                0.0,
                150.0,
                250.0,
                400.0,
                700.0,
                1000.0,
                6000.0,
            ],
            [
                C["a43"],
                C["a43"],
                C["a44"],
                C["a45"],
                C["a46"],
                C["a46"],
                C["a46"],
            ],
        )

        return f2 * factor

    def _regional_term(
        self,
        C: dict,
        imt: str,
        vs30: _np.ndarray,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return regional ASK14 correction.
        """
        return self._regional_term_static(
            self._REGION,
            C,
            imt,
            vs30,
            rrup,
        )

    @classmethod
    def _regional_term_static(
        cls,
        region,
        C: dict,
        imt: str,
        vs30: _np.ndarray,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return regional correction for the supplied region.
        """
        if region == "TWN":
            vs30star = cls._vs30_star(
                vs30,
                imt,
            )

            return (
                C["a31"]
                * _np.log(
                    vs30star / C["vlin"]
                )
                + C["a25"] * rrup
            )

        if region == "CHN":
            return (
                C["a28"] * rrup
            )

        if region == "JPN":
            return (
                cls._linear_interp_extrapolate(
                    vs30,
                    [
                        150.0,
                        250.0,
                        350.0,
                        450.0,
                        600.0,
                        850.0,
                        1150.0,
                        2000.0,
                    ],
                    [
                        C["a36"],
                        C["a37"],
                        C["a38"],
                        C["a39"],
                        C["a40"],
                        C["a41"],
                        C["a42"],
                        C["a42"],
                    ],
                )
                + C["a29"] * rrup
            )

        return _np.zeros_like(
            rrup,
            dtype=float,
        )

    def _z1pt0_reference(
        self,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return regional reference Z1.0 in km.
        """
        if self._REGION == "JPN":
            return (
                _np.exp(
                    -5.23
                    / 2.0
                    * _np.log(
                        (
                            vs30**2
                            + 412.0**2
                        )
                        / (
                            1360.0**2
                            + 412.0**2
                        )
                    )
                )
                / self._METRES_PER_KM
            )

        return (
            _np.exp(
                -7.67
                / 4.0
                * _np.log(
                    (
                        vs30**4
                        + 610.0**4
                    )
                    / (
                        1360.0**4
                        + 610.0**4
                    )
                )
            )
            / self._METRES_PER_KM
        )

    @classmethod
    def _vs30_star(
        cls,
        vs30: _np.ndarray,
        imt: str,
    ) -> _np.ndarray:
        """
        Return capped Vs30* from ASK14 equations 8 and 9.
        """
        period = cls._imt_period(
            imt,
        )

        if period is None:
            v1 = 1500.0
        elif period <= 0.50:
            v1 = 1500.0
        elif period < 3.0:
            v1 = _np.exp(
                -0.35
                * _np.log(
                    period / 0.5
                )
                + _np.log(1500.0)
            )
        else:
            v1 = 800.0

        return _np.minimum(
            vs30,
            v1,
        )

    @classmethod
    def _derivative(
        cls,
        C: dict,
        sa1180: _np.ndarray,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return nonlinear amplification derivative.
        """
        derivative = _np.zeros_like(
            vs30,
            dtype=float,
        )

        idx = vs30 < C["vlin"]

        derivative[idx] = (
            C["b"]
            * sa1180[idx]
            * (
                -1.0
                / (
                    sa1180[idx]
                    + C["c"]
                )
                + 1.0
                / (
                    sa1180[idx]
                    + C["c"]
                    * (
                        vs30[idx]
                        / C["vlin"]
                    )**cls._N
                )
            )
        )

        return derivative

    @classmethod
    def _inter_event_std(
        cls,
        C: dict,
        mag: float,
        sa1180: _np.ndarray,
        vs30: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return ASK14 inter-event standard deviation.
        """
        if mag < 5.0:
            tau_al = C["s3"]
        elif mag >= 7.0:
            tau_al = C["s4"]
        else:
            tau_al = (
                C["s3"]
                + (
                    C["s4"]
                    - C["s3"]
                )
                / 2.0
                * (
                    mag - 5.0
                )
            )

        return (
            tau_al
            * (
                1.0
                + cls._derivative(
                    C,
                    sa1180,
                    vs30,
                )
            )
        )

    def _intra_event_std(
        self,
        C: dict,
        mag: float,
        sa1180: _np.ndarray,
        vs30: _np.ndarray,
        vs30measured: _np.ndarray,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return ASK14 intra-event standard deviation.
        """
        if self._REGION == "JPN":
            phi_al = _np.empty_like(
                rrup,
                dtype=float,
            )

            idx = rrup < 30.0
            phi_al[idx] = C["s5"]

            idx = (
                (rrup >= 30.0)
                & (rrup <= 80.0)
            )

            phi_al[idx] = (
                C["s5"]
                + (
                    C["s6"]
                    - C["s5"]
                )
                / 50.0
                * (
                    rrup[idx] - 30.0
                )
            )

            idx = rrup > 80.0
            phi_al[idx] = C["s6"]

        else:
            s1 = _np.full_like(
                rrup,
                C["s1e"],
                dtype=float,
            )
            s2 = _np.full_like(
                rrup,
                C["s2e"],
                dtype=float,
            )

            s1[vs30measured] = C["s1m"]
            s2[vs30measured] = C["s2m"]

            if mag < 4.0:
                phi_al = s1

            elif mag >= 6.0:
                phi_al = s2

            else:
                phi_al = (
                    s1
                    + (
                        s2 - s1
                    )
                    / 2.0
                    * (
                        mag - 4.0
                    )
                )

        derivative = self._derivative(
            C,
            sa1180,
            vs30,
        )

        phi_amp = _np.full_like(
            rrup,
            0.4,
            dtype=float,
        )

        idx = phi_al < phi_amp
        phi_amp[idx] = (
            0.99 * phi_al[idx]
        )

        phi_b = _np.sqrt(
            phi_al**2
            - phi_amp**2
        )

        return _np.sqrt(
            phi_b**2
            * (
                1.0 + derivative
            )**2
            + phi_amp**2
        )

    @staticmethod
    def _epistemic_sigma(
        mag: float,
        rrup: _np.ndarray,
    ) -> _np.ndarray:
        """
        Return the OpenQuake USGS-2014 epistemic adjustment.
        """
        n = 2.0

        if mag < 6.0 and mag >= 5.0:
            return _np.where(
                rrup <= 10.0,
                0.4 * _np.sqrt(n / 11.0),
                _np.where(
                    rrup < 30.0,
                    0.4 * _np.sqrt(n / 38.0),
                    0.4 * _np.sqrt(n / 94.0),
                ),
            )

        if mag < 7.0 and mag >= 6.0:
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
    def _linear_interp_extrapolate(
        x: _np.ndarray,
        xp,
        fp,
    ) -> _np.ndarray:
        """
        Return one-dimensional linear interpolation with extrapolation.
        """
        x = _np.asarray(
            x,
            dtype=float,
        )
        xp = _np.asarray(
            xp,
            dtype=float,
        )
        fp = _np.asarray(
            fp,
            dtype=float,
        )

        out = _np.interp(
            x,
            xp,
            fp,
        )

        left = x < xp[0]
        if _np.any(left):
            slope = (
                fp[1] - fp[0]
            ) / (
                xp[1] - xp[0]
            )
            out[left] = (
                fp[0]
                + slope
                * (
                    x[left] - xp[0]
                )
            )

        right = x > xp[-1]
        if _np.any(right):
            slope = (
                fp[-1] - fp[-2]
            ) / (
                xp[-1] - xp[-2]
            )
            out[right] = (
                fp[-1]
                + slope
                * (
                    x[right] - xp[-1]
                )
            )

        return out

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
        rrup: _np.ndarray,
        rjb: _np.ndarray,
        ry0: _np.ndarray,
        vs30: _np.ndarray,
    ) -> None:
        """
        Validate physical ASK14 inputs.
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

        if _np.any(rrup < 0.0):
            raise ValueError(
                "rrup must be non-negative."
            )

        if _np.any(rjb < 0.0):
            raise ValueError(
                "rjb must be non-negative."
            )

        if _np.any(ry0 < 0.0):
            raise ValueError(
                "ry0 must be non-negative."
            )

        if _np.any(vs30 <= 0.0):
            raise ValueError(
                "vs30 must be strictly positive."
            )

        if _np.any(vs30 > 6000.0):
            raise ValueError(
                "vs30 must not exceed 6000 m/s."
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


class AbrahamsonEtAl2014RegCHN(AbrahamsonEtAl2014):
    """
    ASK14 regionalization for China.
    """

    _REGION = "CHN"


class AbrahamsonEtAl2014RegJPN(AbrahamsonEtAl2014):
    """
    ASK14 regionalization for Japan.
    """

    _REGION = "JPN"


class AbrahamsonEtAl2014RegTWN(AbrahamsonEtAl2014):
    """
    ASK14 regionalization for Taiwan.
    """

    _REGION = "TWN"


__all__ = [
    "AbrahamsonEtAl2014",
    "AbrahamsonEtAl2014RegCHN",
    "AbrahamsonEtAl2014RegJPN",
    "AbrahamsonEtAl2014RegTWN",
]
