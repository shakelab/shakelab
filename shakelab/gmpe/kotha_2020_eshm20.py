# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
Implementation of the Kotha et al. (2020) GMPE as adapted in ESHM20.

Kotha, S. R., Weatherill, G., Bindi, D., & Cotton, F. (2020).
A regionally-adaptable ground-motion model for shallow crustal
earthquakes in Europe. Bull. Earthq. Eng., 18, 4091–4125.
"""

from __future__ import annotations

import numpy as _np

from .base import GMPE as _GMPE


class KothaEtAl2020ESHM20(_GMPE):
    """
    Kotha et al. (2020) GMPE (ESHM20 variant).

    Median = magnitude piece + distance terms (with random effect on c3)
    + VS30 site term (observed vs inferred). Sigma is heteroscedastic:
      tau = hypot(tau_event_0, tau_l2l)
      phi = hypot(phi_0, phi_s2s) if ergodic else phi_0
      sigma = hypot(tau, phi)
    """

    # --- required by ShakeLab base -----------------------------------------
    REFERENCE_VELOCITY = 800.0
    DISTANCE_METRIC = "joyner-boore"
    MAGNITUDE_TYPE = "MW"

    _COEFF_FILE = "kotha_2020_eshm20.json"
    _COEFF_SET = "default"

    # Model constants (ESHM20)
    _MREF = 4.5
    _RREF = 30.0
    _MH = 5.7
    _H_D10 = 4.0
    _H_10D20 = 8.0
    _H_D20 = 12.0

    def ground_motion(  # noqa: D401
        self,
        imt: str,
        mag: float | _np.ndarray,
        dist: float | _np.ndarray,
        *,
        vs30: float | _np.ndarray,
        hypo_depth: float | _np.ndarray = 10.0,
        vs30measured: bool | _np.ndarray = True,
        c3_epsilon: float = 0.0,
        ergodic: bool = True,
    ) -> tuple[_np.ndarray, _np.ndarray]:
        """
        Compute (mean_ln, sigma_ln) for the requested IMT.

        Parameters
        ----------
        imt
            "PGA", "PGV", or "SA-<T>" (e.g., "SA-1.00").
        mag
            Moment magnitude Mw.
        dist
            Joyner–Boore distance Rjb [km].
        vs30
            Site VS30 [m/s].
        hypo_depth
            Hypocentral depth [km]; sets pseudo-depth h in {4, 8, 12} km.
        vs30measured
            True if VS30 is measured (observed), False if inferred.
        c3_epsilon
            Random effect on anelastic attenuation: c3 += eps * tau_c3.
        ergodic
            If True, include site-to-site in phi; else drop it.

        Returns
        -------
        mean_ln, sigma_ln
            Natural-log mean and total standard deviation.
        """
        C = self.get_coefficients(imt)

        mw = _np.asarray(mag, dtype=float)
        rjb = _np.asarray(dist, dtype=float)
        vs = _np.asarray(vs30, dtype=float)
        dep = _np.asarray(hypo_depth, dtype=float)
        is_obs = _np.asarray(vs30measured, dtype=bool)

        mw, rjb, vs, dep, is_obs = _np.broadcast_arrays(
            mw, rjb, vs, dep, is_obs
        )

        # --- magnitude term (hinge at MH)
        dm = mw - self._MH
        mean = _np.where(
            mw <= self._MH,
            C["e1"] + C["b1"] * dm + C["b2"] * dm**2,
            C["e1"] + C["b3"] * dm,
        )

        # --- effective distance with depth-dependent h
        h = _np.where(
            dep <= 10.0, self._H_D10,
            _np.where(dep > 20.0, self._H_D20, self._H_10D20),
        )
        rval = _np.sqrt(rjb**2 + h**2)
        rref = _np.sqrt(self._RREF**2 + h**2)

        # --- distance term with random effect on c3
        c3 = float(C["c3"]) + float(C.get("tau_c3", 0.0)) * float(c3_epsilon)
        mean += (C["c1"] + C["c2"] * (mw - self._MREF)) * _np.log(rval / rref)
        mean += (c3 * (rval - rref) / 100.0)

        # --- site term: observed vs inferred VS30
        vs_c = _np.clip(vs, None, 1100.0)
        site = _np.empty_like(vs_c, dtype=float)
        site[is_obs] = C["d0_obs"] + C["d1_obs"] * _np.log(vs_c[is_obs])
        site[~is_obs] = C["d0_inf"] + C["d1_inf"] * _np.log(vs_c[~is_obs])
        mean = mean + site

        # --- units: convert cm/s^2 → g for PGA/SA (in log domain)
        if imt.upper().startswith(("PGA", "SA")):
            g0 = 9.80665
            mean = mean - _np.log(100.0 * g0)

        # --- heteroscedastic sigma (hazardlib-ESHM20)
        tau = _np.hypot(C["tau_event_0"], C["tau_l2l"])
        phi0 = C["phi_0"]
        if ergodic:
            # prefer obs/inf-specific S2S if available; else generic phi_s2s
            if "phi_s2s_obs" in C and "phi_s2s_inf" in C:
                phi_s2s = _np.empty_like(mean)
                phi_s2s[is_obs] = C["phi_s2s_obs"]
                phi_s2s[~is_obs] = C["phi_s2s_inf"]
            else:
                phi_s2s = _np.full_like(mean, C.get("phi_s2s", 0.0))
            phi = _np.hypot(phi0, phi_s2s)
        else:
            phi = _np.full_like(mean, phi0)

        sigma = _np.hypot(tau, phi)
        return mean, sigma

