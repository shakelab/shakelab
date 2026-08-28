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
Massa et al. (2008) ground-motion prediction equations for Northern Italy.
"""

from __future__ import annotations

import numpy as _np

from shakelab.groundmotion.gmpe.base import (
    GMPE,
    GroundMotionResult,
    RuptureContext,
    SitesContext,
)


class _MassaEtAl2008Base(GMPE):
    """
    Shared implementation of Massa et al. (2008), equation (3).

    The model predicts the larger of the two horizontal components and
    uses epicentral distance and EC8 site class. The subclasses select
    either local magnitude or moment magnitude coefficients.

    Notes
    -----
    The calibration dataset covers Northern Italy, epicentral distances
    up to 100 km, ML from 3.5 to 6.3, and Mw up to 6.5.

    Site classes follow the paper:
    A: rock
    B: stiff soil
    C: soft soil

    Classes B and C share the same regression coefficient because the
    class-C subset was too small for a separate site term.

    PGA and SA are returned in g. PGV is returned in m/s. Medians and
    standard deviations are expressed in natural-log space.
    """

    REFERENCE_VELOCITY = None

    REQUIRES_RUPTURE_PARAMETERS = {
        "mag",
    }

    REQUIRES_SITE_PARAMETERS = {
        "ec8_class",
    }

    REQUIRES_DISTANCES = {
        "repi",
    }

    _COEFF_FILE = "massa_et_al_2008.json"

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate Massa et al. (2008), equation (3).
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
        repi = self._finite_site_array(
            sites.repi,
            "repi",
        )
        ec8_class = self._site_class_array(
            sites.ec8_class,
        )

        if _np.any(repi < 0.0):
            raise ValueError(
                "repi must be non-negative."
            )

        effective_distance = _np.sqrt(
            repi**2 + C["d"]**2
        )

        site_term = _np.where(
            ec8_class == "A",
            C["s1"],
            C["s2"],
        )

        mean_log10 = (
            C["a"]
            + C["b"] * mag
            + C["c"]
            * _np.log10(
                effective_distance
            )
            + site_term
        )

        mean = self._convert_mean(
            imt,
            mean_log10,
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

    @staticmethod
    def _convert_mean(
        imt: str,
        mean_log10,
    ):
        """
        Convert native log10 units to ShakeLab natural-log units.
        """
        if (
            imt == "PGA"
            or imt.startswith("SA-")
        ):
            # The paper tabulates horizontal PGA and SA directly in g.
            return (
                mean_log10
                * _np.log(10.0)
            )

        if imt == "PGV":
            # The paper tabulates horizontal PGV in cm/s.
            return (
                mean_log10
                * _np.log(10.0)
                - _np.log(100.0)
            )

        raise KeyError(
            f"Unsupported intensity measure type: {imt!r}."
        )

    @staticmethod
    def _site_class_array(
        value,
    ) -> _np.ndarray:
        """
        Return validated EC8 site classes A, B or C.
        """
        array = _np.asarray(value)

        if array.ndim != 1:
            raise ValueError(
                "ec8_class must be one-dimensional."
            )

        normalized = _np.char.upper(
            array.astype(str)
        )

        valid = (
            (normalized == "A")
            | (normalized == "B")
            | (normalized == "C")
        )

        if not _np.all(valid):
            raise ValueError(
                "ec8_class must contain only 'A', 'B' or 'C'."
            )

        return normalized

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


class MassaEtAl2008Ml(_MassaEtAl2008Base):
    """
    Massa et al. (2008) using local magnitude ML.
    """

    MAGNITUDE_TYPE = "ML"
    _COEFF_SET = "ml"


class MassaEtAl2008Mw(_MassaEtAl2008Base):
    """
    Massa et al. (2008) using moment magnitude Mw.
    """

    MAGNITUDE_TYPE = "Mw"
    _COEFF_SET = "mw"


__all__ = [
    "MassaEtAl2008Ml",
    "MassaEtAl2008Mw",
]
