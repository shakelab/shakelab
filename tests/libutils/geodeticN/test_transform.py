#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Unit tests for :mod:`shakelab.libutils.geodeticN.transform`.

These tests validate the core WGS84 <-> Cartesian transforms implemented
in ``transform.py`` using the refactored API based on ``MetricFrame``.

Tested functionality:

- geo_to_ecef / ecef_to_geo (absolute ECEF)
- geo_to_enu / enu_to_geo (local ENU frame)
- wgs_to_metric / metric_to_wgs (high-level helpers)

Run with, e.g.:

    python -m unittest test_transform
or:

    python test_transform.py
"""

from __future__ import annotations

import unittest
from math import isclose

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    MetricPoint,
    geo_to_ecef,
    ecef_to_geo,
    geo_to_enu,
    enu_to_geo,
    wgs_to_metric,
    metric_to_wgs,
)


class TestTransform(unittest.TestCase):
    """Validation tests for geodetic <-> Cartesian transforms."""

    # ------------------------------------------------------------------
    # Absolute ECEF round-trip
    # ------------------------------------------------------------------

    def test_geo_ecef_roundtrip(self) -> None:
        """geo_to_ecef and ecef_to_geo must invert each other."""
        cases = [
            (0.0, 0.0, 0.0),
            (12.0, 45.0, 1000.0),
            (-123.5, 60.0, -500.0),
        ]

        for lon, lat, h in cases:
            x, y, z = geo_to_ecef((lon, lat, h))
            lon2, lat2, h2 = ecef_to_geo((x, y, z))

            self.assertTrue(
                isclose(lon, lon2, abs_tol=1e-8)
                and isclose(lat, lat2, abs_tol=1e-8)
                and isclose(h, h2, abs_tol=1e-3),
                msg=f"ECEF round-trip failed for {(lon, lat, h)}",
            )

    # ------------------------------------------------------------------
    # High-level helpers (absolute ECEF)
    # ------------------------------------------------------------------

    def test_wgs_metric_roundtrip_absolute(self) -> None:
        """wgs_to_metric / metric_to_wgs without frame must round-trip."""
        cases = [
            WgsPoint(0.0, 0.0, 0.0),
            WgsPoint(13.0, 46.0, 200.0),
            WgsPoint(-75.0, -30.0, 1500.0),
        ]

        for p in cases:
            mp = wgs_to_metric(p)
            self.assertIsNone(mp.frame)

            p2 = metric_to_wgs(mp)
            self.assertTrue(
                isclose(p.longitude, p2.longitude, abs_tol=1e-8)
                and isclose(p.latitude, p2.latitude, abs_tol=1e-8)
                and isclose(p.elevation, p2.elevation, abs_tol=1e-3),
                msg=f"Absolute metric round-trip failed for {p}",
            )

    # ------------------------------------------------------------------
    # Local ECEF frame
    # ------------------------------------------------------------------

    def test_local_ecef_reference(self) -> None:
        """Local ECEF frame must map reference point to (0, 0, 0)."""
        ref = WgsPoint(13.0, 46.0, 200.0)
        frame = MetricFrame(ref_geo=ref, orientation="ecef")

        mp0 = wgs_to_metric(ref, frame=frame)
        self.assertTrue(
            isclose(mp0.x, 0.0, abs_tol=1e-6)
            and isclose(mp0.y, 0.0, abs_tol=1e-6)
            and isclose(mp0.z, 0.0, abs_tol=1e-6),
            msg="Reference point not mapped to origin in local ECEF",
        )

        p = WgsPoint(13.01, 46.02, 250.0)
        mp = wgs_to_metric(p, frame=frame)
        p2 = metric_to_wgs(mp)

        self.assertTrue(
            isclose(p.longitude, p2.longitude, abs_tol=1e-8)
            and isclose(p.latitude, p2.latitude, abs_tol=1e-8)
            and isclose(p.elevation, p2.elevation, abs_tol=1e-3),
            msg="Local ECEF round-trip failed",
        )

    # ------------------------------------------------------------------
    # ENU frame (low-level)
    # ------------------------------------------------------------------

    def test_geo_enu_roundtrip(self) -> None:
        """geo_to_enu and enu_to_geo must invert each other."""
        ref = WgsPoint(13.0, 46.0, 200.0)
        frame = MetricFrame(ref_geo=ref, orientation="enu")

        cases = [
            ref,
            WgsPoint(13.01, 46.02, 250.0),
            WgsPoint(12.98, 45.99, 150.0),
        ]

        for p in cases:
            e, n, u = geo_to_enu(p, frame=frame)
            lon2, lat2, h2 = enu_to_geo((e, n, u), frame=frame)

            self.assertTrue(
                isclose(p.longitude, lon2, abs_tol=1e-8)
                and isclose(p.latitude, lat2, abs_tol=1e-8)
                and isclose(p.elevation, h2, abs_tol=1e-3),
                msg=f"ENU round-trip failed for {p}",
            )

    # ------------------------------------------------------------------
    # ENU frame (high-level helpers)
    # ------------------------------------------------------------------

    def test_wgs_metric_roundtrip_enu(self) -> None:
        """wgs_to_metric / metric_to_wgs must round-trip in ENU."""
        ref = WgsPoint(13.0, 46.0, 200.0)
        frame = MetricFrame(ref_geo=ref, orientation="enu")

        cases = [
            ref,
            WgsPoint(13.02, 45.98, 150.0),
            WgsPoint(12.99, 46.03, 300.0),
        ]

        for p in cases:
            mp = wgs_to_metric(p, frame=frame)
            self.assertIsInstance(mp, MetricPoint)
            self.assertEqual(mp.frame, frame)

            p2 = metric_to_wgs(mp)
            self.assertTrue(
                isclose(p.longitude, p2.longitude, abs_tol=1e-8)
                and isclose(p.latitude, p2.latitude, abs_tol=1e-8)
                and isclose(p.elevation, p2.elevation, abs_tol=1e-3),
                msg=f"High-level ENU round-trip failed for {p}",
            )

    def test_wgs_metric_enu_consistency(self) -> None:
        """wgs_to_metric(frame=ENU) must match geo_to_enu."""
        ref = WgsPoint(13.0, 46.0, 200.0)
        frame = MetricFrame(ref_geo=ref, orientation="enu")
        p = WgsPoint(13.015, 46.005, 220.0)

        e, n, u = geo_to_enu(p, frame=frame)
        mp = wgs_to_metric(p, frame=frame)

        self.assertTrue(
            isclose(mp.x, e, abs_tol=1e-6)
            and isclose(mp.y, n, abs_tol=1e-6)
            and isclose(mp.z, u, abs_tol=1e-6),
            msg="ENU inconsistency between low- and high-level helpers",
        )


if __name__ == "__main__":
    unittest.main()
