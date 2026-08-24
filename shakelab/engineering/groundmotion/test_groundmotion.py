# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
# ****************************************************************************
"""
Tests for the engineering ground-motion configuration layer.

Run from the ShakeLab repository root with:

    python -m unittest test_groundmotion_v0.py -v

or:

    python test_groundmotion_v0.py

The tests use a small provider registered only for the test process. No
changes to shakelab.gmmodel are required.
"""

from __future__ import annotations

import json
import tempfile
import unittest

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from shakelab.engineering.groundmotion.groundmotion import (
    GroundMotionModel,
)
from shakelab.gmmodel.groundmotion import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint


_TEST_PROVIDER_ID = "_groundmotion_test_provider_v0"


class _DummyProvider:
    """Minimal deterministic provider used only by this test module."""

    instances = 0
    calls = 0

    def __init__(
        self,
        base_value: float = 1.0,
        **kwargs: Any,
    ) -> None:
        _DummyProvider.instances += 1
        self.base_value = float(base_value)
        self.kwargs = dict(kwargs)

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> tuple[list[float], list[float]]:
        _DummyProvider.calls += 1

        station_code = kwargs.get("station_code")
        offset = 0.0

        if station_code == "UD01":
            offset = 10.0
        elif station_code == "UD02":
            offset = 20.0

        values = [
            self.base_value + offset
            for _site in sites
        ]
        sigma = [0.2 for _site in sites]

        return values, sigma


if _TEST_PROVIDER_ID not in GroundMotionProvider.available_ids():
    GroundMotionProvider.register(_TEST_PROVIDER_ID)(_DummyProvider)


@dataclass
class _Asset:
    id: str
    reference_location: dict[str, float]


@dataclass
class _Exposure:
    assets: list[_Asset]


def _make_exposure() -> _Exposure:
    return _Exposure(
        assets=[
            _Asset(
                id="A001",
                reference_location={
                    "longitude": 13.0,
                    "latitude": 46.0,
                    "elevation": 100.0,
                },
            ),
            _Asset(
                id="A002",
                reference_location={
                    "longitude": 13.1,
                    "latitude": 46.1,
                    "elevation": 120.0,
                },
            ),
            _Asset(
                id="A003",
                reference_location={
                    "longitude": 13.2,
                    "latitude": 46.2,
                    "elevation": 140.0,
                },
            ),
        ]
    )


def _event() -> ScenarioEvent:
    return ScenarioEvent(
        hypocentre=WgsPoint(
            longitude=13.0,
            latitude=46.0,
            elevation=-10000.0,
        ),
        magnitude=5.5,
    )


def _base_config() -> dict[str, Any]:
    return {
        "type": "ShakeLabGroundMotion",
        "schema_version": "1.0.0",
        "providers": [
            {
                "id": "default_test",
                "provider": _TEST_PROVIDER_ID,
                "default": True,
                "config": {
                    "base_value": 1.0,
                },
            },
            {
                "id": "observed_test",
                "provider": _TEST_PROVIDER_ID,
                "config": {
                    "base_value": 2.0,
                },
                "assignments": [
                    {
                        "assets": [
                            "A002",
                        ],
                        "parameters": {
                            "station_code": "UD01",
                        },
                    },
                    {
                        "assets": [
                            "A003",
                        ],
                        "parameters": {
                            "station_code": "UD02",
                        },
                    },
                ],
            },
        ],
    }


def _write_json(
    directory: Path,
    payload: dict[str, Any],
    name: str = "ground_motion.json",
) -> Path:
    path = directory / name
    path.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    return path


class GroundMotionModelTest(unittest.TestCase):
    def setUp(self) -> None:
        _DummyProvider.instances = 0
        _DummyProvider.calls = 0
        self.exposure = _make_exposure()

    def test_load_and_resolve_default_and_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), _base_config())

            model = GroundMotionModel.from_json(
                path,
                exposure_model=self.exposure,
                validate=True,
            )

        cfg, params, explicit = model.resolve_config("A001")
        self.assertEqual(cfg.id, "default_test")
        self.assertEqual(params, {})
        self.assertFalse(explicit)

        cfg, params, explicit = model.resolve_config("A002")
        self.assertEqual(cfg.id, "observed_test")
        self.assertEqual(params["station_code"], "UD01")
        self.assertTrue(explicit)

        cfg, params, explicit = model.resolve_config("A003")
        self.assertEqual(cfg.id, "observed_test")
        self.assertEqual(params["station_code"], "UD02")
        self.assertTrue(explicit)

    def test_runtime_instantiates_each_configured_provider_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), _base_config())

            model = GroundMotionModel.from_json(
                path,
                exposure_model=self.exposure,
                validate=True,
            )

            runtime = model.runtime(_event())

        self.assertEqual(_DummyProvider.instances, 2)

        self.assertIs(
            runtime.provider("default_test"),
            runtime.provider("default_test"),
        )
        self.assertIs(
            runtime.provider("observed_test"),
            runtime.provider("observed_test"),
        )
        self.assertIs(
            runtime.context("observed_test"),
            runtime.context("observed_test"),
        )

    def test_assignment_parameters_reach_provider(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), _base_config())

            model = GroundMotionModel.from_json(
                path,
                exposure_model=self.exposure,
                validate=True,
            )
            runtime = model.runtime(_event())

        a001 = self.exposure.assets[0]
        a002 = self.exposure.assets[1]
        a003 = self.exposure.assets[2]

        median, sigma = runtime.evaluate_at_asset(a001, "PGA")
        self.assertEqual(median, 1.0)
        self.assertEqual(sigma, 0.2)

        median, sigma = runtime.evaluate_at_asset(a002, "PGA")
        self.assertEqual(median, 12.0)
        self.assertEqual(sigma, 0.2)

        median, sigma = runtime.evaluate_at_asset(a003, "PGA")
        self.assertEqual(median, 22.0)
        self.assertEqual(sigma, 0.2)

    def test_runtime_kwargs_override_assignment_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), _base_config())

            model = GroundMotionModel.from_json(
                path,
                exposure_model=self.exposure,
                validate=True,
            )
            runtime = model.runtime(_event())

        a002 = self.exposure.assets[1]

        median, _sigma = runtime.evaluate_at_asset(
            a002,
            "PGA",
            station_code="UD02",
        )

        self.assertEqual(median, 22.0)

    def test_duplicate_asset_assignment_is_rejected(self) -> None:
        payload = _base_config()
        payload["providers"][0]["assignments"] = [
            {
                "assets": ["A002"],
            }
        ]

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            with self.assertRaisesRegex(
                ValueError,
                "assigned more than once",
            ):
                GroundMotionModel.from_json(
                    path,
                    exposure_model=self.exposure,
                    validate=True,
                )

    def test_multiple_defaults_are_rejected(self) -> None:
        payload = _base_config()
        payload["providers"][1]["default"] = True

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            with self.assertRaisesRegex(
                ValueError,
                "Only one ground-motion provider",
            ):
                GroundMotionModel.from_json(
                    path,
                    exposure_model=self.exposure,
                    validate=True,
                )

    def test_unknown_assignment_asset_is_rejected(self) -> None:
        payload = _base_config()
        payload["providers"][1]["assignments"].append(
            {
                "assets": ["UNKNOWN"],
            }
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            with self.assertRaisesRegex(
                ValueError,
                "unknown exposure assets",
            ):
                GroundMotionModel.from_json(
                    path,
                    exposure_model=self.exposure,
                    validate=True,
                )

    def test_missing_coverage_without_default_is_rejected(self) -> None:
        payload = _base_config()
        payload["providers"][0].pop("default")

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            with self.assertRaisesRegex(
                ValueError,
                "without a ground-motion assignment",
            ):
                GroundMotionModel.from_json(
                    path,
                    exposure_model=self.exposure,
                    validate=True,
                )

    def test_unknown_backend_is_rejected(self) -> None:
        payload = _base_config()
        payload["providers"][0]["provider"] = "not_registered"

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            with self.assertRaisesRegex(
                ValueError,
                "Unknown ground-motion provider backend",
            ):
                GroundMotionModel.from_json(
                    path,
                    exposure_model=self.exposure,
                    validate=True,
                )

    def test_resolution_without_default_fails_at_runtime(self) -> None:
        payload = {
            "type": "ShakeLabGroundMotion",
            "schema_version": "1.0.0",
            "providers": [
                {
                    "id": "explicit_only",
                    "provider": _TEST_PROVIDER_ID,
                    "assignments": [
                        {
                            "assets": ["A001"],
                        }
                    ],
                }
            ],
        }

        with tempfile.TemporaryDirectory() as tmp:
            path = _write_json(Path(tmp), payload)

            model = GroundMotionModel.from_json(
                path,
                validate=True,
            )

        with self.assertRaisesRegex(
            KeyError,
            "no default provider",
        ):
            model.resolve_config("A002")


if __name__ == "__main__":
    unittest.main(verbosity=2)
