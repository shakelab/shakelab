"""Reference Vs30 fallback and asset precedence in the impact pipeline."""

from __future__ import annotations

import json
from math import isfinite

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
    load_ground_motion_model,
)
from shakelab.engineering.impact import ImpactConfig, compute_impact_scenario
from shakelab.groundmotion.providers import ScenarioEvent
from tests.engineering.test_impact_integration import _inputs


def _model(reference_vs30):
    return GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        reference_vs30=reference_vs30,
        providers=[
            GroundMotionProviderConfig(
                id="akkar",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "AkkarEtAlRepi2014"},
            ),
        ],
    )


def _runtime(reference_vs30):
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=5.5,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism={"rake": 90.0},
    )
    return GroundMotionRuntime(_model(reference_vs30), event)


def _asset(vs30):
    _, exposure, tree, collection = _inputs(0.0)
    exposure.assets[0].reference_location = {
        "longitude": 13.15,
        "latitude": 46.10,
    }
    exposure.assets[0].reference_geology = (
        None if vs30 is None else {"vs30": vs30}
    )
    return exposure, tree, collection


def _impact(runtime, exposure, tree, collection):
    result = compute_impact_scenario(
        runtime, exposure, tree, collection,
        config=ImpactConfig(uncertainty_mode="lognormal"),
    )
    return result.assets["asset_1"]


def test_reference_vs30_used_when_asset_has_no_vs30():
    runtime = _runtime(760.0)
    exposure, tree, collection = _asset(None)
    median, sigma = runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    assert median > 0.0 and isfinite(sigma)
    impact = _impact(runtime, exposure, tree, collection)
    assert impact.ground_motion_by_imt["PGA"].median == pytest.approx(median)


def test_asset_vs30_overrides_global_reference():
    runtime = _runtime(760.0)
    exposure, tree, collection = _asset(350.0)
    asset_median, _ = runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    impact = _impact(runtime, exposure, tree, collection)
    assert impact.ground_motion_by_imt["PGA"].median == pytest.approx(
        asset_median,
    )
    exposure.assets[0].reference_geology = None
    reference_median, _ = runtime.evaluate_at_asset(
        exposure.assets[0], "PGA",
    )
    assert asset_median != pytest.approx(reference_median, rel=1e-8)


def test_missing_vs30_is_not_silently_defaulted():
    runtime = _runtime(None)
    exposure, tree, collection = _asset(None)
    with pytest.raises((ValueError, AttributeError), match="vs30"):
        _impact(runtime, exposure, tree, collection)


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_reference_vs30_rejected(value):
    with pytest.raises(ValueError, match="reference_vs30"):
        _model(value).validate()


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan")])
def test_invalid_asset_vs30_rejected(value):
    runtime = _runtime(760.0)
    exposure, _, _ = _asset(value)
    with pytest.raises(ValueError, match="vs30"):
        runtime.evaluate_at_asset(exposure.assets[0], "PGA")


def test_reference_vs30_loaded_from_json(tmp_path):
    path = tmp_path / "ground_motion.json"
    path.write_text(json.dumps({
        "type": "ShakeLabGroundMotion",
        "schema_version": "1.0.0",
        "reference_vs30": 800.0,
        "providers": [{
            "id": "akkar",
            "provider": "gmpe",
            "default": True,
            "config": {"gmpe_name": "AkkarEtAlRepi2014"},
        }],
    }), encoding="utf-8")
    model = load_ground_motion_model(path)
    assert model.reference_vs30 == pytest.approx(800.0)
