"""Tests for optional reference-site and exposure geology properties."""

import json

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    load_ground_motion_model,
)
from tests.engineering.test_impact_reference_vs30 import _asset, _runtime


def _model(site):
    return GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        reference_site=site,
        providers=[
            GroundMotionProviderConfig(
                id="gmpe",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "AkkarEtAlRepi2014"},
            ),
        ],
    )


def test_reference_site_validates_all_optional_properties():
    _model({
        "vs30": 760.0,
        "is_measured": False,
        "z1pt0": 250.0,
        "z2pt5": 1.2,
    }).validate()


@pytest.mark.parametrize("site", [
    {"vs30": 0},
    {"is_measured": 1},
    {"z1pt0": -1},
    {"z2pt5": float("nan")},
    {"unknown": 1},
])
def test_invalid_reference_site_rejected(site):
    with pytest.raises(ValueError):
        _model(site).validate()


def test_reference_site_json_round_trip(tmp_path):
    path = tmp_path / "gm.json"
    path.write_text(json.dumps({
        "type": "ShakeLabGroundMotion",
        "schema_version": "1.0.0",
        "reference_site": {"vs30": 760, "is_measured": False},
        "providers": [{
            "id": "gmpe", "provider": "gmpe", "default": True,
            "config": {"gmpe_name": "AkkarEtAlRepi2014"},
        }],
    }), encoding="utf-8")
    model = load_ground_motion_model(path)
    assert model.reference_site["vs30"] == 760
    assert model.reference_site["is_measured"] is False


def test_asset_optional_geology_validated():
    exposure, _, _ = _asset(350.0)
    exposure.assets[0].reference_geology.update({
        "is_measured": True,
        "z1pt0": 150.0,
        "z2pt5": 0.8,
    })
    exposure.assets[0].validate()


@pytest.mark.parametrize("geology", [
    {"is_measured": "yes"},
    {"z1pt0": -1.0},
    {"z2pt5": float("inf")},
    {"vs30": 0.0},
])
def test_invalid_asset_geology_rejected(geology):
    exposure, _, _ = _asset(350.0)
    exposure.assets[0].reference_geology = geology
    with pytest.raises(ValueError):
        exposure.assets[0].validate()


def test_reference_site_vs30_fallback_in_runtime():
    runtime = _runtime(760.0)
    exposure, _, _ = _asset(None)
    median, _ = runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    assert median > 0


def test_independent_site_parameter_precedence(monkeypatch):
    """Asset values override runtime, assignment, and model defaults."""
    from shakelab.engineering.gmdata.associator import GroundMotionRuntime
    from shakelab.groundmotion.providers import ScenarioEvent

    model = _model({
        "vs30": 760.0,
        "is_measured": False,
        "z1pt0": 300.0,
        "z2pt5": 1.5,
    })
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=5.5,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism={"rake": 90.0},
    )
    runtime = GroundMotionRuntime(model, event)
    exposure, _, _ = _asset(350.0)
    asset = exposure.assets[0]
    asset.reference_geology.update({
        "is_measured": True,
        "z1pt0": 120.0,
    })
    observed = {}

    def capture_site(**kwargs):
        observed.update(kwargs["site_parameters"])
        return 0.1, 0.5

    monkeypatch.setattr(runtime, "evaluate_at_site", capture_site)
    runtime.evaluate_at_asset(
        asset,
        "PGA",
        site_parameters={"vs30": 500.0, "z2pt5": 2.0},
    )
    assert observed == {
        "vs30": 350.0,
        "vs30measured": True,
        "z1pt0": 120.0,
        "z2pt5": 2.0,
    }


def test_reference_site_parameters_passed_to_gmpe(monkeypatch):
    """Missing asset values inherit independent model defaults."""
    from shakelab.engineering.gmdata.associator import GroundMotionRuntime
    from shakelab.groundmotion.providers import ScenarioEvent

    model = _model({
        "vs30": 760.0,
        "is_measured": False,
        "z1pt0": 300.0,
        "z2pt5": 1.5,
    })
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=5.5,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism={"rake": 90.0},
    )
    runtime = GroundMotionRuntime(model, event)
    exposure, _, _ = _asset(None)
    observed = {}

    def capture_site(**kwargs):
        observed.update(kwargs["site_parameters"])
        return 0.1, 0.5

    monkeypatch.setattr(runtime, "evaluate_at_site", capture_site)
    runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    assert observed == {
        "vs30": 760.0,
        "vs30measured": False,
        "z1pt0": 300.0,
        "z2pt5": 1.5,
    }


def test_additional_geology_metadata_is_preserved(monkeypatch):
    """Unrelated geology metadata must not reach the GMPE context."""
    exposure, _, _ = _asset(350.0)
    asset = exposure.assets[0]
    asset.reference_geology.update({
        "geological_unit": "alluvium",
        "description": "Local geological survey",
    })
    asset.validate()

    runtime = _runtime(760.0)
    observed = {}

    def capture_site(**kwargs):
        observed.update(kwargs["site_parameters"])
        return 0.1, 0.5

    monkeypatch.setattr(runtime, "evaluate_at_site", capture_site)
    runtime.evaluate_at_asset(asset, "PGA")

    assert observed["vs30"] == 350.0
    assert "geological_unit" not in observed
    assert "description" not in observed
    assert asset.reference_geology["geological_unit"] == "alluvium"
