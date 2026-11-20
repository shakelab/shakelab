"""Regression checks for provider-level GMPE parameters."""

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionAssignment,
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint
from tests.engineering.test_gmpe_configuration_integration import _asset


def _event(source=None):
    return ScenarioEvent.from_lonlat_depth(
        magnitude=6.0, longitude=13.0, latitude=46.0, depth_km=10.0,
        source_parameters=source,
    )


def _provider(source=None, site=None):
    return GroundMotionProvider.gmpe(
        "ChiouYoungs2014",
        config={
            "point_source_distances": True,
            "source_parameters": source,
            "site_parameters": site,
        },
    )


SOURCE = {"strike": 45.0, "rake": 90.0, "dip": 35.0, "ztor": 2.0}
SITE = {"vs30": 500.0, "is_measured": True, "z1pt0": 300.0}


def test_provider_defaults_match_explicit_event_and_site():
    sites = [WgsPoint(13.15, 46.1)]
    configured = _provider(SOURCE, SITE).evaluate("PGA", sites, _event())
    explicit = _provider().evaluate(
        "PGA", sites, _event(SOURCE), site_parameters=SITE,
    )
    for actual, expected in zip(configured, explicit):
        assert actual == pytest.approx(expected)


def test_event_source_overrides_provider_in_context_and_distances():
    sites = [WgsPoint(13.15, 46.1)]
    event = _event({"strike": 120.0, "rake": -90.0})
    configured = _provider(SOURCE, SITE).evaluate("PGA", sites, event)
    explicit = _provider().evaluate(
        "PGA", sites, _event({**SOURCE, **event.source_parameters}),
        site_parameters=SITE,
    )
    for actual, expected in zip(configured, explicit):
        assert actual == pytest.approx(expected)


@pytest.mark.parametrize("name", ["mag", "hypo_depth"])
def test_reserved_source_overrides_rejected(name):
    with pytest.raises(ValueError, match="reserved"):
        _provider({name: 5.0})
    with pytest.raises(ValueError, match="reserved"):
        _provider(SOURCE, SITE).evaluate(
            "PGA", [WgsPoint(13.15, 46.1)], _event(),
            rupture_parameters={name: 5.0},
        )


def test_missing_source_parameter_still_fails():
    source = {key: value for key, value in SOURCE.items() if key != "dip"}
    with pytest.raises((ValueError, AttributeError), match="dip"):
        _provider(source, SITE).evaluate(
            "PGA", [WgsPoint(13.15, 46.1)], _event(),
        )


def test_site_precedence_and_provider_isolation():
    model = GroundMotionModel(
        type="ShakeLabGroundMotion", schema_version="1.0.0",
        reference_site={"vs30": 760.0, "z1pt0": 300.0},
        providers=[
            GroundMotionProviderConfig(
                id="a", provider="gmpe", default=True,
                config={
                    "gmpe_name": "ChiouYoungs2014",
                    "source_parameters": SOURCE,
                    "site_parameters": {"vs30": 600, "is_measured": False},
                    "point_source_distances": True,
                },
                assignments=(GroundMotionAssignment(
                    assets=("asset_1",),
                    parameters={"site_parameters": {"vs30": 500}},
                ),),
            ),
            GroundMotionProviderConfig(
                id="b", provider="gmpe",
                config={
                    "gmpe_name": "ChiouYoungs2014",
                    "source_parameters": SOURCE,
                    "site_parameters": {"vs30": 800, "is_measured": True},
                    "point_source_distances": True,
                },
                assignments=(GroundMotionAssignment(
                    assets=("other",), parameters={},
                ),),
            ),
        ],
    )
    model.validate()
    runtime = GroundMotionRuntime(model, _event())
    exposure, _, _ = _asset()
    asset = exposure.assets[0]
    asset.reference_geology = {"vs30": 350.0}
    loc = asset.reference_location

    def evaluate(asset_id, **site):
        return runtime.evaluate_at_site(
            asset_id, "PGA", loc["longitude"], loc["latitude"],
            site_parameters=site,
        )

    assert runtime.evaluate_at_asset(
        asset, "PGA", site_parameters={"vs30": 400},
    ) == pytest.approx(evaluate("asset_1", vs30=350))
    assert evaluate("asset_1") == pytest.approx(evaluate("asset_1", vs30=500))
    assert evaluate("other") == pytest.approx(
        evaluate("other", vs30=800, vs30measured=True),
    )
    assert evaluate("other")[0] != pytest.approx(evaluate("asset_1")[0])


def test_invalid_provider_config_rejected_during_model_validation():
    config = GroundMotionProviderConfig(
        id="invalid", provider="gmpe",
        config={
            "gmpe_name": "ChiouYoungs2014",
            "site_parameters": {"vs30": -1.0},
        },
    )
    with pytest.raises(ValueError, match="vs30"):
        config.validate()
