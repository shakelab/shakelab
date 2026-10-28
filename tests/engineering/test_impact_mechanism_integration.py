"""Verify focal-mechanism propagation through a real GMPE and impact."""

from __future__ import annotations

from math import isfinite

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.engineering.impact import ImpactConfig, compute_impact_scenario
from shakelab.groundmotion.providers import ScenarioEvent

from tests.engineering.test_impact_integration import _inputs


def _runtime(rake):
    mechanism = None if rake is None else {
        "strike": 30.0,
        "dip": 60.0,
        "rake": rake,
    }
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=5.5,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism=mechanism,
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="akkar",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "AkkarEtAlRepi2014"},
            ),
        ],
    )
    return GroundMotionRuntime(model, event)


def _evaluate(rake):
    runtime = _runtime(rake)
    # This variant derives repi from the epicentre and site coordinates.
    # Vs30 is the only required site parameter.
    return runtime.evaluate_at_site(
        asset_id="asset_1",
        imt="PGA",
        lon=13.15,
        lat=46.10,
        site_parameters={"vs30": 500.0},
    )


@pytest.mark.parametrize("rake", [-90.0, 0.0, 90.0])
def test_mechanism_propagates_to_real_gmpe(rake):
    median, sigma = _evaluate(rake)
    assert isfinite(median) and median > 0.0
    assert isfinite(sigma) and sigma >= 0.0


def test_faulting_style_changes_ground_motion():
    normal, _ = _evaluate(-90.0)
    strike_slip, _ = _evaluate(0.0)
    reverse, _ = _evaluate(90.0)
    assert len({round(normal, 12), round(strike_slip, 12),
                round(reverse, 12)}) > 1


def test_missing_required_rake_fails_explicitly():
    with pytest.raises((ValueError, AttributeError, TypeError),
                       match="rake"):
        _evaluate(None)


def test_mechanism_propagates_through_impact():
    _, exposure, tree, collection = _inputs(0.35)
    # Place the exposure asset at the same site used for the GMPE check.
    exposure.assets[0].reference_location = {
        "longitude": 13.15,
        "latitude": 46.10,
    }

    from shakelab.engineering.gmdata.associator import GroundMotionAssignment

    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="akkar",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "AkkarEtAlRepi2014"},
                assignments=(
                    GroundMotionAssignment(
                        assets=("asset_1",),
                        parameters={"site_parameters": {"vs30": 500.0}},
                    ),
                ),
            ),
        ],
    )

    def calculate(rake):
        event = _runtime(rake).event
        runtime = GroundMotionRuntime(model, event)
        result = compute_impact_scenario(
            runtime, exposure, tree, collection,
            config=ImpactConfig(uncertainty_mode="lognormal"),
        )
        return result.assets["asset_1"]

    normal = calculate(-90.0)
    reverse = calculate(90.0)
    assert normal.ground_motion_by_imt["PGA"].median != pytest.approx(
        reverse.ground_motion_by_imt["PGA"].median, rel=1e-8,
    )
    assert normal.probabilities != reverse.probabilities
