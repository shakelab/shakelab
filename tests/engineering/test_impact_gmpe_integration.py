"""Integration test: real Kotha GMPE through ShakeScenario impact runtime.

Finite-fault distance Rjb is provided explicitly; it is not inferred
from a point-source hypocentre. Fragility parameters are controlled.
"""

from __future__ import annotations

from math import isfinite, log, sqrt

import pytest
from scipy.stats import norm

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.engineering.impact import ImpactConfig, compute_impact_scenario
from shakelab.groundmotion.providers import ScenarioEvent

from tests.engineering.test_impact_integration import (
    BETA,
    COUNT,
    THETA,
    _inputs,
)


def _runtime():
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=5.5,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        event_id="integration-gmpe",
        mechanism={"strike": 45.0, "dip": 30.0, "rake": 90.0},
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="kotha",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "KothaEtAl2020ESHM20"},
            ),
        ],
    )
    runtime = GroundMotionRuntime(model, event)
    return runtime, event


def test_real_gmpe_to_impact_probabilities():
    _, exposure, tree, collection = _inputs(0.35)
    runtime, event = _runtime()
    assert runtime.event.mechanism["rake"] == 90.0

    # Rjb is a finite-fault metric: specify it, never approximate it.
    parameters = {
        "site_parameters": {
            "vs30": 500.0,
            "vs30measured": True,
            "eshm20_region": 0,
        },
        "distances": {"rjb": 15.0},
    }
    # Assign the site and rupture geometry to the configured provider.
    runtime.model.providers[0].assignments  # no mutation of config
    from shakelab.engineering.gmdata.associator import GroundMotionAssignment

    provider_cfg = GroundMotionProviderConfig(
        id="kotha",
        provider="gmpe",
        default=True,
        config={"gmpe_name": "KothaEtAl2020ESHM20"},
        assignments=(
            GroundMotionAssignment(
                assets=("asset_1",),
                parameters=parameters,
            ),
        ),
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[provider_cfg],
    )
    runtime = GroundMotionRuntime(model, event)
    result = compute_impact_scenario(
        runtime, exposure, tree, collection,
        config=ImpactConfig(uncertainty_mode="lognormal"),
    )
    asset = result.assets["asset_1"]
    gm = asset.ground_motion_by_imt["PGA"]

    assert gm.provider_id == "kotha"
    assert isfinite(gm.median) and gm.median > 0.0
    assert isfinite(gm.sigma_ln) and gm.sigma_ln >= 0.0

    # Independent analytical convolution of the actual GMPE output.
    spread = sqrt(BETA ** 2 + gm.sigma_ln ** 2)
    exceed = [
        norm.cdf(log(gm.median / theta) / spread)
        for theta in THETA.values()
    ]
    expected = {
        "D0": 1.0 - exceed[0],
        "D1": exceed[0] - exceed[1],
        "D2": exceed[1] - exceed[2],
        "D3": exceed[2],
    }
    assert sum(asset.probabilities.values()) == pytest.approx(1.0)
    for state, probability in expected.items():
        assert asset.probabilities[state] == pytest.approx(
            probability, abs=1e-11,
        )
        assert asset.expected_counts[state] == pytest.approx(
            COUNT * probability, abs=1e-9,
        )
    assert sum(asset.expected_counts.values()) == pytest.approx(COUNT)
