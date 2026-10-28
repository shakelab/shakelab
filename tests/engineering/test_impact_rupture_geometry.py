"""Finite-fault rupture parameters in the real CY14 -> impact path.

Strike is retained as event metadata but CY14 does not request it.
Finite-fault distances are prescribed, not inferred from strike/dip.
"""

from __future__ import annotations

from math import isfinite

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionAssignment,
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.engineering.impact import ImpactConfig, compute_impact_scenario
from shakelab.groundmotion.providers import ScenarioEvent
from tests.engineering.test_impact_integration import _inputs


SITE = {
    "vs30": 500.0,
    "vs30measured": True,
    "z1pt0": 300.0,
}
DISTANCES = {
    "rrup": 18.0,
    "rjb": 15.0,
    "rx": 8.0,
}
GEOMETRY = {
    "rake": 90.0,
    "dip": 35.0,
    "strike": 40.0,
    "ztor": 2.0,
}


def _runtime(mechanism, *, distances=None):
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=6.0,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism=mechanism,
    )
    assignment = GroundMotionAssignment(
        assets=("asset_1",),
        parameters={
            "site_parameters": SITE,
            "distances": DISTANCES if distances is None else distances,
        },
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="cy14",
                provider="gmpe",
                default=True,
                config={"gmpe_name": "ChiouYoungs2014"},
                assignments=(assignment,),
            ),
        ],
    )
    return GroundMotionRuntime(model, event)


def _evaluate(mechanism, *, distances=None):
    return _runtime(mechanism, distances=distances).evaluate_at_site(
        asset_id="asset_1",
        imt="PGA",
        lon=13.15,
        lat=46.10,
    )


def test_dip_affects_real_gmpe():
    shallow, sigma_shallow = _evaluate({**GEOMETRY, "dip": 30.0})
    steep, sigma_steep = _evaluate({**GEOMETRY, "dip": 75.0})
    assert all(isfinite(v) and v > 0.0 for v in (shallow, steep))
    assert isfinite(sigma_shallow) and sigma_shallow >= 0.0
    assert isfinite(sigma_steep) and sigma_steep >= 0.0
    assert shallow != pytest.approx(steep, rel=1e-8)


def test_strike_does_not_change_cy14_with_fixed_distances():
    first, _ = _evaluate({**GEOMETRY, "strike": 0.0})
    second, _ = _evaluate({**GEOMETRY, "strike": 120.0})
    assert first == pytest.approx(second, rel=1e-12)


@pytest.mark.parametrize("missing", ["dip", "ztor", "rake"])
def test_required_rupture_parameters_cannot_be_omitted(missing):
    mechanism = dict(GEOMETRY)
    del mechanism[missing]
    with pytest.raises((ValueError, AttributeError, TypeError),
                       match=missing):
        _evaluate(mechanism)


def test_finite_fault_distance_is_not_invented():
    distances = {key: val for key, val in DISTANCES.items()
                 if key != "rrup"}
    with pytest.raises((ValueError, AttributeError, TypeError),
                       match="rrup"):
        _evaluate(GEOMETRY, distances=distances)


def test_dip_changes_impact_with_real_gmpe():
    _, exposure, tree, collection = _inputs(0.35)
    exposure.assets[0].reference_location = {
        "longitude": 13.15,
        "latitude": 46.10,
    }

    def calculate(dip):
        runtime = _runtime({**GEOMETRY, "dip": dip})
        result = compute_impact_scenario(
            runtime, exposure, tree, collection,
            config=ImpactConfig(uncertainty_mode="lognormal"),
        )
        return result.assets["asset_1"]

    shallow = calculate(30.0)
    steep = calculate(75.0)
    assert shallow.ground_motion_by_imt["PGA"].median != pytest.approx(
        steep.ground_motion_by_imt["PGA"].median, rel=1e-8,
    )
    assert shallow.probabilities != steep.probabilities
    for asset in (shallow, steep):
        assert sum(asset.probabilities.values()) == pytest.approx(1.0)
        assert sum(asset.expected_counts.values()) == pytest.approx(100.0)
