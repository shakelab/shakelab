"""End-to-end GMPE configuration checks through the impact pipeline."""

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


RUPTURE = {"rake": 90.0, "dip": 35.0, "ztor": 2.0}
DISTANCES = {"rrup": 18.0, "rjb": 15.0, "rx": 8.0}


def _runtime(
    *,
    reference_site=None,
    mechanism=None,
    distances=None,
):
    event = ScenarioEvent.from_lonlat_depth(
        magnitude=6.0,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        mechanism=RUPTURE if mechanism is None else mechanism,
    )
    assignment = GroundMotionAssignment(
        assets=("asset_1",),
        parameters={
            "distances": DISTANCES if distances is None else distances,
        },
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        reference_site=reference_site or {},
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
    model.validate()
    return GroundMotionRuntime(model, event)


def _asset():
    _, exposure, tree, collection = _inputs(0.35)
    exposure.assets[0].reference_location = {
        "longitude": 13.15,
        "latitude": 46.10,
    }
    exposure.assets[0].reference_geology = None
    return exposure, tree, collection


def _impact(runtime, exposure, tree, collection):
    return compute_impact_scenario(
        runtime,
        exposure,
        tree,
        collection,
        config=ImpactConfig(uncertainty_mode="lognormal"),
    ).assets["asset_1"]


def test_reference_site_reaches_real_gmpe_and_impact():
    runtime = _runtime(reference_site={
        "vs30": 500.0,
        "is_measured": True,
        "z1pt0": 300.0,
    })
    exposure, tree, collection = _asset()
    median, sigma = runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    impact = _impact(runtime, exposure, tree, collection)

    assert isfinite(median) and median > 0.0
    assert isfinite(sigma) and sigma >= 0.0
    assert impact.ground_motion_by_imt["PGA"].median == pytest.approx(median)
    assert sum(impact.probabilities.values()) == pytest.approx(1.0)


def test_asset_vs30_overrides_reference_in_real_gmpe():
    runtime = _runtime(reference_site={
        "vs30": 760.0,
        "is_measured": True,
        "z1pt0": 300.0,
    })
    exposure, _, _ = _asset()
    baseline, _ = runtime.evaluate_at_asset(exposure.assets[0], "PGA")
    exposure.assets[0].reference_geology = {"vs30": 350.0}
    modified, _ = runtime.evaluate_at_asset(exposure.assets[0], "PGA")

    assert modified != pytest.approx(baseline, rel=1e-8)


@pytest.mark.parametrize(
    ("reference_site", "missing"),
    [
        ({"vs30": 500.0, "z1pt0": 300.0}, "vs30measured"),
        ({"vs30": 500.0, "is_measured": True}, "z1pt0"),
        ({"is_measured": True, "z1pt0": 300.0}, "vs30"),
    ],
)
def test_missing_site_requirement_fails_in_impact(reference_site, missing):
    runtime = _runtime(reference_site=reference_site)
    exposure, tree, collection = _asset()

    with pytest.raises((ValueError, AttributeError), match=missing):
        _impact(runtime, exposure, tree, collection)


def test_missing_rupture_requirement_fails_in_impact():
    runtime = _runtime(
        reference_site={
            "vs30": 500.0,
            "is_measured": True,
            "z1pt0": 300.0,
        },
        mechanism={"rake": 90.0, "ztor": 2.0},
    )
    exposure, tree, collection = _asset()

    with pytest.raises((ValueError, AttributeError), match="dip"):
        _impact(runtime, exposure, tree, collection)


def test_missing_finite_fault_distance_fails_in_impact():
    runtime = _runtime(
        reference_site={
            "vs30": 500.0,
            "is_measured": True,
            "z1pt0": 300.0,
        },
        distances={"rjb": 15.0, "rx": 8.0},
    )
    exposure, tree, collection = _asset()

    with pytest.raises((ValueError, AttributeError), match="rrup"):
        _impact(runtime, exposure, tree, collection)
