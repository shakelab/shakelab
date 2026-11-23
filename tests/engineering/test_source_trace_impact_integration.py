"""Regression of opt-in source traces through the real impact pipeline.

These tests require no ShakeScenario server, job database or WebUI.
They verify scientific equivalence and the exported trace data contract.
"""

from __future__ import annotations

import json

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
)
from shakelab.engineering.impact import (
    ImpactConfig,
    compute_impact_scenario,
    save_impact_result,
    save_impact_summary,
)
from shakelab.groundmotion.providers import ScenarioEvent
from tests.engineering.test_impact_integration import _inputs


SOURCE = {
    "strike": 120.0,
    "dip": 35.0,
    "rake": 90.0,
    "aspect_ratio": 2.0,
    "lower_seismogenic_depth": 15.0,
}


def _runtime(record_source_trace):
    event = ScenarioEvent.from_lonlat_depth(
        6.2, 13.0, 46.0, 7.0,
        magnitude_type="ML",
        source_parameters={"strike": 90.0},
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="scaled_gmpe",
                provider="gmpe",
                default=True,
                config={
                    "gmpe_name": "AbrahamsonEtAl2014",
                    "source_geometry": "finite",
                    "source_parameters": SOURCE,
                    "source_scaling": {
                        "model": "ThingbaijamEtAl2017",
                        "style": "auto",
                    },
                    "site_parameters": {
                        "vs30": 500.0,
                        "vs30measured": True,
                        "z1pt0": 300.0,
                    },
                },
            ),
        ],
    )
    model.validate()
    return model.runtime(
        event, record_source_trace=record_source_trace,
    )


def _impact(runtime):
    _, exposure, taxonomy, fragilities = _inputs(0.35)
    return compute_impact_scenario(
        ground_motion=runtime,
        exposure_model=exposure,
        taxonomy_tree=taxonomy,
        fragility_collection=fragilities,
        config=ImpactConfig(uncertainty_mode="lognormal"),
    )


def test_trace_does_not_change_ground_motion_or_damage():
    without = _runtime(False)
    with_trace = _runtime(True)
    baseline = _impact(without)
    traced = _impact(with_trace)

    first = baseline.assets["asset_1"]
    second = traced.assets["asset_1"]
    first_gm = first.ground_motion_by_imt["PGA"]
    second_gm = second.ground_motion_by_imt["PGA"]

    assert second_gm.median == pytest.approx(first_gm.median)
    assert second_gm.sigma_ln == pytest.approx(first_gm.sigma_ln)
    assert second.probabilities == pytest.approx(first.probabilities)
    assert second.expected_counts == pytest.approx(first.expected_counts)
    assert without.source_trace()["providers"] == {}

    trace = with_trace.source_trace()
    assert list(trace["providers"]) == ["scaled_gmpe"]
    evaluation, = trace["providers"]["scaled_gmpe"]["evaluations"]
    assert evaluation["imt"] == "PGA"
    assert evaluation["source_geometry"] == "finite"
    assert evaluation["magnitude_type"] == "ML"
    assert evaluation["source_parameters"]["strike"] == 90.0
    assert evaluation["input_layers"]["strike"] == "event"
    assert evaluation["builder_origins"]["area"] == "scaling"
    assert evaluation["scaling"]["used_model"] == "ThingbaijamEtAl2017"
    assert evaluation["rupture"]["area_km2"] > 0
    assert len(evaluation["rupture"]["vertices"]) == 4
    assert evaluation["distances"]["rrup"]["origin"] == "finite_rupture"
    assert evaluation["assumptions"][
        "magnitude_conversion_applied"
    ] is False
    json.dumps(trace, allow_nan=False)


def test_trace_and_impact_artifacts_are_independently_serializable(tmp_path):
    runtime = _runtime(True)
    impact = _impact(runtime)
    config = ImpactConfig(uncertainty_mode="lognormal")
    results = tmp_path / "results"
    results.mkdir()
    impact_path = results / "impact_assets.json"
    summary_path = results / "impact_summary.json"
    trace_path = results / "source_trace.json"

    save_impact_result(impact, str(impact_path), config)
    save_impact_summary(impact, str(summary_path), config)
    trace_path.write_text(
        json.dumps(
            runtime.source_trace(),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        ),
        encoding="utf-8",
    )

    assets = json.loads(impact_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    trace = json.loads(trace_path.read_text(encoding="utf-8"))

    assert assets["type"] == "ShakeScenarioImpactAssets"
    assert "source_trace" not in assets
    assert "source_trace" not in summary
    assert trace["type"] == "ShakeScenarioSourceTrace"
    assert len(trace["providers"]["scaled_gmpe"]["evaluations"]) == 1


def test_trace_is_separate_per_runtime():
    first = _runtime(True)
    second = _runtime(True)
    _impact(first)
    assert len(first.source_trace()["providers"][
        "scaled_gmpe"
    ]["evaluations"]) == 1
    assert second.source_trace()["providers"] == {}
    _impact(second)
    assert len(second.source_trace()["providers"][
        "scaled_gmpe"
    ]["evaluations"]) == 1
