"""Traceability of resolved sources and GMPE distance inputs."""

import json

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.groundmotion.providers import GroundMotionProvider
from tests.engineering.test_finite_source_scaling_provider import (
    SITE,
    SOURCE,
    SITES,
    event,
    provider,
)


def test_trace_contains_final_geometry_and_scaling_provenance():
    gmpe = provider()
    ev = event()
    gmpe.evaluate("PGA", SITES, ev, site_parameters=SITE)
    trace = gmpe.source_trace()
    assert len(trace) == 1
    record = trace[0]
    source = gmpe._resolve_source(ev, {})
    assert record["rupture"]["width_km"] == pytest.approx(
        source.rupture.width
    )
    assert record["rupture"]["ztor_km"] == pytest.approx(
        source.rupture.ztor
    )
    assert len(record["rupture"]["vertices"]) == 4
    assert record["builder_origins"]["area"] == "scaling"
    assert record["scaling"]["used_model"] == "ThingbaijamEtAl2017"
    assert record["scaling"]["estimate"]["median_area_km2"] > 0
    assert record["distances"]["rrup"]["origin"] == "finite_rupture"
    assert record["assumptions"]["magnitude_conversion_applied"] is False
    json.dumps(record, allow_nan=False)


def test_distance_override_is_recorded_without_changing_geometry():
    gmpe = provider()
    gmpe.evaluate(
        "PGA", SITES, event(), site_parameters=SITE,
        distances={"rrup": 42.0},
    )
    record = gmpe.source_trace()[0]
    assert record["distances"]["rrup"]["origin"] == "explicit"
    assert record["distances"]["rrup"]["values_km"] == [42.0] * len(SITES)
    assert record["distances"]["rjb"]["origin"] == "finite_rupture"


def test_layer_provenance_retains_runtime_override():
    gmpe = provider(source={**SOURCE, "strike": 45.0})
    ev = event(source={"strike": 90.0})
    gmpe.evaluate(
        "PGA", SITES, ev, site_parameters=SITE,
        rupture_parameters={"strike": 120.0},
    )
    trace = gmpe.source_trace()[0]
    assert trace["source_parameters"]["strike"] == 120.0
    assert trace["input_layers"]["strike"] == "runtime"
    assert trace["input_layers"]["dip"] == "provider"


def test_trace_records_are_detached():
    gmpe = provider()
    gmpe.evaluate("PGA", SITES, event(), site_parameters=SITE)
    first = gmpe.source_trace()
    first[0]["rupture"]["width_km"] = -1.0
    assert gmpe.source_trace()[0]["rupture"]["width_km"] > 0


def test_explicit_rupture_does_not_claim_empirical_estimate():
    gmpe = provider(source={
        "strike": 120.0, "dip": 35.0, "rake": 90.0,
        "length": 20.0, "width": 12.0,
        "ztor": 2.0, "along_fraction": 0.5,
    })
    gmpe.evaluate("PGA", SITES, event(), site_parameters=SITE)
    scaling = gmpe.source_trace()[0]["scaling"]
    assert scaling["used_model"] is None
    assert scaling["estimate"] is None


def test_nonfinite_geometry_has_no_finite_plane():
    gmpe = GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "point",
            "source_parameters": {
                "strike": 120.0, "dip": 35.0, "rake": 90.0,
                "ztor": 2.0,
            },
        },
    )
    gmpe.evaluate("PGA", SITES, event(), site_parameters=SITE)
    record = gmpe.source_trace()[0]
    assert record["rupture"] is None
    assert record["source_geometry"] == "point"


def test_runtime_source_trace_is_provider_scoped():
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[GroundMotionProviderConfig(
            id="scaled",
            provider="gmpe",
            default=True,
            config={
                "gmpe_name": "AbrahamsonEtAl2014",
                "source_geometry": "finite",
                "source_parameters": SOURCE,
                "source_scaling": {"model": "ThingbaijamEtAl2017"},
            },
        )],
    )
    model.validate()
    runtime = GroundMotionRuntime(model, event())
    runtime.provider("scaled").evaluate(
        "PGA", SITES, runtime.event, site_parameters=SITE,
    )
    trace = runtime.source_trace()
    assert trace["type"] == "ShakeScenarioSourceTrace"
    assert len(trace["providers"]["scaled"]["evaluations"]) == 1
    assert trace["providers"]["scaled"]["backend"] == "gmpe"
    json.dumps(trace, allow_nan=False)
