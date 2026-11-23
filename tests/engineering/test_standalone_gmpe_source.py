"""The standalone GMPE API does not depend on scenario trace output."""

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
)
from shakelab.groundmotion.providers import GroundMotionProvider
from tests.engineering.test_finite_source_scaling_provider import (
    SITE,
    SOURCE,
    SITES,
    event,
)


def make_provider(trace=False):
    return GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "finite",
            "source_parameters": SOURCE,
            "source_scaling": {"model": "ThingbaijamEtAl2017"},
            "record_source_trace": trace,
        },
    )


def test_generic_provider_skips_diagnostic_trace(monkeypatch):
    gmpe = make_provider()

    def unexpected(*args, **kwargs):
        raise AssertionError("Diagnostics must not run by default")

    monkeypatch.setattr(gmpe, "_record_source_trace", unexpected)
    values, sigmas = gmpe.evaluate(
        "PGA", SITES, event(), site_parameters=SITE,
    )
    assert all(value > 0.0 for value in values)
    assert all(sigma >= 0.0 for sigma in sigmas)
    assert gmpe.source_trace() == []


def test_public_resolve_source_preserves_physical_geometry():
    gmpe = make_provider()
    source = gmpe.resolve_source(event())
    assert source.geometry == "finite"
    assert source.build_result.scaling_model == "ThingbaijamEtAl2017"
    overridden = gmpe.resolve_source(event(), {"dip": 40.0})
    assert overridden.rupture.dip == pytest.approx(40.0)
    assert source.rupture.dip == pytest.approx(35.0)


def test_explicit_trace_remains_available():
    gmpe = make_provider(trace=True)
    gmpe.evaluate("PGA", SITES, event(), site_parameters=SITE)
    assert len(gmpe.source_trace()) == 1


def test_invalid_trace_option_rejected():
    with pytest.raises(TypeError, match="record_source_trace"):
        GroundMotionProvider.gmpe(
            "AbrahamsonEtAl2014",
            config={"record_source_trace": 1},
        )


def test_runtime_trace_is_optional_and_independent():
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
    without_trace = model.runtime(event())
    with_trace = model.runtime(event(), record_source_trace=True)
    for runtime in (without_trace, with_trace):
        runtime.provider("scaled").evaluate(
            "PGA", SITES, runtime.event, site_parameters=SITE,
        )
    assert without_trace.source_trace()["providers"] == {}
    assert len(with_trace.source_trace()["providers"]["scaled"][
        "evaluations"
    ]) == 1
