"""Regression tests for single-source GMPE resolution."""

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint
from tests.engineering.test_finite_source_provider import (
    SITE,
    SOURCE,
    SITES,
)


def _provider(geometry="finite", source=None):
    return GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": geometry,
            "source_parameters": source,
        },
    )


def _event(source=None, geometry=None):
    return ScenarioEvent.from_lonlat_depth(
        6.2, 13.0, 46.0, 7.0,
        source_parameters=source,
        source_geometry=geometry,
    )


def test_explicit_rupture_constructed_once_per_evaluation(monkeypatch):
    provider = _provider(source=SOURCE)
    original = provider._build_finite_rupture
    calls = []

    def observed(event, parameters):
        calls.append(dict(parameters))
        return original(event, parameters)

    monkeypatch.setattr(provider, "_build_finite_rupture", observed)
    provider.evaluate("PGA", SITES, _event(), site_parameters=SITE)
    assert len(calls) == 1


def test_resolved_geometry_and_rupture_context_agree():
    provider = _provider(source=SOURCE)
    event = _event()
    resolved = provider._resolve_source(event, {})
    context = provider._build_rupture_context(event, resolved)

    assert resolved.geometry == "finite"
    for name in ("width", "ztor", "dip"):
        if hasattr(context, name):
            assert getattr(context, name) == pytest.approx(
                getattr(resolved.rupture, name)
            )
    with pytest.raises(TypeError):
        resolved.parameters["dip"] = 0.0


def test_runtime_override_flows_to_geometry_and_gmpe_context():
    provider = _provider(source=SOURCE)
    event = _event({"strike": 45.0})
    source = provider._resolve_source(event, {"dip": 40.0})
    context = provider._build_rupture_context(event, source)
    assert source.rupture.strike == pytest.approx(45.0)
    assert source.rupture.dip == pytest.approx(40.0)
    assert context.dip == pytest.approx(40.0)


def test_distance_overrides_do_not_mutate_source():
    provider = _provider(source=SOURCE)
    event = _event()
    source = provider._resolve_source(event, {})
    site = provider._build_sites_context(
        SITES, event, SITE, {"rrup": 42.0}, {},
        resolved_source=source,
    )
    assert site.rrup == pytest.approx([42.0, 42.0])
    assert source.rupture.distances(SITES)["rrup"] != pytest.approx(
        [42.0, 42.0]
    )


@pytest.mark.parametrize("mode", ["native", "point"])
def test_nonfinite_modes_do_not_construct_rupture(mode, monkeypatch):
    provider = _provider(mode, SOURCE)

    def fail(*args):
        raise AssertionError("Unexpected finite source construction")

    monkeypatch.setattr(provider, "_build_finite_rupture", fail)
    source = provider._resolve_source(_event(), {})
    assert source.geometry == mode
    assert source.rupture is None


def test_separate_evaluations_resolve_separate_snapshots():
    provider = _provider(source=SOURCE)
    event = _event()
    first = provider._resolve_source(event, {"dip": 40.0})
    second = provider._resolve_source(event, {"dip": 50.0})
    assert first.rupture.dip == 40.0
    assert second.rupture.dip == 50.0
    assert first.parameters["dip"] == 40.0
