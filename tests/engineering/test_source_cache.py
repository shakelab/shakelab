"""Source resolution caching at GMPE provider and runtime boundaries."""

from concurrent.futures import ThreadPoolExecutor

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
    GroundMotionRuntime,
)
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from tests.engineering.test_finite_source_scaling_provider import (
    SITE,
    SOURCE,
    SITES,
)


def provider(source=None):
    return GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "finite",
            "source_parameters": SOURCE if source is None else source,
            "source_scaling": {"model": "ThingbaijamEtAl2017"},
        },
    )


def event(magnitude=6.2, source=None, longitude=13.0):
    return ScenarioEvent.from_lonlat_depth(
        magnitude, longitude, 46.0, 7.0,
        source_parameters=source,
    )


def count_builds(monkeypatch, gmpe):
    original = gmpe._build_scaled_rupture
    calls = []

    def observed(ev, params):
        calls.append(1)
        return original(ev, params)

    monkeypatch.setattr(gmpe, "_build_scaled_rupture", observed)
    return calls


def test_same_source_reused_across_imts_and_sites(monkeypatch):
    gmpe = provider()
    calls = count_builds(monkeypatch, gmpe)
    ev = event()
    gmpe.evaluate("PGA", SITES, ev, site_parameters=SITE)
    gmpe.evaluate("SA-0.20", SITES[:1], ev, site_parameters=SITE)
    gmpe.evaluate("PGA", SITES, ev,
                  site_parameters={**SITE, "vs30": 600.0})
    assert len(calls) == 1


def test_distance_overrides_do_not_change_cache_key(monkeypatch):
    gmpe = provider()
    calls = count_builds(monkeypatch, gmpe)
    ev = event()
    initial = gmpe._resolve_source(ev, {})
    altered = gmpe._resolve_source(ev, {})
    assert initial is altered
    gmpe.evaluate("PGA", SITES, ev, site_parameters=SITE,
                  distances={"rrup": 42.0})
    assert len(calls) == 1


def test_overrides_and_event_data_invalidate_cache(monkeypatch):
    gmpe = provider()
    calls = count_builds(monkeypatch, gmpe)
    ev = event()
    base = gmpe._resolve_source(ev, {})
    other = gmpe._resolve_source(ev, {"dip": 40.0})
    same = gmpe._resolve_source(ev, {"dip": 40.0})
    assert base is not other
    assert other is same
    assert other.rupture.dip == pytest.approx(40.0)
    gmpe._resolve_source(event(magnitude=6.3), {})
    gmpe._resolve_source(event(longitude=13.1), {})
    gmpe._resolve_source(event(source={"rake": -90.0}), {})
    assert len(calls) == 5


def test_equivalent_effective_parameters_share_cache(monkeypatch):
    gmpe = provider()
    calls = count_builds(monkeypatch, gmpe)
    ev = event()
    first = gmpe._resolve_source(ev, {"dip": 40.0})
    second = gmpe._resolve_source(event(source={"dip": 40.0}), {})
    assert first is second
    assert len(calls) == 1


def test_concurrent_requests_build_only_once(monkeypatch):
    gmpe = provider()
    calls = count_builds(monkeypatch, gmpe)
    ev = event()
    with ThreadPoolExecutor(max_workers=6) as executor:
        sources = list(executor.map(
            lambda _: gmpe._resolve_source(ev, {}),
            range(12),
        ))
    assert all(source is sources[0] for source in sources)
    assert len(calls) == 1


def test_cache_capacity_is_bounded():
    gmpe = provider()
    for index in range(gmpe._SOURCE_CACHE_SIZE + 5):
        gmpe._resolve_source(
            event(),
            {"strike": float(index)},
        )
    assert len(gmpe._source_cache) == gmpe._SOURCE_CACHE_SIZE


def test_provider_instances_are_isolated(monkeypatch):
    first = provider()
    second = provider()
    calls_first = count_builds(monkeypatch, first)
    calls_second = count_builds(monkeypatch, second)
    first._resolve_source(event(), {})
    second._resolve_source(event(), {})
    assert len(calls_first) == len(calls_second) == 1
    assert first._resolve_source(event(), {}) is not (
        second._resolve_source(event(), {})
    )


def test_shakescenario_runtime_owns_provider_cache(monkeypatch):
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[GroundMotionProviderConfig(
            id="a",
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
    one = GroundMotionRuntime(model, event())
    two = GroundMotionRuntime(model, event())
    first_provider = one.provider("a")
    second_provider = two.provider("a")
    assert first_provider is not second_provider
    first_provider._resolve_source(one.event, {})
    assert len(first_provider._source_cache) == 1
    assert len(second_provider._source_cache) == 0


def test_native_and_point_geometry_not_cached():
    for mode in ("native", "point"):
        gmpe = GroundMotionProvider.gmpe(
            "AbrahamsonEtAl2014",
            config={"source_geometry": mode, "source_parameters": SOURCE},
        )
        gmpe._resolve_source(event(), {})
        assert len(gmpe._source_cache) == 0
