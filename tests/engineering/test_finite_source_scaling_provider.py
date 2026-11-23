"""GMPE integration checks for automatic finite rupture scaling."""

from math import isfinite

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.groundmotion.rupture import build_rectangular_rupture
from tests.engineering.test_finite_source_provider import (
    SITE,
    SITES,
)


SOURCE = {
    "strike": 120.0,
    "dip": 35.0,
    "rake": 90.0,
    "aspect_ratio": 2.0,
    "lower_seismogenic_depth": 15.0,
}


def provider(model="ThingbaijamEtAl2017", source=None, **scaling):
    return GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "finite",
            "source_parameters": SOURCE if source is None else source,
            "source_scaling": {"model": model, **scaling},
        },
    )


def event(source=None, geometry=None, magnitude_type=None):
    return ScenarioEvent.from_lonlat_depth(
        6.2, 13.0, 46.0, 7.0, source_parameters=source,
        source_geometry=geometry, magnitude_type=magnitude_type,
    )


@pytest.mark.parametrize("model", [
    "WellsCoppersmith1994",
    "ThingbaijamEtAl2017",
    "Leonard2014_Interplate",
    "Leonard2014_StableContinental",
    "HanksBakun2002",
])
def test_five_models_build_finite_gmpe_geometry(model):
    gmpe = provider(model)
    resolved = gmpe._resolve_source(event(), {})
    assert resolved.build_result is not None
    assert resolved.build_result.scaling_model == model
    assert resolved.build_result.parameter_origins["area"] == "scaling"
    assert resolved.rupture.area == pytest.approx(
        resolved.build_result.scaling_estimate.median,
    )
    values, sigmas = gmpe.evaluate("PGA", SITES, event(),
                                   site_parameters=SITE)
    assert len(values) == len(SITES)
    assert all(isfinite(value) and value > 0 for value in values)
    assert all(isfinite(sigma) and sigma >= 0 for sigma in sigmas)


def test_shared_source_matches_builder_and_rupture_context():
    gmpe = provider()
    ev = event()
    resolved = gmpe._resolve_source(ev, {})
    expected = build_rectangular_rupture(
        ev.hypocentre, magnitude=ev.magnitude, **SOURCE,
        scaling_model="ThingbaijamEtAl2017",
    )
    context = gmpe._build_rupture_context(ev, resolved)
    for key in ("width", "ztor", "dip"):
        assert resolved.parameters[key] == pytest.approx(
            getattr(expected.rupture, key),
        )
        assert getattr(context, key) == pytest.approx(
            getattr(expected.rupture, key),
        )
    distances = gmpe._build_sites_context(
        SITES, ev, SITE, {}, {}, resolved_source=resolved,
    )
    expected_distances = expected.rupture.distances(SITES)
    for metric in ("rrup", "rjb", "rx", "ry0"):
        if hasattr(distances, metric):
            assert getattr(distances, metric) == pytest.approx(
                expected_distances[metric],
            )


def test_seismogenic_width_cap_preserves_area_and_context():
    source = {
        **SOURCE,
        "upper_seismogenic_depth": 6.0,
        "lower_seismogenic_depth": 8.0,
    }
    gmpe = provider(source=source)
    resolved = gmpe._resolve_source(event(), {})
    assert resolved.rupture.area == pytest.approx(
        resolved.build_result.scaling_estimate.median,
    )
    assert resolved.build_result.parameter_origins["width"] == "constrained"
    context = gmpe._build_rupture_context(event(), resolved)
    assert context.width == pytest.approx(resolved.rupture.width)
    assert context.ztor == pytest.approx(resolved.rupture.ztor)


def test_source_precedence_with_scaling_and_distance_override():
    defaults = {**SOURCE, "strike": 45.0}
    gmpe = provider(source=defaults)
    ev = event({"strike": 90.0})
    resolved = gmpe._resolve_source(ev, {"strike": 120.0})
    assert resolved.rupture.strike == pytest.approx(120.0)
    site = gmpe._build_sites_context(
        SITES, ev, SITE, {"rrup": 42.0},
        {"strike": 120.0}, resolved_source=resolved,
    )
    assert site.rrup == pytest.approx([42.0] * len(SITES))


def test_explicit_dimensions_bypass_scaling_estimate():
    source = {
        "strike": 120.0, "dip": 35.0, "rake": 90.0,
        "length": 20.0, "width": 12.0,
        "ztor": 2.0, "along_fraction": 0.5,
    }
    gmpe = provider(source=source)
    resolved = gmpe._resolve_source(event(), {})
    assert resolved.rupture.length == pytest.approx(20.0)
    assert resolved.rupture.width == pytest.approx(12.0)
    assert resolved.build_result.scaling_estimate is None


def test_conflicting_explicit_geometry_fails():
    source = {**SOURCE, "length": 20.0, "width": 12.0}
    with pytest.raises(ValueError, match="aspect_ratio"):
        provider(source=source)._resolve_source(event(), {})


def test_missing_strike_or_dip_fails():
    for missing in ("strike", "dip"):
        source = {k: v for k, v in SOURCE.items() if k != missing}
        with pytest.raises(ValueError, match=missing):
            provider(source=source)._resolve_source(event(), {})


def test_thingbaijam_auto_without_rake_fails():
    source = {k: v for k, v in SOURCE.items() if k != "rake"}
    with pytest.raises(ValueError, match="Unknown fault style"):
        provider(source=source)._resolve_source(event(), {})


def test_explicit_interface_style_does_not_require_rake():
    source = {k: v for k, v in SOURCE.items() if k != "rake"}
    resolved = provider(source=source, style="interface")._resolve_source(
        event(), {},
    )
    assert resolved.build_result.scaling_style == "interface"


def test_event_geometry_override_disables_scaling():
    gmpe = provider()
    for geometry in ("native", "point"):
        resolved = gmpe._resolve_source(event(geometry=geometry), {})
        assert resolved.rupture is None
        assert resolved.build_result is None


def test_magnitude_labels_use_same_numerical_value():
    gmpe = provider()
    mw = gmpe.evaluate(
        "PGA", SITES, event(magnitude_type="Mw"),
        site_parameters=SITE,
    )
    ml = gmpe.evaluate(
        "PGA", SITES, event(magnitude_type="ML"),
        site_parameters=SITE,
    )
    for observed, reference in zip(ml, mw):
        assert observed == pytest.approx(reference)


def test_scaling_rupture_built_once_per_evaluation(monkeypatch):
    gmpe = provider()
    original = gmpe._build_scaled_rupture
    calls = []

    def observed(event_, source):
        calls.append(dict(source))
        return original(event_, source)

    monkeypatch.setattr(gmpe, "_build_scaled_rupture", observed)
    gmpe.evaluate("PGA", SITES, event(), site_parameters=SITE)
    assert len(calls) == 1
