"""Integration tests for finite rectangular rupture providers."""

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.groundmotion.rupture import RectangularRupture
from shakelab.libutils.geodeticN.primitives import WgsPoint
from tests.engineering.test_impact_rupture_geometry import SITE


SOURCE = {
    "strike": 120.0, "dip": 35.0, "rake": 90.0,
    "length": 20.0, "width": 12.0, "ztor": 2.0,
    "along_fraction": 0.5,
}
SITES = [WgsPoint(13.15, 46.1), WgsPoint(13.0, 46.0)]


def _provider(geometry="finite", source=None):
    return GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={"source_geometry": geometry, "source_parameters": source},
    )


def _event(source=None, geometry=None):
    return ScenarioEvent.from_lonlat_depth(
        6.2, 13.0, 46.0, 7.0,
        source_parameters=source, source_geometry=geometry,
    )


def _explicit_distances(event, source):
    rupture = RectangularRupture.from_hypocentre(
        event.hypocentre,
        **{k: v for k, v in source.items() if k != "rake"},
    )
    return rupture.distances(SITES)


@pytest.mark.parametrize("geometry", [None, "finite"])
def test_finite_provider_matches_explicit_distances(geometry):
    event = _event(SOURCE, geometry)
    provider = _provider("finite" if geometry is None else "native")
    distances = _explicit_distances(event, SOURCE)
    context = provider._build_sites_context(SITES, event, SITE, {}, {})
    for key, expected in distances.items():
        assert getattr(context, key) == pytest.approx(expected)
    calculated = provider.evaluate("PGA", SITES, event, site_parameters=SITE)
    explicit = _provider("native").evaluate(
        "PGA", SITES, _event(SOURCE), site_parameters=SITE,
        distances=distances,
    )
    for actual, expected in zip(calculated, explicit):
        assert actual == pytest.approx(expected)


def test_provider_and_event_source_merge():
    event = _event({"strike": 45.0})
    provider = _provider(source=SOURCE)
    context = provider._build_sites_context(SITES, event, SITE, {}, {})
    expected = _explicit_distances(event, {**SOURCE, "strike": 45.0})
    for key, values in expected.items():
        assert getattr(context, key) == pytest.approx(values)


def test_explicit_distance_precedence():
    context = _provider(source=SOURCE)._build_sites_context(
        SITES, _event(), SITE, {"rrup": 42.0}, {},
    )
    assert context.rrup == pytest.approx([42.0, 42.0])


@pytest.mark.parametrize("missing", [
    "strike", "dip", "length", "width", "ztor", "along_fraction",
])
def test_missing_geometry_parameter_fails(missing):
    source = {key: value for key, value in SOURCE.items() if key != missing}
    with pytest.raises(ValueError, match=missing):
        _provider().evaluate("PGA", SITES, _event(source),
                             site_parameters=SITE)


def test_inconsistent_hypocentre_rejected():
    with pytest.raises(ValueError, match="outside"):
        _provider().evaluate(
            "PGA", SITES, _event({**SOURCE, "width": 1}),
            site_parameters=SITE,
        )


def test_point_override_does_not_require_finite_parameters():
    provider = _provider()
    event = _event({"strike": 120.0}, "point")
    context = provider._build_sites_context(SITES, event, SITE, {}, {})
    assert context.rjb.shape == (2,)
