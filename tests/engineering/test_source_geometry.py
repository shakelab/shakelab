"""Test geometry selection without changing distance algorithms."""

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint
from tests.engineering.test_impact_rupture_geometry import SITE, GEOMETRY


def _context(config, geometry=None, distances=None):
    provider = GroundMotionProvider.from_id("gmpe", config={
        "gmpe_name": "ChiouYoungs2014",
        **config,
    })
    event = ScenarioEvent.from_lonlat_depth(
        6.0, 13.0, 46.0, 10.0,
        source_parameters=GEOMETRY,
        source_geometry=geometry,
    )
    return provider._build_sites_context(
        [WgsPoint(13.15, 46.1)], event, SITE, distances or {}, {},
    )


@pytest.mark.parametrize("config", [
    {"source_geometry": "point"},
    {"point_source_distances": True},
    {"source_geometry": "native"},
    {},
])
def test_event_point_overrides_model(config):
    context = _context(config, "point")
    expected = _context({"source_geometry": "point"})
    for metric in ("rrup", "rjb", "rx"):
        assert getattr(context, metric) == pytest.approx(
            getattr(expected, metric),
        )


@pytest.mark.parametrize("config", [
    {"source_geometry": "point"},
    {"point_source_distances": True},
])
def test_event_native_disables_approximation(config):
    with pytest.raises(ValueError, match="distance"):
        _context(config, "native")


def test_explicit_distances_still_win():
    context = _context({}, "point", {"rrup": 42.0})
    assert context.rrup == pytest.approx([42.0])


@pytest.mark.parametrize("config", [
    {"source_geometry": "finite"},
    {"source_geometry": "unknown"},
    {"source_geometry": "point", "point_source_distances": True},
])
def test_invalid_model_geometry_rejected(config):
    with pytest.raises(ValueError):
        _context(config)
