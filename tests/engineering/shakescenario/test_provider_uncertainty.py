"""Regression tests for ShakeMap uncertainty semantics."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint


@pytest.fixture
def event():
    return ScenarioEvent.from_lonlat_depth(
        magnitude=5.0,
        longitude=13.0,
        latitude=46.0,
        depth_km=10.0,
        event_id="test-event",
    )


@pytest.fixture
def site():
    return WgsPoint(longitude=13.0, latitude=46.0, elevation=0.0)


def _fake_shakemap(sigma):
    """Return an in-memory ShakeMap exposing the required grid interface."""
    fields = {"PGA"}
    if sigma is not None:
        fields.add("STDPGA")

    class Grid:
        field_names = sorted(fields)

        @staticmethod
        def has_param(name):
            return name in fields

    def get_ground_motion(coordinates, parameters, method, outside):
        row = {"PGA": 20.0}
        if sigma is not None:
            row["STDPGA"] = sigma
        return [dict(row) for _ in coordinates]

    return SimpleNamespace(
        gm_data=Grid(),
        units={"PGA": "%g"},
        get_ground_motion=get_ground_motion,
    )


@pytest.mark.parametrize(
    ("sigma", "expected"),
    [(None, "unknown"), (0.0, 0.0), (0.6, 0.6)],
)
def test_shakemap_sigma_semantics(tmp_path, monkeypatch, event, site,
                                  sigma, expected):
    provider = GroundMotionProvider.from_id(
        "shakemap", root_path=str(tmp_path),
    )
    monkeypatch.setattr(
        provider, "_load_event", lambda event_id: _fake_shakemap(sigma),
    )
    values, sigmas = provider.evaluate("PGA", [site], event)
    assert values[0] == pytest.approx(0.2)
    if expected == "unknown":
        assert sigmas[0] != sigmas[0]
    else:
        assert sigmas[0] == pytest.approx(expected)


def test_shakemap_required_sigma_rejects_missing(
    tmp_path, monkeypatch, event, site,
):
    provider = GroundMotionProvider.from_id(
        "shakemap", root_path=str(tmp_path), require_uncertainty=True,
    )
    monkeypatch.setattr(
        provider, "_load_event", lambda event_id: _fake_shakemap(None),
    )
    with pytest.raises(KeyError, match="uncertainty field"):
        provider.evaluate("PGA", [site], event)
