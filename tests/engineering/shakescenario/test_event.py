"""Regression tests for submitted ShakeScenario event parameters."""

from __future__ import annotations

import pytest

from shakelab.engineering.shakescenario.shakeserver import (
    _build_scenario_event,
)


@pytest.fixture
def event_data():
    return {
        "magnitude": 5.2,
        "hypocentre": {
            "longitude": 12.34,
            "latitude": 45.67,
            "elevation": -10000.0,
        },
        "event_id": "test-event",
        "origin_time": "2026-10-28T10:00:00Z",
    }


def test_focal_mechanism_is_preserved(event_data):
    mechanism = {"rake": 90.0, "strike": 45.0, "dip": 30.0}
    event_data["mechanism"] = mechanism
    event = _build_scenario_event(event_data)
    assert event.mechanism == mechanism
    assert event.magnitude == 5.2
    assert event.event_id == "test-event"
    assert event.origin_time == "2026-10-28T10:00:00Z"
    assert event.depth_km == pytest.approx(10.0)


def test_missing_mechanism_remains_none(event_data):
    event = _build_scenario_event(event_data)
    assert event.mechanism is None


def test_null_mechanism_remains_none(event_data):
    event_data["mechanism"] = None
    assert _build_scenario_event(event_data).mechanism is None


@pytest.mark.parametrize("mechanism", [42, "reverse", [90.0]])
def test_invalid_mechanism_is_rejected(event_data, mechanism):
    event_data["mechanism"] = mechanism
    with pytest.raises(ValueError, match="scenario.event.mechanism"):
        _build_scenario_event(event_data)
