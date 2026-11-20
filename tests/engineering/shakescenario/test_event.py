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
    source_parameters = {"rake": 90.0, "strike": 45.0, "dip": 30.0}
    event_data["source_parameters"] = source_parameters
    event = _build_scenario_event(event_data)
    assert event.source_parameters == source_parameters
    assert event.magnitude == 5.2
    assert event.event_id == "test-event"
    assert event.origin_time == "2026-10-28T10:00:00Z"
    assert event.depth_km == pytest.approx(10.0)


def test_missing_mechanism_remains_none(event_data):
    event = _build_scenario_event(event_data)
    assert event.source_parameters is None


def test_null_mechanism_remains_none(event_data):
    event_data["source_parameters"] = None
    assert _build_scenario_event(event_data).source_parameters is None


@pytest.mark.parametrize("source_parameters", [42, "reverse", [90.0]])
def test_invalid_mechanism_is_rejected(event_data, source_parameters):
    event_data["source_parameters"] = source_parameters
    with pytest.raises(ValueError, match="scenario.event.source_parameters"):
        _build_scenario_event(event_data)


@pytest.mark.parametrize("parameters", [
    {"mag": 6.0}, {"hypo_depth": 10.0}, {"strike": 361.0},
    {"dip": -1.0}, {"rake": 181.0}, {"width": float("nan")},
    {"width": True}, {"width": "12"},
])
def test_invalid_source_parameter_values(event_data, parameters):
    event_data["source_parameters"] = parameters
    with pytest.raises(ValueError):
        _build_scenario_event(event_data)
