"""Regression tests for ShakeScenario submission entry points."""

from __future__ import annotations

import json
import sys

import pytest

from shakelab.engineering.shakescenario import shakeclient


def test_cli_json_and_python_submit_equivalent_requests(
    tmp_path, monkeypatch,
):
    requests = []

    def request(self, op, payload, req_id):
        assert op == "submit"
        requests.append(payload)
        return {"job_id": 1, "status": "queued"}

    monkeypatch.setattr(shakeclient.ProtocolClient, "request", request)
    monkeypatch.setattr(sys, "argv", [
        "shakeclient", "--format", "json", "submit",
        "--mag", "6.2", "--lon", "13.25", "--lat", "46.25",
        "--depth", "10", "--rake", "-90", "--dip", "35",
        "--source", "ztor=2", "width=12.5",
    ])
    shakeclient.main()

    expected = {"scenario": {"event": {
        "magnitude": 6.2,
        "hypocentre": {
            "longitude": 13.25,
            "latitude": 46.25,
            "elevation": -10000.0,
        },
        "source_parameters": {
            "rake": -90.0, "dip": 35.0, "ztor": 2.0, "width": 12.5,
        },
    }}}
    path = tmp_path / "request.json"
    path.write_text(json.dumps(expected), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", [
        "shakeclient", "--format", "json", "submit",
        "--request-json", str(path),
    ])
    shakeclient.main()

    client = shakeclient.ShakeScenarioClient("localhost", 6000)
    client.submit_scenario(
        magnitude=6.2,
        longitude=13.25,
        latitude=46.25,
        elevation_m=-10000.0,
        source_parameters={"rake": -90.0, "dip": 35.0, "ztor": 2.0, "width": 12.5},
    )

    assert requests == [expected, expected, expected]


def test_cli_overrides_individual_json_fields(tmp_path):
    path = tmp_path / "request.json"
    path.write_text(json.dumps({
        "models": {"model_id": "test-model"},
        "scenario": {"event": {
            "magnitude": 6.2,
            "hypocentre": {"longitude": 13.25, "elevation": -10000.0},
            "source_parameters": {"rake": 90, "dip": 35, "ztor": 2},
        }},
    }), encoding="utf-8")
    args = shakeclient._build_parser().parse_args([
        "submit", "--request-json", str(path), "--rake", "0",
        "--elevation", "0", "--source", "ztor=3",
    ])
    result = shakeclient._merge_submit_payload(args)

    assert result["models"]["model_id"] == "test-model"
    event = result["scenario"]["event"]
    assert event["magnitude"] == 6.2
    assert event["hypocentre"] == {"longitude": 13.25, "elevation": 0.0}
    assert event["source_parameters"] == {"rake": 0.0, "dip": 35, "ztor": 3.0}


@pytest.mark.parametrize(("options", "message"), [
    (["--rake", "90", "--source", "rake=45"], "Duplicate.*rake"),
    (["--strike", "0", "--strike", "90"], "Duplicate.*strike"),
    (["--source", "ztor=2", "ztor=3"], "Duplicate.*ztor"),
    (["--source", "ztor=2", "--source", "width=3"], "once"),
    (["--source", "90,35"], "NAME=VALUE"),
    (["--source", "ztor=abc"], "finite"),
    (["--source", "ztor=nan"], "finite"),
])
def test_cli_rejects_invalid_mechanism_options(options, message):
    args = shakeclient._build_parser().parse_args(["submit", *options])
    with pytest.raises(ValueError, match=message):
        shakeclient._merge_submit_payload(args)


def test_cli_and_python_reject_both_vertical_coordinates():
    with pytest.raises(SystemExit) as error:
        shakeclient._build_parser().parse_args([
            "submit", "--depth", "10", "--elevation", "-10000",
        ])
    assert error.value.code == 2

    client = shakeclient.ShakeScenarioClient("localhost", 6000)
    with pytest.raises(ValueError, match="only one"):
        client.submit_scenario(
            magnitude=6.2, longitude=13.25, latitude=46.25,
            depth_km=10, elevation_m=-10000,
        )
