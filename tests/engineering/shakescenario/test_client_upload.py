"""Client-side tests for one-product ground-motion uploads."""

import hashlib
import socket
import threading
from types import SimpleNamespace
import zipfile

import pytest

from shakelab.engineering.shakescenario.protocol import (
    ProtocolClient,
    recv_message,
    send_message,
)
from shakelab.engineering.shakescenario.shakeclient import (
    ShakeScenarioClient,
    _build_parser,
)


def test_upload_flags_are_mutually_exclusive():
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([
            "submit", "--shakemap", "map.zip",
            "--imt-file", "groundmotion.json",
        ])
    args = parser.parse_args([
        "submit", "--shakemap", "map.zip",
    ])
    assert args.shakemap == "map.zip"
    assert args.imt_file is None


def test_api_rejects_two_products():
    client = ShakeScenarioClient("127.0.0.1", 6000)
    with pytest.raises(ValueError, match="either"):
        client.submit({}, shakemap="map.zip", imt_file="data.json")


def test_client_upload_imt_descriptor(tmp_path, monkeypatch):
    source = tmp_path / "custom.json"
    content = b'{"type":"ShakeScenarioIMT"}'
    source.write_bytes(content)
    client = ShakeScenarioClient("127.0.0.1", 6000)
    calls = []

    def capture(payload, req_id, path, kind, fmt, filename=None):
        calls.append((payload, path.read_bytes(), kind, fmt, filename))
        return {"job_id": 3, "status": "queued"}

    client._client = SimpleNamespace(request_attachment=capture)
    result = client.submit({"tag": "x"}, imt_file=source)
    assert result["job_id"] == 3
    assert calls == [
        ({"tag": "x"}, content, "imt", "json", None)
    ]


def test_client_zips_shakemap_directory(tmp_path, monkeypatch):
    directory = tmp_path / "event"
    directory.mkdir()
    (directory / "grid.xml").write_text("<grid/>")
    (directory / "uncertainty.xml").write_text("<uncertainty/>")
    client = ShakeScenarioClient("127.0.0.1", 6000)
    seen = []

    def capture(payload, req_id, path, kind, fmt, filename=None):
        with zipfile.ZipFile(path) as archive:
            seen.append((sorted(archive.namelist()), kind, fmt, filename))
        return {"job_id": 4, "status": "queued"}

    client._client = SimpleNamespace(request_attachment=capture)
    result = client.submit({}, shakemap=directory)
    assert result["job_id"] == 4
    assert seen == [
        (
            ["grid.xml", "uncertainty.xml"],
            "shakemap", "zip", "event.zip",
        )
    ]


def test_client_reuses_legacy_json_submit(monkeypatch):
    client = ShakeScenarioClient("127.0.0.1", 6000)
    seen = []

    def capture(op, payload, req_id):
        seen.append((op, payload))
        return {"job_id": 9, "status": "queued"}

    client._client = SimpleNamespace(request=capture)
    assert client.submit({"tag": "old"})["job_id"] == 9
    assert seen == [("submit", {"tag": "old"})]


def test_protocol_streams_binary_after_json(tmp_path, monkeypatch):
    content = b"shakemap-binary-data" * 3000
    source = tmp_path / "map.zip"
    source.write_bytes(content)
    client = ProtocolClient("127.0.0.1", 6000)
    receiver, sender = socket.socketpair()
    observed = {}

    def handle():
        try:
            request = recv_message(receiver)
            descriptor = request["attachment"]
            remaining = descriptor["size"]
            digest = hashlib.sha256()
            while remaining:
                block = receiver.recv(min(4096, remaining))
                if not block:
                    raise AssertionError("Unexpected EOF")
                digest.update(block)
                remaining -= len(block)
            observed["request"] = request
            observed["digest"] = digest.hexdigest()
            send_message(receiver, {
                "v": 1, "req_id": request["req_id"],
                "ok": True,
                "result": {"job_id": 12, "status": "queued"},
            })
        finally:
            receiver.close()

    worker = threading.Thread(target=handle)
    worker.start()
    monkeypatch.setattr(
        socket, "create_connection",
        lambda address, timeout: sender,
    )
    try:
        result = client.request_attachment(
            {}, "request-1", source, "shakemap", "zip",
        )
    finally:
        worker.join(timeout=5)
    assert not worker.is_alive()
    assert result["job_id"] == 12
    assert observed["request"]["attachment"]["sha256"] == (
        hashlib.sha256(content).hexdigest()
    )
    assert observed["digest"] == hashlib.sha256(content).hexdigest()


def test_python_scenario_api_accepts_uploaded_file(tmp_path, monkeypatch):
    source = tmp_path / "data.json"
    source.write_text("{}")
    client = ShakeScenarioClient("127.0.0.1", 6000)
    seen = []

    def capture(payload, *, shakemap=None, imt_file=None):
        seen.append((payload, shakemap, imt_file))
        return {"job_id": 14}

    monkeypatch.setattr(client, "submit", capture)
    result = client.submit_scenario(
        magnitude=6.2, longitude=13.25, latitude=46.25,
        depth_km=8.0, model_id="nei", imt_file=source,
    )
    assert result["job_id"] == 14
    assert seen[0][2] == source
    assert seen[0][0]["scenario"]["event"]["magnitude"] == 6.2
