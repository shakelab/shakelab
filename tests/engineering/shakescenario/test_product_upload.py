"""Protocol and staging regression tests for single-product submits."""

import hashlib
import io
import json
import socket
import zipfile
from types import SimpleNamespace

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
)
from shakelab.engineering.shakescenario.database import JobDatabase
from shakelab.engineering.shakescenario.product_upload import (
    Attachment,
    MAX_ATTACHMENT_BYTES,
    extract_shakemap,
    receive_attachment,
)
from shakelab.engineering.shakescenario.shakeserver import (
    ShakeScenarioServer,
)


def _descriptor(content, kind="imt", fmt="json"):
    return {
        "kind": kind,
        "format": fmt,
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "filename": "groundmotion.json" if kind == "imt" else "map.zip",
    }


def _imt_data():
    return json.dumps({
        "type": "ShakeScenarioIMT",
        "schema_version": "1.0.0",
        "imts": {"PGA": {"units": "g"}},
        "stations": [{
            "station_code": "A",
            "values": {"PGA": 0.1},
        }],
    }).encode("utf-8")


def _payload():
    return {
        "models": {"model_id": "test"},
        "scenario": {"event": {
            "magnitude": 5.0,
            "hypocentre": {"longitude": 13.0, "latitude": 46.0},
        }},
    }


def _server(tmp_path, monkeypatch, kind="imt", **provider_options):
    from shakelab.engineering.shakescenario import shakeserver

    root = tmp_path / "root"
    root.mkdir()
    jobs = tmp_path / "runs"
    db = JobDatabase(tmp_path / "jobs.db")
    db.initialize()
    cfg = SimpleNamespace(
        defaults=SimpleNamespace(model_id="test", impact_config={}),
        paths=SimpleNamespace(model_root=root),
    )
    provider = GroundMotionProviderConfig(
        id="uploaded",
        provider=kind,
        config={"root_path": str(root), **provider_options},
    )
    ground_motion = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[provider],
    )
    monkeypatch.setattr(
        shakeserver, "load_model",
        lambda *args: SimpleNamespace(ground_motion=ground_motion),
    )
    monkeypatch.setattr(shakeserver, "model_paths", lambda *args: {})
    server = ShakeScenarioServer(
        "127.0.0.1", 0, db, jobs, cfg
    )
    monkeypatch.setattr(server._executor, "submit", lambda *args: None)
    return server, db


def _transmit(server, payload, descriptor, content):
    receiver, sender = socket.socketpair()
    try:
        sender.sendall(content)
        sender.shutdown(socket.SHUT_WR)
        return server._op_submit_attachment(
            payload, descriptor, receiver
        )
    finally:
        receiver.close()
        sender.close()


def test_submit_imt_without_event_id(tmp_path, monkeypatch):
    server, db = _server(tmp_path, monkeypatch)
    data = _imt_data()
    result = _transmit(
        server, _payload(), _descriptor(data), data
    )
    job_id = result["job_id"]
    folder = tmp_path / "runs" / f"job_{job_id:06d}"
    assert (folder / "inputs/imt/groundmotion.json").read_bytes() == data
    request = json.loads((folder / "request.json").read_text())
    assert request["_product"]["provider_id"] == "uploaded"
    assert request["_product"]["path"] == "imt/groundmotion.json"
    assert request["scenario"]["event"].get("event_id") is None
    assert db.get_job(job_id)["status"] == "queued"
    assert not list((tmp_path / "runs").glob(".upload_*"))
    server._executor.shutdown(wait=False)


@pytest.mark.parametrize("failure", ["checksum", "interrupted"])
def test_failed_transfer_never_queues(tmp_path, monkeypatch, failure):
    server, db = _server(tmp_path, monkeypatch)
    data = _imt_data()
    descriptor = _descriptor(data)
    if failure == "checksum":
        descriptor["sha256"] = "0" * 64
    else:
        descriptor["size"] += 10
    with pytest.raises((ValueError, EOFError, RuntimeError)):
        _transmit(server, _payload(), descriptor, data)
    assert db.list_jobs() == []
    assert not list((tmp_path / "runs").glob(".upload_*"))
    server._executor.shutdown(wait=False)


def test_invalid_imt_json_never_queues(tmp_path, monkeypatch):
    server, db = _server(tmp_path, monkeypatch)
    data = b'{"type":'
    with pytest.raises((ValueError, json.JSONDecodeError)):
        _transmit(server, _payload(), _descriptor(data), data)
    assert db.list_jobs() == []
    server._executor.shutdown(wait=False)


def test_descriptor_rejects_oversize_and_invalid_kind():
    descriptor = _descriptor(b"x")
    descriptor["size"] = MAX_ATTACHMENT_BYTES + 1
    with pytest.raises(ValueError, match="size"):
        Attachment.parse(descriptor)
    descriptor["size"] = 1
    descriptor["kind"] = "gmpe"
    with pytest.raises(ValueError, match="kind"):
        Attachment.parse(descriptor)


@pytest.mark.parametrize("name", ["../grid.xml", "/grid.xml",
                                    "folder/../../grid.xml",
                                    "folder\\grid.xml"])
def test_zip_rejects_dangerous_paths(tmp_path, name):
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as output:
        output.writestr(name, "unsafe")
    with pytest.raises(ValueError, match="Unsafe"):
        extract_shakemap(archive, tmp_path / "unpacked")


def test_zip_accepts_single_enclosing_directory(tmp_path):
    archive = tmp_path / "map.zip"
    with zipfile.ZipFile(archive, "w") as output:
        output.writestr("event/grid.xml", "<grid/>")
        output.writestr("event/uncertainty.xml", "<grid/>")
    dest = tmp_path / "unpacked"
    extract_shakemap(archive, dest)
    assert (dest / "grid.xml").is_file()
    assert (dest / "uncertainty.xml").is_file()


def test_missing_grid_rejected_before_queue(tmp_path, monkeypatch):
    server, db = _server(tmp_path, monkeypatch, kind="shakemap")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("info.json", "{}")
    data = buffer.getvalue()
    with pytest.raises(FileNotFoundError, match="grid"):
        _transmit(
            server, _payload(), _descriptor(data, "shakemap", "zip"), data
        )
    assert db.list_jobs() == []
    server._executor.shutdown(wait=False)


def test_required_uncertainty_rejected(tmp_path, monkeypatch):
    server, db = _server(
        tmp_path, monkeypatch, kind="shakemap",
        require_uncertainty=True
    )
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("grid.xml", "<grid/>")
    data = buffer.getvalue()
    with pytest.raises(ValueError, match="uncertainty"):
        _transmit(
            server, _payload(), _descriptor(data, "shakemap", "zip"), data
        )
    assert db.list_jobs() == []
    server._executor.shutdown(wait=False)


def test_framed_tcp_upload_submit(tmp_path, monkeypatch):
    """The binary body follows the existing framed JSON request."""
    import threading

    from shakelab.engineering.shakescenario.protocol import (
        recv_message, send_message,
    )

    server, db = _server(tmp_path, monkeypatch)
    data = _imt_data()
    receiver, sender = socket.socketpair()
    worker = threading.Thread(
        target=server._handle_connection,
        args=(receiver, ("local", 0)),
    )
    worker.start()
    try:
        send_message(sender, {
            "v": 1,
            "op": "submit",
            "req_id": "upload-test",
            "payload": _payload(),
            "attachment": _descriptor(data),
        })
        sender.sendall(data)
        response = recv_message(sender)
        assert response["ok"] is True
        assert response["req_id"] == "upload-test"
        job_id = response["result"]["job_id"]
        assert db.get_job(job_id)["status"] == "queued"
    finally:
        sender.close()
        worker.join(timeout=5)
        server._executor.shutdown(wait=False)
    assert not worker.is_alive()


def test_framed_tcp_legacy_submit(tmp_path, monkeypatch):
    """Requests without attachments retain the old JSON framing."""
    import threading

    from shakelab.engineering.shakescenario.protocol import (
        recv_message, send_message,
    )

    server, db = _server(tmp_path, monkeypatch)
    receiver, sender = socket.socketpair()
    worker = threading.Thread(
        target=server._handle_connection,
        args=(receiver, ("local", 0)),
    )
    worker.start()
    try:
        send_message(sender, {
            "v": 1,
            "op": "submit",
            "req_id": "legacy-test",
            "payload": _payload(),
        })
        response = recv_message(sender)
        assert response["ok"] is True
        assert db.get_job(response["result"]["job_id"]) is not None
    finally:
        sender.close()
        worker.join(timeout=5)
        server._executor.shutdown(wait=False)
    assert not worker.is_alive()


def test_upload_rejects_ambiguous_provider(tmp_path, monkeypatch):
    from shakelab.engineering.shakescenario import shakeserver

    server, db = _server(tmp_path, monkeypatch)
    cfg = GroundMotionProviderConfig(
        id="second", provider="imt",
        config={"root_path": str(tmp_path / "root")},
    )
    model = GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="first", provider="imt",
                config={"root_path": str(tmp_path / "root")},
            ),
            cfg,
        ],
    )
    monkeypatch.setattr(
        shakeserver, "load_model",
        lambda *args: SimpleNamespace(ground_motion=model),
    )
    data = _imt_data()
    with pytest.raises(ValueError, match="exactly one matching"):
        _transmit(server, _payload(), _descriptor(data), data)
    assert db.list_jobs() == []
    server._executor.shutdown(wait=False)
