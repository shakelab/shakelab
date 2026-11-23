# product_upload.py
# -*- coding: utf-8 -*-
"""Validation and staging of ShakeScenario ground-motion attachments."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import stat
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from shakelab.engineering.shakescenario.protocol import ProtocolError

MAX_ATTACHMENT_BYTES = 128 * 1024 * 1024
MAX_UNPACKED_BYTES = 512 * 1024 * 1024
MAX_ZIP_ENTRIES = 256
CHUNK_SIZE = 1024 * 1024
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class Attachment:
    """Metadata preceding one raw binary product in a submit request."""

    kind: str
    format: str
    size: int
    sha256: str
    filename: str

    @classmethod
    def parse(cls, value: Any) -> "Attachment":
        """Validate an attachment descriptor before reading its payload."""
        if not isinstance(value, dict):
            raise ValueError("attachment must be a JSON object.")
        if set(value) != {"kind", "format", "size", "sha256", "filename"}:
            raise ValueError("Invalid attachment descriptor fields.")
        kind, fmt = value["kind"], value["format"]
        if (kind, fmt) not in {("shakemap", "zip"), ("imt", "json")}:
            raise ValueError("Unsupported attachment kind or format.")
        size = value["size"]
        if type(size) is not int or not 0 < size <= MAX_ATTACHMENT_BYTES:
            raise ValueError("Attachment size exceeds limit or is invalid.")
        digest = value["sha256"]
        if not isinstance(digest, str) or not _SHA256_RE.fullmatch(digest):
            raise ValueError("attachment.sha256 must be lowercase SHA-256.")
        filename = value["filename"]
        if (
            not isinstance(filename, str)
            or not filename
            or len(filename) > 255
            or filename in {".", ".."}
            or "/" in filename
            or "\\" in filename
            or any(ord(c) < 32 for c in filename)
        ):
            raise ValueError("Invalid attachment filename.")
        return cls(kind, fmt, size, digest, filename)


def receive_attachment(sock: Any, dest: Path, attachment: Attachment) -> None:
    """Stream exactly the declared bytes into a temporary regular file."""
    digest = hashlib.sha256()
    remaining = attachment.size
    deadline = time.monotonic() + 600.0
    with dest.open("xb") as output:
        while remaining:
            remaining_time = deadline - time.monotonic()
            if remaining_time <= 0:
                raise ProtocolError("Attachment transfer deadline exceeded.")
            sock.settimeout(min(60.0, remaining_time))
            chunk = sock.recv(min(CHUNK_SIZE, remaining))
            if not chunk:
                raise ProtocolError("Attachment transfer interrupted.")
            output.write(chunk)
            digest.update(chunk)
            remaining -= len(chunk)
    if digest.hexdigest() != attachment.sha256:
        raise ValueError("Attachment SHA-256 checksum mismatch.")


def extract_shakemap(archive: Path, destination: Path) -> None:
    """Safely extract a flat ShakeMap ZIP, bounding expanded content."""
    destination.mkdir()
    try:
        with zipfile.ZipFile(archive) as source:
            entries = source.infolist()
            if not entries or len(entries) > MAX_ZIP_ENTRIES:
                raise ValueError("Invalid ShakeMap ZIP entry count.")
            unpacked = 0
            seen = set()
            file_names = [
                e.filename for e in entries if not e.is_dir()
            ]
            if not file_names:
                raise ValueError("ShakeMap ZIP is empty.")
            # Archives may have one enclosing event directory.
            parts = [PurePosixPath(name).parts for name in file_names]
            wrapped = (
                all(len(p) == 2 for p in parts)
                and len({p[0] for p in parts}) == 1
            )
            for entry in entries:
                name = entry.filename
                path = PurePosixPath(name)
                mode = (entry.external_attr >> 16) & 0xFFFF
                if (
                    not name
                    or name.startswith("/")
                    or "\\" in name
                    or path.is_absolute()
                    or ".." in path.parts
                    or any(":" in part for part in path.parts)
                    or entry.flag_bits & 1
                ):
                    raise ValueError("Unsafe ShakeMap ZIP entry.")
                if entry.is_dir():
                    if wrapped and name == parts[0][0] + "/":
                        continue
                    raise ValueError("Unexpected ShakeMap ZIP directory.")
                if (
                    len(path.parts) != (2 if wrapped else 1)
                    or (mode and stat.S_IFMT(mode) not in (0, stat.S_IFREG))
                ):
                    raise ValueError("Unsafe ShakeMap ZIP entry.")
                flat_name = path.parts[-1]
                if flat_name in seen:
                    raise ValueError("Duplicate ShakeMap ZIP filename.")
                seen.add(flat_name)
                if (
                    entry.file_size > 200 * max(1, entry.compress_size)
                    and entry.file_size > 1024 * 1024
                ):
                    raise ValueError("Excessive ZIP compression ratio.")
                unpacked += entry.file_size
                if unpacked > MAX_UNPACKED_BYTES:
                    raise ValueError("ShakeMap ZIP exceeds expanded size limit.")
                if entry.file_size > MAX_UNPACKED_BYTES:
                    raise ValueError("ShakeMap ZIP member too large.")
            for entry in entries:
                if entry.is_dir():
                    continue
                target = destination / PurePosixPath(entry.filename).name
                with source.open(entry) as reader, target.open("xb") as writer:
                    copied = 0
                    while True:
                        block = reader.read(CHUNK_SIZE)
                        if not block:
                            break
                        copied += len(block)
                        if copied > entry.file_size:
                            raise ValueError("ShakeMap ZIP member size mismatch.")
                        writer.write(block)
                    if copied != entry.file_size:
                        raise ValueError("ShakeMap ZIP member size mismatch.")
    except (zipfile.BadZipFile, RuntimeError, EOFError) as exc:
        raise ValueError("Invalid or corrupted ShakeMap ZIP.") from exc


def materialize_product(
    attachment: Attachment, source: Path, inputs: Path
) -> Path:
    """Create a validated product path inside the temporary inputs tree."""
    if attachment.kind == "shakemap":
        folder = inputs / "shakemap"
        extract_shakemap(source, folder)
        return folder

    folder = inputs / "imt"
    folder.mkdir()
    target = folder / "groundmotion.json"
    shutil.copyfile(source, target)
    return target


def validate_product(
    path: Path,
    attachment: Attachment,
    provider_config: Any,
    event_id: str | None,
) -> None:
    """Validate scientific input with the existing ShakeLab readers."""
    if attachment.kind == "shakemap":
        from shakelab.groundmotion.usgs import ShakeMapEvent

        config = provider_config.config
        if config.get("require_uncertainty") and not (
            path / "uncertainty.xml"
        ).is_file():
            raise ValueError("Required ShakeMap uncertainty.xml is missing.")
        product = ShakeMapEvent(
            folder=path,
            load_stations=False,
            load_uncertainty=config.get("load_uncertainty", True),
            require_info=config.get("require_info", False),
        )
        # Event metadata in info.json must not contradict the submitted id.
        product_id = product.event_id
        if (
            event_id is not None
            and product_id != "N/A"
            and product_id != event_id
        ):
            raise ValueError("ShakeMap event_id does not match the request.")
        return

    from shakelab.groundmotion.providers import GroundMotionProvider

    provider = GroundMotionProvider.imt(
        root_path=str(path.parent), config={"product_file": str(path)}
    )
    provider._load_event(None)
    data = json.loads(path.read_text(encoding="utf-8"))
    product_id = data.get("event_id")
    if (
        event_id is not None
        and product_id is not None
        and str(product_id) != event_id
    ):
        raise ValueError("IMT event_id does not match the request.")


def select_provider(model: Any, attachment: Attachment) -> Any:
    """Select exactly one compatible provider; never choose arbitrarily."""
    providers = [
        provider for provider in model.ground_motion.providers
        if provider.provider == attachment.kind
    ]
    if len(providers) != 1:
        raise ValueError(
            "Uploaded product requires exactly one matching "
            f"{attachment.kind} provider; found {len(providers)}."
        )
    return providers[0]
