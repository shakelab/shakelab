"""Regression tests for explicit ground-motion product resolution."""

from pathlib import Path

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionModel,
    GroundMotionProviderConfig,
)
from shakelab.groundmotion.providers import ScenarioEvent


def _event():
    return ScenarioEvent.from_lonlat_depth(
        magnitude=5.0,
        longitude=13.0,
        latitude=46.0,
        depth_km=8.0,
    )


def _model(provider, root):
    return GroundMotionModel(
        type="ShakeLabGroundMotion",
        schema_version="1.0.0",
        providers=[
            GroundMotionProviderConfig(
                id="product",
                provider=provider,
                config={"root_path": str(root)},
            )
        ],
    )


@pytest.mark.parametrize(
    "backend,attribute,subpath",
    [
        ("shakemap", "product_directory", "inputs/shakemap"),
        ("imt", "product_file", "inputs/imt/groundmotion.json"),
    ],
)
def test_explicit_product_without_event_id(
    tmp_path, backend, attribute, subpath
):
    target = tmp_path / subpath
    target.parent.mkdir(parents=True, exist_ok=True)
    if backend == "shakemap":
        target.mkdir()
    else:
        target.write_text("{}", encoding="utf-8")

    runtime = _model(backend, tmp_path).runtime(
        _event(), product_paths={"product": target}
    )
    provider = runtime._providers["product"]
    assert getattr(provider, attribute) == target.resolve()
    assert runtime.event.event_id is None


@pytest.mark.parametrize("backend", ["shakemap", "imt"])
def test_legacy_path_resolution_unchanged(tmp_path, backend):
    runtime = _model(backend, tmp_path).runtime(_event())
    provider = runtime._providers["product"]

    if backend == "shakemap":
        assert provider.product_directory is None
        assert provider.event_directory("123") == tmp_path / "123"
    else:
        assert provider.product_file is None
        assert provider.event_file("123") == tmp_path / "123.json"


def test_unknown_provider_id_rejected(tmp_path):
    with pytest.raises(ValueError, match="Unknown product provider"):
        _model("shakemap", tmp_path).runtime(
            _event(), product_paths={"missing": tmp_path}
        )


def test_incompatible_provider_rejected(tmp_path):
    with pytest.raises(ValueError, match="only for shakemap and imt"):
        _model("gmpe", tmp_path).runtime(
            _event(), product_paths={"product": tmp_path}
        )
