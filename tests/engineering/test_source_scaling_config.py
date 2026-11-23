"""Configuration-only tests for opt-in GMPE source scaling."""

import pytest

from shakelab.engineering.gmdata.associator import (
    GroundMotionProviderConfig,
)
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    validate_source_scaling,
)


def config(scaling, geometry="finite"):
    return GroundMotionProviderConfig(
        id="gmpe",
        provider="gmpe",
        config={
            "gmpe_name": "ChiouYoungs2014",
            "source_geometry": geometry,
            "source_scaling": scaling,
        },
    )


@pytest.mark.parametrize("model", [
    "WellsCoppersmith1994",
    "ThingbaijamEtAl2017",
    "Leonard2014_Interplate",
    "Leonard2014_StableContinental",
    "HanksBakun2002",
])
def test_all_models_accept_auto_without_rake_at_configuration(model):
    normalized = validate_source_scaling({"model": model})
    assert normalized["style"] == "auto"
    assert normalized["strict"] is False
    config({"model": model}).validate()


@pytest.mark.parametrize("model,style", [
    ("ThingbaijamEtAl2017", "interface"),
    ("ThingbaijamEtAl2017", "reverse"),
    ("WellsCoppersmith1994", "all"),
    ("Leonard2014_StableContinental", "normal"),
    ("HanksBakun2002", "all"),
])
def test_supported_styles(model, style):
    config({"model": model, "style": style}).validate()


@pytest.mark.parametrize("scaling,exception", [
    ({}, ValueError),
    ({"model": "missing"}, ValueError),
    ({"model": "HanksBakun2002", "style": "reverse"}, ValueError),
    ({"model": "ThingbaijamEtAl2017", "style": "all"}, ValueError),
    ({"model": "WellsCoppersmith1994", "style": "interface"},
     ValueError),
    ({"model": "WellsCoppersmith1994", "unknown": 1}, ValueError),
    ({"model": "WellsCoppersmith1994", "strict": 1}, TypeError),
    ({"model": "WellsCoppersmith1994", "epsilon": True}, ValueError),
    ({"model": "WellsCoppersmith1994", "epsilon": float("nan")},
     ValueError),
    ({"model": "HanksBakun2002", "strict": True}, ValueError),
    ({"model": "Leonard2014_Interplate", "strict": True}, ValueError),
    ({"model": "HanksBakun2002", "epsilon": 1}, ValueError),
])
def test_invalid_scaling_fails_before_runtime(scaling, exception):
    with pytest.raises(exception):
        config(scaling).validate()


def test_strict_supported_for_calibrated_models():
    config({"model": "ThingbaijamEtAl2017", "strict": True}).validate()
    config({"model": "WellsCoppersmith1994", "strict": True}).validate()


@pytest.mark.parametrize("geometry", ["native", "point", None])
def test_scaling_requires_finite_provider_geometry(geometry):
    with pytest.raises(ValueError, match="source_geometry='finite'"):
        config({"model": "WellsCoppersmith1994"},
               geometry=geometry).validate()


def test_existing_provider_without_scaling_is_unchanged():
    config(None, geometry="native").validate()
    assert validate_source_scaling(None) is None


def test_constructor_rejects_incompatible_geometry():
    with pytest.raises(ValueError, match="source_geometry='finite'"):
        GroundMotionProvider.from_id("gmpe", config={
            "gmpe_name": "ChiouYoungs2014",
            "source_geometry": "point",
            "source_scaling": {"model": "WellsCoppersmith1994"},
        })
