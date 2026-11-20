"""Tests of automatic finite-rupture construction."""

import pytest

from shakelab.groundmotion.rupture import build_rectangular_rupture
from shakelab.libutils.geodeticN.primitives import WgsPoint


def build(**kwargs):
    options = dict(strike=120, dip=90)
    options.update(kwargs)
    depth = options.pop("depth", 8)
    return build_rectangular_rupture(
        WgsPoint(13, 46, -1000 * depth), **options,
    )


@pytest.mark.parametrize("shape", [
    {"length": 20, "width": 10},
    {"area": 200, "length": 20},
    {"area": 200, "width": 10},
    {"area": 200, "aspect_ratio": 2},
    {"length": 20, "aspect_ratio": 2},
    {"width": 10, "aspect_ratio": 2},
])
def test_complete_dimensions(shape):
    result = build(**shape)
    assert result.rupture.length == pytest.approx(20)
    assert result.rupture.width == pytest.approx(10)
    assert result.rupture.area == pytest.approx(200)
    assert result.rupture.ztor == pytest.approx(3)
    assert result.rupture.down_fraction == pytest.approx(.5)
    assert result.rupture.along_fraction == .5


@pytest.mark.parametrize("model,style", [
    ("WellsCoppersmith1994", "auto"),
    ("ThingbaijamEtAl2017", "auto"),
    ("Leonard2014_Interplate", "auto"),
    ("Leonard2014_StableContinental", "auto"),
    ("HanksBakun2002", "auto"),
])
def test_all_scaling_models(model, style):
    result = build(
        magnitude=6.5, rake=90, aspect_ratio=2,
        scaling_model=model, scaling_style=style,
    )
    assert result.rupture.area == pytest.approx(
        result.scaling_estimate.median,
    )
    assert result.scaling_model == model
    assert result.parameter_origins["area"] == "scaling"
    assert result.scaling_style == (
        "all" if model == "HanksBakun2002" else "reverse"
    )


def test_explicit_dimensions_bypass_scaling():
    result = build(length=20, width=10, scaling_model="unknown")
    assert result.scaling_estimate is None
    assert result.parameter_origins["length"] == "explicit"
    assert result.parameter_origins["area"] == "derived"


def test_width_cap_preserves_area():
    result = build(
        depth=5, area=200, aspect_ratio=2,
        upper_seismogenic_depth=2, lower_seismogenic_depth=8,
    )
    assert result.rupture.width == pytest.approx(6)
    assert result.rupture.area == pytest.approx(200)
    assert result.rupture.ztor == pytest.approx(2)
    assert result.parameter_origins["width"] == "constrained"
    assert result.parameter_origins["aspect_ratio"] == "constrained"


@pytest.mark.parametrize("depth,top", [(1, 0), (14, 5)])
def test_shift_central_nucleation(depth, top):
    result = build(depth=depth, length=20, width=10,
                   lower_seismogenic_depth=15)
    assert result.rupture.ztor == pytest.approx(top)
    assert result.parameter_origins["ztor"] == "constrained"


def test_min_width_for_explicit_ztor():
    result = build(area=100, aspect_ratio=4, ztor=0)
    assert result.rupture.width == pytest.approx(8)
    assert result.rupture.area == pytest.approx(100)
    assert result.rupture.down_fraction == pytest.approx(1)


def test_ztor_and_fraction_determine_width():
    result = build(area=200, ztor=3, down_fraction=.5)
    assert result.rupture.width == pytest.approx(10)
    assert result.rupture.length == pytest.approx(20)


def test_explicit_fraction_limits_width():
    result = build(depth=2, area=200, aspect_ratio=2,
                   down_fraction=.5)
    assert result.rupture.width == pytest.approx(4)
    assert result.rupture.ztor == pytest.approx(0)
    assert result.rupture.area == pytest.approx(200)


def test_horizontal_rupture():
    result = build(dip=0, area=200, aspect_ratio=2,
                   lower_seismogenic_depth=9)
    assert result.rupture.ztor == pytest.approx(8)
    assert result.rupture.down_fraction == .5


@pytest.mark.parametrize("options", [
    {"length": 20, "width": 10, "area": 100},
    {"length": 20, "width": 10, "aspect_ratio": 1},
    {"length": 20, "width": 10, "lower_seismogenic_depth": 9},
    {"length": 20, "width": 10, "depth": 20,
     "lower_seismogenic_depth": 15},
    {"length": 20, "width": 10, "ztor": 4, "down_fraction": .5},
    {"length": 20, "width": 10, "down_fraction": 1},
    {"length": 20, "width": 10, "dip": 0, "ztor": 2},
    {"area": 200},
    {"aspect_ratio": 2, "magnitude": 6},
    {"area": 200, "aspect_ratio": 0},
    {"area": 200, "aspect_ratio": 2, "ztor": 9},
    {"area": 200, "aspect_ratio": 2, "down_fraction": 2},
])
def test_incompatible_constraints(options):
    with pytest.raises(ValueError):
        build(**options)


def test_unknown_style_without_rake():
    with pytest.raises(ValueError, match="Unknown fault style"):
        build(magnitude=6.5, scaling_model="ThingbaijamEtAl2017",
              aspect_ratio=2)


def test_scaling_residual_is_used():
    result = build(magnitude=6.5, scaling_model="WellsCoppersmith1994",
                   aspect_ratio=2, scaling_epsilon=1)
    assert result.rupture.area == pytest.approx(
        result.scaling_estimate.at_epsilon(1),
    )


def test_explicit_area_takes_precedence():
    result = build(area=200, aspect_ratio=2, magnitude=6.5,
                   scaling_model="WellsCoppersmith1994")
    assert result.scaling_estimate is None
    assert result.parameter_origins["area"] == "explicit"


@pytest.mark.parametrize("depth,fraction", [(0, 0), (15, 1)])
def test_hypocentre_at_layer_boundary(depth, fraction):
    result = build(depth=depth, area=200, aspect_ratio=2,
                   lower_seismogenic_depth=15)
    assert result.rupture.down_fraction == pytest.approx(fraction)
    assert result.rupture.area == pytest.approx(200)


def test_explicit_length_cannot_be_changed_by_width_cap():
    with pytest.raises(ValueError, match="width"):
        build(area=200, length=20, lower_seismogenic_depth=9)


def test_inclined_rupture_respects_depth_constraints():
    result = build(depth=6, dip=30, area=400, aspect_ratio=1,
                   upper_seismogenic_depth=3, lower_seismogenic_depth=9)
    assert result.rupture.width == pytest.approx(12)
    assert result.rupture.ztor == pytest.approx(3)
    assert result.rupture.down_fraction == pytest.approx(.5)
    assert result.rupture.area == pytest.approx(400)
