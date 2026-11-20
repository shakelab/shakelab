"""Regression tests against WC1994 Table 2A and geometric identities."""

import pytest

from shakelab.groundmotion.scaling import (
    RuptureDimensions,
    WellsCoppersmith1994,
    fault_style_from_rake,
)


@pytest.mark.parametrize("style,area,length,width,sigma", [
    ("all", 1.97, 1.10, .91, .24),
    ("strike-slip", 1.98, 1.15, .86, .22),
    ("reverse", 1.89, 1.06, .85, .26),
    ("normal", 2.05, 1.12, .96, .22),
])
def test_forward_regressions_at_mw6(style, area, length, width, sigma):
    model = WellsCoppersmith1994(style)
    result = model.magnitude_to_area(6)
    assert result.median == pytest.approx(10 ** area)
    assert result.sigma == sigma
    assert result.sigma_space == "log10"
    assert result.at_epsilon(1) == pytest.approx(10 ** (area + sigma))
    assert model.magnitude_to_length(6).median == pytest.approx(10 ** length)
    assert model.magnitude_to_width(6).median == pytest.approx(10 ** width)


@pytest.mark.parametrize("style,expected,sigma", [
    ("all", 6.03, .24), ("strike-slip", 6.02, .23),
    ("reverse", 6.13, .25), ("normal", 5.97, .25),
])
def test_area_reverse_regression(style, expected, sigma):
    estimate = WellsCoppersmith1994(style).area_to_magnitude(100)
    assert estimate.median == pytest.approx(expected)
    assert estimate.sigma == sigma
    assert estimate.sigma_space == "Mw"
    assert estimate.at_epsilon(-1) == pytest.approx(expected - sigma)


@pytest.mark.parametrize("rake,style", [
    (None, "all"), (0, "strike-slip"), (90, "reverse"),
    (-90, "normal"), (180, "strike-slip"), (-180, "strike-slip"),
    (45, "strike-slip"), (-135, "strike-slip"),
])
def test_rake_classification(rake, style):
    assert fault_style_from_rake(rake) == style


def test_extrapolation_policy():
    with pytest.warns(UserWarning, match="outside"):
        result = WellsCoppersmith1994().magnitude_to_area(4.0)
    assert result.extrapolated
    with pytest.raises(ValueError, match="outside"):
        WellsCoppersmith1994(strict=True).magnitude_to_area(4.0)


@pytest.mark.parametrize("kwargs", [
    {"aspect_ratio": 2}, {"length": 20}, {"width": 10},
])
def test_rectangular_area_preserved(kwargs):
    dimensions = RuptureDimensions.from_area(200, **kwargs)
    assert dimensions.area == pytest.approx(200)
    assert dimensions.length == pytest.approx(20)
    assert dimensions.width == pytest.approx(10)


def test_width_cap_preserves_area():
    result = RuptureDimensions.from_area(
        200, aspect_ratio=2, max_width=5,
    )
    assert result.width == 5
    assert result.length == 40


@pytest.mark.parametrize("kwargs", [
    {}, {"width": 10, "length": 20}, {"aspect_ratio": 0},
    {"width": 10, "max_width": 5}, {"length": 20, "max_width": 5},
])
def test_invalid_dimensions(kwargs):
    with pytest.raises(ValueError):
        RuptureDimensions.from_area(200, **kwargs)


def test_invalid_scaling_inputs():
    with pytest.raises(ValueError):
        WellsCoppersmith1994("unknown")
    with pytest.raises(ValueError):
        WellsCoppersmith1994().area_to_magnitude(0)
    with pytest.raises(ValueError):
        WellsCoppersmith1994().magnitude_to_area(float("nan"))
    with pytest.raises(ValueError):
        fault_style_from_rake(181)


@pytest.mark.parametrize("style,length,width,sigma_l,sigma_w", [
    ("all", 5.87, 6.31, .26, .41),
    ("strike-slip", 5.82, 6.39, .24, .45),
    ("reverse", 5.98, 6.32, .26, .32),
    ("normal", 5.88, 6.15, .31, .31),
])
def test_reverse_length_and_width(style, length, width, sigma_l, sigma_w):
    model = WellsCoppersmith1994(style)
    result_l = model.length_to_magnitude(10)
    result_w = model.width_to_magnitude(10)
    assert result_l.median == pytest.approx(length)
    assert result_w.median == pytest.approx(width)
    assert result_l.sigma == sigma_l
    assert result_w.sigma == sigma_w
    assert not result_l.extrapolated
    assert not result_w.extrapolated


def test_legacy_interface_delegates_to_canonical_model():
    from shakelab.hazard.magnitude_scaling import wells_coppersmith_1994

    legacy = wells_coppersmith_1994()
    model = WellsCoppersmith1994("normal")
    expected = model.magnitude_to_area(6)
    assert legacy.magnitude_to_area(6, "n") == pytest.approx(
        (expected.median, expected.sigma),
    )
    expected = model.area_to_magnitude(100)
    assert legacy.area_to_magnitude(100, "normal") == pytest.approx(
        (expected.median, expected.sigma),
    )
