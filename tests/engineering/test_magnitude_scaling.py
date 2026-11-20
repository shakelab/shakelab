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


@pytest.mark.parametrize("style,rows", [
    ("reverse", [( -2.693, .614, .083), (-1.669, .435, .087),
                 (-4.362, 1.049, .121)]),
    ("normal", [(-1.722, .485, .128), (-.829, .323, .128),
                (-2.551, .808, .181)]),
    ("strike-slip", [(-2.943, .681, .151), (-.543, .261, .105),
                     (-3.486, .942, .184)]),
    ("interface", [(-2.412, .583, .107), (-.880, .366, .099),
                   (-3.292, .949, .150)]),
])
def test_tmg2017_table1(style, rows):
    from shakelab.groundmotion.scaling import create_scaling_relation

    model = create_scaling_relation("ThingbaijamEtAl2017", style=style)
    for quantity, (a, b, sigma) in zip(
        ("length", "width", "area"), rows,
    ):
        estimate = model.estimate(7.0, quantity)
        assert estimate.median == pytest.approx(10 ** (a + b * 7))
        assert estimate.sigma == sigma
        reverse = model.estimate_magnitude(estimate.median, quantity)
        assert reverse.median == pytest.approx(7)
        assert reverse.sigma == pytest.approx(sigma / b)


def test_scaling_factory():
    from shakelab.groundmotion.scaling import (
        available_scaling_relations, create_scaling_relation,
    )

    assert available_scaling_relations() == (
        "HanksBakun2002", "Leonard2014_Interplate",
        "Leonard2014_StableContinental", "ThingbaijamEtAl2017",
        "WellsCoppersmith1994",
    )
    assert isinstance(
        create_scaling_relation("WellsCoppersmith1994"),
        WellsCoppersmith1994,
    )
    with pytest.raises(ValueError, match="Unknown scaling relation"):
        create_scaling_relation("unknown")
    with pytest.raises(ValueError, match="Unknown fault style"):
        create_scaling_relation("ThingbaijamEtAl2017", style="all")
    with pytest.raises(ValueError, match="outside"):
        create_scaling_relation(
            "ThingbaijamEtAl2017", style="normal", strict=True,
        ).magnitude_to_area(5)


@pytest.mark.parametrize("name,style,constant", [
    ("Leonard2014_Interplate", "all", 3.995),
    ("Leonard2014_Interplate", "strike-slip", 3.99),
    ("Leonard2014_Interplate", "normal", 4.00),
    ("Leonard2014_Interplate", "reverse", 4.00),
    ("Leonard2014_StableContinental", "all", 4.185),
    ("Leonard2014_StableContinental", "strike-slip", 4.18),
    ("Leonard2014_StableContinental", "normal", 4.19),
    ("Leonard2014_StableContinental", "reverse", 4.19),
])
def test_leonard_area(name, style, constant):
    from shakelab.groundmotion.scaling import create_scaling_relation

    model = create_scaling_relation(name, style=style)
    estimate = model.magnitude_to_area(6)
    assert estimate.median == pytest.approx(10 ** (6 - constant))
    assert model.area_to_magnitude(100).median == pytest.approx(2 + constant)
    assert model.area_to_magnitude(estimate.median).median == pytest.approx(6)
    assert estimate.sigma is None
    assert estimate.extrapolated is None


@pytest.mark.parametrize("area,expected", [
    (100, 5.98), (1000, 7.07),
])
def test_hanks_bakun_branches(area, expected):
    from shakelab.groundmotion.scaling import create_scaling_relation

    model = create_scaling_relation("HanksBakun2002")
    assert model.area_to_magnitude(area).median == pytest.approx(expected)
    assert model.magnitude_to_area(expected).median == pytest.approx(area)


def test_hanks_bakun_crossover():
    from math import log10
    from shakelab.groundmotion.scaling import HanksBakun2002

    model = HanksBakun2002()
    assert model.area_to_magnitude(537).median == pytest.approx(
        3.98 + log10(537),
    )
    assert model.magnitude_to_area(
        model.crossover_magnitude,
    ).median == pytest.approx(537)
    assert model.area_to_magnitude(538).median == pytest.approx(
        3.07 + (4 / 3) * log10(538),
    )


@pytest.mark.parametrize("name", [
    "HanksBakun2002", "Leonard2014_Interplate", "Leonard2014_StableContinental",
])
def test_area_only_model_contract(name):
    from shakelab.groundmotion.scaling import create_scaling_relation

    model = create_scaling_relation(name)
    estimate = model.estimate(6, "area")
    assert estimate.at_epsilon(0) == estimate.median
    with pytest.raises(ValueError, match="dispersion"):
        estimate.at_epsilon(1)
    with pytest.raises(ValueError, match="Unsupported"):
        model.estimate(6, "length")
    with pytest.raises(ValueError, match="Unsupported"):
        model.estimate_magnitude(10, "width")
    with pytest.raises(ValueError, match="calibration"):
        create_scaling_relation(name, strict=True)
    with pytest.raises(ValueError, match="positive"):
        model.area_to_magnitude(0)
    with pytest.raises(ValueError, match="finite"):
        model.magnitude_to_area(float("nan"))
