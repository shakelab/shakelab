"""Regression tests for ground-motion uncertainty in impact calculations."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from shakelab.engineering.impact import damage_probabilities


@pytest.fixture
def resolved():
    """Provide a minimal fragility model with observable method calls."""
    model = SimpleNamespace(
        imt="PGA",
        damage_scale=SimpleNamespace(levels=["D1"]),
        poe_all=Mock(return_value={"D1": [0.2]}),
        poe_lognormal_im_all=Mock(return_value={"D1": 0.4}),
    )
    return [(model, 1.0)]


@pytest.mark.parametrize("sigma", [None, float("nan"), float("inf"),
                                    -float("inf"), -0.1])
def test_lognormal_rejects_missing_or_invalid_sigma(resolved, sigma):
    with pytest.raises(ValueError, match="sigma_ln"):
        damage_probabilities(
            resolved, im_median=0.2, sigma_ln=sigma,
            mode="lognormal", output="exceed",
        )
    resolved[0][0].poe_lognormal_im_all.assert_not_called()
    resolved[0][0].poe_all.assert_not_called()


@pytest.mark.parametrize("sigma", [0.0, 0.5])
def test_lognormal_uses_probabilistic_method(resolved, sigma):
    result = damage_probabilities(
        resolved, im_median=0.2, sigma_ln=sigma,
        mode="lognormal", output="exceed",
    )
    assert result["D1"] == pytest.approx(0.4)
    resolved[0][0].poe_lognormal_im_all.assert_called_once_with(0.2, sigma)
    resolved[0][0].poe_all.assert_not_called()


@pytest.mark.parametrize("sigma", [None, float("nan"), float("inf"), -0.1])
def test_median_only_ignores_sigma(resolved, sigma):
    result = damage_probabilities(
        resolved, im_median=0.2, sigma_ln=sigma,
        mode="median_only", output="exceed",
    )
    assert result["D1"] == pytest.approx(0.2)
    resolved[0][0].poe_lognormal_im_all.assert_not_called()
