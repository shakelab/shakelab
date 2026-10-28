"""Numerical integration of real fragility models and impact calculation.

The ground-motion evaluator is controlled; the fragility and impact
implementations are real. Parameters are synthetic, not a calibrated
building-class model.
"""

from __future__ import annotations

from math import log, sqrt
from types import SimpleNamespace

import pytest
from scipy.stats import norm

from shakelab.engineering.exposure.exposure import (
    Asset,
    ExposureModel,
    Typology,
)
from shakelab.engineering.fragility.fragility import (
    DamageScale,
    FragilityCollection,
    FragilityModel,
    IMBounds,
    LognormalDamageModel,
)
from shakelab.engineering.impact import ImpactConfig, compute_impact_scenario
from shakelab.engineering.taxonomy.taxonomy_tree import (
    FragilityRef,
    TaxonomyTree,
)


THETA = {"D1": 0.1, "D2": 0.2, "D3": 0.4}
BETA = 0.5
MEDIAN = 0.25
COUNT = 100


def _inputs(sigma):
    levels = tuple(THETA)
    model = FragilityModel(
        id="synthetic_rc",
        taxonomy="RC",
        imt="PGA",
        model_type="lognormal_continuous",
        damage_scale=DamageScale(id="test", levels=levels),
        im_bounds=IMBounds(min=0.0, max=2.0),
        damage_models={
            level: LognormalDamageModel(theta=theta, beta=BETA)
            for level, theta in THETA.items()
        },
    )
    collection = FragilityCollection(models={model.id: model})
    tree = TaxonomyTree(
        mappings={"RC": [FragilityRef(id=model.id, weight=1.0)]},
    )
    exposure = ExposureModel(
        type="ShakeLabExposure",
        schema_version="1.0.0",
        metadata={"name": "Synthetic integration fixture"},
        assets=[
            Asset(
                id="asset_1",
                aggregated=True,
                reference_location={"longitude": 13.0, "latitude": 46.0},
                typologies=[Typology(taxonomy="RC", count=COUNT)],
            ),
        ],
    )

    class Context:
        def evaluate_at_site(self, imt, lon, lat, elevation_m, **kwargs):
            assert imt == "PGA"
            return MEDIAN, sigma

    class Runtime:
        def resolve(self, asset_id):
            assert asset_id == "asset_1"
            return SimpleNamespace(
                provider_id="controlled",
                context=Context(),
                parameters={},
            )

    return Runtime(), exposure, tree, collection


def _expected_states(sigma):
    spread = sqrt(BETA ** 2 + sigma ** 2)
    exceed = [
        norm.cdf(log(MEDIAN / theta) / spread)
        for theta in THETA.values()
    ]
    return dict(zip(
        ("D0", "D1", "D2", "D3"),
        (1.0 - exceed[0],
         exceed[0] - exceed[1],
         exceed[1] - exceed[2],
         exceed[2]),
    ))


@pytest.mark.parametrize("sigma", [0.0, 0.35])
def test_real_fragility_impact_matches_analytic_solution(sigma):
    runtime, exposure, tree, collection = _inputs(sigma)
    result = compute_impact_scenario(
        runtime, exposure, tree, collection,
        config=ImpactConfig(uncertainty_mode="lognormal"),
    )
    asset = result.assets["asset_1"]
    expected = _expected_states(sigma)

    assert asset.n_units == COUNT
    assert asset.ground_motion_by_imt["PGA"].sigma_ln == sigma
    assert sum(asset.probabilities.values()) == pytest.approx(1.0)
    for state, probability in expected.items():
        assert asset.probabilities[state] == pytest.approx(
            probability, abs=1e-12,
        )
        assert asset.expected_counts[state] == pytest.approx(
            COUNT * probability, abs=1e-10,
        )
    assert sum(asset.expected_counts.values()) == pytest.approx(COUNT)


@pytest.mark.parametrize("sigma", [None, float("nan")])
def test_missing_uncertainty_rejected_in_lognormal_mode(sigma):
    runtime, exposure, tree, collection = _inputs(sigma)
    with pytest.raises(ValueError, match="sigma_ln"):
        compute_impact_scenario(
            runtime, exposure, tree, collection,
            config=ImpactConfig(uncertainty_mode="lognormal"),
        )


def test_median_only_ignores_missing_uncertainty():
    runtime, exposure, tree, collection = _inputs(None)
    result = compute_impact_scenario(
        runtime, exposure, tree, collection,
        config=ImpactConfig(uncertainty_mode="median_only"),
    )
    expected = _expected_states(0.0)
    actual = result.assets["asset_1"].probabilities
    for state, probability in expected.items():
        assert actual[state] == pytest.approx(probability, abs=1e-12)
