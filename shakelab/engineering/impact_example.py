# ---------------------------------------------------------------------------
# Impact example with corrected damage-state semantics
# ---------------------------------------------------------------------------

from shakelab.engineering.impact import (
    ImpactConfig,
    compute_impact_scenario,
    save_impact_result,
    save_impact_summary,
)

from shakelab.engineering.exposure.exposure import ExposureModel
from shakelab.engineering.fragility.fragility import FragilityCollection
from shakelab.engineering.groundmotion.groundmotion import GroundMotionModel
from shakelab.engineering.taxonomy.taxonomy_tree import TaxonomyTree

from shakelab.gmmodel.groundmotion import ScenarioEvent
from shakelab.libutils.geodeticN.primitives import WgsPoint


# ---------------------------------------------------------------------------
# Input models
# ---------------------------------------------------------------------------

exposure = ExposureModel.from_json(
    "model/exposure_example.json",
    validate=True,
)

taxonomy_tree = TaxonomyTree.from_json(
    "model/taxonomy_tree_example.json"
)

fragility = FragilityCollection.from_json(
    "model/fragility_example.json"
)

ground_motion_model = GroundMotionModel.from_json(
    "model/groundmotion_example.json",
    exposure_model=exposure,
    validate=True,
)


# ---------------------------------------------------------------------------
# Scenario event
# ---------------------------------------------------------------------------

event = ScenarioEvent(
    hypocentre=WgsPoint(
        longitude=13.0,
        latitude=46.0,
        elevation=-10000.0,
    ),
    magnitude=5.5,
)

ground_motion = ground_motion_model.runtime(event)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def print_result(title, result):
    print(f"\n{title}")

    for asset_id, asset_result in result.assets.items():
        print(f"\nAsset: {asset_id}")

        for imt, gm_value in asset_result.ground_motion_by_imt.items():
            print(
                f"  {imt}: "
                f"median={gm_value.median:.6g}, "
                f"sigma_ln={gm_value.sigma_ln:.6g}, "
                f"provider={gm_value.provider_id}"
            )

        print("  probabilities:")
        for key, value in asset_result.probabilities.items():
            print(f"    {key}: {value:.8f}")

        print("  expected_counts:")
        for key, value in asset_result.expected_counts.items():
            print(f"    {key}: {value:.8f}")


# ---------------------------------------------------------------------------
# State probabilities: D0..Dn
# ---------------------------------------------------------------------------

state_config = ImpactConfig(
    uncertainty_mode="lognormal",
    output="state",
    typology_weighting="count",
    normalize_asset_probabilities=True,
    include_typology_breakdown=True,
)

state_result = compute_impact_scenario(
    ground_motion=ground_motion,
    exposure_model=exposure,
    taxonomy_tree=taxonomy_tree,
    fragility_collection=fragility,
    config=state_config,
)

print_result(
    "MUTUALLY EXCLUSIVE DAMAGE STATES",
    state_result,
)


# ---------------------------------------------------------------------------
# Exceedance probabilities: P(DS >= Dk)
# ---------------------------------------------------------------------------

exceed_config = ImpactConfig(
    uncertainty_mode="lognormal",
    output="exceed",
    typology_weighting="count",
    normalize_asset_probabilities=True,
    include_typology_breakdown=True,
)

exceed_result = compute_impact_scenario(
    ground_motion=ground_motion,
    exposure_model=exposure,
    taxonomy_tree=taxonomy_tree,
    fragility_collection=fragility,
    config=exceed_config,
)

print_result(
    "DAMAGE-STATE EXCEEDANCE PROBABILITIES",
    exceed_result,
)


# ---------------------------------------------------------------------------
# Consistency checks
# ---------------------------------------------------------------------------

for asset_id in state_result.assets:
    states = state_result.assets[asset_id].probabilities
    exceed = exceed_result.assets[asset_id].probabilities

    state_sum = sum(states.values())
    if abs(state_sum - 1.0) > 1e-10:
        raise RuntimeError(
            f"State probabilities for {asset_id} sum to {state_sum}."
        )

    # Reconstruct exceedance from mutually exclusive states.
    damage_levels = [
        key for key in states
        if key != state_config.no_damage_key
    ]

    for i, level in enumerate(damage_levels):
        reconstructed = sum(
            states[key]
            for key in damage_levels[i:]
        )

        if abs(reconstructed - exceed[level]) > 1e-10:
            raise RuntimeError(
                f"Inconsistent state/exceedance conversion for "
                f"{asset_id}, {level}: "
                f"{reconstructed} != {exceed[level]}"
            )

print("\nState/exceedance consistency checks passed.")


# ---------------------------------------------------------------------------
# Save state-based canonical smoke-test products
# ---------------------------------------------------------------------------

save_impact_result(
    result=state_result,
    output_path="output/impact_assets.json",
    config=state_config,
)

save_impact_summary(
    result=state_result,
    output_path="output/impact_summary.json",
    config=state_config,
)
