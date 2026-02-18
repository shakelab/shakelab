# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
ShakeScenario: impact scenario calculator.

This module computes damage probabilities (or exceedance) over an exposure
model by combining:
- ground-motion evaluation (GroundMotionContext)
- explicit taxonomy->fragility mapping (TaxonomyTree)
- fragility model database (FragilityCollection)

Design choices
--------------
- The impact layer depends on the *ground-motion context* (already binding
  event + provider) rather than re-creating it internally. This keeps the
  impact module decoupled from the ground-motion factory and provider types.
- Damage is computed per typology and then aggregated to the asset level.
  Asset aggregation supports both:
    1) normalized damage *probabilities* (mixture over typologies)
    2) expected *counts* (count-weighted, not normalized)

Conventions
-----------
- GroundMotionContext returns (im_median_linear, sigma_ln).
- FragilityModel is queried either deterministically (poe_all) or with
  lognormal IM uncertainty (poe_lognormal_im_all).
- Output can be exceedance ("exceed") or mutually-exclusive states ("state").
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any, Dict, List, Literal, Mapping, Optional, Sequence, Tuple
import json, os

from shakelab.engineering.groundmotion import GroundMotionContext
from shakelab.engineering.exposure.exposure import Asset, ExposureModel
from shakelab.engineering.taxonomy.taxonomy_tree import TaxonomyTree
from shakelab.engineering.fragility.fragility import FragilityCollection


DamageOutput = Literal["state", "exceed"]
UncertaintyMode = Literal["lognormal", "median_only"]
MissingTaxonomyPolicy = Literal["raise", "skip"]
TypologyWeighting = Literal["count", "uniform"]


__all__ = [
    "ImpactConfig",
    "ImpactResult",
    "compute_impact_scenario",
    "damage_probabilities",
]


@dataclass(frozen=True)
class ImpactConfig:
    """
    Impact scenario configuration.

    Attributes
    ----------
    uncertainty_mode
        - "lognormal": convolve fragility PoE with IM lognormal uncertainty
          using (median IM, sigma_ln).
        - "median_only": deterministic computation using the median IM
          (sigma_ln ignored).
    output
        - "exceed": return exceedance probabilities P(DS >= level)
        - "state": return mutually exclusive probabilities including D0 and
          a tail state.
    typology_weighting
        How typologies are mixed at asset level:
        - "count": weights proportional to typology.count (default)
        - "uniform": each typology contributes equally
    normalize_asset_probabilities
        If True, the asset-level probability output is normalized so that the
        returned probabilities sum to ~1 (for output="state") or remain in
        [0, 1] (for output="exceed") as a convex mixture across typologies.
        If False, typology mixtures are summed without normalization.
        (Expected counts are always computed using typology.count.)
    missing_taxonomy
        Behavior when a typology.taxonomy is not found in the TaxonomyTree.
    no_damage_key
        Output key for the no-damage state (only for output="state").
    tail_key
        Output key for the tail state beyond the last level (only for
        output="state").
    """

    uncertainty_mode: UncertaintyMode = "lognormal"
    output: DamageOutput = "state"

    typology_weighting: TypologyWeighting = "count"
    normalize_asset_probabilities: bool = True

    missing_taxonomy: MissingTaxonomyPolicy = "raise"

    no_damage_key: str = "D0"
    tail_key: str = "GT_LAST"


@dataclass
class ImpactResult:
    """
    Scenario outputs.

    Attributes
    ----------
    damage_prob_by_asset
        Asset-level damage probabilities (keys depend on config.output).
    expected_count_by_asset
        Asset-level expected counts per damage state (only meaningful if the
        exposure typologies provide `count`). Keys follow the "state" output
        convention (including D0 and tail).
    """

    damage_prob_by_asset: Dict[str, Dict[str, float]]
    expected_count_by_asset: Dict[str, Dict[str, float]]


def compute_impact_scenario(
    gm_context: GroundMotionContext,
    exposure_model: ExposureModel,
    taxonomy_tree: TaxonomyTree,
    fragility_collection: FragilityCollection,
    config: Optional[ImpactConfig] = None,
) -> ImpactResult:
    """
    Compute a damage scenario for an exposure model.

    Parameters
    ----------
    gm_context
        Ground-motion context (event + provider).
    exposure_model
        ExposureModel instance.
    taxonomy_tree
        Explicit taxonomy -> fragility mapping.
    fragility_collection
        Fragility model database.
    config
        Optional ImpactConfig.

    Returns
    -------
    ImpactResult
        Asset-level probability outputs + expected counts.

    Notes
    -----
    - Probability output is a *mixture* over typologies, optionally normalized.
    - Expected counts are computed as sum(count_typ * P_state_typ), and are
      always returned in "state" convention (D0..Dn + tail).
    """
    cfg = config or ImpactConfig()

    damage_prob_by_asset: Dict[str, Dict[str, float]] = {}
    expected_count_by_asset: Dict[str, Dict[str, float]] = {}

    for asset in exposure_model.assets:
        asset_id = str(asset.id)

        is_aggregated = bool(getattr(asset, "aggregated", False))
        
        if not is_aggregated and len(asset.typologies) > 1:
            # Non-aggregated asset with multiple typologies is
            # ambiguous (composition vs epistemic uncertainty).
            raise ValueError(
                f"Asset {asset_id!r} is not aggregated but has multiple "
                "typologies. Ambiguous composition vs epistemic "
                "uncertainty. Set aggregated=True or provide a single "
                "typology."
            )

        lon, lat, elev = _asset_reference_lonlat(asset)

        # 1) Per-typology probabilities (asset mixture)
        typ_prob_list: List[Tuple[Dict[str, float], float]] = []

        # 2) Per-typology expected counts (always "state")
        asset_expected: Dict[str, float] = {}

        for typ in asset.typologies:
            tax = str(typ.taxonomy)

            if tax not in taxonomy_tree:
                if cfg.missing_taxonomy == "skip":
                    continue
                raise KeyError(
                    f"Exposure taxonomy not found in taxonomy tree: {tax!r}"
                )

            resolved = taxonomy_tree.resolve(tax, fragility_collection)

            # Damage probabilities at this site for this typology (cfg.output)
            probs = damage_probabilities(
                resolved=resolved,
                gm_context=gm_context,
                lon=lon,
                lat=lat,
                elevation_m=elev,
                output=cfg.output,
                mode=cfg.uncertainty_mode,
                no_damage_key=cfg.no_damage_key,
                tail_key=cfg.tail_key,
            )

            # Weight for asset probability mixture
            if cfg.typology_weighting == "uniform":
                w_typ = 1.0
            else:
                w_typ = float(getattr(typ, "count", 0.0) or 0.0)

            typ_prob_list.append((probs, w_typ))

            # Expected counts use "state" convention regardless of cfg.output
            probs_state = probs
            if cfg.output == "exceed":
                # Convert exceedance -> state using model levels from
                # the first resolved model (all are assumed compatible).
                levels = resolved[0][0].damage_scale.levels
            
                # Scientific guard: ensure exceedance keys match damage levels
                missing = [lv for lv in levels if lv not in probs]
                extra = [k for k in probs.keys() if k not in set(levels)]
                if missing or extra:
                    raise ValueError(
                        "Cannot convert exceedance to state probabilities: "
                        f"missing levels={missing}, extra keys={extra}."
                    )
            
                probs_state = _exceed_to_state(
                    probs,
                    levels,
                    no_damage_key=cfg.no_damage_key,
                    tail_key=cfg.tail_key,
                )

            cnt = float(getattr(typ, "count", 0.0) or 0.0)
            for k, p in probs_state.items():
                asset_expected[k] = asset_expected.get(k, 0.0) + cnt * float(p)

        damage_prob_by_asset[asset_id] = _mix_probabilities(
            typ_prob_list,
            normalize=cfg.normalize_asset_probabilities,
        )
        expected_count_by_asset[asset_id] = asset_expected

    return ImpactResult(
        damage_prob_by_asset=damage_prob_by_asset,
        expected_count_by_asset=expected_count_by_asset,
    )


def damage_probabilities(
    resolved: List[Tuple[Any, float]],
    gm_context: GroundMotionContext,
    lon: float,
    lat: float,
    elevation_m: float = 0.0,
    output: DamageOutput = "state",
    mode: UncertaintyMode = "lognormal",
    no_damage_key: str = "D0",
    tail_key: str = "GT_LAST",
) -> Dict[str, float]:
    """
    Compute weighted damage probabilities for one site and one taxonomy.

    Parameters
    ----------
    resolved
        Output of TaxonomyTree.resolve(): list of (FragilityModel, weight).
        Weights represent epistemic logic-tree weights.
    gm_context
        Ground-motion context returning:
            evaluate_at_site(imt, lon, lat, elevation_m) -> (im_median, sigma_ln)
    lon, lat, elevation_m
        Site coordinates.
    output
        "exceed" -> exceedance probabilities P(DS >= level).
        "state"  -> mutually exclusive probabilities including D0 and tail.
    mode
        "lognormal" uses (median IM, sigma_ln) and poe_lognormal_im_all.
        "median_only" uses median IM deterministically via poe_all.
    no_damage_key, tail_key
        Keys used only for output="state".

    Returns
    -------
    dict
        Damage probabilities.
    """
    if not resolved:
        raise ValueError("resolved must be a non-empty list.")

    model0 = resolved[0][0]
    levels_ref = list(model0.damage_scale.levels)
    imt_ref = str(model0.imt).strip()

    # Scientific guards: ensure the mapping is internally consistent.
    for m, _w in resolved[1:]:
        if list(m.damage_scale.levels) != levels_ref:
            raise ValueError(
                "Inconsistent damage scale in taxonomy mapping: "
                f"{levels_ref!r} vs {list(m.damage_scale.levels)!r}."
            )
        if str(m.imt).strip() != imt_ref:
            raise ValueError(
                "Inconsistent IMT in taxonomy mapping: "
                f"{imt_ref!r} vs {str(m.imt).strip()!r}. "
                "Do not mix fragility models with different IMT under "
                "the same taxonomy mapping."
            )

    combined: Dict[str, float] = {}
    wsum = 0.0

    for model, w in resolved:
        w_f = float(w)
        if not isfinite(w_f) or w_f <= 0.0:
            continue

        imt = str(model.imt).strip()
        im_med, sigma_ln = gm_context.evaluate_at_site(
            imt=imt,
            lon=float(lon),
            lat=float(lat),
            elevation_m=float(elevation_m),
        )

        # Scientific guards: avoid silently producing nonsense.
        if not isfinite(float(im_med)) or float(im_med) <= 0.0:
            raise ValueError(
                f"Non-positive IM returned by GM context for IMT={imt}: {im_med}"
            )
        if isfinite(float(sigma_ln)) and float(sigma_ln) < 0.0:
            raise ValueError(
                f"Negative sigma_ln returned by GM context for IMT={imt}: "
                f"{sigma_ln}"
            )

        # Compute exceedance probabilities for this model.
        if mode == "lognormal":
            sig = float(sigma_ln)
            if not isfinite(sig):
                poe = _poe_det(model, float(im_med), levels_ref)
            else:
                poe = model.poe_lognormal_im_all(float(im_med), sig)
        else:
            poe = _poe_det(model, float(im_med), levels_ref)

        # Enforce PoE monotonicity before mixing/conversion.
        _enforce_poe_monotone(poe, levels_ref)

        if output == "exceed":
            for lv in levels_ref:
                combined[lv] = combined.get(lv, 0.0) + w_f * _clip01(poe[lv])
        else:
            state = _exceed_to_state(
                poe,
                levels_ref,
                no_damage_key=no_damage_key,
                tail_key=tail_key,
            )
            for key, p in state.items():
                combined[key] = combined.get(key, 0.0) + w_f * _clip01(p)

        wsum += w_f

    if wsum > 0.0:
        for k in list(combined.keys()):
            combined[k] = float(combined[k]) / wsum
    else:
        raise ValueError(
            "No valid fragility weights found in resolved mapping "
            "(wsum <= 0)."
        )

    return combined


def save_impact_result(
    result: ImpactResult,
    output_path: str,
    config: ImpactConfig,
    metadata: Optional[Dict[str, Any]] = None,
    create_dirs: bool = True,
) -> None:
    """
    Save ImpactResult to a JSON file.

    Parameters
    ----------
    result
        ImpactResult instance.
    output_path
        Path to the output JSON file.
    config
        ImpactConfig used to compute the result.
    metadata
        Optional additional metadata.
    create_dirs
        If True, create parent directories if needed.
    """
    if create_dirs:
        directory = os.path.dirname(os.path.abspath(output_path))
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

    payload: Dict[str, Any] = {
        "type": "ShakeLabImpactResult",
        "schema_version": "1.0.0",
        "config": {
            "uncertainty_mode": config.uncertainty_mode,
            "output": config.output,
            "typology_weighting": config.typology_weighting,
            "normalize_asset_probabilities":
                config.normalize_asset_probabilities,
            "no_damage_key": config.no_damage_key,
            "tail_key": config.tail_key,
        },
        "damage_prob_by_asset": result.damage_prob_by_asset,
        "expected_count_by_asset":
            result.expected_count_by_asset,
    }

    if metadata:
        payload["metadata"] = metadata

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _asset_reference_lonlat(asset: Asset) -> Tuple[float, float, float]:
    """
    Extract reference lon/lat/elevation from an ExposureAsset.

    The exposure schema uses `reference_location` with keys:
    - longitude, latitude, elevation
    """
    loc = getattr(asset, "reference_location", None)
    if not isinstance(loc, Mapping):
        raise ValueError(
            f"Asset {asset.id!r} has no valid reference_location mapping."
        )

    lon = float(loc.get("longitude"))
    lat = float(loc.get("latitude"))
    elev = float(loc.get("elevation", 0.0))

    if not (isfinite(lon) and isfinite(lat) and isfinite(elev)):
        raise ValueError(
            f"Asset {asset.id!r} reference_location has non-finite values."
        )

    return lon, lat, elev


def _mix_probabilities(
    prob_list: Sequence[Tuple[Dict[str, float], float]],
    normalize: bool = True,
) -> Dict[str, float]:
    """
    Mix a list of probability dictionaries with associated weights.

    Scientific guard
    ---------------
    This function enforces that all non-zero-weight probability dictionaries
    share the same set of keys. Mixing distributions with different damage
    keys (e.g. D1D3 vs D1D5) is not scientifically interpretable without an
    explicit remapping policy, so we fail fast.
    """
    out: Dict[str, float] = {}
    wsum = 0.0

    ref_keys: Optional[set[str]] = None

    for probs, w in prob_list:
        w_f = float(w)
        if not isfinite(w_f) or w_f <= 0.0:
            continue

        keys = set(probs.keys())
        if ref_keys is None:
            ref_keys = keys
        elif keys != ref_keys:
            missing = sorted(ref_keys - keys)
            extra = sorted(keys - ref_keys)
            raise ValueError(
                "Inconsistent damage keys while mixing typologies. "
                f"Missing: {missing}. Extra: {extra}."
            )

        for k, p in probs.items():
            out[k] = out.get(k, 0.0) + w_f * float(p)

        wsum += w_f

    if normalize and wsum > 0.0:
        for k in list(out.keys()):
            out[k] = float(out[k]) / wsum

    return out


def _poe_det(model: Any, im: float, levels: Sequence[str]) -> Dict[str, float]:
    """
    Deterministic PoE for all levels; returns scalar floats.
    """
    poe_arr = model.poe_all(float(im))
    out: Dict[str, float] = {}
    for lv in levels:
        # poe_all returns 1D arrays (length 1 for scalar input).
        out[lv] = float(poe_arr[lv][0])
    return out


def _clip01(x: float) -> float:
    """Clip a probability to [0, 1]."""
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return float(x)


def _exceed_to_state(
    poe: Mapping[str, float],
    levels: Sequence[str],
    no_damage_key: str = "D0",
    tail_key: str = "GT_LAST",
) -> Dict[str, float]:
    """
    Convert exceedance P(DS>=k) to mutually exclusive P(DS=k).

    Assumes levels are ordered increasingly (D1, D2, ...).
    """
    out: Dict[str, float] = {}

    p_prev = 1.0
    for lv in levels:
        p_exc = _clip01(float(poe[lv]))
        out[lv] = _clip01(p_prev - p_exc)
        p_prev = p_exc

    out[tail_key] = _clip01(p_prev)
    out[no_damage_key] = _clip01(1.0 - sum(out.values()))

    # Renormalization
    ssum = sum(out.values())
    if isfinite(ssum) and ssum > 0.0:
        # Renormalize only if drifting noticeably (numerical safety)
        if abs(ssum - 1.0) > 1e-6:
            for k in list(out.keys()):
                out[k] = _clip01(out[k] / ssum)

    return out


def _enforce_poe_monotone(poe: Dict[str, float], levels: Sequence[str]) -> None:
    # Enforce PoE(L1) >= PoE(L2) >= ... >= PoE(Ln)
    prev = 1.0
    for lv in levels:
        val = _clip01(float(poe[lv]))
        if val > prev:
            val = prev
        poe[lv] = val
        prev = val


# ---------------------------------------------------------------------------
# Minimal smoke example
# ---------------------------------------------------------------------------

def main() -> None:
    exposure_file = "model/exposure_example.json"
    fragility_file = "model/fragility_example.json"
    taxonomy_tree_file = "model/taxonomy_tree_example.json"

    exposure_model = ExposureModel.from_json(exposure_file, validate=True)
    taxonomy_tree = TaxonomyTree.from_json(taxonomy_tree_file)
    fragility_collection = FragilityCollection.from_json(fragility_file)

    # Ground motion layer
    from shakelab.engineering.groundmotion import (
        GroundMotionProvider,
        ScenarioEvent,
    )
    from shakelab.libutils.geodeticN.primitives import WgsPoint

    event = ScenarioEvent(
        hypocentre=WgsPoint(longitude=13.0, latitude=46.0, elevation=-1e4),
        magnitude=5.5,
    )

    provider = GroundMotionProvider.gmpe(
        gmpe_name="BragatoSlejko2005",
        distance_approx="ellipsoid",
    )

    gm_context = GroundMotionContext(event=event, provider=provider)

    config = ImpactConfig(
        uncertainty_mode="lognormal",
        output="state",
        typology_weighting="count",
        normalize_asset_probabilities=True,
    )

    res = compute_impact_scenario(
        gm_context=gm_context,
        exposure_model=exposure_model,
        taxonomy_tree=taxonomy_tree,
        fragility_collection=fragility_collection,
        config=config,
    )

    print(res)

    save_impact_result(res, "impact_result.json", config)

if __name__ == "__main__":
    main()