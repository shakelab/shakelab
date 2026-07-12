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
Layer builders.

This module contains helpers that transform ShakeScenario outputs into
application-oriented geospatial layers.

The functions in this module operate on plain JSON-like dictionaries. This
keeps them independent from the engineering calculation classes and allows
them to be reused by both the ShakeScenario server and the WebUI backend.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


__all__ = [
    "build_impact_geojson",
]


DAMAGE_LEVELS = (
    "D0",
    "D1",
    "D2",
    "D3",
    "D4",
    "D5",
    "GT_LAST",
)


def build_impact_geojson(
    impact_assets: dict[str, Any],
    geometry_geojson: dict[str, Any],
    missing_geometry: str = "skip",
    geometry_id_keys: tuple[str, ...] = (
        "id",
        "asset_id"
    ),
) -> dict[str, Any]:
    """
    Build a GeoJSON impact layer by joining impact results and geometries.

    Parameters
    ----------
    impact_assets
        Dictionary loaded from ``impact_assets.json``.
    geometry_geojson
        GeoJSON FeatureCollection containing asset geometries. Each feature
        must expose an asset identifier in one of the
        ``geometry_id_keys`` properties, or in the top-level feature ``id``.
    missing_geometry
        Policy for assets without matching geometry:
        - ``"skip"``: omit assets without geometry.
        - ``"raise"``: raise a ValueError.
    geometry_id_keys
        Property names searched in geometry features to find the asset id.

    Returns
    -------
    dict
        GeoJSON FeatureCollection. Feature geometries are copied from
        ``geometry_geojson`` and feature properties are populated with the
        corresponding impact result fields. Common scalar fields are also
        exposed directly for WebUI styling, filtering, and popups.

    Raises
    ------
    ValueError
        If inputs are malformed, if duplicate geometries are found, or if
        ``missing_geometry`` is invalid.
    """
    if missing_geometry not in ("skip", "raise"):
        raise ValueError(
            "missing_geometry must be either 'skip' or 'raise'."
        )

    assets = _impact_assets_list(impact_assets)
    geom_index = _geometry_index(
        geometry_geojson,
        geometry_id_keys=geometry_id_keys,
    )

    features: list[dict[str, Any]] = []
    missing: list[str] = []

    for asset in assets:
        asset_id = _asset_result_id(asset)
        geom_feature = _find_geometry(asset_id, geom_index)

        if geom_feature is None:
            missing.append(asset_id)
            continue

        geometry = geom_feature.get("geometry")
        if not isinstance(geometry, dict):
            raise ValueError(
                f"Geometry feature for asset {asset_id!r} has no geometry."
            )

        feature = {
            "type": "Feature",
            "id": asset_id,
            "geometry": deepcopy(geometry),
            "properties": _impact_properties(asset, geom_feature),
        }
        features.append(feature)

    if missing and missing_geometry == "raise":
        raise ValueError(
            "Missing geometry for impact assets: "
            + ", ".join(sorted(missing))
        )

    if not features and assets:
        raise ValueError(
            "No matching geometries found for any impact asset. "
            "Check that asset identifiers match the geometry layer."
        )

    return {
        "type": "FeatureCollection",
        "features": features,
    }


def _impact_assets_list(
    impact_assets: dict[str, Any],
) -> list[dict[str, Any]]:
    """Return and validate the asset list from an impact result object."""
    if not isinstance(impact_assets, dict):
        raise ValueError("impact_assets must be a JSON object.")

    assets = impact_assets.get("assets")
    if not isinstance(assets, list):
        raise ValueError("impact_assets must contain an 'assets' list.")

    out: list[dict[str, Any]] = []
    for i, asset in enumerate(assets):
        if not isinstance(asset, dict):
            raise ValueError(f"impact asset at index {i} is not an object.")
        out.append(asset)

    return out


def _geometry_index(
    geometry_geojson: dict[str, Any],
    geometry_id_keys: tuple[str, ...],
) -> dict[str, dict[str, Any]]:
    """Build an asset-id index from a GeoJSON FeatureCollection."""
    if not isinstance(geometry_geojson, dict):
        raise ValueError("geometry_geojson must be a JSON object.")

    if geometry_geojson.get("type") != "FeatureCollection":
        raise ValueError("geometry_geojson must be a FeatureCollection.")

    features = geometry_geojson.get("features")
    if not isinstance(features, list):
        raise ValueError("geometry_geojson must contain a 'features' list.")

    index: dict[str, dict[str, Any]] = {}

    for i, feature in enumerate(features):
        if not isinstance(feature, dict):
            raise ValueError(f"geometry feature at index {i} is not an object.")

        if feature.get("type") != "Feature":
            raise ValueError(f"geometry feature at index {i} is not a Feature.")

        asset_id = _geometry_feature_id(
            feature,
            index=i,
            geometry_id_keys=geometry_id_keys,
        )

        if asset_id in index:
            raise ValueError(f"Duplicate geometry for asset id: {asset_id!r}")

        index[asset_id] = feature

    return index


def _geometry_feature_id(
    feature: dict[str, Any],
    index: int,
    geometry_id_keys: tuple[str, ...],
) -> str:
    """Extract the asset id from a geometry feature."""
    props = feature.get("properties")
    if props is not None and not isinstance(props, dict):
        raise ValueError(
            f"geometry feature at index {index} has invalid properties."
        )

    if isinstance(props, dict):
        for key in geometry_id_keys:
            if key in props and props[key] is not None:
                return _normalize_id(props[key])

    if feature.get("id") is not None:
        return _normalize_id(feature["id"])

    raise ValueError(
        "geometry feature at index "
        f"{index} has no usable geometry identifier."
    )


def _asset_result_id(asset: dict[str, Any]) -> str:
    """Extract the asset id from a serialized impact asset."""
    asset_id = asset.get("id")
    if asset_id is None:
        raise ValueError("impact asset without id.")
    return str(asset_id)


def _find_geometry(
    asset_id: str,
    geom_index: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    """Find geometry using the exact normalized asset id."""
    key = _normalize_id(asset_id)
    return geom_index.get(key)


def _normalize_id(value: Any) -> str:
    """Normalize ids loaded from JSON/GeoJSON properties."""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _impact_properties(
    asset: dict[str, Any],
    geom_feature: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Build GeoJSON properties for one impact asset.

    The output keeps the impact result nested enough to avoid losing the
    complete calculation result, while exposing common scalar fields directly
    for styling and filtering in the WebUI.
    """
    asset_id = _asset_result_id(asset)

    properties: dict[str, Any] = {
        "id": asset_id,
        "asset_id": asset_id,
    }

    geom_props = None
    if isinstance(geom_feature, dict):
        geom_props = geom_feature.get("properties")

    if isinstance(geom_props, dict):
        properties["geometry_properties"] = deepcopy(geom_props)
        _copy_first(
            geom_props,
            properties,
            ("name", "asset_name"),
            "name",
        )

    if "n_units" in asset:
        properties["n_units"] = asset["n_units"]

    if "reference_location" in asset:
        properties["reference_location"] = deepcopy(
            asset["reference_location"]
        )

    ground_motion = asset.get("ground_motion")
    if isinstance(ground_motion, dict):
        properties["ground_motion"] = deepcopy(ground_motion)
        _add_ground_motion_scalars(properties, ground_motion)

    damage = asset.get("damage")
    if isinstance(damage, dict):
        properties["damage"] = deepcopy(damage)
        _add_damage_scalars(properties, damage)

    if "typologies" in asset:
        properties["typologies"] = deepcopy(asset["typologies"])

    return properties


def _add_ground_motion_scalars(
    properties: dict[str, Any],
    ground_motion: dict[str, Any],
) -> None:
    """Expose common ground-motion scalar values."""
    pga = ground_motion.get("PGA")

    if isinstance(pga, dict):
        pga_median = _as_number(pga.get("median"))
        pga_sigma = _as_number(pga.get("sigma_ln"))
    else:
        pga_median = _as_number(pga)
        pga_sigma = None

    if pga_median is not None:
        properties["PGA"] = pga_median
        properties["pga"] = pga_median
        properties["pga_median"] = pga_median

    if pga_sigma is not None:
        properties["pga_sigma_ln"] = pga_sigma


def _add_damage_scalars(
    properties: dict[str, Any],
    damage: dict[str, Any],
) -> None:
    """Expose common damage scalar values."""
    expected = damage.get("expected_counts")
    probabilities = damage.get("probabilities")

    if isinstance(expected, dict):
        for level in DAMAGE_LEVELS:
            value = _as_number(expected.get(level))
            if value is not None:
                properties[level] = value
                properties[f"expected_{level}"] = value

        d4 = _as_number(expected.get("D4")) or 0.0
        d5 = _as_number(expected.get("D5")) or 0.0
        gt_last = _as_number(expected.get("GT_LAST")) or 0.0

        properties["damage_d4_d5"] = d4 + d5
        properties["damage_gt_last"] = gt_last
        properties["damage_total"] = d4 + d5 + gt_last
        properties["damage"] = d4 + d5 + gt_last

    if isinstance(probabilities, dict):
        for level in DAMAGE_LEVELS:
            value = _as_number(probabilities.get(level))
            if value is not None:
                properties[f"prob_{level}"] = value

        prob_d4 = _as_number(probabilities.get("D4")) or 0.0
        prob_d5 = _as_number(probabilities.get("D5")) or 0.0
        prob_gt_last = _as_number(probabilities.get("GT_LAST")) or 0.0

        properties["prob_damage_d4_d5"] = prob_d4 + prob_d5
        properties["prob_damage_total"] = prob_d4 + prob_d5 + prob_gt_last


def _copy_first(
    source: dict[str, Any],
    target: dict[str, Any],
    source_keys: tuple[str, ...],
    target_key: str,
) -> None:
    """Copy the first available source key to a target key."""
    for key in source_keys:
        if key in source and source[key] is not None:
            target[target_key] = source[key]
            return


def _as_number(value: Any) -> int | float | None:
    """Convert a JSON value to int or float when possible."""
    if value is None or value == "":
        return None

    try:
        number = float(value)
    except (TypeError, ValueError):
        return None

    if number.is_integer():
        return int(number)

    return number
