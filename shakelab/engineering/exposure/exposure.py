# ****************************************************************************
# Copyright (C) 2019-2025, ShakeLab Developers.
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
Conversion functions between ShakeLab JSON exposure models and GeoJSON
format compatible with QGIS.
"""
import json
from pathlib import Path
from typing import Union, List, Optional, Dict


class Subtype:
    def __init__(
        self,
        taxonomy: str,
        count: int,
        avg_stories: float,
        avg_area_m2: float,
        total_value_eur: float,
        avg_occupants: float,
        fragility_set: str
    ):
        self.taxonomy = taxonomy
        self.count = count
        self.avg_stories = avg_stories
        self.avg_area_m2 = avg_area_m2
        self.total_value_eur = total_value_eur
        self.avg_occupants = avg_occupants
        self.fragility_set = fragility_set

    def to_dict(self) -> dict:
        return self.__dict__


class Asset:
    def __init__(
        self,
        asset_id: str,
        aggregated: bool,
        name: Optional[str] = None,
        taxonomy: Optional[str] = None,
        usage: Optional[str] = None,
        year_built: Optional[int] = None,
        stories: Optional[int] = None,
        area_m2: Optional[float] = None,
        value_eur: Optional[float] = None,
        occupants: Optional[Dict[str, float]] = None,
        vulnerability_class: Optional[str] = None,
        fragility_set: Optional[str] = None,
        critical: Optional[bool] = False,
        geometry_shape: Optional[Dict] = None,
        reference_site_location: Optional[Dict] = None,
        subtypes: Optional[List[Subtype]] = None
    ):
        self.id = asset_id
        self.aggregated = aggregated
        self.name = name
        self.taxonomy = taxonomy
        self.usage = usage
        self.year_built = year_built
        self.stories = stories
        self.area_m2 = area_m2
        self.value_eur = value_eur
        self.occupants = occupants
        self.vulnerability_class = vulnerability_class
        self.fragility_set = fragility_set
        self.critical = critical
        self.geometry_shape = geometry_shape
        self.reference_site_location = reference_site_location
        self.subtypes = subtypes or []

    def to_dict(self) -> dict:
        data = self.__dict__.copy()
        data["subtypes"] = [s.to_dict() for s in self.subtypes]
        return data


class Exposure:
    def __init__(self, metadata: Dict, assets: List[Asset]):
        self.metadata = metadata
        self.assets = assets

    def to_dict(self) -> dict:
        return {
            "metadata": self.metadata,
            "assets": [a.to_dict() for a in self.assets]
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'Exposure':
        assets = []
        for a in data.get("assets", []):
            subtypes = [Subtype(**s) for s in a.get("subtypes", [])]
            asset = Asset(
                asset_id=a["id"],
                aggregated=a["aggregated"],
                name=a.get("name"),
                taxonomy=a.get("taxonomy"),
                usage=a.get("usage"),
                year_built=a.get("year_built"),
                stories=a.get("stories"),
                area_m2=a.get("area_m2"),
                value_eur=a.get("value_eur"),
                occupants=a.get("occupants"),
                vulnerability_class=a.get("vulnerability_class"),
                fragility_set=a.get("fragility_set"),
                critical=a.get("critical", False),
                geometry_shape=a.get("geometry_shape"),
                reference_site_location=a.get("reference_site_location"),
                subtypes=subtypes
            )
            assets.append(asset)
        return cls(metadata=data.get("metadata", {}), assets=assets)


def convert_exposure_file_to_geojson(
    input_path: Union[str, Path],
    output_path: Union[str, Path] = None,
    feature: str = None,
    explode_typologies: bool = True,
    save_metadata: bool = True
    ) -> None:
    """
    Convert a ShakeLab exposure model to a GeoJSON file for QGIS use.

    Parameters
    ----------
    input_path : str or Path
        Path to the input ShakeLab JSON exposure model.
    output_path : str or Path, optional
        Path to save the GeoJSON. Defaults to .geojson extension.
    feature : str, optional
        If 'point' or 'polygon', filter by geometry type.
    explode_typologies : bool, default True
        If True, creates a feature per typology.
    save_metadata : bool, default True
        If True, saves metadata in a separate .meta.json file.
    """
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_suffix(".geojson")
    else:
        output_path = Path(output_path)

    with open(input_path, "r", encoding="utf-8") as f:
        exposure = json.load(f)

    features = []

    for asset in exposure.get("assets", []):
        if "geometry_shape" in asset:
            geometry = asset["geometry_shape"]
        elif "reference_site_location" in asset:
            rsl = asset["reference_site_location"]
            geometry = {
                "type": "Point",
                "coordinates": [
                    rsl["longitude"],
                    rsl["latitude"]
                ]
            }
        else:
            continue

        geom_type = geometry.get("type", "").lower()
        if feature:
            if feature == "point" and geom_type != "point":
                continue
            elif feature == "polygon" and geom_type != "polygon":
                continue

        rsl = asset.get("reference_site_location", {})
        rsl_flat = {
            f"reference_{k}": v for k, v in rsl.items()
            if k in ("longitude", "latitude", "elevation_m")
        }

        if explode_typologies:
            for i, typ in enumerate(asset.get("typologies", [])):
                props = {
                    "parent_id": asset["id"],
                    "typology_index": i,
                    "name": asset.get("name"),
                    "aggregated": asset.get("aggregated", False),
                    "usage": asset.get("usage"),
                    "critical": asset.get("critical", False),
                    **typ,
                    **rsl_flat
                }
                features.append({
                    "type": "Feature",
                    "geometry": geometry,
                    "properties": props
                })
        else:
            props = {
                k: v for k, v in asset.items()
                if k not in ("geometry_shape", "reference_site_location")
            }
            props.update(rsl_flat)
            features.append({
                "type": "Feature",
                "geometry": geometry,
                "properties": props
            })

    geojson = {
        "type": "FeatureCollection",
        "features": features
    }

    # Salva GeoJSON con indentazione normale
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(geojson, f, indent=2)

    # Salva file .meta.json con i metadati
    if save_metadata and "metadata" in exposure:
        meta_path = output_path.with_suffix(".meta.json")
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(exposure["metadata"], f, indent=2)


def convert_geojson_to_exposure_file(
    input_paths: Union[str, Path, List[Union[str, Path]]],
    output_path: Union[str, Path],
    metadata_path: Union[str, Path, None] = None
    ) -> None:
    """
    Convert one or more GeoJSON files into a ShakeLab-compatible
    JSON exposure model. Features with the same parent_id are grouped
    under a single asset with typologies.

    Parameters
    ----------
    input_paths : str, Path or list of str/Path
        One or more GeoJSON input files.
    output_path : str or Path
        Path where the resulting ShakeLab JSON file will be saved.
    metadata_path : str, Path or None
        Optional path to a .meta.json file with metadata. If None or missing,
        generic metadata fields are used.
    """
    if isinstance(input_paths, (str, Path)):
        input_paths = [input_paths]

    combined_assets = []
    grouped_typologies = {}

    for path in input_paths:
        path = Path(path)
        with open(path, "r", encoding="utf-8") as f:
            geojson = json.load(f)

        if geojson.get("type") != "FeatureCollection":
            raise ValueError(f"Invalid GeoJSON format in {path.name}")

        for feature in geojson.get("features", []):
            geometry = feature.get("geometry", {})
            props = feature.get("properties", {}).copy()

            geom_type = geometry.get("type", "").lower()
            if geom_type not in ("point", "polygon"):
                continue

            rsl = {}
            for k in ("longitude", "latitude", "elevation_m"):
                val = props.pop(f"reference_{k}", None)
                if val is not None:
                    rsl[k] = val

            if "parent_id" in props:
                pid = props.pop("parent_id")
                typology = {
                    k: v for k, v in props.items()
                    if k not in (
                        "typology_index", "name", "geometry_shape",
                        "aggregated", "usage", "critical"
                    )
                }

                if pid not in grouped_typologies:
                    grouped_typologies[pid] = {
                        "id": pid,
                        "aggregated": True,
                        "name": props.get("name"),
                        "usage": props.get("usage"),
                        "critical": props.get("critical", False),
                        "geometry_shape": geometry,
                        "reference_site_location": rsl,
                        "typologies": []
                    }

                grouped_typologies[pid]["typologies"].append(typology)
            else:
                asset = {
                    "id": props.get("id"),
                    "aggregated": props.get("aggregated", False),
                    "name": props.get("name"),
                    "usage": props.get("usage"),
                    "critical": props.get("critical", False),
                    "geometry_shape": geometry,
                    "reference_site_location": rsl,
                    "typologies": props.get("typologies", [])
                }
                if isinstance(asset["typologies"], dict):
                    asset["typologies"] = [asset["typologies"]]
                combined_assets.append(asset)

    combined_assets.extend(grouped_typologies.values())

    # Gestione metadati opzionale
    if metadata_path:
        try:
            with open(metadata_path, "r", encoding="utf-8") as f:
                metadata = json.load(f)
        except FileNotFoundError:
            metadata = None
    else:
        metadata = None

    metadata_default = {
        "name": "MergedExposureModel",
        "description": "Imported from GeoJSON",
        "region": "Unknown",
        "crs": "EPSG:4326",
        "source": "GeoJSON",
        "date": "YYYY-MM-DD",
        "version": "2.1"
    }

    model = {
        "metadata": metadata or metadata_default,
        "assets": combined_assets
    }

    output_path = Path(output_path)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(model, f, indent=2)

