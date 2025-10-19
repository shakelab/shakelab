from exposure import (
    convert_exposure_file_to_geojson,
    convert_geojson_to_exposure_file
)

# === 1. JSON → GeoJSON (explode typologies, save metadata)
convert_exposure_file_to_geojson(
    input_path="exposure_model_fvg_2025.json",
    output_path="converted/exposure_model_fvg_2025.geojson",
    explode_typologies=True,
    save_metadata=True
)

# === 2. JSON → GeoJSON (no explosion, no metadata)
convert_exposure_file_to_geojson(
    input_path="exposure_model_fvg_2025.json",
    output_path="converted/exposure_model_flat.geojson",
    explode_typologies=False,
    save_metadata=False
)

# === 3. GeoJSON + meta → JSON
convert_geojson_to_exposure_file(
    input_paths="converted/exposure_model_fvg_2025.geojson",
    output_path="converted/exposure_model_fvg_2025_restored.json",
    metadata_path="exposure_model_fvg_2025.meta.json"
)

# === 4. GeoJSON → JSON (no metadata)
convert_geojson_to_exposure_file(
    input_paths="converted/exposure_model_flat.geojson",
    output_path="converted/restored_without_meta.json"
)

# === 5. Multi-file GeoJSON (points + polygons) → JSON
convert_exposure_file_to_geojson(
    input_path="exposure_model_fvg_2025.json",
    output_path="converted/only_points.geojson",
    feature="point"
)

convert_exposure_file_to_geojson(
    input_path="exposure_model_fvg_2025.json",
    output_path="converted/only_polygons.geojson",
    feature="polygon"
)

# Poi unisce entrambi in un file JSON
convert_geojson_to_exposure_file(
    input_paths=[
        "converted/only_points.geojson",
        "converted/only_polygons.geojson"
    ],
    output_path="converted/merged_from_types.json"
)
