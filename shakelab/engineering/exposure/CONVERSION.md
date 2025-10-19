
# ShakeLab GeoJSON Conversion Guide

This document describes how to convert exposure models from the ShakeLab JSON
format to GeoJSON and back. It covers all supported options, features, and
provides concrete examples.

## Purpose

The conversion to GeoJSON allows:
- Direct visualization in GIS software (e.g., QGIS)
- Manual editing of typologies and asset attributes
- Exporting maps or thematic styling

The conversion back to ShakeLab JSON ensures:
- Round-trip consistency
- Compatibility with ShakeLab's computational modules

---

## Conversion: JSON → GeoJSON

Function:
```python
convert_exposure_file_to_geojson(...)
```

### Options

| Parameter            | Description                                       |
|----------------------|---------------------------------------------------|
| `input_path`         | Path to ShakeLab JSON exposure model              |
| `output_path`        | Output GeoJSON path (default: same name + .geojson) |
| `explode_typologies` | If True, creates one feature per typology         |
| `feature`            | Filter by geometry type: 'point', 'polygon', or None |
| `save_metadata`      | If True, saves a `.meta.json` with model metadata |

### Examples

Convert full model with exploded typologies and metadata:

```python
convert_exposure_file_to_geojson(
    "exposure_model_fvg_2025.json",
    "exposure_model_fvg_2025.geojson",
    explode_typologies=True,
    save_metadata=True
)
```

Convert without exploding typologies:

```python
convert_exposure_file_to_geojson(
    "exposure_model_fvg_2025.json",
    "exposure_model_flat.geojson",
    explode_typologies=False,
    save_metadata=False
)
```

Convert only point geometries:

```python
convert_exposure_file_to_geojson(
    "exposure_model_fvg_2025.json",
    "only_points.geojson",
    feature="point"
)
```

---

## Conversion: GeoJSON → JSON

Function:
```python
convert_geojson_to_exposure_file(...)
```

### Options

| Parameter       | Description                                      |
|-----------------|--------------------------------------------------|
| `input_paths`   | One or more GeoJSON files                        |
| `output_path`   | Output path for the ShakeLab JSON file           |
| `metadata_path` | Optional `.meta.json` file containing metadata   |

### Examples

Convert back using the metadata file:

```python
convert_geojson_to_exposure_file(
    "exposure_model_fvg_2025.geojson",
    "restored_model.json",
    metadata_path="exposure_model_fvg_2025.meta.json"
)
```

Convert back without a metadata file (default metadata will be used):

```python
convert_geojson_to_exposure_file(
    "exposure_model_flat.geojson",
    "restored_without_meta.json"
)
```

Merge multiple GeoJSON files:

```python
convert_geojson_to_exposure_file(
    ["only_points.geojson", "only_polygons.geojson"],
    "merged_model.json"
)
```

---

## File Extensions

| File                     | Description                             |
|--------------------------|-----------------------------------------|
| `.json`                  | ShakeLab exposure model                 |
| `.geojson`               | GeoJSON file for GIS                    |
| `.meta.json`             | Optional metadata companion file        |

---

## Notes

- `typologies` are always represented as a list in JSON and GeoJSON.
- If `explode_typologies=True`, each feature corresponds to one typology.
- The metadata file is required only to restore original metadata on import.
- Geometry is encoded in standard GeoJSON format.

---

© 2025 ShakeLab Developers. Licensed under GPL v3.
