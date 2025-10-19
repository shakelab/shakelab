
# Exposure Model Format for ShakeLab

This document describes the JSON-based exposure model format used in the
ShakeLab framework for seismic impact and risk assessment. The format is
designed to be flexible, human-readable, and easily convertible to GeoJSON
for use in GIS software like QGIS.

## Structure

The model consists of two main sections:

- `metadata`: General information about the dataset.
- `assets`: A list of exposed elements, either individual structures or
  aggregated entities.

## Metadata Fields

| Field         | Type   | Description                                    |
|---------------|--------|------------------------------------------------|
| `name`        | string | Model name                                     |
| `description` | string | Description of the model                       |
| `region`      | string | Area covered                                   |
| `crs`         | string | Coordinate reference system (e.g., EPSG:4326)  |
| `source`      | string | Data source or provenance                      |
| `date`        | string | Creation or update date in YYYY-MM-DD format   |
| `version`     | string | Version of the model                           |

## Asset Fields

Each element in the `assets` list represents either an individual structure
or an aggregated asset.

| Field             | Type     | Description |
|------------------|----------|-------------|
| `id`             | string   | Unique identifier |
| `aggregated`     | boolean  | `false` for single elements, `true` for aggregates |
| `name`           | string   | Optional name or label |
| `usage`          | string   | Functional classification (e.g., residential, school) |
| `critical`       | boolean  | Indicates if the asset is critical infrastructure |
| `geometry_shape` | object   | GeoJSON geometry (Point or Polygon) |
| `reference_site_location` | object | Location used for ground motion calculation |

### `reference_site_location` Example

```json
"reference_site_location": {
  "longitude": 13.7663,
  "latitude": 45.6489,
  "elevation_m": 115.6
}
```

### Typologies

All asset typologies are listed under a single key `typologies`, which is
an array even for single buildings (`count = 1`).

Each typology describes a structural class, occupancy, and vulnerability
model.

| Field              | Type    | Description |
|-------------------|---------|-------------|
| `taxonomy`        | string  | Structural type classification (e.g., RC-LRS) |
| `count`           | integer | Number of buildings represented |
| `avg_stories`     | float   | Average number of above-ground stories |
| `avg_area_m2`     | float   | Average area per building (m²) |
| `total_value_eur` | float   | Total economic value in euros |
| `avg_occupants`   | float   | Average number of occupants |
| `fragility_set`   | string  | Identifier for the associated fragility curve |
| `retrofitted`     | boolean | Whether the buildings have been retrofitted |
| `building_type`   | string  | Descriptive structural category |
| `code_level`      | string  | Seismic design level (e.g., none, moderate, high) |

#### Example

```json
"typologies": [
  {
    "taxonomy": "RC-LRS",
    "count": 1,
    "avg_stories": 2,
    "avg_area_m2": 880,
    "total_value_eur": 1750000,
    "avg_occupants": 2.0,
    "fragility_set": "EMS98_B",
    "retrofitted": false,
    "building_type": "reinforced concrete",
    "code_level": "moderate"
  }
]
```

## Compatibility

This format is compatible with:

- ShakeLab computational modules
- GeoJSON conversion (for QGIS or other GIS software)
- Python-based parsers and validators

## Notes

- All coordinates use WGS84 (EPSG:4326).
- Asset geometries may be simplified to points or polygons.
- Asset typologies are always represented as an array.
- Typologies for single buildings must use `count = 1`.

## Example Assets

See `exposure_model_complete.json` for an example including both
single buildings and aggregated areas with multiple structural classes.

© 2025 ShakeLab Developers. Licensed under GPL v3.
