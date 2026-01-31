# ShakeLab Exposure Model Format

This document describes the **JSON-based exposure model format** used by the
ShakeLab framework for seismic impact and scenario analysis.

The goal of this format is to provide a representation of exposed buildings
that is:

- **clear and explicit**, even for non-developers,
- **easy to validate programmatically**,
- **flexible**, supporting both minimal and richly described datasets,
- **compatible with GIS tools** (e.g. QGIS),
- **stable**, so it can act as a reference interchange format.

The current specification corresponds to:

- **type**: `ShakeLabExposure`
- **schema_version**: `1.0.0`

This README is intended as a **didactic guide** to the format and its usage,
not as a changelog of individual schema updates.

---

## 1. Conceptual Overview

The exposure model is organized in three hierarchical levels:

1. **Exposure** (root object)
2. **Assets** (georeferenced exposure units)
3. **Typologies** (building classes inside each asset)

This separation is intentional:

- **Assets** describe *where* exposed elements are located.
- **Typologies** describe *what* those elements are (construction class,
  usage, vulnerability-relevant attributes).

This structure allows representing:

- individual buildings,
- aggregated building stock (e.g. census zones, grid cells),
- mixed-use and mixed-construction areas,
- partially known exposure inventories.

The model supports **progressive enrichment**: a minimal exposure can be
valid, and additional information can be added incrementally as it becomes
available.

---

## 2. Root Object: `ShakeLabExposure`

The root object wraps the entire exposure model.

### Required Fields

| Field            | Type    | Description |
|------------------|---------|-------------|
| `type`           | string  | Must be exactly `"ShakeLabExposure"` |
| `schema_version` | string  | Must be `"1.0.0"` |
| `metadata`       | object  | Descriptive metadata |
| `assets`         | array   | List of exposure assets (non-empty) |

### Minimal Example

```json
{
  "type": "ShakeLabExposure",
  "schema_version": "1.0.0",
  "metadata": {
    "name": "Example exposure",
    "date": "2025-01-15"
  },
  "assets": []
}
```

---

## 3. Metadata

The `metadata` object stores general information about the dataset.
Only `name` and `date` are mandatory; all other fields are optional but
strongly recommended for clarity and reproducibility.

### Common Metadata Fields

| Field | Type | Description |
|------|------|-------------|
| `name` | string | Model name |
| `description` | string | Free-text description |
| `region` | string | Geographic coverage |
| `source` | string | Data provenance |
| `date` | string | ISO date (`YYYY-MM-DD`) |
| `version` | string | Dataset-level version |
| `crs` | string | Coordinate reference system (e.g. `EPSG:4326`) |
| `currency` | string | ISO 4217 code (e.g. `EUR`) |
| `units` | object | Units for numeric quantities |
| `occupants_unit` | string | Recommended: `persons_per_building` |
| `license` | string | Data license |

### Example

```json
"metadata": {
  "name": "NEI25 Exposure Model",
  "description": "Exposure model for North-East Italy",
  "region": "Friuli Venezia Giulia / Veneto",
  "source": "ISTAT 2021, OGS elaborations",
  "date": "2025-01-15",
  "crs": "EPSG:4326",
  "currency": "EUR",
  "units": {
    "aggregation_area": "m2",
    "elevation": "m"
  },
  "occupants_unit": "persons_per_building"
}
```

---

## 4. Assets

Each element of the `assets` array represents an **exposure unit**.
An asset can be either:

- a **single building**, or
- an **aggregated asset** (polygon representing many buildings).

### Asset Fields

| Field | Type | Description |
|------|------|-------------|
| `id` | string | Unique identifier |
| `aggregated` | boolean | `false` = single building, `true` = aggregate |
| `reference_location` | object | Representative point (mandatory) |
| `typologies` | array | List of building typologies (non-empty) |
| `name` | string \| null | Optional human-readable label |
| `aggregation_area` | number \| null | Area of aggregate |
| `critical` | boolean \| null | Optional critical-infrastructure flag |
| `geometry` | object \| null | GeoJSON-like geometry |
| `reference_geology` | object \| null | Optional site/geology attributes |

Only the first four fields are strictly required. All other fields may be
omitted or set to `null`.

---

## 5. Geometry

Geometry follows a GeoJSON-like structure and is **optional**.

When provided:

- If `aggregated = false` → geometry must be `Point`
- If `aggregated = true` → geometry must be `Polygon`

Coordinates are always `[longitude, latitude]`.

### Point Example (Single Building)

```json
"geometry": {
  "type": "Point",
  "coordinates": [13.7663, 45.6489]
}
```

### Polygon Example (Aggregated Asset)

```json
"geometry": {
  "type": "Polygon",
  "coordinates": [
    [
      [13.4521, 45.9479],
      [13.4538, 45.9479],
      [13.4538, 45.9490],
      [13.4521, 45.9490],
      [13.4521, 45.9479]
    ]
  ]
}
```

If geometry is not provided, spatial analyses rely exclusively on
`reference_location`.

---

## 6. Reference Location

`reference_location` defines a representative point for the asset and is
**always required**.

It is used for:

- site-condition assignment,
- ground-motion extraction,
- labeling and plotting,
- fallback spatial reference when geometry is absent.

```json
"reference_location": {
  "longitude": 13.7663,
  "latitude": 45.6489,
  "elevation": 115.6
}
```

For single buildings, this usually coincides with the Point geometry.
For polygons, it typically represents a centroid or meaningful location.

---

## 7. Typologies

Typologies describe **homogeneous groups of buildings** within an asset.

All assets must contain **at least one typology**.
Single buildings are represented as typologies with `count = 1`.

### Typology Fields

| Field | Type | Description |
|------|------|-------------|
| `taxonomy` | string | Construction class identifier |
| `count` | integer | Number of buildings (≥ 1) |
| `usage` | string \| null | Functional use |
| `building_type` | string \| null | Descriptive structural category |
| `code_level` | string \| null | Seismic code compliance level |
| `occupants` | object \| null | Day/night occupancy |
| `period` | object \| null | Construction period |
| `replacement_cost` | number \| null | Replacement cost |
| `stories` | integer \| null | Number of stories |
| `damage_state` | string \| null | Optional damage label |

Only `taxonomy` and `count` are mandatory.

---

## 8. Minimal Valid Exposure Example

The following illustrates the **minimum valid exposure model**:

- one asset,
- no geometry,
- no geology,
- one typology with taxonomy and count.

```json
{
  "type": "ShakeLabExposure",
  "schema_version": "1.0.0",
  "metadata": {
    "name": "Minimal example",
    "date": "2025-01-15"
  },
  "assets": [
    {
      "id": "A001",
      "aggregated": false,
      "reference_location": {
        "longitude": 13.5,
        "latitude": 45.9
      },
      "typologies": [
        {
          "taxonomy": "RC",
          "count": 1
        }
      ]
    }
  ]
}
```

This model is fully valid and can be progressively enriched.

---

## 9. Validation and Usage

This format is:

- validated by the ShakeLab exposure validator,
- suitable for conversion to GeoJSON for GIS visualization,
- designed for scenario-based damage and impact calculations.

The exposure model is **neutral with respect to hazard, vulnerability and
loss models**. Typology attributes are intended to be mapped to fragility
or vulnerability models by downstream components.

---

## 10. Versioning Policy

- `schema_version = 1.0.0` identifies the structure described here.
- Backward-incompatible changes require a schema version bump.
- Metadata `version` may be used for dataset-level revisions.

---

© 2026 ShakeLab Developers  
Licensed under the GNU General Public License v3
