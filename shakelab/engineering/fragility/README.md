# ShakeLab Fragility Model Format

This document describes the **JSON-based fragility model format** used by the
ShakeLab framework for seismic damage and impact analysis.

The goal of this format is to provide a representation of fragility models that
is:

- **clear and explicit**, separating conceptual roles,
- **easy to validate programmatically**,
- **flexible**, supporting continuous and discrete formulations,
- **extensible**, allowing multiple damage scales and authors,
- **stable**, so it can act as a reference interchange format.

The current specification corresponds to:

- **type**: `ShakeLabFragility`
- **schema_version**: `1.0.0`

This README is intended as a **didactic guide** to the format and its usage,
not as a changelog of individual schema updates.

---

## 1. Conceptual Overview

The fragility model format is organized in three conceptual layers:

1. **FragilityCollection** (root object)
2. **FragilityModel** (one model per taxonomy / author / IMT)
3. **Damage models** (one per damage state)

This separation is intentional:

- **FragilityCollection** acts as a database or catalogue.
- **FragilityModel** represents a published or calibrated fragility model
  applicable to a specific building typology and intensity measure.
- **Damage models** describe the probability of exceeding a *single* damage
  state as a function of the intensity measure.

The format supports both:

- **continuous (lognormal)** fragility formulations,
- **discrete (tabulated)** fragility formulations.

---

## 2. Root Object: `ShakeLabFragility`

The root object wraps the entire fragility database.

### Required Fields

| Field            | Type    | Description |
|------------------|---------|-------------|
| `type`           | string  | Must be exactly `"ShakeLabFragility"` |
| `schema_version` | string  | Must be `"1.0.0"` |
| `metadata`       | object  | Descriptive metadata |
| `models`         | array   | List of fragility models (non-empty) |

### Minimal Example

```json
{
  "type": "ShakeLabFragility",
  "schema_version": "1.0.0",
  "metadata": {
    "name": "Example fragility database",
    "date": "2026-02-07"
  },
  "models": []
}
```

---

## 3. Metadata

The `metadata` object stores general information about the fragility database.
Only `name` and `date` are mandatory.

### Common Metadata Fields

| Field | Type | Description |
|------|------|-------------|
| `name` | string | Database name |
| `description` | string | Free-text description |
| `source` | string | Data provenance |
| `date` | string | ISO date (`YYYY-MM-DD`) |
| `version` | string | Dataset-level version |
| `license` | string | Data license |

---

## 4. Fragility Models

Each element of the `models` array defines a **FragilityModel**.

A FragilityModel represents:

- one author or source,
- one construction taxonomy,
- one intensity measure type,
- one damage scale.

### Required Fields

| Field | Type | Description |
|------|------|-------------|
| `id` | string | Unique model identifier |
| `taxonomy` | string | Building taxonomy |
| `imt` | string | Intensity measure type (e.g. `PGA`, `SA(0.3)`) |
| `model_type` | string | Model family |
| `damage_scale` | object | Damage scale definition |
| `im_bounds` | object | Valid IM range |
| `parameters` / `tables` | object | Model parameters |

Additional fields (e.g. `author`, `year`, `notes`) are allowed and stored as
model metadata.

---

## 5. Damage Scale

The `damage_scale` object defines the damage states used by the model.

```json
"damage_scale": {
  "id": "D1D5",
  "levels": ["D1", "D2", "D3", "D4", "D5"]
}
```

- `levels` must be **ordered** from least to most severe.
- Damage state identifiers are arbitrary strings, but must be unique.

---

## 6. IM Bounds

`im_bounds` defines the valid range of the intensity measure.

```json
"im_bounds": {
  "min": 0.01,
  "max": 5.0
}
```

These bounds are informative and may be used for validation or clipping.

---

## 7. Continuous Lognormal Models

For `model_type = "lognormal_continuous"`, fragility is defined by a lognormal
distribution for each damage state.

```json
"model_type": "lognormal_continuous",
"parameters": {
  "D1": { "theta": 0.06, "beta": 0.60 },
  "D2": { "theta": 0.12, "beta": 0.60 },
  "D3": { "theta": 0.20, "beta": 0.60 }
}
```

Where:

- `theta` is the median capacity (in IM units),
- `beta` is the logarithmic standard deviation (natural log).

Each damage state is treated independently.

---

## 8. Discrete Models

For `model_type = "discrete"`, fragility is defined by tabulated PoE values.

```json
"model_type": "discrete",
"tables": {
  "D1": {
    "im": [0.05, 0.10, 0.20, 0.50],
    "poe": [0.01, 0.10, 0.40, 0.80],
    "log_im": true
  }
}
```

Where:

- `im` is a strictly increasing IM grid,
- `poe` are probabilities of exceedance,
- `log_im` specifies interpolation in log(IM) space.

Each damage state has its own table.

---

## 9. Probabilities of Exceedance vs State Probabilities

Fragility models define **probabilities of exceedance**:

> P(damage ≥ Di | IM)

These are cumulative and monotonic.

From them, ShakeLab can derive **state probabilities**:

> P(damage = Di | IM)

An additional state `GT_LAST` represents damage exceeding the highest defined
state.

---

## 10. Validation and Usage

This format is:

- validated by the ShakeLab fragility validator,
- designed to integrate with exposure and hazard components,
- suitable for deterministic and probabilistic damage calculations.

Fragility models are **independent of exposure geometry** and are linked to
assets via the `taxonomy` field.

---

## 11. Versioning Policy

- `schema_version = 1.0.0` identifies the structure described here.
- Backward-incompatible changes require a schema version bump.
- Dataset-level changes should use metadata `version`.

---

© 2026 ShakeLab Developers  
Licensed under the GNU General Public License v3
