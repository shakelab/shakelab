# Empirical earthquake magnitude scaling

`shakelab.groundmotion.scaling` provides independent empirical
magnitude–rupture-size relations and helpers for consistent
rectangular dimensions. It can be used without a GMPE provider,
rupture surface, scenario model or server.

**Units and calibration:** magnitude is **Mw**, length and width are
in km, and rupture area is in km². For geometry construction and
`source_geometry` modes, see [RUPTURE.md](RUPTURE.md).

## 1. Available models

Use the exact identifiers returned by
`available_scaling_relations()` with
`create_scaling_relation(name, style=..., strict=...)`.

| Model identifier | Forward prediction | Reverse prediction | Styles | Residual dispersion / calibration checks |
| --- | --- | --- | --- | --- |
| `WellsCoppersmith1994` | Mw → area, length, width | Area, length, width → Mw | `all`, `strike-slip`, `reverse`, `normal` | Available; published calibration bounds |
| `ThingbaijamEtAl2017` | Mw → area, length, width | Area, length, width → Mw | `strike-slip`, `reverse`, `normal`, `interface` | Available; style-dependent bounds |
| `Leonard2014_Interplate` | Mw → area | Area → Mw | `all`, `strike-slip`, `reverse`, `normal` | Not supplied in this implementation |
| `Leonard2014_StableContinental` | Mw → area | Area → Mw | `all`, `strike-slip`, `reverse`, `normal` | Not supplied in this implementation |
| `HanksBakun2002` | Mw → area | Area → Mw | `all` only | Not supplied in this implementation |

The tabulated availability is about **implemented functionality**, not
a statement that all models share the same tectonic calibration or
statistical interpretation. `Leonard2014_Interplate` means active
crustal faulting in this implementation, **not** a specialised
subduction-interface law.

### Model references and implementation details

**Wells & Coppersmith (1994).** *Bulletin of the Seismological Society
of America*, 84(4), 974–1002, Table 2A.
https://doi.org/10.1785/BSSA0840040974.
Independent forward and reverse regressions are implemented for
area, length and width, as defined by the selected coefficient rows. Forward coefficients predict log10(size); reverse
coefficients predict Mw. These are separate empirical fits, **not
algebraic inverses**.

The forward calibration interval depends on style and size measure.
Reverse length and area diagnostics use the published dimension ranges.
The width bounds printed in Table 2A repeat the subsurface length
ranges; therefore reverse-width diagnostics instead use the predicted
Mw limits, avoiding a false claim of a published width interval.

**Thingbaijam, Mai & Goda (2017).** *Bulletin of the Seismological
Society of America*, 107, 2225–2246, Table 1.
https://doi.org/10.1785/0120170017.
Relations are mechanism/tectonic-setting dependent. No pooled `all`
fit is provided. Its `interface` category is a tectonic setting and
must be selected explicitly, not inferred from rake. The area
coefficients in this implementation are derived from length and
width regressions. Reverse values invert the fitted lines, using
residual scale `sigma / abs(slope)`; they should **not** be
interpreted as independent conditional Mw regressions.

**Leonard (2014).** https://doi.org/10.1785/0120140087.
`Leonard2014_Interplate` and
`Leonard2014_StableContinental` implement the magnitude–area
subset, with median constants following the OpenQuake hazardlib
convention. Normal and reverse use the same dip-slip constant; the
`all` constant averages the relevant fault-style constants. No
length/width law or calibrated residual sigma is asserted here.

**Hanks & Bakun (2002).** *Bulletin of the Seismological Society of
America*, 92, 1841–1846.
https://doi.org/10.1785/0120010148.
The implemented bilinear area relation is

```text
Mw = log10(A) + 3.98             for A <= 537 km²
Mw = (4/3) log10(A) + 3.07       for A >  537 km²
```

The magnitude-to-area inverse switches at the **lower-branch Mw
value** for 537 km². The rounded published coefficients give a small
discontinuity around the transition. `HanksBakun2002` has only
`style="all"`; no mechanism-specific variants are invented.

The older `hazard/magnitude_scaling.py` module has been removed;
use this independent `groundmotion.scaling` module.

## 2. Standalone model selection and evaluation

```python
from shakelab.groundmotion.scaling import (
    RuptureDimensions,
    available_scaling_relations,
    create_scaling_relation,
    fault_style_from_rake,
)

print(available_scaling_relations())

style = fault_style_from_rake(90.0)  # "reverse"
model = create_scaling_relation(
    "WellsCoppersmith1994",
    style=style,
)

area_estimate = model.magnitude_to_area(6.2)
print(area_estimate.median, area_estimate.sigma)

dimensions = RuptureDimensions.from_area(
    area_estimate.median,
    aspect_ratio=2.0,
    max_width=15.0,
)
print(dimensions.length, dimensions.width, dimensions.area)

magnitude_estimate = model.area_to_magnitude(200.0)
print(magnitude_estimate.median)
```

`RuptureDimensions.from_area()` requires **exactly one** of
`aspect_ratio`, `length`, or `width`. If a maximum width is given
and an aspect ratio defined the shape, the width may be capped and
length increased to preserve area. When a dimension was explicit,
a conflicting cap raises an error rather than changing it.

This helper deals with **dimensions only**. It does not define the
hypocentre, dip, fault vertices, seismogenic depth or rupture-to-site
distances. For those, use
`build_rectangular_rupture()` or `RectangularRupture` as described
in [RUPTURE.md](RUPTURE.md).

## 3. Fault style and tectonic domain

`fault_style_from_rake(rake)` uses these conventions:

| Rake | Automatic crustal style |
| --- | --- |
| Strictly between +45° and +135° | `reverse` |
| Strictly between -135° and -45° | `normal` |
| All other valid angles, including boundaries | `strike-slip` |
| `None` | `all` |

The valid range is [-180°, 180°]. The ±45° and ±135° boundaries are
assigned to strike-slip. The function is a **classification
convention**, not an additional empirical magnitude-scaling law.

For `ThingbaijamEtAl2017`, automatic selection with an unknown rake
would yield `all`, which is unsupported: provide a rake or choose
an explicit supported style. For `interface`, always choose
`style="interface"` explicitly; a reverse rake alone cannot
establish the interface tectonic domain.

`HanksBakun2002` uses `all` regardless of rake when a provider's
`style="auto"` is selected. Other models use the rake-derived style
in their automatic provider path.

## 4. Predictions, dispersion and validity

The returned `ScalingEstimate` contains:

| Field | Meaning |
| --- | --- |
| `median` | Central predicted size or magnitude |
| `sigma` | Residual standard deviation, or `None` if unavailable |
| `sigma_space` | Regression-space label (`log10` for size, `Mw` for reverse magnitude) |
| `extrapolated` | `True` / `False` when calibration bounds are implemented; otherwise `None` |

For published magnitude intervals, forward predictions outside
calibration **warn by default**, while `strict=True` rejects them.
Reverse diagnostics use their implemented range checks. The
absence of calibration bounds is encoded as `None`, **not** as
`extrapolated=False`.

`ScalingEstimate.at_epsilon(epsilon)` applies a signed number of
residual standard deviations, deterministically. For log10 size:

```text
size(epsilon) = median * 10 ** (epsilon * sigma)
```

For reverse Mw estimates it adds `epsilon * sigma` in magnitude
units. `epsilon=0` returns the median. A nonzero epsilon requires
an available residual dispersion; no sigma is invented for Leonard
or Hanks–Bakun.

```python
estimate = create_scaling_relation(
    "WellsCoppersmith1994",
    style="reverse",
).magnitude_to_area(6.5)

upper_area = estimate.at_epsilon(1.0)
```

A residual offset is **not** a joint stochastic draw of area, length,
width and earthquake magnitude. In particular, fitting separate
length and width medians does not guarantee that their product equals
the independently estimated median area.

For Leonard and Hanks–Bakun, the implemented calibration intervals
and residual prediction sigmas are unknown. Thus `strict=True` and
nonzero `epsilon` are rejected, rather than implying validated
ranges or zero uncertainty. Standard errors of fitted coefficients
are not residual prediction standard deviations.

## 5. Rupture area versus finite-rupture geometry

Magnitude–size relations supply empirical dimensions, **not fault
location or orientation**. The constrained rupture builder uses
estimated **area** when needed, with a caller-supplied shape
constraint (`length`, `width`, or `aspect_ratio`) or sufficient
depth/nucleation constraints.

```python
from shakelab.groundmotion.rupture import build_rectangular_rupture
from shakelab.libutils.geodeticN.primitives import WgsPoint

result = build_rectangular_rupture(
    WgsPoint(13.0, 46.0, -8000.0),
    strike=120.0,
    dip=45.0,
    rake=90.0,
    magnitude=6.5,
    scaling_model="WellsCoppersmith1994",
    aspect_ratio=2.0,
    lower_seismogenic_depth=15.0,
)

assert result.scaling_estimate is not None
print(result.rupture.area, result.parameter_origins["area"])
```

If sufficient dimensions or an explicit area are supplied, the
builder does not need to invoke magnitude scaling. Explicit and
redundant constraints are validated for consistency. When the width
is derived from area and aspect ratio, the seismogenic layer can cap
it while preserving area by extending length; **explicit dimensions
are never silently altered**.

The builder does not infer strike, dip, rake, tectonic domain or a
universal aspect ratio. A central hypocentre position is the initial
default only when nucleation fractions are not specified. See
[RUPTURE.md](RUPTURE.md) for its depth/nucleation rules and the
meaning of `RuptureBuildResult.parameter_origins`.

## 6. Optional scaling in standalone GMPE providers

A GMPE provider declares a scaling relation in `source_scaling`;
orientation, dimensions and seismogenic depth bounds remain in
`source_parameters`.

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)

event = ScenarioEvent.from_lonlat_depth(
    6.5, 13.0, 46.0, 8.0, magnitude_type="Mw",
)

provider = GroundMotionProvider.gmpe(
    "AbrahamsonEtAl2014",
    config={
        "source_geometry": "finite",
        "source_parameters": {
            "strike": 120.0,
            "dip": 45.0,
            "rake": 90.0,
            "aspect_ratio": 2.0,
            "lower_seismogenic_depth": 15.0,
        },
        "source_scaling": {
            "model": "ThingbaijamEtAl2017",
            "style": "auto",
            "strict": False,
            "epsilon": 0.0,
        },
    },
)

resolved = provider.resolve_source(event)
print(resolved.rupture.length, resolved.rupture.width)
print(resolved.build_result.scaling_model)
```

The accepted `source_scaling` keys are:

| Key | Default | Meaning |
| --- | --- | --- |
| `model` | Required | Exact relation identifier |
| `style` | `"auto"` | Resolve supported fault style using rake, or provide one explicitly |
| `strict` | `False` | Reject extrapolation when calibration bounds exist |
| `epsilon` | `0.0` | Deterministic residual offset for estimated area |

Unknown keys and incompatible model/style combinations are rejected
at configuration validation. `source_scaling` requires
`source_geometry="finite"`. The scaling style is resolved at
**evaluation time** when `"auto"` is selected because rake can
come from the event or a runtime override.

The provider merges source parameters with precedence
**provider < event < runtime**, resolves a single consistent finite
rupture, and uses its final parameters in both the GMPE rupture
context and distance calculation. An explicit distance override to
`evaluate(distances=...)` still wins over a calculated distance.
Its source resolver can be called independently as
`provider.resolve_source(event, rupture_parameters=None)`.

If dimensions or area are already sufficient, scaling may be
configured but **not actually used**. Inspect
`resolved.build_result.scaling_model` and
`resolved.build_result.scaling_estimate` to distinguish configured
from applied scaling.


## 7. Practical scaling recipes and comparisons

The examples here operate entirely on the **independent scaling
module**. They do not require ShakeScenario, a GMPE, or an exposure
model. Use physically appropriate Mw and tectonic setting for any
scientific application.

### 7.1. Compare predicted area across implemented relations

Models with different calibration domains should not be interpreted
as interchangeable estimates of uncertainty. For a numerical
comparison, first select valid styles explicitly:

```python
from shakelab.groundmotion.scaling import create_scaling_relation

models_and_styles = {
    "WellsCoppersmith1994": "reverse",
    "ThingbaijamEtAl2017": "reverse",
    "Leonard2014_Interplate": "reverse",
    "Leonard2014_StableContinental": "reverse",
    "HanksBakun2002": "all",
}

for model_name, style in models_and_styles.items():
    relation = create_scaling_relation(model_name, style=style)
    estimate = relation.magnitude_to_area(6.5)
    print(
        model_name,
        f"area={estimate.median:.1f} km²",
        f"sigma={estimate.sigma}",
        f"extrapolated={estimate.extrapolated}",
    )
```

A `None` sigma means that this implementation does **not** supply
a residual prediction dispersion; it does not indicate a
deterministic or perfectly constrained relationship.

### 7.2. Different mechanisms and explicit interface selection

Rake supports an approximate crustal mechanism classification,
but cannot identify a subduction-interface tectonic environment:

```python
from shakelab.groundmotion.scaling import (
    create_scaling_relation, fault_style_from_rake,
)

for rake in (-90.0, 0.0, 90.0):
    style = fault_style_from_rake(rake)
    relation = create_scaling_relation(
        "ThingbaijamEtAl2017", style=style,
    )
    estimate = relation.magnitude_to_area(7.0)
    print(rake, style, estimate.median)

interface = create_scaling_relation(
    "ThingbaijamEtAl2017", style="interface",
)
print("Interface:", interface.magnitude_to_area(7.0).median)
```

The interface style is an **explicit tectonic decision**, not the
output of `fault_style_from_rake()`.

### 7.3. Residual offsets are deterministic, not Monte Carlo draws

```python
from shakelab.groundmotion.scaling import create_scaling_relation

relation = create_scaling_relation(
    "WellsCoppersmith1994", style="reverse",
)
estimate = relation.magnitude_to_area(6.5)
low = estimate.at_epsilon(-1.0)
central = estimate.at_epsilon(0.0)
high = estimate.at_epsilon(1.0)
print(low, central, high)
assert low < central < high
```

These are three points along the **same area regression**. They are
not sampled lengths, widths, hypocentral positions or a joint rupture
distribution. A nonzero epsilon for `HanksBakun2002` and the
two Leonard variants is rejected because their residual sigmas
are not implemented.

### 7.4. Forward and reverse regressions are not always inverses

For Wells–Coppersmith, the regressions from Mw to area and from
area to Mw are independently fitted:

```python
from shakelab.groundmotion.scaling import create_scaling_relation

relation = create_scaling_relation(
    "WellsCoppersmith1994", style="reverse",
)
area = relation.magnitude_to_area(6.5).median
estimated_mw = relation.area_to_magnitude(area).median
print("Input Mw:", 6.5, "Reverse estimate:", estimated_mw)
```

The two numbers need not agree. For Thingbaijam, the implemented
reverse path instead algebraically inverts the forward fitted line
and propagates the associated residual scale accordingly.

### 7.5. Convert one area into alternative rectangle shapes

```python
from shakelab.groundmotion.scaling import (
    RuptureDimensions, create_scaling_relation,
)

relation = create_scaling_relation(
    "WellsCoppersmith1994", style="reverse",
)
area = relation.magnitude_to_area(6.5).median

by_ratio = RuptureDimensions.from_area(area, aspect_ratio=2.0)
by_length = RuptureDimensions.from_area(area, length=20.0)
by_width = RuptureDimensions.from_area(area, width=10.0)

for shape in (by_ratio, by_length, by_width):
    print(shape.length, shape.width, shape.area)
```

All three preserve the area estimate while making **different
assumptions about rupture shape**. The empirical area relationship
alone cannot select among them.

### 7.6. Reject unsupported styles and calibration assumptions

The API validates model/style compatibility; unknown uncertainty
information is not silently invented:

```python
from shakelab.groundmotion.scaling import create_scaling_relation

try:
    create_scaling_relation(
        "ThingbaijamEtAl2017", style="all",
    )
except ValueError as exc:
    print("Unsupported style:", exc)

try:
    create_scaling_relation(
        "HanksBakun2002", style="all", strict=True,
    )
except ValueError as exc:
    print("Calibration bounds unavailable:", exc)
```

When supported, `strict=True` rejects extrapolation. For other
models, `extrapolated=None` means calibration limits are unknown
in this implementation.

### Temporary magnitude-type policy

The empirical models are calibrated against **Mw**. Currently the
GMPE provider uses the event's magnitude **numerically as given**,
even if `magnitude_type` says `ML` or another scale. It does not
perform a physical conversion or reject a magnitude-type mismatch.
This is a provisional operational assumption, **not a scientific
equivalence between magnitude scales**. Standalone empirical model
calls should use Mw, and `strict` checks still apply to the
supplied numerical magnitude.

No magnitude conversion, stochastic finite-fault sampling, fault
orientation inference or joint dimension uncertainty modelling is
implemented by this module.
