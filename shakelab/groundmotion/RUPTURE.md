# Earthquake source representation and rupture geometry

This guide describes the scientific source representation used by
`shakelab.groundmotion.rupture` and the GMPE backend of
`shakelab.groundmotion.providers`. It applies to standalone Python
applications; no impact model, server, or application-specific
configuration is needed.

For empirical magnitude–size relations, their calibration, and the
meaning of `epsilon` and `strict`, see [SCALING.md](SCALING.md).

## 1. Source representation modes

A GMPE provider accepts exactly three values of `source_geometry`:

| Mode | Representation | Distance handling | Magnitude scaling |
| --- | --- | --- | --- |
| `native` (default) | No rupture surface constructed | Explicit finite-fault metrics; `repi` and `rhypo` may be computed from the event and sites | Not used |
| `point` | Hypocentre, optionally with strike | Point-source approximations for supported metrics | Not used |
| `finite`, explicit | Rectangular plane specified by geometric inputs | Analytic distances to the plane | Not used without `source_scaling` |
| `finite`, constructed | Rectangular plane completed by `build_rectangular_rupture()` | Analytic distances to the resolved plane | Opt-in, used only when needed |

The last two rows are **two construction paths for the same `finite`
mode**, not additional `source_geometry` values.

`source_geometry=None` resolves to `native` unless the legacy
`point_source_distances` flag is used; `True` selects `point` and
`False` selects `native`. Do not provide both configuration options.

The event may specify `source_geometry` to override the provider's
selected mode. A provider configured with `source_scaling` must have
`source_geometry="finite"` at configuration time, but an explicit
event-level override to `native` or `point` bypasses finite
construction for that evaluation.

## 2. Conventions and source parameters

`ScenarioEvent` is the general event container used by the
ground-motion provider. It holds a numerical magnitude and WGS84
hypocentre; optional metadata include `magnitude_type`, `event_id`,
`source_parameters`, and `source_geometry`.

```python
from shakelab.groundmotion.providers import ScenarioEvent

event = ScenarioEvent.from_lonlat_depth(
    magnitude=6.2,
    longitude=13.0,
    latitude=46.0,
    depth_km=7.0,
    magnitude_type="Mw",
)
```

The location and orientation conventions are:

| Quantity | Convention |
| --- | --- |
| Longitude / latitude | WGS84 degrees |
| `WgsPoint.elevation` | Metres, positive upward |
| Event `depth_km` | Kilometres, positive downward |
| `strike` | Degrees clockwise from North, in [0, 360] |
| `dip` | Degrees from horizontal, in [0, 90], dipping right of strike |
| `rake` | Degrees in [-180, 180] |
| `length`, `width`, `ztor` | Kilometres |
| `area` | Square kilometres |
| `aspect_ratio` | Dimensionless, `length / width` |
| `along_fraction` | Hypocentre position along strike, from 0 to 1 |
| `down_fraction` | Hypocentre position down dip, from 0 to 1 |
| `upper_seismogenic_depth`, `lower_seismogenic_depth` | Kilometres, positive downward |

`ztor` is the depth to the top edge of the rectangular rupture. The
epicentre is the surface projection of the hypocentre.

The provider accepts additional finite numerical values in
`source_parameters` for GMPE-specific rupture context fields.
`mag`, `magnitude`, `hypo_depth`, `depth`, and `depth_km` are
reserved and cannot be supplied as source overrides. The event
magnitude and hypocentral depth come from `ScenarioEvent`.

### Parameter precedence

For source parameters, a later layer overrides an earlier one:

`provider.source_parameters < event.source_parameters < evaluate(rupture_parameters=...)`

Runtime `None` values are ignored as overrides. The final merged
parameters are used consistently to construct the source, build the
GMPE rupture context, and compute finite-fault distances.

This input hierarchy is separate from **geometric constraints**:
explicit physical dimensions and area take precedence over magnitude
scaling, but mutually inconsistent explicit constraints are rejected
rather than silently overwritten.

Site properties, such as `vs30`, belong in `site_parameters`, not
`source_parameters`. Site-specific distance overrides belong in
`distances`, not in either source mapping.

## 3. Native and point-source modes

### Native: no invented finite-fault geometry

`native` is the default. The provider does not construct a rupture
plane. Supply `rrup`, `rjb`, `rx`, `ry0`, or other finite-fault
metrics explicitly when required by the chosen GMPE. The elementary
`repi` and `rhypo` distances can be computed from the event and
sites when those metrics are requested; their availability does not
imply a finite-rupture approximation.

For a GMPE that requires `rjb`, for example, an explicit distance can
be supplied as follows:

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint

event = ScenarioEvent.from_lonlat_depth(6.2, 13.0, 46.0, 7.0)
provider = GroundMotionProvider.gmpe(
    "KothaEtAl2020ESHM20",
    config={"source_geometry": "native"},
)
sites = [WgsPoint(13.1, 46.1)]

median, sigma_ln = provider.evaluate(
    "PGA",
    sites,
    event,
    site_parameters={"vs30": 500.0, "eshm20_region": 0},
    distances={"rjb": 15.0},
)
```

The required site, rupture, and distance fields depend on the GMPE.
An explicit distance takes precedence over any distance otherwise
derived by the provider.

### Point: explicitly approximate finite-fault metrics

With `source_geometry="point"`, the provider uses:

| Requested metric | Point approximation |
| --- | --- |
| `rjb` | `repi` |
| `rrup` | `rhypo` |
| `rx` | Signed cross-strike epicentral offset, if strike is available |
| `ry0` | Absolute along-strike epicentral offset, if strike is available |

The oriented offsets are projections about the **epicentre**, not
distances to a finite plane. A required `rx` or `ry0` cannot be
derived in point mode without strike or an explicit distance. These
approximations should not be interpreted as finite-fault geometry.

`point` still passes required rupture parameters (for example
`width`, `ztor`, `dip`, or `rake`) to the GMPE when the selected
model requires them. The mode only changes **distance construction**;
it does not automatically invent missing GMPE rupture parameters.

`point_source_distances=True` is a supported legacy alias for
`source_geometry="point"`. Prefer the explicit new option.

## 4. Explicit rectangular rupture geometry

`RectangularRupture` represents an oriented planar rectangle in a
local WGS84 East–North–Up (ENU) tangent frame. It is independent of
GMPEs, empirical scaling, and application frameworks.

### Construct from the upper-edge reference point

```python
from shakelab.groundmotion.rupture import RectangularRupture
from shakelab.libutils.geodeticN.primitives import WgsPoint

rupture = RectangularRupture(
    WgsPoint(13.0, 46.0, -2000.0),
    strike=120.0,
    dip=35.0,
    length=20.0,
    width=12.0,
)

vertices = rupture.vertices     # Property: tuple of four WgsPoint values
distances = rupture.distances([WgsPoint(13.1, 46.1)])
```

The reference is at `along_fraction=0.5` and `down_fraction=0.0`
by default: the midpoint of the upper edge. Vertices follow the order
upper-start, upper-end, lower-end, lower-start.

### Construct a rectangle containing an existing hypocentre

```python
hypocentre = WgsPoint(13.0, 46.0, -7000.0)

rupture = RectangularRupture.from_hypocentre(
    hypocentre,
    strike=120.0,
    dip=35.0,
    length=20.0,
    width=12.0,
    ztor=2.0,
    along_fraction=0.5,
)
```

When `ztor` and dip constrain the hypocentral depth, the down-dip
fraction is derived. For a horizontal plane, `ztor` must equal the
hypocentral depth and an appropriate `down_fraction` is needed.
Inconsistent hypocentre placement or a rupture extending above the
local surface is rejected.

### Finite-fault distance definitions

`rupture.distances(sites)` returns analytic arrays in kilometres:

| Key | Definition |
| --- | --- |
| `rrup` | Shortest three-dimensional distance to the rectangle |
| `rjb` | Shortest horizontal distance to its surface projection |
| `rx` | Signed cross-strike offset from the infinite upper-edge line; positive to the right of strike |
| `ry0` | Along-strike distance outside the upper-edge interval |

Site elevations are included in the ENU coordinate transformation.
The geometry uses one local tangent plane, not a polygon draped over a
curved Earth. Earth curvature, topography, nonplanar faults and
multi-segment rupture surfaces are not modelled by this class.
It is most appropriate for regional rectangular-fault approximations.

### Explicit finite geometry in a GMPE provider

The legacy explicit-construction path is selected by
`source_geometry="finite"` **without** `source_scaling`. The
provider requires all six entries:

`strike`, `dip`, `length`, `width`, `ztor`,
`along_fraction`.

`down_fraction` is optional and must be compatible if supplied.
This path does **not** automatically estimate missing dimensions.
The following provider can be used directly in Python:

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint

event = ScenarioEvent.from_lonlat_depth(6.2, 13.0, 46.0, 7.0)
sites = [WgsPoint(13.1, 46.1)]

provider = GroundMotionProvider.gmpe(
    "AbrahamsonEtAl2014",
    config={
        "source_geometry": "finite",
        "source_parameters": {
            "strike": 120.0,
            "dip": 35.0,
            "rake": 90.0,
            "length": 20.0,
            "width": 12.0,
            "ztor": 2.0,
            "along_fraction": 0.5,
        },
    },
)

median, sigma_ln = provider.evaluate(
    "PGA",
    sites,
    event,
    site_parameters={
        "vs30": 500.0,
        "vs30measured": True,
        "z1pt0": 300.0,
    },
)
```

The event used in this example has a 7 km hypocentral depth.
Its finite geometry must physically contain that hypocentre.

## 5. Constrained rupture construction

The independent `build_rectangular_rupture()` function can complete
an incomplete rectangular geometry from area, shape information, and
optional magnitude scaling. It returns `RuptureBuildResult`, not
just a `RectangularRupture`.

```python
from shakelab.groundmotion.rupture import build_rectangular_rupture

result = build_rectangular_rupture(
    WgsPoint(13.0, 46.0, -8000.0),
    strike=120.0,
    dip=45.0,
    rake=90.0,
    magnitude=6.5,
    aspect_ratio=2.0,
    scaling_model="WellsCoppersmith1994",
    lower_seismogenic_depth=15.0,
)
rupture = result.rupture
print(rupture.length, rupture.width, rupture.ztor)
print(result.parameter_origins)
```

No focal angles, tectonic domain, or aspect ratio are inferred. At
least one shape constraint (`length`, `width`, or `aspect_ratio`)
is normally needed; alternatively, depth/nucleation constraints may
determine width. The function accepts an explicit `area`, so a
scaling model is not needed if area or the dimensions are sufficient.

### Dimension and constraint rules

1. If both length and width are explicit, their product defines area;
   any additional area or aspect ratio must be consistent.
2. If area and one dimension are explicit, the other is derived.
3. If area and aspect ratio are supplied, the initial dimensions are
   obtained from their ratio and product.
4. When area is absent and dimensions are incomplete, a configured
   scaling relation may estimate area from the numerical magnitude.
5. If seismogenic depth constraints limit a *derived* width and only
   aspect ratio fixes the shape, width can be reduced while length
   increases to preserve area; the resulting aspect ratio is reported
   as constrained.
6. An explicitly fixed dimension is **never silently modified** to
   fit the seismogenic layer. Contradictions raise `ValueError`.

The `upper_seismogenic_depth` defaults to 0 km;
`lower_seismogenic_depth` is unbounded by default. The hypocentre
must lie in the permitted layer. An explicitly specified `ztor` or
`down_fraction` is a constraint, not a preference.

By default, `along_fraction=0.5`. When neither `ztor` nor
`down_fraction` is supplied, the nucleation position is initially
centred down dip and can move to fit the layer. An explicit nucleation
fraction is never silently changed. Horizontal planes have no
vertical width constraint.

The result exposes:

- `rupture`: the final `RectangularRupture`.
- `source_parameters`: resolved geometry and available rake.
- `parameter_origins`: `explicit`, `scaling`, `derived`,
  `default`, or `constrained` for applicable parameters.
- `scaling_model`, `scaling_style`, and `scaling_estimate`:
  populated only when a relation actually supplied area.

An explicitly complete rupture can therefore be built even when a
scaling model is configured: the model is **not used** unless needed
to determine area.

## 6. Automatic construction in GMPE providers

Magnitude scaling is **opt-in**. The provider requires
`source_geometry="finite"` when `source_scaling` is configured:

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)

event = ScenarioEvent.from_lonlat_depth(6.2, 13.0, 46.0, 7.0)
provider = GroundMotionProvider.gmpe(
    "AbrahamsonEtAl2014",
    config={
        "source_geometry": "finite",
        "source_parameters": {
            "strike": 120.0,
            "dip": 35.0,
            "rake": 90.0,
            "aspect_ratio": 2.0,
            "lower_seismogenic_depth": 15.0,
        },
        "source_scaling": {
            "model": "ThingbaijamEtAl2017",
            "style": "auto",
        },
    },
)

resolved = provider.resolve_source(event)
rupture = resolved.rupture
print(rupture.vertices)
print(resolved.build_result.parameter_origins)
```

The provider passes the event's numerical magnitude to the builder,
merges source inputs once, and uses the **final** width, `ztor`, dip
and other resolved fields in both the rupture context and distance
calculation. Explicit distance overrides to `evaluate(distances=...)`
still take precedence; they do not alter the physical rupture.

When `source_scaling` is omitted, the explicit finite path described
in section 4 is used instead. When it is present, the builder may be
used even if scaling is not needed because dimensions are complete.
See [SCALING.md](SCALING.md) for models, fault styles and diagnostics.

## 7. Standalone source resolution, caching and trace recording

`provider.resolve_source(event, rupture_parameters=None)` is a
public scientific API and does **not** evaluate a GMPE. It returns a
`ResolvedSource` with:

- `geometry`: effective `native`, `point`, or `finite` mode.
- `parameters`: merged and resolved rupture parameters.
- `rupture`: a finite `RectangularRupture`, or `None`.
- `build_result`: a `RuptureBuildResult` when the automatic
  builder was used, otherwise `None`.

The same resolution logic is used internally by
`provider.evaluate(imt, sites, event, ...)`. The returned snapshot
should be treated as read-only.

Finite sources are cached **per provider instance** using a bounded
64-entry LRU protected by a reentrant lock. Keys include the
hypocentral coordinates and elevation, numerical magnitude,
effective source mode, merged source parameters, and scaling options.
Site parameters, IMT and explicit distance overrides do not affect
rupture geometry and are not part of the key. Point and native source
resolutions are not cached. Distances remain site-dependent and are
calculated for each evaluation as needed.

### Optional source trace

Generic GMPE users do not need trace collection. It is disabled by
default; explicitly set `record_source_trace=True` when needed:

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint

event = ScenarioEvent.from_lonlat_depth(6.2, 13.0, 46.0, 7.0)
sites = [WgsPoint(13.1, 46.1)]

provider = GroundMotionProvider.gmpe(
    "AbrahamsonEtAl2014",
    config={
        "source_geometry": "finite",
        "source_parameters": {
            "strike": 120.0, "dip": 35.0, "rake": 90.0,
            "length": 20.0, "width": 12.0,
            "ztor": 2.0, "along_fraction": 0.5,
        },
        "record_source_trace": True,
    },
)
provider.evaluate(
    "PGA", sites, event,
    site_parameters={"vs30": 500.0,
                     "vs30measured": True, "z1pt0": 300.0},
)
records = provider.source_trace()
```

`source_trace()` returns detached, JSON-compatible records of
successful GMPE evaluations. They include resolved parameters, finite
geometry when present, input-layer precedence, builder origins (if
available), scaling estimates (if actually used), site coordinates,
and the distance values supplied to the GMPE together with their
provenance. Trace recording is separate from the numerical algorithm:
enabling it must not change GMPE results. A trace-generation error
is not silently suppressed when recording is explicitly enabled.
No output file format or file path is prescribed by `groundmotion`.

## 8. Scope and scientific assumptions

- The finite source is a **single planar rectangle**, not a
  multi-fault, nonplanar or dynamically propagating rupture.
- Point-source finite-distance approximations are useful only when
  their limitations are acceptable for the chosen GMPE and event.
- Required GMPE context fields remain model-specific. Geometry mode
  does not fill unrelated missing rupture/site parameters.
- Empirical scaling models were calibrated for **Mw**. Currently the
  GMPE provider uses the event's numerical magnitude without
  converting between Mw, ML or other labels, and does not reject a
  magnitude-type mismatch. This is a **provisional operational
  approximation**, not a physical equivalence of magnitude scales.
  Direct scientific use of the scaling relations should supply Mw.
- Seismogenic depth bounds, fault style and aspect ratio are not
  inferred from geographical location. They must be chosen by the
  caller from appropriate scientific evidence.
