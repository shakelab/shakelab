# Finite rupture geometry

`shakelab.groundmotion.rupture.RectangularRupture` is independent of
GMPEs and scenario providers. It represents a rectangular plane in a
WGS84 local East-North-Up frame using the existing geodeticN transforms.

Lengths, depths, local coordinates exposed by this module and distances
are in kilometres. WgsPoint elevation is in metres, positive upward.
Strike is measured clockwise from North; dip is measured downward
from horizontal, to the right of strike.

Construct a rupture from its upper-edge midpoint:

```python
from shakelab.groundmotion.rupture import RectangularRupture
from shakelab.libutils.geodeticN.primitives import WgsPoint

rupture = RectangularRupture(
    WgsPoint(13.0, 46.0, -2000.0),
    strike=120.0, dip=35.0, length=20.0, width=12.0,
)
distances = rupture.distances([WgsPoint(13.1, 46.1)])
```

Or define a plane containing a given hypocentre:

```python
rupture = RectangularRupture.from_hypocentre(
    WgsPoint(13.0, 46.0, -7000.0),
    strike=120.0, dip=35.0, length=20.0, width=12.0,
    ztor=2.0, along_fraction=0.5,
)
```

The along-strike hypocentral position is explicit. The down-dip
position is constrained by depth and ztor, except for a horizontal
rupture, which also requires an explicit down_fraction. The constructor
rejects an inconsistent hypocentre or a rupture reaching above the
local surface.

Vertices are ordered upper-start, upper-end, lower-end, lower-start.
Distances are analytic in the local Cartesian representation:

- rrup: shortest three-dimensional distance to the finite plane;
- rjb: shortest horizontal distance to the projected rectangle;
- rx: signed cross-strike distance to the infinite upper-edge line;
- ry0: along-strike distance outside the interval spanned by the edge.

Site elevations are retained in the ENU transform. This representation
uses one local tangent plane: depths and horizontal distances are local,
not distances to a curved ellipsoidal polygon. It is intended for
regional finite-rupture modelling; Earth curvature and surface
topography need separate consideration for large extents.
The RectangularRupture class itself performs no magnitude scaling or
fault selection. The separate builder below can use magnitude scaling.
The GMPE provider uses this class when source_geometry is finite; the
geometry itself remains independent of the provider.


## Automatic rectangular rupture construction

```python
from shakelab.groundmotion.rupture import build_rectangular_rupture

result = build_rectangular_rupture(
    WgsPoint(13.0, 46.0, -8000.0),
    strike=120, dip=45, rake=90, magnitude=6.5,
    scaling_model="WellsCoppersmith1994", aspect_ratio=2,
    lower_seismogenic_depth=15,
)
rupture = result.rupture
```

Scaling is opt-in and was calibrated using Mw. At present the GMPE
provider passes the numerical scenario magnitude unchanged, regardless
of its magnitude label; this is a provisional assumption, not a
magnitude conversion. Explicit dimensions or area take
precedence. At least one shape constraint is necessary, except where
the depth and nucleation constraints fix width. Redundant explicit
constraints must agree. When only aspect ratio defines shape, width
may be adjusted to fit the depth bounds, preserving area. Dimensions
fixed by other explicit constraints are not changed.

The default along-strike fraction is 0.5. Where neither ztor nor
down_fraction is supplied, the rupture is initially centred in depth,
then shifted to fit the permitted layer. An explicit fraction is
never shifted. Upper depth defaults to zero; lower depth is optional.
Horizontal planes have no vertical width constraint.

The result includes resolved source parameters and parameter origins,
with scaling model, style and estimate recorded only when used.
With `source_scaling` configured, GMPE providers call the builder
and retain its `RuptureBuildResult` in the resolved source snapshot.

## Shared GMPE source resolution

The GMPE provider resolves provider, event and runtime source parameters
into a single `ResolvedSource` snapshot per `evaluate()` call. In finite
mode the explicit rectangular rupture is built once and its final
geometry is used for both rupture parameters and site distances.
Explicit distance overrides still take precedence, and point/native
modes retain their existing numerical behavior. Standalone calls to
`_build_sites_context()` remain supported for existing integrations.
Finite sources are now retained in an instance-local, bounded LRU
cache (64 entries) for repeated evaluations of the same effective
source. The key includes hypocentral longitude, latitude and elevation,
numerical magnitude, finite geometry mode, merged source parameters
and scaling options. Site parameters, IMT and distance overrides are
excluded because they do not determine rupture geometry. Source
construction and LRU updates are protected by a reentrant lock.
A changed event, source override or scaling configuration yields a
different key; no cache state is global or shared across provider
instances. Native and point modes are never cached. Cached finite
sources may be reused across assets and IMTs, while distances are
computed for each site and explicit distance overrides still apply.
The cache is an internal provider optimization; resolved source
objects should be treated as read-only.

No scaling is performed for native or point geometry, even when
an event overrides a finite provider's geometry selection.

Magnitude labels are provisionally treated as equivalent in provider
evaluation; no magnitude conversion is performed. GMPE-specific
magnitude-type mismatches no longer stop evaluation.
