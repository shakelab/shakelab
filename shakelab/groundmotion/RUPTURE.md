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

Scaling is opt-in and assumes Mw. Explicit dimensions or area take
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
The builder is not yet connected to the GMPE providers.
