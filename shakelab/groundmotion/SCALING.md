# Magnitude scaling

The independent module `shakelab.groundmotion.scaling` implements
Wells and Coppersmith (1994), BSSA 84(4), 974–1002, Table 2A,
https://doi.org/10.1785/BSSA0840040974.

Use moment magnitude (Mw), lengths and widths in km and area in km².

```python
from shakelab.groundmotion.scaling import (
    WellsCoppersmith1994, RuptureDimensions, fault_style_from_rake,
)

model = WellsCoppersmith1994(fault_style_from_rake(90))
estimate = model.magnitude_to_area(6.2)
dimensions = RuptureDimensions.from_area(
    estimate.median, aspect_ratio=2.0, max_width=15.0,
)
```

The default fault style is all. Forward predictions return a median
and a residual sigma in log10 units; reverse predictions return a
central magnitude and sigma in magnitude units. The two regressions
are independent, not inverses. Use at_epsilon() for deterministic
residual offsets, not for sampling a joint dimension distribution.

Forward estimates flag and warn about extrapolation outside the
published Mw range; strict=True rejects it. Reverse length and area
use the published size ranges. The printed width range in Table 2A
repeats the subsurface length range, so the reverse width diagnostic
uses the predicted Mw range instead, not a claimed width domain.

Independent median length, width and area are not necessarily
geometrically consistent. RuptureDimensions.from_area preserves
area with one explicit shape constraint: aspect ratio, length or
width. A width cap may change the aspect ratio, but cannot override
an explicitly specified length or width. No seismogenic thickness
or nucleation position is assumed. Scaling is not yet activated in
the GMPE providers.

Rake classification uses conventional 45-degree boundaries;
boundary values select strike-slip. This classification is not a
separate empirical relationship from the original paper.

The historical `shakelab.hazard.magnitude_scaling` module delegates its
scalar area/magnitude methods to this implementation. New applications
should use the class in `groundmotion.scaling` directly.
