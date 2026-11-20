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

## Model selection

Use available_scaling_relations() to discover supported names and
create_scaling_relation(name, **kwargs) to select a model.

ThingbaijamEtAl2017 implements Table 1 of Thingbaijam, Mai and Goda,
BSSA 107, 2225-2246, doi:10.1785/0120170017. Specify style as
reverse (shallow crustal), normal, strike-slip or interface.
No pooled all-style fit is supplied. Interface is a tectonic setting
and must not be inferred solely from rake.

The forward magnitude ranges are style specific; reverse estimates
use the observed dimension ranges in Table 1. Reverse values invert
the fitted line, with residual scale sigma/abs(slope), rather than
claiming a separately fitted conditional magnitude distribution.
Unlike WC1994, TMG2017 area coefficients are derived from the length
and width regressions. No joint stochastic sampling is implemented.

The obsolete hazard/magnitude_scaling.py module has been removed.


## Additional magnitude-area models

Leonard2014_Interplate and Leonard2014_StableContinental implement the magnitude-area
subset of Leonard (2014), doi:10.1785/0120140087, using median constants
as in OpenQuake hazardlib. Interplate means active crustal faults here,
not a dedicated subduction interface model. Normal and reverse share
the dip-slip constant; all uses the mean of the two style constants.

HanksBakun2002 (doi:10.1785/0120010148) uses Mw = log10(A) + 3.98 up
to 537 km² and Mw = (4/3) log10(A) + 3.07 above that area. The inverse
branches switch at the lower-branch magnitude for 537 km². Rounded
published coefficients produce a small discontinuity around the break.

These models support area only. Use RuptureDimensions to supply shape
constraints. Residual sigma and calibration status are None (unknown),
not zero and not confirmed in-range. Nonzero epsilon and strict=True
are rejected rather than silently assuming uncertainty or calibration
limits. No mechanism-specific variants are invented for Hanks–Bakun.
