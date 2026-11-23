# Ground-motion prediction equations (GMPEs)

This guide describes the ground-motion models registered in
`shakelab.groundmotion.gmpe`, how to evaluate them independently of
application frameworks, and how to implement and validate an additional
model. The primary sources of truth are
[`registry.json`](registry.json), [`base.py`](base.py), the individual
model modules, and their coefficient files under [`data/`](data/).

For earthquake source representations, rupture geometry, and automatic
rupture construction see [`../RUPTURE.md`](../RUPTURE.md). For empirical
magnitude–rupture-size relations see
[`../SCALING.md`](../SCALING.md). These are separate from GMPE
coefficients and their aleatory variability.

## 1. Scientific and numerical conventions

A GMPE evaluates an **intensity measure type (IMT)** for one rupture
and a collection of sites. The low-level contract is:

- `RuptureContext`: scalar rupture properties (e.g., `mag`, `rake`).
- `SitesContext`: site properties and distances as one-dimensional
  read-only NumPy arrays; scalars are broadcast to all sites.
- `GroundMotionResult.mean`: **natural logarithm of the median**
  intensity measure, *not* the linear median.
- `GroundMotionResult.sigma`: total aleatory standard deviation in
  natural-log space, or `None` when not supplied.
- `GroundMotionResult.tau` and `phi`: inter-event and intra-event
  natural-log standard deviations when implemented; do not infer them
  from total sigma if they are absent.

Within ShakeLab's current ground-motion conventions, **PGA and SA
are in g**, and **PGV is in m/s**. Each GMPE is responsible for
converting any original coefficient or publication units to these
outputs. The median in linear space is `numpy.exp(result.mean)`;
the approximate `median * exp(±sigma)` bounds are appropriate only
when total lognormal dispersion is available and interpreted as such.
Do not confuse a residual sigma with epistemic model uncertainty,
rupture-scaling uncertainty, or a confidence interval on the median.

`MAGNITUDE_TYPE` documents a model's intended magnitude scale. This
low-level API expects the caller to supply the **appropriate magnitude**.
The higher-level provider currently passes a numerical event magnitude
through without converting differing magnitude labels; **that is a
temporary operational assumption, not a physical equivalence of Mw and
ML**. The implementation does not automatically enforce published
magnitude/distance validity ranges for every GMPE. Users must check the
original reference and the specific implementation before extrapolating.

Individual model files describe significant implementation-specific
choices and differences from their published or OpenQuake counterparts.
The registry denotes implemented code availability, **not scientific
validation or universal applicability**.

## 2. Available models and variants

All identifiers below come from the package's current
[`registry.json`](registry.json). Use
`list_gmpes()` to obtain the authoritative installed list; aliases
are accepted through `resolve_name()`. Regional and distance
variants are **distinct registered classes**, not interchangeable
names of a single coefficient set.

### Model families

| Family / implementation | Registered identifiers or variants | Required distance(s) | Required source parameters in addition to `mag` | Required site parameters |
| --- | --- | --- | --- | --- |
| Abrahamson et al. (2014), ASK14 | `AbrahamsonEtAl2014`, `AbrahamsonEtAl2014RegCHN`, `AbrahamsonEtAl2014RegJPN`, `AbrahamsonEtAl2014RegTWN` | `rrup`, `rjb`, `rx`, `ry0` | `rake`, `dip`, `ztor`, `width` | `vs30`, `z1pt0`, `vs30measured` |
| Akkar et al. (2014) | `AkkarEtAlRepi2014`, `AkkarEtAlRhypo2014`, `AkkarEtAlRjb2014` | Respectively `repi`, `rhypo`, `rjb` | `rake` | `vs30` |
| Atkinson (2015), induced-seismicity model | `Atkinson2015`, `Atkinson2015AltDistSat` | `rhypo` | None | None |
| Bindi et al. (2011) | `BindiEtAl2011` | `rjb` | `rake` | `vs30` |
| Bindi et al. (2014) | `BindiEtAl2014Rhyp`, `BindiEtAl2014Rjb` | Respectively `rhypo`, `rjb` | `rake` | `vs30` |
| Boore et al. (2014), BSSA14 | `BooreEtAl2014` | `rjb` | `rake` | `vs30` |
| Bragato & Slejko (2005) | `BragatoSlejko2005` | `repi` | None | None |
| Campbell & Bozorgnia (2014), CB14 | `CampbellBozorgnia2014`, `CampbellBozorgnia2014HighQ`, `CampbellBozorgnia2014LowQ`, `CampbellBozorgnia2014JapanSite`, `CampbellBozorgnia2014HighQJapanSite`, `CampbellBozorgnia2014LowQJapanSite` | `rrup`, `rjb`, `rx` | `rake`, `dip`, `ztor`, `width`, `hypo_depth` | `vs30`, `z2pt5` |
| Chiou & Youngs (2014), CY14 | `ChiouYoungs2014`, `ChiouYoungs2014Italy`, `ChiouYoungs2014Japan`, `ChiouYoungs2014Wenchuan`, `ChiouYoungs2014NearFaultEffect` | `rrup`, `rjb`, `rx`; near-fault variant also `rcdpp` | `rake`, `dip`, `ztor` | `vs30`, `vs30measured`, `z1pt0` |
| Kotha et al. (2020), ESHM20 adaptation | `KothaEtAl2020ESHM20` | `rjb` | `hypo_depth` | `vs30`, `vs30measured`, `eshm20_region` |
| Lanzano et al. (2019) | `LanzanoEtAl2019Rjb`, `LanzanoEtAl2019Rrup` | Respectively `rjb`, `rrup` | `rake` | `vs30` |
| Massa et al. (2008) | `MassaEtAl2008Ml`, `MassaEtAl2008Mw` | `repi` | None | `ec8_class` |

The table reflects the **declared context requirements** in the
current Python classes, not a claim that these are all the physical
inputs needed for responsible scientific use. In particular,
`rcdpp` is declared as a distance input in the CY14 near-fault
implementation; model-specific quantities must be examined in its
source. Some GMPEs also support optional inputs or constructor
options. The definitions in the class remain authoritative.

The core registry does **not** include the `BragatoSlejko2005JB`
class present in the source module: it is a separate implementation
but is not selectable by its name through the packaged registry.
Adding it to the registry would be a distinct code change.

### Magnitude and applicability notes

The `MAGNITUDE_TYPE` declarations are `Mw` for the registered
families other than `BragatoSlejko2005` (`ML`) and
`MassaEtAl2008Ml` (`ML`); `MassaEtAl2008Mw` uses `Mw`.
This is metadata, not a conversion routine.

Useful source-specific calibration statements already documented in
the implementations include:

- **Bragato & Slejko (2005):** Eastern Alps, local magnitude,
  approximately ML 2.5–6.3, with epicentral distance for the
  registered version.
- **Massa et al. (2008):** Northern Italy; the class documentation
  identifies epicentral distances up to 100 km, ML 3.5–6.3 and
  Mw up to 6.5 for the corresponding variants.
- **Kotha et al. (2020), ESHM20:** shallow European crustal
  application; `eshm20_region` is an integer 0–5 (0 is generic).
- **ASK14, BSSA14, CB14 and CY14:** NGA-West2 family implementations
  with nontrivial source, site, regional and near-fault terms; consult
  each class and its original paper before selecting a variant.
- **Atkinson (2015):** induced-seismicity formulation; it should not
  be treated as an unrestricted replacement for tectonic-earthquake
  models.

**A complete, uniformly verified table of magnitude, distance and
period calibration limits is not currently encoded in the registry.**
Do not infer one from the list of available IMTs. Prefer the original
publication, the model-specific class documentation, and a dedicated
scientific validation suite when assessing extrapolation.

### IMT discovery: do not guess coefficient keys

An exact IMT such as `SA-0.20` may be stored under a different
textual key from `SA-0.2`. The low-level `get_coefficients()`
performs **exact string matching**. Inspect `list_imts()` for each
model and coefficient set rather than assuming interpolation or
automatic period normalization.

```python
from shakelab.groundmotion.gmpe.registry import (
    create_gmpe, list_gmpes, resolve_name,
)

print(list_gmpes())
print(resolve_name("ASK14"))  # "AbrahamsonEtAl2014"

gmpe = create_gmpe("AbrahamsonEtAl2014")
print(gmpe.list_imts()[:10])
print(gmpe.required_rupture_parameters)
print(gmpe.required_site_parameters)
print(gmpe.required_distances)
```

Models can expose additional nonstandard measures such as IA or
Fourier-amplitude keys; the registry and coefficients decide what
`list_imts()` reports. Those measures do not automatically have
the same units as PGA, PGV or SA.

## 3. Evaluating a GMPE directly

The low-level API deliberately depends on **neither** an earthquake
catalogue object **nor** a geometry or application-specific provider.
For example, this computes PGAs for one synthetic Bragato–Slejko
source at several epicentral distances:

```python
import numpy as np

from shakelab.groundmotion.gmpe.base import (
    RuptureContext, SitesContext,
)
from shakelab.groundmotion.gmpe.registry import create_gmpe

gmpe = create_gmpe("BragatoSlejko2005")

rupture = RuptureContext(mag=5.5)  # Local magnitude (ML)
sites = SitesContext(repi=[5.0, 10.0, 20.0, 50.0])

result = gmpe.ground_motion("PGA", rupture, sites)

median_pga_g = np.exp(result.mean)
sigma_ln = result.sigma

print("PGA median (g):", median_pga_g)
print("Total sigma, ln(PGA):", sigma_ln)
```

The GMPE itself receives distances and site fields that have already
been defined. `SitesContext(size=n, vs30=500.0, rjb=[...])`
broadcasts a scalar `vs30` to `n` evaluation sites and rejects
incompatible array lengths. Inputs are stored as read-only NumPy
arrays. `RuptureContext` and `SitesContext` support attribute access
and mapping access.

### Attenuation curves and response spectra

Evaluate a vector of distances in one call, then plot
`np.exp(result.mean)` against the **distance measure required by
the selected model**, rather than silently substituting another
distance definition. For spectra, iterate over the **exact available
SA keys** and parse their periods for the x-axis; do not interpolate
between periods unless that interpolation has been intentionally
implemented and scientifically justified. Periods and damping must
be consistent with the selected model's definition.

The general-purpose graphical demonstration
[`../groundmotion_example.py`](../groundmotion_example.py) compares
native, point and finite-source provider calculations, and displays
rupture geometry and actual GMPE distance inputs. Matplotlib is
optional and not required for core GMPE evaluation.

## 4. The higher-level GroundMotionProvider

The `GroundMotionProvider` adapts a `ScenarioEvent` and WGS84
site locations to the numerical contexts required by an individual
GMPE. It can derive elementary point-source distances, construct
finite rectangular rupture distances when requested, merge source
overrides and optionally collect trace records. These capabilities
are general ShakeLab APIs, **not** obligations of a GMPE class.

```python
from shakelab.groundmotion.providers import (
    GroundMotionProvider, ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint

event = ScenarioEvent.from_lonlat_depth(
    5.5, 13.0, 46.0, 8.0, magnitude_type="ML",
)
sites = [WgsPoint(13.1, 46.1), WgsPoint(13.2, 46.2)]

provider = GroundMotionProvider.gmpe(
    "BragatoSlejko2005",
    config={"source_geometry": "native"},
)
median_g, sigma_ln = provider.evaluate("PGA", sites, event)
```

Unlike the low-level `GMPE.ground_motion()` result, this provider
returns the **linear-space median** and natural-log sigma. With
`source_geometry="native"`, `repi` and `rhypo` can be derived
when needed; `rjb`, `rrup`, `rx` and `ry0` generally need
explicit values unless the user selects a point approximation or
finite rupture. Source parameters follow
**provider < event < runtime** precedence; explicitly passed
distances override computed distances. See
[`../RUPTURE.md`](../RUPTURE.md) for examples of each mode.

If a GMPE has required site fields such as `z1pt0`,
`vs30measured`, `z2pt5` or `ec8_class`, they must be supplied
according to its implementation; a geometry selection does not fill
missing site properties.

## 5. Implementing a new GMPE

New scientific formulations belong in
`shakelab/groundmotion/gmpe/<model_name>.py`; their coefficients
belong under `gmpe/data/`. Keep imports of engineering applications,
server logic, scenario models, and ObsPy out of the numerical model.
Follow the repository's Python style and documented public APIs.

### Step 1 — Identify the scientific specification

Record the exact publication, equation numbers, coefficient table,
horizontal component definition, IMTs, damping, magnitude scale,
distance definition, site classification, units, applicability
bounds, regional variants and aleatory variability decomposition.
Identify which features are implemented, deliberately omitted or
approximated. Check licensing and attribution for any reused
coefficients.

Do not invent sigma, `tau`/`phi` or calibration limits if the
publication does not supply enough information to support them.

### Step 2 — Declare the context contract

Define a subclass of `GMPE` with class-level declarations:

- `REQUIRES_RUPTURE_PARAMETERS` and
  `OPTIONAL_RUPTURE_PARAMETERS`;
- `REQUIRES_SITE_PARAMETERS` and
  `OPTIONAL_SITE_PARAMETERS`;
- `REQUIRES_DISTANCES` and `OPTIONAL_DISTANCES`;
- `MAGNITUDE_TYPE` and `REFERENCE_VELOCITY`, where appropriate;
- `_COEFF_FILE` and `_COEFF_SET` for the packaged coefficients.

Names must be non-empty strings with no overlap between conflicting
categories. The base class validates declarations and can validate
the required fields of both contexts.

### Step 3 — Provide the coefficients

The loader in `GMPE.import_coeff_from_json()` expects this
**existing ShakeLab schema**:

```json
{
  "coefficients": {
    "default": {
      "keys": ["a", "b", "sigma"],
      "type": {
        "PGA": [1.0, -1.0, 0.5],
        "SA-0.20": [1.2, -1.1, 0.6]
      }
    }
  }
}
```

These are **invented demonstration values, not scientific
coefficients**. Each `type` row follows the same order as `keys`.
The file is loaded from `gmpe/data/` relative to `base.py`,
and entries are converted to floats. The implementation does not
automatically verify physical consistency, regression-space units,
that every row has the same length as `keys`, or whether every
published period is present. Add explicit regression checks for
those properties in the new model's tests.

For variants with different coefficient sets, use a distinct
`_COEFF_SET` referring to another object under `coefficients`.
The existing Akkar and Massa implementations are useful patterns.

### Step 4 — Implement the numerical method

Below is a minimal **structural example**; the coefficients and
functional form are placeholders and must **not** be used as a
scientific prediction equation.

```python
import numpy as np

from shakelab.groundmotion.gmpe.base import (
    GMPE, GroundMotionResult, RuptureContext, SitesContext,
)


class ExampleGMPE(GMPE):
    """Illustrative API shape only; not a calibrated GMPE."""

    MAGNITUDE_TYPE = "Mw"
    REQUIRES_RUPTURE_PARAMETERS = {"mag"}
    REQUIRES_SITE_PARAMETERS = {"vs30"}
    REQUIRES_DISTANCES = {"rjb"}
    _COEFF_FILE = "example_gmpe.json"
    _COEFF_SET = "default"

    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        self.validate_contexts(rupture, sites)
        coeff = self.get_coefficients(imt)

        magnitude = float(rupture.mag)
        rjb = np.asarray(sites.rjb, dtype=float)
        vs30 = np.asarray(sites.vs30, dtype=float)

        if not np.isfinite(magnitude):
            raise ValueError("Magnitude must be finite.")
        if not np.isfinite(rjb).all() or (rjb < 0).any():
            raise ValueError("rjb must be finite and non-negative.")
        if not np.isfinite(vs30).all() or (vs30 <= 0).any():
            raise ValueError("vs30 must be finite and positive.")

        # DEMONSTRATION ONLY: not a scientifically calibrated equation.
        mean_ln = (
            coeff["a"]
            + magnitude
            + coeff["b"] * np.log(rjb + 1.0)
            + 0.1 * np.log(vs30 / 760.0)
        )
        sigma_ln = np.full(sites.size, coeff["sigma"])

        return GroundMotionResult(
            mean=mean_ln,
            sigma=sigma_ln,
        )
```

This illustrates vectorization, validation, and the return contract
only. In a real model, all terms, clipping and saturation, distance
dependencies, unit conversions, variability and special cases must
follow the scientific reference exactly. Published base-10 log
regressions require conversion to natural-log space, including their
sigmas; PGV published in cm/s must be converted to m/s when that is
the intended output. Do not double-convert outputs.

If a subclass supplies its own `__init__`, call `super().__init__()`
once after setting any instance options used during coefficient
loading. If `tau` and `phi` are available, return them separately
and verify that they are consistent with the total sigma for the
particular model's specified decomposition.

### Step 5 — Register the class and test discovery

Add an entry to `gmpe/registry.json`:

```json
{
  "ExampleGMPE": {
    "pointer": "shakelab.groundmotion.gmpe.example_gmpe:ExampleGMPE",
    "alias": ["ExampleModel"]
  }
}
```

This is an **entry inside** the existing
`{"gmpemap": {...}}` object, not a replacement for the entire
registry. The module path, class name and aliases must be valid.
The registry uses explicit, lazy import; adding a Python file alone
does **not** register the new model. Keep alias names unique.

Test discovery and construction:

```python
from shakelab.groundmotion.gmpe.registry import (
    create_gmpe, list_gmpes, resolve_name,
)

assert "ExampleGMPE" in list_gmpes()
assert resolve_name("ExampleModel") == "ExampleGMPE"

gmpe = create_gmpe("ExampleGMPE")
assert "PGA" in gmpe.list_imts()
```

The environment variable `SHAKELAB_GMPE_REGISTRY` can override the
packaged registry with a JSON file at an **absolute path**. This is
useful for controlled development, but test the final packaged
registry before release. `clear_caches()` resets registry and class
lookups in tests.

### Step 6 — Numerical and scientific validation

Place pytest regression tests under `tests/groundmotion/` (or a
dedicated GMPE test subdirectory when the suite grows). At minimum:

1. **Published reference points:** fixed source/site/IMT examples
   matched to independent reference tables or a validated external
   implementation; specify tolerances and any known differences.
2. **IMTs and coefficient integrity:** exact available periods,
   coefficient row lengths, finite values and correct coefficient set.
3. **Units and log base:** PGA/SA in g, PGV in m/s, `mean` as ln(IM),
   `sigma`, `tau`, `phi` in natural-log space.
4. **Input validation:** missing or overlapping declared context
   fields, incompatible site-array lengths, invalid negative
   distances, nonfinite values and model-specific domain constraints.
5. **Vectorization:** compare multi-site evaluations against
   equivalent one-site calls; ensure input arrays are not mutated.
6. **Rupture/site dependence:** test expected responses to magnitude,
   distance, mechanism, site conditions, regions and near-fault
   effects without presuming monotonicity where the model does not
   require it.
7. **Provider consistency:** compare direct `ground_motion()`
   results with `GroundMotionProvider.gmpe(...).evaluate()`
   after matching the **same actual** context fields and distances.
8. **Applicability and provenance:** document reference and
   implementation limits, extrapolation decisions and any missing
   terms or unsupported variants.

For example, a basic direct-model vectorization regression:

```python
import numpy as np

from shakelab.groundmotion.gmpe.base import RuptureContext, SitesContext
from shakelab.groundmotion.gmpe.registry import create_gmpe

gmpe = create_gmpe("BragatoSlejko2005")
rupture = RuptureContext(mag=5.5)
distances = [5.0, 10.0, 20.0]

all_sites = gmpe.ground_motion(
    "PGA", rupture, SitesContext(repi=distances),
)
for index, distance in enumerate(distances):
    single = gmpe.ground_motion(
        "PGA", rupture, SitesContext(repi=[distance]),
    )
    np.testing.assert_allclose(
        np.asarray(all_sites.mean)[index],
        np.asarray(single.mean)[0],
        rtol=1e-12, atol=1e-12,
    )
```

A passing smoke test or successful provider evaluation is **not**
sufficient scientific validation. Preserve the published reference
data and expected numerical results in regression tests rather than
checking only that outputs are positive or finite.

## 6. Developer checklist

Before submitting a new GMPE:

- [ ] Original reference, applicable tectonic regime, magnitude scale,
      distance definition, component, units and ranges documented.
- [ ] Exact coefficient source and sets reviewed and attributed.
- [ ] Required and optional rupture/site/distance fields declared.
- [ ] All outputs converted into ShakeLab's logarithmic/physical units.
- [ ] Sigma availability and decomposition handled correctly.
- [ ] Registry entry and aliases added without collisions.
- [ ] Independent scientific benchmarks and edge cases tested.
- [ ] Standalone and provider-level use tested.
- [ ] Example and documentation updated without binding core GMPE
      code to a particular application.

The guidance above describes the current API; it intentionally does
not prescribe how any particular application must save results,
choose a rupture surface, or assemble hazard/impact workflows.
