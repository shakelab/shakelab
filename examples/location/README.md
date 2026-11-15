# Location examples and tests

## Examples

These scripts show how to read seismic picks, locate events and calculate
travel times in flat layers. For scripts that use the bundled Socchieve
data, run from this directory or set the data path in the script.

| Example | Task |
| --- | --- |
| [01](Examples/01_picks_and_stations.py) | Read station coordinates, the PUN origin and OBS picks. Try DAT, QuakeML and HPL through the same reader. |
| [02](Examples/02_locate_homogeneous.py) | Relocate Socchieve with constant P and S velocities. |
| [03](Examples/03_locate_layered.py) | Relocate the same event using the supplied 1D model. |
| [04](Examples/04_locate_relative.py) | Recover a synthetic hypocenter from inter-station time differences. |
| [05](Examples/05_estimate_origin_time.py) | Estimate a weighted origin time for a known hypocenter. |
| [06](Examples/06_first_arrival_curves.py) | Plot the crossover from direct to critically refracted first arrivals. |
| [07](Examples/07_reflected_ray.py) | Trace a same-phase reflection from a chosen interface. |
| [08](Examples/08_compare_socchieve.py) | Compare P and S residuals at the catalog origin and after relocation. Uncomment the `fig.savefig` line to save the figure. |
| [09](Examples/09_velocity_profile.py) | Plot the Friuli P and S velocities as functions of depth. |
| [10](Examples/10_reference_profiles.py) | Browse the built-in profiles and calculate P travel times without input files. |

Install ShakeLab from the checkout, then run the examples:

```sh
python3 -m pip install -e '.[plot]'
cd examples/location
python3 01_picks_and_stations.py
python3 10_reference_profiles.py
```

The `plot` extra provides Matplotlib for examples with figures.

## Reference velocity profiles

`shakelab.location.profiles` provides offline access to published layer
values. Use `list_models()` to discover identifiers, `load_model(name)`
to create a `LayeredModel1D`, and `model_info(name)` to inspect the
source, region and depth reference. The loaded model exposes layer tops
as `z_top` and velocities as `vp` and `vs`, in metres and metres per
second. The catalog is discovered from the JSON files in
[`shakelab/location/models/`](../../shakelab/location/models/).

| Identifier | Source and scope |
| --- | --- |
| `ogs-friuli` | Friuli Venezia Giulia bulletin model documented by [Saraò et al. (2021), §3.2](https://doi.org/10.5194/essd-13-2245-2021). Vs is derived from Vp with the published ratio of 1.78. |
| `iasp91-crust-halfspace` | Flat crustal approximation to [Kennett and Engdahl (1991)](https://doi.org/10.1111/j.1365-246X.1991.tb06724.x). The velocities just below the 35 km Moho are held constant in the half-space. |
| `ingv-rsn1` | Historical Italian reference profile from [Peruzza et al. (2021), Table 2](https://doi.org/10.5194/se-12-2021-2021). Vs is derived from Vp with the published ratio of 1.73. |

The [AlpsLocPS dataset](https://doi.org/10.35097/1965) is available
separately under CC BY-NC-SA 4.0.

### Read a user-supplied model

```python
from shakelab.location.profiles import read_model

model = read_model("my-model.csv")
```

Use `read_model("my-model.json")` for JSON or pass `format="csv"` or
`format="json"` if the extension differs. List layers from the surface
down, with thickness in metres and velocities in metres per second. A
half-space can appear only as the last layer. It is written as `inf` in
CSV and `null` in JSON. A finite last layer is also supported.

```csv
thickness (m),vp (m/s),vs (m/s)
2000,4000,2200
inf,6000,3200
```

```json
{
  "layers": [
    {"thickness_m": 2000, "vp_m_s": 4000, "vs_m_s": 2200},
    {"thickness_m": null, "vp_m_s": 6000, "vs_m_s": 3200}
  ]
}
```

The CSV density column is optional for `read_model()`. Use the bundled
files as examples of the optional JSON `metadata` object.

## Tests and data

From the checkout root, run `python3 -m pytest -q tests/location`.
The Socchieve tests and examples share the files under `data/`. These
include the origin and pick exports for the local earthquake of
27 March 2024, the station list and the experimental 20-layer
velocity profile used in the comparison scripts.
