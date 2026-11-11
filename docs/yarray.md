# YArray: seismic-array processing

`shakelab.exploration.yarray` measures the waves crossing a seismic
array by frequency-wavenumber (f-k) analysis: phase velocity and
propagation azimuth as functions of frequency and time, and the signed
ellipticity of Rayleigh waves. It also tells, before or after a
survey, what an array can measure: its geometric response, the
statistical performance of the estimators and the frequency-velocity
domain where the analysis is reliable.

The same components are used in three ways:

- as a **library** from Python (`beamform`, `Beamformer`, results,
  statistics, planning);
- from the **command line** (`slab-yarray`);
- in **real time**, as a function of the SeedLink client
  (`ArrayMonitor`).

This guide gives the theory (Part I), the use (Part II) and the
limits (Part III). The modules, the design choices and the state of
the work are in the [package README](../shakelab/exploration/yarray/README.md);
runnable scripts of every way of use are in
[`examples/exploration/`](../examples/README.md#exploration). The guide
derives from the four documents of the YArray development code
(August 2026: theoretical background, software architecture, user
manual, array response theory), updated to the implementation in
ShakeLab.

## Contents

Part I. Theory

1. [Fourier coefficients and cross-spectral matrices](#1-fourier-coefficients-and-cross-spectral-matrices)
2. [Plane waves and steering vectors](#2-plane-waves-and-steering-vectors)
3. [Estimators](#3-estimators)
4. [Components and Rayleigh three-component beamforming](#4-components-and-rayleigh-three-component-beamforming)
5. [Uncertainty of the power](#5-uncertainty-of-the-power)
6. [Spectral maxima and their density](#6-spectral-maxima-and-their-density)
7. [Geometric array response](#7-geometric-array-response)
8. [Statistical performance and operational limits](#8-statistical-performance-and-operational-limits)

Part II. Use

9. [Data](#9-data)
10. [Configuration](#10-configuration)
11. [Processing and results](#11-processing-and-results)
12. [Command line](#12-command-line)
13. [Real time](#13-real-time)
14. [Survey planning](#14-survey-planning)
15. [Practical guidance and processing times](#15-practical-guidance-and-processing-times)

Part III. Limits

16. [Current limitations and planned extensions](#16-current-limitations-and-planned-extensions)
17. [References](#17-references)

# Part I. Theory

The wavefield is assumed locally plane over the array and propagating
horizontally. Coordinates are east (x), north (y) and up, in metres
(or in any distance unit, consistently with velocities and
wavenumbers); azimuths are clockwise from north and give the direction
towards which the wave travels.

## 1. Fourier coefficients and cross-spectral matrices

Each record is cut into elementary **windows**, tapered (Tukey) and
transformed at the analysis frequencies. The windows have a fixed
duration, or a duration of `n_c` cycles of each frequency,
`T_w = n_c / f`; they may overlap. Consecutive windows form
**covariance blocks**, which may also overlap: each block gives one
estimate of the cross-spectral matrix (CSM) and therefore one
spectrum, so that a long record gives a time series of results.

For window `m` the coefficients of the `N` channels at frequency `f`
form a vector `x_m(f)`. The CSM of a block of `L` windows is

```
C(f) = (1 / L) sum_m x_m(f) x_m(f)^H
```

The coefficients are scaled as a one-sided power spectral density,

```
x(f) = sqrt(2 delta / sum_n w_n**2) sum_n w_n x_n exp(-2 pi i f n delta)
```

(`w` taper, `delta` sampling interval), so that the diagonal of `C` is
the PSD of each channel and Bartlett and Capon give the density of a
plane wave in the units of the data squared per Hz. Data not corrected
for the instrument work the same way: the power is then in counts
squared per Hz, or `(unit)**2/Hz` when the units are unknown.

Two transforms give the coefficients
(`shakelab.signals.crossspectral.coefficients`):

- **direct**: the transform evaluated at each requested frequency; any
  frequency, fixed or cycle windows;
- **RFFT**: one real FFT per window (fixed windows only), at the
  **nearest** bin (repeated bins are processed once, so the
  frequencies are the bins) or **interpolated** at the requested
  frequencies (zero padding and local cubic interpolation of the
  complex coefficients). Zero padding densifies the sampled transform
  but does not improve the resolution of the window.

**Effective number of windows.** Overlapping tapered windows are
correlated, so `L` windows carry the information of fewer independent
ones (Percival and Walden, 1993, Section 6.17):

```
L_eff = L / (1 + 2 sum_q (1 - q / L) rho_q**2)
```

with `rho_q` the overlap of the taper with itself shifted by `q`
windows. `L_eff` is used by the standard errors of the power.

## 2. Plane waves and steering vectors

A plane wave of angular wavenumber `k` and azimuth `theta` has the
wavevector

```
kx = k sin(theta),   ky = k cos(theta),   k = 2 pi f / c = 2 pi f s
```

(`c` phase velocity, `s` slowness), and at station `i` the phase factor

```
a_i(k) = exp(-i (kx x_i + ky y_i))
```

The steering vectors `a` are evaluated on a polar **search grid**,
whose radial coordinate is velocity, slowness or wavenumber, sampled
linearly or logarithmically, with uniform or explicit azimuths.

## 3. Estimators

**Bartlett** (conventional beamforming):

```
P_B = a^H C a / N**2
```

robust and inexpensive; its resolution is that of the array response
(Section 7).

**Capon** (minimum variance, Capon 1969), with trace-scaled diagonal
loading `lambda` for numerical stability:

```
R = C + lambda (tr(C) / N) I,      P_C = 1 / (a^H R^-1 a)
```

It suppresses the power not coherent with the tested wave and resolves
closer waves than Bartlett. A CSM of `L` windows underestimates the
Capon power by about `(L - N + 1) / L` (Capon and Goodman, 1970): the
bias does not move the maxima. With `strict=True` (default) Capon and
MUSIC require more windows per block than the dimension of the matrix.

**MUSIC** (Schmidt, 1986): with `E_n` the eigenvectors of the noise
subspace (all but `music_sources`),

```
P_M = 1 / (a^H E_n E_n^H a)
```

a pseudo-spectrum without physical units: very sharp when the
subspace model holds, sensitive to the assumed number of sources.

## 4. Components and Rayleigh three-component beamforming

For each tested azimuth the horizontal motion is projected on the
radial and transverse directions,

```
R = E sin(theta) + N cos(theta),    T = E cos(theta) - N sin(theta)
```

(transverse: radial rotated by +90 degrees, clockwise). The component
modes are: `vertical`, `radial`, `transverse`,
`three_component_separate` (the three spectra on the same grid) and
`rayleigh_joint`.

**Rayleigh three-component beamforming** (RTBF, Wathelet et al.,
2018) combines, for each azimuth, the radial and vertical coefficients
in one matrix of dimension `2N`, steered by

```
E = [-i e a, a]
```

where `e` is the signed ellipticity (radial over vertical amplitude,
positive for radial motion `e sin(wt)` with vertical `cos(wt)`). The
Capon power `P` of `E` is computed for every trial ellipticity and the
objective `e**2 P**2` is maximized: the spectrum is that maximum, and
the ellipticity at the maximum is kept with each spectral maximum. The
ellipticity is sampled uniformly in its angle `chi = arctan(e)`
strictly inside (-90, 90) degrees (default -85 to 85, 61 values), so
that large ellipticities are reached without an arbitrary bound on
`e`; zero is excluded. Its power has the square of the units of the
power.

Wathelet (2024) shows that the 2018 formulation is distorted by
incoherent noise for vanishing or very large ellipticities and
proposes a normalized steering in the angular ellipticity,
`[-i sin(xi) a, cos(xi) a]`, with the Capon power maximized directly
over `xi` (ARTBF, and ARDS on the raw three components), a fast Newton
search of `xi` and the correction of the bias due to the difference of
incoherent noise between vertical and horizontal components. These are
not implemented yet (Section 16).

## 5. Uncertainty of the power

`statistical_uncertainty="empirical"` gives the standard error of the
power at every grid point (single components):

- Bartlett: the power of each window is computed separately and the
  variance of their mean is divided by `L_eff`;
- Capon: first-order propagation (delta method) of the variance of the
  windows through `dP = P**2 b^H dR b`, `b = R^-1 a`.

`geometric_uncertainty="expected_coherence"` (Bartlett) gives the
power expected over independent Gaussian errors of the station
coordinates, with covariance `S_i`: the steering coefficients are
damped by the expected coherence `g_i = exp(-k^T S_i k / 2)` and the
diagonal is kept exactly,

```
P = (a g)^H C (a g) / N**2 + sum_i (1 - g_i**2) C_ii / N**2
```

These are uncertainties of the power, not of the position of a
maximum.

## 6. Spectral maxima and their density

At each frequency, block and component the local maxima of the
spectrum (azimuth x radial neighbourhood, default 3 x 3, optional
minimum power relative to the largest) are ranked by power; rank 1 is
the strongest, not a mode. Each maximum keeps frequency, block and its
time (centre of the block), wavenumber and azimuth (velocity,
slowness, `kx`, `ky` derived), power, its standard error and, for
RTBF, the ellipticity.

The maxima of many blocks (a long record, or the analyses of the
monitor joined) scatter around the dispersion curves. Their density at
each frequency is estimated as a weighted mixture of Gaussian kernels
in wavenumber, log-slowness or log-velocity:

```
D(x | f) = sum_i w_i N(x; x_i, sigma_i)
```

with uniform, power or relative-power weights and kernels as wide as
the resolution of the search grid (default), Silverman's rule or a
fixed width. Its summary gives mode, median and a quantile interval;
its local maxima, with their prominence, are candidate modes at each
frequency (not yet connected across frequencies).

## 7. Geometric array response

The response of an array depends only on the station positions. It
tells which wavevectors the array resolves and which it confuses,
independently of the data and of the estimator.

**Array factor and ambiguity.** For stations at `r_n`,

```
B(dk) = (1 / N) sum_n exp(-i dk . r_n),     A(dk) = |B(dk)|**2
```

`A(0) = 1` and `A` does not change with a translation of the array.
With normalized steering vectors `a`, `|a(k1)^H a(k2)|**2 = A(k2 - k1)`:
`A` is the similarity of two wavevectors as a function of their
difference. Its central region describes the **resolution**; remote
maxima, or connected regions of high response, are **ambiguities**
(aliases, grating lobes of regular arrays, strong sidelobes of
irregular ones). For an irregular array both depend on the direction,
so the full two-dimensional response is the primary result; the
classical rules `2 pi / d_max` (Rayleigh), `2 pi / d_min` (alias
period) and `pi / d_min` (Nyquist) are diagnostic proxies only.

**Main lobe.** For a level `T_R` (default 0.5) the main lobe is the
connected region of `A >= T_R` containing the origin; its radius
`r_R(theta)` is the first crossing of the level along each axis
azimuth, `2 r_R` the full width. The conservative scalar is the
largest radius, `resolution_radius`: below it two waves are not
resolved (a single wave can still be measured at smaller `k`, with a
precision limited by noise and averaging; Section 8). Near the origin
`A(k) ~ 1 - k^T C_r k`, with `C_r` the covariance of the centred
coordinates: its principal directions and the anisotropy ratio
describe the geometry.

**Ambiguity.** For a level `T_A` (default 0.5) the regions of
`A >= T_A` other than the central one are the remote ambiguity
regions; their local maxima are the ambiguity peaks (displacements,
not propagation directions). Two distances must not be confused: the
**peak distance**, from the origin to the nearest significant peak,
and the **onset distance**, from the origin to the nearest point of a
remote region; the onset is shorter, because a region crosses the
level before its peak. The peak sidelobe level and the ambiguity
margin `1 - max A_remote` (a squared correlation, not a probability)
summarize the strength of the ambiguities.

**Search domains.** A search over a domain `D` of wavevectors is free
of ambiguity at level `T_A` if no difference `k2 - k1` of two of its
wavevectors falls into a remote region. For a full-azimuth disk
`|k| <= K` the differences fill the disk `|dk| <= 2K`, whence the
conservative limit

```
aliasing_radius = onset_distance / 2
```

the rigorous form of the factor one half of the classical aliasing
rules (for a regular 3 x 3 array of 10 m: peak at `2 pi / 10 =
0.628` rad/m, onset 0.531 rad/m, aliasing radius 0.266 rad/m; the
limit tends to `pi / d` as `T_A` tends to one).

**Velocity limits.** At frequency `f` the wavenumber limits become

```
c_max(f) = 2 pi f / resolution_radius,   c_min(f) = 2 pi f / aliasing_radius
```

geometric bounds at the chosen levels, not confidence limits. When
the aliasing radius is smaller than the resolution radius no
wavenumber is both resolved and unambiguous (e.g. three stations).

**Coordinate uncertainty.** With independent Gaussian errors of
covariance `S_n` and `q_n(k) = exp(-k^T S_n k / 2)`, the expected
response is

```
E[A(k)] = (|sum_n q_n exp(-i k . r_n)|**2 + sum_n (1 - q_n**2)) / N**2
```

where the second term keeps the auto terms, in which the same error
appears twice and cancels. Uncertain positions damp the coherent
structure increasingly with the wavenumber.

A linear array has a valid response but no closed two-dimensional
main lobe (it does not resolve the wavevector across the line): the
limits are refused for it. The numerical grid of the response is
chosen from the proxies only for sampling, and an excessive grid is an
error, never silently coarsened.

## 8. Statistical performance and operational limits

The geometric response is the point-spread function of Bartlett for a
single coherent wave; it is not that of Capon, MUSIC or RTBF, whose
behaviour depends on the SNR, on the number of windows, on the
loading and on the subspace model. `A = 0.8` is not a probability of
aliasing. The performance of an estimator is therefore measured by
Monte Carlo simulation with the estimators themselves:

```
x_m = s_m a(k0) + n_m,    s_m ~ CN(0, P_s),    n_m ~ CN(0, P_n I)
```

(independent realizations, SNR `P_s / P_n` per station). Each trial
builds a sample CSM from `L` realizations, computes the spectrum on
the search grid and classifies its maximum by the error
`dk = k_hat - k0` with the geometric response:

| Outcome | Meaning |
|---|---|
| `BOUNDARY` | maximum on the radial boundary of the search grid |
| `LOCAL` | `dk` in the main lobe (geometrically local, not "accurate") |
| `ALIAS` | `dk` in a remote ambiguity region (an alias family) |
| `OTHER` | anywhere else (noise-driven maxima) |

The probabilities (with Wilson score intervals), the bias and
covariance of the local errors (radial and transverse) and the exact
velocity statistics (each estimate converted, `c = 2 pi f / k`)
quantify the performance; a sweep over SNR and windows tells how much
data a survey needs.

The **operational assessment** repeats the experiment over a grid of
frequencies, velocities and azimuths and marks a point usable when it
meets explicit criteria (minimum LOCAL probability, maximum ALIAS,
OTHER and BOUNDARY probabilities, maximum relative velocity RMSE,
grid quantization small relative to the resolution radius); otherwise
the first failed criterion is its limiting factor (search boundary,
grid quantization, aliasing, other failures, low SNR, resolution).
There is no opaque quality score. Recommendations explain the limits
and suggest processing changes; the geometry itself is not optimized.

# Part II. Use

```python
from shakelab.exploration.yarray import ArrayData, Beamformer
```

The package exposes the main classes and functions; the modules are
`data`, `config`, `estimators`, `beamforming`, `statistics`,
`response`, `performance`, `operational`, `io`, `monitor` and `cli`,
the figures `shakelab.graphics.yarray`. Options with fixed values are
enumerations whose string values are accepted everywhere
(`method="capon"`).

## 9. Data

An `ArrayData` holds the stations (`StationData`: records, ENU
coordinates, orientation of the components, optional coordinate
covariance and geographic position). It is filled in one go or
progressively, from any source:

| Source | Use |
|---|---|
| Array file (JSON: local coordinates, one waveform file per component, reader options) | `ArrayData("array.json")`, or `array.read("array.json")` |
| Waveforms read by ShakeLab (StreamCollection, Stream, records: any format and source), with metadata or a table of coordinates | `ArrayData(collection, select=..., coordinates=..., coordinate_covariance=...)`, or `array.add(...)` |
| Stations already built | `array.add(station)` |

```python
from shakelab.exploration.yarray import ArrayData

array = ArrayData()
array.read("array.json")
array.add(records, coordinates={"XX.A10.": (120.0, -35.0)})
array.remove("XX.A03.")
array.info()
```

The stations need not share their time axis: each keeps the interval
covered by all its channels (or `starttime`, `endtime`), and the
window common to all of them (`starttime`, `endtime`, `nsamp` of the
array) is taken when the array is processed (`array.window(...)`,
used by the beamformer). Stations with gaps are excluded with an
`ArrayDataWarning` (`gaps="drop"`) or raise (`gaps="raise"`); samples
that are not simultaneous are an error (no implicit resampling);
sampling interval and units must be the same for all the stations.

With geographic positions in the metadata the coordinates are
computed in a local ENU frame (`ArrayData.frame`). Its origin is fixed:
the `reference` given to `ArrayData` or written in the array file
(`"reference": [lon, lat, elevation]`), otherwise the position of the
first station added; adding stations does not move it.
`set_reference(point)` (or `set_reference()`, the mean of the current
stations) changes it and recomputes the coordinates.

## 10. Configuration

`BeamformingConfig` collects every parameter. Each one can be given
singly, with the name of the corresponding command-line option
(`--grid-max` is `grid_max`), or grouped in the objects
`FrequencyConfig`, `WindowConfig`, `SearchGrid` and `EllipticityGrid`
(immutable, also as dictionaries with the fields of the JSON
sections); single parameters apply on top of the objects and of the
file. The configuration changes with `update` (all the values are
validated together) or by assigning a parameter.

```python
from shakelab.exploration.yarray import BeamformingConfig

config = BeamformingConfig(fmin=2.0, fmax=12.0, nfreq=30,
                           grid_min=100.0, grid_max=2000.0, grid_size=100,
                           window_type="fixed", window_length=5.0)
config = BeamformingConfig("capon.json", method="music")  # file + changes
config.update(components="radial", maxima=3)
config.window_length = 4.0
config.write("radial.json")
```

| Parameter | Default | Meaning |
|---|---|---|
| `frequencies` | required | `FrequencyConfig`, dictionary or explicit values; or `fmin`, `fmax`, `nfreq`, `frequency_scale` |
| `grid` | required | `SearchGrid` or dictionary; or `grid_min`, `grid_max`, `grid_size`, `grid_quantity`, `grid_scale`, `azimuth_size`, `azimuths` |
| `window` | `WindowConfig()` | windows and covariance blocks; or `window_type`, `window_length`, `window_taper`, `window_overlap`, `windows_per_block`, `block_overlap` |
| `method` | `capon` | `bartlett`, `capon`, `music` |
| `components` | `vertical` | `vertical`, `radial`, `transverse`, `three_component_separate`, `rayleigh_joint` |
| `diagonal_loading` | 0.001 | relative loading of Capon |
| `music_sources` | 1 | sources excluded from the MUSIC noise subspace |
| `ellipticity_grid` | -85 to 85 degrees, 61 | `EllipticityGrid` of RTBF; or `ellipticity_angle_min`, `ellipticity_angle_max`, `ellipticity_angle_size` |
| `spectra` | False | keep also the complete spectra of every block (the maxima are always kept); they can take much memory |
| `maxima` | 1 | maxima kept per spectrum |
| `peak_neighborhood` | (3, 3) | neighbourhood of a maximum (azimuth, radial) |
| `minimum_relative_power` | 0 | smallest power of a maximum, relative to the largest |
| `statistical_uncertainty` | `none` | `empirical` (Bartlett, Capon; single components) |
| `geometric_uncertainty` | `none` | `expected_coherence` (Bartlett; needs coordinate covariances) |
| `strict` | True | more windows than the matrix dimension for Capon and MUSIC |
| `transform` | `direct` | `direct`, `rfft` (fixed windows) |
| `fft_selection` | `nearest` | `nearest`, `interp` (with `rfft`) |
| `fft_padding` | None | zero padding (1 for nearest, 8 for interp) |

| `FrequencyConfig` | Default | Meaning |
|---|---|---|
| `values` | None | explicit frequencies (Hz) |
| `minimum`, `maximum`, `size` | None | range of `size` frequencies; with the RFFT a range without `size` takes every bin |
| `scale` | `log` | `log`, `linear` |

| `WindowConfig` | Default | Meaning |
|---|---|---|
| `type` | `cycles` | `cycles` (length in cycles of each frequency) or `fixed` (seconds) |
| `length` | 20 | cycles or seconds |
| `taper` | 0.1 | Tukey fraction |
| `overlap` | 0.5 | overlap of the windows |
| `windows_per_block` | 40 | windows of a covariance block |
| `block_overlap` | 0.5 | overlap of the blocks |

| `SearchGrid` | Default | Meaning |
|---|---|---|
| `minimum`, `maximum`, `size` | required | radial range and samples |
| `quantity` | `velocity` | `velocity`, `slowness`, `wavenumber` |
| `scale` | log (linear for wavenumbers) | `log`, `linear` |
| `azimuth_size` | 180 | uniform azimuths over 360 degrees |
| `azimuths` | None | explicit azimuths |

The JSON file has the sections `frequencies`, `spectral` (`transform`,
`fft_selection`, `fft_padding`), `method` (`name`, `diagonal_loading`,
`music_sources`), `components`, `grid`, `window`, `ellipticity`
(`minimum_angle`, `maximum_angle`, `size`), `output` (`spectra`,
`maxima`, `neighborhood`, `minimum_relative_power`; the `mode` of the
August files is read as `spectra`), `uncertainty` (`statistical`,
`geometric`) and `strict`:

```json
{
  "version": 1,
  "frequencies": {"minimum": 2.0, "maximum": 12.0, "size": 100,
                  "scale": "log"},
  "spectral": {"transform": "rfft", "fft_selection": "nearest"},
  "method": {"name": "capon", "diagonal_loading": 0.001},
  "components": "vertical",
  "grid": {"quantity": "velocity", "scale": "log", "minimum": 100.0,
           "maximum": 2000.0, "size": 100, "azimuth_size": 180},
  "window": {"type": "fixed", "length": 5.0, "taper": 0.1,
             "overlap": 0.5, "windows_per_block": 40,
             "block_overlap": 0.5},
  "output": {"spectra": true, "maxima": 3},
  "uncertainty": {"statistical": "none", "geometric": "none"},
  "graphics": {"dispersion": true, "ranks": [1, 2, 3]}
}
```

Missing fields take the defaults; unknown sections are ignored by the
processing, so the section `graphics` of the command line can live in
the same file.

## 11. Processing and results

```python
from shakelab.exploration.yarray import ArrayData, Beamformer

array = ArrayData("array.json")
beamformer = Beamformer(array, "capon.json")
result = beamformer.process()                       # common window
part = beamformer.process("2026-11-01T10:00", "2026-11-01T10:10")
result.info()
result.dispersion().info()

for maximum in result.maxima_at(rank=1):
    print(maximum.time, maximum.frequency, maximum.velocity,
          maximum.azimuth, maximum.power)
```

`Beamformer(array, config)` keeps the array by reference (stations
added later are processed too) and a copy of the configuration
(`beamformer.config`, which can be changed between runs).
`cache=True` keeps the steering vectors while station coordinates and
grid do not change (repeated analyses, real time). Frequencies whose
windows do not fill one block give no result and a
`BeamformingWarning`, not an error.

`BeamformingResult` holds the maxima (`SpectralMaximum`), the spectra
if requested (`spectrum(frequency, component)` returns a
`FrequencySpectrum` with the power of every block on the grid, and
its standard error or ellipticity), the configuration, the stations
(`sids`, `coordinates`, and `stations_at(frequency, block)` for the
stations of each block), the block times and the units of the power
(`power_units`). `dispersion(rank=1)` gives, at each frequency, the
median velocity of the maxima with its 16-84 % interval, their mean
azimuth and their number (`DispersionCurve`). `append` adds the
results of later windows (blocks renumbered; results with the same
configuration, frequencies and units).

```python
result.write("daily.pkl")
beamformer.process(start, end).write("daily.pkl", append=True)
daily = BeamformingResult("daily.pkl")              # joined
```

A result file holds one result, possibly written in several pieces:
`write(append=True)` checks that the configuration is the same.
Pickle is a provisional format for trusted local files: never open one
from an untrusted source.

The density of the maxima:

```python
density = result.density(coordinate="log_velocity",
                          weighting="relative_power")
summary = density.summary()              # mode, median, quantiles
peaks = density.peaks(3, minimum_relative_density=0.1)
```

Figures (`shakelab.graphics.yarray`, matplotlib imported on use):
`plot_dispersion`, `plot_ellipticity`, `plot_spectrum` (power,
standard error or ellipticity; wavenumber plane or azimuth and grid
coordinates), `plot_density`. See the examples `yarray_01_array.py`
to `yarray_10_density.py` in
[`examples/exploration/`](../examples/README.md#exploration).

## 12. Command line

`slab-yarray` (installed with ShakeLab; `python -m
shakelab.exploration.yarray` is equivalent):

```sh
slab-yarray --input array.json --config capon.json -o result.pkl \
    --plot-dispersion --save-figures figures
slab-yarray --input array.json --config capon.json --method music \
    --fmin 2 --fmax 8 --nfreq 20 -o music.pkl     # options replace the JSON
slab-yarray --data "data/*.mseed" --metadata array.xml \
    --select "XX.*..HH?" --config rtbf.json -o rtbf.pkl --append
slab-yarray --result rtbf.pkl --plot-spectrum --frequency 4.5 --show
slab-yarray response --input array.json --fmin 1 --fmax 20 --plot-all
```

| Group | Options |
|---|---|
| Input (one) | `--input` array file; `--data` waveform files with `--metadata`, `--select`, `--start`, `--end`; `--result` saved result to plot |
| Configuration and output | `--config`, `--write-config`, `-o/--output`, `--append`, `--owrite` |
| Processing (replace the JSON) | `--frequencies` or `--fmin --fmax --nfreq`, `--frequency-scale`, `--method`, `--components`, `--diagonal-loading`, `--music-sources`, `--grid-quantity`, `--grid-scale`, `--grid-min`, `--grid-max`, `--grid-size`, `--azimuth-size`, `--window-type`, `--window-length`, `--window-taper`, `--window-overlap`, `--windows-per-block`, `--block-overlap`, `--ellipticity-angle-min/max/size`, `--spectra`, `--maxima`, `--minimum-relative-power`, `--statistical-uncertainty`, `--geometric-uncertainty`, `--strict`, `--transform`, `--fft-selection`, `--fft-padding` |
| Figures | `--plot-dispersion`, `--plot-ellipticity`, `--plot-spectrum` (`--frequency`, `--block`, `--spectrum-view`, `--spectrum-quantity`), `--plot-density` (`--density-coordinate`), `--component`, `--ranks`, `--show`, `--save-figures`, `--figure-format`, `--dpi` |

The figures can also be requested in the section `graphics` of the
configuration file (`dispersion`, `ellipticity`, `component`, `ranks`,
`spectrum` {`enabled`, `frequency`, `block`, `view`, `quantity`},
`density` {`enabled`, `coordinate`}, `show`, `save` {`directory`,
`format`, `dpi`}); the options of the command line prevail.

`slab-yarray response` analyses the geometry only (Section 7): input
`--input`, `--data` or `--coordinates` (a text file of lines `x y [z]`
or `name x y [z]`, to check a layout before deployment); levels
`--resolution-level`, `--ambiguity-level`; `--expected` (coordinate
uncertainty of the array); grid `--grid-max --grid-size`; velocity
limits at `--frequencies` or `--fmin --fmax --nfreq`; `--search-max K`
to check a search domain; figures `--plot-geometry`,
`--plot-response`, `--plot-resolution`, `--plot-velocity`,
`--plot-search-domain` or `--plot-all`. The Monte Carlo analysis
needs too many choices for a command line and is used from Python.
See [`yarray_11_command_line.sh`](../examples/exploration/yarray_11_command_line.sh).

## 13. Real time

```python
from shakelab.exploration.yarray import ArrayMonitor
from shakelab.signals.io import metadata_reader
from shakelab.signals.seedlink import Client

client = Client("server", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=900)
monitor = ArrayMonitor(config, window=600, output="monitor.pkl")
monitor.on_result(lambda monitor, result: print(result.summary()))
client.every(600, monitor)        # clock="data" to replay an archive
client.start()
...                               # from any thread, while it runs:
monitor.status(); monitor.latest(); monitor.result(start, end)
client.stop(); monitor.close()
```

The monitor is called by the client in its receiving thread: it only
copies the latest `window` seconds of the buffer (`client.get_data`)
and queues them for its worker thread, so the reception is never
delayed. The queue is short (`queue_size=1`): when the processing is
slower than the calls the oldest pending window is dropped and
reported (`on_skip`), so the monitor never accumulates delay. Windows
that cannot be processed (buffer shorter than the window, too few
stations, no complete block) are skipped and reported; an unexpected
error stops the worker and is raised by `close()`. `status()`,
`latest()`, `results(start, end)`, `result(start, end)` (joined) and
`array()` can be called at any time and return copies; each result can
also be appended
to a file (`output`). `background=False` processes in the call;
`monitor.process(collection, start, end)` processes given data. Any
object with `get_data(select)` returning a `StreamCollection` can feed
the monitor. See
[`yarray_12_realtime.py`](../examples/exploration/yarray_12_realtime.py).

## 14. Survey planning

```python
from shakelab.exploration.yarray import (
    ArrayData, ArrayResponse, PerformanceConfig,
)

response = ArrayResponse(ArrayData("array.json"))  # or coordinates
analysis = response.analyze()                      # levels 0.5
analysis.info()
limits = analysis.velocity_limits([1, 2, 5, 10, 20])
response.search_domain(0.1).info()

config = PerformanceConfig(snr_db=0.0, method="capon", grid_min=0.0,
                           grid_max=0.2, grid_size=41, trials=500,
                           realizations=40, seed=1)
wave = {"velocity": 400.0, "frequency": 8.0, "azimuth": 30.0}
response.performance(config, **wave).info()
response.performance_sweep(config, snr_db=[-10, 0, 10],
                           realizations=[10, 20, 40], **wave).info()
assessment = response.assess(config, frequencies=[2, 4, 8, 16],
                             velocities=[150, 300, 600, 1200])
assessment.info()
```

`ArrayResponse` accepts an `ArrayData` (and then its coordinate
covariances, for the expected response: `analyze(expected=True)`), an
`ArrayGeometry` or coordinates, and is computed once: the response does
not follow later changes of the array. The wavenumber grid is chosen
from the geometry (extent 1.25 times `2 pi / minimum_spacing`, 12
samples per `2 pi / maximum_aperture`); `grid_max` and `grid_size`
replace either choice.

`PerformanceConfig` works as `BeamformingConfig` (single parameters,
JSON file in the constructor, `update`, `copy(**changes)`, `read`,
`write`): `snr_db` and the search grid (`grid_min`, `grid_max`,
`grid_size`, `azimuth_size` or `azimuths`; always of wavenumbers, the
experiment does not depend on the frequency) are required; `method`
(Capon by default), `diagonal_loading`, `music_sources`, `strict`;
`trials`, `realizations`, `seed` (reproducible runs; point `i` of a
sweep or assessment uses a derived seed); `resolution_level`,
`ambiguity_level` and `boundary_margin` of the classification. It is
not derived from a beamforming configuration, whose grid is usually of
velocities and whose overlapping windows are fewer independent
realizations than `windows_per_block`.

The plane wave is given by `wavenumber`, or by `velocity` and
`frequency`, with `azimuth` (0) and `signal_power` (1); a known
frequency is kept in `result.source` and used by
`result.velocity_statistics()`. The criteria of `assess` are
parameters: `min_local` (0.9), `max_alias` (0.05), `max_other`
(0.05), `max_boundary` (0.01), `max_velocity_rmse` (0.1, relative) and
`max_grid_error` (0.1, relative to the resolution radius), with
`conservative_azimuth` (all the azimuths, or at least half). Figures:
`plot_array_geometry`, `plot_array_response`, `plot_resolution`,
`plot_search_domain`, `plot_velocity_limits`, `plot_performance`,
`plot_performance_sweep`, `plot_operational`. See
[`yarray_13_response.py`](../examples/exploration/yarray_13_response.py),
[`yarray_14_performance.py`](../examples/exploration/yarray_14_performance.py)
and [`yarray_15_operational.py`](../examples/exploration/yarray_15_operational.py).

## 15. Practical guidance and processing times

1. Check the array first (`slab-yarray response`): the velocity limits
   tell which frequencies and velocities the grid should cover.
2. Keep the configuration in a version-controlled JSON file and start
   with a compact run (`spectra=True`, few frequencies) to inspect
   spectra, number of blocks and windows per block.
3. For Capon and MUSIC use more windows per block than stations
   (`2N` for RTBF); prefer more windows to more loading.
4. Use the direct transform for cycle windows or arbitrary
   frequencies; the RFFT with nearest bins for many frequencies on
   fixed windows (the frequencies become the bins); the interpolated
   RFFT keeps the frequencies but is not faster than the direct
   transform when beamforming dominates.
5. RTBF is the most expensive analysis (a `2N` system for every
   frequency, block, azimuth and ellipticity): use moderate grids.
6. Save the results and plot or post-process them from the file,
   rather than repeating the processing.

Times on the nine SESAME stations (three components, 60 s;
[`yarray_16_profiling.py`](../examples/exploration/yarray_16_profiling.py),
one core of a laptop, indicative only): vertical Capon on
180 azimuths x 100 radial values, fixed 2 s windows, takes about
0.5 s for 50 frequencies with the direct transform and 0.13 s with the
nearest-bin RFFT (which keeps 15 distinct bins), and grows linearly
with the frequencies for the direct and the interpolated transforms;
RTBF with 61 ellipticities takes about 5.5 s for 10 frequencies; the
cache of the steering vectors halves the time of repeated analyses of
the same array. With the 38 stations of the full simulation the August
measurements were of the order of 15 s (Bartlett, 80 spectra) to
1000 s (RTBF, 60 spectra), equivalent in ShakeLab.

# Part III. Limits

## 16. Current limitations and planned extensions

- Uncertainty of the position of the maxima (wavenumber, azimuth,
  velocity, ellipticity), to be used also as kernel widths of the
  density; uncertainty of RTBF.
- Wathelet (2024) formulation of the three-component beamforming
  (ARTBF, ARDS, Newton search of the ellipticity, correction of the
  vertical-horizontal noise difference).
- Bias correction of Capon, `(L - N + 1) / L`.
- Multitaper spectral estimates with a rigorous number of degrees of
  freedom; a `CrossSpectra` class shared with frequency domain
  decomposition (FDD) and spatial autocorrelation (SPAC).
- Results in `.npz` and JSON instead of pickle.
- Mode tracking across frequencies (manual, automatic, hybrid) from
  the density peaks.
- Planning with correlated windows, coloured noise, several sources,
  sampled coordinate errors and RTBF; recommendations on aperture and
  spacing; a constrained geometry optimizer.
- Interactive graphics.

## 17. References

Capon, J. (1969). High-resolution frequency-wavenumber spectrum
analysis. Proceedings of the IEEE, 57(8), 1408-1418.

Capon, J., and Goodman, N. R. (1970). Probability distributions for
estimators of the frequency-wavenumber spectrum. Proceedings of the
IEEE, 58(10), 1785-1786.

Percival, D. B., and Walden, A. T. (1993). Spectral Analysis for
Physical Applications. Cambridge University Press.

Schmidt, R. O. (1986). Multiple emitter location and signal parameter
estimation. IEEE Transactions on Antennas and Propagation, 34(3),
276-280.

Wathelet, M. (2024). Incoherent noise-induced distortions of Rayleigh
wave ellipticity measurements obtained with three-component
beamforming. Geophysical Journal International, 236(3), 1804-1827.

Wathelet, M., Guillier, B., Roux, P., Cornou, C., and Ohrnberger, M.
(2018). Rayleigh wave three-component beamforming: signed ellipticity
assessment from high-resolution frequency-wavenumber processing of
ambient vibration arrays. Geophysical Journal International, 215(1),
507-523.

Welch, P. D. (1967). The use of fast Fourier transform for the
estimation of power spectra: a method based on time averaging over
short, modified periodograms. IEEE Transactions on Audio and
Electroacoustics, 15(2), 70-73.
