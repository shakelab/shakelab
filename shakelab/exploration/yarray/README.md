# YArray: seismic-array processing

`shakelab.exploration.yarray` performs frequency-wavenumber (f-k)
analysis of synchronized seismic arrays: phase velocity, propagation
azimuth and, for Rayleigh waves, signed ellipticity of the waves
crossing the array, as functions of frequency and time. It is used
from Python, from the command line and within the callbacks of the
SeedLink client (real-time monitoring), with the same components.

The package derives from the YArray development code written in
August 2026 outside the repository (project MODERN), brought into
ShakeLab in blocks. This file records its design, its conventions and
the state of the work; the theory and the use are in the guide
[`docs/yarray.md`](../../../docs/yarray.md), the runnable examples in
[`examples/exploration/`](../../../examples/README.md#exploration).

## Status

| Block | Content | State |
|---|---|---|
| 1 | `signals.windowing`, `data` | done (commit `4574d38`) |
| 2 | `fourier.dft`, `signals.crossspectral`, `config`, `estimators`, `beamforming` | done (commit `e7a8246`) |
| 3 | `io` (array files, results), `monitor` (`ArrayMonitor`), `cli` (`slab-yarray`), `shakelab.graphics.yarray`, example | done (commit `14ecf0b`) |
| 4 | `statistics` (density of maxima), `plot_density`, `--plot-density` | done (commit `d1c1f9c`) |
| 5 | `response`, `performance`, `operational` (array design and survey planning), their figures, `slab-yarray response` | done (commit `21b3c8e`) |
| 6 | Nine-station three-component test data and references; examples of every way of use (library, command line, real time, planning, profiling); `docs/yarray.md` | done (commit `f711fb0`) |
| A | Interface: `ArrayData` filled progressively (time window taken at processing, fixed ENU origin), `BeamformingConfig` with single parameters, `Beamformer(array, config).process(start, end)`, results with `read`/`write`/`append`, stations of each block, `dispersion`, `density`; shorter examples | done (commit `02dee03`) |
| B | Interface of the survey planning: `ArrayResponse` with methods (`analyze`, `search_domain`, `performance`, `performance_sweep`, `assess`), flat `PerformanceConfig`, criteria as parameters; examples 13-15 | in progress |

## Modules

| Module | Content |
|---|---|
| `data` | `StationData` (immutable: records, ENU coordinates, component orientations, geographic position), `ArrayData` (filled from array files, ShakeLab waveforms of any source, also a SeedLink buffer, or stations: `read`, `add`, `remove`; common window with `window`; fixed ENU origin, `set_reference`) |
| `config` | `BeamformingConfig` (single parameters named as the command-line options, or the objects `FrequencyConfig`, `WindowConfig`, `SearchGrid`, `EllipticityGrid`; JSON file in the constructor, `read`, `write`, `update`, `copy`), `PerformanceConfig` (same mechanics: SNR, estimator, wavenumber grid, trials, realizations, seed, classification levels) |
| `estimators` | Bartlett, Capon, MUSIC and their statistics, on matrices and steering vectors (shared with the future performance analysis) |
| `beamforming` | `Beamformer(array, config)` (`process(starttime, endtime)`), results `BeamformingResult` (`append`, `read`, `write`, `dispersion`, `density`, `stations_at`), `FrequencySpectrum`, `SpectralMaximum`, `DispersionCurve` |
| `statistics` | `maximum_density` (density of the maxima versus frequency in wavenumber, log-slowness or log-velocity; kernels of the grid resolution, Silverman or fixed width; uniform or power weights), `MaximumDensity.summary` (mode, median, quantiles), `MaximumDensity.peaks` (candidate modes, prominences of `scipy.signal`) |
| `response` | `ArrayResponse(array or coordinates)`: geometric response (also expected under coordinate uncertainty) and entry point of survey planning: `analyze` (`ResponseAnalysis`: main lobe, ambiguities, limits, `velocity_limits`), `search_domain`, `performance`, `performance_sweep`, `assess`; `ArrayGeometry`, `ResponseGrid` |
| `performance` | Monte Carlo performance of Bartlett, Capon and MUSIC for one plane wave in white noise: LOCAL/ALIAS/OTHER/BOUNDARY outcomes, Wilson intervals, local errors, exact velocity statistics (`PerformanceResult`, `PerformanceSweepResult`, `PlaneWaveSource`) |
| `operational` | frequency-velocity-azimuth envelope with explicit usability criteria, limiting factors and recommendations (`OperationalAssessment`) |
| `io` | file formats: array files (local coordinates, a waveform file per component, reader options) and result files (pickle), used by `ArrayData.read` and `BeamformingResult.read`/`write` |
| `monitor` | `ArrayMonitor`: periodic processing of the data of the SeedLink client (`client.every(interval, monitor)`) in a worker thread, queryable while running |
| `cli` | command line `slab-yarray` (`python -m shakelab.exploration.yarray`) |

Figures: `shakelab.graphics.yarray` (`plot_dispersion`,
`plot_ellipticity`, `plot_spectrum` of power, standard error or
ellipticity, in the wavenumber plane or in azimuth and grid
coordinates, `plot_density` with peaks or summary;
`plot_array_geometry`, `plot_array_response`, `plot_resolution`,
`plot_search_domain`, `plot_velocity_limits`, `plot_performance`,
`plot_performance_sweep`, `plot_operational`); matplotlib is imported
when a function is called.

Generic tools live in `shakelab.signals`, for reuse by FDD, SPAC and
other methods: `windowing` (windows, blocks, tapers, effective number
of windows), `crossspectral` (Fourier coefficients of windows, direct
or RFFT; cross-spectral matrices with variance), `fourier.dft`
(transform at arbitrary frequencies).

## Usage

```python
from shakelab.signals.io import reader, metadata_reader
from shakelab.exploration.yarray import (
    ArrayData, Beamformer, BeamformingConfig, BeamformingResult,
)

data = reader("data/*.mseed", metadata=metadata_reader("array.xml"))
array = ArrayData(data, select="XX.*..HH?")

config = BeamformingConfig("capon_rtbf.json", window_length=4.0)
beamformer = Beamformer(array, config)
result = beamformer.process()                       # common window
result.info()
result.dispersion().info()

result.append(beamformer.process("2026-11-01T10:00", "2026-11-01T10:10"))
result.write("result.pkl")
result = BeamformingResult("result.pkl")
```

Real time:

```python
from shakelab.signals.seedlink import Client
from shakelab.exploration.yarray import ArrayMonitor

client = Client("server", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=900)
monitor = ArrayMonitor(config, window=600, output="monitor.pkl")
monitor.on_result(lambda monitor, result: print(result.summary()))
client.every(600, monitor)          # clock="data" to replay an archive
client.start()
...                                 # from any thread, while it runs:
monitor.status(); monitor.latest(); monitor.result(start, end)
monitor.array(); client.get_data(...)
client.stop(); monitor.close()
```

The callback, run in the receiving thread of the client, only copies
the data of the buffer (`client.get_data`); the processing runs in the
worker thread of the monitor. The queue is short: when the processing
is slower than the calls the oldest pending window is dropped and
reported (`on_skip`), so the monitor never accumulates delay. Windows
that cannot be processed (no data, too few stations, too short) are
skipped and reported; unexpected errors stop the worker and are raised
by `close()`. `background=False` processes in the call;
`monitor.process(collection, start, end)` processes given data.

Command line:

```sh
slab-yarray --input array.json --config capon.json -o result.pkl \
    --plot-dispersion --save-figures figures
slab-yarray --input array.json --config capon.json --fmin 2 --fmax 8 \
    --nfreq 20 --method music -o music.pkl      # options replace the JSON
slab-yarray --data "data/*.mseed" --metadata array.xml \
    --select "XX.*..HH?" --config rtbf.json -o rtbf.pkl --append
slab-yarray --result rtbf.pkl --plot-spectrum --frequency 4.5 --show
```

The figures can also be requested in the section `graphics` of the
configuration file (`dispersion`, `ellipticity`, `component`, `ranks`,
`spectrum` {`enabled`, `frequency`, `block`, `view`, `quantity`},
`show`, `save` {`directory`, `format`, `dpi`}), as in August.

Array design and survey planning (geometry only, then estimator
performance):

```python
from shakelab.exploration.yarray import (
    ArrayData, ArrayResponse, PerformanceConfig,
)

response = ArrayResponse(ArrayData("array.json"))   # or coordinates
analysis = response.analyze()                # levels 0.5 by default
analysis.info()                              # resolution, aliasing
limits = analysis.velocity_limits([1, 2, 5, 10, 20])
response.search_domain(0.05).info()

config = PerformanceConfig(snr_db=0.0, method="capon", grid_min=0.0,
                           grid_max=0.2, grid_size=41, trials=500,
                           realizations=40, seed=1)  # or a JSON file
wave = {"velocity": 400.0, "frequency": 8.0, "azimuth": 30.0}
response.performance(config, **wave).info()  # or wavenumber=...
response.performance_sweep(config, snr_db=[-10, 0, 10],
                           realizations=[10, 20, 40], **wave).info()
assessment = response.assess(config, frequencies=[2, 4, 8, 16],
                             velocities=[150, 300, 600, 1200],
                             min_local=0.8)  # criteria: min_*, max_*
assessment.info()                            # usable ranges, advice
```

From the command line, the geometric part only (the Monte Carlo
analysis needs too many choices for a command line):

```sh
slab-yarray response --input array.json --plot-all --show
slab-yarray response --coordinates stations.txt --fmin 1 --fmax 20 \
    --search-max 0.05 --plot-all --save-figures figures
```

`--coordinates` reads a text file of lines `x y [z]` or
`name x y [z]`, so that a layout can be checked before deployment.

## Conventions

- Coordinates East, North, Up in metres; when computed from
  geographic positions, in the local ENU frame `ArrayData.frame`
  (default: the first station added; see `ArrayData`).
- Azimuths in degrees clockwise from north; propagation azimuth (the
  direction towards which the wave travels).
- Wavevector `k = 2 pi f / c`, `kx = k sin(az)`, `ky = k cos(az)`;
  steering vectors `exp(-1j (kx x + ky y))`; Fourier kernel
  `exp(-2j pi f t)`.
- Radial = horizontal motion along the propagation azimuth;
  transverse = radial rotated by +90 degrees (clockwise).
- Signed ellipticity (RTBF, Wathelet et al., 2018): radial coefficient
  `-1j e` times the vertical one; `e > 0` for radial motion
  `e sin(wt)` with vertical `cos(wt)`.
- Units: the cross-spectral matrices are one-sided cross power
  spectral densities, so that Bartlett and Capon give the density of a
  plane wave in the units of the data squared per Hz. Data not
  corrected for the instrument work the same way (counts squared per
  Hz); `BeamformingResult.power_units` reports the units (`unit` when
  unknown). MUSIC is a pseudo-power (no units); the RTBF objective
  `e**2 P**2` has the square of the units of the power.
- Times: `Date`; block times are the centres of the covariance blocks.

## Design choices

- Layers: `signals` primitives (no array concepts) < `data` <
  `estimators` < `config` < `beamforming` < applications (command
  line, monitor, graphics). Numerical modules have no file I/O, no
  plotting and no global state.
- `ArrayData` is a container filled progressively (decision of
  2026-11-11): the stations keep their own time axis and the common
  window is taken when processing (`window`, used by `Beamformer`), so
  that stations from different sources can be added one at a time.
  `StationData` stays immutable and the builders create new records:
  later changes of the input (e.g. a SeedLink buffer) do not affect
  the array.
- The origin of the ENU frame is fixed at the first station added
  with a geographic position, or given (`reference`, or the
  `reference` of the array file); adding
  stations does not move it; `set_reference` recomputes the
  coordinates.
- Stations not covering the window without gaps are excluded with an
  `ArrayDataWarning` (`gaps="drop"`, default) or raise an error
  (`gaps="raise"`); samples not simultaneous between stations are an
  error (no implicit resampling).
- `BeamformingConfig`: one name per parameter, the same as the
  command-line option; the objects (`WindowConfig`...) stay immutable,
  the configuration changes with `update`, validated as a whole.
  `Beamformer` and `ArrayMonitor` keep a copy, the results their own.
- A result file holds one result, possibly written in pieces
  (`append=True`, same configuration checked); the stations and
  coordinates are recorded for each block, so that results of windows
  with different stations (real time with `gaps="drop"`) join
  correctly.
- Frequencies whose windows do not fill a covariance block (data too
  short) give no result and a `BeamformingWarning` instead of an
  error, as needed by real-time processing of short windows.
- `Beamformer(cache=True)` keeps the steering vectors while the
  coordinates and the grid do not change (about 15-30% faster for
  repeated calls).
- Options with fixed values are enumerations; their string values are
  accepted everywhere.
- Survey planning starts from `ArrayResponse`, computed once from the
  coordinates of the moment (an `ArrayData` changed later does not
  change it); its methods replace free functions. `PerformanceConfig`
  shares the mechanics of `BeamformingConfig` (a private base in
  `config`), but not its search grid: the Monte Carlo grid is of
  wavenumbers, independent of the frequency, while beamforming grids
  are usually of velocities; and overlapping windows are fewer
  independent realizations than `windows_per_block`. So a performance
  configuration is not derived from a beamforming one. The plane wave
  is given by parameters (`wavenumber`, or `velocity` and `frequency`;
  `azimuth`); the frequency, if known, is kept in `result.source` and
  used by `velocity_statistics()`. The criteria of `assess` are single
  parameters (`min_local`, `max_alias`, `max_other`, `max_boundary`,
  `max_velocity_rmse`, `max_grid_error`) recorded in
  `assessment.criteria`. Results keep a copy of the configuration.

## Changes from the August 2026 version

- Same algorithms: on the SESAME simulation (3 and 38 stations) the
  maxima (wavenumber, azimuth, rank, block, ellipticity) are identical
  for every method and component mode (test against
  `tests/exploration/yarray/data/reference.json`); processing
  times are equivalent (38 stations, Capon, 50 frequencies: 10 s
  direct, 7 s RFFT).
- Power in physical units (see Conventions): the old coefficients were
  `2 / N sum w x e` with the taper normalized to unit mean square, so
  that the new power is `N delta / 2` times the old one (squared for
  RTBF; maxima and relative thresholds unchanged).
- Defaults of `WindowConfig` unified with those of the JSON files
  (overlap 0.5, 40 windows per block, block overlap 0.5; the Python
  classes had 0, 20, 0).
- Spectral code moved to `signals` (`windowing`, `crossspectral`,
  `fourier.dft`), with the conventions of `fourier`.
- `Yarray` (runner) replaced by `Beamformer`; `read_processing_config`
  by `BeamformingConfig(file)`; `npts` by `nsamp`; result times are
  `Date`; `BeamformingResult` gains `block_times`, `sids`,
  `coordinates`, `units`, `maxima_at`, `spectrum`, `append`,
  `stations_at`, `dispersion`, `density`, `info`.
- The SAC reader accepts the SESAME files (local coordinates in
  STLA/STLO).
- Block 5 (`response`, `performance`, `operational`): same results as
  August, identical to the last bit with the same seeds (responses,
  limits, Monte Carlo draws and outcomes, envelopes; regression tests
  in `test_yarray_performance.py` and `test_yarray_operational.py`).
  `print_*_summary` and `*_summary` became `info()` methods;
  `analyze_response` returns a `ResponseAnalysis` instead of a tuple
  `(main_lobe, ambiguity, limits)`; `array_response` and
  `ResponseGrid.from_geometry` also accept an `ArrayData` (with its
  coordinate covariances) or an `ArrayGeometry`; `use_expected_response`
  is `expected`; methods and outcomes accept strings;
  `PlaneWaveSource.from_velocity`; the sample matrices use
  `signals.crossspectral.cross_spectral_matrix`. The classification
  maps and the steering vectors of the search grid are computed once
  per sweep or assessment (about 30% faster). Fixes: the limiting
  factor converted the velocity at 1 Hz instead of the frequency of the
  point (no effect on the result, the relative RMSE does not depend on
  it); the recommendation on grid quantization compared the error with
  the Rayleigh proxy, the point criterion with the resolution radius
  of the response: both use the resolution radius now.
- `SearchGrid` scale defaults to linear for wavenumbers (the only
  spacing allowed) and log otherwise; before, a wavenumber grid needed
  `scale="linear"` explicitly.
- Interface of the planning (November 2026, block B): `ArrayResponse`
  with methods instead of `array_response`, `analyze_response`,
  `analyze_search_domain`, `velocity_limits`, `evaluate_performance`,
  `performance_sweep` and `assess_array`; flat `PerformanceConfig`
  instead of `EstimatorConfig`, `MonteCarloConfig` and
  `ClassificationConfig`; criteria as parameters instead of
  `OperationalConfig` and `OperationalCriteria`; default method Capon
  (Bartlett before), as in `BeamformingConfig`. Same results, bit for
  bit, with the same seeds.

## Notes and open points

- Capon estimated from L windows underestimates the power by a factor
  of about `(L - n + 1) / L` for n stations (Capon and Goodman, 1970):
  -8% with 40 windows and 10 stations in the tests. A bias correction
  could be offered as an option.
- RTBF follows Wathelet et al. (2018): steering vectors
  `[-1j e a, a]` and objective `e**2 P**2`. Only the search grid is
  angular (`e = tan(chi)`, `chi` within (-90, 90) degrees, zero
  excluded), so vanishing and very large ellipticities cannot be
  reached and are distorted by incoherent noise. Wathelet (2024)
  reformulates the problem in the angular ellipticity, with a
  different noise model and a new steering matrix: candidate
  extension, not implemented (nor in the August code).
- Empirical uncertainty: Bartlett and Capon on single components (not
  MUSIC, not RTBF); coordinate uncertainty: Bartlett only.
- Recovered in block 3: processing options of the command line as
  overrides of the JSON (runner of 10 August), reader options per
  component in the array file, appended results and `--result`, maps
  of the power uncertainty (`plot_spectrum(quantity=
  "standard_error")`). The profiling suite of August (shell script and
  configurations) became `examples/exploration/yarray_16_profiling.py`
  (block 6).
- Block 6: the test data hold nine three-component SESAME stations
  (`sesame_array.json`); the maxima of the August code on them
  (vertical, radial, transverse, three components, Rayleigh joint) are
  identical (`reference.json`, which replaces the three-station
  reference). `ArrayMonitor` skips, without warnings about stations, a
  buffer shorter than the window (the start of an acquisition).
- `statistics` (block 4) gives the same densities, summaries and peak
  positions as August (1e-14); two fixes: `numpy.trapz` (removed in
  NumPy 2.4) replaced by `trapezoid`, and the peak prominence, which
  was measured from the lowest density on each side instead of the
  saddle towards a higher peak (secondary modes on the flank of a
  stronger one were overrated: 0.251 instead of 0.153 in the test
  case), now computed by `scipy.signal.find_peaks`. The option
  `use_uncertainty` was removed: maxima have no positional errors yet
  (candidate extension).
- The local SeedLink test server sends the stations one after the
  other: the integration test checks the call in the receiving thread
  and then the analysis of the complete buffer.
- Planning (block 5): the experiments model one scalar plane wave in
  white noise with independent realizations; correlated windows,
  coloured noise, several sources, sampled coordinate errors and RTBF
  are not modelled. The recommendation categories `ARRAY_APERTURE`,
  `ARRAY_SPACING`, `AZIMUTH_SAMPLING` and the limiting factor
  `INSUFFICIENT_REALIZATIONS` are defined but never assigned (as in
  August): a low LOCAL probability is reported as `LOW_SNR`. A
  geometry optimizer (constrained station placement) is a candidate
  extension.
- On the 3-station SESAME excerpt the main lobe (radius 0.060 rad/m)
  is wider than the alias-free domain (0.057 rad/m): no wavenumber is
  both resolved and unambiguous at level 0.5, as `info()` reports.
- Documentation: the four documents of August (theory, architecture,
  user manual, array response) are merged and updated in
  `docs/yarray.md`.

## Test data

`tests/exploration/yarray/data/`: excerpt of the SESAME array
simulation (9 three-component stations, 60 s, with their array file;
see its README) and the maxima of the August implementation on it.

## References

Capon, J. (1969). High-resolution frequency-wavenumber spectrum
analysis. Proceedings of the IEEE, 57(8), 1408-1418.

Capon, J., and Goodman, N. R. (1970). Probability distributions for
estimators of the frequency-wavenumber spectrum. Proceedings of the
IEEE, 58(10), 1785-1786.

Wathelet, M. (2024). Incoherent noise-induced distortions of Rayleigh
wave ellipticity measurements obtained with three-component
beamforming. Geophysical Journal International, 236(3), 1804-1827.

Wathelet, M., Guillier, B., Roux, P., Cornou, C., and Ohrnberger, M.
(2018). Rayleigh wave three-component beamforming: signed ellipticity
assessment from high-resolution frequency-wavenumber processing of
ambient vibration arrays. Geophysical Journal International, 215(1),
507-523.
