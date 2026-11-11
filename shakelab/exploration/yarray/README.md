# YArray: seismic-array processing

`shakelab.exploration.yarray` performs frequency-wavenumber (f-k)
analysis of synchronized seismic arrays: phase velocity, propagation
azimuth and, for Rayleigh waves, signed ellipticity of the waves
crossing the array, as functions of frequency and time. It is used
from Python, from the command line and within the callbacks of the
SeedLink client (real-time monitoring), with the same components.

The package derives from the YArray development code written in
August 2026 outside the repository (project MODERN) and is being
brought into ShakeLab in blocks. This file records its design, its
conventions and the state of the work.

## Status

| Block | Content | State |
|---|---|---|
| 1 | `signals.windowing`, `data` | done (commit `4574d38`) |
| 2 | `fourier.dft`, `signals.crossspectral`, `config`, `estimators`, `beamforming` | done (commit `e7a8246`) |
| 3 | `io` (array files, results), `monitor` (`ArrayMonitor`), `cli` (`slab-yarray`), `shakelab.graphics.yarray`, example | done (commit `14ecf0b`) |
| 4 | `statistics` (density of maxima), `plot_density`, `--plot-density` | done |
| 5 | `response`, `performance`, `operational` (array design and survey planning) | planned |
| 6 | Exhaustive examples of every way of use (library, command line, real time, planning), `docs/yarray.md` | planned |

## Modules

| Module | Content |
|---|---|
| `data` | `StationData`, `ArrayData`: immutable snapshots of the waveforms of the stations on a common time axis, ENU coordinates, component orientations; `ArrayData.from_collection` / `from_records` build them from ShakeLab waveforms (any source, also a SeedLink buffer) |
| `config` | `BeamformingConfig` and its sections (`FrequencyConfig`, `WindowConfig`, `SearchGrid`, `EllipticityGrid`), with `read`/`write` (JSON), `to_dict`/`from_dict`, `replace` |
| `estimators` | Bartlett, Capon, MUSIC and their statistics, on matrices and steering vectors (shared with the future performance analysis) |
| `beamforming` | `Beamformer` (`process`), `beamform`, results `BeamformingResult`, `FrequencySpectrum`, `SpectralMaximum` |
| `statistics` | `maximum_density` (density of the maxima versus frequency in wavenumber, log-slowness or log-velocity; kernels of the grid resolution, Silverman or fixed width; uniform or power weights), `MaximumDensity.summary` (mode, median, quantiles), `MaximumDensity.peaks` (candidate modes, prominences of `scipy.signal`) |
| `io` | `read_array` (array file: local coordinates, a waveform file per component, reader options), `write_result` (also `append`), `read_results`, `read_result` (pickle) |
| `monitor` | `ArrayMonitor`: periodic processing of the data of the SeedLink client (`client.every(interval, monitor)`) in a worker thread, queryable while running |
| `cli` | command line `slab-yarray` (`python -m shakelab.exploration.yarray`) |

Figures: `shakelab.graphics.yarray` (`plot_dispersion`,
`plot_ellipticity`, `plot_spectrum` of power, standard error or
ellipticity, in the wavenumber plane or in azimuth and grid
coordinates, `plot_density` with peaks or summary); matplotlib is
imported when a function is called.

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
array = ArrayData.from_collection(data, select="XX.*..HH?")

config = BeamformingConfig.read("capon_rtbf.json")
result = Beamformer(config).process(array)
result.info()

for maximum in result.maxima_at(rank=1):
    print(maximum.time, maximum.frequency, maximum.velocity,
          maximum.azimuth, maximum.ellipticity)

total = BeamformingResult.concatenate([result_1, result_2])
```

Real time:

```python
from shakelab.signals import seedlink
from shakelab.exploration.yarray import ArrayMonitor

client = seedlink.Client("server", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=900)
monitor = ArrayMonitor(Beamformer(config, cache=True), window=600,
                       output="monitor.pkl")
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

## Conventions

- Coordinates East, North, Up in metres; when computed from
  geographic positions, in the local ENU frame `ArrayData.frame`
  (default: mean position of the stations).
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
- `ArrayData` is immutable and built with new records cut on the
  common window: later changes of the input (e.g. a SeedLink buffer)
  do not affect a snapshot being processed.
- Stations not covering the window without gaps are excluded with an
  `ArrayDataWarning` (`gaps="drop"`, default) or raise an error
  (`gaps="raise"`); samples not simultaneous between stations are an
  error (no implicit resampling).
- Frequencies whose windows do not fill a covariance block (data too
  short) give no result and a `BeamformingWarning` instead of an
  error, as needed by real-time processing of short windows.
- `Beamformer(cache=True)` keeps the steering vectors while the
  coordinates do not change (about 15% faster for repeated calls).
- Options with fixed values are enumerations; their string values are
  accepted everywhere.

## Changes from the August 2026 version

- Same algorithms: on the SESAME simulation (3 and 38 stations) the
  maxima (wavenumber, azimuth, rank, block, ellipticity) are identical
  for every method and component mode (test against
  `tests/exploration/yarray/data/reference_20260812.json`); processing
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
  by `BeamformingConfig.read`; `npts` by `nsamp`; result times are
  `Date`; `BeamformingResult` gains `block_times`, `sids`,
  `coordinates`, `units`, `maxima_at`, `spectrum`, `concatenate`,
  `info`.
- The SAC reader accepts the SESAME files (local coordinates in
  STLA/STLO).

## Notes and open points

- Capon estimated from L windows underestimates the power by a factor
  of about `(L - n + 1) / L` for n stations (Capon and Goodman, 1970):
  -8% with 40 windows and 10 stations in the tests. A bias correction
  could be offered as an option.
- RTBF is unstable for ellipticities near zero or very large with
  incoherent noise; Wathelet (2024) solves it with the angular
  ellipticity and an all-component steering matrix: candidate
  extension.
- Empirical uncertainty: Bartlett and Capon on single components (not
  MUSIC, not RTBF); coordinate uncertainty: Bartlett only.
- Recovered in block 3: processing options of the command line as
  overrides of the JSON (runner of 10 August), reader options per
  component in the array file, appended results and `--result`, maps
  of the power uncertainty (`plot_spectrum(quantity=
  "standard_error")`). Still to do: profiling suite as an example
  (block 6).
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
- Documentation: the four documents of August (theory, architecture,
  user manual, array response) will become `docs/yarray.md`, updated.

## Test data

`tests/exploration/yarray/data/`: excerpt of the SESAME array
simulation (3 stations, 60 s; see its README) and the maxima of the
August implementation on it.

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
