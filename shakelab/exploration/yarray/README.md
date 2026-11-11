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
| 2 | `fourier.dft`, `signals.crossspectral`, `config`, `estimators`, `beamforming` | done |
| 3 | JSON array input, results I/O, command line `slab-yarray`, `ArrayMonitor` (SeedLink), graphics in `shakelab/graphics` | planned |
| 4 | `statistics` (density of maxima) | planned |
| 5 | `response`, `performance`, `operational` (array design and survey planning) | planned |

## Modules

| Module | Content |
|---|---|
| `data` | `StationData`, `ArrayData`: immutable snapshots of the waveforms of the stations on a common time axis, ENU coordinates, component orientations; `ArrayData.from_collection` / `from_records` build them from ShakeLab waveforms (any source, also a SeedLink buffer) |
| `config` | `BeamformingConfig` and its sections (`FrequencyConfig`, `WindowConfig`, `SearchGrid`, `EllipticityGrid`), with `read`/`write` (JSON), `to_dict`/`from_dict`, `replace` |
| `estimators` | Bartlett, Capon, MUSIC and their statistics, on matrices and steering vectors (shared with the future performance analysis) |
| `beamforming` | `Beamformer` (`process`), `beamform`, results `BeamformingResult`, `FrequencySpectrum`, `SpectralMaximum` |

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

Real time (block 3): `client.every(600, ArrayMonitor(Beamformer(config,
cache=True), window=600, ...))`; the callback only copies the data from
the buffer (`client.get_data`), the processing runs in a worker thread
of the monitor, which can be queried from any thread (`status()`,
`latest()`, `result(start, end)`, `array()`).

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
- Recovered functionality to implement in block 3: processing options
  of the command line as overrides of the JSON (runner of 10 August),
  reader options per component in the array JSON, appended results
  and `--result` (plot without processing), `plot_density` and maps of
  the power uncertainty (only in scripts in August), profiling suite
  as an example.
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
