# YARRAY: seismic array processing

`shakelab.exploration.yarray` analyses the wavefield recorded by a
seismic array, that is a group of stations recording at the same time
over an area. The waves that cross the array reach the stations with
small delays, which depend on their speed and direction: by comparing
the recordings in the frequency domain, YARRAY estimates the phase
velocity and the propagation azimuth of these waves as functions of
frequency and time. For Rayleigh waves recorded on three components it
also measures the signed ellipticity of the particle motion.

These measurements are the basis of passive surface-wave methods
(dispersion curves for site characterization, ambient-vibration
arrays), of the monitoring of the wavefield in time, and of the study
of the sources of seismic noise. Since the quality of the results
depends first of all on the geometry of the array, YARRAY also helps
to design an array before a survey, from the station layout alone.

The same objects are used in three ways: as a Python library, from the
command line (`slab-yarray`), and in real time within the callbacks of
the ShakeLab SeedLink client.

## 1. Features

- Array data from any source read by ShakeLab (MiniSEED, SAC, SDS
  archives, FDSN services, a SeedLink buffer...), with coordinates from
  the station metadata or from a simple array file.
- Frequency-wavenumber (f-k) beamforming with the Bartlett, Capon and
  MUSIC estimators.
- Vertical, radial, transverse and three-component processing;
  Rayleigh three-component beamforming with signed ellipticity
  (RTBF, ARTBF, ARDS).
- Spectral maxima, dispersion curves and the density of the maxima,
  to identify the modes of propagation.
- Uncertainty of the beam power and of the station coordinates.
- Real-time monitoring of an array with the SeedLink client.
- Survey planning: array response, resolution and aliasing limits,
  Monte Carlo performance of the estimators and operational envelopes.
- Figures for every result (`shakelab.graphics.yarray`).

## 2. Quick start

A typical analysis needs three objects: the array data (`ArrayData`,
the stations with their recordings and positions), a configuration
(`BeamformingConfig`, which says what to compute and how) and a
`Beamformer`, which applies the configuration to the data and returns
a `BeamformingResult`.

In this example the array is described by an array file (Section 3.1)
and the vertical component is processed with the Capon estimator,
searching phase velocities between 100 and 1500 m/s at 20 frequencies
between 2 and 10 Hz:

```python
import matplotlib.pyplot as plt

from shakelab.exploration.yarray import (
    ArrayData,
    Beamformer,
    BeamformingConfig,
)
from shakelab.graphics.yarray import plot_dispersion

array = ArrayData("array.json")
array.info()

config = BeamformingConfig(
    fmin=2.0, fmax=10.0, nfreq=20,                 # frequencies (Hz)
    grid_min=100.0, grid_max=1500.0, grid_size=60,  # velocities (m/s)
)
result = Beamformer(array, config).process()       # common time window
result.dispersion().info()

plot_dispersion(result)
plt.show()
```

The recordings are divided into short windows, grouped into covariance
blocks. For each block and frequency the beamformer computes the beam
power over a grid of velocities and azimuths and keeps its maxima:
velocity, azimuth, power and, for Rayleigh processing, ellipticity.
The dispersion curve summarizes the maxima of all the blocks at each
frequency.

## 3. Common tasks

### 3.1 Build the array data

An `ArrayData` holds the stations of the array, each with its
recordings (one to three components), its position and the
orientation of its components. It can be built in two ways, depending
on how the data are available.

**From waveforms read by ShakeLab.** When the recordings come with
station metadata (StationXML, an FDSN service, a SeedLink server with
metadata), the waveforms are read with `shakelab.signals.io` or any
other ShakeLab tool and passed to `ArrayData`. Positions and
orientations are taken from the metadata, and the geographic positions
are converted to local East-North-Up coordinates in metres. A channel
pattern selects the stations and components to use:

```python
from shakelab.signals.io import reader, metadata_reader
from shakelab.exploration.yarray import ArrayData

data = reader("data/*.mseed", metadata=metadata_reader("array.xml"))
array = ArrayData(data, select="XX.*..HH?")
```

Without metadata, the local coordinates can be given directly, in
metres, for each station (`NET.STA.LOC`):

```python
array = ArrayData(data, coordinates={"XX.A01.": (0.0, 0.0),
                                     "XX.A02.": (25.0, -10.0)})
```

**From an array file.** For data without usable metadata, such as SAC
files of a temporary survey, a JSON array file lists the stations,
their local coordinates and the waveform file of each component:

```json
{
  "version": 1,
  "coordinate_system": {"type": "local_cartesian", "unit": "m"},
  "data_root": "sac",
  "stations": [
    {"id": "A01",
     "coordinates": {"x": 0.0, "y": 0.0, "z": 0.0},
     "components": {"Z": "A01.Z.sac", "N": "A01.N.sac",
                    "E": "A01.E.sac"}},
    {"id": "A02",
     "coordinates": {"x": 25.0, "y": -10.0, "z": 0.0},
     "components": {"Z": "A02.Z.sac", "N": "A02.N.sac",
                    "E": "A02.E.sac"}}
  ]
}
```

```python
array = ArrayData("array.json")
```

The file paths are relative to `data_root` (by default, the folder of
the array file). Components named Z, N and E need no orientation;
other names take an explicit one. Optional fields give the geographic
position of the origin (`reference`), the reader options of each file
and the uncertainty of the coordinates; the complete format is
described in the [guide](../../../docs/exploration/yarray.md#array-file).

An `ArrayData` can also be filled one station at a time (`add`,
`remove`), mixing different sources. The stations keep their own time
axis: the window common to all of them is taken only when the array is
processed, and stations with gaps in that window are left out with a
warning.

### 3.2 Configure the processing

A `BeamformingConfig` collects all the choices of the processing.
Each parameter has the same name as the corresponding command-line
option (`--grid-max` is `grid_max`), so that a configuration can be
written in Python, in a JSON file or on the command line in the same
way. Parameters not given take their default values; frequencies and
search grid are required.

```python
from shakelab.exploration.yarray import BeamformingConfig

config = BeamformingConfig("capon.json", window_length=4.0)
config.update(method="music", maxima=2)
config.write("music.json")
```

Here a configuration is read from a file with one value replaced,
then changed and saved. The changes are validated as a whole, so that
a configuration is never left in an inconsistent state.

The main choices are the estimator (`method`: `"capon"`, the default,
`"bartlett"` or `"music"`), the components to process (`components`),
the frequencies (`fmin`, `fmax`, `nfreq`), the search grid of
velocities, slownesses or wavenumbers (`grid_min`, `grid_max`,
`grid_size`, `grid_quantity`) and the time windows (`window_type`,
`window_length`, `windows_per_block`). Windows can have a fixed
length in seconds or a fixed number of cycles at each frequency.

### 3.3 Rayleigh waves and ellipticity

With three-component stations, the radial and vertical motions of
Rayleigh waves can be processed together. Besides velocity and
azimuth, this gives the ellipticity of the particle motion, a property
of the subsoil that complements the dispersion curve:

```python
config = BeamformingConfig(
    "capon.json",
    components="rayleigh_angular",   # or "rayleigh_joint", "rayleigh_direct"
)
result = Beamformer(array, config).process()
```

Three formulations are available: `"rayleigh_joint"` (RTBF, Wathelet
et al., 2018), which searches the ellipticity on a grid, and
`"rayleigh_angular"` (ARTBF) and `"rayleigh_direct"` (ARDS) of
Wathelet (2024), which find the ellipticity angle exactly, also for
purely vertical or purely horizontal motion. The ellipticity is the
ratio of radial to vertical amplitude, positive for prograde and
negative for retrograde particle motion (Section 6).

Single components can also be processed separately (`"vertical"`,
`"radial"`, `"transverse"`, `"three_component_separate"`), for
example to compare Rayleigh and Love waves.

### 3.4 Results: maxima, dispersion, density, files

A `BeamformingResult` keeps the maxima of every block and frequency,
the stations used in each block and the configuration that produced
it. The maxima can be selected by rank, frequency and component, or
summarized in a dispersion curve:

```python
from shakelab.exploration.yarray import BeamformingResult

best = result.maxima_at(rank=1, frequency=5.0)
curve = result.dispersion()

density = result.density(coordinate="log_velocity")
density.summary()                       # mode, median, quantiles
density.peaks(2)                        # candidate modes

result.write("result.pkl")
result = BeamformingResult("result.pkl")
```

With many blocks the maxima scatter around the dispersion curves of
the different modes. Their density at each frequency, in velocity,
slowness, wavenumber or ellipticity angle, shows these curves more
clearly than single maxima; its summary and its peaks give the
candidate modes.

Results of successive time windows, as in long recordings or in real
time, can be joined with `append`, also in the same file
(`write(..., append=True)`). The complete beam-power spectra are kept
only on request (`spectra=True` in the configuration), since they can
take much memory.

### 3.5 Real-time monitoring

During an acquisition, the SeedLink client receives the data in a
background thread and calls an `ArrayMonitor` at regular intervals.
At each call the monitor copies the latest window of data from the
client buffer and processes it in its own thread, so that the
reception is never delayed; the results can be queried at any time
while the acquisition runs, or written to a file:

```python
from shakelab.signals.io import metadata_reader
from shakelab.signals.seedlink import Client
from shakelab.exploration.yarray import ArrayMonitor

client = Client("rtserve.example.org", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=600)

monitor = ArrayMonitor(config, window=300, output="monitor.pkl")
client.every(60, monitor)
client.start()

monitor.latest()                        # from any thread
```

Here the last 300 s of data are processed every 60 s. Functions can be
called on each result (`on_result`) or when a window cannot be
processed (`on_skip`), for example because too few stations have data.
With `clock="data"` the calls follow the time of the data instead of
the wall clock, which allows to replay an archive through a SeedLink
server.

### 3.6 Plan a survey

The ability of an array to resolve waves of a given wavelength, and
to tell them apart from aliases, depends only on the positions of the
stations. The array response can therefore be computed from a planned
layout, before the survey, to check which velocities can be measured
at each frequency:

```python
from shakelab.exploration.yarray import ArrayResponse, PerformanceConfig

response = ArrayResponse(array)         # or coordinates in m
analysis = response.analyze()           # resolution and aliasing
analysis.velocity_limits([2.0, 5.0, 10.0])

config = PerformanceConfig(snr_db=0.0, grid_min=0.0, grid_max=0.2,
                           grid_size=41, trials=200, seed=1)
response.performance(config, velocity=400.0, frequency=8.0).info()
```

The geometric analysis gives the resolution and the aliasing limits.
The performance experiments go further: they simulate a plane wave in
noise many times (Monte Carlo) and count how often the estimator finds
it, for a given signal-to-noise ratio and number of windows.
`response.assess(...)` combines these results into the
frequency-velocity range in which the array is usable, with the
factors that limit it.

## 4. Command line

The `slab-yarray` command runs the same processing without writing
Python. The data are given as an array file (`--input`) or as
waveforms with metadata (`--data`, `--metadata`, `--select`); the
configuration comes from a JSON file, and any option on the command
line replaces the corresponding value of the file. Figures can be
shown or saved, also from a result written before (`--result`):

```sh
slab-yarray --input array.json --config capon.json -o result.pkl \
    --plot-dispersion --show
slab-yarray --data "data/*.mseed" --metadata array.xml \
    --select "XX.*..HH?" --config capon.json --components rayleigh_angular
slab-yarray --result result.pkl --plot-density
slab-yarray response --input array.json --plot-all --show
```

The subcommand `slab-yarray response` computes the array response
and its limits; the Monte Carlo experiments, which need many more
choices, are available from Python only. `slab-yarray --help` lists
all the options.

## 5. Structure

The package is organized in layers: the data and the configuration
at the bottom, then the estimators and the beamforming, and on top the
applications (command line, real-time monitor, figures). Each module
has a single task:

| Module | Main objects | Role |
|---|---|---|
| `data` | `ArrayData`, `StationData` | stations, coordinates, common time window |
| `config` | `BeamformingConfig`, `PerformanceConfig` | processing and planning parameters |
| `estimators` | | Bartlett, Capon and MUSIC on cross-spectral matrices |
| `beamforming` | `Beamformer`, `BeamformingResult`, `DispersionCurve` | f-k processing and its results |
| `statistics` | | density of the maxima, summaries, peaks |
| `response` | `ArrayResponse`, `ResponseAnalysis` | array response and limits; entry point of planning |
| `performance` | `PerformanceResult` | Monte Carlo performance of the estimators |
| `operational` | `OperationalAssessment` | usable frequency-velocity envelope |
| `monitor` | `ArrayMonitor` | real-time processing with the SeedLink client |
| `io` | | array files and result files |
| `cli` | | the `slab-yarray` command |

The numerical modules do no file input/output and no plotting, so
they can be used in any context. The generic spectral tools on which
they rely (windows, Fourier coefficients, cross-spectral matrices) are
in `shakelab.signals` (`windowing`, `crossspectral`, `fourier.dft`),
where they are shared with other array methods. The figures are in
`shakelab.graphics.yarray`; matplotlib is needed only to draw them.

## 6. Conventions

These conventions are specific to array processing; the general ones
of ShakeLab are in [`docs/conventions.md`](../../../docs/conventions.md).

- **Coordinates** are East, North, Up in metres. When they are
  computed from geographic positions, the origin of the local frame is
  the first station added, or a `reference` point given by the user or
  by the array file. Adding stations does not move the origin.
- **Azimuths** are in degrees clockwise from north and give the
  direction towards which the wave travels (not where it comes from).
- **Radial** is the horizontal motion along the propagation azimuth;
  **transverse** is the radial direction rotated by +90 degrees.
- **Ellipticity** (Wathelet et al., 2018) is the ratio of radial to
  vertical amplitude, with the radial positive along the propagation
  and the vertical positive upwards: $`e > 0`$ for prograde and $`e < 0`$
  for retrograde particle motion. The fundamental Rayleigh mode of a
  homogeneous half-space is retrograde. Data with the vertical
  positive downwards give the opposite sign.
- **Beam power** is a power spectral density, in the units of the data
  squared per hertz (Bartlett and Capon); data not corrected for the
  instrument response give counts squared per hertz. MUSIC gives a
  pseudo-power without units.
- **Times** are ShakeLab `Date` objects; the time of a block is the
  centre of its covariance window.

## 7. See also

- The [guide](../../../docs/exploration/yarray.md)
  (`docs/exploration/yarray.md`): theory of the methods, complete
  description of the options and of the file formats, practical
  guidance and current limitations.
- [`examples/exploration/yarray/`](../../../examples/README.md#exploration):
  one script per topic (array data, configuration, estimators,
  components, Rayleigh waves, transforms, uncertainty, results,
  density, command line, real time, planning).
- `shakelab.signals.io`, `shakelab.signals.fdsnws` and
  `shakelab.signals.seedlink`: reading, downloading and receiving the
  data in real time.

## 8. References

Capon, J. (1969). High-resolution frequency-wavenumber spectrum
analysis. Proceedings of the IEEE, 57(8), 1408-1418.

Wathelet, M. (2024). Incoherent noise-induced distortions of Rayleigh
wave ellipticity measurements obtained with three-component
beamforming. Geophysical Journal International, 236(3), 1804-1827.

Wathelet, M., Guillier, B., Roux, P., Cornou, C., and Ohrnberger, M.
(2018). Rayleigh wave three-component beamforming: signed ellipticity
assessment from high-resolution frequency-wavenumber processing of
ambient vibration arrays. Geophysical Journal International, 215(1),
507-523.
