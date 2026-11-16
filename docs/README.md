# ShakeLab documentation

ShakeLab is an open-source Python library for seismology, engineering
seismology and the processing of geophysical data. This page presents
the ideas behind it and its main parts, and is the starting point of
the documentation: every document is linked from the
[map](#map-of-the-documentation) below. Installation and dependencies
are in the [README of the repository](../README.md).

## What ShakeLab is

ShakeLab gives seismologists, geophysicists and engineers a set of
basic tools that work together: waveforms and their metadata,
instrument responses and spectra, array processing, ground-motion
models, site response, seismicity and hazard, exposure and impact.
Each tool is designed as a building block of a common ecosystem rather
than as the solution of a single project, so that the output of one is
the natural input of the next: the records read from a file or a data
centre carry the description of their channel through filtering,
response correction and spectral analysis, into the array processing
or the computation of ground-motion parameters.

The library follows a few principles, which the user meets everywhere:

- **Science first.** The scientific meaning of every operation is
  explicit in its interface: units and conventions are declared,
  relevant choices are parameters rather than hidden defaults, and
  methods are validated against references, published results and
  synthetic data. The [conventions](conventions.md) shared by all the
  library (time, coordinates, units) are written down.
- **Objects and primitives.** Every module works on the objects of
  ShakeLab (records, streams, metadata, spectra, array data) and
  returns objects that are convenient to inspect and to pass on; the
  numerical functions on plain arrays that they are built on remain
  available.
- **One implementation, many ways in.** Each capability exists first
  as a Python interface; the command line tools (`slab-fdsn`,
  `slab-seedlink`, `slab-yarray`), the real-time callbacks of the
  SeedLink client and future services use the same code.
- **Separate concerns.** Data objects hold data and the operations
  that belong to them; file formats are read and written in one place,
  figures are made by `shakelab.graphics`, and the layers of the
  library depend only on those below them.
- **Few dependencies.** NumPy, SciPy and the MiniSEED codec
  `cymseed3`; matplotlib only for figures.

ShakeLab is distributed under the GNU Lesser General Public License,
version 3 or later.

## The subpackages

| Subpackage | Content |
|---|---|
| [`signals`](signals/signals.md) | the foundation for time series: records, streams and collections, channel metadata and instrument responses, processing, spectra and spectrograms, file formats, FDSN web services and SeedLink |
| [`exploration`](exploration/exploration.md) | geophysical exploration from ambient vibrations and seismic arrays: array processing, frequency domain decomposition, survey planning (YARRAY) |
| `groundmotion` | ground-motion models (GMPEs) and ground-motion providers (e.g. ShakeMap grids) |
| `site` | soil models, one-dimensional site response, quarter-wavelength methods, surface waves |
| `hazard` | magnitude-frequency distributions, completeness, magnitude scaling |
| `seismicity` | earthquake catalogues, magnitudes, QuakeML and ISC tools |
| `engineering` | exposure, fragility and taxonomy models, impact scenarios, engineering parameters of records |
| `graphics` | figures for the objects of the other subpackages |
| `libutils` | shared utilities: time, geodesy, JSON and ASCII input and output, time zones |

## Map of the documentation

**Library-wide**

- [Conventions](conventions.md): time, coordinates, units and
  numerical representation adopted across the library.
- [Development guidelines](development.md): architecture, scientific
  semantics and engineering practice, for contributors.
- [Changelog](../CHANGELOG.md): notable changes, deprecations and
  scientific corrections.

**Signals** ([overview](signals/signals.md),
[package README](../shakelab/signals/README.md))

- [Waveforms](signals/waveforms.md): records, streams and collections,
  tolerances, views and copies, processing.
- [Metadata and instrument response](signals/response.md): channel
  metadata, response stages, correction.
- [Spectra and spectrograms](signals/spectral.md): Fourier spectra,
  time-frequency transforms, windows and cross-spectra.
- [Input and output](signals/io.md): file formats, recognition, name
  templates.
- [FDSN web services](signals/fdsnws.md): waveforms, metadata, events
  and availability from data centres.
- [SeedLink](signals/seedlink.md): real-time acquisition.

**Exploration** ([overview](exploration/exploration.md))

- [YARRAY](exploration/yarray.md): seismic-array processing (theory,
  use from Python, the command line and real time, survey planning).

**Examples**

- [Examples](../examples/README.md): short scripts illustrating the use
  of each subpackage.

Every document starts with a line that leads back to this page and to
the overview of its subpackage.
