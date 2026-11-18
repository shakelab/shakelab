# Signals: waveforms, metadata and spectra

`shakelab.signals` holds the seismic recordings and everything needed
to work with them: the waveform objects, the description of the
recording channels (positions, orientations, instrument responses),
the reading and writing of the common file formats, the basic
processing, Fourier spectra and time-frequency representations, and
the clients that download data from FDSN web services or receive them
in real time from SeedLink servers.

A recording is a `Record`: a continuous series of samples with its
start time, sampling interval and source identifier
(`NET.STA.LOC.CHA`). The records of one channel form a `Stream`, and
the streams of many channels a `StreamCollection`, which is what the
readers and the clients return. Each record also carries the
metadata of its channel, so that coordinates, orientation and
instrument response travel with the data and are used by the
processing without being passed every time.

The other parts of ShakeLab build on these objects: array processing
(`shakelab.exploration.yarray`), ground-motion and site analyses, and
the real-time monitors.

## 1. Features

- Waveform containers (`Record`, `Stream`, `StreamCollection`) with
  selection by channel patterns and time windows, merging of
  contiguous records and in-place processing at every level.
- Reading and writing of MiniSEED 2 and 3 (libmseed or pure Python),
  SAC, DYNA 1.2 (ESM, ITACA), Simple ASCII, text tables and pickle;
  reading of NI TDMS. Formats recognized from the file content,
  folders, patterns, SDS archives and file name templates.
- Channel metadata (`MetadataCollection`) from StationXML, RESP, the
  text format of the FDSN station service and JSON, attached to the
  records and selected by channel and time.
- Instrument response model following the SEED manual and FDSN
  StationXML: removal and application of the complete response, of
  chosen stages or of the sensitivity only, with water level and
  prefilter, and validation of the metadata.
- Processing: detrend, taper, Butterworth filters, resampling,
  integration and differentiation, time shifts, cuts, convolution and
  correlation.
- Fourier spectra (`Spectrum`) and time-frequency representations
  (`Spectrogram`) with the windowed Fourier, Morlet wavelet and
  Stockwell transforms, calibrated in the units of the data.
- FDSN web services client (`fdsnws`): waveforms, also corrected for
  the response, metadata, events and data availability from any data
  center of the FDSN registry.
- SeedLink client (`seedlink`): real-time data in a buffer, in
  MiniSEED archives and in user callbacks, with automatic resumption
  after interruptions.
- Command-line tools `slab-fdsn` and `slab-seedlink`.

## 2. Quick start

This example reads the three components of a strong-motion station,
attaches the channel metadata, removes the instrument response to get
accelerations in m/s², and writes the result to SAC files:

```python
from shakelab.signals.io import metadata_reader, reader, writer

metadata = metadata_reader("station.xml")
waveforms = reader("data/*.mseed", metadata=metadata)
waveforms.info()

waveforms.detrend("linear")
waveforms.taper(0.05)                                 # seconds
waveforms.deconvolve_response(prefilter=(0.05, 0.1, 30.0, 40.0))
waveforms.filter(highpass=0.1, lowpass=25.0)

record = waveforms["IV.MODE..HNZ"][0]                 # stream, record
print(record.head.units, record.peak)

writer(waveforms, "corrected/{net}.{sta}.{loc}.{cha}.sac")
```

`reader` returns a `StreamCollection`; indexing it by source
identifier gives the `Stream` of that channel, and indexing the stream
gives its records. The processing methods work in place on a record,
a stream or a whole collection. `deconvolve_response` uses the
response of each channel, selected at the time of each record; the
prefilter (four corner frequencies in Hz) limits the band in which the
response is inverted.

## 3. Common tasks

### 3.1 Read and write waveforms

`reader` accepts a file, a wildcard pattern, a folder, a list of them
or a text file listing them. The format of each file is taken from its
extension or, when the extension is missing or ambiguous, from its
first bytes; `format=` imposes one. Channels and a time window can be
selected while reading, and archives organized by channel and day
(the SeisComP Data Structure) are read opening only the files needed:

```python
from shakelab.signals.io import reader, writer

data = reader("data/", select="IV.MODE..HN?",
              starttime="2012-05-20T02:04:00",
              endtime="2012-05-20T02:05:00")
archive = reader("sds/", layout="sds", select="OX.*..HHZ",
                 starttime="2026-11-05", endtime="2026-11-06")

writer(data, "mode.mseed")                       # one file
writer(data, "sac/", format="sac")               # one file per record
writer(data, "{sta}/{sta}.{cha}.{start}.sac")    # template
```

`writer` writes a record, a stream or a collection to one file, to a
folder or to files named by a template with the codes of the channel
and the start time (`{net}`, `{sta}`, `{loc}`, `{cha}`, `{start}`,
`{year}`, `{jday}`...); no file is written if one exists and
`owrite=True` is not given. `formats()` lists the formats with their
extensions and capabilities (Section 5).

### 3.2 Records, streams and collections

A `Record` can also be built directly from an array of samples:

```python
import numpy as np
from shakelab.signals.base import Record, StreamCollection

record = Record(data=np.zeros(1000), delta=0.01,
                time="2026-01-01T00:00:00", sid="XX.TEST..HHZ",
                units="m/s")
record.info()

collection = StreamCollection()
collection.append(record)              # stored as a copy
stream = collection["XX.TEST..HHZ"]
part = record.extract("2026-01-01T00:00:02", "2026-01-01T00:00:05")
first = record[:100]                   # Record of the first 100 samples
vertical = collection.select("*.*.*.HHZ",
                             starttime="2026-01-01T00:00:02")
```

Records of the same channel are kept in their stream in the order
they are appended, and a record contiguous with the last one is joined
to it (within half a sample); `extract` returns a new object with a
time window, `cut` trims in place, and a slice of a record is a new
record of consecutive samples (`record.data` holds the samples as a
NumPy array). `select` returns the channels of a collection that match
a pattern, within an optional time window, as the reader does. The
time axis (`time_axis()`), peak value (`peak`),
duration and end time are available on every record.

### 3.3 Channel metadata

The description of the channels is read into a `MetadataCollection`,
a flat list of channel epochs (one channel in a validity interval)
with their station:

```python
from shakelab.signals.io import metadata_reader

metadata = metadata_reader("station.xml")   # or RESP, FDSN text, JSON
metadata.info()
epoch = metadata.get("IV.MODE..HNZ", "2012-05-20T02:04:00")
print(epoch.location, epoch.orientation, epoch.sample_rate)

waveforms.set_metadata(metadata)            # if not given to the reader
problems = metadata.validate()
```

Each record keeps the epoch of its channel (`record.head.meta`): its
`location`, `orientation` and `response` are attributes of the record.
`validate` reports inconsistencies (overlapping epochs, responses whose
gains or units do not agree, sampling rates different from those of
the response) without correcting them.

### 3.4 Instrument response

The response of a channel is a set of stages (poles and zeros, digital
filters, gains, response tables) evaluated as in the SEED manual and
in FDSN StationXML. Three corrections are available on records,
streams and collections:

```python
waveforms.deconvolve_response()                  # complete response
waveforms.deconvolve_response(waterlevel=40.0,
                              prefilter=(0.05, 0.1, 30.0, 40.0))
waveforms.deconvolve_response(stages=("gain", "paz"))  # analog part
waveforms.remove_sensitivity()                   # flat band only
```

The water level is in dB below the maximum of the response, and limits
the amplification where the response is small. `convolve_response`
applies a response instead (for example to simulate an instrument).
The units of the data follow the output units of the correction.

### 3.5 Processing

```python
record.remove_mean()
record.detrend("linear")
record.taper(1.0)                          # seconds at each end
record.filter(highpass=0.5, lowpass=10.0)  # Butterworth, zero phase
record.resample(0.02)                      # new sampling interval
record.integrate()                         # e.g. velocity from acceleration
record.shift(0.5)                          # seconds
```

All the operations change the data in place and are available on
streams and collections as well, which apply them to every record.
Filters are zero-phase by default (`minphase=True` for a causal
filter); integration and differentiation work in the frequency domain
by default.

### 3.6 Spectra and spectrograms

```python
spectrum = record.to_spectrum()
frequency, amplitude = spectrum.frequency, spectrum.amplitude

spectrogram = record.spectrogram("cwt", fmin=0.2, fmax=20.0)
spectrogram.info()
```

A `Spectrum` is the one-sided Fourier transform of a record, with its
metadata, scaled to approximate the continuous transform (a record in
m/s gives a spectrum in m/s·s); it can be filtered, smoothed, corrected
for the response and transformed back (`to_record`). A `Spectrogram`
holds a time-frequency representation computed with the windowed
Fourier transform (`"wft"`), the Morlet wavelet transform (`"cwt"`) or
the Stockwell transform (`"st"`), with amplitude, power, phase and the
cone of influence of the edge effects; the three methods are
calibrated in the units of the data, so that they can be compared.

### 3.7 Data from FDSN web services

```python
from shakelab.signals.fdsnws import Client

client = Client("INGV")
waveforms = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:03:23",
                                 "2012-05-20T02:06:23", correct=True)
metadata = client.get_metadata("IV.MO*", level="channel")
```

The client reaches any data center of the FDSN registry by name, or a
server by its address; it downloads waveforms (also corrected for the
response), metadata, events and data availability. See the
[`fdsnws` README](fdsnws/README.md).

### 3.8 Real-time data from SeedLink

```python
from shakelab.signals.seedlink import Client

client = Client("rtserve.example.org")
client.select("IV.MODE..HN?")
client.set_buffer(max_duration=600)
client.every(60, lambda client: print(client.buffer))
client.start()
```

The client receives the data in a background thread and sends them to
a buffer, to MiniSEED archives and to user functions, and resumes the
stream after interruptions. See the
[`seedlink` README](seedlink/README.md).

## 4. Command line

Two commands are installed with ShakeLab:

```sh
slab-fdsn waveforms INGV "IV.MODE..HN?" \
    --start 2012-05-20T02:03:23 --end 2012-05-20T02:06:23 \
    -o mode.mseed --metadata mode.xml
slab-fdsn metadata INGV "IV.MO*" --level station

slab-seedlink info rtserve.example.org --streams IV.MODE
slab-seedlink archive rtserve.example.org "IV.*..HH?" \
    --sds archive --state iv.state
```

`slab-fdsn` downloads waveforms, metadata, events and availability;
`slab-seedlink` lists the streams of a server, monitors the latency of
the channels, collects a time window or archives the data
continuously. `--help` after a subcommand lists its options.

## 5. Formats

| Format | Read | Write | Module | Notes |
|---|---|---|---|---|
| MiniSEED 2 and 3 (`mseed`) | yes | yes | `libio.cymseed`, `libio.mseed` | libmseed by default; pure-Python backend for MiniSEED 2 |
| SAC (`sac`) | yes | yes | `libio.sac` | one record per file; header fields kept |
| DYNA 1.2 (`dyna`) | yes | yes | `libio.dyna` | ESM and ITACA; one record per file; values in cm |
| Simple ASCII (`simpleascii`) | yes | yes | `libio.simpleascii` | format of `mseed2ascii` and of the EarthScope services |
| Text tables (`ascii`) | yes | yes | `libio.ascii` | any layout, CSV; parameters in the file or as arguments |
| NI TDMS (`tdms`) | yes | no | `libio.tdms` | DAQ systems and DAS interrogators |
| Pickle (`pickle`) | yes | yes | `io` | intermediate results; trusted files only |
| StationXML | yes | no | `stationxml` | metadata |
| RESP | yes | no | `respfile` | metadata (response only) |
| FDSN station text (`fdsntext`) | yes | no | `stationtext` | metadata without response stages |
| JSON | yes | yes | `metadata` | metadata (`MetadataCollection.write`) |

The common scheme of the format modules, and how to add a format, are
described in the [`libio` README](libio/README.md).

## 6. Structure

The modules are organized in layers, each depending only on those
below it:

| Module | Main objects | Role |
|---|---|---|
| `fourier` | | Fourier primitives on arrays: frequency axes, scaled transforms, Butterworth filters, smoothing |
| `codes` | `SidPattern` | source identifiers and channel patterns |
| `response` | `StageSet`, `StagePoleZero`, `StageCoefficients`... | instrument response model |
| `metadata` | `MetadataCollection`, `ChannelEpoch`, `StationEpoch` | channel and station metadata |
| `header` | `Header` | metadata of a record: sampling, times, identifiers, channel epoch |
| `base` | `Record`, `Stream`, `StreamCollection` | waveform containers and processing |
| `spectrum` | `Spectrum` | Fourier spectrum of a record |
| `spectrogram` | `Spectrogram` | time-frequency representations |
| `windowing` | | division of series into windows and blocks |
| `crossspectral` | `CrossSpectra` | Fourier coefficients of windows and cross-spectral matrices |
| `wavenumber` | | frequency-wavenumber transforms of DAS data |
| `circular` | `CircularStream` | time-limited streams for real-time buffers (experimental) |
| `libio` | | file formats (see its README) |
| `stationxml`, `respfile`, `stationtext` | | metadata formats |
| `io` | `reader`, `writer`, `metadata_reader` | single entry point for files of every format |
| `fdsnws` | `Client` | FDSN web services (see its README) |
| `seedlink` | `Client` | SeedLink client (see its README) |

`windowing` and `crossspectral` are the building blocks of the methods
that average many windows of several channels: they divide the
records into windows and blocks and estimate cross-spectral matrices,
on which the array processing and the frequency domain decomposition
of `shakelab.exploration.yarray` are built, and future methods such as
the spatial autocorrelation will be. They are part of `signals`
because they know nothing about arrays.

Modules import only from lower layers; where a lower object returns a
higher one (`Record.to_spectrum`, `Record.spectrogram`), the import is
made inside the method. Figures are in `shakelab.graphics`; the
objects have no plotting methods.

## 7. Conventions

The general conventions of ShakeLab are in
[`docs/conventions.md`](../../docs/conventions.md); these are specific
to waveforms and metadata.

- **Source identifiers** are `NET.STA.LOC.CHA`, with an empty location
  allowed (`IV.MODE..HNZ`). Patterns select them code by code: `*`,
  `?`, `[ZNE]`, alternatives `{A,B}`, and codes missing at the end
  match anything (`"IV.MODE"` is every channel of MODE).
- **Times** are ShakeLab `Date` objects in UTC, with nanosecond
  resolution; they are given as ISO 8601 strings. Time windows include
  both ends.
- **Sampling** is given by the interval `delta` in seconds (`rate` in
  Hz is its inverse); records are evenly sampled.
- **Units** are declared as written in the files or metadata (`M/S`,
  `CM/S**2`), or left undefined when the format has none (raw counts
  of MiniSEED); response corrections set the units of their output.
- **Coordinates** are `(longitude, latitude, elevation)` with the
  elevation in m positive upwards; the burial depth of a sensor is in
  km positive downwards; orientations are azimuth (clockwise from
  north) and dip (positive downwards) in degrees, as in SEED.
- **Instrument response**: definitions of the SEED manual (v2.4) and
  of FDSN StationXML; the declared sensitivity is authoritative; a
  positive time correction of the datalogger is a time advance; the
  water level is relative to the maximum of the response, in dB.
- **Spectra** are one-sided and scaled by the sampling interval
  (`delta * DFT`), with the kernel `exp(-2πift)`.

## 8. See also

- [Signals overview](../../docs/signals/signals.md): how the parts of
  the package work together, with the map of its documentation.
- Guides in [`docs/signals/`](../../docs/signals/): waveforms and
  processing ([`waveforms.md`](../../docs/signals/waveforms.md)),
  metadata and instrument response
  ([`response.md`](../../docs/signals/response.md)), spectra and
  time-frequency analysis
  ([`spectral.md`](../../docs/signals/spectral.md)), file formats
  ([`io.md`](../../docs/signals/io.md)), FDSN web services
  ([`fdsnws.md`](../../docs/signals/fdsnws.md)) and SeedLink
  ([`seedlink.md`](../../docs/signals/seedlink.md)).
- [`examples/signals/`](../../examples/README.md#signals): scripts on
  reading and writing, metadata and response correction, processing,
  spectra and the network clients.
- READMEs of the sub-packages: [`libio`](libio/README.md),
  [`fdsnws`](fdsnws/README.md), [`seedlink`](seedlink/README.md).
- [`shakelab.exploration.yarray`](../exploration/yarray/README.md):
  array processing built on these objects.

## 9. References

FDSN (2019). FDSN StationXML schema, version 1.1. International
Federation of Digital Seismograph Networks.

IRIS (2012). SEED Reference Manual, version 2.4. Incorporated Research
Institutions for Seismology.

Goldstein, P., and Snoke, A. (2005). SAC availability for the IRIS
community. Incorporated Research Institutions for Seismology Data
Management Center Electronic Newsletter, 7(1).

Torrence, C., and Compo, G. P. (1998). A practical guide to wavelet
analysis. Bulletin of the American Meteorological Society, 79(1),
61-78.

Stockwell, R. G., Mansinha, L., and Lowe, R. P. (1996). Localization of
the complex spectrum: the S transform. IEEE Transactions on Signal
Processing, 44(4), 998-1001.
