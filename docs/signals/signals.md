[ShakeLab documentation](../README.md) › Signals

# Signals

`shakelab.signals` is the foundation of ShakeLab for time series. It
reads seismic and geophysical recordings from files, data centres and
real-time servers, keeps them together with the description of the
channels that recorded them, processes them, and turns them into
spectra and time-frequency representations. Every other part of the
library that works with recordings (array processing, ground-motion
parameters, site response) starts from its objects.

This page tells how the parts of the package work together and where
each of them is documented. A quick start with examples is in the
[README of the package](../../shakelab/signals/README.md).

## How the parts work together

**Data.** A recording is held in a `Record`: a continuous segment of
uniformly sampled samples, with a header that gives its sampling, its
start time and its channel. The records of one channel form a
`Stream`, which keeps them in order of arrival and joins those that
follow each other without gaps; the streams of many channels form a
`StreamCollection`. These three containers are what every reader
returns and what every processing method works on: detrending,
tapering, filtering, resampling, integration, time shifts. They are
described in the [waveforms guide](waveforms.md).

**Metadata and response.** Each record carries the description of its
channel: the station position, the orientation of the sensor, the
sampling rate declared by the network and the instrument response,
read from StationXML, RESP or the text format of the FDSN services.
The description travels with the data through every operation, so that
correcting for the instrument or placing stations on a map needs
nothing more than the record itself. The response model
follows the conventions of SEED and StationXML, stage by stage, and
the correction regularizes the division with a water level and an
optional band-pass window ([response guide](response.md)).

**Spectral representations.** A record is turned into a `Spectrum`,
whose coefficients are scaled to approximate the continuous Fourier
transform (a record in m/s gives a spectrum in m/s·s), or into a `Spectrogram` by a windowed
Fourier, a continuous wavelet or a Stockwell transform. Underneath,
the division of records into windows and the estimation of
cross-spectral matrices between channels are the building blocks of
the multichannel methods, such as the array processing of YARRAY
([spectral guide](spectral.md)).

**Formats.** Files are read and written through one entry point,
`shakelab.signals.io`, which recognizes the format from the content
(MiniSEED, SAC, the DYNA format of ESM and ITACA, Simple ASCII, text
tables, TDMS), reads folders, patterns and SDS archives, selects
channels and time windows, and writes one file per channel with name
templates. Each format is a module of `libio` with the same scheme,
which also allows the direct use of its full header
([input and output guide](io.md)).

**Services.** The FDSN client downloads waveforms, metadata, events
and availability from the data centres, by channel patterns and
geographic regions, correcting the response on request
([FDSN guide](fdsnws.md)). The SeedLink client receives data in real
time into circular buffers, calls user functions on every record or at
regular intervals, and archives the stream; the real-time analyses of
ShakeLab are built on it ([SeedLink guide](seedlink.md)). Both are
also available from the command line, as `slab-fdsn` and
`slab-seedlink`.

## Where to read

| To... | Read |
|---|---|
| build, join, cut, select and process records, streams and collections | [waveforms](waveforms.md) |
| read metadata, correct the instrument response, check a response | [response](response.md) |
| compute spectra, spectrograms, windows and cross-spectra | [spectral](spectral.md) |
| read and write files, use a format directly, add a format | [input and output](io.md), [`libio` README](../../shakelab/signals/libio/README.md) |
| download from data centres | [FDSN web services](fdsnws.md), [`fdsnws` README](../../shakelab/signals/fdsnws/README.md) |
| acquire data in real time | [SeedLink](seedlink.md), [`seedlink` README](../../shakelab/signals/seedlink/README.md) |
| see the package at a glance | [package README](../../shakelab/signals/README.md) |
| run a complete example | [examples](../../examples/README.md) |

The modules of the package, organized in layers, are listed in the
Structure section of the [package README](../../shakelab/signals/README.md#6-structure);
the conventions it shares with the rest of the library are in the
[conventions](../conventions.md).
