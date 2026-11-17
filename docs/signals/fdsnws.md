[ShakeLab documentation](../README.md) › [Signals](signals.md) ›
FDSN web services

# FDSN web services

Almost every seismological data center publishes its data through the
web services defined by the International Federation of Digital
Seismograph Networks (FDSN): waveforms, station and channel metadata,
earthquake catalogues and the time spans of the available data, each
behind a standard address and with standard parameters. The package
`shakelab.signals.fdsnws` is a client of these services: a `Client`
talks to one data center and returns ShakeLab objects, so that
downloaded data are used exactly as data read from files. Waveforms
can be downloaded with the metadata of their channels and already
corrected for the instrument response. The client uses only the
Python standard library.

This guide describes the services and the concepts of the client
(Part I), the use of every method and of the command `slab-fdsn`
(Part II), and the design choices, validation and limits (Part III).
The objects returned are described in [`waveforms.md`](waveforms.md)
(records and collections), [`response.md`](response.md) (metadata,
responses and their correction) and [`io.md`](io.md) (writing the
downloaded data to files).

The examples download from INGV the three components of the
strong-motion station IV.MODE during the earthquake of 20 May 2012 in
Emilia (raw counts at 100 Hz, from 02:03:23 to 02:06:23 UTC), its
metadata, and the catalogue of the sequence from USGS. The numbers in
the comments are those of these data.

## Contents

Part I. Services and concepts

1. [The FDSN web services](#1-the-fdsn-web-services)
2. [Data centers and addresses](#2-data-centers-and-addresses)
3. [From a ShakeLab selection to a request](#3-from-a-shakelab-selection-to-a-request)
4. [Corrected downloads](#4-corrected-downloads)
5. [Long windows](#5-long-windows)
6. [Transport](#6-transport)
7. [Conventions](#7-conventions)

Part II. Use

8. [Client](#8-client)
9. [Waveforms](#9-waveforms)
10. [Metadata](#10-metadata)
11. [Events](#11-events)
12. [Availability](#12-availability)
13. [Errors and warnings](#13-errors-and-warnings)
14. [Command line](#14-command-line)
15. [Workflows](#15-workflows)

Part III. Design, validation and limits

16. [Design choices](#16-design-choices)
17. [Validation](#17-validation)
18. [Current limitations](#18-current-limitations)
19. [References](#19-references)

# Part I. Services and concepts

## 1. The FDSN web services

The FDSN specifications (version 1.2) define four services used by
the client, each under a standard path of the host of the data
center:

| Service | Path | Method of the client | Reply | ShakeLab object |
|---|---|---|---|---|
| `dataselect` | `/fdsnws/dataselect/1/` | `get_waveforms` | MiniSEED | `StreamCollection` |
| `station` | `/fdsnws/station/1/` | `get_metadata` | StationXML, or text | `MetadataCollection` |
| `event` | `/fdsnws/event/1/` | `get_events` | QuakeML | `EqDatabase` (`shakelab.seismicity`) |
| `availability` | `/fdsnws/availability/1/` | `get_availability` | JSON | list of `TimeSpan` |

Each service has the method `query` (the request itself, with GET
parameters or a POST body), `version` and documentation pages; the
availability service has also `extent`. A request selects data by
codes (network, station, location, channel, with the wildcards `*`
and `?` and lists separated by commas), by time, and for metadata and
events by geography.

Two replies need a word. Status **204** (or 404, if the request asks
for it) means *no data*: the request was valid but nothing matches,
and the client returns an empty object. Errors (400 for a bad request,
413 for a request too large, 5xx for a failure of the server) carry a
short text from the server, which the client reports. Servers limit
the size of a request (the number of channels times the duration) and
refuse larger ones; how large depends on the data center.

## 2. Data centers and addresses

A client is created from the name of a data center. The FDSN keeps a
registry of the data centers and of the addresses of their services
(https://www.fdsn.org/datacenters/); ShakeLab has a built-in copy of
it, with 34 data centers, so that no request is needed to find an
address:

- names are those of the registry, matched without regard to case;
- former names are accepted as aliases: `IRIS` and `IRISDMC` for
  EarthScope, `ORFEUS` for ODC, `RESIF` for EPOSFR;
- the services of a data center may live on different hosts or under
  non-standard paths: each service has its own address;
- a few services missing from the registry are added: the event
  service of INGV and the servers of OGS (`OGS-ANT`, `OGS-SCP`).

`datacenters(refresh=True)` reads the registry again and uses it for
the rest of the session.

Any other server can be used, not only those of the registry: the
FDSN services of an observatory, of a temporary network, or the local
server of an acquisition system (SeisComP, for instance, offers them
on its own port). It is given by the `http://` or `https://` address
of its host, with an optional port and path prefix
(`http://192.168.1.10:8080`, `https://example.org/seismo`), when it
uses the standard paths `/fdsnws/<service>/1/`; otherwise, or for
services on different hosts, by the address of each service. The OGS
servers of the table are themselves addresses with a port.

## 3. From a ShakeLab selection to a request

Channels are selected with the patterns of ShakeLab,
`NET.STA.LOC.CHA` with `*`, `?`, sets `[...]` and alternatives
`{A,B}`, where codes missing at the end match anything (see
[`waveforms.md`](waveforms.md) §8.3). The services know only `*`, `?`
and lists, so each pattern is translated:

| ShakeLab pattern | Codes sent to the server |
|---|---|
| `IV.MODE..HN?` | `IV`, `MODE`, `--`, `HN?` (an empty location is written `--`) |
| `IV.MODE..HN[ZNE]` | `IV`, `MODE`, `--`, `HNZ,HNN,HNE` (a set becomes a list) |
| `IV.MODE..HN[!Z]` | `IV`, `MODE`, `--`, `HN?` (a negated set becomes `?`) |
| `IV.{MODE,CAVE}..HN?` | two requests, one per alternative |
| `IV.MODE` | `IV`, `MODE`, `*`, `*` |

Sets that would expand into more than 100 words also become `?`. The
server may then return more than asked; the client keeps only the
channels, stations and spans that match the original patterns.

A selection with one line is sent with GET, the parameters in the
address. Several lines, or patterns with different time windows, are
sent with POST, as a body with the parameters first and then one line
per combination of codes:

```text
quality=D
IV MODE -- HN? 2012-05-20T02:03:00.000000 2012-05-20T02:07:00.000000
IV CAVE -- HNZ 2012-05-29T07:00:00.000000 2012-05-29T07:05:00.000000
```

Parameter values are formatted as the services expect: times in UTC
with microseconds and without zone designator, real numbers in
decimal notation (never exponential), booleans as `true` or `false`;
parameters left None are not sent, so that the defaults of the server
apply.

## 4. Corrected downloads

`get_waveforms(..., correct=True)` returns records corrected for the
instrument response. The correction is a division in the frequency
domain, which spreads the effects of the edges of a record (the
transient of the filter, the wrap-around of the discrete transform)
over its first and last seconds ([`response.md`](response.md) §5).
The client places those effects outside the requested window:

1. it requests a margin `pad` before and after each window, by default
   10% of the longest window and at least 30 s;
2. it requests the metadata of the channels, with their responses, and
   attaches them;
3. on each padded record it removes a linear trend, tapers the part
   of the record outside the window (a half cosine at each end, at
   most as long as the margin), corrects the response with the water
   level (default 60 dB) and the optional `prefilter`, with zero
   padding to avoid the wrap-around, and applies the optional
   `highpass` and `lowpass` filters;
4. it cuts the records to the requested window.

The same margin, detrend and taper are used when only filters are
asked. Records whose channel has no response are left out, with one
`FdsnWarning`.

**Margins without data.** The margin protects the window only where
the server has data in it. At the beginning or end of the available
data, and at the edges of a gap inside the window, a record has an
end inside the window, or close to it. The taper never covers more
of a record than the part outside the window, so that the data of the
window are not attenuated, but an end needs at least a short taper:
after the detrend a record does not begin at zero, and the correction
of the digital stages, whose response falls near the Nyquist
frequency, amplifies the step up to the water level. The taper is
therefore never shorter than 10 samples (0.1 s at 100 Hz), and where
the record has less than that outside the window it enters the window
by the difference. On the vertical acceleration of IV.MODE, corrected
without prefilter in a window that begins 4 s before the peak
acceleration, with no data before the window (amplitude relative to
the same window corrected with a full margin of 30 s):

| Taper at the beginning of the window | 0-0.1 s | 0.1-0.2 s | 1-2 s | 2-3 s | Peak (reference 0.3264 m/s²) |
|---|---|---|---|---|---|
| as long as the margin, inside the window | 0.00 | 0.00 | 0.01 | 0.02 | 0.085 |
| none | 23.8 | 34.9 | 1.75 | 1.42 | 1.260 |
| 10 samples (the client) | 0.67 | 1.04 | 1.00 | 1.00 | 0.3261 |

With the prefilter `(0.05, 0.1, 30, 40)` the step is not amplified:
without taper the window is correct from the first sample, with 10
samples from the first tenth of a second, while the long taper still
attenuates it for 30 s. The client reports the records
without the whole margin in one `FdsnWarning`, with the side
(`before`, `after`): the first or last seconds of the window may still
carry the transient of the correction at low frequencies, which no
short taper removes (sensors whose response falls at low frequencies,
§18).

```python
from shakelab.signals.fdsnws import Client

client = Client("INGV")
data = client.get_waveforms("IV.MODE..HNZ", "2012-05-20T02:03:23",
                            "2012-05-20T02:05:30", correct=True)
```

On the test data, which hold only these three minutes, the call warns
that there are no data "in the whole margin of the window for
IV.MODE..HNZ (before and after)".

A longer `pad` widens the protection where the server has data, at the
cost of a larger download; it cannot add data before the first sample
or inside a gap.

## 5. Long windows

Servers refuse requests larger than their limit, and a long download in
one reply is lost entirely if the connection breaks. `chunk` splits
each window into requests of at most `chunk` seconds, sent one after
the other.

The servers return whole MiniSEED data records, a few seconds long
each, that intersect the requested window, so the replies of
consecutive chunks overlap at their boundaries. The client cuts each
reply to its chunk (excluding the sample at the end, which belongs to
the next chunk) and joins the pieces again: the result is the same as
with one request, contiguous records of a channel joined into one
record.

With `target`, the MiniSEED replies of the chunks are written one after
the other into one file, through a temporary file renamed at the end,
so that no partial file is left if a request fails.

## 6. Transport

The requests use `urllib` of the standard library:

- replies are requested compressed (gzip) and decompressed while they
  are read, also when written to a file block by block;
- **retries**: timeouts, connection failures and the statuses 429
  (too many requests), 500, 502, 503 and 504 are retried, `retries`
  times (default 3), after waits of 1, 2, 4... seconds (at most 30),
  or after the time given by the server in `Retry-After`;
- other statuses raise `FdsnError` at once, with the status, the
  address and the first line of the message of the server;
- `timeout` (default 120 s) is the time without data before a request
  is abandoned, not the duration of the whole download;
- a POST redirected with status 307 or 308 is sent again, body
  included, to the new address;
- certificates are verified with the files named by the environment
  variables `SSL_CERT_FILE` or `SSL_CERT_DIR` (useful behind a proxy
  that inspects the traffic), otherwise with the bundle of `certifi` if
  installed, otherwise with those of the system; verification is never
  disabled;
- requests identify themselves as `ShakeLab/<version>
  (https://shakelab.org)`.

## 7. Conventions

- **Times** are ISO 8601 strings or `Date` objects, in UTC; window
  bounds are included.
- **Geography** follows ShakeLab, not the order of the parameters of
  the services: `region=(lonmin, lonmax, latmin, latmax)` and
  `center=(lon, lat)` in degrees, `radius=(rmin, rmax)` in degrees from
  the center, `depth=(min, max)` in km positive downwards,
  `magnitude=(min, max)`. Either bound of a range may be None.
- **Selections** are patterns `NET.STA.LOC.CHA`, a list of them, or
  `(pattern, starttime, endtime)` tuples with their own windows.
- **Files** are not overwritten unless `owrite=True`.

# Part II. Use

## 8. Client

```python
from shakelab.signals.fdsnws import Client, datacenters

client = Client("INGV")
print(client)
# Client('INGV')
print(client.services["station"])
# https://webservices.ingv.it/fdsnws/station/1/
print(Client("IRIS"))                      # a former name
# Client('EarthScope')
print(sorted(Client("USGS").services))     # an event service only
# ['event']
Client("XYZ")                              # ValueError: unknown data center
```

Any server, also one not in the registry, is given by the address of
its host, with the port and a path prefix if needed, when it uses the
standard paths; otherwise service by service (the address of a service
may include the method, which is removed). Explicit services replace
those of a named data center. The scheme `http://` or `https://` is
required, and addresses without a host or with an invalid port are
rejected when the client is created; a name that looks like an address
(with `.`, `:` or `/`) gives an error that suggests the scheme:

```python
local = Client("http://192.168.1.10:8080")      # e.g. a SeisComP server
print(local.services["dataselect"])
# http://192.168.1.10:8080/fdsnws/dataselect/1/
print(sorted(local.services))
# ['availability', 'dataselect', 'event', 'station']
prefixed = Client("https://example.org/seismo")
print(prefixed.services["station"])
# https://example.org/seismo/fdsnws/station/1/
mixed = Client(station="https://meta.example.org/fdsnws/station/1/query",
               dataselect="http://192.168.1.10:8080/custom/dataselect/1/")
print(mixed)
# Client(services=['dataselect', 'station'])
inherited = Client("INGV", station="http://localhost:8080/fdsnws/station/1/")
print(inherited.services["station"], inherited.services["dataselect"])
# http://localhost:8080/fdsnws/station/1/
# https://webservices.ingv.it/fdsnws/dataselect/1/
Client("192.168.1.10:8080")         # ValueError: unknown data center; write 'http://192.168.1.10:8080'
Client("http://")                   # ValueError: service address without host
patient = Client("INGV", timeout=300, retries=5)
```

`timeout` is the time without data before a request is abandoned and
`retries` the number of further attempts (§6). `version(service)`
returns the version of a service, and `datacenters()` the table of the
registry:

```python
print(client.version("station"))
# 1.1.0                       the version reported by the server
centers = datacenters()
print(len(centers), centers["GEOFON"].services["dataselect"])
# 34 https://geofon.gfz.de/fdsnws/dataselect/1/
centers = datacenters(refresh=True)        # needs network: the registry
```

A method of a service that the data center does not offer raises a
`ValueError` before any request.

## 9. Waveforms

### 9.1 A time window

```python
data = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:04:00",
                            "2012-05-20T02:05:00")
print(data.sid)
# ['IV.MODE..HNE', 'IV.MODE..HNN', 'IV.MODE..HNZ']
print(data["IV.MODE..HNZ"][0])
# starttime 2012-05-20T02:04:00.000000Z, duration 60.0s,
# sampling rate 100.0 Hz, 6001 samples
```

The result is a `StreamCollection` (empty if the server has no data),
with the records cut to the window: the servers return whole MiniSEED
data records, which usually begin before the window and end after it.
Every pattern needs a time window.

### 9.2 Several windows

A selection can give each pattern its own window, as a tuple; patterns
without a window take that of the arguments. Several windows are sent
in one POST request:

```python
event = ("2012-05-20T02:03:40", "2012-05-20T02:04:40")
data = client.get_waveforms(
    [("IV.MODE..HNZ", *event), "IV.MODE..HNE"],
    "2012-05-20T02:05:00", "2012-05-20T02:05:30")
for record in data.records():
    print(record.sid, record.time, record.duration)
# IV.MODE..HNE 2012-05-20T02:05:00.000000Z 30.0
# IV.MODE..HNZ 2012-05-20T02:03:40.000000Z 60.0
client.get_waveforms("IV.MODE..HNZ")       # ValueError: no time window
```

### 9.3 With their metadata

`metadata=True` also requests the metadata of the channels, at
response level, and attaches them to the records, as `metadata=` does
when reading files:

```python
data = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:03:23",
                            "2012-05-20T02:06:23", metadata=True)
record = data["IV.MODE..HNZ"][0]
print(record.location)
# WgsPoint(longitude=10.949167, latitude=44.629722, elevation=41.0)
data.detrend()                       # later, as for data read from files
data.taper(5.0)
data.deconvolve_response(waterlevel=40.0)
```

Channels without metadata keep their records, with one `FdsnWarning`.

`metadata` also takes metadata already at hand, as a
`MetadataCollection` or StationXML files: they are attached instead of
the reply of the server, and used by `correct`, without a request to
the station service. This corrects the data with metadata of your own,
for example a StationXML fixed by hand, or with metadata downloaded
once for several requests:

```python
from shakelab.signals.io import metadata_reader

network = metadata_reader("station.xml")
data = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:04:00",
                            "2012-05-20T02:05:00", metadata=network,
                            correct=True)
print(round(data["IV.MODE..HNZ"][0].peak, 4))
# 0.3264                           without prefilter
```

Channels missing from the given metadata behave as channels that the
server has no metadata for: kept with a warning, or left out with a
warning when they are to be corrected.

### 9.4 Corrected and filtered

```python
acceleration = client.get_waveforms(
    "IV.MODE..HNZ", "2012-05-20T02:04:00", "2012-05-20T02:05:00",
    correct=True, prefilter=(0.05, 0.1, 30.0, 40.0))
record = acceleration.records()[0]
print(record.head.units, round(record.peak, 4))
# m/s**2 0.3258
print(record.duration)
# 60.0                             the margins have been cut
```

The steps are described in §4. The options:

| Option | Meaning | Default |
|---|---|---|
| `correct` | correct the instrument response (implies `metadata=True` when no metadata are given) | False |
| `pad` | margin before and after each window, s | 10% of the window, at least 30 s |
| `prefilter` | `(f1, f2, f3, f4)`, frequency window of the correction, Hz (with `correct` only) | none |
| `waterlevel` | water level of the correction, dB below the maximum | 60 |
| `highpass`, `lowpass` | corners of a Butterworth filter applied after the correction, Hz, as `Record.filter` | none |

Filters can also be applied without the correction, to the raw counts:

```python
window = ("2012-05-20T02:04:00", "2012-05-20T02:05:00")
filtered = client.get_waveforms("IV.MODE..HNZ", *window, highpass=1.0)
client.get_waveforms("IV.MODE..HNZ", *window, prefilter=(0.05, 0.1, 30, 40))  # ValueError
```

`prefilter` without `correct`, and `pad` without correction or
filters, raise a `ValueError`.

### 9.5 Long windows

```python
day = client.get_waveforms("IV.MODE..HNZ", "2012-05-20T02:03:23",
                           "2012-05-20T02:06:23", chunk=60.0)
print(len(day.records()), day.records()[0].nsamp)
# 1 18001                     three requests of 60 s, one record
```

`chunk` is the longest window of one request, in seconds (§5): one
hour (`chunk=3600`) is a reasonable choice for a day of 100 Hz data
of a few channels. The result does not depend on it.

### 9.6 Options of the server

`quality` asks for a data quality code (`D`, `R`, `Q`, `M`, or `B` for
the best available), `minimumlength` drops continuous segments shorter
than the given seconds, and `longestonly=True` keeps only the longest
segment of each channel; the server applies them. `mseed_backend`
chooses the MiniSEED reader, as in `io.reader`.

```python
data = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:04:00",
                            "2012-05-20T02:05:00", quality="D",
                            minimumlength=30.0, mseed_backend="python")
```

### 9.7 Saving the MiniSEED as received

`target` writes the MiniSEED records to a file exactly as the server
sends them, without decoding them, and returns the number of bytes
written; it cannot be combined with metadata, corrections or filters,
which need the records decoded:

```python
size = client.get_waveforms("IV.MODE..HN?", "2012-05-20T02:03:23",
                            "2012-05-20T02:06:23", target="mode.mseed")
print(size)
# 83456
whole = ("2012-05-20T02:03:23", "2012-05-20T02:06:23")
client.get_waveforms("IV.MODE..HN?", *whole, target="mode.mseed")  # FileExistsError
```

The file holds whole data records, so it usually covers a little more
than the window. Corrected or filtered data are written with
`io.writer` in any format ([`io.md`](io.md) §9).

### 9.8 Adding to a collection

```python
from shakelab.signals.base import StreamCollection

collection = StreamCollection()
for window in (("2012-05-20T02:04:00", "2012-05-20T02:04:30"),
               ("2012-05-20T02:05:00", "2012-05-20T02:05:30")):
    client.get_waveforms("IV.MODE..HNZ", *window,
                         stream_collection=collection)
print(len(collection["IV.MODE..HNZ"]))
# 2                                  two records, with a gap
```

Metadata and corrections are applied only to the new records.

## 10. Metadata

`get_metadata` returns a `MetadataCollection`
([`response.md`](response.md) §8). `level` sets the detail: station
epochs, channel epochs (the default), or channel epochs with their
responses:

```python
stations = client.get_metadata("IV.MO*", level="station")
print(stations)
# MetadataCollection(stations=1, channels=0, epochs=0)
channels = client.get_metadata("IV.MODE..HN?",
                               "2012-05-20", "2012-05-21")
print(channels.sid)
# ['IV.MODE..HNZ', 'IV.MODE..HNN', 'IV.MODE..HNE']    order of the reply
responses = client.get_metadata("IV.MODE..HN?", "2012-05-20",
                                "2012-05-21", level="response")
```

With a time window, only the epochs active in it are returned.
`text=True` asks for the text format of the station service, much
lighter than StationXML: coordinates, orientation, sampling rate and
overall sensitivity of each channel, without the response stages (and
not at response level):

```python
light = client.get_metadata("IV.*..HN?", text=True)
print(light.sid)
# ['IV.MODE..HNE', 'IV.MODE..HNN', 'IV.MODE..HNZ']
client.get_metadata("IV.MODE", level="response", text=True)  # ValueError
```

Stations and channels can be selected by geography, by the dates of
their epochs, and by other constraints of the service:

```python
around = client.get_metadata("IV.*..HN?", level="station",
                             center=(11.26, 44.90), radius=(None, 0.5))
inside = client.get_metadata("*", level="station",
                             region=(9.0, 13.0, 43.5, 46.0))
recent = client.get_metadata("IV.*", level="station",
                             startafter="2020-01-01",
                             matchtimeseries=True)
client.get_metadata("IV.*", radius=(None, 0.5))   # ValueError: needs center
```

| Option | Meaning |
|---|---|
| `region` | `(lonmin, lonmax, latmin, latmax)`, degrees |
| `center`, `radius` | `(lon, lat)` and `(rmin, rmax)` in degrees from it |
| `startbefore`, `startafter`, `endbefore`, `endafter` | constraints on the start and end of the epochs |
| `includerestricted` | whether restricted channels are listed (server default: yes) |
| `matchtimeseries` | only channels with data in the archive |
| `updatedafter` | only metadata updated after a time |

`target` writes the reply (StationXML or text) to a file, as received,
and returns its size; the file is read later with `metadata_reader`
([`io.md`](io.md) §11.3):

```python
client.get_metadata("IV.MODE..HN?", level="response", target="mode.xml")
```

## 11. Events

`get_events` returns an earthquake catalogue as an `EqDatabase` of
`shakelab.seismicity`:

```python
usgs = Client("USGS")
catalogue = usgs.get_events("2012-05-01", "2012-06-30",
                            region=(9.0, 13.0, 43.5, 46.0),
                            magnitude=(5.5, None))
for event in catalogue.event:
    origin = event.location.prime
    size = event.magnitude.prime
    print(origin.date, origin["Longitude"], origin["Latitude"],
          origin["Depth"], size["MagSize"], size["MagType"])
```

Each event has its origins (`event.location`) and magnitudes
(`event.magnitude`), the preferred one in `prime`; the loop prints
the time, the epicentre, the depth in km, the magnitude and its type
of each event.

The constraints:

| Option | Meaning |
|---|---|
| `starttime`, `endtime` | origin time window |
| `region`, `center`, `radius` | geography, as for metadata |
| `depth` | `(min, max)` in km, positive downwards |
| `magnitude`, `magnitudetype` | `(min, max)` magnitude, and its type (`"ML"`, `"Mw"`...) |
| `eventtype`, `catalog`, `contributor`, `eventid` | type of event, catalogue, contributor, one event |
| `limit`, `offset` | number of events, and index of the first (from 1) |
| `orderby` | `"time"` (default of the services), `"time-asc"`, `"magnitude"`, `"magnitude-asc"` |
| `includeallorigins`, `includeallmagnitudes`, `includearrivals` | whether to include every origin, magnitude and arrival |
| `updatedafter` | only events updated after a time |

```python
latest = usgs.get_events(center=(11.26, 44.90), radius=(None, 1.0),
                         magnitude=(4.0, None), orderby="time", limit=5)
usgs.get_events(depth=(20.0, 10.0))         # ValueError: max below min
usgs.get_events("2012-05-01", "2012-06-30", target="emilia.xml")
```

`target` writes the QuakeML reply to a file. A catalogue without
events is an empty `EqDatabase`.

## 12. Availability

`get_availability` returns the time spans of the data held by the
data center, as `TimeSpan` tuples (`sid`, `starttime`, `endtime`,
`quality`, `sample_rate`):

```python
spans = client.get_availability("IV.MODE..HN?", "2012-05-20",
                                "2012-05-21")
for span in spans:
    print(span.sid, span.starttime, span.endtime, span.quality)

extent = client.get_availability("IV.MODE..HN?", "2012-05-20",
                                 "2012-05-21", extent=True)
# one span per channel: its first and last sample
```

`extent=True` uses the method `extent`: one span per channel, from the
first to the last sample, gaps ignored. `merge` merges spans across
`"samplerate"`, `"quality"` and (not with `extent`) `"overlap"`;
`mergegaps` ignores gaps shorter than the given seconds; `quality` and
`includerestricted` select the data. Spans are returned in the order
of the server. Not every data center offers the service.

## 13. Errors and warnings

| Situation | Result |
|---|---|
| unknown data center, invalid address | `ValueError` when creating the client |
| service not offered by the data center | `ValueError` before the request |
| wrong argument types, invalid ranges, incompatible options | `TypeError` or `ValueError` before the request |
| no data for the request (status 204) | empty object |
| server error (400, 413...) | `FdsnError` with `status`, `url` and `detail` (the message of the server); for 413, a hint to use `chunk` |
| timeout, connection failure, busy server (429, 5xx) | retried; `FdsnError` after the last attempt |
| certificate not verified | `FdsnError` explaining how to name the certificates |
| `target` exists and `owrite` is False | `FileExistsError`, before the request |
| channels without metadata (`metadata=True`) | records kept, `FdsnWarning` |
| channels without response (`correct=True`) | records left out, `FdsnWarning` |

```python
from shakelab.signals.fdsnws import FdsnError

try:
    data = client.get_waveforms("IV.MODE..HNZ", "2012-05-20T02:04:00",
                                "2012-05-20T02:05:00")
except FdsnError as error:
    print(error.status, error.detail)
```

An empty result is not an error: check `len(data)` (or the size
returned with `target`, which is 0).

## 14. Command line

`slab-fdsn` (also `python -m shakelab.signals.fdsnws`) gives the four
services and the table of data centers from a terminal. The data
center is a name of the registry or the address of a host with the
standard paths (`slab-fdsn waveforms http://192.168.1.10:8080 ...`;
services with other paths need the library); patterns and times are
those of the library:

```sh
slab-fdsn waveforms INGV "IV.MODE..HN?" \
    --start 2012-05-20T02:03:23 --end 2012-05-20T02:06:23 \
    -o mode.mseed --metadata mode.xml
slab-fdsn waveforms INGV "IV.MODE..HN?" \
    --start 2012-05-20T02:04:00 --end 2012-05-20T02:05:30 \
    --correct --prefilter 0.05 0.1 30 40 --lowpass 20 \
    -o "{net}.{sta}.{loc}.{cha}.sac"
slab-fdsn metadata INGV "IV.MO*" --level station
slab-fdsn metadata INGV "IV.MODE..HN?" --level response -o mode.xml
slab-fdsn events USGS --start 2012-05-01 --end 2012-06-30 \
    --minmag 5.5 --region 9 13 43.5 46
slab-fdsn availability INGV "IV.MODE..HN?" --start 2012-05-20 \
    --end 2012-05-21 --extent
slab-fdsn datacenters --service availability --urls
```

| Subcommand | Options | Output |
|---|---|---|
| `waveforms` | patterns, `--start`, `--end`, `-o` (required), `--metadata FILE`, `--correct`, `--pad`, `--prefilter F1 F2 F3 F4`, `--waterlevel`, `--highpass`, `--lowpass`, `--chunk`, `--quality`, `--format`, `--owrite` | MiniSEED as received; with `--correct` or filters, any format of `io.writer`, also with file name templates and folders |
| `metadata` | patterns (default all), `--start`, `--end`, `--level`, `--text`, `-o`, `--owrite` | the reply written to a file, or a table of stations or channels |
| `events` | `--start`, `--end`, `--minmag`, `--maxmag`, `--magtype`, `--mindepth`, `--maxdepth`, `--region`, `--center`, `--maxradius`, `--limit`, `--orderby`, `-o`, `--owrite` | QuakeML written to a file, or a table of the events |
| `availability` | patterns, `--start`, `--end`, `--extent`, `--merge`, `--mergegaps` | a table of the spans |
| `datacenters` | `--refresh`, `--service`, `--urls` | the data centers and their services |

With `--metadata FILE` and `--correct`, the StationXML written to the
file is also the one used for the correction (one request to the
station service). All subcommands except `datacenters` accept
`--timeout`. The summary
and the warnings go to the standard error; the exit status is 0 on
success, 1 on errors and when there are no data, 130 when
interrupted.

## 15. Workflows

### 15.1 The stations around an epicentre, corrected

```python
epicentre = (11.26, 44.90)
stations = client.get_metadata("IV.*..HN?", "2012-05-20", "2012-05-21",
                               level="station", center=epicentre,
                               radius=(None, 0.5))
codes = [f"{station.code}..HN?" for station in stations.stations()]
print(codes)
# ['IV.MODE..HN?']               stations within 0.5° (about 55 km)

data = client.get_waveforms(codes, "2012-05-20T02:03:40",
                            "2012-05-20T02:05:40", correct=True,
                            prefilter=(0.05, 0.1, 30.0, 40.0))
for record in data.records():
    print(record.sid, round(record.peak, 3), record.head.units)
# IV.MODE..HNE 0.27 m/s**2
# IV.MODE..HNN 0.418 m/s**2
# IV.MODE..HNZ 0.326 m/s**2
```

The peaks agree within 0.5% with those of the accelerograms processed
by ESM (0.270, 0.419 and 0.327 m/s²).

### 15.2 A day of data into an archive

```python
from shakelab.signals.io import writer

day = client.get_waveforms("IV.MODE..HN?", "2012-05-20",
                           "2012-05-20T23:59:59.99", chunk=3600)
writer(day, "archive/{year}/{net}/{sta}/{cha}.D/"
            "{net}.{sta}.{loc}.{cha}.D.{year}.{jday}", format="mseed")
```

The archive is then read with `reader(..., layout="sds")`
([`io.md`](io.md) §8.5).

### 15.3 The metadata of a network, kept as JSON

```python
network = client.get_metadata("IV.MODE..HN?", level="response")
network.write("mode.json")
```

The JSON of ShakeLab keeps the responses and is read back with
`metadata_reader` or `MetadataCollection`, without a new request.

### 15.4 Checking the availability before downloading

```python
spans = client.get_availability("IV.MODE..HN?", "2012-05-20",
                                 "2012-05-21")
for span in spans:
    data = client.get_waveforms(span.sid, span.starttime, span.endtime,
                                chunk=3600)
```

# Part III. Design, validation and limits

## 16. Design choices

- **ShakeLab objects, ShakeLab conventions.** Waveforms, metadata and
  events come back as the objects used for files, so that downloaded
  data need no special treatment; selections use the patterns of
  ShakeLab and geography the order (longitude, latitude) of ShakeLab,
  not the parameter names of the services.
- **The client filters what it receives.** Patterns that the services
  cannot express are widened in the request and narrowed again on the
  reply, so that the result is always what the pattern says.
- **Corrections where the data are complete.** The margin, the taper
  and the zero padding place the edge effects of the correction
  outside the window, which the client cuts at the end; the parameters
  are those of `Record.deconvolve_response`. The taper does not enter
  the window beyond the 10 samples that an end without margin needs
  (§4), and missing margins are reported.
- **Raw data untouched when asked.** `target` writes the bytes of the
  server, without decoding and encoding them again.
- **A built-in registry.** Creating a client needs no request; the
  table can be refreshed from the FDSN, and any server can be given by
  its address.
- **Standard library only.** No dependency for HTTP; compression,
  retries, redirects and certificates handled explicitly, verification
  never disabled.
- **Errors that explain.** The message of the server is kept in
  `FdsnError`; arguments are checked before any request.
- **One core, several entry points.** `slab-fdsn` calls the same
  methods as the library.

## 17. Validation

The package is covered by the tests of `tests/signals/fdsnws/` (94
tests), which run against a fake FDSN server on a local socket: the
server answers the four services with the test data of ShakeLab
(IV.MODE in MiniSEED, StationXML and text, a QuakeML catalogue and
availability replies), records every request, and can be made to
fail, compress, redirect or answer with chosen statuses.

| Check | Result |
|---|---|
| Translation of patterns (empty location, sets, negated and large sets, alternatives) | the codes of §3; GET for one line, POST bodies line by line |
| Formatting of parameters (times, small and large reals, booleans) | as the FDSN specifications require |
| Waveforms: window, several windows, chunks joined, `target` with chunks | records cut to the windows; chunked and single downloads equal |
| Corrected downloads, with and without prefilter and filters | margins requested; records in m/s², cut to the window; peak acceleration of IV.MODE..HNZ within 0.5% of the accelerogram processed by ESM |
| Corrected downloads without data in the margins: a window from the start of the data, a gap inside the window | peak acceleration within 0.5% of ESM; amplitude next to the gap within 20% of the download without gap; one `FdsnWarning` with the sides; no warning with whole margins |
| Metadata given as a collection or a file, also with channels missing; `slab-fdsn waveforms --metadata FILE --correct` | no request to the station service (one with the command line); the PGA of ESM; missing channels reported |
| Metadata at the three levels, text format, geography | the expected epochs; parameters in the order of the services |
| Events, availability (`query` and `extent`) | the expected catalogue and spans |
| Statuses 204, 400, 413, 429 and 5xx, `Retry-After`, timeouts, gzip, 307 redirects of POST requests | empty results, `FdsnError` with the message of the server, retries |
| Registry: names, aliases, addresses on several hosts, refresh | the expected services |
| Command line: every subcommand, outputs and exit statuses | as documented |

## 18. Current limitations

- **Margins without data.** Where the server has no data in a margin
  (start or end of the available data, gaps), the window is protected
  only by a taper of 10 samples (§4), with an `FdsnWarning`. For
  sensors whose response falls at low frequencies (broadband
  velocimeters, geophones) the low-frequency transient of the
  correction can last longer than the taper: use a `prefilter`, or ask
  for a window that the data exceed by the margin.
- **One data center at a time.** There is no federation (EIDA routing,
  the FDSN federator), no parallel downloads and no authentication for
  restricted data.
- **Size of the requests.** The client does not know the limits of the
  servers: a request refused as too large (413) must be split with
  `chunk`, or by channels.
- **Event catalogue.** `get_events` returns the `EqDatabase` of
  `shakelab.seismicity`, whose interface has not been revised yet; the
  text format of the event service is not read.
- **Interface asymmetries** (under review): `target` exists for
  waveforms, metadata and events, not for availability, and makes the
  methods return a number of bytes instead of an object;
  `stream_collection` exists for waveforms only; the command line
  writes `--minmag`, `--maxmag`, `--mindepth`, `--maxdepth`,
  `--maxradius` for the ranges of the library and has no geographic
  constraints for metadata.

## 19. References

FDSN (2019). FDSN Web Service Specifications, version 1.2.
International Federation of Digital Seismograph Networks.
https://www.fdsn.org/webservices/

FDSN (2023). FDSN Data Availability web service, version 1.0.
International Federation of Digital Seismograph Networks.

FDSN (2016). FDSN StationXML, version 1.1. International Federation of
Digital Seismograph Networks.

Schorlemmer, D., Euchner, F., Kästli, P., and Saul, J. (2011).
QuakeML: status of the XML-based seismological data exchange format.
Annals of Geophysics, 54(1), 59-65.
