# FDSNWS: data from FDSN web services

`shakelab.signals.fdsnws` downloads data from the web services defined
by the International Federation of Digital Seismograph Networks
(FDSN), which almost every seismological data center offers:
waveforms (`dataselect`), station and channel metadata (`station`),
earthquake catalogues (`event`) and the time spans of the available
data (`availability`).

A `Client` talks to one data center, named as in the FDSN registry
(`"INGV"`, `"GEOFON"`, `"EarthScope"`...) or given by its address. Its
methods return ShakeLab objects: waveforms as a `StreamCollection`,
metadata as a `MetadataCollection`, events as a catalogue of
`shakelab.seismicity`. Waveforms can be downloaded with the metadata
of their channels and already corrected for the instrument response.
The client uses only the Python standard library.

## 1. Features

- All the data centers of the FDSN registry by name (with the former
  names `IRIS`, `ORFEUS`, `RESIF` as aliases), or any server by its
  address; services on different hosts or paths handled.
- Waveforms as ShakeLab objects or as MiniSEED files written as
  received; long windows split into several requests.
- Corrected downloads: margins before and after the window, detrend,
  taper, response correction with water level and prefilter, optional
  filter, cut to the window.
- Metadata at station, channel or response level, as
  StationXML or in the light text format.
- Earthquake catalogues with geographic, depth and magnitude
  constraints; data availability.
- Channel selection with the patterns of ShakeLab; geographic
  constraints in the conventions of ShakeLab.
- Retries on timeouts and on temporary server errors, with clear error
  messages.
- Command `slab-fdsn`.

## 2. Quick start

Download the three components of a strong-motion station during an
earthquake, corrected for the instrument response:

```python
from shakelab.signals.fdsnws import Client

client = Client("INGV")
waveforms = client.get_waveforms(
    "IV.MODE..HN?", "2012-05-20T02:03:23", "2012-05-20T02:06:23",
    correct=True, prefilter=(0.05, 0.1, 30.0, 40.0),
)
waveforms.info()                    # accelerations in m/s**2
```

With `correct=True` the client requests the metadata of the channels
and a margin before and after the window (`pad`), removes a linear
trend, tapers the margins, corrects the response in the frequency
domain and cuts the records to the window, so that the edge effects
fall in the discarded margins. Where the server has no data in a
margin (start of the data, gaps) a `FdsnWarning` says so.

## 3. Common tasks

### 3.1 Waveforms

```python
data = client.get_waveforms("IV.MODE..HN?", start, end, metadata=True)
data.deconvolve_response(waterlevel=40.0)  # later, as for files

client.get_waveforms("IV.MODE..HN?", start, end, target="mode.mseed")
day = client.get_waveforms("IV.MODE..HNZ", "2026-11-05", "2026-11-06",
                           chunk=3600)      # one request per hour
```

`metadata=True` attaches the channel metadata to the records, as
`metadata=` does for files; `target` writes the MiniSEED records to a
file unchanged; `highpass` and `lowpass` filter the corrected data.

### 3.2 Metadata

```python
metadata = client.get_metadata("IV.MO*", level="channel")
stations = client.get_metadata("IV.*", level="station",
                               region=(9.0, 13.0, 43.5, 46.0))
client.get_metadata("IV.MODE..HN?", level="response",
                    target="mode.xml")
light = client.get_metadata("IV.*..HN?", text=True)
```

The result is a `MetadataCollection` (see the
[`signals` README](../README.md#33-channel-metadata)). The text format
(`text=True`) is much lighter than StationXML and gives coordinates,
orientations, sampling rates and sensitivities, without response
stages.

### 3.3 Events and availability

```python
usgs = Client("USGS")
catalogue = usgs.get_events("2012-05-01", "2012-06-30",
                            region=(9.0, 13.0, 43.5, 46.0),
                            magnitude=(5.5, None))

spans = client.get_availability("IV.MODE..HN?", "2026-11-01",
                                "2026-11-05")
for span in spans:
    print(span.sid, span.starttime, span.endtime)
```

### 3.4 Data centers

```python
from shakelab.signals.fdsnws import datacenters

for name, center in datacenters().items():
    print(name, center)
```

`datacenters()` returns the built-in table of the registry
(`refresh=True` reads it again from the FDSN). A server not in the
registry is given by its address, `Client("example.org:8080")`
(`https` if the scheme is omitted; servers that serve plain HTTP only
need `http://`, and `port=` replaces the port), or service by service
(`station=`, `dataselect=`, `event=`, `availability=`).

## 4. Command line

```sh
slab-fdsn waveforms INGV "IV.MODE..HN?" \
    --start 2012-05-20T02:03:23 --end 2012-05-20T02:06:23 \
    -o mode.mseed --metadata mode.xml
slab-fdsn waveforms INGV "IV.MODE..HN?" \
    --start 2012-05-20T02:04:00 --end 2012-05-20T02:05:30 \
    --correct --lowpass 40 -o "{net}.{sta}.{loc}.{cha}.sac"
slab-fdsn metadata INGV "IV.MO*" --level station
slab-fdsn events USGS --start 2012-05-01T00:00:00 \
    --end 2012-06-30T00:00:00 --minmag 5.5 --region 9 13 43.5 46
slab-fdsn availability INGV "IV.MODE..HN?" --start 2026-11-01 \
    --end 2026-11-05
slab-fdsn datacenters
```

Waveforms are written as received (MiniSEED) or, when corrected or
filtered, in any format of `shakelab.signals.io` with its file name
templates; metadata and events are written as files or printed as
tables. `python -m shakelab.signals.fdsnws` is the same command.

## 5. Structure

| Module | Main objects | Role |
|---|---|---|
| `transport` | `FdsnError`, `FdsnWarning` | HTTP requests, retries, errors and parameter formatting |
| `registry` | `DataCenter`, `datacenters` | data centers and the addresses of their services |
| `selection` | | channel patterns translated into requests |
| `client` | `Client`, `TimeSpan` | the four services, as ShakeLab objects |
| `cli` | | the `slab-fdsn` command |

The MiniSEED replies are decoded by
`shakelab.signals.io.read_mseed_bytes`, the StationXML and text replies
by `shakelab.signals.stationxml` and `shakelab.signals.stationtext`,
the QuakeML replies by `shakelab.seismicity`.

## 6. Conventions

- Channels are patterns `NET.STA.LOC.CHA` of `shakelab.signals.codes`.
  The services know only `*`, `?` and lists, so the patterns are
  translated (an empty location becomes `--`, `HH[ZNE]` becomes
  `HHZ,HHN,HHE`, alternatives `{A,B}` separate requests), and the client
  keeps only the channels that match the original patterns.
- Times are ISO 8601 strings or ShakeLab `Date` objects, in UTC.
- Geographic constraints follow ShakeLab, not the order of the
  services: `region=(lonmin, lonmax, latmin, latmax)`, `center=(lon,
  lat)`, `radius=(rmin, rmax)` in degrees, `depth=(min, max)` in km
  positive downwards, `magnitude=(min, max)`; either end may be None.
- Nothing is written over an existing file unless `owrite=True`.

## 7. See also

- The guide of the FDSN client
  ([`docs/signals/fdsnws.md`](../../../docs/signals/fdsnws.md)).
- [`examples/signals/`](../../../examples/README.md#signals).
- The [`signals` README](../README.md): processing and response
  correction of the downloaded data.

## 8. References

FDSN (2019). FDSN Web Service Specifications, version 1.2.
International Federation of Digital Seismograph Networks,
https://www.fdsn.org/webservices/.

FDSN. Data Availability web service specification, version 1.0.
International Federation of Digital Seismograph Networks.
