# SeedLink: real-time waveform data

`shakelab.signals.seedlink` receives waveform data from SeedLink
servers, the protocol with which seismic networks distribute their
data in real time (SeisComP, ringserver and other servers). The data
arrive as MiniSEED records a few seconds after they are recorded, and
the servers usually keep the last hours or days, which can also be
requested.

The `Client` connects to a server, requests the channels selected and
sends every record it receives to one or more destinations: a buffer
in memory, MiniSEED archive files, and functions written by the user
that run at every record, at regular intervals or at connection
events. It reconnects after network failures and resumes each station
from the packet after the last one received, so that no data are lost
while the server keeps them. Three functions cover the simplest uses:
`read` collects a time window, `archive` archives continuously and
`records` iterates over the records as they arrive.

## 1. Features

- Protocol versions 3.1 and 4, negotiated with the server; MiniSEED 2
  and 3 payloads.
- Channel selection with the patterns of ShakeLab (`NET.STA.LOC.CHA`
  with `*`, `?`, `[...]` and `{A,B}`), from a list or a file.
- Destinations combined freely: a buffer (unlimited, or circular with
  a maximum duration), MiniSEED archives (one file, files by channel
  or station, rotation in time, the SeisComP Data Structure) and
  callbacks.
- Real-time and dial-up modes, requests of a time window of the data
  held by the server.
- Automatic reconnection and resumption with the sequence numbers,
  also after a restart with a state file.
- Background thread with thread-safe copies of the buffer, the status
  and latency of every channel, and the stations and streams offered
  by the server.
- Command `slab-seedlink` for information, monitoring, downloads and
  archives.

## 2. Quick start

Collect one minute of data of three channels:

```python
from shakelab.signals import seedlink

data = seedlink.read("rtserve.example.org", "IV.MODE..HH?", duration=60)
data.info()
```

`read` returns a `StreamCollection` with the contiguous records of each
channel joined, ready for the processing of `shakelab.signals`.

## 3. Common tasks

### 3.1 A client with a buffer and periodic functions

```python
from shakelab.signals.seedlink import Client

client = Client("rtserve.example.org")
client.select(["IV.MODE..HH?", "OX.*..HHZ"])
client.set_buffer(max_duration=600)          # last 10 minutes

def check(client):
    data = client.get_data("IV.MODE..HHZ")   # copy of the buffer
    data.info()

client.every(30, check)
client.start()                               # background thread
...
client.stop()
```

The functions registered with `on_record` (called as
`func(client, record)`), `every` (`func(client)`) and `on_event`
(`func(client, event)`) run in the thread that receives the data, so
they should be quick; `every` takes `clock="data"` to follow the time
of the data instead of the wall clock, for example when an archive is
replayed through a server. Other threads read the buffer with
`get_data`, which returns an independent copy (in physical units with
`sensitivity=True`, when the client has the channel metadata), and the
reception with `status` (last sample and latency of every channel). `run` does the
same as `start` in the calling thread, until a duration, an end time
or `stop`.

### 3.2 Archives

```python
client = Client("rtserve.example.org", state_file="iv.state")
client.select("IV.*..HH?")
client.set_archive("archive", layout="sds")  # one file per channel and day
client.run()
```

Records are appended to the files as they arrive, unchanged. A path
template separates the files by network, station or channel
(`"{net}/{sta}/{sta}.{cha}.mseed"`), and `rotate` starts a new file at
regular intervals aligned with the clock (`rotate="1h"`, with
`{time:%Y%m%dT%H}` in the template). With a state file the client
resumes after a restart from the packet after the last one archived.
`archive(host, select, path, ...)` does the same in one call.

### 3.3 Data held by the server

```python
data = seedlink.read("rtserve.example.org", "IV.MODE..HH?",
                     start="2026-11-05T12:00:00",
                     end="2026-11-05T12:05:00")
```

A start time asks the server for the data it still holds from that
time; with an end time in the past the client stops when the server has
sent them (dial-up mode, `fetch=True`).

### 3.4 Records as they arrive

```python
for record in seedlink.records("rtserve.example.org", "IV.MODE..HHZ"):
    print(record.sid, record.starttime, record.nsamp)
```

### 3.5 Metadata and events

With `metadata=` (a `MetadataCollection`, for example from
`shakelab.signals.io.metadata_reader`) the records receive the
metadata of their channel as they arrive, so that the buffer can be
corrected for the instrument response. `on_event` reports connections,
disconnections, gaps and overlaps between records, and the end of a
dial-up transfer (`SeedLinkEvent`). `client.info()`, `stations()` and
`streams()` describe what the server offers.

## 4. Command line

```sh
slab-seedlink info rtserve.example.org --streams IV.MODE
slab-seedlink monitor rtserve.example.org "IV.*..HHZ" --status 10
slab-seedlink get rtserve.example.org "IV.MODE..HH?" \
    --start 2026-11-05T12:00:00 --end 2026-11-05T12:05:00 -o mode.mseed
slab-seedlink archive rtserve.example.org "IV.*..HH?" \
    --sds archive --state iv.state
```

`info` lists the stations or the streams of a server, `monitor` prints
the records received or a table of the latencies, `get` writes a time
window to a MiniSEED file and `archive` archives continuously. The
server is `HOST[:PORT]` (port 18000 by default).
`python -m shakelab.signals.seedlink` is the same command.

## 5. Structure

| Module | Main objects | Role |
|---|---|---|
| `protocol` | `Connection`, `Packet`, `SeedLinkError` | connection, negotiation of the protocol version, commands and packets |
| `archiver` | `Archive` | MiniSEED archive files |
| `client` | `Client`, `read`, `archive`, `records`, `SeedLinkEvent`, `ChannelStatus` | requests and destinations of the records |
| `cli` | | the `slab-seedlink` command |

The records are decoded by `shakelab.signals.io.read_mseed_bytes`, with
the MiniSEED backend of `shakelab.signals.io`. A circular buffer holds
time-limited streams (`shakelab.signals.circular`, experimental).

## 6. Conventions

- Channels are selected with the patterns of `shakelab.signals.codes`;
  the selection is also applied by the client, so that only the
  channels requested reach the destinations even when the server
  supports coarser selections only.
- Times are ShakeLab `Date` objects in UTC; latencies are in seconds,
  from the last sample received to the time of the query.
- A record is archived and buffered whole, also when it overlaps the
  end of a requested time window.

## 7. See also

- The guide of the SeedLink client
  ([`docs/signals/seedlink.md`](../../../docs/signals/seedlink.md)).
- [`examples/signals/`](../../../examples/README.md#signals).
- The [`signals` README](../README.md), and the real-time monitors of
  array processing in
  [`shakelab.exploration.yarray`](../../exploration/yarray/README.md).

## 8. References

FDSN. SeedLink protocol, version 4. International Federation of
Digital Seismograph Networks, https://docs.fdsn.org/projects/seedlink.

SeisComP. SeedLink, in the documentation of SeisComP (GFZ Potsdam and
gempa GmbH): version 3.1 of the protocol, the de facto standard of
SeisComP and RingServer.
