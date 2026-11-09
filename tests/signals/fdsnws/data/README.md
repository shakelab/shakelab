# Test data of `signals.fdsnws`

Replies served by the fake FDSN server of the tests
(`tests/signals/fdsnws/fdsn_server.py`); the MiniSEED and StationXML
replies are `tests/signals/data/MODE.mseed` and `MODE.xml`.

| File | Content |
|---|---|
| `channels.txt` | Station service, text format, channel level: the three channels of `MODE.xml` and a fictitious `IV.MODE.00.HHZ` |
| `stations.txt` | Station service, text format, station level: `IV.MODE` and a fictitious `IV.MOXX` |
| `events.xml` | Event service, QuakeML 1.2: two events with the times, epicentres and magnitudes of the 20 and 29 May 2012 Emilia earthquakes (rounded values, written for these tests) |
| `availability.json` | Availability service, `query` method, JSON format (fictitious spans) |
| `extent.json` | Availability service, `extent` method, JSON format (fictitious) |

All files were written by hand for these tests (ShakeLab, LGPL).
