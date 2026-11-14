# ShakeLab conventions

*Version 1.0, 14 November 2026*

This document lists the scientific and technical conventions adopted
across the ShakeLab library: how time, coordinates, units and numbers
are represented. They apply to every module and to the values exchanged
through the public API.

The conventions are **normative**. A module that intentionally departs
from one of them declares it in the "Conventions" section of its module
docstring. Code that does not yet comply is considered technical debt
and is being corrected; see the [development guidelines](development.md)
for the principles behind these choices.

## Contents

1. [Time](#1-time)
2. [Space and coordinates](#2-space-and-coordinates)
3. [Units](#3-units)
4. [Numerical representation](#4-numerical-representation)
5. [Random numbers and reproducibility](#5-random-numbers-and-reproducibility)
6. [Supported versions](#6-supported-versions)

## 1. Time

Absolute times are handled by `Date` (`shakelab.libutils.timeN`).

| Topic | Convention |
|---|---|
| Absolute time | `Date`: an absolute UTC instant stored as one integer number of nanoseconds from the internal epoch `0001-01-01T00:00:00Z`. |
| Calendar | Proleptic Gregorian calendar, astronomical year numbering (year 0 = 1 BC), ISO 8601 text representation. |
| Time scale | UTC. ISO 8601 offsets are accepted on input and converted to UTC; the offset is not retained. |
| Time zones | Local times are produced only for output, by passing a numeric offset to `Date.to_iso8601(timezone="+02:00")`. Offsets of named zones, daylight saving time included, are obtained with `shakelab.libutils.timezones.offset_string` (codes of the `TIMEZONE` table or IANA names such as `"Europe/Rome"`). |
| Leap seconds | Not represented: seconds lie in [0, 60). Differences between dates spanning a leap second differ from the true elapsed SI time by one second per leap second, as in POSIX time. Applications needing exact elapsed time across leap seconds must handle it explicitly. |
| Resolution | 1 ns; finer inputs are rounded to the nearest nanosecond, ties to even. |
| Unix time | Exchange format only, through `Date.from_unix`, `Date.to_unix` and their nanosecond variants. |
| Bare numbers | A bare number used as a date means seconds from the internal epoch. |
| Durations and shifts | Real numbers in seconds (int, float, `Decimal` or NumPy scalars). Non-finite values and booleans are rejected. |
| Sampling | Sampling interval `delta` in s and sampling rate `rate` in Hz, kept consistent (`rate = 1 / delta`). Sample `i` of a record starts at `starttime + i * delta`, computed from the start time and the index, never by accumulating increments. |
| Timing error | Metadata of a date: preserved by copies and shifts, ignored by comparison and hashing. |
| Time tolerances | Tolerances on time alignment are fractions of the sampling interval, defined once in `shakelab.signals.base`: `MERGE_TOLERANCE` (0.5) decides whether two records are contiguous, absorbing the time resolution of data formats; `SNAP_TOLERANCE` (1e-6) decides whether a time coincides with a sample, absorbing floating-point rounding. |

Example:

```python
from shakelab.libutils.timeN import Date
from shakelab.libutils.timezones import offset_string

date = Date("2026-07-15T12:00:00Z")
date.to_iso8601(precision=0, timezone=offset_string("Europe/Rome", date))
# '2026-07-15T14:00:00+02:00'
```

## 2. Space and coordinates

Geographic objects and transformations are handled by
`shakelab.libutils.geodeticN`.

| Topic | Convention |
|---|---|
| Datum | WGS84 geodetic coordinates. |
| Order | Always `(longitude, latitude[, elevation])`: in constructors (`WgsPoint(longitude, latitude, elevation)`), tuples, function arguments and array columns (column 0 longitude, column 1 latitude). This matches GeoJSON and the x/y order of projected coordinates. Keyword arguments are recommended at call sites. |
| Angles | Degrees in the public API, radians only internally. |
| Longitude | Accepted in either `[-180, 180)` or `[0, 360)` and kept as given. Algorithms work on longitude differences, so results do not depend on the convention and remain valid across the antimeridian. Generated longitudes follow the convention of the input (e.g. of `lon_min` for a range). `normalize_longitude` converts to `[-180, 180)`. A range with `lon_max < lon_min` crosses the antimeridian. |
| Latitude | In `[-90, 90]`. |
| Elevation | Metres above the WGS84 ellipsoid (ellipsoidal height), positive upward. Heights referenced to mean sea level may be used as ellipsoidal heights only as an explicit, documented approximation made at the I/O boundary (the geoid undulation can reach about 100 m). |
| Depth | Seismological domain objects (hypocentres, catalogues, rupture contexts) use depth in km, positive downward. The only conversion to geometry is `elevation_m = -1000 * depth_km`, performed in one identified place per domain object. |
| Distances | Metres in the geodetic primitives (great-circle, epicentral, hypocentral, chord distances); areas in square metres. Ground-motion distance metrics (Repi, Rhypo, Rjb, Rrup, Rx) are in km, as in the published models (see [Units](#3-units)). |
| Azimuth | Degrees clockwise from geographic north, in `[0, 360)`. |
| Local frames | `MetricFrame` and `MetricPoint`, in metres: ECEF, or ENU with x = East, y = North, z = Up. A metric coordinate is always associated with an explicit frame. Map projections such as UTM are intentionally not provided. |
| Components | Ground-motion components follow SEED orientation codes: Z positive up, N positive north, E positive east; other codes (1, 2, ...) carry an explicit azimuth and dip. |

Example:

```python
from shakelab.libutils.geodeticN import WgsPoint

hypocentre = WgsPoint(longitude=13.1, latitude=46.3, elevation=-1000.0 * 8.0)
station = WgsPoint(longitude=13.25, latitude=46.06)
hypocentre.hypocentral_distance_to(station)  # metres
```

## 3. Units

Primitives and scientific computations use SI units: m, s, Hz, kg, Pa,
m/s, m/s². The following exceptions are adopted because they are
established in seismological and engineering practice. Any other
non-SI quantity is declared in the conventions of its module and
carries its unit in the name (for example `depth_km`, `distance_km`).

| Quantity | Unit |
|---|---|
| Hypocentral depth | km, positive downward |
| Ground-motion distance metrics | km (Repi, Rhypo, Rjb, Rrup, Rx) |
| PGA, spectral acceleration | g (standard gravity, 9.80665 m/s², `shakelab.libutils.constants.GRAVITY`) |
| PGV | m/s |
| PGD, spectral displacement | m |
| Vs30 and shear-wave velocities | m/s |
| Basin depths (Z1.0, Z2.5) | as defined by each model, stated in the name or docstring |
| Magnitude, intensity | dimensionless; the scale (Mw, ML, MMI, EMS-98, ...) is always explicit |
| Geographic angles | degrees |
| Signal amplitudes | stored in the `units` field of the record header; never assumed |

Ground-motion models return medians in the units above and standard
deviations in natural-log units. Conversions to other units (cm/s, %g,
gal) happen only in I/O and presentation.

## 4. Numerical representation

- Floating-point data are `float64` unless a format or a memory
  constraint justifies another dtype; integer samples read from a
  format keep their dtype until a computation requires conversion.
- Absolute times are never stored as floating-point seconds over long
  spans; `Date` (integer nanoseconds) is used instead.
- Missing numerical values are `numpy.nan` in float arrays and `None`
  for scalar metadata; sentinel values from file formats are converted
  at the I/O boundary.

## 5. Random numbers and reproducibility

Stochastic functions and configurations (sampling, simulations, Monte
Carlo, random meshes) take a parameter named `seed`:

- an integer makes the result reproducible;
- `None` (the default) uses a new generator with fresh entropy
  (`numpy.random.default_rng()`): the result is not reproducible;
- functions may also accept a `numpy.random.Generator`, to continue
  the stream of the caller.

The global random states of NumPy (`numpy.random.seed`,
`numpy.random.rand`) and of Python (`random`) are never used. In
configurations stored as files (JSON) the seed is an integer or
`null`. Library code never depends on wall-clock time, environment
variables or the iteration order of unordered containers in ways that
change scientific results.

## 6. Supported versions

| Item | Requirement |
|---|---|
| Python | >= 3.10 (tested on 3.10 to 3.14) |
| NumPy | >= 1.24, compatible with NumPy 2.x |
| SciPy | >= 1.10 |
| Optional extras | `plot` (matplotlib >= 3.7, for the figures), `webui` (Flask), `dev` (pytest); `all` installs every run-time extra, including the wxPython GUI dependency |

`pyproject.toml` is the single source of the version requirements and
of the package version (`shakelab.__version__`).
