# Changelog

All notable changes to ShakeLab are recorded in this file: new
features, changes of behavior, deprecations, removals and scientific
corrections. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **License**: ShakeLab is now released under the GNU Lesser General
  Public License v3.0 or later (previously AGPL v3 in the source
  headers and GPL v3 in the license file). The license texts are in
  `COPYING.LESSER` and `COPYING`.
- Python 3.10 or later is required.
- `pyproject.toml` is the single source of package metadata, version
  and dependencies; `setup.py`, `requirements.txt` and `version.txt`
  were removed. The installed version is available as
  `shakelab.__version__`.
- `libutils.timeN.Date` stores an absolute UTC instant as integer
  nanoseconds from the internal epoch `0001-01-01T00:00:00Z`; it
  supports years before 1 (proleptic Gregorian calendar), NumPy
  scalars, stricter validation, and numeric offsets on ISO 8601 output.
- `libutils.geodeticN`: vectorized distance and transform functions,
  validated `WgsPoint`, longitudes accepted in `[-180, 180)` or
  `[0, 360)`, antimeridian-safe polygons and meshes, explicit `rng`
  for random sampling, CSV and 3-D GeoJSON I/O.
- `libutils.ascii.AsciiTable` and `libutils.json` were rewritten with
  the same API; errors are reported with exceptions instead of
  printed messages.
- `libutils.utils.a_round` returns a rounded copy instead of modifying
  its input.
- `libutils.timezones.TIMEZONE`: outdated offsets updated (e.g.
  Caracas, Amman, Damascus, Casablanca, Samoa, Minsk, Novosibirsk,
  Astana, Greenland, Asuncion) and an IANA name added to each entry.

### Added

- `docs/conventions.md` (library-wide conventions) and
  `docs/development.md` (development and architecture guidelines,
  version 2.1), replacing the former Word document; the normative
  conventions, previously an appendix of the guidelines, are now in
  `docs/conventions.md`.
- `libutils.timeN.Date`: Unix time conversions (`from_unix`,
  `to_unix` and nanosecond variants), `replace`, `nanosecond`,
  `is_set`.
- `libutils.timezones`: `utc_offset` and `offset_string` (daylight
  saving time through the standard-library `zoneinfo`),
  `format_offset`, `find_timezone`.
- `libutils.ascii.AsciiTable`: quoted fields, `add_row`, `to_array`,
  `unique`, `filter`, `select`, `sort`.
- `libutils.json`: `write_str`, `JsonParseError` with line and column,
  `NaN`/`Infinity` support, `allow_comments`, `float_format`,
  `sort_keys`, `inline_numbers`, `allow_nan`.
- `libutils.utils.pick`.
- Test suite for all `libutils` modules; continuous integration with
  GitHub Actions on Linux and macOS, Python 3.10 to 3.14.

### Fixed

- `signals.base.Header.location`: a `(longitude, latitude)` sequence
  was stored with longitude and latitude swapped.
- `seismicity.catalogue` hypocentre: longitude and latitude were
  swapped, and the depth in km was used as an elevation in metres; the
  elevation is now `-1000 * depth_km`.
- Run-time data files (GMPE coefficients, web interface assets) were
  missing from non-editable installations.
- `libutils.utils.cast_value` and `is_empty` failed on NumPy arrays;
  `serialize_ndarray` did not convert NumPy integer scalars.
- `libutils.ascii.AsciiTable.write`, `rename_key` and `concatenate`
  always failed; reading a file with unnamed columns misaligned the
  values.
- `libutils.json` could write invalid JSON (unescaped backslashes,
  newlines and control characters, `nan`, tuples) and did not decode
  all JSON escapes.
- `libutils.timeN`: ISO 8601 formatting with zero precision, rounding
  to 60 seconds, `now()` resolution and hashing of equal dates.

### Deprecated

- `libutils.utils.slice`: use `pick` (the old name shadows the
  built-in `slice`).

### Removed

- `libutils.time`, replaced by `libutils.timeN`.
- `libutils.geodetic`, replaced by `libutils.geodeticN` (note the
  `(longitude, latitude)` argument order of `WgsPoint`).
- `libutils.timeN` aliases `MSEC`, `HSEC`, `DSEC`, `YDAYS`.
- Unused dependency `shapely`.
