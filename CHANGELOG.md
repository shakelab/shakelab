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
- `signals.Spectrum` moved from `signals.fourier` to the new module
  `signals.spectrum`; `signals.fourier` now holds only array functions.
  The constructor takes spectral values (`Spectrum(data, delta,
  nsamp)`); spectra of time series are built with
  `Spectrum.from_record`, `Spectrum.from_samples` or
  `Record.to_spectrum`. `nsamp` is stored explicitly. Removed the
  `norm` option and `set_data`/`fft`/`ifft` methods (replaced by the
  constructor, `from_samples` and `to_samples`); `resample` renamed
  `interpolate`; `filter(filt_type=...)` renamed `filter(method=...)`
  with `minphase`.
- `signals.fourier.butterworth(freq, corner, order, btype, minphase)`:
  explicit filter type instead of negative orders.
- Instrument response correction: the water level is now relative, in
  dB below the maximum of the response amplitude (default 60 dB), and
  set once for the whole response with
  `deconvolve_response(resp, waterlevel=..., method=...)` on `Record`,
  `Stream`, `StreamCollection` and `Spectrum`. The correction uses the
  total response of the active stages instead of stage by stage, and
  the units of the data become the input (deconvolution) or output
  (convolution) units of the response. `fourier.inverse_spectrum`
  follows the same relative definition, with the new default method
  `"water"` (amplitude clipped to the water level, phase kept).
- `signals.response` no longer depends on waveform classes: stages
  provide `transfer_function(freq, delta=None)` and stage sets the
  product over the active stages. Reading a missing channel or epoch
  raises `KeyError` (it printed a message and returned None).
  `ResponseCollection(ftype=None)` infers the format from the file
  extension; `ResponseCollection.append` accepts stream responses only
  and merges the epochs of an existing channel.
- `signals.Header` moved to the new module `signals.header` (it
  depends only on time and geodetic primitives); it is still available
  from `signals.base`.
- `signals.base.Header.time` accepts every format accepted by `Date`
  (strings, seconds, sequences, dictionaries) and stores its own copy,
  so later changes to the assigned `Date` no longer move the record;
  `Header.location` also accepts NumPy arrays.
- `signals.base`: contiguous records are joined when the start time is
  within half a sample of the expected one (previously 1e-6 samples),
  absorbing the time resolution of data formats.
- `signals.base`: the `tolerance` parameter is replaced by
  `merge_tolerance` (joining records: `append`, `merge`) and
  `snap_tolerance` (matching times to samples: `cut`, `extract`);
  `Stream.extract` and `StreamCollection.extract` accept both.
- `signals.base.Stream`: slicing returns a `Stream` (same source
  identifier, the stored records) instead of a list; numerical times
  in `extract` are seconds from the start of the stream (previously
  from the start of each record).
- `signals.base.Stream`: `filter`, `resample`, `taper`, `remove_mean`,
  `shift`, `integrate` and `differentiate` apply the operation to
  every record; the same methods are available on `StreamCollection`.
- `signals.base.StreamCollection`: slicing returns a collection;
  `extract` without `sid` extracts the window from every stream; `*`
  and `/` accept every response container.
- `signals.base`: a stream response is applied only to records with
  its own source identifier (previously, e.g., `collection / rc[sid]`
  silently applied one channel's response to all channels).
- `signals.base.Record.cut`/`extract` accept windows not aligned with
  the sampling grid and keep the samples inside the window.
- `signals.base.Record` operators `*` and `/` accept every response
  container (stage, stage set, stream response, response collection).
- `libutils.timezones.TIMEZONE`: outdated offsets updated (e.g.
  Caracas, Amman, Damascus, Casablanca, Samoa, Minsk, Novosibirsk,
  Astana, Greenland, Asuncion) and an IANA name added to each entry.

### Added

- `signals.base`: `MERGE_TOLERANCE` (0.5 samples, used when joining
  records) and `SNAP_TOLERANCE` (snapping times to samples);
  `Record(units=...)`;
  `Record.time_axis(reference=date)` for seconds relative to any date.
- `signals.fourier`: `fft`/`ifft` with the `delta * DFT` scaling,
  `smooth` (Gaussian in log-frequency and Konno-Ohmachi windows, at the
  input or at chosen frequencies), `smooth_amplitude` and
  `smooth_amplitude_phase`; `inverse_spectrum` moved here from
  `signals.response`. `Spectrum.logsmooth` gains `method`, `bandwidth`
  and `phase`.
- `signals.response.select_response` (response applicable to a channel
  and time), `StageSet.input_units`/`output_units`,
  `Spectrum.convolve_response`/`deconvolve_response`. Poles and zeros
  given in Hz (StationXML `LAPLACE (HERTZ)`) are supported.
- `signals.circular.CircularStream` (experimental): time-limited
  buffer of waveform blocks for real-time acquisition, used by
  `StreamCollection(max_duration=...)`. It copies appended records and
  refuses block-by-block processing.
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

- `signals.base.Record.resample`: Fourier resampling produced a
  sampling interval different from the declared one; linear
  interpolation could assign the last value to a time beyond the end
  of the data.
- `signals.fourier.butterworth` failed with NumPy 2 (`np.complex_`)
  and returned |H|^2 instead of |H| (-6 dB at the corner); for odd
  orders it was infinite at the corner and negative above it. The
  minimum-phase option now gives the causal analog response.
- `signals.Spectrum`: the constructor with an array failed; `resample`
  used a missing attribute; `filter` ignored unknown types; the header
  copied from a record held a full copy of the record; `logsmooth`
  averaged the wrapped phase. `site.soil.Model1D.soil_response`, which
  failed through the constructor, works again.
- `signals.base.Record.shift` failed on complex data; `filter`
  rejected NumPy integer orders.
- `signals.response`: `ResponseCollection(file)` failed with the
  default file type; reading JSON with a `coefficients` stage raised
  `NameError`, and several epochs of one channel were dropped; the
  poles-and-zeros response modified the stage when zeros were missing.
- `signals.base.Stream.append` stored the caller's record and extended
  it in place when later records were merged, so the caller's record
  silently changed; it also overwrote its source identifier. The
  stream now stores a copy.
- `signals.base.Stream`: NumPy integers are accepted as indices;
  `True`/`False` are no longer taken as positions 1 and 0.
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
- `signals.response`: `convolve_record`/`deconvolve_record` of stages
  and stage sets, `to_spectrum` and `response_function` of stages
  (use `transfer_function` and the `*_response` methods of the data
  classes), and `response.inverse_spectrum` (see
  `fourier.inverse_spectrum`).
- Unused dependency `shapely`.
