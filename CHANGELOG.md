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
  complete response (all stages, digital filters included, each at its
  own sampling rate) instead of stage by stage; a subset is chosen with
  `stages=...`, e.g. `stages=("gain", "paz")` for the analog part. The
  units of the data become the input (deconvolution) or output
  (convolution) units of the response. `fourier.inverse_spectrum`
  follows the same relative definition, with the new default method
  `"water"` (amplitude clipped to the water level, phase kept).
- `signals.response` no longer depends on waveform classes: stages
  provide `transfer_function(freq, delta=None)` and stage sets the
  product over the selected stages (`ACTIVE_STAGES` replaced by the
  `stages` argument). Digital stages (FIR, digital coefficients,
  digital poles and zeros) are evaluated at their own input sampling
  rate and include the time correction declared by the datalogger
  (conventions of SEED v2.4, Appendix C). `StageFIR` renamed
  `StageCoefficients`, with `cf_type` (`DIGITAL`,
  `ANALOG (RADIANS/SECOND)`, `ANALOG (HERTZ)`) as `pz_type` in
  `StagePoleZero`. `StagePolynomial` now represents a non-linear sensor
  (`x = sum(a_k y**k)`: the input as a function of the output, as in
  SEED blockette 62 and StationXML `Polynomial`) with `evaluate`
  applied to output values; `transfer_function` gives the gain
  `1 / a_1` of a linear polynomial. It was evaluated as a rational
  function of the frequency.
  A missing channel or epoch raises `KeyError` (it printed a message
  and returned None).
- **Channel metadata**: `signals.metadata` was rewritten (the previous
  module, based on nested dictionaries, could not be imported). Channel
  epochs (`ChannelEpoch`: coordinates, burial depth in km,
  `Orientation`, sampling rate, sensor and datalogger `Equipment`,
  instrument response, station) and station epochs (`StationEpoch`)
  are held by `MetadataCollection`, read from StationXML, RESP
  (response only) or JSON with `MetadataCollection(file)` or `read`;
  `get(sid, time)` returns the epoch valid at a time, `select(sid,
  time)` a sub-collection (shell wildcards), `stations()` the station
  epochs, `validate()` the inconsistencies, `write()` a JSON file.
  `MetadataCollection` replaces `signals.response.ResponseCollection`
  and `StreamResponse`, which were removed: the response model no
  longer holds identifiers and epochs, so `StageSet(sample_rate=None)`
  has no `starttime`, `endtime` and `match`, and gains
  `to_dict`/`from_dict`. `select_response` and `select_sensitivity`
  moved to `signals.metadata`; response corrections accept stages,
  stage sets, channel epochs and metadata collections. The JSON format
  changed accordingly (stations and channel epochs). Without a time,
  `get` requires a channel with a single epoch (the stage set valid
  "now" was used before). `StageSet.append` raises `TypeError` for
  objects that are not stages (it raised `ValueError`), and indexing
  a stage set with a non-integer raises `TypeError` (it raised
  `NotImplementedError`).
- **Metadata of the records**: the source identifier and the
  description of the channel are held by `Header.meta`, a
  `ChannelEpoch` always present (empty by default), shared by the
  records of the same epoch. `Record.location`, `orientation` and
  `response` read it; assigning them (or `sid`) replaces the epoch of
  that record with an updated copy (`ChannelEpoch.replace`), so copies
  are never affected. `set_metadata(source)` on `Record`, `Stream` and
  `StreamCollection` assigns the epoch valid at the start of each
  record from a `MetadataCollection` (or a given `ChannelEpoch`), one
  copy per epoch, with a `MetadataWarning` when different metadata are
  replaced, when the declared sampling rate differs, or when a record
  extends beyond its epoch; streams and collections leave channels
  missing from the collection unchanged, with a single warning listing
  them. The response
  corrections (`deconvolve_response`, `convolve_response`,
  `remove_sensitivity`) of `Record`, `Stream`, `StreamCollection` and
  `Spectrum` use the response of the metadata when called without a
  response, and raise `ValueError` when the declared units of the data
  differ from those expected by the response (e.g. data already
  corrected). Records with different channel metadata are not joined:
  `Record.append` raises `ValueError`, `Stream.append` keeps them as
  separate records; empty metadata are compatible with any. The source
  identifier must be an FDSN code `NET.STA.LOC.CHA` (or None).
  The free dictionary of a record is `Header.extra` (it was
  `Header.meta`); free information about the channel goes in
  `ChannelEpoch.extra`. `ChannelEpoch` accepts `sid=None` (empty epoch
  of a record without identifier) and provides `empty`.
- `signals.response.StageSet.transfer_function`: poles-and-zeros and
  coefficient stages are normalized to unit amplitude at the frequency
  of the gain of the same stage, so that the stage amplitude there is
  the declared gain (SEED blockette 58, StationXML `StageGain`), as in
  evalresp. Nothing changes when the gain and normalization frequencies
  coincide and the metadata are consistent; otherwise the response
  previously differed from the declared gain (e.g. by 1.5% for
  IU.ANMO.10.BHZ, normalized at 0.1 Hz with the gain at 0.02 Hz).
- `signals.stationxml.parse_sxml` renamed `read_stationxml`; it accepts
  paths, XML text with or without declaration, bytes and open files,
  and returns a `MetadataCollection` with stations, channel epochs
  (coordinates, depth, orientation, equipment) and responses; channels
  without response are kept. `parse_response(element, sample_rate)`
  no longer takes the epoch.
  `stationxml_to_dict` returns `{"tag", "attrib", "text", "children"}`
  and keeps values as strings (codes such as `"00"` became floats).
- `signals.respfile.parse_resp` renamed `read_resp`; it accepts paths,
  RESP text, bytes and open files, and returns a `MetadataCollection`
  with one channel epoch per RESP epoch. Lines are identified by their
  `BxxxFyy` code, so RESP files of rdseed, of the Nominal Response
  Library and of evalresp are read alike. The sampling rate of an
  epoch is the output rate of the last decimation stage, as in
  evalresp.
- RESP files are recognized by `MetadataCollection.read`
  (`ftype="resp"`, inferred for file names starting with `RESP.` or
  ending with `.resp`).
- `signals.Header` moved to the new module `signals.header` (it
  depends only on time primitives and on `signals.metadata`); it is
  still available from `signals.base`.
- `signals.base.Header.time` accepts every format accepted by `Date`
  (strings, seconds, sequences, dictionaries) and stores its own copy,
  so later changes to the assigned `Date` no longer move the record;
  `Record.location` also accepts NumPy arrays.
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
- `signals.libio.mseed` (pure-Python MiniSEED 2 backend) rewritten:
  byte order detected from each record (the `byte_order` argument
  defaults to None), float64 data (encoding 5), little-endian STEIM
  frames, record lengths from 128 to 65536 bytes, blockette 100
  (actual sampling rate), and sources given as paths, bytes or binary
  file objects. STEIM decoding is vectorized and the records of each
  channel are sorted and joined into one `Record` per contiguous
  segment before building the waveforms: reading is about 5 times
  faster for STEIM and 20-30 times for uncompressed data, and no longer
  grows quadratically with the file length. Records without waveform
  (log, detection, timing) are skipped; records with unsupported
  encodings are skipped with one `MSeedWarning`; a mismatch of the
  last STEIM sample gives a warning instead of an exception; trailing
  bytes after the last record are ignored (with a warning unless they
  are zeros). Integer samples are returned in native byte order.
- MiniSEED writing (both backends): with `encoding=None` (new default)
  the encoding follows the data type: STEIM2 for integers (STEIM1, or
  int32, when the differences between samples exceed its words), 4 for
  float32, 5 for float64. Start times keep the microseconds in
  blockette 1001 when needed.
- `signals.libio.mseed.MSRecord`: built with `MSRecord.from_bytes` and
  encoded with `to_bytes`; `code` renamed `sid`; new `rate`,
  `start_ns`, `encoding` and `is_waveform`; `delta` is None for records
  without sampling.
- `signals.libio.cymseed` no longer depends on `signals.fdsnws`;
  integer samples are int32. Until `cymseed3` packs non-integer
  samples correctly, MiniSEED 2 files needing encodings 1, 4 or 5 are
  written by the pure-Python backend, and MiniSEED 3 output is limited
  to encodings 3, 10 and 11 (`NotImplementedError` otherwise).
- `libutils.timezones.TIMEZONE`: outdated offsets updated (e.g.
  Caracas, Amman, Damascus, Casablanca, Samoa, Minsk, Novosibirsk,
  Astana, Greenland, Asuncion) and an IANA name added to each entry.

### Added

- `signals.io.read_metadata(source, ftype=None)`: reads channel
  metadata (StationXML, RESP, JSON) from a file, a wildcard pattern, a
  directory or a sequence of them into one `MetadataCollection`.
  `io.reader(..., metadata=...)` and `StreamCollection.read(...,
  metadata=...)` assign the metadata (a collection, a channel epoch, or
  files read with `read_metadata`) to the records read; records already
  in a given collection are not modified.
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
- `signals.response`: `Sensitivity` (overall sensitivity, read from
  StationXML `InstrumentSensitivity`), `select_sensitivity` and
  `remove_sensitivity` on `Record`, `Stream`, `StreamCollection` and
  `Spectrum`; digital poles and zeros (`DIGITAL (Z-TRANSFORM)`) and
  analog coefficients;
  `StageResponseList` (amplitude and phase tables); `StageUnsupported`
  (elements of unknown type raise an error when used);
  `StageSet.validate()` (gains, sensitivity, normalization, units,
  decimation chain, FIR symmetry); `ResponseWarning`; `same_units`.
- `signals.respfile`: transfer function types (blockettes 53 and 54),
  denominators of blockette 54, response lists (55), decimation and
  time correction (57), overall sensitivity (58 of stage 0), FIR
  filters with symmetry (61), polynomials (62, stage 0 as
  `instrument_polynomial`), generic responses and response references
  (56, 60) as `StageUnsupported`; continued blockettes are merged;
  `parse_resp_time` for SEED times, also truncated.
- `signals.response.expand_fir_symmetry`, shared by the StationXML and
  RESP readers.
- `signals.stationxml`: reading of `Coefficients`, `ResponseList`,
  `InstrumentSensitivity`, `InstrumentPolynomial` and `Decimation`;
  `iter_channels`.
- `examples/`: runnable examples organized by subpackage, run by the
  test suite; first example: instrument response correction of the
  IV.MODE records compared with the ESM accelerations.
- `libutils.xml`: reading of XML documents from any source with
  namespaces removed (`read_xml`), conversion to dictionaries
  (`to_dict`), `find_text`, `find_float`.
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

- `signals.libio.mseed`: a negative sampling-rate multiplier gave a
  negative rate (it multiplied by 1/m instead of dividing by |m|); the
  time correction of the fixed header was ignored when not yet applied;
  the microseconds of blockette 1001 were read without sign and
  ignored; files written by libmseed in little-endian order, with
  float64 data, text or detection records could not be read; writing
  floating-point data failed with the default encoding.
- `signals.libio.cymseed`: input files were not closed; codes padded
  with spaces by libmseed (e.g. channel `Z`) gave source identifiers
  with trailing spaces.
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
- `signals.response`: reading a response file failed with the default
  file type; reading JSON with a `coefficients` stage raised
  `NameError`, and several epochs of one channel were dropped; the
  poles-and-zeros response modified the stage when zeros were missing.
- `signals.libio.dyna`: the header was assumed to have 64 lines; with
  the 65-line header of ESM 3.0 files the first sample was lost and the
  series was shifted by one sampling interval.
- `signals.respfile`: the zeros of poles-and-zeros stages were never
  read (e.g. the response of a velocity sensor lacked the factor
  `s**2`); all channels and epochs of a file were merged into one
  epoch; truncated SEED times (e.g. `1998,299,20`) and `No Ending Time`
  were replaced by the current time; the overall sensitivity was read
  as a stage; the frequency of the gains was set to 1 Hz; transfer
  function types and units were ignored.
- `signals.response.StagePoleZero`: a stage without poles raised
  `ValueError`; SEED and StationXML allow it (e.g. a pure gain).
- `signals.stationxml`: `Coefficients` stages were skipped, so the
  output units of the response were lost; `FIR` elements were read as
  a unit filter (`NumeratorCoefficient` ignored) and symmetric filters
  were not expanded; digital poles and zeros were evaluated as analog;
  `ResponseList` stages were dropped without notice; a channel without
  response raised `AttributeError`; joining lines could merge words of
  multi-line texts.
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
- `signals.Header`: `location` (constructor argument and attribute)
  and `response`; they are part of the channel metadata
  (`Record.location`, `Record.response`, `Header.meta`).
- `signals.response`: `convolve_record`/`deconvolve_record` of stages
  and stage sets, `to_spectrum` and `response_function` of stages
  (use `transfer_function` and the `*_response` methods of the data
  classes), and `response.inverse_spectrum` (see
  `fourier.inverse_spectrum`).
- Unused dependency `shapely`.
- `signals.libio.mseed`: `MSRecord.append`, `MSRecord.read`/`write`
  on `ByteStream` objects (replaced by `from_bytes`/`to_bytes`) and
  the `ByteStream` arguments of the module functions; the module no
  longer uses `signals.binutils`.
