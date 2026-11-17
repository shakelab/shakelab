# libio: waveform file formats

`shakelab.signals.libio` holds one module for each waveform file
format that ShakeLab reads or writes. The usual way to use them is
not to call them directly but through `shakelab.signals.io`, whose
`reader` and `writer` recognize the format of each file, select
channels and time windows, and handle folders, patterns, archives and
file name templates (see the [`signals` README](../README.md)).

The modules are useful on their own when a file must be examined or
changed in the terms of its format: all the header fields of a SAC
file, the spectra of an ESM record, the raw records of a MiniSEED file,
the groups and properties of a TDMS file. Every module follows the
same scheme (Section 4), so that what is learned for one format holds
for the others, and a new format can be added in the same way
(Section 5).

## 1. Features

- MiniSEED 2 and 3 through libmseed (`cymseed`, the default), and
  MiniSEED 2 in pure Python (`mseed`), with the same results; raw
  access to the records of a file.
- SAC binary (`sac`), header versions 6 and 7, both byte orders; all
  the header fields kept through reading and writing.
- DYNA 1.2 (`dyna`), the ASCII format of ESM and ITACA: waveforms and,
  through the format class, response and Fourier spectra.
- Simple ASCII (`simpleascii`), the format of `mseed2ascii` and of the
  EarthScope data services, as lists or time-sample pairs.
- Text tables of any layout (`ascii`): delimiters, decimal commas,
  column names and `key: value` parameters recognized automatically.
- NI TDMS (`tdms`), the format of LabVIEW and of DAS interrogators,
  read with memory mapping and NI scales applied.

## 2. Quick start

A SAC file read through `io` becomes a `Record` with the coordinates
and orientation of the station in its channel metadata; the other
header fields are kept as they are:

```python
from shakelab.signals.io import reader

record = reader("event.sac")["XX.STA..HHZ"][0]
print(record.location, record.orientation)
print(record.head.extra["sac"])          # the other header fields
```

The format class gives the whole content of the file in the terms of
the format, and can be changed and written back:

```python
from shakelab.signals.libio.sac import Sac

sac = Sac.read("event.sac")
sac.info()                               # every defined header field
sac.head["KUSER0"] = "checked"
sac.write("event_checked.sac", owrite=True)

record = sac.to_shakelab()               # the same Record as above
```

## 3. Common tasks

### 3.1 Read and write a format directly

Each module has a reader and a writer named after the format, with
the same arguments for every format:

```python
from shakelab.signals.libio.dyna import dynaread, dynawrite
from shakelab.signals.libio.mseed import msread, mswrite

record = dynaread("XX.STA..HNE.ACC.ASC")
dynawrite("copy.asc", record, owrite=True)

collection = msread("data.mseed")
mswrite("steim1.mseed", collection, encoding=10, owrite=True)
```

Formats holding one waveform per file (SAC, DYNA) read and write a
`Record`; the others a `StreamCollection`. The same options are
accepted by `io.reader` and `io.writer` (`writer(data, "x.mseed",
encoding=10)`), which check them against the signature of the format
function.

### 3.2 Choose the MiniSEED backend

MiniSEED is read and written by libmseed through the `cymseed3`
package, or by the pure-Python module. The two give the same records;
libmseed is faster, reads and writes MiniSEED 3 and supports every
encoding, the Python module needs no compiled code:

```python
from shakelab.signals import io

collection = io.reader("data.mseed", mseed_backend="python")
io.USE_LIBMSEED = False                  # default for every call
```

`msrawread` and `msrawiter` of the Python module return the MiniSEED
records as they are in the file (`MSRecord`, with their sequence
numbers, encodings and record lengths), also those that are not
waveforms.

### 3.3 Spectra of ESM and ITACA

DYNA files also hold acceleration response spectra, displacement
response spectra and Fourier amplitude spectra. They are not time
series, so they are read by the format class only:

```python
from shakelab.signals.libio.dyna import Dyna

spectrum = Dyna.read("XX.STA..HNE.SA.ASC")
spectrum.info()
```

### 3.4 Text tables

`asciiread` reads tables written by spreadsheets or simple programs,
one record per data column. The sampling rate, start time, codes and
units can be written in the file as `key: value` lines or given as
arguments, which take precedence; a column of times gives the
sampling interval and the start time:

```python
from shakelab.signals.io import reader

data = reader("logger.csv", rate=200.0, start="2026-03-01T10:00:00",
              network="XX", station="LOG", units="m/s**2")
```

### 3.5 TDMS files

`tdmsread` reads the numeric channels of one group of a TDMS file as
records, taking sampling rate, start time and units from the
properties of the channels, of the group or of the file. `TdmsFile`
gives all the groups, channels and properties:

```python
from shakelab.signals.libio.tdms import TdmsFile

tdms = TdmsFile.read("acquisition.tdms")
tdms.info()
```

## 4. The common scheme

Every format module has the same parts, named after the format (`xxx`
below stands for `sac`, `dyna`, `ms`...):

| Part | Content |
|---|---|
| `xxxread(source[, stream_collection], **options)` | read a file as a `Record` (one waveform per file) or a `StreamCollection`, optionally added to a given collection |
| `xxxwrite(target, data, ..., owrite=False)` | write a `Record` or a collection; existing files are not overwritten unless `owrite=True` |
| `is_xxx(head)` | whether the first bytes of a file belong to the format, used by `io` to recognize files without a known extension |
| format class (`Sac`, `Dyna`, `SimpleAscii`, `MSRecord`; `TdmsFile` for reading) | the content of a file, or of one MiniSEED record, in the terms of the format |
| `XxxWarning` | warning class of the format |

The format class has:

- the header as a dictionary with the names of the format, where
  undefined fields are None (`sac.head["STLA"]`), and the data;
- `read(source)` and `write(target, owrite=False)`, and the conversion
  from and to bytes or text (`from_bytes`, `to_bytes`, `from_text`,
  `to_text`);
- `to_shakelab()` and `from_shakelab(record)`, the conversion to and
  from ShakeLab objects (for the formats that can be written);
- `info()`, and the properties `sid`, `start_ns` and `is_timeseries`
  (False for content that is not a waveform, such as spectra or text);
  `to_shakelab()` raises `NotWaveformError` (from `common.py`, a
  `ValueError`) for such content, so that `io.reader` can skip the
  file in folders and patterns.

In the conversion, the source identifier, start time, sampling and
units go to the record header; coordinates, burial depth and
orientation go to the channel metadata (`record.head.meta`); every
other field is kept in `record.head.extra["<format>"]` and written back
by `from_shakelab`, so that a file read and written again keeps its
information. Values are kept as written (units are declared, not
converted) unless the format requires a conversion, which is then
documented in the module.

A source can be a path, the content of a file (bytes or text), or a
file already opened. Modules are named after the format, not after the
organization that defined it, and are imported only when a file of
that format is read or written.

## 5. Adding a format

A new format follows the scheme of Section 4, taking `sac.py` as the
model:

1. Write `libio/<format>.py` with the header of ShakeLab and a module
   docstring describing the format (with its reference document) and
   the sections "Reading" and "Writing": what is read, what is kept
   where, what is converted and what is not supported. Define
   `__all__`, the reader, the writer if the format can be written,
   `is_<format>`, the format class and the warning class.
2. Register the format in `shakelab/signals/io.py`, in the table
   `_FORMATS`: description, module, names of the reader and writer
   functions, extensions, suffix of the files written, `per_record`
   (one waveform per file) and the name of the detection function;
   aliases of the name go to `_ALIASES`. The options of the reader and
   writer become options of `io.reader` and `io.writer` automatically.
3. Add the tests in `tests/signals/libio/test_<format>.py`, with small
   test files in `tests/signals/libio/data/<format>/` and their
   provenance and license in the README of the data folder. Compare
   with a reference implementation or with files written by the
   reference software of the format, and test reading and writing
   again (round trip), the detection and the errors.
4. Document the format: the table of formats in the
   [`signals` README](../README.md#5-formats), the guide of the file
   formats, an example if useful, and the CHANGELOG.

## 6. Structure

| Module | Format | Main objects |
|---|---|---|
| `cymseed` | MiniSEED 2 and 3 (libmseed) | `msread`, `mswrite` |
| `mseed` | MiniSEED 2 (pure Python) | `msread`, `mswrite`, `msrawread`, `msrawiter`, `msrawwrite`, `MSRecord`, `is_mseed` |
| `sac` | SAC binary | `sacread`, `sacwrite`, `Sac`, `is_sac` |
| `dyna` | DYNA 1.2 | `dynaread`, `dynawrite`, `Dyna`, `is_dyna` |
| `simpleascii` | Simple ASCII | `simpleasciiread`, `simpleasciiwrite`, `SimpleAscii`, `is_simpleascii` |
| `ascii` | text tables | `asciiread`, `asciiwrite` |
| `tdms` | NI TDMS | `tdmsread`, `TdmsFile`, `is_tdms` |
| `common` | definitions shared by the formats | `NotWaveformError` |

The pickle format of `io` and the metadata formats (StationXML, RESP,
the text format of the FDSN station service) are in
`shakelab.signals`, outside this folder.

## 7. See also

- The [`signals` README](../README.md): `reader`, `writer` and the
  table of formats.
- The guide of the file formats
  ([`docs/signals/io.md`](../../../docs/signals/io.md)).
- [`examples/signals/`](../../../examples/README.md#signals).

## 8. References

Goldstein, P., and Snoke, A. (2005). SAC availability for the IRIS
community. Incorporated Research Institutions for Seismology Data
Management Center Electronic Newsletter, 7(1).

IRIS (2012). SEED Reference Manual, version 2.4, chapter 8 (MiniSEED).
Incorporated Research Institutions for Seismology.

FDSN (2023). MiniSEED 3 specification. International Federation of
Digital Seismograph Networks.
