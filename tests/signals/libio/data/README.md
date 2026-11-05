# Test data for `shakelab.signals.libio`

| File | Content |
|---|---|
| `mseed/libmseed/` | Unmodified copy of the folder `test/data` of libmseed (EarthScope, [github.com/EarthScope/libmseed](https://github.com/EarthScope/libmseed), commit `d6e6ad3`, 5 August 2026), licensed under the Apache License 2.0 (`mseed/libmseed/LICENSE`): MiniSEED 2 and 3 files with the same series in every encoding and byte order (`reference-testdata-*`), text and detection records, historical encodings (CDSN, DWWSSN, GEOSCOPE, SRO), invalid blockette offsets, records without blockette 1000, records of mixed lengths in mixed order, an unapplied time correction, odd sampling rates, dates before 1970 and nanosecond times |
| `mseed/libmseed_reference.json` | What libmseed (same commit, `ms3_readtracelist`) reads from the MiniSEED 2 files of `mseed/libmseed` and of `tests/signals/data`: source identifier, start time, rate, number of samples and SHA-256 of the samples as big-endian float64, for each segment. Files that ShakeLab does not support (historical encodings, invalid blockettes, no blockette 1000) and segments without samples or with text are left out |
| `sac/obspy/` | SAC files from the test data of ObsPy ([github.com/obspy/obspy](https://github.com/obspy/obspy), `obspy/io/sac/tests/data`, commit `50033e1`, 28 August 2026), licensed under the GNU LGPL v3 (`sac/obspy/LICENSE.txt`), unmodified: `seism.sac` (the example distributed with SAC), `test.sac` and `test.sac.swap` (the same file in both byte orders), `LMOW.BHE.SAC`, `dis.G.SCZ.__.BHE_short`, `null_terminated.sac`, `non_ascii.sac`, `test_encode.sac` (text fields holding integers) and `testxy.sac` (alphanumeric SAC) |
| `simpleascii/` | Simple ASCII files written by `mseed2ascii` 2.6 (EarthScope, [github.com/EarthScope/mseed2ascii](https://github.com/EarthScope/mseed2ascii), commit `6e54686`) from the MiniSEED 2 files of `mseed/libmseed`: `int32` (sample list in one and six columns, time-sample pairs), `float32`, `float64` (sample list, time-sample pairs), `oddrate` (1080 sps, start with microseconds) and `3channel` (three channels of IU.COLA), with `mseed2ascii -f 1 [-c 6]` or `-f 2` and `-o <file>`. The program is GPL-3.0; the files are data generated from the libmseed test files |
| `tdms/nptdms/` | Example files of npTDMS ([github.com/adamreeve/npTDMS](https://github.com/adamreeve/npTDMS), `nptdms/test/data`, commit `ab68cac`, 4 November 2026), licensed under the GNU LGPL v3 (`tdms/nptdms/LICENSE.txt`), unmodified: `big_endian.tdms`, `Digital_Input.tdms` (three groups), `raw_timestamps.tdms` (LabVIEW waveform) and `raw1.tdms` (DAQmx raw data, not supported) |
| `tdms/scenarios/` | One file per test scenario of npTDMS (same commit, `nptdms/test/scenarios.py`: object lists reused and extended across segments, chunked, interleaved and incomplete data, booleans, floats with units, complex numbers, linear scales, timestamps), with the expected channel data in `expected.npz` (key `<scenario>/<group>/<channel>`) |
| `tdms/generated/` | Files written with the npTDMS writer: text and numeric channels in two segments (`mixed`) and the same file truncated by 10 bytes (`truncated`), LabVIEW waveforms with start offset and units (`waveform`), and NI scales (`scaled`: polynomial, linear followed by polynomial, subtraction, thermocouple, already scaled, scale of the group) |
| `tdms/das/` | Excerpt of an iDAS acquisition (Silixa iDAS v3, 200 Hz, 2 m spatial resolution, 10 m gauge length) of 4 May 2026: channels 0-15 of the first 2 s, written as two consecutive files of 1 s with the properties of the original file (start times of the second file advanced by 1 s); the identifying fields (file and host names, serial numbers, GPS coordinates of the interrogator) are replaced by fictitious values |

The MiniSEED reference is generated with a small libmseed program and a Python
script kept outside the repository
(`AiSpace/Claude/Shakelab/tools/libmseed_reference/`); the script
docstring gives the commands. The TDMS files of `nptdms`, `scenarios`
and `generated` are written by `tdms_testdata.py`, kept in the same
folder.

Real waveforms used by these tests (station IV.MODE and the ESM
records of the 2012 Emilia earthquake) are in `tests/signals/data`,
with their sources and licenses.
