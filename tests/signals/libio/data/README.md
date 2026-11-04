# Test data for `shakelab.signals.libio`

| File | Content |
|---|---|
| `mseed/libmseed/` | Unmodified copy of the folder `test/data` of libmseed (EarthScope, [github.com/EarthScope/libmseed](https://github.com/EarthScope/libmseed), commit `d6e6ad3`, 5 August 2026), licensed under the Apache License 2.0 (`mseed/libmseed/LICENSE`): MiniSEED 2 and 3 files with the same series in every encoding and byte order (`reference-testdata-*`), text and detection records, historical encodings (CDSN, DWWSSN, GEOSCOPE, SRO), invalid blockette offsets, records without blockette 1000, records of mixed lengths in mixed order, an unapplied time correction, odd sampling rates, dates before 1970 and nanosecond times |
| `mseed/libmseed_reference.json` | What libmseed (same commit, `ms3_readtracelist`) reads from the MiniSEED 2 files of `mseed/libmseed` and of `tests/signals/data`: source identifier, start time, rate, number of samples and SHA-256 of the samples as big-endian float64, for each segment. Files that ShakeLab does not support (historical encodings, invalid blockettes, no blockette 1000) and segments without samples or with text are left out |
| `sac/obspy/` | SAC files from the test data of ObsPy ([github.com/obspy/obspy](https://github.com/obspy/obspy), `obspy/io/sac/tests/data`, commit `50033e1`, 28 August 2026), licensed under the GNU LGPL v3 (`sac/obspy/LICENSE.txt`), unmodified: `seism.sac` (the example distributed with SAC), `test.sac` and `test.sac.swap` (the same file in both byte orders), `LMOW.BHE.SAC`, `dis.G.SCZ.__.BHE_short`, `null_terminated.sac`, `non_ascii.sac`, `test_encode.sac` (text fields holding integers) and `testxy.sac` (alphanumeric SAC) |

The MiniSEED reference is generated with a small libmseed program and a Python
script kept outside the repository
(`AiSpace/Claude/Shakelab/tools/libmseed_reference/`); the script
docstring gives the commands.

Real waveforms used by these tests (station IV.MODE and the ESM
records of the 2012 Emilia earthquake) are in `tests/signals/data`,
with their sources and licenses.
