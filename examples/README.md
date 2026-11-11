# ShakeLab examples

Scripts showing how to use ShakeLab. The folder mirrors the structure
of the package and of the tests: the examples of `shakelab.signals`
are in `examples/signals/`, and so on.

Run the examples from the root of the repository, e.g.

```bash
python examples/signals/response_correction_esm.py
```

They use the data in `tests/<subpackage>/data/` and require
`matplotlib` for the figures. Every example is also run by the test
suite (`tests/test_examples.py`, with a non-interactive backend), so
that it stays in line with the library.

## Signals

| Example | Content |
|---|---|
| [`signals/response_correction_esm.py`](signals/response_correction_esm.py) | Instrument response correction of the IV.MODE strong-motion records (StationXML), compared with the ESM processed accelerations (DYNA ASCII files, with the ESM filter band of each component) |

## Exploration

| Example | Content |
|---|---|
| [`exploration/yarray_sesame.py`](exploration/yarray_sesame.py) | Quick tour of the array processing on nine three-component stations of the SESAME simulation: vertical Capon beamforming with the RFFT, Rayleigh three-component beamforming with signed ellipticity, the same analysis window by window with the real-time `ArrayMonitor`, density of the maxima |
| [`exploration/yarray_processing.py`](exploration/yarray_processing.py) | The array processing as a library: three ways of building the array, configuration in Python and JSON, Bartlett/Capon/MUSIC, every component mode (radial, transverse, three components, RTBF), empirical and coordinate uncertainties, results (maxima, spectra, joining, files, density) |
| [`exploration/yarray_command_line.py`](exploration/yarray_command_line.py) | The command `slab-yarray`: processing with a configuration file and its `graphics` section, options replacing the file, RTBF with `--write-config`, appended results, figures of saved results, `slab-yarray response` on the array and on a planned layout |
| [`exploration/yarray_realtime.py`](exploration/yarray_realtime.py) | Real-time monitoring: an `ArrayMonitor` called by the client every 10 s of data, processing in its worker thread while the main thread queries it; replay of the SESAME records, or a live SeedLink server with `--server` |
| [`exploration/yarray_planning.py`](exploration/yarray_planning.py) | Array design and survey planning: geometric response and velocity limits of two layouts, coordinate uncertainty, search domain, Monte Carlo performance of the estimators, SNR-windows sweep, frequency-velocity envelope with criteria and recommendations |
| [`exploration/yarray_profiling.py`](exploration/yarray_profiling.py) | Processing times of the choices of a configuration (estimators, components, direct transform versus RFFT, RTBF, cache, planning); `--full` for realistic sizes |
