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
| [`exploration/yarray_sesame.py`](exploration/yarray_sesame.py) | Frequency-wavenumber analysis of three stations of the SESAME array simulation (array file with local coordinates): vertical Capon beamforming with the RFFT, Rayleigh three-component beamforming with signed ellipticity, and the same analysis window by window with the real-time `ArrayMonitor` |
