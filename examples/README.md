# ShakeLab examples

Scripts showing how to use ShakeLab. The folder mirrors the structure
of the package and of the tests: the examples of `shakelab.signals`
are in `examples/signals/`, and so on.

Run the examples from the root of the repository, e.g.

```bash
python examples/signals/response_correction_esm.py
```

They use the data in `tests/<subpackage>/data/` and require
`matplotlib` for the figures. They are references for the user, short
and commented step by step; some need a server or data of the user
and do not run as they are. Files written by the examples go to
`output/` (ignored by git).

## Signals

| Example | Content |
|---|---|
| [`signals/response_correction_esm.py`](signals/response_correction_esm.py) | Instrument response correction of the IV.MODE strong-motion records (StationXML), compared with the ESM processed accelerations (DYNA ASCII files, with the ESM filter band of each component) |

## Exploration

Array processing (`shakelab.exploration.yarray`, guide in
[`docs/yarray.md`](../docs/yarray.md)), on nine three-component
stations of the SESAME simulation:

| Example | Content |
|---|---|
| [`exploration/yarray_01_array.py`](exploration/yarray_01_array.py) | Building the array data: array file, waveforms with a table of coordinates (channel and time selection), records in memory |
| [`exploration/yarray_02_config.py`](exploration/yarray_02_config.py) | The configuration: `BeamformingConfig` from Python, other frequency, window and grid forms, `replace`, JSON |
| [`exploration/yarray_03_beamforming.py`](exploration/yarray_03_beamforming.py) | A vertical Capon analysis: `beamform` and `Beamformer`, the maxima, dispersion and spectrum figures |
| [`exploration/yarray_04_methods.py`](exploration/yarray_04_methods.py) | Bartlett, Capon and MUSIC on the same data, with their options |
| [`exploration/yarray_05_components.py`](exploration/yarray_05_components.py) | Radial, transverse and the three components separately |
| [`exploration/yarray_06_rtbf.py`](exploration/yarray_06_rtbf.py) | Rayleigh three-component beamforming with the signed ellipticity |
| [`exploration/yarray_07_transforms.py`](exploration/yarray_07_transforms.py) | Direct transform and RFFT (nearest bins, interpolation, all the bins of a range), windows of cycles |
| [`exploration/yarray_08_uncertainty.py`](exploration/yarray_08_uncertainty.py) | Empirical standard error of the power, uncertain station coordinates |
| [`exploration/yarray_09_results.py`](exploration/yarray_09_results.py) | Results: maxima, spectra, joining consecutive results, files |
| [`exploration/yarray_10_density.py`](exploration/yarray_10_density.py) | Density of the maxima: summary, peaks, figure |
| [`exploration/yarray_11_command_line.sh`](exploration/yarray_11_command_line.sh) | `slab-yarray` and `slab-yarray response`, command by command (shell script) |
| [`exploration/yarray_12_realtime.py`](exploration/yarray_12_realtime.py) | Real-time monitoring of an array with the SeedLink client and `ArrayMonitor` (needs a server) |
| [`exploration/yarray_13_response.py`](exploration/yarray_13_response.py) | Geometric response of a planned layout and of the test array: limits, velocities, search domain, uncertain positions |
| [`exploration/yarray_14_performance.py`](exploration/yarray_14_performance.py) | Monte Carlo performance of the estimators; SNR-windows sweep |
| [`exploration/yarray_15_operational.py`](exploration/yarray_15_operational.py) | Frequency-velocity envelope with criteria, limiting factors and recommendations |
| [`exploration/yarray_16_profiling.py`](exploration/yarray_16_profiling.py) | Processing times of the main choices of a configuration |
