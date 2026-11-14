# ShakeLab examples

Scripts showing how to use ShakeLab. The folder mirrors the structure
of the package and of the tests: the examples of `shakelab.signals`
are in `examples/signals/`, and so on.

The examples illustrate the interfaces, the configurations and the
functions available: they are references for the user, short and
direct, with comments on the options. File names (`array.json`,
`data/*.mseed`, `array.xml`...) and servers are generic: replace them
with your own data to run an example. The figures need `matplotlib`.

## Signals

| Example | Content |
|---|---|
| [`signals/response_correction_esm.py`](signals/response_correction_esm.py) | Instrument response correction of the IV.MODE strong-motion records (StationXML), compared with the ESM processed accelerations (DYNA ASCII files, with the ESM filter band of each component) |

## Exploration

Array processing (`shakelab.exploration.yarray`, guide in
[`docs/exploration/yarray.md`](../docs/exploration/yarray.md)):

| Example | Content |
|---|---|
| [`exploration/yarray/01_array.py`](exploration/yarray/01_array.py) | Building the array data: array file, waveforms with metadata (channel and time selection) or with local coordinates, one station at a time |
| [`exploration/yarray/02_config.py`](exploration/yarray/02_config.py) | The configuration: single parameters, `update` and attributes, objects and dictionaries, JSON files with changes |
| [`exploration/yarray/03_beamforming.py`](exploration/yarray/03_beamforming.py) | A vertical Capon analysis: `Beamformer`, a time window, the maxima, the dispersion curve, figures |
| [`exploration/yarray/04_methods.py`](exploration/yarray/04_methods.py) | Bartlett, Capon and MUSIC on the same data, with their options |
| [`exploration/yarray/05_components.py`](exploration/yarray/05_components.py) | Radial, transverse and the three components separately |
| [`exploration/yarray/06_rtbf.py`](exploration/yarray/06_rtbf.py) | Rayleigh three-component beamforming with the signed ellipticity |
| [`exploration/yarray/07_transforms.py`](exploration/yarray/07_transforms.py) | Direct transform and RFFT (nearest bins, interpolation, all the bins of a range), windows of cycles |
| [`exploration/yarray/08_uncertainty.py`](exploration/yarray/08_uncertainty.py) | Empirical standard error of the power, uncertain station coordinates |
| [`exploration/yarray/09_results.py`](exploration/yarray/09_results.py) | Results: maxima, spectra, dispersion curve, joining successive windows, files |
| [`exploration/yarray/10_density.py`](exploration/yarray/10_density.py) | Density of the maxima: summary, peaks, figure |
| [`exploration/yarray/11_command_line.sh`](exploration/yarray/11_command_line.sh) | `slab-yarray` and `slab-yarray response`, command by command (shell script) |
| [`exploration/yarray/12_realtime.py`](exploration/yarray/12_realtime.py) | Real-time monitoring of an array with the SeedLink client and `ArrayMonitor` |
| [`exploration/yarray/13_response.py`](exploration/yarray/13_response.py) | Geometric response of an array or of a planned layout: limits, velocities, search domain, uncertain positions |
| [`exploration/yarray/14_performance.py`](exploration/yarray/14_performance.py) | Monte Carlo performance of the estimators; SNR-windows sweep |
| [`exploration/yarray/15_operational.py`](exploration/yarray/15_operational.py) | Frequency-velocity envelope with criteria, limiting factors and recommendations |
| [`exploration/yarray/16_profiling.py`](exploration/yarray/16_profiling.py) | Processing times of the main choices of a configuration |
