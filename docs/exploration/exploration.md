[ShakeLab documentation](../README.md) › Exploration

# Exploration

`shakelab.exploration` collects the methods of geophysical exploration
that work on recordings of ambient vibrations and earthquakes: what
the wavefield crossing a group of stations reveals about the ground
beneath it and about the structures standing on it.

Its main component is **YARRAY**, the processing of seismic arrays. By
comparing simultaneous recordings in the frequency domain, it measures
the phase velocity and the direction of the waves crossing the array
(frequency-wavenumber beamforming, with the Bartlett, Capon and MUSIC
estimators) and the ellipticity of Rayleigh waves on three
components: the dispersion curves on which passive surface-wave
methods for site characterization are built. From the same
cross-spectral matrices, the frequency domain decomposition identifies
the resonances recorded by the stations, with their shapes and
damping. Since the quality of every result depends first on the
geometry of the array, YARRAY also helps to design an array before a
survey. It is used from Python, from the command line (`slab-yarray`)
and in real time with the SeedLink client of `signals`.

The package also contains spectral ratios of ambient-noise recordings
(`ambientnoise`) and multichannel analysis of surface waves (`masw`).

## Where to read

| To... | Read |
|---|---|
| process array data, plan a survey, identify resonances | [YARRAY guide](yarray.md) |
| see YARRAY at a glance | [YARRAY README](../../shakelab/exploration/yarray/README.md) |
| run complete examples | [examples](../../examples/README.md) |

YARRAY reads its data with `shakelab.signals` (see the
[signals overview](../signals/signals.md)) and draws its figures with
`shakelab.graphics.yarray`.
