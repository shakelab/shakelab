# Test data for `shakelab.exploration.yarray`

| File | Content |
|---|---|
| `sesame/M02.1_3001_*_0000.{1,2,3}.sac` | Nine stations of a synthetic ambient-noise array simulation of the SESAME project (model `M2.1`), components Z, N, E (files `.1`, `.2`, `.3`), 114.29 Hz (`DELTA` 0.00875 s); excerpt of 60 s (samples 11429 to 18286 of the original 405 s records, `B` 100 s) |
| `sesame_array.json` | Array file of the nine stations (local coordinates, one file per component) |
| `reference.json` | Maxima of the YArray code of 12 August 2026 on the nine stations, for every component mode (vertical, radial, transverse, three components, Rayleigh joint), with their configurations |

The nine stations are a subset of the 38 of the simulation, chosen
with `yarray.response` among random subsets containing the central
station and the first three: `S1001` (2048, 1964), `S1002` (2028,
2000), `S1003` (2060, 2008), `S1006` (2000, 2028), `S1017` (2036,
2048), `S1019` (2048, 2048, centre), `S1021` (2060, 2048), `S1023`
(2116, 2048), `S1028` (2032, 2056), in metres. At level 0.5 the
resolved and alias-free wavenumbers span 0.030 to 0.240 rad/m, while
no wavenumber is both resolved and alias-free with the first three
stations alone (the excerpt of the first version of these data).

The files keep the header of the original simulation (binary SAC,
`KSTNM` = station, `KCMPNM` = component, no network code), apart from
`NPTS`, `B`, `E`, `DEPMIN`, `DEPMAX` and `DEPMEN`, updated for the
excerpt, and `KEVNM`, written as one undefined 16-character field.
The station coordinates are local, in metres: East in
`USER7`, North in `USER8`, elevation in `USER9`; the same values are
also written in `STLO` and `STLA`, which therefore do not contain
geographic coordinates (the SAC reader ignores them with a warning).

The polarity of the components is not documented (no `CMPAZ` or
`CMPINC`). The signed Rayleigh ellipticity of the joint processing
(positive for prograde motion with the vertical positive upwards) has
the opposite of the expected sequence of a soft layer over a stiff
half-space on the 38 stations of the simulation (retrograde between
the H/V peak and the trough near 4 Hz, prograde above): the vertical
component is probably positive downwards (or the horizontal ones
reversed). The tests do not depend on it.

The simulation was produced within the European project SESAME (Site
EffectS assessment using AMbient Excitations, 2001-2004) and
distributed with the array-processing material of ETH Zurich
(ArrayTool3C file list, `NAME = M2.1`, `LOCATION = ETH Zurich`).
The files are distributed with ShakeLab for testing and examples only,
with acknowledgement of the source.
