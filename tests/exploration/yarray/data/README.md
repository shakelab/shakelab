# Test data for `shakelab.exploration.yarray`

| File | Content |
|---|---|
| `sesame/M02.1_3001_*_0000.{1,2,3}.sac` | Three stations (`S1001`, `S1002`, `S1003`) of a synthetic ambient-noise array simulation of the SESAME project (model `M2.1`), components Z, N, E (files `.1`, `.2`, `.3`), 114.29 Hz (`DELTA` 0.00875 s); excerpt of 60 s (samples 11429 to 18286 of the original 405 s records, `B` 100 s) |

The files keep the header of the original simulation (binary SAC,
`KSTNM` = station, `KCMPNM` = component, no network code), apart from
`NPTS`, `B`, `E`, `DEPMIN`, `DEPMAX` and `DEPMEN`, updated for the
excerpt. The station coordinates are local, in metres: East in
`USER7`, North in `USER8`, elevation in `USER9`; the same values are
also written in `STLO` and `STLA`, which therefore do not contain
geographic coordinates (the SAC reader ignores them with a warning).

The simulation was produced within the European project SESAME (Site
EffectS assessment using AMbient Excitations, 2001-2004) and
distributed with the array-processing material of ETH Zurich
(ArrayTool3C file list, `NAME = M2.1`, `LOCATION = ETH Zurich`).
The files are distributed with ShakeLab for testing and examples only,
with acknowledgement of the source.
