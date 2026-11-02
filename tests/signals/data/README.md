# Test data for `shakelab.signals`

| File | Content |
|---|---|
| `MODE.mseed` | Strong-motion records (channels `IV.MODE..HNE`, `HNN`, `HNZ`, 100 Hz, counts) of the Mw 6.1 Emilia earthquake of 20 May 2012, time window 02:03:23 to 02:06:23 UTC |
| `MODE.xml` | FDSN StationXML 1.2 (level `response`) of the same channels, as distributed by the INGV fdsnws-station service on 2 November 2026 |
| `ESM/ASCII/IV.MODE..HN?.D.IT-2012-0008.{ACC,VEL,DIS}.MP.ASC` | Acceleration (cm/s**2), velocity (cm/s) and displacement (cm) of the same records processed by ESM (manual processing `MP`, 200 Hz), DYNA 1.2 format, ESM version 3.0, retrieved on 2 November 2026; the header gives the processing parameters (e.g. Butterworth band-pass of order 2, 0.03-40 Hz for HNE and HNZ, 0.02-40 Hz for HNN) |
| `ESM/MSEED/IV.MODE..HN?.IT-2012-0008.{ACC,VEL,DIS}.MP.mseed` | The same quantities in MiniSEED, retrieved in December 2023 from an earlier ESM version (differences from the ASCII files below 1 %) |

Station MODE (Modena), Italian Seismic Network (network code IV,
INGV, doi:[10.13127/SD/X0FXNH7QFY](https://doi.org/10.13127/SD/X0FXNH7QFY)).
Raw data and metadata retrieved from INGV through the FDSN web
services (EIDA node).

Processed data from the Engineering Strong-Motion database (ESM),
licensed CC BY 4.0: Mascandola C., Felicetta C., Russo E., Luzi L.,
Lanzano G., Sgobba S., Brunelli G., Ramadan F., Pacor F., ORFEUS
Strong-Motion SMC (2026). Engineering Strong Motion Database (ESM),
version 3.0. Istituto Nazionale di Geofisica e Vulcanologia (INGV).
[doi:10.13127/esm.3](https://doi.org/10.13127/esm.3).

All files are distributed with ShakeLab for testing and examples only,
with acknowledgement of the sources.
