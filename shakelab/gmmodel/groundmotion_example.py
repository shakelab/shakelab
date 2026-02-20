"""
Example: GMPE decay with the GMPE distance metric (log-log + ±1σ).

- Sites are generated along +East in a local ENU frame centered at the
  event epicentre, then converted back to WGS84.
- X-axis distance is computed using the *same* distance metric declared by
  the GMPE (provider._gmpe.DISTANCE_METRIC).
- The provider returns IM median (linear) and sigma_ln; ±1σ band is
  IM * exp(±sigma_ln).
"""

from __future__ import annotations

from math import exp

import matplotlib.pyplot as plt

from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    MetricPoint,
    metric_to_wgs,
)

from shakelab.engineering.groundmotion import (
    GroundMotionContext,
    GroundMotionProvider,
    ScenarioEvent
)


def build_sites_enu_east_line(
    origin: WgsPoint,
    distances_km: list[float],
) -> list[WgsPoint]:
    """
    Build WGS84 sites at given distances along the +East axis in a local ENU
    frame centered at `origin`.
    """
    frame = MetricFrame(ref_geo=origin, orientation="enu")
    sites: list[WgsPoint] = []

    for d_km in distances_km:
        mp = MetricPoint(float(d_km) * 1000.0, 0.0, 0.0, frame=frame)
        p = metric_to_wgs(mp, frame=frame)
        sites.append(
            WgsPoint(longitude=p.longitude, latitude=p.latitude, elevation=0.0)
        )
    return sites


def logspace_km(dmin_km: float, dmax_km: float, n: int) -> list[float]:
    """
    Log-spaced distances in km (no numpy).

    distances[i] = dmin * (dmax/dmin)^(i/(n-1))
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    if dmin_km <= 0.0 or dmax_km <= 0.0:
        raise ValueError("dmin_km and dmax_km must be > 0")
    if dmax_km <= dmin_km:
        raise ValueError("dmax_km must be > dmin_km")

    ratio = dmax_km / dmin_km
    return [dmin_km * (ratio ** (i / (n - 1))) for i in range(n)]


def main() -> None:
    # --------------------------------------------------------------
    # Event definition (hypocentre elevation negative = depth in meters)
    # --------------------------------------------------------------
    event = ScenarioEvent(
        hypocentre=WgsPoint(
            longitude=13.10,
            latitude=46.05,
            elevation=-10e4,  # 10 km depth
        ),
        magnitude=6.5,
    )

    # --------------------------------------------------------------
    # Provider (GMPE resolved via registry.json)
    # --------------------------------------------------------------
    distance_approx = "ellipsoid"

    provider = GroundMotionProvider.gmpe(
        gmpe_name="BragatoSlejko2005",
        distance_approx="ellipsoid",
    )

    gm = GroundMotionContext(event=event, provider=provider)

    # --------------------------------------------------------------
    # Build sites (log-spaced epicentral distances)
    # Note: if GMPE metric is hypocentral, distances will be bounded
    # below by depth (e.g. ~10 km here). That's expected.
    # --------------------------------------------------------------
    n = 50
    dmin_km = 1
    dmax_km = 200.0
    distances_km = logspace_km(dmin_km, dmax_km, n)

    sites = build_sites_enu_east_line(event.epicentre, distances_km)

    # --------------------------------------------------------------
    # X-axis: use the GMPE-declared metric
    # --------------------------------------------------------------
    metric = getattr(provider, "_gmpe").DISTANCE_METRIC

    if metric in {"rhypo", "hypocentral"}:
        x_km = [
            event.hypocentre.hypocentral_distance_to(p) / 1000.0
            for p in sites
        ]
        x_label = f"Distance (km) [hypocentral, {distance_approx}]"

    elif metric in {"repi", "epicentral"}:
        x_km = [
            event.epicentre.epicentral_distance_to(p) / 1000.0
            for p in sites
        ]
        x_label = "Distance (km) [epicentral]"

    else:
        raise ValueError(f"Unsupported GMPE distance metric: {metric!r}")

    # --------------------------------------------------------------
    # Ground motion evaluation
    # --------------------------------------------------------------
    imt = "PGA"
    im_med, sigma_ln = gm.evaluate(imt=imt, sites=sites)

    im_lo = [m * exp(-s) for m, s in zip(im_med, sigma_ln)]
    im_hi = [m * exp(+s) for m, s in zip(im_med, sigma_ln)]

    # --------------------------------------------------------------
    # Plot (log-log) + sigma band
    # --------------------------------------------------------------
    plt.figure()
    plt.loglog(x_km, im_med, marker="o", linestyle="-", label=f"{imt} median")
    plt.fill_between(x_km, im_lo, im_hi, alpha=0.25, label=f"{imt} ±1σ")

    plt.xlabel(x_label)
    plt.ylabel(f"{imt} (linear)")
    plt.title(f"GMPE decay: {imt} vs distance")
    plt.grid(True, which="both")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()