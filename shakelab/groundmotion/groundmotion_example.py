"""
Example: ground-motion evaluation with the ShakeLab GMPE provider.

This example demonstrates two complementary cases supported by the
context-based GMPE architecture:

1. BragatoSlejko2005
   - requires ML and epicentral distance (Repi);
   - Repi is computed automatically by the GMPE provider from the event
     epicentre and the site coordinates.

2. BindiEtAl2011
   - requires Mw, rake, Joyner-Boore distance (Rjb), and Vs30;
   - rake is read from ScenarioEvent.mechanism;
   - Rjb and Vs30 are supplied explicitly to the provider as site-dependent
     data.

The provider returns:
    (median intensity measure in linear space, sigma_ln)

The plotted +/- 1 sigma interval is:
    IM * exp(+/- sigma_ln)

Notes
-----
For the Bindi example, Rjb is assigned explicitly from the synthetic
distance vector used to construct the sites. This is intentional: the
current ScenarioEvent contains a point-source hypocentre but no finite-fault
geometry, so Rjb must not be silently approximated from epicentral distance.
"""

from __future__ import annotations

from math import exp

import matplotlib.pyplot as plt

from shakelab.groundmotion.providers import (
    GroundMotionContext,
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    MetricPoint,
    metric_to_wgs,
)


def build_sites_enu_east_line(
    origin: WgsPoint,
    distances_km: list[float],
) -> list[WgsPoint]:
    """
    Build WGS84 sites along the positive East direction.

    Parameters
    ----------
    origin
        Reference WGS84 point.
    distances_km
        Horizontal offsets from the origin in km.

    Returns
    -------
    list of WgsPoint
        Site coordinates converted back to WGS84.
    """
    frame = MetricFrame(
        ref_geo=origin,
        orientation="enu",
    )

    sites: list[WgsPoint] = []

    for distance_km in distances_km:
        point_metric = MetricPoint(
            float(distance_km) * 1000.0,
            0.0,
            0.0,
            frame=frame,
        )

        point_wgs = metric_to_wgs(
            point_metric,
            frame=frame,
        )

        sites.append(
            WgsPoint(
                longitude=point_wgs.longitude,
                latitude=point_wgs.latitude,
                elevation=0.0,
            )
        )

    return sites


def logspace_km(
    dmin_km: float,
    dmax_km: float,
    n: int,
) -> list[float]:
    """
    Return logarithmically spaced distances in km without NumPy.
    """
    if n < 2:
        raise ValueError("n must be >= 2.")

    if dmin_km <= 0.0 or dmax_km <= 0.0:
        raise ValueError(
            "dmin_km and dmax_km must be > 0."
        )

    if dmax_km <= dmin_km:
        raise ValueError(
            "dmax_km must be greater than dmin_km."
        )

    ratio = dmax_km / dmin_km

    return [
        dmin_km * ratio ** (index / (n - 1))
        for index in range(n)
    ]


def sigma_band(
    median: list[float],
    sigma_ln: list[float],
) -> tuple[list[float], list[float]]:
    """
    Return the lower and upper +/- 1 sigma curves in linear space.
    """
    lower = [
        value * exp(-sigma)
        for value, sigma in zip(
            median,
            sigma_ln,
        )
    ]

    upper = [
        value * exp(sigma)
        for value, sigma in zip(
            median,
            sigma_ln,
        )
    ]

    return lower, upper


def run_bragato_example(
    sites: list[WgsPoint],
) -> None:
    """
    Evaluate Bragato-Slejko (2005) using automatic Repi calculation.
    """
    event = ScenarioEvent(
        hypocentre=WgsPoint(
            longitude=13.10,
            latitude=46.05,
            elevation=-10000.0,
        ),
        magnitude=5.5,
        magnitude_type="ML",
    )

    provider = GroundMotionProvider.gmpe(
        gmpe_name="BragatoSlejko2005",
        distance_approx="ellipsoid",
    )

    ground_motion = GroundMotionContext(
        event=event,
        provider=provider,
    )

    imt = "PGA"

    median, sigma_ln = ground_motion.evaluate(
        imt=imt,
        sites=sites,
    )

    x_km = [
        event.epicentre.epicentral_distance_to(site)
        / 1000.0
        for site in sites
    ]

    lower, upper = sigma_band(
        median,
        sigma_ln,
    )

    plt.figure()
    plt.loglog(
        x_km,
        median,
        marker="o",
        linestyle="-",
        label=f"{imt} median",
    )
    plt.fill_between(
        x_km,
        lower,
        upper,
        alpha=0.25,
        label=f"{imt} +/- 1 sigma",
    )

    plt.xlabel("Epicentral distance Repi (km)")
    plt.ylabel(f"{imt} (g)")
    plt.title("Bragato & Slejko (2005)")
    plt.grid(
        True,
        which="both",
    )
    plt.legend()


def run_bindi_example(
    sites: list[WgsPoint],
    rjb_km: list[float],
) -> None:
    """
    Evaluate Bindi et al. (2011) using explicit Rjb and Vs30.
    """
    event = ScenarioEvent(
        hypocentre=WgsPoint(
            longitude=13.10,
            latitude=46.05,
            elevation=-10000.0,
        ),
        magnitude=6.0,
        magnitude_type="Mw",
        mechanism={
            "rake": 90.0,
        },
    )

    provider = GroundMotionProvider.gmpe(
        gmpe_name="BindiEtAl2011",
    )

    ground_motion = GroundMotionContext(
        event=event,
        provider=provider,
    )

    imt = "PGA"

    # One scalar Vs30 is intentionally used here. SitesContext broadcasts
    # it internally to all evaluation sites.
    site_parameters = {
        "vs30": 600.0,
    }

    # Rjb cannot be derived from the current point-source ScenarioEvent.
    # It is therefore supplied explicitly as one value per site.
    distances = {
        "rjb": rjb_km,
    }

    median, sigma_ln = ground_motion.evaluate(
        imt=imt,
        sites=sites,
        site_parameters=site_parameters,
        distances=distances,
    )

    lower, upper = sigma_band(
        median,
        sigma_ln,
    )

    plt.figure()
    plt.loglog(
        rjb_km,
        median,
        marker="o",
        linestyle="-",
        label=f"{imt} median",
    )
    plt.fill_between(
        rjb_km,
        lower,
        upper,
        alpha=0.25,
        label=f"{imt} +/- 1 sigma",
    )

    plt.xlabel("Joyner-Boore distance Rjb (km)")
    plt.ylabel(f"{imt} (g)")
    plt.title(
        "Bindi et al. (2011), "
        "Vs30 = 600 m/s, rake = 90 deg"
    )
    plt.grid(
        True,
        which="both",
    )
    plt.legend()


def main() -> None:
    """
    Run the Bragato-Slejko and Bindi examples.
    """
    origin = WgsPoint(
        longitude=13.10,
        latitude=46.05,
        elevation=0.0,
    )

    distances_km = logspace_km(
        dmin_km=1.0,
        dmax_km=200.0,
        n=50,
    )

    sites = build_sites_enu_east_line(
        origin,
        distances_km,
    )

    run_bragato_example(
        sites,
    )

    run_bindi_example(
        sites,
        rjb_km=distances_km,
    )

    plt.show()


if __name__ == "__main__":
    main()
