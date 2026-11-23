"""Standalone GMPE examples for the principal ShakeLab source modes.

Run from the repository root with Matplotlib installed:

    python -m shakelab.groundmotion.groundmotion_example

Three independent plots demonstrate:

1. Native mode: BragatoSlejko2005 with automatically derived Repi.
2. Native mode: BindiEtAl2011 with an explicitly supplied Rjb.
3. The same AbrahamsonEtAl2014 GMPE with point, finite explicit,
   and magnitude-scaled finite rupture geometry.

All numbers are synthetic illustrations, not earthquake-specific
source calibrations or recommendations about fault geometry.

The provider returns linear-space median intensity measures and
sigma_ln, the standard deviation in natural-log space. Plotting
uses median * exp(+/- sigma_ln) when sigma is available.
"""

from __future__ import annotations

from math import exp, isfinite

import matplotlib.pyplot as plt

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint
from shakelab.libutils.geodeticN.transform import (
    MetricFrame,
    MetricPoint,
    metric_to_wgs,
)


def logspace_km(start: float, stop: float, count: int) -> list[float]:
    """Return positive logarithmically spaced distances in kilometres."""
    if count < 2 or not 0 < start < stop:
        raise ValueError("Require count >= 2 and 0 < start < stop.")
    return [
        start * (stop / start) ** (index / (count - 1))
        for index in range(count)
    ]


def eastward_sites(
    origin: WgsPoint,
    offsets_km: list[float],
) -> list[WgsPoint]:
    """Generate surface sites at positive East offsets from an epicentre."""
    frame = MetricFrame(ref_geo=origin, orientation="enu")
    sites = []
    for offset in offsets_km:
        point = MetricPoint(
            offset * 1000.0, 0.0, 0.0, frame=frame,
        )
        geographic = metric_to_wgs(point, frame=frame)
        sites.append(WgsPoint(
            geographic.longitude, geographic.latitude, 0.0,
        ))
    return sites


def draw_curve(
    ax,
    x_km: list[float],
    median: list[float],
    sigma_ln: list[float],
    label: str,
) -> None:
    """Draw median and (when available) its natural-log sigma band."""
    (line,) = ax.plot(x_km, median, label=label)
    if all(isfinite(sigma) and sigma >= 0 for sigma in sigma_ln):
        lower = [
            value * exp(-sigma)
            for value, sigma in zip(median, sigma_ln)
        ]
        upper = [
            value * exp(sigma)
            for value, sigma in zip(median, sigma_ln)
        ]
        ax.fill_between(
            x_km, lower, upper, alpha=0.15, color=line.get_color(),
        )


def native_repi_example(ax, offsets_km, sites) -> None:
    """Use the automatic epicentral distance in native mode."""
    event = ScenarioEvent.from_lonlat_depth(
        5.5, 13.10, 46.05, 10.0, magnitude_type="ML",
    )
    provider = GroundMotionProvider.gmpe(
        "BragatoSlejko2005",
        config={"source_geometry": "native"},
    )
    median, sigma_ln = provider.evaluate("PGA", sites, event)
    draw_curve(
        ax, offsets_km, median, sigma_ln,
        "Bragato–Slejko 2005 (native Repi)",
    )
    ax.set_title("Native: derived Repi")


def native_explicit_example(ax, offsets_km, sites) -> None:
    """Supply synthetic Rjb values explicitly in native mode."""
    event = ScenarioEvent.from_lonlat_depth(
        6.0, 13.10, 46.05, 10.0,
        magnitude_type="Mw",
        source_parameters={"rake": 90.0},
    )
    provider = GroundMotionProvider.gmpe(
        "BindiEtAl2011",
        config={"source_geometry": "native"},
    )
    # An illustrative supplied Rjb vector; it is not inferred from the
    # hypocentre and should not be treated as a real rupture distance.
    rjb_km = offsets_km
    median, sigma_ln = provider.evaluate(
        "PGA", sites, event,
        site_parameters={"vs30": 600.0},
        distances={"rjb": rjb_km},
    )
    draw_curve(
        ax, rjb_km, median, sigma_ln,
        "Bindi et al. 2011 (explicit Rjb)",
    )
    ax.set_title("Native: supplied Rjb")


def source_geometry_comparison(ax, offsets_km, sites) -> None:
    """Evaluate one GMPE and event with three source geometries."""
    event = ScenarioEvent.from_lonlat_depth(
        6.2, 13.10, 46.05, 7.0, magnitude_type="Mw",
    )
    site_parameters = {
        "vs30": 500.0,
        "vs30measured": True,
        "z1pt0": 300.0,
    }
    common_source = {
        "strike": 120.0,
        "dip": 35.0,
        "rake": 90.0,
    }

    point = GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "point",
            "source_parameters": {
                **common_source,
                "width": 12.0,
                "ztor": 2.0,
            },
        },
    )
    explicit = GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "finite",
            "source_parameters": {
                **common_source,
                "length": 20.0,
                "width": 12.0,
                "ztor": 2.0,
                "along_fraction": 0.5,
            },
        },
    )
    scaled = GroundMotionProvider.gmpe(
        "AbrahamsonEtAl2014",
        config={
            "source_geometry": "finite",
            "source_parameters": {
                **common_source,
                "aspect_ratio": 2.0,
                "lower_seismogenic_depth": 15.0,
            },
            "source_scaling": {
                "model": "ThingbaijamEtAl2017",
                "style": "auto",
            },
        },
    )

    for provider, label in [
        (point, "Point approximation"),
        (explicit, "Finite explicit"),
        (scaled, "Finite with scaling"),
    ]:
        median, sigma_ln = provider.evaluate(
            "PGA", sites, event, site_parameters=site_parameters,
        )
        draw_curve(ax, offsets_km, median, sigma_ln, label)

    source = scaled.resolve_source(event)
    print(
        "Scaled rupture: length={:.2f} km, width={:.2f} km, "
        "ztor={:.2f} km".format(
            source.rupture.length,
            source.rupture.width,
            source.rupture.ztor,
        )
    )
    ax.set_title("Same event and GMPE, three geometries")


def main() -> None:
    """Plot generic examples of ShakeLab source handling."""
    origin = WgsPoint(13.10, 46.05, 0.0)
    offsets_km = logspace_km(1.0, 200.0, 50)
    sites = eastward_sites(origin, offsets_km)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    native_repi_example(axes[0], offsets_km, sites)
    native_explicit_example(axes[1], offsets_km, sites)
    source_geometry_comparison(axes[2], offsets_km, sites)

    for index, ax in enumerate(axes):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel("PGA (g)")
        ax.set_xlabel(
            "Supplied Rjb (km)" if index == 1
            else "Epicentral East offset (km)"
        )
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize="small")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
