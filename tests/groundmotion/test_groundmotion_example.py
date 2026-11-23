"""Smoke test for the standalone ground-motion plotting example."""

import pytest


def test_groundmotion_example_runs_with_optional_matplotlib(monkeypatch):
    """Evaluate the five illustrated modes and construct three panels."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    from shakelab.groundmotion import groundmotion_example

    monkeypatch.setattr(plt, "show", lambda: None)
    try:
        groundmotion_example.main()
        figures = plt.get_fignums()
        assert len(figures) == 2
        plots = [plt.figure(number) for number in figures]
        assert sorted(len(fig.axes) for fig in plots) == [3, 6]
        comparison = next(fig for fig in plots if len(fig.axes) == 3)
        diagnostic = next(fig for fig in plots if len(fig.axes) == 6)
        assert all(axis.lines for axis in comparison.axes)
        assert all(axis.lines for axis in diagnostic.axes[:5])
    finally:
        plt.close("all")


def test_diagnostic_distances_match_finite_geometry():
    """Recorded distance inputs match the independently exposed rupture."""
    pytest.importorskip("matplotlib")
    from matplotlib import pyplot as plt
    from shakelab.groundmotion import groundmotion_example as example
    from shakelab.libutils.geodeticN.primitives import WgsPoint

    offsets = [1.0, 3.0, 8.0, 20.0, 50.0]
    sites = example.eastward_sites(
        WgsPoint(13.10, 46.05, 0.0), offsets,
    )
    fig, ax = plt.subplots()
    try:
        result = example.source_geometry_comparison(ax, offsets, sites)
        cases = result["cases"]
        for kind in ("explicit", "scaled"):
            geometry = cases[kind]["source"].rupture
            distances = geometry.distances(sites)
            for metric in ("rrup", "rjb", "rx", "ry0"):
                record = cases[kind]["trace"]["distances"][metric]
                assert record["origin"] == "finite_rupture"
                assert record["values_km"] == pytest.approx(
                    distances[metric],
                )
        for metric in ("rrup", "rjb", "rx", "ry0"):
            assert cases["point"]["trace"]["distances"][metric][
                "origin"
            ] == "point_approximation"
        assert cases["point"]["source"].rupture is None
        assert cases["explicit"]["source"].rupture.width == pytest.approx(
            12.0,
        )
    finally:
        plt.close(fig)
