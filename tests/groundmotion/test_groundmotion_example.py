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
        assert len(figures) == 1
        figure = plt.figure(figures[0])
        assert len(figure.axes) == 3
        assert all(axis.lines for axis in figure.axes)
    finally:
        plt.close("all")
