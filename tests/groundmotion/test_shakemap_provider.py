"""Integration tests for the USGS ShakeMap XML-backed provider.

The XML products are deliberately small and synthetic; no remote data or
external ShakeMap service is required.
"""

from __future__ import annotations

import json
from math import isnan

import pytest

from shakelab.groundmotion.providers import (
    GroundMotionProvider,
    ScenarioEvent,
)
from shakelab.libutils.geodeticN.primitives import WgsPoint


GRID_XML = """<?xml version="1.0" encoding="US-ASCII"?>
<shakemap_grid xmlns="http://earthquake.usgs.gov/eqcenter/shakemap">
  <grid_specification
      lon_min="13.0"
      lon_max="14.0"
      lat_min="46.0"
      lat_max="47.0"
      nominal_lon_spacing="1.0"
      nominal_lat_spacing="1.0"
      nlon="2"
      nlat="2" />
  <grid_field index="1" name="LON" units="dd" />
  <grid_field index="2" name="LAT" units="dd" />
  <grid_field index="3" name="PGA" units="%g" />
  <grid_field index="4" name="PGV" units="cm/s" />
  <grid_field index="5" name="PSA03" units="%g" />
  <grid_field index="6" name="PSA10" units="%g" />
  <grid_field index="7" name="PSA30" units="%g" />
  <grid_data>
13.0 46.0 10.0 1.0 20.0 15.0 5.0
14.0 46.0 20.0 2.0 30.0 25.0 10.0
13.0 47.0 30.0 3.0 40.0 35.0 15.0
14.0 47.0 40.0 4.0 50.0 45.0 20.0
  </grid_data>
</shakemap_grid>
"""

UNCERTAINTY_XML = """<?xml version="1.0" encoding="US-ASCII"?>
<shakemap_grid xmlns="http://earthquake.usgs.gov/eqcenter/shakemap">
  <grid_specification
      lon_min="13.0"
      lon_max="14.0"
      lat_min="46.0"
      lat_max="47.0"
      nominal_lon_spacing="1.0"
      nominal_lat_spacing="1.0"
      nlon="2"
      nlat="2" />
  <grid_field index="1" name="LON" units="dd" />
  <grid_field index="2" name="LAT" units="dd" />
  <grid_field index="3" name="STDPGA" units="ln(g)" />
  <grid_field index="4" name="STDPGV" units="ln(cm/s)" />
  <grid_field index="5" name="STDPSA03" units="ln(g)" />
  <grid_field index="6" name="STDPSA10" units="ln(g)" />
  <grid_field index="7" name="STDPSA30" units="ln(g)" />
  <grid_data>
13.0 46.0 0.4 0.5 0.6 0.7 0.8
14.0 46.0 0.4 0.5 0.6 0.7 0.8
13.0 47.0 0.4 0.5 0.6 0.7 0.8
14.0 47.0 0.4 0.5 0.6 0.7 0.8
  </grid_data>
</shakemap_grid>
"""


@pytest.fixture
def shakemap_root(tmp_path):
    """Build one synthetic four-node ShakeMap product set."""
    event_dir = tmp_path / "TEST-1"
    event_dir.mkdir()
    (event_dir / "grid.xml").write_text(
        GRID_XML, encoding="utf-8",
    )
    (event_dir / "uncertainty.xml").write_text(
        UNCERTAINTY_XML, encoding="utf-8",
    )
    (event_dir / "info.json").write_text(
        json.dumps({
            "input": {
                "event_information": {
                    "event_id": "TEST-1",
                    "magnitude": 5.5,
                    "longitude": 13.5,
                    "latitude": 46.5,
                    "depth": 10.0,
                },
            },
        }),
        encoding="utf-8",
    )
    return tmp_path


@pytest.fixture
def event():
    """Synthetic event associated with the XML files."""
    return ScenarioEvent.from_lonlat_depth(
        5.5, 13.5, 46.5, 10.0, event_id="TEST-1",
    )


@pytest.fixture
def provider(shakemap_root):
    """Use the public factory and unit conventions."""
    return GroundMotionProvider.shakemap(
        config={
            "root_path": str(shakemap_root),
            "interp": "linear",
            "outside": "raise",
            "acceleration_unit": "g",
            "velocity_unit": "cm/s",
        },
    )


@pytest.mark.parametrize(
    "imt, expected, expected_sigma",
    [
        ("PGA", 0.25, 0.4),
        ("PGV", 2.5, 0.5),
        ("SA(0.3)", 0.35, 0.6),
        ("SA(1.0)", 0.25, 0.7),
        ("SA(3.0)", 0.125, 0.8),
    ],
)
def test_bilinear_interpolation_and_units(
    provider, event, imt, expected, expected_sigma,
):
    """Interpolation at the centre averages all four grid corners."""
    centre = WgsPoint(13.5, 46.5, 0.0)
    values, sigmas = provider.evaluate(imt, [centre], event)
    assert values == pytest.approx([expected])
    assert sigmas == pytest.approx([expected_sigma])


def test_multiple_sites_and_nearest_interpolation(provider, event):
    """Nearest interpolation returns the nearest-node values."""
    sites = [
        WgsPoint(13.1, 46.1, 0.0),
        WgsPoint(13.9, 46.9, 0.0),
    ]
    values, sigmas = provider.evaluate(
        "PGA", sites, event, interp="nearest",
    )
    assert values == pytest.approx([0.1, 0.4])
    assert sigmas == pytest.approx([0.4, 0.4])


def test_event_cache_and_public_invalidation(provider, event):
    """The reader loads each event once unless cache is cleared."""
    centre = WgsPoint(13.5, 46.5, 0.0)
    provider.evaluate("PGA", [centre], event)
    provider.evaluate("PGV", [centre], event)
    assert len(provider._cache) == 1
    cached = provider._cache["TEST-1"]
    provider.evaluate("PGA", [centre], event)
    assert provider._cache["TEST-1"] is cached

    provider.clear_cache("TEST-1")
    assert provider._cache == {}
    provider.evaluate("PGA", [centre], event)
    assert len(provider._cache) == 1


def test_outside_grid_policy(provider, event):
    outside = WgsPoint(20.0, 50.0, 0.0)
    with pytest.raises(ValueError):
        provider.evaluate("PGA", [outside], event)

    values, sigmas = provider.evaluate(
        "PGA", [outside], event, outside="nan",
    )
    assert isnan(values[0])
    # The interpolation policy also applies to uncertainty fields.
    assert isnan(sigmas[0])


def test_event_id_is_required(provider):
    anonymous = ScenarioEvent.from_lonlat_depth(
        5.5, 13.5, 46.5, 10.0,
    )
    with pytest.raises(ValueError, match="event_id"):
        provider.evaluate(
            "PGA", [WgsPoint(13.5, 46.5, 0.0)], anonymous,
        )


def test_required_uncertainty_is_enforced(
    shakemap_root, event,
):
    (shakemap_root / "TEST-1" / "uncertainty.xml").unlink()
    strict_provider = GroundMotionProvider.shakemap(
        config={
            "root_path": str(shakemap_root),
            "require_uncertainty": True,
        },
    )
    with pytest.raises(KeyError, match="STDPGA"):
        strict_provider.evaluate(
            "PGA", [WgsPoint(13.5, 46.5, 0.0)], event,
        )


def test_unknown_imt_is_rejected(provider, event):
    with pytest.raises(ValueError, match="Unsupported ShakeMap IMT"):
        provider.evaluate(
            "SA(0.5)", [WgsPoint(13.5, 46.5, 0.0)], event,
        )
