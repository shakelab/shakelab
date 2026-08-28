"""
Synthetic integration test for the ShakeMap reader and provider.

Prerequisites
-------------
Save:
- usgs_v1.py as shakelab/gmmodel/usgs.py
- groundmotion_v1.py as shakelab/gmmodel/groundmotion.py

Run:
    python3 test_shakemap_provider_v1.py
"""

from __future__ import annotations

import json
import tempfile

from pathlib import Path

from shakelab.gmmodel.groundmotion import (
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


def assert_close(value, expected, tolerance=1e-12):
    if abs(value - expected) > tolerance:
        raise AssertionError(
            f"Got {value}, expected {expected}."
        )


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        event_dir = root / "TEST-1"
        event_dir.mkdir()

        (event_dir / "grid.xml").write_text(
            GRID_XML,
            encoding="utf-8",
        )
        (event_dir / "uncertainty.xml").write_text(
            UNCERTAINTY_XML,
            encoding="utf-8",
        )
        (event_dir / "info.json").write_text(
            json.dumps(
                {
                    "input": {
                        "event_information": {
                            "event_id": "TEST-1",
                            "magnitude": 5.5,
                            "longitude": 13.5,
                            "latitude": 46.5,
                            "depth": 10.0,
                        }
                    }
                }
            ),
            encoding="utf-8",
        )

        provider = GroundMotionProvider.shakemap(
            config={
                "root_path": str(root),
                "interp": "linear",
                "outside": "raise",
                "acceleration_unit": "g",
                "velocity_unit": "cm/s",
            }
        )

        event = ScenarioEvent(
            magnitude=5.5,
            hypocentre=WgsPoint(
                longitude=13.5,
                latitude=46.5,
                elevation=-10000.0,
            ),
            event_id="TEST-1",
        )

        centre = WgsPoint(
            longitude=13.5,
            latitude=46.5,
            elevation=0.0,
        )

        pga, pga_sigma = provider.evaluate(
            "PGA",
            [centre],
            event,
        )
        assert_close(pga[0], 0.25)
        assert_close(pga_sigma[0], 0.4)

        pgv, pgv_sigma = provider.evaluate(
            "PGV",
            [centre],
            event,
        )
        assert_close(pgv[0], 2.5)
        assert_close(pgv_sigma[0], 0.5)

        sa03, sa03_sigma = provider.evaluate(
            "SA(0.3)",
            [centre],
            event,
        )
        assert_close(sa03[0], 0.35)
        assert_close(sa03_sigma[0], 0.6)

        # Cache must contain exactly one loaded event.
        if len(provider._cache) != 1:
            raise AssertionError(
                "ShakeMap event caching failed."
            )

        # Outside-domain behavior must fail explicitly.
        outside_site = WgsPoint(
            longitude=20.0,
            latitude=50.0,
            elevation=0.0,
        )

        try:
            provider.evaluate(
                "PGA",
                [outside_site],
                event,
            )
        except ValueError:
            pass
        else:
            raise AssertionError(
                "Outside-domain evaluation did not raise."
            )

        # ShakeMap provider requires event_id.
        no_id_event = ScenarioEvent(
            magnitude=5.5,
            hypocentre=WgsPoint(
                longitude=13.5,
                latitude=46.5,
                elevation=-10000.0,
            ),
        )

        try:
            provider.evaluate(
                "PGA",
                [centre],
                no_id_event,
            )
        except ValueError:
            pass
        else:
            raise AssertionError(
                "Missing event_id did not raise."
            )

    print("ShakeMap provider synthetic test passed.")


if __name__ == "__main__":
    main()
