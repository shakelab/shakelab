"""Required GMPE site parameters remain explicit at evaluation time."""

from __future__ import annotations

import pytest

from shakelab.groundmotion.gmpe.abrahamson_et_al_2014 import (
    AbrahamsonEtAl2014,
)
from shakelab.groundmotion.gmpe.base import RuptureContext, SitesContext
from shakelab.groundmotion.gmpe.campbell_bozorgnia_2014 import (
    CampbellBozorgnia2014,
)
from shakelab.groundmotion.gmpe.chiou_youngs_2014 import ChiouYoungs2014


@pytest.mark.parametrize(
    ("gmpe_type", "required"),
    [
        (AbrahamsonEtAl2014, {"vs30", "vs30measured", "z1pt0"}),
        (ChiouYoungs2014, {"vs30", "vs30measured", "z1pt0"}),
        (CampbellBozorgnia2014, {"vs30", "z2pt5"}),
    ],
)
def test_missing_required_site_parameters_are_reported(
    gmpe_type,
    required,
):
    """No generic default is injected for absent basin parameters."""
    gmpe = gmpe_type()
    rupture = RuptureContext()
    sites = SitesContext(size=1, vs30=760.0)
    with pytest.raises(ValueError) as error:
        gmpe.validate_contexts(rupture, sites)

    message = str(error.value)
    assert gmpe_type.__name__ in message
    assert "sites:" in message
    for name in required - {"vs30"}:
        assert name in message


@pytest.mark.parametrize(
    ("gmpe_type", "basin"),
    [
        (ChiouYoungs2014, "z1pt0"),
        (CampbellBozorgnia2014, "z2pt5"),
    ],
)
def test_explicit_unavailable_basin_depth_is_not_missing(
    gmpe_type,
    basin,
):
    """A model-specific sentinel is distinct from an omitted field."""
    gmpe = gmpe_type()
    parameters = {"vs30": 760.0, basin: -999.0}
    if "vs30measured" in gmpe.required_site_parameters:
        parameters["vs30measured"] = False

    sites = SitesContext(size=1, **parameters)
    missing = gmpe.required_site_parameters - sites.keys()
    assert not missing
