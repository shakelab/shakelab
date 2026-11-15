# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""Check catalog values against published model parameters."""

import json
from math import hypot, inf

import numpy as np
import pytest

from shakelab.location.profiles import (
    list_models,
    load_model,
    model_info,
    read_model,
)
from shakelab.location.ttsolver import LayeredModel1D, solve_traveltime


def test_ogs_friuli_profile():
    """Reproduce Saraò et al. (2021), section 3.2."""
    model = load_model("ogs-friuli")
    assert isinstance(model, LayeredModel1D)
    np.testing.assert_array_equal(model.thickness, [22000, 17500, inf])
    np.testing.assert_array_equal(model.vp, [5850, 6800, 8000])
    np.testing.assert_allclose(
        model.vs, np.array([5850, 6800, 8000]) / 1.78, rtol=0, atol=1e-12,
    )


def test_iasp91_crust_and_moho_velocities():
    """Use the two IASP91 crustal layers and velocities at the Moho base."""
    model = load_model("iasp91-crust-halfspace")
    np.testing.assert_array_equal(model.z_top, [0, 20000, 35000])
    np.testing.assert_array_equal(model.vp, [5800, 6500, 8040])
    np.testing.assert_array_equal(model.vs, [3360, 3750, 4470])
    assert model.total_thickness() == inf


def test_ingv_rsn1_profile():
    """Reproduce model 5 from Peruzza et al. (2021), Table 2."""
    model = load_model("ingv-rsn1")
    np.testing.assert_array_equal(model.z_top, [0, 10000, 30000])
    np.testing.assert_array_equal(model.vp, [5000, 6000, 8100])
    np.testing.assert_allclose(
        model.vs, np.array([5000, 6000, 8100]) / 1.73, rtol=0, atol=1e-12,
    )


def test_profiles_are_discoverable_and_models_are_independent():
    """Expose provenance without sharing writable model arrays."""
    assert {"ogs-friuli", "iasp91-crust-halfspace", "ingv-rsn1"} <= set(
        list_models(),
    )
    name = "ogs-friuli"
    model = load_model(name)
    model.vp[0] = 1000
    assert load_model(name).vp[0] == 5850
    info = model_info(name)
    assert info["source_url"] == "https://doi.org/10.5194/essd-13-2245-2021"
    info["source_url"] = ""
    assert model_info(name)["source_url"]


def test_unknown_profile_lists_available_models():
    """Report a bad identifier along with the available choices."""
    with pytest.raises(ValueError, match="Unknown model.*ogs-friuli"):
        load_model("unknown")
    with pytest.raises(TypeError, match="string"):
        model_info(None)


def test_reference_model_passes_directly_to_solver():
    """A short ray above the first OGS interface is straight."""
    model = load_model("ogs-friuli")
    time, _ = solve_traveltime(
        model, 10000, 9000, phase="P", include_diving=False,
    )
    assert float(time) == pytest.approx(
        hypot(10000, 9000) / 5850, rel=0, abs=1e-12,
    )


@pytest.mark.parametrize("suffix", ["json", "csv"])
def test_read_custom_model(tmp_path, suffix):
    """Read equivalent models with a half-space and no density column."""
    path = tmp_path / f"my-model.{suffix.upper()}"
    if suffix == "json":
        contents = json.dumps({"layers": [
            {"thickness_m": 2000, "vp_m_s": 4000, "vs_m_s": 2200},
            {"thickness_m": None, "vp_m_s": 6000, "vs_m_s": 3200},
        ]})
    else:
        contents = (
            "thickness (m),vp (m/s),vs (m/s)\n"
            "2000,4000,2200\ninf,6000,3200\n"
        )
    path.write_text(contents, encoding="utf-8")

    model = read_model(path)
    np.testing.assert_array_equal(model.thickness, [2000, inf])
    np.testing.assert_array_equal(model.vp, [4000, 6000])
    np.testing.assert_array_equal(model.vs, [2200, 3200])


def test_read_model_uses_explicit_format(tmp_path):
    """Allow CSV files without a conventional file extension."""
    path = tmp_path / "model.txt"
    path.write_text(
        "thickness (m),vp (m/s),vs (m/s),density (kg/m3)\n"
        "1000,3000,1800,2200\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Unsupported model format"):
        read_model(path)
    assert read_model(str(path), format="csv").total_thickness() == 1000


@pytest.mark.parametrize(("contents", "message"), [
    ('{"layers": []}', "at least one layer"),
    ('{"layers": [{"thickness_m": Infinity}]}', "Invalid JSON number"),
    ('{"layers": [{"thickness_m": 1000, "vp_m_s": 0, "vs_m_s": 1800}]}',
     "finite and positive"),
    ('{"layers": [{"thickness_m": null, "vp_m_s": 4000, "vs_m_s": 2300},'
     '{"thickness_m": 1000, "vp_m_s": 4200, "vs_m_s": 2400}]}',
     "Layer 1: invalid thickness"),
])
def test_read_model_rejects_invalid_layers(tmp_path, contents, message):
    """Reject malformed or physically invalid JSON layer data."""
    path = tmp_path / "bad-model.json"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_model(path)
