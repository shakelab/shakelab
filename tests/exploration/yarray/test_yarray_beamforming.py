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
"""
Unit tests for :mod:`shakelab.exploration.yarray.beamforming`.

Synthetic plane waves of known velocity, azimuth, polarization and
ellipticity are recovered by every estimator and component mode; the
beam power of a random plane wave equals its power spectral density;
the results on the SESAME excerpt are compared with the independent
YArray implementation of August 2026 (same maxima, power scaled by the
new spectral-density normalization).
"""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
    BeamformingWarning,
    beamform,
)
from shakelab.exploration.yarray.config import (
    BeamformingConfig,
    SearchGrid,
    WindowConfig,
)
from shakelab.exploration.yarray.data import ArrayData, StationData
from shakelab.exploration.yarray.io import read_array
from shakelab.signals import fourier
from shakelab.signals.base import Record
from shakelab.signals.crossspectral import coefficients
from shakelab.signals.io import reader

DATA = Path(__file__).parent / "data"
DELTA = 0.01
START = "2026-01-01T00:00:00"
VELOCITY = 400.0
AZIMUTH = 60.0


def _coordinates(nsta=10, aperture=200.0, seed=1):
    """Random station positions."""
    return np.random.default_rng(seed).uniform(
        -aperture / 2, aperture / 2, (nsta, 2)
    )


def _array(motion, coordinates, units=None, comps="ZNE"):
    """ArrayData from motion[i] = {"Z": ..., "N": ..., "E": ...}."""
    stations = []

    for i, (xy, channels) in enumerate(zip(coordinates, motion)):
        records = [
            Record(data=channels[code], delta=DELTA, time=START,
                   sid=f"XX.S{i:02d}..HH{code}", units=units)
            for code in comps
        ]
        stations.append(StationData(f"XX.S{i:02d}.", records, xy))

    return ArrayData(stations)


def _delays(coordinates, velocity=VELOCITY, azimuth=AZIMUTH):
    """Plane-wave delays x . s at the stations."""
    angle = math.radians(azimuth)
    slowness = np.array([math.sin(angle), math.cos(angle)]) / velocity
    return coordinates @ slowness


def _harmonic_array(freqs=(4.0, 6.0), kind="rayleigh", ellipticity=0.5,
                    noise=0.3, nsamp=12000, seed=2):
    """
    Harmonic plane waves with noise.

    ``rayleigh``: vertical cos, radial ``ellipticity * sin`` (positive
    ellipticity in the RTBF convention); ``love``: transverse cos.
    """
    rng = np.random.default_rng(seed)
    coordinates = _coordinates()
    time = np.arange(nsamp) * DELTA
    angle = math.radians(AZIMUTH)
    motion = []

    for tau in _delays(coordinates):
        z = sum(np.cos(2 * np.pi * f * (time - tau)) for f in freqs)
        h = sum(np.sin(2 * np.pi * f * (time - tau)) for f in freqs)

        if kind == "rayleigh":
            vertical, east, north = (
                z, ellipticity * h * math.sin(angle),
                ellipticity * h * math.cos(angle),
            )
        else:
            vertical = np.zeros(nsamp)
            east, north = z * math.cos(angle), -z * math.sin(angle)

        motion.append({
            code: value + noise * rng.standard_normal(nsamp)
            for code, value in (("Z", vertical), ("N", north), ("E", east))
        })

    return _array(motion, coordinates)


def _config(**kwargs):
    """Configuration with a fine grid around the synthetic waves."""
    options = dict(
        frequencies=[4.0, 6.0],
        grid=SearchGrid(200.0, 800.0, 61, scale="linear",
                        azimuth_size=72),
        window=WindowConfig(type="fixed", length=4.0, windows_per_block=40),
        output="both",
    )
    options.update(kwargs)
    return BeamformingConfig(**options)


# ----------------------------------------------------------------------------
# Plane waves
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["bartlett", "capon", "music"])
def test_vertical_plane_wave(method):
    """Every estimator finds velocity and azimuth of the wave."""
    result = beamform(_harmonic_array(), _config(method=method))
    best = result.maxima_at(rank=1)

    assert len(best) == len(result.maxima) > 0
    for item in best:
        assert item.velocity == pytest.approx(VELOCITY, abs=10.0)
        assert item.azimuth == pytest.approx(AZIMUTH)


def test_radial_and_transverse():
    """
    Rayleigh motion is radial, Love motion transverse.

    The polarized component peaks at the wave; the other one is
    weaker and peaks elsewhere (side lobes of the projected motion).
    """
    config = _config(components="three_component_separate",
                     method="bartlett", frequencies=[4.0])
    rayleigh = beamform(
        _harmonic_array(freqs=(4.0,), ellipticity=1.0), config
    )
    love = beamform(_harmonic_array(freqs=(4.0,), kind="love"), config)

    for result, polarized, other in (
        (rayleigh, "radial", "transverse"), (love, "transverse", "radial"),
    ):
        best = result.maxima_at(rank=1, component=polarized)[0]
        weak = result.maxima_at(rank=1, component=other)[0]
        assert best.azimuth == pytest.approx(AZIMUTH)
        assert best.velocity == pytest.approx(VELOCITY, abs=10.0)
        assert best.power > 2 * weak.power
        assert weak.azimuth != AZIMUTH


@pytest.mark.parametrize("ellipticity", [0.5, -1.5])
def test_rayleigh_joint_ellipticity(ellipticity):
    """
    RTBF recovers velocity, azimuth and signed ellipticity.

    The ellipticity grid has steps of 2.83 degrees: the angle must be
    within 3 degrees of atan(e).
    """
    array = _harmonic_array(freqs=(5.0,), ellipticity=ellipticity)
    result = beamform(array, _config(components="rayleigh_joint",
                                     frequencies=[5.0]))
    best = result.maxima_at(rank=1)[0]

    assert best.component == "rayleigh_joint"
    assert best.velocity == pytest.approx(VELOCITY, abs=10.0)
    assert best.azimuth == pytest.approx(AZIMUTH)
    assert best.ellipticity_angle == pytest.approx(
        math.degrees(math.atan(ellipticity)), abs=3.0
    )
    assert result.spectrum(5.0).ellipticity.shape == (
        result.spectrum(5.0).power.shape
    )


def test_power_is_the_density_of_the_wave():
    """
    Bartlett gives the PSD of a random plane wave; Capon slightly less.

    The wave is white noise of standard deviation 2 (one-sided PSD
    2 * 4 * delta = 0.08), delayed at each station; the grid contains
    the true wavevector. The Bartlett power at the wavevector equals
    the mean PSD measured at the stations (same windows) within 2%,
    the small loss being the incoherence of windows shifted by the
    delays (at most 0.12 s in 4 s windows). The measured PSD equals
    0.08 within its statistical error (about 2% over 16 frequencies):
    tolerance 8%. Capon estimated from L windows is biased low, by a
    factor of about (L - n + 1) / L for n stations (Capon and Goodman,
    1970): it lies between that fraction and the Bartlett power.
    """
    rng = np.random.default_rng(5)
    coordinates = _coordinates(aperture=50.0)
    source = 2.0 * rng.standard_normal(30000)
    motion = [
        {"Z": fourier.shift_time(source, DELTA, tau)}
        for tau in _delays(coordinates)
    ]
    array = _array(motion, coordinates, comps="Z")
    frequencies = np.linspace(3.0, 10.5, 16)
    grid = SearchGrid(200.0, 600.0, 5, scale="linear",
                      azimuths=[AZIMUTH, AZIMUTH + 90.0])
    _, values = coefficients(array.vertical(), DELTA, frequencies, 400,
                             overlap=0.5)
    measured = np.mean(np.abs(values) ** 2)
    power = {}

    for method in ("bartlett", "capon"):
        result = beamform(array, _config(
            method=method, frequencies=frequencies, grid=grid,
        ))
        power[method] = np.mean([
            spectrum.power[:, 0, 2] for spectrum in result.spectra
        ])

    assert measured == pytest.approx(0.08, rel=0.08)
    assert 0.98 < power["bartlett"] / measured <= 1.0
    assert (40 - 10 + 1) / 40 < power["capon"] / power["bartlett"] <= 1.0


def test_units_of_the_power():
    """Units follow the data; MUSIC and RTBF are flagged."""
    array = _harmonic_array(freqs=(4.0,))
    counts = ArrayData([
        StationData(station.sid, [
            Record(data=rec.data, delta=DELTA, time=START, sid=rec.sid,
                   units="counts")
            for rec in station.records
        ], station.coordinates)
        for station in array
    ])

    assert beamform(counts, _config()).power_units == "(counts)**2/Hz"
    assert beamform(array, _config()).power_units == "(unit)**2/Hz"
    assert beamform(array, _config(method="music")).power_units is None
    assert beamform(
        array, _config(components="rayleigh_joint", frequencies=[4.0])
    ).power_units == "((unit)**2/Hz)**2"


# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------


def test_rfft_nearest_equals_direct_at_bins():
    """With fixed windows the RFFT gives the direct result at bins."""
    array = _harmonic_array()
    direct = beamform(array, _config(frequencies=[4.0, 6.0]))
    rfft = beamform(array, _config(frequencies=[4.03, 5.98],
                                   transform="rfft"))

    np.testing.assert_allclose(rfft.frequencies, [4.0, 6.0])
    for first, second in zip(direct.spectra, rfft.spectra):
        np.testing.assert_allclose(second.power, first.power, rtol=1e-9)


def test_rfft_all_bins_in_range():
    """A range without size processes every FFT bin in it."""
    from shakelab.exploration.yarray.config import FrequencyConfig

    result = beamform(_harmonic_array(), _config(
        frequencies=FrequencyConfig(minimum=3.9, maximum=5.1),
        transform="rfft", output="maxima",
    ))

    np.testing.assert_allclose(result.frequencies, [4.0, 4.25, 4.5, 4.75,
                                                    5.0])
    assert result.spectra is None


def test_cycle_windows_and_block_times():
    """Windows in cycles; block times are the centres of the blocks."""
    array = _harmonic_array(nsamp=6000)
    config = _config(window=WindowConfig(length=8, windows_per_block=20,
                                         overlap=0.5, block_overlap=0.0))
    result = beamform(array, config)

    # 4 Hz: windows of 200 samples, step 100; block 0 covers samples
    # 0-2099, centre 1049.5 samples = 10.495 s.
    times = result.block_times[4.0]
    assert times[0] == array.starttime + 10.495
    assert len(times) == 2
    assert len(result.block_times[6.0]) == 4


def test_statistical_uncertainty():
    """Empirical standard errors are positive and in the spectra."""
    result = beamform(_harmonic_array(), _config(
        method="capon", statistical_uncertainty="empirical",
        components="radial", frequencies=[4.0],
    ))
    best = result.maxima_at(rank=1)[0]

    assert 0.0 < best.power_standard_error < best.power
    assert result.spectrum(4.0).power_standard_error.shape == (
        result.spectrum(4.0).power.shape
    )


def test_geometric_uncertainty():
    """
    Coordinate uncertainty lowers the peak and needs covariances.

    With zero covariance the result equals the plain Bartlett power.
    """
    array = _harmonic_array(freqs=(6.0,))
    config = _config(method="bartlett", frequencies=[6.0],
                     geometric_uncertainty="expected_coherence")

    with pytest.raises(ValueError, match="covariances"):
        beamform(array, config)

    def with_covariance(value):
        return ArrayData([
            StationData(st.sid, st.records, st.coordinates,
                        coordinate_covariance=value)
            for st in array
        ])

    plain = beamform(array, config.replace(geometric_uncertainty="none"))
    zero = beamform(with_covariance(0.0), config)
    blurred = beamform(with_covariance(25.0), config)

    np.testing.assert_allclose(zero.spectra[0].power, plain.spectra[0].power)
    assert np.max(blurred.spectra[0].power) < np.max(plain.spectra[0].power)


def test_output_modes():
    """maxima, spectra or both are kept."""
    array = _harmonic_array()

    assert beamform(array, _config(output="maxima")).spectra is None
    spectra = beamform(array, _config(output="spectra"))
    assert spectra.maxima == () and spectra.spectra

    several = beamform(array, _config(maxima=3))
    ranks = {item.rank for item in several.maxima}
    assert ranks <= {1, 2, 3} and 1 in ranks


def test_strict_realizations():
    """Capon with too few windows per block is an error."""
    config = _config(window=WindowConfig(type="fixed", length=4.0,
                                         windows_per_block=10))

    with pytest.raises(ValueError, match="realizations"):
        beamform(_harmonic_array(), config)

    beamform(_harmonic_array(), config.replace(strict=False))


def test_short_data_gives_empty_frequencies():
    """Frequencies without a complete block are warned and empty."""
    config = _config(window=WindowConfig(length=40, windows_per_block=40))

    with pytest.warns(BeamformingWarning, match="2 frequencies"):
        result = beamform(_harmonic_array(nsamp=3000), config)

    assert result.maxima == ()
    assert result.block_times == {4.0: (), 6.0: ()}
    assert result.spectra[0].power.shape == (0, 72, 61)


# ----------------------------------------------------------------------------
# Beamformer and results
# ----------------------------------------------------------------------------


def test_beamformer_cache_and_config(tmp_path):
    """The cache does not change the result; configs from files."""
    array = _harmonic_array()
    config = _config()
    reference = beamform(array, config)
    beamformer = Beamformer(config, cache=True)

    for _ in range(2):
        result = beamformer.process(array)
        np.testing.assert_array_equal(result.spectra[0].power,
                                      reference.spectra[0].power)

    assert set(beamformer._cache) == {4.0, 6.0}

    path = tmp_path / "config.json"
    config.replace(method="bartlett").write(path)
    beamformer.set_config(path)

    assert beamformer.config.method.value == "bartlett"
    assert beamformer._cache == {}

    with pytest.raises(TypeError):
        Beamformer(42)


def test_concatenate():
    """Results of consecutive data join with renumbered blocks."""
    array = _harmonic_array(nsamp=24000)
    config = _config()
    whole = beamform(array, config)
    first, second = (
        beamform(ArrayData([
            StationData(st.sid, [rec.extract(start, end)
                                 for rec in st.records], st.coordinates)
            for st in array
        ]), config)
        for start, end in ((None, 119.99), (120.0, None))
    )
    joined = BeamformingResult.concatenate([first, second])

    assert len(joined.block_times[4.0]) == (
        len(first.block_times[4.0]) + len(second.block_times[4.0])
    )
    blocks = [item.block for item in joined.maxima_at(frequency=4.0)]
    assert blocks == list(range(len(blocks)))
    assert joined.spectrum(4.0).power.shape[0] == len(blocks)
    assert len(whole.maxima) >= len(joined.maxima)

    with pytest.raises(ValueError):
        BeamformingResult.concatenate(
            [first, beamform(array, config.replace(method="bartlett"))]
        )


def test_result_access_and_summary(capsys):
    """Selections, nearest-frequency lookup and info."""
    result = beamform(_harmonic_array(), _config(maxima=2))

    assert result.spectrum(4.1).frequency == 4.0
    assert all(item.frequency == 6.0
               for item in result.maxima_at(frequency=5.9))
    assert result.components == ("vertical",)
    assert len(result.sids) == 10 and result.coordinates.shape == (10, 3)

    best = result.maxima_at(rank=1)[0]
    assert best.slowness == pytest.approx(1.0 / best.velocity)
    assert math.hypot(best.kx, best.ky) == pytest.approx(best.wavenumber)

    result.info()
    assert "median rank-1 velocity" in capsys.readouterr().out

    with pytest.raises(KeyError):
        result.spectrum(4.0, "radial")


# ----------------------------------------------------------------------------
# Comparison with the implementation of August 2026
# ----------------------------------------------------------------------------


def _sesame_array():
    """The nine SESAME stations (array file)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_array(DATA / "sesame_array.json")


@pytest.mark.parametrize(
    "case",
    ["vertical_capon", "radial_bartlett", "transverse_music",
     "three_component_capon", "rayleigh_joint"],
)
def test_reference_august_2026(case):
    """
    Maxima equal those of the YArray implementation of August 2026.

    The reference (``data/reference.json``) was computed by the
    original code of 12 August 2026 on the same data. Wavevectors,
    ranks, blocks, components and ellipticities are identical; the
    power differs by the change of normalization: the old coefficients
    were ``2 / N * sum w x e`` with the taper normalized to unit mean
    square, so that the new spectral density is ``N * delta / 2`` times
    the old power (squared for the RTBF objective; MUSIC is unchanged
    up to the same factor, inverse, being a reciprocal of a projection
    without the matrix scale).
    """
    reference = json.loads(
        (DATA / "reference.json").read_text()
    )["cases"][case]
    config = BeamformingConfig.from_dict(reference["config"])
    array = _sesame_array()
    result = beamform(array, config)

    assert len(result.maxima) == len(reference["maxima"])

    for new, old in zip(result.maxima, reference["maxima"]):
        assert new.frequency == old["frequency"]
        assert (new.block, new.rank) == (old["block"], old["rank"])
        assert new.component == old["component"]
        assert new.wavenumber == pytest.approx(old["wavenumber"], rel=1e-12)
        assert new.azimuth == old["azimuth"]

        if old["ellipticity"] is not None:
            assert new.ellipticity == pytest.approx(old["ellipticity"])

        size = config.window.size(array.delta, new.frequency)
        factor = size * array.delta / 2.0

        if config.method.value == "music":
            ratio = 1.0
        elif config.components.value == "rayleigh_joint":
            ratio = factor ** 2
        else:
            ratio = factor

        assert new.power == pytest.approx(old["power"] * ratio, rel=1e-8)

        if old["power_standard_error"] is not None:
            assert new.power_standard_error == pytest.approx(
                old["power_standard_error"] * ratio, rel=1e-6
            )
