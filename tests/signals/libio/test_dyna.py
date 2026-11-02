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
Unit tests for :mod:`shakelab.signals.libio.dyna`.

Regression test for reading DYNA 1.2 files of ESM, whose header length
depends on the database version.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from shakelab.signals.libio.dyna import Dyna, dynaread

ESM = (
    Path(__file__).resolve().parents[1] / "data" / "ESM" / "ASCII"
    / "IV.MODE..HNZ.D.IT-2012-0008.ACC.MP.ASC"
)


def test_read_esm_file_keeps_all_samples():
    """
    All samples are read, starting from the first one.

    The ESM 3.0 header has 65 lines; the reader assumed 64 and lost the
    first sample, shifting the series by one sampling interval.
    """
    dyna = Dyna()
    dyna.read(str(ESM))

    assert dyna.header["NDATA"] == 40000
    assert len(dyna.data) == 40000
    assert dyna.data[:3] == [0.0, 0.000006, 0.000013]
    assert dyna.header["LOW_CUT_FREQUENCY_HZ"] == 0.03
    assert dyna.header["HIGH_CUT_FREQUENCY_HZ"] == 40.0


def test_dynaread_record():
    """The record carries identifier, start time and sampling."""
    record = dynaread(str(ESM))

    assert record.sid == "IV.MODE..HNZ"
    assert str(record.time).startswith("2012-05-20T02:03:23")
    assert record.delta == 0.005
    assert record.nsamp == 40000
    assert np.abs(record.data).max() == 32.733094
