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
Fixtures of the tests of :mod:`shakelab.signals.fdsnws`.
"""

from __future__ import annotations

import pytest

from fdsn_server import FakeFdsnServer
from shakelab.signals.fdsnws import transport


# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------


@pytest.fixture
def server():
    """Fake FDSN server, closed at the end of the test."""
    fake = FakeFdsnServer()
    yield fake
    fake.close()


@pytest.fixture(autouse=True)
def local_requests(monkeypatch):
    """Retries without waiting, local requests without proxy."""
    monkeypatch.setattr(transport, "RETRY_WAIT", 0.0)
    monkeypatch.setattr(transport, "RETRY_MAX_WAIT", 0.0)
    for name in ("NO_PROXY", "no_proxy"):
        monkeypatch.setenv(name, "127.0.0.1,localhost")
