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
ShakeLab: Python tools for engineering seismology.

The version of the installed package is available as
``shakelab.__version__``; it is defined once in ``pyproject.toml``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

name = "shakelab"

try:
    __version__ = version("shakelab")
except PackageNotFoundError:
    # Source tree used without installation.
    __version__ = "unknown"
