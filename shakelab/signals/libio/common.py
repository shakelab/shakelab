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
Definitions shared by the waveform file formats.

:class:`NotWaveformError` is raised by the conversion of a valid file
that holds no waveform (a response spectrum in DYNA, a spectrum or x-y
data in SAC), so that :func:`shakelab.signals.io.reader` can skip such
files when it reads directories and patterns.
"""

from __future__ import annotations

__all__ = [
    "NotWaveformError",
]


# ----------------------------------------------------------------------------
# Exceptions
# ----------------------------------------------------------------------------

class NotWaveformError(ValueError):
    """
    A valid file of a waveform format that holds no waveform.

    Subclass of :class:`ValueError`, raised when such a file is
    converted to a ShakeLab Record; the format class still reads it.
    """
