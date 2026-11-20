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
Waveform file formats.

One module per format: MiniSEED with libmseed (:mod:`.cymseed`) or in
pure Python (:mod:`.mseed`), SAC (:mod:`.sac`), DYNA 1.2
(:mod:`.dyna`), Simple ASCII (:mod:`.simpleascii`), text tables
(:mod:`.ascii`) and NI TDMS (:mod:`.tdms`). Every module follows the
same scheme: a reader ``<format>read``, a writer ``<format>write``, a
detection function ``is_<format>`` and a class with the content of a
file in the terms of the format.

The formats are usually read and written through
:func:`shakelab.signals.io.reader` and
:func:`shakelab.signals.io.writer`, which recognize them and import
each module only when needed. See ``README.md`` in this folder for
the common scheme and for adding a format.
"""

__all__: list = []
