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
Waveforms, channel metadata, spectra and data services.

The package holds the seismic recordings and what is needed to work
with them, in layers that depend only on those below:

- primitives: Fourier transforms (:mod:`.fourier`), source identifiers
  and channel patterns (:mod:`.codes`);
- the instrument response (:mod:`.response`) and the channel metadata
  (:mod:`.metadata`, :mod:`.header`);
- waveforms (:mod:`.base`: ``Record``, ``Stream``,
  ``StreamCollection``), spectra (:mod:`.spectrum`), time-frequency
  representations (:mod:`.spectrogram`), windows and cross-spectral
  matrices for multichannel methods (:mod:`.windowing`,
  :mod:`.crossspectral`), frequency-wavenumber transforms of DAS data
  (:mod:`.wavenumber`) and time-limited streams (:mod:`.circular`,
  experimental);
- file formats (:mod:`.libio`, :mod:`.stationxml`, :mod:`.respfile`,
  :mod:`.stationtext`) and the single entry point for files
  (:mod:`.io`);
- data services: FDSN web services (:mod:`.fdsnws`) and SeedLink
  (:mod:`.seedlink`).

The modules are imported explicitly (nothing is imported with the
package). See ``README.md`` in this folder.

.. code-block:: python

    from shakelab.signals.io import metadata_reader, reader

    waveforms = reader("data/*.mseed", metadata=metadata_reader("st.xml"))
    waveforms.deconvolve_response()
"""

__all__: list = []
