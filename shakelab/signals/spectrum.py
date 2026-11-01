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
Frequency-domain representation of waveforms.

This module provides :class:`Spectrum`, the frequency-domain
counterpart of :class:`~shakelab.signals.base.Record`: a one-sided
complex spectrum together with the metadata of the corresponding time
series.

Design choices
--------------
- The constructor receives spectral values (e.g. a transfer function);
  spectra of time series are built with :meth:`Spectrum.from_record`
  or :meth:`Spectrum.from_samples`. ``Record.to_spectrum`` and
  ``Record.from_spectrum`` delegate to this module, which holds the
  canonical conversion between the two domains.
- The number of samples of the time series is stored explicitly: it is
  needed to distinguish even and odd lengths in the inverse transform.

Conventions
-----------
- One-sided spectra scaled as ``delta * DFT`` (see
  :mod:`shakelab.signals.fourier`); frequencies in Hz.
- The header keeps the metadata of the time series (sampling interval,
  start time, identifiers, units).
- The ``info()`` method prints a summary for interactive inspection.
"""

from __future__ import annotations

from numbers import Number
from typing import Any, Optional

import numpy as np

from shakelab.signals import fourier
from shakelab.signals.base import Record
from shakelab.signals.header import Header

__all__ = ["Spectrum"]


# ----------------------------------------------------------------------------
# Spectrum
# ----------------------------------------------------------------------------


class Spectrum:
    """
    One-sided complex spectrum of a real time series.

    Parameters
    ----------
    data
        Complex spectral values at the frequencies
        ``k / (nsamp * delta)``, ``k = 0, ..., nsamp // 2``.
    delta
        Sampling interval of the time series, in seconds.
    nsamp
        Number of samples of the time series. Default is
        ``2 * (len(data) - 1)`` (even length).
    head
        Header with the metadata of the time series; a copy is stored.
        If given, ``delta`` defaults to ``head.delta``.

    Attributes
    ----------
    head : Header
        Metadata of the time series.
    data : numpy.ndarray
        Complex spectral values.

    Raises
    ------
    ValueError
        If ``nsamp`` does not match the number of spectral values.
    """

    def __init__(
        self,
        data: Any = None,
        delta: Optional[float] = None,
        nsamp: Optional[int] = None,
        head: Optional[Header] = None,
    ) -> None:
        self.head = Header() if head is None else head.copy()

        if delta is not None:
            self.head.delta = delta

        self.data = np.array([], dtype=complex)
        self._nsamp = 0

        if data is not None:
            self.set(data, nsamp)

    # ------------------------------------------------------------------------
    # Construction and conversion
    # ------------------------------------------------------------------------

    def set(self, data: Any, nsamp: Optional[int] = None) -> None:
        """
        Set the spectral values.

        Parameters
        ----------
        data
            Complex spectral values.
        nsamp
            Number of samples of the time series; default
            ``2 * (len(data) - 1)``.

        Raises
        ------
        ValueError
            If ``data`` is not one-dimensional or ``nsamp`` does not
            match its length.
        """
        data = np.array(data, dtype=complex)

        if data.ndim == 0:
            data = data.reshape(1)

        if data.ndim != 1:
            raise ValueError("spectrum data must be one-dimensional")

        if nsamp is None:
            nsamp = max(2 * (data.size - 1), 1)

        nsamp = int(nsamp)

        if fourier.rfft_length(nsamp) != data.size:
            raise ValueError(
                f"{data.size} spectral values do not match a time "
                f"series of {nsamp} samples"
            )

        self.data = data
        self._nsamp = nsamp

    @classmethod
    def from_samples(
        cls,
        data: Any,
        delta: float,
        nsamp: Optional[int] = None,
    ) -> Spectrum:
        """
        Return the spectrum of a real time series.

        Parameters
        ----------
        data
            Real samples.
        delta
            Sampling interval in seconds.
        nsamp
            Transform length (zero padding or truncation). Default is
            the number of samples.

        Returns
        -------
        Spectrum
            Spectrum scaled as ``delta * DFT``.
        """
        data = np.asarray(data)

        if nsamp is None:
            nsamp = data.size

        spectrum = cls(delta=delta)
        spectrum.set(fourier.fft(data, delta, nsamp), nsamp)

        return spectrum

    @classmethod
    def from_record(cls, record: Record) -> Spectrum:
        """
        Return the spectrum of a Record, with a copy of its metadata.

        Parameters
        ----------
        record
            Record with a sampling interval.

        Returns
        -------
        Spectrum
            Spectrum of the record samples.

        Raises
        ------
        TypeError
            If ``record`` is not a Record.
        ValueError
            If the record has no sampling interval.
        """
        if not isinstance(record, Record):
            raise TypeError("record must be a Record instance")

        if record.delta is None:
            raise ValueError(
                "sampling interval is required to compute a spectrum"
            )

        spectrum = cls.from_samples(record.data, record.delta)
        spectrum.head = record.head.copy()

        return spectrum

    def to_samples(self) -> np.ndarray:
        """
        Return the real time series of the spectrum.

        Returns
        -------
        numpy.ndarray
            ``nsamp`` real samples (inverse of :meth:`from_samples`).
        """
        self._require_delta()
        return fourier.ifft(self.data, self.delta, self.nsamp)

    def to_record(self) -> Record:
        """
        Return the Record of the spectrum, with a copy of its metadata.

        Returns
        -------
        Record
            Time series with the header of the spectrum.
        """
        record = Record()
        record.head = self.head.copy(parent=record)
        record.data = self.to_samples()

        return record

    def copy(self) -> Spectrum:
        """Return an independent copy of the spectrum."""
        spectrum = Spectrum(head=self.head)
        spectrum.data = self.data.copy()
        spectrum._nsamp = self._nsamp

        return spectrum

    # ------------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------------

    def __len__(self) -> int:
        return self.nfreq

    def __getitem__(self, item: Any) -> Any:
        return self.data[item]

    def __repr__(self) -> str:
        return (
            f"Spectrum(sid={self.head.sid!r}, nfreq={self.nfreq}, "
            f"dfreq={self.dfreq!r})"
        )

    def __str__(self) -> str:
        return (
            f"source identifier {self.head.sid}, {self.nfreq} "
            f"frequencies, frequency step {self.dfreq} Hz, time series "
            f"of {self.nsamp} samples at {self.delta} s"
        )

    def info(self) -> None:
        """Print a human-readable description of the spectrum."""
        print(self)

    @property
    def delta(self) -> Optional[float]:
        """Sampling interval of the time series, in seconds."""
        return self.head.delta

    @property
    def nsamp(self) -> int:
        """Number of samples of the time series."""
        return self._nsamp

    @property
    def nfreq(self) -> int:
        """Number of frequencies."""
        return self.data.size

    @property
    def dfreq(self) -> Optional[float]:
        """Frequency step in Hz."""
        if self.delta is None or self.nsamp == 0:
            return None
        return 1.0 / (self.nsamp * self.delta)

    @property
    def frequency(self) -> np.ndarray:
        """Frequencies of the spectral values, in Hz."""
        self._require_delta()
        return fourier.frequency_axis(self.delta, self.nsamp)

    @property
    def amplitude(self) -> np.ndarray:
        """Amplitude (modulus) of the spectrum."""
        return np.abs(self.data)

    @property
    def phase(self) -> np.ndarray:
        """Phase of the spectrum, in radians, in (-pi, pi]."""
        return np.angle(self.data)

    @property
    def phase_unwrap(self) -> np.ndarray:
        """Phase of the spectrum unwrapped along frequency, in radians."""
        return np.unwrap(self.phase)

    # ------------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------------

    def _operand(self, value: Any, name: str) -> Any:
        """Return the array or scalar to combine with the data."""
        if isinstance(value, Number) and not isinstance(value, bool):
            return value

        if isinstance(value, Spectrum):
            if value.nsamp != self.nsamp or not np.isclose(
                value.delta, self.delta, rtol=1e-9, atol=0.0
            ):
                raise ValueError(
                    "spectra have different frequency axes"
                )
            return value.data

        raise TypeError(
            f"unsupported operand type for Spectrum {name}: "
            f"{type(value).__name__}"
        )

    def __mul__(self, value: Any) -> Spectrum:
        result = self.copy()
        result.data = self.data * self._operand(value, "multiplication")
        return result

    __rmul__ = __mul__

    def __truediv__(self, value: Any) -> Spectrum:
        operand = self._operand(value, "division")

        if isinstance(operand, Number) and operand == 0:
            raise ZeroDivisionError("cannot divide a Spectrum by zero")

        result = self.copy()
        result.data = self.data / operand
        return result

    # ------------------------------------------------------------------------
    # Processing
    # ------------------------------------------------------------------------

    def filter(
        self,
        highpass: Optional[float] = None,
        lowpass: Optional[float] = None,
        order: int = 4,
        minphase: bool = False,
        method: str = "butterworth",
    ) -> None:
        """
        Filter the spectrum in place.

        Parameters
        ----------
        highpass, lowpass
            Corner frequencies in Hz; at least one is required.
        order
            Butterworth order.
        minphase
            Use the causal (minimum-phase) Butterworth response instead
            of the zero-phase one (see
            :func:`shakelab.signals.fourier.butterworth`).
        method
            ``"butterworth"`` (default) or ``"sharp"``, which sets to
            zero the values outside the pass band.

        Raises
        ------
        ValueError
            If no corner is given or the method is unknown.
        """
        if highpass is None and lowpass is None:
            raise ValueError("at least one corner frequency is required")

        freq = self.frequency

        if method == "sharp":
            if highpass is not None:
                self.data[freq < highpass] = 0.0
            if lowpass is not None:
                self.data[freq > lowpass] = 0.0
            return

        if method != "butterworth":
            raise ValueError("method must be 'butterworth' or 'sharp'")

        if highpass is not None:
            self.data = self.data * fourier.butterworth(
                freq, highpass, order, "highpass", minphase
            )

        if lowpass is not None:
            self.data = self.data * fourier.butterworth(
                freq, lowpass, order, "lowpass", minphase
            )

    def interpolate(self, frequency: Any) -> np.ndarray:
        """
        Return the spectrum at given frequencies by linear interpolation.

        Real and imaginary parts are interpolated separately.

        Parameters
        ----------
        frequency
            Frequencies in Hz, within the frequency axis.

        Returns
        -------
        numpy.ndarray
            Complex values at ``frequency``.

        Raises
        ------
        ValueError
            If a frequency is outside the frequency axis.
        """
        freq = self.frequency
        frequency = np.asarray(frequency, dtype=float)

        if np.any(frequency < freq[0]) or np.any(frequency > freq[-1]):
            raise ValueError("frequency outside the spectrum range")

        return (
            np.interp(frequency, freq, self.data.real)
            + 1j * np.interp(frequency, freq, self.data.imag)
        )

    def invert(self, waterlevel: float = 100.0, method: str = "smooth"):
        """
        Replace the spectrum with its regularized inverse, in place.

        See :func:`shakelab.signals.fourier.inverse_spectrum`.
        """
        self.data = fourier.inverse_spectrum(
            self.data, waterlevel=waterlevel, method=method
        )

    def logsmooth(
        self,
        sigma: float = 0.2,
        truncation: float = 3.0,
        method: str = "gaussian",
        bandwidth: float = 40.0,
        phase: bool = True,
    ) -> None:
        """
        Smooth the spectrum in logarithmic frequency, in place.

        Parameters
        ----------
        sigma, truncation
            Width and truncation of the Gaussian kernel (natural-log
            frequency units).
        method
            ``"gaussian"`` (default) or ``"konno-ohmachi"``.
        bandwidth
            Bandwidth coefficient of the Konno-Ohmachi window.
        phase
            If True (default), smooth log-amplitude and unwrapped phase
            (:func:`~shakelab.signals.fourier.smooth_amplitude_phase`);
            otherwise smooth only the log-amplitude and keep the phase
            (:func:`~shakelab.signals.fourier.smooth_amplitude`).

        Notes
        -----
        The zero-frequency value is not smoothed.
        """
        options = dict(
            method=method,
            sigma=sigma,
            truncation=truncation,
            bandwidth=bandwidth,
        )

        if phase:
            self.data = fourier.smooth_amplitude_phase(
                self.frequency, self.data, **options
            )
        else:
            self.data = fourier.smooth_amplitude(
                self.frequency, self.data, **options
            )

    def _require_delta(self) -> None:
        if self.delta is None:
            raise ValueError("the sampling interval is not defined")
