# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
"""
import numpy as np
from scipy import interpolate
from copy import deepcopy

import shakelab.signals.base as base
from shakelab.signals import response

numeric_type = (int, float, complex,
                np.int8, np.int16, np.int32, np.int64,
                np.uint8, np.uint16, np.uint32, np.uint64)


class Spectrum():
    """
    A class for Discrete Fourier spectrum in RFFT format (one-sided positive).

    Attributes:
        head (base.Header): The ShakeLab signal's header object.
        data (array-like): The RFFT complex spectrum
        dfreq (float): The frequency spacing.

    Methods:
        __init__(self, record=None):
            Initialize a Spectrum object.
        __len__(self):
            Get the actual length of the spectrum.
        __getitem__(self, sliced):
            Get a slice of the spectrum.
        fft(self, data, delta=None, norm=False):
            Perform the FFT on the input data.
        ifft(self, norm=False):
            Perform the IFFT on the spectrum.
        from_record(self, record, norm=False):
            Initialize the Spectrum from a Record.
        to_record(self, norm=False):
            Create a Record object from the Spectrum.
        filter(self, highpass=None, lowpass=None):
            Apply highpass and lowpass filtering.
        resample(self, frequency):
            Extract spectrum samples at specific frequencies.
        logsmooth(self, sigma=0.2, memsafe=False):
            Logarithm smoothing of (complex) spectra.
    """

    def __init__(self, data=None, delta=None, nsamp=None):
        """
        Initialize a Spectrum object.

        Args:
            record (base.Record, optional): A record object to initialize the
                spectrum from.
        """
        self.head = base.Header()
        self.data = np.array([], dtype="complex")
        self.dfreq = None

        if isinstance(data, base.Record):
            self.from_record(data)

        elif isinstance(data, (list, np.ndarray)):
            if delta is None:
                raise ValueError('Delta must be provided')
            else:
                self._fft(data, delta)

    def __len__(self):
        return self.nfreq

    def __getitem__(self, sliced):
        return self.data[sliced]

    def __truediv__(self, value):
        """
        """
        spec_mod = self.copy()
        if isinstance(value, numeric_type):
            spec_mod.data = spec_mod.data / value
        elif isinstance(value, Spectrum):
            spec_mod.data /= value.data
        else:
            raise TypeError('unsupported operand type(s) for /')

        return spec_mod

    @property
    def nfreq(self):
        """
        Return the actual length of the spectrum.
        """
        return len(self.data)

    @property
    def nsamp(self):
        """
        Return the length of the corresponding time series.
        """
        #return 2*(self.nfreq-1) + self._parity
        return _nsamp(self.head.delta, self.dfreq)

    @property
    def amplitude(self):
        """
        Property to compute the amplitude (modulus) of the spectrum.
        """
        return np.abs(self.data)

    @property
    def phase(self):
        """
        Property to compute the phase of the spectrum.
        """
        return np.angle(self.data)

    @property
    def phase_unwrap(self):
        """
        Property to compute the unwrapped phase of the spectrum.
        """
        return np.unwrap(self.phase)

    @property
    def frequency_axis(self):
        """
        Property to compute the frequency axis of the spectrum.
        """
        return frequency_axis(self.head.delta, self.nsamp)

    def set_data(self, data, delta, nsamp=None):
        """
        """
        if nsamp is None:
            nsamp = (len(data)-1) * 2

        self.data = data
        self.head.delta = delta
        self.dfreq = _dfreq(delta, nsamp)

    def fft(self, data, delta=None, nsamp=None, norm=False):
        """
        Compute the FFT spectrum from real value input data.

        Args:
            data (array-like): The input time-series data (must be real)
            delta (float, optional): The time step. If provided, it updates
                the internal time step.
            norm (bool, optional): If True, normalizes the result.

        Returns:
            None
        """
        if nsamp is None:
            nsamp = len(data)

        if delta is not None:
            self.head.delta = delta

        self.dfreq = _dfreq(self.head.delta, nsamp)

        self.data = _fft(data, nsamp)

        if norm:
            # to check
            self.data =  self.data / self.head.delta / self.nsamp

    def ifft(self, nsamp=None, norm=False):
        """
        Compute the IFFT of the spectrum.

        Args:
            norm (bool, optional): If True, normalizes the result.

        Returns:
            array-like: The inverse FFT result.
        """
        if nsamp is None:
            nsamp = self.nsamp

        data = _ifft(self.data, nsamp)

        if norm:
            # to check
            data = data / self.dfreq * self.nfreq

        return data

    def invert(self, waterlevel=100, method='smooth'):
        """
        """
        self.data = response.inverse_spectrum(self.data,
                                              waterlevel=waterlevel,
                                              method=method)

    def from_record(self, record, norm=False):
        """
        Initialize the spectrum from a Record.

        Args:
            record (base.Record): The Record object to initialize from.
            norm (bool, optional): If True, normalizes the spectrum.
        """
        if isinstance(record, base.Record):
            self.head = deepcopy(record.head)
            self.fft(record.data, norm=norm)
        else:
            pass

    def to_record(self, norm=False):
        """
        Create a Record object from the Spectrum.

        Args:
            norm (bool, optional): If True, normalizes the data in the Record.

        Returns:
            base.Record: The created Record object.
        """
        record = base.Record()
        record.head = deepcopy(self.head)
        record.data = self.ifft(norm=norm)

        return record

    def copy(self):
        """
        """
        return deepcopy(self)

    def filter(self, highpass=None, lowpass=None, order=4, filt_type='bw'):
        """
        In-place highpass and lowpass filtering to the spectrum.

        Args:
            highpass (float, optional): The highpass frequency cutoff.
            lowpass (float, optional): The lowpass frequency cutoff.
        """
        freq = self.frequency_axis

        if filt_type == 'sharp':
            if (highpass is not None):
                self.data[freq < highpass] = 0.

            if (lowpass is not None):
                self.data[freq > lowpass] = 0.

        elif filt_type == 'bw':
            if (highpass is not None):
                self.data *= butterworth(freq, highpass, order=-order)

            if (lowpass is not None):
                self.data *= butterworth(freq, lowpass, order=order)

        else:
            assert TypeError('Not implemented filter type')

    def resample(self, frequency):
        """
        Extract spectrum samples at specific frequencies.
        """
        f = interpolate.interp1d(self.frequency, self.data)
        return f(frequency)

    def logsmooth(self, sigma=0.2, memsafe=False):
        """
        Logarithm smoothing of (complex) spectra.

        Args:
            sigma (float, optional): The smoothing parameter.
            memsafe (bool, optional): If True, a memory-saving algorithm is
            used.

        Returns:
            array-like: The smoothed spectrum.

        Notes:
            - 0-frequency is preserved
            - algorithm is slow and memory consuming, must be optimized.
        """

        slen = len(self) - 1
        freq = np.log(self.frequency[1:]) / (np.sqrt(2) * sigma)
        data = np.log(self.data[1:])
        s0 = self.data[0]

        if not memsafe:
            # Fast vectorial version, although memory consuming

            # Gaussian weighting window
            g = np.exp(-(np.tile(freq, slen) - np.repeat(freq, slen))**2)
            g = g.reshape(slen, slen)
            
            sumg = np.matmul(g, np.ones(slen))
            smat = g / sumg

            return np.insert(np.exp(np.matmul(data, smat)), 0, s0)

        else:
            # Slower version, but memory saving for large time series.
            sdata = np.zeros(slen, dtype='complex')

            for i,f0 in enumerate(freq):
                # Gaussian weighting window
                g = np.exp(-(freq - f0)**2)

                sdata[i]= np.matmul(data, g/sum(g))

            return np.insert(np.exp(sdata), 0, s0)


def rfft_length(nsamp):
    """
    Calculate the length of the Real Fast Fourier Transform (RFFT) positive
    frequency axis, accounting for odd/even numbers of samples.

    Parameters:
    - nsamp (int): The number of samples in the input signal.

    Returns:
    - int: The length of the positive frequency axis in the RFFT.

    Example:
    >>> rfft_length(10)
    6
    """
    if nsamp % 2 == 1:
        return (nsamp + 1) // 2
    else:
        return nsamp // 2 + 1

def irfft_length(delta, dfreq):
    """
    Calculate the length of the time series corresponding to a rfft
    spectrum of given time spacing (delta) and frequency spacing (dfreq).

    Parameters:
    - delta (float): Time spacing between samples.
    - dfreq (float): Frequency spacing between samples.

    Returns:
    - int: The number of samples needed to satisfy the specified spacing.

    Example:
    >>> irfft_length(0.01, 10)
    100
    """
    return round(1 / (dfreq * delta))

def _fnum(nsamp):
    """
    Wrapper function to compute the length of the expected RFFT spectrum.
    """
    return rfft_length(nsamp)

def _nsamp(delta, dfreq):
    """
    Wrapper function to comnpute the length of the time series from a RFFT spectrum.
    """
    return irfft_length(delta, dfreq)

def _dfreq(delta, nsamp):
    """
    """
    return 1 / (delta * nsamp)

def _delta(dfreq, nsamp):
    """
    """
    return 1 / (dfreq * nsamp)

def _fft(data, nsamp=None):
    """
    Wrapper to numpy's rfft'.
    """
    return np.fft.rfft(data, nsamp)

def _ifft(data, nsamp=None):
    """
    Wrapper to numpy's irfft'.
    """
    return np.fft.irfft(data, nsamp)

def frequency_axis(delta, nsamp):
    """
    """
    return np.fft.rfftfreq(nsamp, delta)

def frequency_range(fmin, fmax, fnum, log=True):
    """
    Compute a linear or logarithmic frequency range

    :param float fmin:
        Minimum frequency

    :param float fmax:
        Maximum frequency

    :param int fnum:
        Number of frequencies

    :param boolean log:
        Select between linear or logarithmic spacing
        (default is logarithmic)

    :return numpy.array freq:
        The frequency axis
    """

    if log:
        freq = np.logspace(np.log10(fmin), np.log10(fmax), fnum)
    else:
        freq = np.linspace(fmin, fmax, fnum)

    return freq

def _to_complex(amplitude, phase):
    """
    """
    return amplitude * np.exp(1j * phase)

def shift_time(signal, delta, shift):
    """
    Shift a signal in time by using fft-based circular convolution.
    No zero-padding is assumed.
    """
    nsamp = len(signal)
    freq = frequency_axis(delta, len(signal))
    expt = np.exp(-2*1j*np.pi*shift*freq)

    return _ifft(_fft(signal, nsamp)*expt, nsamp)

def butterworth(freq, corner_freq, order=4, minimum_phase=False):
    """
    Note: if order is negative, the filter is highpass
    """
    hf = np.zeros(len(freq), dtype=np.complex_)

    for i, f in enumerate(freq):
        if f > 0 or order > 0:
            hf_mag = 1/(1+(1j*f/corner_freq)**(2*order))

            if minimum_phase:
                # Convert to minimum phase by taking
                # the negative square root
                hf_phase = -np.sqrt(1-np.abs(hf_mag)**2)
                hf[i] = np.exp(1j*hf_phase) * np.abs(hf_mag)
            else:
                hf[i] = hf_mag

    return hf
