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
Instrument response model.

This module represents the response of a recording instrument: stages
(gain, poles and zeros, digital filters) grouped into stage sets valid
in a time interval (epochs), stream responses holding the epochs of one
channel, and collections of stream responses indexed by source
identifier.

Design choices
--------------
- The response is a model of the instrument: each stage provides its
  transfer function ``H(f)`` (:meth:`StageResponse.transfer_function`)
  and a stage set the product over its active stages. The model does
  not depend on waveform classes: responses are applied to data by
  :meth:`shakelab.signals.spectrum.Spectrum.deconvolve_response` and
  the corresponding methods of records, streams and collections.
- :func:`select_response` chooses the stage or stage set applicable to
  a source identifier and a time.
- Only gain and poles-and-zeros stages are applied (``ACTIVE_STAGES``).
  Digital stages (FIR, coefficients, polynomial) are read and stored
  but not applied: FIR filters act at the sampling rate of their own
  stage, and StationXML polynomials describe non-linear sensors rather
  than frequency responses.

Conventions
-----------
- Frequencies in Hz; poles and zeros in rad/s (Laplace transform,
  ``s = 2j * pi * f``).
- The input units of a response are those of the ground motion (e.g.
  ``m/s**2``), the output units those of the recorded data.
- The ``info()`` methods print a summary for interactive inspection.
"""

from __future__ import annotations

import json as _json
import os as _os
from copy import deepcopy
from typing import Any, Optional

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import (
    cast_value,
    deserialize_complex,
    serialize_ndarray,
)

__all__ = [
    "ACTIVE_STAGES",
    "ResponseCollection",
    "StreamResponse",
    "StageSet",
    "StageResponse",
    "StageGain",
    "StagePoleZero",
    "StagePolynomial",
    "StageFIR",
    "select_response",
    "paz_transfer_function",
    "polynomial_transfer_function",
    "fir_transfer_function",
]

ACTIVE_STAGES = {
    "gain": True,
    "paz": True,
    "coefficients": False,
    "polynomial": False,
    "fir": False,
}
"""Stage types applied when a response is applied to data."""


# ----------------------------------------------------------------------------
# Collections and epochs
# ----------------------------------------------------------------------------


class ResponseCollection:
    """
    Collection of stream responses indexed by source identifier.

    Parameters
    ----------
    byte_stream
        Optional file to read (see :meth:`read`).
    ftype
        File type (``"xml"`` or ``"json"``); inferred from the file
        extension if None.
    """

    def __init__(self, byte_stream: Any = None, ftype: Optional[str] = None):
        self.response: list[StreamResponse] = []
        if byte_stream is not None:
            self.read(byte_stream, ftype=ftype)

    def __len__(self) -> int:
        return len(self.response)

    def __iter__(self):
        return iter(self.response)

    def __getitem__(self, key: Any) -> StreamResponse:
        if isinstance(key, str):
            return self.response[self._idx(key)]
        return self.response[key]

    def __repr__(self) -> str:
        return f"ResponseCollection(responses={len(self.response)})"

    def _idx(self, sid: str) -> int:
        """
        Return the position of a source identifier.

        Raises
        ------
        KeyError
            If the identifier is not in the collection.
        """
        sid_list = self.sid
        if sid not in sid_list:
            raise KeyError(f"source identifier not found: {sid}")
        return sid_list.index(sid)

    @property
    def sid(self) -> list:
        """Source identifiers of the stream responses."""
        return [resp.sid for resp in self.response]

    def append(self, obj: Any) -> None:
        """
        Append a stream response.

        If a stream response with the same identifier is already in the
        collection, the stage sets of ``obj`` are added to it.

        Raises
        ------
        ValueError
            If the object is not a StreamResponse.
        """
        if not isinstance(obj, StreamResponse):
            raise ValueError("Unsupported type to append")
        if obj.sid in self.sid:
            self[obj.sid].append(obj.stageset)
        else:
            self.response.append(obj)

    def extract(
        self,
        sid: str,
        time: Any = None,
        stage_number: Optional[int] = None,
    ) -> Any:
        """
        Return the stage set (or one of its stages) of a channel.

        Parameters
        ----------
        sid
            Source identifier.
        time
            Time of validity (default: now).
        stage_number
            If given, return that stage of the stage set.

        Raises
        ------
        KeyError
            If the channel or the epoch is not found.
        """
        result = self[sid].extract(time)
        if stage_number is not None:
            return result[stage_number]
        return result

    def read(
        self,
        byte_stream: Any,
        ftype: Optional[str] = None,
        append: bool = False,
    ) -> None:
        """
        Read responses from a StationXML or JSON file.

        Parameters
        ----------
        byte_stream
            File path (or, for StationXML, an open stream).
        ftype
            ``"xml"`` (also ``"sxml"``) or ``"json"``; inferred from the
            file extension if None.
        append
            If True, add to the current content; otherwise replace it.

        Raises
        ------
        ValueError
            If the file type is unknown or a stage type is unsupported.
        """
        if ftype is None:
            _, ext = _os.path.splitext(str(byte_stream))
            ftype = ext.lstrip(".").lower()

        if ftype in ("xml", "sxml"):
            # Local import: the StationXML parser builds response
            # objects, so it sits above this module.
            from shakelab.signals.stationxml import parse_sxml

            new = parse_sxml(byte_stream)
            if append:
                for resp in new.response:
                    self.append(resp)
            else:
                self.response = new.response

        elif ftype == "json":
            with open(byte_stream, "r") as fobj:
                content = _json.load(fobj)

            if not append:
                self.response = []

            for entry in content:
                sid = entry["sid"]
                stage_set = StageSet(entry["starttime"], entry["endtime"])

                for item in entry["stages"]:
                    data = deserialize_complex(item["data"])
                    data["stage_number"] = item["stage_number"]
                    stage_set.append(_stage_from_type(item["stage_type"],
                                                      data))

                if sid in self.sid:
                    self[sid].append(stage_set)
                else:
                    stream_response = StreamResponse(sid)
                    stream_response.append(stage_set)
                    self.append(stream_response)

        else:
            raise ValueError(f"Unknown file type: {ftype}")

    def write(self, file: Any, ftype: str = "json") -> None:
        """
        Write the responses to a JSON file.

        Raises
        ------
        NotImplementedError
            For file types other than JSON.
        """
        if ftype != "json":
            raise NotImplementedError("Only the json format is supported")

        out = []

        for resp in self.response:
            for stage_set in resp.stageset:
                stages = [
                    {
                        "stage_number": stage.stage_number,
                        "stage_type": stage.stage_type,
                        "data": serialize_ndarray(stage.to_dict()),
                    }
                    for stage in stage_set.stage
                ]
                out.append({
                    "sid": resp.sid,
                    "starttime": (
                        str(stage_set.starttime)
                        if stage_set.starttime else None
                    ),
                    "endtime": (
                        str(stage_set.endtime)
                        if stage_set.endtime else None
                    ),
                    "stages": stages,
                })

        with open(file, "w") as fobj:
            _json.dump(out, fobj, indent=2)

    def info(self) -> None:
        """Print a summary of all stream responses."""
        for resp in self.response:
            resp.info()


class StreamResponse:
    """
    Stage sets (epochs) of one channel.

    Parameters
    ----------
    sid
        Source identifier of the channel.
    """

    def __init__(self, sid: Optional[str]):
        self.sid = sid
        self.stageset: list[StageSet] = []

    def __len__(self) -> int:
        return len(self.stageset)

    def __getitem__(self, time: Any) -> StageSet:
        return self.extract(time)

    def __repr__(self) -> str:
        return f"StreamResponse(sid={self.sid!r}, epochs={len(self)})"

    def append(self, stage_set: Any) -> None:
        """
        Append one or more stage sets.

        Raises
        ------
        ValueError
            If an item is not a StageSet.
        """
        if not isinstance(stage_set, list):
            stage_set = [stage_set]
        for item in stage_set:
            if not isinstance(item, StageSet):
                raise ValueError("Not a valid stage set")
            self.stageset.append(item)

    def extract(self, time: Any = None) -> StageSet:
        """
        Return the stage set valid at a time.

        Parameters
        ----------
        time
            Any value accepted by :class:`Date`, or ``"now"``/None.

        Raises
        ------
        KeyError
            If no stage set is valid at that time.
        """
        for stage_set in self.stageset:
            if stage_set.match(time):
                return stage_set
        raise KeyError(
            f"no response epoch of {self.sid} is valid at {time}"
        )

    def info(self) -> None:
        """Print a summary of the stage sets."""
        print(f"Stream ID: {self.sid}")
        for stage_set in self.stageset:
            stage_set.info()


class StageSet:
    """
    Stages of a response, valid between two times.

    Parameters
    ----------
    starttime, endtime
        Validity interval (``endtime`` excluded); None means unbounded.
    """

    def __init__(self, starttime: Any = None, endtime: Any = None):
        if starttime is not None and not isinstance(starttime, Date):
            starttime = Date(starttime)
        if endtime is not None and not isinstance(endtime, Date):
            endtime = Date(endtime)

        self.starttime = starttime
        self.endtime = endtime
        self.stage: list[StageResponse] = []

    def __getitem__(self, idx: int) -> StageResponse:
        """Return the stage with the given stage number."""
        if isinstance(idx, int):
            return self.get_stage(idx)
        raise NotImplementedError("Not a valid stage index")

    def __repr__(self) -> str:
        return (
            f"StageSet(starttime={self.starttime}, endtime={self.endtime}, "
            f"stages={len(self.stage)})"
        )

    def info(self) -> None:
        """Print the validity interval and the list of stages."""
        start_str = str(self.starttime) if self.starttime else "unknown"
        end_str = str(self.endtime) if self.endtime else "now"
        print(f"Response valid from: {start_str} to: {end_str}")

        if not self.stage:
            print("Stages: none")
            return

        stage_str = "; ".join(
            f"({s.stage_number}) {s.stage_type}" for s in self.stage
        )
        print(f"Stages: {stage_str}")

    def append(self, stage_response: Any) -> None:
        """
        Append one or more stages.

        Raises
        ------
        ValueError
            If an item is not a StageResponse.
        """
        if not isinstance(stage_response, list):
            stage_response = [stage_response]

        for stage in stage_response:
            if not isinstance(stage, StageResponse):
                raise ValueError("Not a valid stage response")
            self.stage.append(stage)

    def match(self, time: Any) -> bool:
        """Return whether a time falls within the validity interval."""
        if not isinstance(time, Date):
            time = Date("now") if time in ("now", None) else Date(time)

        return (
            (self.starttime is None or time >= self.starttime)
            and (self.endtime is None or time < self.endtime)
        )

    def get_stage(self, stage_number: int) -> StageResponse:
        """
        Return the stage with the given stage number.

        Raises
        ------
        IndexError
            If the stage is not found.
        """
        for stage in self.stage:
            if stage.stage_number == stage_number:
                return stage
        raise IndexError(f"Stage {stage_number} not found")

    @property
    def active_stages(self) -> list:
        """Stages applied to data (see ``ACTIVE_STAGES``)."""
        return [s for s in self.stage if ACTIVE_STAGES[s.stage_type]]

    @property
    def input_units(self) -> Optional[str]:
        """
        Units of the ground motion: input units of the first stage.

        None if the first stage does not declare them.
        """
        if not self.stage:
            return None
        return getattr(self.stage[0], "input_units", None) or None

    @property
    def output_units(self) -> Optional[str]:
        """
        Units of the recorded data: output units of the last stage.

        None if the last stage does not declare them (e.g. a gain-only
        stage), so that units are never guessed from inner stages.
        """
        if not self.stage:
            return None
        return getattr(self.stage[-1], "output_units", None) or None

    def transfer_function(
        self,
        frequency: Any,
        delta: Optional[float] = None,
    ) -> np.ndarray:
        """
        Return the transfer function of the active stages.

        Parameters
        ----------
        frequency
            Frequencies in Hz.
        delta
            Sampling interval of the data (needed only by digital
            stages, which are not active by default).

        Returns
        -------
        numpy.ndarray
            Product of the transfer functions of the active stages.
        """
        result = np.ones(np.shape(frequency), dtype=complex)
        for stage in self.active_stages:
            result = result * stage.transfer_function(frequency, delta)
        return result

    def copy(self) -> StageSet:
        """Return a deep copy of the stage set."""
        return deepcopy(self)


# ----------------------------------------------------------------------------
# Stages
# ----------------------------------------------------------------------------


class StageResponse:
    """
    Base class of response stages.

    Parameters are defined in ``_KEYMAP`` as
    ``{"name": (type, default)}`` and can be given as a dictionary, as
    keyword arguments, or read from a JSON file.
    """

    _KEYMAP: dict = {}
    stage_type: Optional[str] = None

    def __init__(self, data: Optional[dict] = None, **kwargs):
        for key in self._KEYMAP:
            setattr(self, key, self._KEYMAP[key][1])

        if data is not None:
            self.from_dict(data)

        for key in kwargs:
            self[key] = kwargs[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a parameter, converting it to the type of ``_KEYMAP``."""
        if key in self._KEYMAP:
            dtype = self._KEYMAP[key][0]
            setattr(self, key, cast_value(value, dtype))
        else:
            raise KeyError(f"Invalid key: {key}")

    def __getitem__(self, key: str) -> Any:
        """Get a parameter by key."""
        if key in self._KEYMAP:
            return getattr(self, key)
        raise KeyError(f"Invalid key: {key}")

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(stage_number="
            f"{getattr(self, 'stage_number', None)!r})"
        )

    def from_dict(self, data: dict) -> None:
        """
        Set several parameters from a dictionary.

        Raises
        ------
        ValueError
            If ``data`` is not a dictionary.
        """
        if isinstance(data, dict):
            for key, value in data.items():
                self[key] = value
        else:
            raise ValueError("Input data must be a dictionary.")

    def to_dict(self) -> dict:
        """Return the parameters as a dictionary."""
        return {key: getattr(self, key) for key in self._KEYMAP}

    def from_json(self, file_name: Any) -> None:
        """
        Load parameters from a JSON file.

        Raises
        ------
        ValueError
            If the file cannot be read or is not valid JSON.
        """
        try:
            with open(file_name, "r") as fobj:
                data = _json.load(fobj)
            self.from_dict(data)
        except (IOError, _json.JSONDecodeError) as error:
            raise ValueError(f"Failed to load JSON file: {error}")

    def to_json(self, file_name: Any) -> None:
        """Save the parameters to a JSON file."""
        with open(file_name, "w") as fobj:
            _json.dump(self.to_dict(), fobj, indent=2)

    def transfer_function(
        self,
        frequency: Any,
        delta: Optional[float] = None,
    ) -> np.ndarray:
        """
        Return the complex transfer function at given frequencies.

        Parameters
        ----------
        frequency
            Frequencies in Hz.
        delta
            Sampling interval of the data, for digital stages.
        """
        raise NotImplementedError(
            f"transfer function not available for {type(self).__name__}"
        )


class StageGain(StageResponse):
    """Gain stage: a constant, real sensitivity."""

    _KEYMAP = {
        "description": (str, None),
        "sensitivity": (float, 1.0),
        "frequency": (float, None),
        "stage_number": (int, None),
    }
    stage_type = "gain"

    def transfer_function(self, frequency, delta=None):
        """
        Return the sensitivity at every frequency.

        Raises
        ------
        ValueError
            If the sensitivity is not set.
        """
        if self.sensitivity is None:
            raise ValueError("Sensitivity is not set.")
        return np.full(np.shape(frequency), self.sensitivity, dtype=complex)


class StagePoleZero(StageResponse):
    """Analog stage defined by poles, zeros and normalization factor."""

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "pz_type": (str, None),
        "normalization_factor": (float, None),
        "normalization_frequency": (float, None),
        "poles": (complex, None),
        "zeros": (complex, None),
        "stage_number": (int, None),
    }
    stage_type = "paz"

    def transfer_function(self, frequency, delta=None):
        """
        Return ``A0 * prod(s - z) / prod(s - p)`` at ``s = 2j*pi*f``
        (``s = j*f`` if the poles and zeros are given in Hz).

        Raises
        ------
        ValueError
            If the poles are not defined.
        """
        if self.poles is None:
            raise ValueError("Poles must be defined.")

        zeros = [] if self.zeros is None else self.zeros
        a0 = 1.0 if self.normalization_factor is None else (
            self.normalization_factor
        )

        # Poles and zeros given in Hz (StationXML "LAPLACE (HERTZ)"):
        # evaluate at s = j*f, i.e. with the frequency scaled by 1/2pi.
        frequency = np.asarray(frequency, dtype=float)
        if self.pz_type and "HERTZ" in self.pz_type.upper():
            frequency = frequency / (2.0 * np.pi)

        return paz_transfer_function(a0, self.poles, zeros, frequency)


class StagePolynomial(StageResponse):
    """
    Polynomial stage.

    Evaluated as ``H(s) = sum(n_k s**k) / sum(d_k s**k)``. Note that a
    StationXML ``<Polynomial>`` describes a non-linear sensor, so this
    stage is not active by default.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "numerator": (float, [1.0]),
        "denominator": (float, [1.0]),
        "stage_number": (int, None),
    }
    stage_type = "polynomial"

    def transfer_function(self, frequency, delta=None):
        """Return the rational function in ``s = 2j*pi*f``."""
        return polynomial_transfer_function(
            self.numerator, self.denominator, frequency
        )


class StageFIR(StageResponse):
    """
    Digital FIR stage.

    Not active by default: a FIR acts at the sampling rate of its own
    stage, which is not yet taken into account.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "simmetry": (str, None),
        "normalization": (str, None),
        "numerator": (float, [1.0]),
        "denominator": (float, [1.0]),
        "stage_number": (int, None),
    }
    stage_type = "fir"

    def transfer_function(self, frequency, delta=None):
        """
        Return the response at given frequencies for sampling ``delta``.

        Raises
        ------
        ValueError
            If ``delta`` is not given.
        """
        if delta is None:
            raise ValueError("the sampling interval is required")

        omega = 2.0 * np.pi * np.asarray(frequency, dtype=float) * delta
        numerator = np.atleast_1d(
            self.numerator if self.numerator is not None else 1.0
        )
        denominator = np.atleast_1d(
            self.denominator if self.denominator is not None else 1.0
        )

        num = np.exp(-1j * np.outer(omega, np.arange(numerator.size)))
        den = np.exp(-1j * np.outer(omega, np.arange(denominator.size)))
        num = num @ numerator
        den = den @ denominator
        den = np.where(np.abs(den) < 1e-12, 1e-12, den)

        return num / den


_STAGE_TYPES = {
    "gain": StageGain,
    "paz": StagePoleZero,
    "fir": StageFIR,
    # StationXML <Coefficients> stages are stored as FIR stages.
    "coefficients": StageFIR,
    "polynomial": StagePolynomial,
}


def _stage_from_type(stage_type: str, data: dict) -> StageResponse:
    """
    Build a stage from its type name and parameters.

    Raises
    ------
    ValueError
        If the stage type is unsupported.
    """
    try:
        return _STAGE_TYPES[stage_type](data)
    except KeyError:
        raise ValueError(f"Unsupported stage_type: {stage_type}") from None


# ----------------------------------------------------------------------------
# Selection
# ----------------------------------------------------------------------------


def select_response(
    resp: Any,
    sid: Optional[str] = None,
    time: Any = None,
) -> Any:
    """
    Return the response applicable to a channel at a time.

    Parameters
    ----------
    resp
        Stage, stage set, stream response or response collection.
    sid
        Source identifier of the data.
    time
        Start time of the data.

    Returns
    -------
    StageResponse or StageSet
        Object providing ``transfer_function``. Stages and stage sets
        are returned unchanged.

    Raises
    ------
    TypeError
        If the response type is unsupported.
    ValueError
        If the identifier or the time needed for the selection are
        missing, or if a stream response belongs to another channel.
    KeyError
        If the channel or the epoch are not available.
    """
    if isinstance(resp, (StageResponse, StageSet)):
        return resp

    if isinstance(resp, StreamResponse):
        if resp.sid is not None and sid is not None and resp.sid != sid:
            raise ValueError(
                f"the response of {resp.sid} cannot be applied to {sid}"
            )
        if time is None:
            raise ValueError(
                "start time is required to select a response epoch"
            )
        return resp.extract(time)

    if isinstance(resp, ResponseCollection):
        if sid is None:
            raise ValueError(
                "source identifier is required to select a response"
            )
        return select_response(resp[sid], sid=sid, time=time)

    raise TypeError(f"unsupported response object: {type(resp).__name__}")


# ----------------------------------------------------------------------------
# Transfer functions
# ----------------------------------------------------------------------------


def paz_transfer_function(a0, poles, zeros, freq):
    """
    Return the transfer function of a poles-and-zeros system.

    Parameters
    ----------
    a0 : float
        Normalization factor.
    poles, zeros : array-like of complex
        Poles and zeros in rad/s.
    freq : array-like
        Frequencies in Hz.

    Returns
    -------
    numpy.ndarray
        ``a0 * prod(s - z) / prod(s - p)`` at ``s = 2j*pi*f``.
    """
    s = 2j * np.pi * np.asarray(freq, dtype=float)

    numerator = np.ones_like(s, dtype=np.complex128)
    denominator = np.ones_like(s, dtype=np.complex128)

    for zero in (zeros if zeros is not None else []):
        numerator *= s - zero

    for pole in (poles if poles is not None else []):
        denominator *= s - pole

    denominator = np.where(np.abs(denominator) < 1e-12, 1e-12, denominator)

    return a0 * numerator / denominator


def polynomial_transfer_function(ncoeff, dcoeff, freq):
    """
    Return a rational function of ``s = 2j*pi*f``.

    Parameters
    ----------
    ncoeff, dcoeff : array-like
        Numerator and denominator coefficients, ascending order.
    freq : array-like
        Frequencies in Hz.

    Returns
    -------
    numpy.ndarray
        ``sum(n_k s**k) / sum(d_k s**k)``.
    """
    s = 2j * np.pi * np.asarray(freq, dtype=float)
    ncoeff = np.atleast_1d(ncoeff)
    dcoeff = np.atleast_1d(dcoeff)

    num = sum(c * s ** k for k, c in enumerate(ncoeff))
    den = sum(c * s ** k for k, c in enumerate(dcoeff))
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)

    return num / den


def fir_transfer_function(ncoeff, dcoeff, nsamp, delta=1.0):
    """
    Return the response of a digital filter on ``nsamp`` frequencies.

    Parameters
    ----------
    ncoeff, dcoeff : array-like
        Numerator (impulse response) and denominator coefficients.
    nsamp : int
        Number of frequencies, from 0 to the Nyquist frequency.
    delta : float, optional
        Sampling interval in seconds.

    Returns
    -------
    tuple of numpy.ndarray
        Complex response and frequencies in Hz.
    """
    if ncoeff is None or np.isscalar(ncoeff):
        ncoeff = [ncoeff or 1.0]
    if dcoeff is None or np.isscalar(dcoeff):
        dcoeff = [dcoeff or 1.0]

    omega = np.linspace(0, np.pi, nsamp)
    jn = np.arange(len(ncoeff)).reshape(-1, 1)
    jd = np.arange(len(dcoeff)).reshape(-1, 1)

    num = np.sum(np.array(ncoeff)[:, None] * np.exp(-1j * omega * jn), axis=0)
    den = np.sum(np.array(dcoeff)[:, None] * np.exp(-1j * omega * jd), axis=0)
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)

    return num / den, omega / (2 * np.pi * delta)
