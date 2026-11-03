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
(gain, poles and zeros, digital filters, response tables) grouped into
stage sets valid in a time interval (epochs), stream responses holding
the epochs of one channel, and collections of stream responses indexed
by source identifier.

Design choices
--------------
- The response is a model of the instrument: each stage provides its
  transfer function ``H(f)`` (:meth:`StageResponse.transfer_function`)
  and a stage set the product over its stages. The model does not
  depend on waveform classes: responses are applied to data by
  :meth:`shakelab.signals.spectrum.Spectrum.deconvolve_response` and
  the corresponding methods of records, streams and collections.
- By default the complete response is used (all stages). A subset is
  chosen with the ``stages`` argument, e.g. ``("gain", "paz")`` for
  the analog part only.
- A response is never computed wrongly in silence: stages that cannot
  be evaluated (non-linear polynomials, unsupported elements) raise an
  explicit error when used. Inconsistent metadata are reported by
  :meth:`StageSet.validate`, never corrected automatically.
- :func:`select_response` chooses the stage or stage set applicable to
  a source identifier and a time.

Conventions
-----------
- Frequencies in Hz. Analog poles and zeros in rad/s
  (``s = 2j*pi*f``) or in Hz (``s = 1j*f``); digital poles and zeros
  in the z-plane, with ``z = exp(2j*pi*f/fs)``.
- Digital stages are evaluated at their own input sampling rate
  ``fs`` (StationXML ``Decimation/InputSampleRate``). A time
  correction ``correction`` applied by the datalogger multiplies the
  response by ``exp(2j*pi*f*correction)``: a positive correction is a
  time advance (SEED blockette 57, field 8), so a correction equal to
  the delay of a linear-phase filter cancels its phase.
- The definitions follow the SEED manual (v2.4, Appendix C) and the
  FDSN StationXML documentation.
- The input units of a response are those of the ground motion (e.g.
  ``m/s**2``), the output units those of the recorded data. Units are
  kept as written in the metadata.
- The ``info()`` methods print a summary for interactive inspection.
"""

from __future__ import annotations

import json as _json
import math
import os as _os
import warnings
from copy import deepcopy
from typing import Any, Optional

import numpy as np
from numpy.polynomial import polynomial as _poly

from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import (
    cast_value,
    deserialize_complex,
    serialize_ndarray,
)

__all__ = [
    "STAGE_TYPES",
    "ResponseWarning",
    "ResponseCollection",
    "StreamResponse",
    "StageSet",
    "Sensitivity",
    "StageResponse",
    "StageGain",
    "StagePoleZero",
    "StageCoefficients",
    "StageResponseList",
    "StagePolynomial",
    "StageUnsupported",
    "select_response",
    "select_sensitivity",
    "response_transfer",
    "same_units",
    "expand_fir_symmetry",
    "paz_transfer_function",
    "polynomial_transfer_function",
    "fir_transfer_function",
]

STAGE_TYPES = (
    "gain",
    "paz",
    "coefficients",
    "response_list",
    "polynomial",
    "unsupported",
)
"""Stage types, as accepted by the ``stages`` arguments."""


class ResponseWarning(UserWarning):
    """Warning about incomplete or unusual response metadata."""


# ----------------------------------------------------------------------------
# Collections and epochs
# ----------------------------------------------------------------------------


class ResponseCollection:
    """
    Collection of stream responses indexed by source identifier.

    Parameters
    ----------
    byte_stream
        Optional source to read (see :meth:`read`).
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

    def __contains__(self, sid: Any) -> bool:
        return sid in self.sid

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
        Append a copy of a stream response.

        If a stream response with the same identifier is already in the
        collection, the stage sets of ``obj`` are added to it. The
        collection stores copies, so later changes to ``obj`` do not
        affect it.

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
            self.response.append(deepcopy(obj))

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

    def validate(self, rtol: float = 1e-3) -> list:
        """
        Check the consistency of all epochs (see :meth:`StageSet.validate`).

        Returns
        -------
        list of str
            Problems found, each prefixed by channel and epoch start.
        """
        issues = []
        for resp in self.response:
            for stage_set in resp.stageset:
                prefix = f"{resp.sid} [{stage_set.starttime}]"
                issues.extend(
                    f"{prefix}: {issue}"
                    for issue in stage_set.validate(rtol=rtol)
                )
        return issues

    def read(
        self,
        byte_stream: Any,
        ftype: Optional[str] = None,
        append: bool = False,
    ) -> None:
        """
        Read responses from a StationXML, RESP or JSON file.

        Parameters
        ----------
        byte_stream
            File path; for StationXML and RESP also text, bytes or an
            open file (see
            :func:`shakelab.signals.stationxml.read_stationxml` and
            :func:`shakelab.signals.respfile.read_resp`).
        ftype
            ``"xml"`` (also ``"sxml"``), ``"resp"`` or ``"json"``. If
            None, it is inferred from the file name: ``"resp"`` for
            names starting with ``RESP.`` (as written by rdseed) or
            ending with ``.resp``, otherwise the extension.
        append
            If True, add to the current content; otherwise replace it.

        Raises
        ------
        ValueError
            If the file type is unknown or a stage type is unsupported.
        """
        if ftype is None:
            ftype = _file_type(byte_stream)
        ftype = ftype.lower()

        if ftype in ("xml", "sxml"):
            # Local import: the StationXML reader builds response
            # objects, so it sits above this module.
            from shakelab.signals.stationxml import read_stationxml

            new = read_stationxml(byte_stream)

        elif ftype == "resp":
            from shakelab.signals.respfile import read_resp

            new = read_resp(byte_stream)

        elif ftype == "json":
            with open(byte_stream, "r") as fobj:
                new = _collection_from_json(_json.load(fobj))

        else:
            raise ValueError(f"Unknown file type: {ftype}")

        if not append:
            self.response = []
        for resp in new.response:
            self.append(resp)

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
                polynomial = stage_set.instrument_polynomial
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
                    "sample_rate": stage_set.sample_rate,
                    "sensitivity": (
                        stage_set.sensitivity.to_dict()
                        if stage_set.sensitivity else None
                    ),
                    "instrument_polynomial": (
                        serialize_ndarray(polynomial.to_dict())
                        if polynomial else None
                    ),
                    "stages": [
                        {
                            "stage_number": stage.stage_number,
                            "stage_type": stage.stage_type,
                            "data": serialize_ndarray(stage.to_dict()),
                        }
                        for stage in stage_set.stage
                    ],
                })

        with open(file, "w") as fobj:
            _json.dump(out, fobj, indent=2)

    def info(self) -> None:
        """Print a summary of all stream responses."""
        for resp in self.response:
            resp.info()


def _file_type(path: Any) -> str:
    """Return the response file type inferred from a file name."""
    name = _os.path.basename(str(path))
    if name.upper().startswith("RESP."):
        return "resp"
    _, ext = _os.path.splitext(name)
    return ext.lstrip(".").lower()


def _collection_from_json(content: list) -> ResponseCollection:
    """Build a collection from the content of a JSON file."""
    collection = ResponseCollection()

    for entry in content:
        stage_set = StageSet(
            entry.get("starttime"),
            entry.get("endtime"),
            sample_rate=entry.get("sample_rate"),
        )
        if entry.get("sensitivity"):
            stage_set.sensitivity = Sensitivity(**entry["sensitivity"])
        if entry.get("instrument_polynomial"):
            stage_set.instrument_polynomial = _stage_from_type(
                "polynomial", entry["instrument_polynomial"]
            )

        for item in entry["stages"]:
            data = deserialize_complex(item["data"])
            data["stage_number"] = item["stage_number"]
            stage_set.append(_stage_from_type(item["stage_type"], data))

        stream_response = StreamResponse(entry["sid"])
        stream_response.append(stage_set)
        collection.append(stream_response)

    return collection


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
        Append copies of one or more stage sets.

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
            self.stageset.append(item.copy())

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


class Sensitivity:
    """
    Overall sensitivity of a response.

    Parameters
    ----------
    value
        Sensitivity (output units per input unit).
    frequency
        Frequency in Hz at which the value applies.
    input_units, output_units
        Units of the ground motion and of the recorded data.
    """

    def __init__(
        self,
        value: Optional[float] = None,
        frequency: Optional[float] = None,
        input_units: Optional[str] = None,
        output_units: Optional[str] = None,
    ):
        self.value = None if value is None else float(value)
        self.frequency = None if frequency is None else float(frequency)
        self.input_units = input_units
        self.output_units = output_units

    def __repr__(self) -> str:
        return (
            f"Sensitivity(value={self.value}, frequency={self.frequency}, "
            f"input_units={self.input_units!r}, "
            f"output_units={self.output_units!r})"
        )

    def to_dict(self) -> dict:
        """Return the attributes as a dictionary."""
        return {
            "value": self.value,
            "frequency": self.frequency,
            "input_units": self.input_units,
            "output_units": self.output_units,
        }


class StageSet:
    """
    Stages of a response, valid between two times.

    Parameters
    ----------
    starttime, endtime
        Validity interval (``endtime`` excluded); None means unbounded.
    sample_rate
        Sampling rate of the recorded data in Hz, if known.

    Attributes
    ----------
    stage : list of StageResponse
        Stages in order of application.
    sensitivity : Sensitivity or None
        Declared overall sensitivity.
    instrument_polynomial : StagePolynomial or None
        Declared overall polynomial of a non-linear instrument.
    """

    def __init__(
        self,
        starttime: Any = None,
        endtime: Any = None,
        sample_rate: Optional[float] = None,
    ):
        if starttime is not None and not isinstance(starttime, Date):
            starttime = Date(starttime)
        if endtime is not None and not isinstance(endtime, Date):
            endtime = Date(endtime)

        self.starttime = starttime
        self.endtime = endtime
        self.sample_rate = (
            None if sample_rate is None else float(sample_rate)
        )
        self.stage: list[StageResponse] = []
        self.sensitivity: Optional[Sensitivity] = None
        self.instrument_polynomial: Optional[StagePolynomial] = None

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

        if self.sensitivity is not None:
            sens = self.sensitivity
            print(
                f"Sensitivity: {sens.value} {sens.output_units} per "
                f"{sens.input_units} at {sens.frequency} Hz"
            )

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
        Return the first stage with the given stage number.

        Raises
        ------
        IndexError
            If the stage is not found.
        """
        for stage in self.stage:
            if stage.stage_number == stage_number:
                return stage
        raise IndexError(f"Stage {stage_number} not found")

    def select_stages(self, stages: Any = None) -> list:
        """
        Return the stages of the given types.

        Parameters
        ----------
        stages
            Stage type or sequence of stage types (see
            ``STAGE_TYPES``); None selects all stages.

        Raises
        ------
        ValueError
            If a stage type is unknown.
        """
        if stages is None:
            return list(self.stage)
        if isinstance(stages, str):
            stages = (stages,)
        unknown = set(stages) - set(STAGE_TYPES)
        if unknown:
            raise ValueError(f"unknown stage types: {sorted(unknown)}")
        return [s for s in self.stage if s.stage_type in stages]

    @property
    def input_units(self) -> Optional[str]:
        """
        Units of the ground motion.

        Input units of the first stage that declares units (gain stages
        never do), or of the overall sensitivity.
        """
        units = _chain_units(self.stage)[0]
        if not units and self.sensitivity is not None:
            units = self.sensitivity.input_units
        return units or None

    @property
    def output_units(self) -> Optional[str]:
        """
        Units of the recorded data.

        Output units of the last stage that declares units (gain stages
        never do), or of the overall sensitivity.
        """
        units = _chain_units(self.stage)[1]
        if not units and self.sensitivity is not None:
            units = self.sensitivity.output_units
        return units or None

    def units(self, stages: Any = None) -> tuple:
        """
        Return the input and output units of the selected stages.

        Parameters
        ----------
        stages
            Stage types (see :meth:`select_stages`); None for the
            complete response, which falls back on the units of the
            overall sensitivity.

        Returns
        -------
        tuple
            ``(input_units, output_units)``, None where not declared.
        """
        if stages is None:
            return self.input_units, self.output_units
        return _chain_units(self.select_stages(stages))

    def transfer_function(
        self,
        frequency: Any,
        delta: Optional[float] = None,
        stages: Any = None,
    ) -> np.ndarray:
        """
        Return the transfer function of the selected stages.

        Parameters
        ----------
        frequency
            Frequencies in Hz.
        delta
            Sampling interval of the data, used by digital stages that
            do not declare their own sampling rate.
        stages
            Stage types to include (see :meth:`select_stages`); None
            for the complete response.

        Returns
        -------
        numpy.ndarray
            Product of the transfer functions of the stages.

        Raises
        ------
        ValueError, NotImplementedError
            If a selected stage cannot be evaluated.

        Notes
        -----
        The gain of a stage is its amplitude at the gain frequency (SEED
        v2.4, blockette 58; StationXML ``StageGain``). Poles-and-zeros
        and coefficient stages are therefore normalized to unit
        amplitude at the frequency of the gain of the same stage,
        ``T(f) / |T(f_gain)|``, so that the stage amplitude at
        ``f_gain`` is the declared gain also when ``f_gain`` differs from
        the normalization frequency of the poles and zeros (or when the
        digital filter does not have unit gain at ``f_gain``), as in
        evalresp. With consistent metadata the factor is 1. If the
        amplitude at ``f_gain`` is zero or not finite, the stage is used
        as declared, with a :class:`ResponseWarning`.
        """
        gain_frequency = self._gain_frequencies()
        result = np.ones(np.shape(frequency), dtype=complex)
        for stage in self.select_stages(stages):
            value = stage.transfer_function(frequency, delta)
            reference = gain_frequency.get(stage.stage_number)
            if reference is not None and isinstance(
                stage, (StagePoleZero, StageCoefficients)
            ):
                value = value / _stage_amplitude(stage, reference, delta)
            result = result * value
        return result

    def _gain_frequencies(self) -> dict:
        """Return the gain frequency of each stage number."""
        frequencies = {}
        for stage in self.stage:
            if (
                isinstance(stage, StageGain)
                and stage.stage_number is not None
                and stage.frequency is not None
            ):
                frequencies.setdefault(stage.stage_number, stage.frequency)
        return frequencies

    def validate(self, rtol: float = 1e-3) -> list:
        """
        Check the consistency of the response metadata.

        The checks are: amplitude of the total response at the
        sensitivity frequency against the declared sensitivity, and the
        product of the stage gains when they refer to that frequency;
        normalization of poles-and-zeros stages; continuity of units
        between stages; decimation chain against the sampling rate;
        symmetry of FIR filters whose delay implies linear phase.

        Parameters
        ----------
        rtol
            Relative tolerance of the numerical checks.

        Returns
        -------
        list of str
            Problems found (empty if none). Nothing is modified.
        """
        issues = []
        sens = self.sensitivity

        if sens is not None and sens.value:
            # The product of the gains equals the sensitivity only if
            # all gains refer to the sensitivity frequency; otherwise
            # the amplitude of the total response is the check.
            gains = [s for s in self.stage if isinstance(s, StageGain)]
            if gains and all(
                g.frequency is None or sens.frequency is None
                or math.isclose(
                    g.frequency, sens.frequency, abs_tol=1e-9
                )
                for g in gains
            ):
                product = float(np.prod([g.sensitivity for g in gains]))
                if not math.isclose(product, sens.value, rel_tol=rtol):
                    issues.append(
                        f"product of stage gains {product:.6g} differs "
                        f"from the sensitivity {sens.value:.6g}"
                    )

            if sens.frequency is not None:
                try:
                    amplitude = abs(self.transfer_function(
                        np.array([sens.frequency]),
                        delta=_delta(self.sample_rate),
                    )[0])
                except (ValueError, NotImplementedError) as error:
                    issues.append(f"total response not computable: {error}")
                else:
                    if not math.isclose(amplitude, sens.value, rel_tol=rtol):
                        issues.append(
                            f"total response amplitude {amplitude:.6g} at "
                            f"{sens.frequency} Hz differs from the "
                            f"sensitivity {sens.value:.6g}"
                        )

        for stage in self.stage:
            issues.extend(stage.validate(rtol=rtol))

        declared = [
            s for s in self.stage
            if getattr(s, "input_units", None)
            and getattr(s, "output_units", None)
        ]
        for first, second in zip(declared, declared[1:]):
            if not same_units(first.output_units, second.input_units):
                issues.append(
                    f"stage {second.stage_number}: input units "
                    f"{second.input_units!r} do not match the output "
                    f"units {first.output_units!r} of stage "
                    f"{first.stage_number}"
                )
        if sens is not None and declared:
            if sens.input_units and not same_units(
                sens.input_units, declared[0].input_units
            ):
                issues.append(
                    f"sensitivity input units {sens.input_units!r} differ "
                    f"from those of the first stage "
                    f"{declared[0].input_units!r}"
                )
            if sens.output_units and not same_units(
                sens.output_units, declared[-1].output_units
            ):
                issues.append(
                    f"sensitivity output units {sens.output_units!r} "
                    f"differ from those of the last stage "
                    f"{declared[-1].output_units!r}"
                )

        rate = None
        for stage in self.stage:
            fs = getattr(stage, "input_sample_rate", None)
            if fs is None:
                continue
            if rate is not None and not math.isclose(fs, rate, rel_tol=rtol):
                issues.append(
                    f"stage {stage.stage_number}: input sampling rate {fs} "
                    f"Hz differs from the output rate {rate} Hz of the "
                    f"previous digital stage"
                )
            rate = fs / (stage.decimation_factor or 1)
        if rate is not None and self.sample_rate is not None:
            if not math.isclose(rate, self.sample_rate, rel_tol=rtol):
                issues.append(
                    f"decimation chain ends at {rate} Hz, but the "
                    f"sampling rate is {self.sample_rate} Hz"
                )

        return issues

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
        KeyError
            If a key is not a parameter of the stage.
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
            _json.dump(serialize_ndarray(self.to_dict()), fobj, indent=2)

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
            Sampling interval of the data, for digital stages that do
            not declare their own sampling rate.
        """
        raise NotImplementedError(
            f"transfer function not available for {type(self).__name__}"
        )

    def validate(self, rtol: float = 1e-3) -> list:
        """Return the problems found in the stage parameters."""
        return []


_DECIMATION_KEYMAP = {
    "input_sample_rate": (float, None),
    "decimation_factor": (int, None),
    "decimation_offset": (int, None),
    "delay": (float, None),
    "correction": (float, None),
}
"""Parameters of digital stages (StationXML ``Decimation``)."""


def _delta(sample_rate: Optional[float]) -> Optional[float]:
    """Return the sampling interval of a sampling rate (None if unknown)."""
    return None if not sample_rate else 1.0 / sample_rate


def _normalized_frequency(stage, frequency, delta) -> np.ndarray:
    """
    Return ``f / fs`` for a digital stage.

    Raises
    ------
    ValueError
        If neither the stage nor ``delta`` give the sampling rate.
    """
    fs = getattr(stage, "input_sample_rate", None)
    if not fs:
        if delta is None:
            raise ValueError(
                f"stage {stage.stage_number}: the sampling rate of the "
                "digital stage is unknown (give delta)"
            )
        fs = 1.0 / delta
    return np.asarray(frequency, dtype=float) / fs


def _time_correction(stage, frequency) -> Any:
    """Return ``exp(2j*pi*f*correction)`` (1 without correction)."""
    correction = getattr(stage, "correction", None)
    if not correction:
        return 1.0
    return np.exp(2j * np.pi * np.asarray(frequency, dtype=float)
                  * correction)


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
    """
    Stage defined by poles, zeros and normalization factor.

    ``pz_type`` follows StationXML: ``"LAPLACE (RADIANS/SECOND)"``
    (default), ``"LAPLACE (HERTZ)"`` or ``"DIGITAL (Z-TRANSFORM)"``.
    Digital stages also carry the decimation parameters.
    """

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
        **_DECIMATION_KEYMAP,
    }
    stage_type = "paz"

    @property
    def is_digital(self) -> bool:
        """Whether poles and zeros are given in the z-plane."""
        return bool(self.pz_type) and "DIGITAL" in self.pz_type.upper()

    def transfer_function(self, frequency, delta=None):
        """
        Return ``A0 * prod(x - z) / prod(x - p)``.

        ``x`` is ``s = 2j*pi*f`` (rad/s), ``s = 1j*f`` (Hz) or
        ``z = exp(2j*pi*f/fs)`` (digital).

        Missing poles or zeros mean none (e.g. a stage with gain only:
        SEED allows zero poles and zeros).

        Raises
        ------
        ValueError
            If the sampling rate of a digital stage is unknown.
        """
        zeros = [] if self.zeros is None else self.zeros
        poles = [] if self.poles is None else self.poles
        a0 = 1.0 if self.normalization_factor is None else (
            self.normalization_factor
        )
        frequency = np.asarray(frequency, dtype=float)

        if self.is_digital:
            z = np.exp(2j * np.pi * _normalized_frequency(
                self, frequency, delta
            ))
            return (
                a0 * _root_product(z, zeros) / _root_product(z, poles)
                * _time_correction(self, frequency)
            )

        # Poles and zeros given in Hz: evaluate at s = j*f, i.e. with
        # the frequency scaled by 1/2pi.
        if self.pz_type and "HERTZ" in self.pz_type.upper():
            frequency = frequency / (2.0 * np.pi)

        return paz_transfer_function(a0, poles, zeros, frequency)

    def validate(self, rtol=1e-3):
        """Check that ``|A0 * H(fn)| = 1`` at the normalization frequency."""
        if self.normalization_frequency is None:
            return []
        try:
            value = abs(self.transfer_function(
                np.array([self.normalization_frequency])
            )[0])
        except ValueError:
            return []
        if math.isclose(value, 1.0, rel_tol=rtol):
            return []
        return [
            f"stage {self.stage_number}: normalized poles-and-zeros "
            f"amplitude is {value:.6f} at {self.normalization_frequency} Hz "
            "(expected 1)"
        ]


class StageCoefficients(StageResponse):
    """
    Stage given as a ratio of polynomials (StationXML ``Coefficients``).

    ``cf_type`` follows StationXML:

    - ``"DIGITAL"`` (default): ``H(z) = N(z**-1) / D(z**-1)`` with
      ``z = exp(2j*pi*f/fs)``, coefficients in forward order
      (``numerator[k]`` multiplies ``z**-k``), evaluated at the input
      sampling rate of the stage (SEED Appendix C, eqs. 9 and 11). It
      also represents StationXML ``FIR`` (after expansion of the
      symmetry). Without coefficients ``H = 1``, e.g. an
      analog-to-digital converter, which only changes the units.
    - ``"ANALOG (RADIANS/SECOND)"`` and ``"ANALOG (HERTZ)"``:
      ``H(s) = sum(n_k s**k) / sum(d_k s**k)`` in ascending powers,
      with ``s = 2j*pi*f`` or ``s = 1j*f`` (SEED Appendix C, eq. 7).

    Digital stages carry the decimation parameters.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "cf_type": (str, "DIGITAL"),
        "numerator": (float, [1.0]),
        "denominator": (float, [1.0]),
        "stage_number": (int, None),
        **_DECIMATION_KEYMAP,
    }
    stage_type = "coefficients"

    @property
    def is_digital(self) -> bool:
        """Whether the coefficients define a digital filter."""
        return not self.cf_type or "DIGITAL" in self.cf_type.upper()

    def transfer_function(self, frequency, delta=None):
        """
        Return the response at given frequencies.

        Raises
        ------
        ValueError
            If the type is unknown, or the sampling rate of a digital
            stage is not declared and ``delta`` is not given.
        """
        numerator = _coefficients(self.numerator)
        denominator = _coefficients(self.denominator)

        if not self.is_digital:
            cf_type = self.cf_type.upper()
            frequency = np.asarray(frequency, dtype=float)
            if "HERTZ" in cf_type:
                frequency = frequency / (2.0 * np.pi)
            elif "RADIANS" not in cf_type:
                raise ValueError(
                    f"stage {self.stage_number}: unknown coefficients "
                    f"type {self.cf_type!r}"
                )
            return polynomial_transfer_function(
                numerator, denominator, frequency
            )

        nu = _normalized_frequency(self, frequency, delta)
        # Polynomials in w = z**-1, evaluated with Horner's scheme.
        w = np.exp(-2j * np.pi * nu)
        num = _poly.polyval(w, numerator)
        den = _poly.polyval(w, denominator)
        den = np.where(np.abs(den) < 1e-12, 1e-12, den)
        return num / den * _time_correction(self, frequency)

    def validate(self, rtol=1e-3):
        """Check the symmetry of filters whose delay implies linear phase."""
        numerator = _coefficients(self.numerator)
        denominator = _coefficients(self.denominator)
        fs = self.input_sample_rate
        if (
            not self.is_digital or numerator.size < 2
            or denominator.size != 1 or not fs or self.delay is None
        ):
            return []
        linear_delay = (numerator.size - 1) / 2.0 / fs
        if not math.isclose(self.delay, linear_delay, rel_tol=1e-6):
            return []
        difference = np.abs(numerator - numerator[::-1]).max()
        if difference <= rtol * np.abs(numerator).max():
            return []
        return [
            f"stage {self.stage_number}: delay {self.delay} s implies a "
            "linear-phase filter, but the coefficients are not symmetric "
            f"(largest difference {difference:.3g})"
        ]


class StageResponseList(StageResponse):
    """
    Stage given as a table of amplitude and phase (degrees).

    The response is interpolated linearly in log-amplitude versus
    log-frequency and in unwrapped phase versus frequency. Outside the
    tabulated range it is zero (no extrapolation), with a warning.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "frequency": (float, None),
        "amplitude": (float, None),
        "phase": (float, None),
        "stage_number": (int, None),
    }
    stage_type = "response_list"

    def transfer_function(self, frequency, delta=None):
        """
        Return the interpolated response.

        Raises
        ------
        ValueError
            If the table is empty.
        """
        if self.frequency is None or self.amplitude is None:
            raise ValueError(
                f"stage {self.stage_number}: empty response list"
            )

        table_f = np.atleast_1d(np.asarray(self.frequency, dtype=float))
        order = np.argsort(table_f)
        table_f = table_f[order]
        table_a = np.atleast_1d(np.asarray(self.amplitude, float))[order]
        table_p = (
            np.zeros(table_f.size) if self.phase is None
            else np.atleast_1d(np.asarray(self.phase, float))[order]
        )
        table_p = np.unwrap(np.deg2rad(table_p))

        frequency = np.asarray(frequency, dtype=float)
        result = np.zeros(frequency.shape, dtype=complex)
        inside = (frequency >= table_f[0]) & (frequency <= table_f[-1])

        if np.all(table_f > 0) and np.all(table_a > 0):
            inside &= frequency > 0
            amplitude = np.exp(np.interp(
                np.log(frequency[inside]), np.log(table_f), np.log(table_a)
            ))
        else:
            amplitude = np.interp(frequency[inside], table_f, table_a)
        phase = np.interp(frequency[inside], table_f, table_p)
        result[inside] = amplitude * np.exp(1j * phase)

        if not np.all(inside):
            warnings.warn(
                f"stage {self.stage_number}: response list defined from "
                f"{table_f[0]} to {table_f[-1]} Hz; set to zero outside",
                ResponseWarning,
                stacklevel=2,
            )
        return result


class StagePolynomial(StageResponse):
    """
    Non-linear stage: ``x = sum(a_k * y**k)`` (MacLaurin series).

    The polynomial gives the input of the stage ``x`` (in
    ``input_units``, e.g. a temperature or a pressure) as a function of
    its output ``y`` (in ``output_units``, e.g. Volts): it is applied to
    the recorded values to recover the physical quantity (SEED v2.4,
    blockette 62, notes 5-7; FDSN StationXML ``Polynomial`` and
    ``InstrumentPolynomial``). It is valid within the approximation
    bounds and is not a frequency response: :meth:`evaluate` applies it
    to output values, while :meth:`transfer_function` exists only for a
    linear polynomial.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "approximation_type": (str, "MACLAURIN"),
        "frequency_lower_bound": (float, None),
        "frequency_upper_bound": (float, None),
        "approximation_lower_bound": (float, None),
        "approximation_upper_bound": (float, None),
        "maximum_error": (float, None),
        "coefficients": (float, None),
        "stage_number": (int, None),
    }
    stage_type = "polynomial"

    def evaluate(self, y: Any) -> np.ndarray:
        """
        Return the input values ``x`` for output values ``y``.

        Parameters
        ----------
        y
            Values in the output units of the stage (e.g. Volts).

        Returns
        -------
        numpy.ndarray
            Values in the input units (e.g. a physical quantity).

        Raises
        ------
        ValueError
            If no coefficients are defined.
        """
        if self.coefficients is None:
            raise ValueError("Polynomial coefficients are not defined.")
        return _poly.polyval(np.asarray(y, dtype=float),
                             _coefficients(self.coefficients))

    def transfer_function(self, frequency, delta=None):
        """
        Return the gain ``1 / a_1`` of a linear polynomial.

        For ``x = a_0 + a_1 * y`` the stage converts input to output with
        gain ``dy/dx = 1 / a_1``, as every other stage; the offset
        ``a_0`` has no frequency response and is ignored.

        Raises
        ------
        ValueError
            If the polynomial is not linear or ``a_1`` is zero.
        """
        coefficients = _coefficients(self.coefficients)
        if coefficients.size > 2 and np.any(coefficients[2:] != 0):
            raise ValueError(
                f"stage {self.stage_number}: non-linear polynomial "
                "response cannot be applied in the frequency domain"
            )
        slope = coefficients[1] if coefficients.size > 1 else 0.0
        if slope == 0:
            raise ValueError(
                f"stage {self.stage_number}: polynomial without linear "
                "term has no gain"
            )
        return np.full(np.shape(frequency), 1.0 / slope, dtype=complex)


class StageUnsupported(StageResponse):
    """
    Placeholder for a response element that cannot be evaluated.

    It keeps the units and position of the element in the response,
    so that the response is known to be incomplete; using it raises
    ``NotImplementedError``.
    """

    _KEYMAP = {
        "description": (str, None),
        "input_units": (str, None),
        "output_units": (str, None),
        "element": (str, None),
        "stage_number": (int, None),
    }
    stage_type = "unsupported"

    def transfer_function(self, frequency, delta=None):
        """Raise ``NotImplementedError``."""
        raise NotImplementedError(
            f"stage {self.stage_number}: response element "
            f"{self.element!r} is not supported"
        )


_STAGE_TYPES = {
    "gain": StageGain,
    "paz": StagePoleZero,
    "coefficients": StageCoefficients,
    # Name used by earlier versions for digital coefficients.
    "fir": StageCoefficients,
    "response_list": StageResponseList,
    "polynomial": StagePolynomial,
    "unsupported": StageUnsupported,
}


def _stage_from_type(stage_type: str, data: dict) -> StageResponse:
    """
    Build a stage from its type name and parameters.

    Parameters not defined for the stage type are ignored, so that
    files written by earlier versions can be read.

    Raises
    ------
    ValueError
        If the stage type is unsupported.
    """
    try:
        cls = _STAGE_TYPES[stage_type]
    except KeyError:
        raise ValueError(f"Unsupported stage_type: {stage_type}") from None
    return cls({k: v for k, v in data.items() if k in cls._KEYMAP})


def _stage_amplitude(stage, frequency: float, delta) -> float:
    """
    Return the amplitude of a stage at the frequency of its gain.

    Returns 1 (no normalization) with a warning if the amplitude is
    zero or not finite, e.g. a gain declared at 0 Hz for a stage with
    zeros at the origin.
    """
    value = abs(stage.transfer_function(np.array([frequency]), delta)[0])
    if value == 0 or not np.isfinite(value):
        warnings.warn(
            f"stage {stage.stage_number}: amplitude {value} at the gain "
            f"frequency {frequency} Hz; the stage is not normalized",
            ResponseWarning,
            stacklevel=3,
        )
        return 1.0
    return value


def _chain_units(stages: list) -> tuple:
    """
    Return the input and output units of a chain of stages.

    Gain stages carry no units: in StationXML they belong to the stage
    whose transfer element declares the units, so they are skipped.
    """
    declared = [s for s in stages if not isinstance(s, StageGain)]
    if not declared:
        return None, None
    return (
        getattr(declared[0], "input_units", None) or None,
        getattr(declared[-1], "output_units", None) or None,
    )


def _coefficients(values: Any) -> np.ndarray:
    """Return coefficients as a 1-D array (``[1.0]`` if missing)."""
    if values is None:
        return np.array([1.0])
    return np.atleast_1d(np.asarray(values, dtype=float))


def _root_product(x: np.ndarray, roots: Any) -> np.ndarray:
    """Return ``prod(x - r)`` over the roots."""
    result = np.ones_like(x, dtype=complex)
    for root in (roots if roots is not None else []):
        result = result * (x - root)
    return result


# ----------------------------------------------------------------------------
# Selection and units
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


def select_sensitivity(
    resp: Any,
    sid: Optional[str] = None,
    time: Any = None,
) -> Sensitivity:
    """
    Return the overall sensitivity applicable to a channel at a time.

    Parameters
    ----------
    resp
        Gain stage, stage set, stream response or response collection.
    sid, time
        Source identifier and start time of the data (see
        :func:`select_response`).

    Raises
    ------
    ValueError
        If the stage set declares no sensitivity.
    TypeError
        If the selected response is not a stage set or a gain stage.
    """
    selected = select_response(resp, sid=sid, time=time)

    if isinstance(selected, StageSet):
        sens = selected.sensitivity
        if sens is None or not sens.value:
            raise ValueError("the response declares no overall sensitivity")
        return sens

    if isinstance(selected, StageGain):
        return Sensitivity(selected.sensitivity, selected.frequency)

    raise TypeError(
        f"no sensitivity available from {type(selected).__name__}"
    )


def response_transfer(
    resp: Any,
    frequency: Any,
    delta: Optional[float] = None,
    stages: Any = None,
) -> tuple:
    """
    Return transfer function and units of a selected response.

    Parameters
    ----------
    resp
        Stage or stage set (see :func:`select_response`).
    frequency
        Frequencies in Hz.
    delta
        Sampling interval of the data.
    stages
        Stage types of a stage set to include (None: all); ignored for
        a single stage.

    Returns
    -------
    tuple
        ``(transfer, input_units, output_units)``.
    """
    if isinstance(resp, StageSet):
        transfer = resp.transfer_function(frequency, delta, stages=stages)
        return (transfer, *resp.units(stages))
    transfer = resp.transfer_function(frequency, delta)
    return (
        transfer,
        getattr(resp, "input_units", None) or None,
        getattr(resp, "output_units", None) or None,
    )


def same_units(first: Optional[str], second: Optional[str]) -> bool:
    """
    Return whether two unit names are equivalent.

    The comparison ignores case and surrounding spaces, and treats
    ``count`` and ``counts`` as equal. Units are never rewritten.
    """
    def normalize(units):
        units = (units or "").strip().lower()
        return "count" if units == "counts" else units

    return normalize(first) == normalize(second)


def expand_fir_symmetry(coefficients: Any, symmetry: str = "NONE") -> list:
    """
    Return all the coefficients of a symmetric FIR filter.

    Symmetric filters are stored by their first half only (FDSN
    StationXML ``FIR/Symmetry``; SEED blockette 61, field 5).

    Parameters
    ----------
    coefficients
        Stored coefficients, in forward order.
    symmetry
        ``"NONE"`` (all coefficients given), ``"ODD"`` (odd total
        number: first half and central coefficient given) or
        ``"EVEN"`` (even total number: first half given). Case is
        ignored.

    Returns
    -------
    list of float
        The complete coefficients.

    Raises
    ------
    ValueError
        If the symmetry is unknown.

    Examples
    --------
    >>> expand_fir_symmetry([1.0, 2.0, 3.0], "ODD")
    [1.0, 2.0, 3.0, 2.0, 1.0]
    >>> expand_fir_symmetry([1.0, 2.0], "EVEN")
    [1.0, 2.0, 2.0, 1.0]
    """
    values = [float(c) for c in np.atleast_1d(coefficients)]
    symmetry = (symmetry or "NONE").strip().upper()
    if symmetry == "NONE":
        return values
    if symmetry == "EVEN":
        return values + values[::-1]
    if symmetry == "ODD":
        return values + values[-2::-1]
    raise ValueError(f"Unknown FIR symmetry: {symmetry!r}")


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
