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
Reader of SEED RESP instrument responses.

RESP is the text form of the SEED response blockettes written by
rdseed, by the IRIS/EarthScope Nominal Response Library and by other
tools, and read by evalresp. Each line holds one field of a blockette,
identified by the code ``BxxxFyy`` (blockette ``xxx``, field ``yy``);
repeated fields (e.g. poles) are written one group per line with a
field range such as ``B053F15-18``. The reader converts every channel
epoch into a channel epoch of :mod:`shakelab.signals.metadata` with
the response model of :mod:`shakelab.signals.response`. RESP files
describe only the response: the channel epochs have no coordinates,
orientation or equipment.

Mapping of the blockettes
-------------------------
=================================  =====================================
SEED (RESP)                        ShakeLab
=================================  =====================================
50, 52 (station, channel, epoch)   ``ChannelEpoch`` and ``StageSet``
53 Poles & Zeros                   ``StagePoleZero``
54 Coefficients                    ``StageCoefficients``
55 Response List                   ``StageResponseList``
57 Decimation                      attributes of the stage
58 Sensitivity/Gain, stage > 0     ``StageGain``
58 Sensitivity/Gain, stage 0       ``StageSet.sensitivity``
61 FIR (symmetry expanded)         ``StageCoefficients``
62 Polynomial, stage > 0           ``StagePolynomial``
62 Polynomial, stage 0             ``StageSet.instrument_polynomial``
56 Generic, 60 Reference           ``StageUnsupported``
=================================  =====================================

Design choices
--------------
- Lines are identified by their ``BxxxFyy`` code only: the labels
  after the code differ between writers (e.g. ``Gain`` or
  ``Sensitivity``, ``Input sample rate`` or ``Input sample rate
  (HZ)``) and are ignored. Comment lines start with ``#``.
- The metadata are read as written; units are reduced to their code
  (``M/S - Velocity in Meters Per Second`` becomes ``M/S``) and are
  not rewritten. Inconsistent values are not corrected (see
  :meth:`shakelab.signals.response.StageSet.validate`).
- A blockette continued in a second blockette of the same stage (SEED
  allows it when a blockette exceeds 9999 characters) is merged with
  the first one.
- Blockettes that cannot be evaluated (generic responses, response
  references, composite types) are kept as ``StageUnsupported`` with
  a :class:`shakelab.signals.response.ResponseWarning`, so that an
  incomplete response cannot be used without notice.
- RESP files do not contain the sampling rate of the channel (SEED
  blockette 52, field 18). As in evalresp, it is taken as the output
  rate of the last decimation stage (input rate divided by the
  decimation factor); it is None if no stage declares a decimation.
- The overall sensitivity (blockette 58, stage 0) has no units in
  RESP: they are those of the response, i.e. the input units of the
  first stage and the output units of the last stage that declare
  them (SEED v2.4, blockette 58: "counts per ground motion").

References
----------
- SEED Reference Manual, version 2.4 (2012), chapter 3 (``TIME``
  fields) and chapter 6 (blockettes 50 to 62); Appendix C.
- evalresp, IRIS/EarthScope (https://github.com/EarthScope/evalresp).
"""

from __future__ import annotations

import math
import os
import re
import warnings
from typing import Any, Iterable, Optional

from shakelab.libutils.timeN import Date
from shakelab.signals import metadata as mdm
from shakelab.signals import response as rspm

__all__ = [
    "read_resp",
    "parse_resp_time",
]

_LINE = re.compile(r"^B(\d{3})F(\d{2})(?:-(\d{2}))?(?:\s+(.*))?$")

_TIME = re.compile(
    r"^(\d{4})"
    r"(?:,(\d{1,3})"
    r"(?:,(\d{1,2})"
    r"(?::(\d{1,2})"
    r"(?::(\d{1,2}(?:\.\d*)?))?)?)?)?$"
)

_NO_TIME = ("", "no ending time", "none", "n/a")

_PZ_TYPES = {
    "A": "LAPLACE (RADIANS/SECOND)",
    "B": "LAPLACE (HERTZ)",
    "D": "DIGITAL (Z-TRANSFORM)",
}
"""Blockette 53, field 3."""

_CF_TYPES = {
    "A": "ANALOG (RADIANS/SECOND)",
    "B": "ANALOG (HERTZ)",
    "D": "DIGITAL",
}
"""Blockette 54, field 3."""

_FIR_SYMMETRY = {"A": "NONE", "B": "ODD", "C": "EVEN"}
"""Blockette 61, field 5."""

_POLYNOMIAL_TYPES = {"M": "MACLAURIN"}
"""Blockette 62, field 7."""

_TRANSFER = ("053", "054", "055", "056", "060", "061", "062")
"""Blockettes that define the transfer element of a stage."""

_CHANNEL = ("050", "052")
"""Blockettes of the channel header."""


# ----------------------------------------------------------------------------
# Reader
# ----------------------------------------------------------------------------


def read_resp(source: Any) -> mdm.MetadataCollection:
    """
    Read the instrument responses of a RESP file.

    Parameters
    ----------
    source
        File path (``str`` or path-like), RESP text (``str`` with more
        than one line), ``bytes`` or an open text or binary file.

    Returns
    -------
    MetadataCollection
        One channel epoch per RESP epoch, with its response.

    Raises
    ------
    FileNotFoundError
        If a path does not exist.
    TypeError
        If the source type is not supported.
    ValueError
        If a value cannot be read, or the text holds no response.
    """
    lines = _read_lines(source)
    collection = mdm.MetadataCollection()

    for header, blockettes in _split_epochs(lines):
        sid = ".".join((
            header.get("network", ""),
            header.get("station", ""),
            header.get("location", ""),
            header.get("channel", ""),
        ))
        stage_set = _build_stage_set(sid, blockettes)
        collection.append(mdm.ChannelEpoch(
            sid,
            starttime=parse_resp_time(header.get("starttime")),
            endtime=parse_resp_time(header.get("endtime")),
            sample_rate=stage_set.sample_rate,
            response=stage_set,
        ))

    if not len(collection):
        raise ValueError("No RESP channel response found")

    return collection


def parse_resp_time(text: Optional[str]) -> Optional[Date]:
    """
    Convert a SEED time field into a date.

    SEED times are written as ``YYYY,DDD,HH:MM:SS.FFFF`` (UTC, day of
    year) and may be truncated after any subfield, e.g. ``1987,023``
    (SEED v2.4, chapter 3).

    Parameters
    ----------
    text
        Time field; empty, None or ``"No Ending Time"`` mean an open
        interval.

    Returns
    -------
    Date or None
        The date, or None for an open interval.

    Raises
    ------
    ValueError
        If the text is not a SEED time.

    Examples
    --------
    >>> str(parse_resp_time("1998,299,20"))
    '1998-10-26T20:00:00.000000Z'
    >>> parse_resp_time("No Ending Time") is None
    True
    """
    text = "" if text is None else text.strip()
    if text.lower() in _NO_TIME:
        return None

    match = _TIME.match(text)
    if match is None:
        raise ValueError(f"Not a SEED time: {text!r}")

    year, doy, hour, minute, second = match.groups()
    return Date([
        int(year),
        int(doy or 1),
        int(hour or 0),
        int(minute or 0),
        float(second or 0.0),
    ])


# ----------------------------------------------------------------------------
# Text and structure
# ----------------------------------------------------------------------------


def _read_lines(source: Any) -> list:
    """Return the lines of a RESP source."""
    if isinstance(source, (bytes, bytearray)):
        text = bytes(source).decode("latin-1")
    elif isinstance(source, str) and "\n" in source:
        text = source
    elif isinstance(source, (str, os.PathLike)):
        if not os.path.isfile(source):
            raise FileNotFoundError(f"File not found: {source}")
        with open(source, "r", encoding="latin-1") as fobj:
            text = fobj.read()
    elif hasattr(source, "read"):
        text = source.read()
        if isinstance(text, (bytes, bytearray)):
            text = bytes(text).decode("latin-1")
    else:
        raise TypeError(
            "source must be a path, RESP text, bytes or a file object"
        )
    return text.splitlines()


class _Blockette:
    """Fields of one blockette as read from RESP lines."""

    def __init__(self, btype: str):
        self.btype = btype
        self.fields: dict = {}
        self.rows: dict = {}

    def value(self, field: str, default: Any = None) -> Any:
        """Return the text of a single field."""
        return self.fields.get(field, default)

    def number(self, field: str) -> Optional[float]:
        """
        Return the leading number of a field (None if missing).

        Raises
        ------
        ValueError
            If the field does not start with a number.
        """
        text = self.fields.get(field)
        if text is None or not text.strip():
            return None
        token = text.split()[0]
        try:
            return float(token)
        except ValueError:
            raise ValueError(
                f"blockette {self.btype}, field {field}: not a number: "
                f"{text!r}"
            ) from None

    def integer(self, field: str) -> Optional[int]:
        """Return a field as an integer (None if missing)."""
        value = self.number(field)
        return None if value is None else int(value)

    def code(self, field: str) -> str:
        """Return the first character of a field, upper case."""
        text = (self.fields.get(field) or "").strip()
        return text[:1].upper()

    def units(self, field: str) -> Optional[str]:
        """Return the unit code of a units field (before ``" - "``)."""
        text = (self.fields.get(field) or "").strip()
        if not text:
            return None
        return text.split(" - ", 1)[0].strip() or None

    @property
    def stage_number(self) -> Optional[int]:
        """Stage sequence number."""
        # Field 3 holds the stage number, except in blockettes 53, 54
        # and 62 (field 3 is the type) and 60 (field 3 is the number of
        # stages).
        field = "04" if self.btype in ("053", "054", "060", "062") else "03"
        return self.integer(field)


def _split_epochs(lines: Iterable[str]) -> list:
    """
    Group the lines into channel epochs.

    Returns
    -------
    list of tuple
        ``(header, blockettes)`` for each epoch, where ``header`` holds
        network, station, location, channel and the time interval.
    """
    epochs = []
    header: dict = {}
    blockettes: list = []
    current: Optional[_Blockette] = None

    def close():
        if blockettes:
            epochs.append((dict(header), list(blockettes)))
            blockettes.clear()

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = _LINE.match(line)
        if match is None:
            continue
        btype, field, last, rest = match.groups()
        rest = rest or ""

        if btype in _CHANNEL:
            # A channel header after response blockettes starts a new
            # epoch; header fields not repeated are kept.
            close()
            current = None
            _set_header(header, btype, field, rest)
            continue

        if last is None and ":" in rest:
            # Single field: "label: value". A new blockette starts at
            # its first field (field 3) or when a field repeats.
            value = rest.split(":", 1)[1].strip()
            if (
                current is None or current.btype != btype
                or field == "03" or field in current.fields
            ):
                current = _Blockette(btype)
                blockettes.append(current)
            current.fields[field] = value
        else:
            # Repeated fields: one group of numbers per line.
            if current is None or current.btype != btype:
                current = _Blockette(btype)
                blockettes.append(current)
            key = field if last is None else f"{field}-{last}"
            current.rows.setdefault(key, []).append(
                [_to_float(token, btype, key) for token in rest.split()]
            )

    close()
    return epochs


def _set_header(header: dict, btype: str, field: str, text: str) -> None:
    """Store a field of blockette 50 or 52."""
    value = text.split(":", 1)[1].strip() if ":" in text else text.strip()
    if btype == "050":
        if field == "03":
            header["station"] = value
        elif field == "16":
            header["network"] = value
    elif field == "03":
        # rdseed writes "??" for an empty location code.
        header["location"] = "" if value in ("??", "--") else value
    elif field == "04":
        header["channel"] = value
    elif field == "22":
        header["starttime"] = value
    elif field == "23":
        header["endtime"] = value


def _to_float(token: str, btype: str, field: str) -> float:
    """Convert a token of a repeated field."""
    try:
        return float(token)
    except ValueError:
        raise ValueError(
            f"blockette {btype}, field {field}: not a number: {token!r}"
        ) from None


# ----------------------------------------------------------------------------
# Stages
# ----------------------------------------------------------------------------


def _build_stage_set(sid: str, blockettes: list) -> rspm.StageSet:
    """Convert the blockettes of one epoch into a stage set."""
    stage_set = rspm.StageSet()

    # Collect the blockettes of each stage, merging continuations.
    stages: dict = {}
    for blockette in blockettes:
        number = blockette.stage_number
        if number is None:
            warnings.warn(
                f"{sid}: blockette {blockette.btype} without stage "
                "number ignored",
                rspm.ResponseWarning,
                stacklevel=4,
            )
            continue
        parts = stages.setdefault(number, {})
        if blockette.btype in parts:
            _merge(parts[blockette.btype], blockette)
        else:
            parts[blockette.btype] = blockette

    for number in sorted(stages):
        parts = stages[number]
        if number == 0:
            _set_overall(stage_set, parts)
            continue
        stage_set.append(_stage(sid, number, parts))

    if stage_set.sensitivity is not None:
        stage_set.sensitivity.input_units = stage_set.input_units
        stage_set.sensitivity.output_units = stage_set.output_units

    stage_set.sample_rate = _output_sample_rate(stage_set.stage)
    return stage_set


def _merge(first: _Blockette, second: _Blockette) -> None:
    """Append the repeated fields of a continuation blockette."""
    for key, rows in second.rows.items():
        first.rows.setdefault(key, []).extend(rows)
    for key, value in second.fields.items():
        first.fields.setdefault(key, value)


def _set_overall(stage_set: rspm.StageSet, parts: dict) -> None:
    """Set the overall sensitivity or polynomial (stage 0)."""
    if "058" in parts:
        b058 = parts["058"]
        stage_set.sensitivity = rspm.Sensitivity(
            value=b058.number("04"),
            frequency=b058.number("05"),
        )
    if "062" in parts:
        polynomial = _polynomial(parts["062"])
        polynomial.stage_number = 0
        stage_set.instrument_polynomial = polynomial


def _stage(sid: str, number: int, parts: dict) -> list:
    """Return the stages (transfer element and gain) of one stage."""
    transfers = [btype for btype in _TRANSFER if btype in parts]
    if len(transfers) > 1:
        warnings.warn(
            f"{sid}: stage {number} has several response blockettes "
            f"({', '.join(transfers)}); only {transfers[0]} is used",
            rspm.ResponseWarning,
            stacklevel=5,
        )

    result = []
    transfer = None
    if transfers:
        transfer = _transfer(parts[transfers[0]])
    elif "057" in parts:
        # Decimation without a filter: unit digital response.
        transfer = rspm.StageCoefficients(cf_type="DIGITAL")

    if transfer is not None:
        transfer.stage_number = number
        if "057" in parts and "input_sample_rate" in transfer._KEYMAP:
            _set_decimation(transfer, parts["057"])
        if isinstance(transfer, rspm.StageUnsupported):
            warnings.warn(
                f"{sid}: stage {number}: {transfer.element} is not "
                "supported; the response cannot be evaluated completely",
                rspm.ResponseWarning,
                stacklevel=5,
            )
        result.append(transfer)

    if "058" in parts:
        b058 = parts["058"]
        gain = rspm.StageGain(
            sensitivity=b058.number("04"),
            frequency=b058.number("05"),
        )
        gain.stage_number = number
        result.append(gain)

    return result


def _transfer(blockette: _Blockette) -> rspm.StageResponse:
    """Convert a response blockette."""
    return _PARSERS[blockette.btype](blockette)


def _units(blockette: _Blockette, first: str) -> dict:
    """Return the input and output units starting at a field."""
    second = f"{int(first) + 1:02d}"
    return {
        "input_units": blockette.units(first),
        "output_units": blockette.units(second),
    }


def _column(rows: list, index: int, ncols: int) -> list:
    """
    Return one column of repeated fields.

    Rows may start with an index (written by most tools) or not:
    ``ncols`` is the number of values without the index.
    """
    values = []
    for row in rows:
        offset = len(row) - ncols
        if offset < 0:
            raise ValueError(
                f"repeated field with {len(row)} values, expected {ncols}"
            )
        values.append(row[offset + index])
    return values


def _complex_rows(rows: list) -> list:
    """Return complex roots from rows ``[i,] re, im, re_err, im_err``."""
    real = _column(rows, 0, 4)
    imag = _column(rows, 1, 4)
    return [complex(r, i) for r, i in zip(real, imag)]


def _parse_paz(blockette: _Blockette) -> rspm.StageResponse:
    """Convert blockette 53 (poles and zeros)."""
    code = blockette.code("03")
    units = _units(blockette, "05")
    if code not in _PZ_TYPES:
        return rspm.StageUnsupported(
            element=f"Poles & Zeros of type {code or 'unknown'!r}",
            **units,
        )
    return rspm.StagePoleZero(
        pz_type=_PZ_TYPES[code],
        normalization_factor=blockette.number("07"),
        normalization_frequency=blockette.number("08"),
        zeros=_complex_rows(blockette.rows.get("10-13", [])),
        poles=_complex_rows(blockette.rows.get("15-18", [])),
        **units,
    )


def _parse_coefficients(blockette: _Blockette) -> rspm.StageResponse:
    """
    Convert blockette 54 (coefficients).

    Without numerators the stage has a unit response (e.g. an
    analog-to-digital converter).
    """
    code = blockette.code("03")
    units = _units(blockette, "05")
    if code not in _CF_TYPES:
        return rspm.StageUnsupported(
            element=f"Coefficients of type {code or 'unknown'!r}",
            **units,
        )
    numerator = _column(blockette.rows.get("08-09", []), 0, 2)
    denominator = _column(blockette.rows.get("11-12", []), 0, 2)
    return rspm.StageCoefficients(
        cf_type=_CF_TYPES[code],
        numerator=numerator or [1.0],
        denominator=denominator or [1.0],
        **units,
    )


def _parse_response_list(blockette: _Blockette) -> rspm.StageResponse:
    """Convert blockette 55 (frequency, amplitude, phase)."""
    rows = blockette.rows.get("07-11", [])
    return rspm.StageResponseList(
        frequency=_column(rows, 0, 5),
        amplitude=_column(rows, 1, 5),
        phase=_column(rows, 3, 5),
        **_units(blockette, "04"),
    )


def _parse_fir(blockette: _Blockette) -> rspm.StageResponse:
    """
    Convert blockette 61 (FIR), expanding the symmetry.

    Raises
    ------
    ValueError
        If the symmetry code is unknown.
    """
    code = blockette.code("05") or "A"
    if code not in _FIR_SYMMETRY:
        raise ValueError(f"Unknown FIR symmetry code: {code!r}")
    coefficients = _column(blockette.rows.get("09", []), 0, 1)
    coefficients = rspm.expand_fir_symmetry(
        coefficients, _FIR_SYMMETRY[code]
    )
    return rspm.StageCoefficients(
        cf_type="DIGITAL",
        numerator=coefficients or [1.0],
        denominator=[1.0],
        description=blockette.value("04"),
        **_units(blockette, "06"),
    )


def _polynomial(blockette: _Blockette) -> rspm.StagePolynomial:
    """Convert blockette 62 (polynomial)."""
    lower = blockette.number("09")
    upper = blockette.number("10")
    if blockette.code("08") == "A":
        # Valid frequency bounds given in rad/s.
        lower = None if lower is None else lower / (2.0 * math.pi)
        upper = None if upper is None else upper / (2.0 * math.pi)
    code = blockette.code("07")
    return rspm.StagePolynomial(
        approximation_type=_POLYNOMIAL_TYPES.get(code, code or None),
        frequency_lower_bound=lower,
        frequency_upper_bound=upper,
        approximation_lower_bound=blockette.number("11"),
        approximation_upper_bound=blockette.number("12"),
        maximum_error=blockette.number("13"),
        coefficients=_column(blockette.rows.get("15-16", []), 0, 2),
        **_units(blockette, "05"),
    )


def _parse_polynomial(blockette: _Blockette) -> rspm.StageResponse:
    """Convert blockette 62 within the cascade."""
    code = blockette.code("03")
    if code and code != "P":
        return rspm.StageUnsupported(
            element=f"Polynomial of type {code!r}",
            **_units(blockette, "05"),
        )
    return _polynomial(blockette)


def _parse_unsupported(blockette: _Blockette) -> rspm.StageResponse:
    """Keep blockettes 56 and 60 as unsupported stages."""
    names = {"056": "Generic Response", "060": "Response Reference"}
    units = _units(blockette, "04") if blockette.btype == "056" else {}
    return rspm.StageUnsupported(
        element=f"{names[blockette.btype]} (blockette {blockette.btype})",
        **units,
    )


def _set_decimation(stage: rspm.StageResponse, blockette: _Blockette) -> None:
    """Copy the fields of blockette 57 to a stage."""
    stage.input_sample_rate = blockette.number("04")
    stage.decimation_factor = blockette.integer("05")
    stage.decimation_offset = blockette.integer("06")
    stage.delay = blockette.number("07")
    stage.correction = blockette.number("08")


def _output_sample_rate(stages: list) -> Optional[float]:
    """Return the output rate of the last decimation stage (evalresp)."""
    rate = None
    for stage in stages:
        fs = getattr(stage, "input_sample_rate", None)
        if fs:
            rate = fs / (getattr(stage, "decimation_factor", None) or 1)
    return rate


_PARSERS = {
    "053": _parse_paz,
    "054": _parse_coefficients,
    "055": _parse_response_list,
    "056": _parse_unsupported,
    "060": _parse_unsupported,
    "061": _parse_fir,
    "062": _parse_polynomial,
}
