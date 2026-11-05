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
SAC binary reader and writer.

SAC (Seismic Analysis Code) files hold one waveform: a 632-byte header,
the samples as 32-bit floats and, in version 7 of the header, a footer
of 22 double-precision values (SAC manual, "SAC data file format").

Reading
-------
- Byte order detected from the header version, or imposed.
- Header versions 6 and 7 (the double-precision footer takes
  precedence over the corresponding header values).
- The start time is the reference time (``NZYEAR`` ... ``NZMSEC``)
  plus ``B``.
- Coordinates, burial depth and orientation (``STLA``, ``STLO``,
  ``STEL``, ``STDP``, ``CMPAZ``, ``CMPINC``) become the channel
  metadata of the record (:attr:`Record.head.meta
  <shakelab.signals.header.Header.meta>`); the other defined header
  fields (event, picks, user fields, ...) are kept in
  ``record.head.extra["sac"]``, with the time markers (``O``, ``A``,
  ``F``, ``T0`` to ``T9``) as absolute dates.
- Files that are not evenly sampled time series (spectra, x-y data)
  are read by :class:`Sac` only.

Writing
-------
- Header version 6, big-endian unless requested otherwise.
- The reference time is the start time truncated to the millisecond,
  and ``B`` holds the rest, so the start time keeps its microseconds.
- Channel metadata and ``head.extra["sac"]`` are written back, so a
  file read and written again keeps its picks and event information.
- Samples are stored as 32-bit floats: float64 data lose precision.
"""

from __future__ import annotations

import os
import struct
import warnings
from datetime import date as _date
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import Record

__all__ = [
    "HEADER_FIELDS",
    "SAC_UNDEFINED",
    "Sac",
    "SacWarning",
    "is_sac",
    "sacread",
    "sacwrite",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

SAC_UNDEFINED = -12345
DEFAULT_BYTE_ORDER = "be"
HEADER_SIZE = 632

_FLOATS = (
    "DELTA", "DEPMIN", "DEPMAX", "SCALE", "ODELTA", "B", "E", "O", "A",
    "INTERNAL1", "T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8",
    "T9", "F", "RESP0", "RESP1", "RESP2", "RESP3", "RESP4", "RESP5",
    "RESP6", "RESP7", "RESP8", "RESP9", "STLA", "STLO", "STEL", "STDP",
    "EVLA", "EVLO", "EVEL", "EVDP", "MAG", "USER0", "USER1", "USER2",
    "USER3", "USER4", "USER5", "USER6", "USER7", "USER8", "USER9",
    "DIST", "AZ", "BAZ", "GCARC", "SB", "SDELTA", "DEPMEN", "CMPAZ",
    "CMPINC", "XMINIMUM", "XMAXIMUM", "YMINIMUM", "YMAXIMUM", "UNUSED1",
    "UNUSED2", "UNUSED3", "UNUSED4", "UNUSED5", "UNUSED6", "UNUSED7",
)
_INTEGERS = (
    "NZYEAR", "NZJDAY", "NZHOUR", "NZMIN", "NZSEC", "NZMSEC", "NVHDR",
    "NORID", "NEVID", "NPTS", "NSNPTS", "NWFID", "NXSIZE", "NYSIZE",
    "UNUSED8", "IFTYPE", "IDEP", "IZTYPE", "UNUSED9", "IINST", "ISTREG",
    "IEVREG", "IEVTYP", "IQUAL", "ISYNTH", "IMAGTYP", "IMAGSRC",
    "IBODY", "UNUSED10", "UNUSED11", "UNUSED12", "UNUSED13", "UNUSED14",
    "UNUSED15", "UNUSED16",
)
_LOGICALS = ("LEVEN", "LPSPOL", "LOVROK", "LCALDA", "UNUSED17")
_STRINGS = (
    "KSTNM", "KEVNM", "KHOLE", "KO", "KA", "KT0", "KT1", "KT2", "KT3",
    "KT4", "KT5", "KT6", "KT7", "KT8", "KT9", "KF", "KUSER0", "KUSER1",
    "KUSER2", "KCMPNM", "KNETWK", "KDATRD", "KINST",
)

HEADER_FIELDS = _FLOATS + _INTEGERS + _LOGICALS + _STRINGS
_HEADER_FORMAT = "70f40i192s"

# Fields of the version-7 footer (double precision), in file order.
FOOTER_FIELDS = (
    "DELTA", "B", "E", "O", "A", "T0", "T1", "T2", "T3", "T4", "T5", "T6",
    "T7", "T8", "T9", "F", "EVLO", "EVLA", "STLO", "STLA", "SB", "SDELTA",
)

# Enumerated values used here.
ITIME = 1
IB = 9

# Time markers relative to the reference time.
TIME_MARKERS = ("O", "A", "F") + tuple(f"T{k}" for k in range(10))

# Fields set from the record (not kept in ``head.extra["sac"]``).
_RECORD_FIELDS = {
    "DELTA", "B", "E", "NPTS", "NZYEAR", "NZJDAY", "NZHOUR", "NZMIN",
    "NZSEC", "NZMSEC", "NVHDR", "IFTYPE", "LEVEN", "IZTYPE", "DEPMIN",
    "DEPMAX", "DEPMEN", "KNETWK", "KSTNM", "KHOLE", "KCMPNM", "STLA",
    "STLO", "STEL", "STDP", "CMPAZ", "CMPINC",
}

_UNIX_DAY0 = _date(1970, 1, 1).toordinal()


class SacWarning(UserWarning):
    """Warning about SAC header values that are ignored."""


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def sacread(source: Any, byte_order: Optional[str] = None) -> Record:
    """
    Read a SAC file as a ShakeLab Record.

    Parameters
    ----------
    source : str, path-like, bytes or binary file object
        SAC file, its content, or a file opened in binary mode.
    byte_order : {"be", "le"}, optional
        Byte order; detected from the header if None.

    Returns
    -------
    Record
        Waveform (float32 samples) with its channel metadata and, in
        ``head.extra["sac"]``, the other defined header fields.

    Raises
    ------
    ValueError
        If the file is invalid or is not an evenly sampled time series.
    """
    return Sac.read(source, byte_order).to_shakelab()


def sacwrite(
    target: Any,
    record: Record,
    byte_order: str = DEFAULT_BYTE_ORDER,
    owrite: bool = False,
) -> None:
    """
    Write a ShakeLab Record to a SAC file.

    Parameters
    ----------
    target : str, path-like or binary file object
        Output file, or a file opened in binary mode.
    record : Record
        Waveform to write.
    byte_order : {"be", "le"}, optional
        Byte order. Default is ``"be"``.
    owrite : bool, optional
        If True, overwrite an existing file.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    """
    Sac.from_shakelab(record).write(target, byte_order, owrite)


def is_sac(head: bytes) -> bool:
    """
    Tell whether the beginning of a file is a binary SAC header.

    Parameters
    ----------
    head : bytes
        First bytes of the file (at least the 632 of the header).

    Returns
    -------
    bool
        True if, in one byte order, the header version (``NVHDR``) is
        6 or 7, the number of points is not negative, the file type
        (``IFTYPE``) is known and ``LEVEN`` is 0 or 1.

    Examples
    --------
    >>> with open("record.sac", "rb") as fid:
    ...     is_sac(fid.read(1024))
    True
    """
    if not isinstance(head, (bytes, bytearray, memoryview)):
        raise TypeError("head must be bytes")
    if len(head) < HEADER_SIZE:
        return False
    buffer = memoryview(bytes(head))

    def integer(name: str, prefix: str) -> int:
        offset = 4 * (70 + _INTEGERS.index(name))
        return struct.unpack_from(prefix + "i", buffer, offset)[0]

    leven = 4 * (70 + len(_INTEGERS))
    for prefix in (">", "<"):
        if (
            integer("NVHDR", prefix) in (6, 7)
            and integer("NPTS", prefix) >= 0
            and integer("IFTYPE", prefix) in (1, 2, 3, 4, 51)
            and struct.unpack_from(prefix + "i", buffer, leven)[0]
            in (0, 1)
        ):
            return True
    return False


# ----------------------------------------------------------------------------
# SAC container
# ----------------------------------------------------------------------------


class Sac:
    """
    Content of a SAC file.

    Parameters
    ----------
    **fields
        Initial header values (names as in :data:`HEADER_FIELDS`).

    Attributes
    ----------
    head : dict
        Header fields by name; undefined fields are None. Logical fields
        are booleans, text fields stripped strings.
    data : numpy.ndarray
        First data block (float32): the samples of a time series.
    data2 : numpy.ndarray or None
        Second data block (imaginary part or phase of a spectrum,
        abscissas of x-y data), None if absent.
    byte_order : str
        Byte order of the file read (``"be"`` or ``"le"``).
    """

    def __init__(self, **fields: Any) -> None:
        self.head: dict = {name: None for name in HEADER_FIELDS}
        self.head.update(NVHDR=6, IFTYPE=ITIME, LEVEN=True)
        self.data = np.array([], dtype=np.float32)
        self.data2: Optional[np.ndarray] = None
        self.byte_order = DEFAULT_BYTE_ORDER

        for key, value in fields.items():
            if key not in self.head:
                raise KeyError(f"unknown SAC header field: {key}")
            self.head[key] = value

    def __repr__(self) -> str:
        return (
            f"Sac({self.sid}, npts={self.head['NPTS']}, "
            f"delta={self.head['DELTA']}, iftype={self.head['IFTYPE']})"
        )

    # -- properties -----------------------------------------------------

    @property
    def is_timeseries(self) -> bool:
        """True for evenly sampled time series."""
        return self.head["IFTYPE"] == ITIME and bool(self.head["LEVEN"])

    @property
    def sid(self) -> Optional[str]:
        """Source identifier ``NET.STA.LOC.CHA`` (None if no codes)."""
        codes = [
            self.head[name] or ""
            for name in ("KNETWK", "KSTNM", "KHOLE", "KCMPNM")
        ]
        return ".".join(codes) if any(codes) else None

    @property
    def reference_ns(self) -> Optional[int]:
        """Reference time in nanoseconds from the Unix epoch."""
        h = self.head
        fields = ("NZYEAR", "NZJDAY", "NZHOUR", "NZMIN", "NZSEC", "NZMSEC")
        if any(h[name] is None for name in fields):
            return None
        days = _date(h["NZYEAR"], 1, 1).toordinal() - _UNIX_DAY0
        seconds = (
            (days + h["NZJDAY"] - 1) * 86400
            + h["NZHOUR"] * 3600 + h["NZMIN"] * 60 + h["NZSEC"]
        )
        return seconds * 10**9 + h["NZMSEC"] * 10**6

    @property
    def start_ns(self) -> Optional[int]:
        """Time of the first sample (reference time + ``B``)."""
        reference = self.reference_ns
        if reference is None:
            return None
        return reference + round((self.head["B"] or 0.0) * 1e9)

    # -- conversion -----------------------------------------------------

    @classmethod
    def read(cls, source: Any, byte_order: Optional[str] = None) -> Sac:
        """
        Read a SAC file.

        Parameters
        ----------
        source : str, path-like, bytes or binary file object
            SAC file, its content, or a file opened in binary mode.
        byte_order : {"be", "le"}, optional
            Byte order; detected from the header if None.

        Raises
        ------
        ValueError
            If the file is invalid.
        """
        return cls.from_bytes(_read_source(source), byte_order)

    @classmethod
    def from_bytes(cls, buffer: Any, byte_order: Optional[str] = None) -> Sac:
        """
        Decode a SAC file.

        Parameters
        ----------
        buffer : bytes-like
            Content of the file.
        byte_order : {"be", "le"}, optional
            Byte order; detected from the header version if None.

        Raises
        ------
        ValueError
            If the header or the length of the file is invalid.
        """
        buffer = memoryview(buffer)
        if len(buffer) < HEADER_SIZE:
            raise ValueError("file shorter than a SAC header")
        if byte_order is None:
            byte_order = _detect_byte_order(buffer)
        prefix = _prefix(byte_order)

        sac = cls()
        sac.byte_order = byte_order
        values = struct.unpack_from(prefix + _HEADER_FORMAT, buffer)
        numbers = values[:110]
        text = values[110]

        head = {}
        for name, value in zip(_FLOATS, numbers[:70]):
            head[name] = None if value == SAC_UNDEFINED else float(value)
        for name, value in zip(_INTEGERS, numbers[70:105]):
            head[name] = None if value == SAC_UNDEFINED else int(value)
        for name, value in zip(_LOGICALS, numbers[105:110]):
            head[name] = None if value == SAC_UNDEFINED else bool(value)
        position = 0
        for name in _STRINGS:
            width = 16 if name == "KEVNM" else 8
            head[name] = _decode_text(text[position:position + width])
            position += width
        sac.head = head

        version = head["NVHDR"]
        npts = head["NPTS"]
        if version not in (6, 7):
            raise ValueError(f"unsupported SAC header version {version}")
        if npts is None or npts < 0:
            raise ValueError(f"invalid number of samples {npts}")

        blocks = 1 if sac.is_timeseries else 2
        size = HEADER_SIZE + 4 * npts * blocks
        footer = 8 * len(FOOTER_FIELDS) if version == 7 else 0
        if len(buffer) < size + footer:
            raise ValueError("SAC file shorter than declared")

        dtype = np.dtype(prefix + "f4")
        sac.data = np.frombuffer(
            buffer, dtype=dtype, count=npts, offset=HEADER_SIZE
        ).astype(np.float32)
        if blocks == 2:
            sac.data2 = np.frombuffer(
                buffer, dtype=dtype, count=npts,
                offset=HEADER_SIZE + 4 * npts,
            ).astype(np.float32)

        if version == 7:
            doubles = struct.unpack_from(
                prefix + f"{len(FOOTER_FIELDS)}d", buffer, size
            )
            for name, value in zip(FOOTER_FIELDS, doubles):
                head[name] = None if value == SAC_UNDEFINED else value

        if head["DELTA"] is None or head["DELTA"] <= 0:
            if sac.is_timeseries:
                raise ValueError(f"invalid sampling interval {head['DELTA']}")
        return sac

    def to_bytes(self, byte_order: str = DEFAULT_BYTE_ORDER) -> bytes:
        """
        Encode the file with header version 6.

        Parameters
        ----------
        byte_order : {"be", "le"}, optional
            Byte order. Default is ``"be"``.
        """
        prefix = _prefix(byte_order)
        head = dict(self.head)
        data = np.asarray(self.data, dtype=np.float32)
        head["NVHDR"] = 6
        head["NPTS"] = int(data.size)
        if data.size:
            head["DEPMIN"] = float(data.min())
            head["DEPMAX"] = float(data.max())
            head["DEPMEN"] = float(data.mean())
        if head["DELTA"] is not None and head["B"] is not None and data.size:
            head["E"] = head["B"] + (data.size - 1) * head["DELTA"]

        numbers = []
        for name in _FLOATS:
            value = head[name]
            numbers.append(SAC_UNDEFINED if value is None else float(value))
        for name in _INTEGERS + _LOGICALS:
            value = head[name]
            numbers.append(SAC_UNDEFINED if value is None else int(value))
        text = b"".join(
            _encode_text(head[name], 16 if name == "KEVNM" else 8)
            for name in _STRINGS
        )

        content = struct.pack(prefix + _HEADER_FORMAT, *numbers, text)
        content += data.astype(prefix + "f4").tobytes()
        if self.data2 is not None:
            content += np.asarray(self.data2).astype(prefix + "f4").tobytes()
        return content

    def write(
        self,
        target: Any,
        byte_order: str = DEFAULT_BYTE_ORDER,
        owrite: bool = False,
    ) -> None:
        """
        Write the SAC file.

        Parameters
        ----------
        target : str, path-like or binary file object
            Output file, or a file opened in binary mode.
        byte_order : {"be", "le"}, optional
            Byte order. Default is ``"be"``.
        owrite : bool, optional
            If True, overwrite an existing file.

        Raises
        ------
        FileExistsError
            If the file exists and ``owrite`` is False.
        """
        content = self.to_bytes(byte_order)
        if isinstance(target, (str, os.PathLike)):
            if os.path.exists(target) and not owrite:
                raise FileExistsError(f"File exists: {target}")
            with open(target, "wb") as fid:
                fid.write(content)
        elif hasattr(target, "write"):
            target.write(content)
        else:
            raise TypeError("target must be a path or a binary file object")

    def to_shakelab(self) -> Record:
        """
        Convert an evenly sampled time series to a ShakeLab Record.

        Raises
        ------
        ValueError
            If the file is not an evenly sampled time series.
        """
        from shakelab.signals.base import Record
        from shakelab.signals.metadata import ChannelEpoch

        if not self.is_timeseries:
            raise ValueError(
                "the SAC file is not an evenly sampled time series "
                f"(IFTYPE={self.head['IFTYPE']}, LEVEN={self.head['LEVEN']});"
                " read it with Sac"
            )

        start = self.start_ns
        record = Record(
            data=self.data.copy(),
            delta=self.head["DELTA"],
            time=None if start is None else Date.from_unix_ns(start),
        )
        record.head.meta = ChannelEpoch(self.sid, **self._channel())
        record.head.extra["sac"] = self._extra()
        return record

    @classmethod
    def from_shakelab(cls, record: Record) -> Sac:
        """
        Build the SAC content of a ShakeLab Record.

        The channel metadata of the record and the fields in
        ``record.head.extra["sac"]`` are written to the header.
        """
        sac = cls()
        head = sac.head
        extra = dict(record.head.extra.get("sac", {}))

        start = (
            None if record.head.time is None
            else record.head.time.to_unix_ns()
        )
        reference = None
        if start is not None:
            reference = start - start % 10**6
            head.update(_time_fields(reference))
            head["B"] = (start - reference) / 1e9
            head["IZTYPE"] = IB
        else:
            head["B"] = 0.0

        for name, value in extra.items():
            if name not in head or name in _RECORD_FIELDS:
                continue
            if name in TIME_MARKERS and isinstance(value, Date):
                if reference is None:
                    continue
                value = (value.to_unix_ns() - reference) / 1e9
            head[name] = value

        head["DELTA"] = record.head.delta
        sid = record.head.sid
        if sid:
            codes = sid.split(".")
            for name, code in zip(
                ("KNETWK", "KSTNM", "KHOLE", "KCMPNM"), codes
            ):
                head[name] = code or None

        meta = record.head.meta
        if meta.location is not None:
            head["STLO"] = meta.location.longitude
            head["STLA"] = meta.location.latitude
            head["STEL"] = meta.location.elevation
        if meta.depth is not None:
            head["STDP"] = meta.depth * 1000.0
        if meta.orientation is not None:
            head["CMPAZ"] = meta.orientation.azimuth
            head["CMPINC"] = meta.orientation.dip + 90.0

        sac.data = np.asarray(record.data, dtype=np.float32)
        return sac

    def info(self) -> None:
        """Print the defined header fields."""
        print(repr(self))
        for name in HEADER_FIELDS:
            value = self.head[name]
            if value is not None and not name.startswith(("UNUSED", "INT")):
                print(f"{name:>12} = {value}")

    # -- helpers --------------------------------------------------------

    def _channel(self) -> dict:
        """Channel metadata from the header."""
        h = self.head
        channel = {}
        if h["STLO"] is not None and h["STLA"] is not None:
            channel["location"] = (
                (h["STLO"], h["STLA"])
                if h["STEL"] is None
                else (h["STLO"], h["STLA"], h["STEL"])
            )
        if h["STDP"] is not None:
            channel["depth"] = h["STDP"] / 1000.0
        if h["CMPAZ"] is not None and h["CMPINC"] is not None:
            if 0.0 <= h["CMPINC"] <= 180.0:
                channel["orientation"] = (h["CMPAZ"], h["CMPINC"] - 90.0)
            else:
                warnings.warn(
                    f"CMPINC {h['CMPINC']} out of [0, 180]: orientation "
                    "ignored",
                    SacWarning,
                    stacklevel=4,
                )
        return channel

    def _extra(self) -> dict:
        """Defined header fields not represented by the record."""
        reference = self.reference_ns
        extra = {}
        for name in HEADER_FIELDS:
            value = self.head[name]
            if (
                value is None
                or name in _RECORD_FIELDS
                or name.startswith(("UNUSED", "INTERNAL"))
                or name in ("NSNPTS", "SB", "SDELTA")
            ):
                continue
            if name in TIME_MARKERS and reference is not None:
                value = Date.from_unix_ns(reference + round(value * 1e9))
            extra[name] = value
        return extra


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _read_source(source: Any) -> bytes:
    """Return the content of a path, bytes-like or binary file."""
    if isinstance(source, (bytes, bytearray, memoryview)):
        return source
    if isinstance(source, (str, os.PathLike)):
        with open(source, "rb") as fid:
            return fid.read()
    if hasattr(source, "read"):
        content = source.read()
        if not isinstance(content, (bytes, bytearray)):
            raise TypeError("the file must be opened in binary mode")
        return content
    raise TypeError("source must be a path, bytes or a binary file object")


def _prefix(byte_order: str) -> str:
    """Return the struct/NumPy prefix of a byte order."""
    if byte_order == "be":
        return ">"
    if byte_order == "le":
        return "<"
    raise ValueError(f"byte order must be 'be' or 'le', not {byte_order!r}")


def _detect_byte_order(buffer: memoryview) -> str:
    """Byte order in which the header version (NVHDR) is 6 or 7."""
    offset = 4 * (70 + _INTEGERS.index("NVHDR"))
    for byte_order in ("be", "le"):
        (version,) = struct.unpack_from(_prefix(byte_order) + "i",
                                        buffer, offset)
        if version in (6, 7):
            return byte_order
    raise ValueError("not a binary SAC file (header version not 6 or 7)")


def _decode_text(raw: bytes) -> Optional[str]:
    """
    Decode a text field.

    Blank fields, the undefined value (also repeated, as in the 16-byte
    ``KEVNM``) and fields that are not ASCII (e.g. the integer -12345
    written in place of text) are undefined (None).
    """
    raw = raw.split(b"\x00", 1)[0]
    if not raw.isascii():
        return None
    text = raw.decode("ascii").strip()
    if not text or set(text.split()) == {str(SAC_UNDEFINED)}:
        return None
    return text


def _encode_text(value: Optional[str], width: int) -> bytes:
    """Encode a text field, padded with spaces."""
    text = str(SAC_UNDEFINED) if value is None else str(value)
    if not text.isascii():
        raise ValueError(f"SAC text fields must be ASCII: {text!r}")
    if len(text) > width:
        raise ValueError(f"{text!r} exceeds the SAC field width {width}")
    return text.ljust(width).encode("ascii")


def _time_fields(unix_ns: int) -> dict:
    """Reference-time fields of an instant (to the millisecond)."""
    days, rest = divmod(unix_ns, 86_400 * 10**9)
    day = _date.fromordinal(days + _UNIX_DAY0)
    seconds, ns = divmod(rest, 10**9)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return {
        "NZYEAR": day.year,
        "NZJDAY": day.timetuple().tm_yday,
        "NZHOUR": hours,
        "NZMIN": minutes,
        "NZSEC": seconds,
        "NZMSEC": ns // 10**6,
    }
