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
DYNA 1.2 ASCII reader and writer.

DYNA 1.2 is the ASCII format of the Italian Accelerometric Archive
(ITACA) and of the Engineering Strong Motion database (ESM): a header
of 64 ``KEY: value`` lines followed by the data, one sample per line
for waveforms (``ACC``, ``VEL``, ``DIS``) and two columns for spectra
(period and ordinate for the response spectra ``SA``, ``SD``,
frequency and amplitude for the Fourier amplitude spectra ``FAS``).
Times are UTC (GMT in the ITACA manual); amplitudes are in cm/s**2,
cm/s or cm, as declared by ``UNITS``.

Reading
-------
- The header ends at the first line that is not ``KEY: value``; files
  of ESM 3.0, which write ``DATABASE_VERSION`` and ``HEADER_FORMAT`` on
  one line, are read correctly.
- Empty values and ``None`` become None; numeric fields are converted.
- The start time keeps its milliseconds; the units are declared in the
  record header (``CM/S**2``, ``CM/S``, ``CM``), the values are not
  converted.
- Station coordinates and sensor depth become the channel metadata of
  the record (:attr:`Record.head.meta
  <shakelab.signals.header.Header.meta>`), ``EVENT_ID`` the event
  identifier; the other fields (event, site, instrument, processing,
  data citation, ...) are kept in ``record.head.extra["dyna"]``.
- Spectra are read by :class:`Dyna` only.

Writing
-------
- The 64 header lines in the order of ESM, with the peak value of the
  quantity written (signed, as in ESM) and its time from the first
  sample, the number of samples and the duration computed from the
  data.
- Amplitudes in m/s**2, m/s or m are converted to cm/s**2, cm/s or cm.
- The start time is written to the millisecond.
- ``head.extra["dyna"]`` is written back, so a file read and written
  again keeps its event, site and processing information.
"""

from __future__ import annotations

import re
import warnings
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from shakelab.libutils.ascii import read_text, write_text
from shakelab.libutils.timeN import Date

if TYPE_CHECKING:
    from shakelab.signals.base import Record

__all__ = [
    "DEFAULT_FORMAT",
    "HEADER_FIELDS",
    "Dyna",
    "DynaWarning",
    "dynaread",
    "dynawrite",
    "is_dyna",
]


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

HEADER_FORMAT = "DYNA 1.2"
DEFAULT_FORMAT = "%.9g"

# Header fields in file order (ESM 3.0). PEAK stands for the pair of
# fields of the peak value of the quantity in the file.
_LAYOUT = (
    "EVENT_NAME", "EVENT_ID", "EVENT_DATE_YYYYMMDD", "EVENT_TIME_HHMMSS",
    "EVENT_LATITUDE_DEGREE", "EVENT_LONGITUDE_DEGREE", "EVENT_DEPTH_KM",
    "HYPOCENTER_REFERENCE", "MAGNITUDE_W", "MAGNITUDE_W_REFERENCE",
    "MAGNITUDE_L", "MAGNITUDE_L_REFERENCE", "FOCAL_MECHANISM",
    "NETWORK", "STATION_CODE", "STATION_NAME", "STATION_LATITUDE_DEGREE",
    "STATION_LONGITUDE_DEGREE", "STATION_ELEVATION_M", "LOCATION",
    "SENSOR_DEPTH_M", "VS30_M/S", "SITE_CLASSIFICATION_EC8",
    "MORPHOLOGIC_CLASSIFICATION", "EPICENTRAL_DISTANCE_KM",
    "EARTHQUAKE_BACKAZIMUTH_DEGREE",
    "DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS",
    "DATE_TIME_FIRST_SAMPLE_PRECISION", "SAMPLING_INTERVAL_S", "NDATA",
    "DURATION_S", "STREAM", "UNITS", "INSTRUMENT",
    "INSTRUMENT_ANALOG/DIGITAL", "INSTRUMENTAL_FREQUENCY_HZ",
    "INSTRUMENTAL_DAMPING", "FULL_SCALE_G", "N_BIT_DIGITAL_CONVERTER",
    "PEAK", "BASELINE_CORRECTION", "FILTER_TYPE", "FILTER_ORDER",
    "LOW_CUT_FREQUENCY_HZ", "HIGH_CUT_FREQUENCY_HZ",
    "LATE/NORMAL_TRIGGERED", "DATABASE_VERSION", "HEADER_FORMAT",
    "DATA_TYPE", "PROCESSING", "DATA_TIMESTAMP_YYYYMMDD_HHMMSS",
    "DATA_LICENSE", "DATA_CITATION", "DATA_CREATOR",
    "ORIGINAL_DATA_MEDIATOR_CITATION", "ORIGINAL_DATA_MEDIATOR",
    "ORIGINAL_DATA_CREATOR_CITATION", "ORIGINAL_DATA_CREATOR",
    "USER1", "USER2", "USER3", "USER4", "USER5",
)

# Quantities: file units, record units, data type, peak fields.
_QUANTITIES = {
    "acceleration": (
        "cm/s^2", "CM/S**2", "ACCELERATION", "PGA_CM/S^2", "TIME_PGA_S",
    ),
    "velocity": ("cm/s", "CM/S", "VELOCITY", "PGV_CM/S", "TIME_PGV_S"),
    "displacement": ("cm", "CM", "DISPLACEMENT", "PGD_CM", "TIME_PGD_S"),
}
_PEAK_FIELDS = tuple(
    name for q in _QUANTITIES.values() for name in q[3:]
)

#: All the known header fields, in file order (peaks of the three
#: quantities in place of the single pair of a file).
HEADER_FIELDS = tuple(
    name
    for field in _LAYOUT
    for name in (_PEAK_FIELDS if field == "PEAK" else (field,))
)

_FLOAT_FIELDS = {
    "EVENT_LATITUDE_DEGREE", "EVENT_LONGITUDE_DEGREE", "EVENT_DEPTH_KM",
    "MAGNITUDE_W", "MAGNITUDE_L", "STATION_LATITUDE_DEGREE",
    "STATION_LONGITUDE_DEGREE", "STATION_ELEVATION_M", "SENSOR_DEPTH_M",
    "VS30_M/S", "EPICENTRAL_DISTANCE_KM", "EARTHQUAKE_BACKAZIMUTH_DEGREE",
    "SAMPLING_INTERVAL_S", "DURATION_S", "INSTRUMENTAL_FREQUENCY_HZ",
    "INSTRUMENTAL_DAMPING", "FULL_SCALE_G", "LOW_CUT_FREQUENCY_HZ",
    "HIGH_CUT_FREQUENCY_HZ",
} | set(_PEAK_FIELDS)
_INT_FIELDS = {"NDATA", "FILTER_ORDER", "N_BIT_DIGITAL_CONVERTER"}

# Fields set from the record (not kept in ``head.extra["dyna"]``).
_RECORD_FIELDS = {
    "EVENT_ID", "NETWORK", "STATION_CODE", "LOCATION", "STREAM",
    "STATION_LATITUDE_DEGREE", "STATION_LONGITUDE_DEGREE",
    "STATION_ELEVATION_M", "SENSOR_DEPTH_M",
    "DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS",
    "DATE_TIME_FIRST_SAMPLE_PRECISION", "SAMPLING_INTERVAL_S", "NDATA",
    "DURATION_S", "UNITS", "DATA_TYPE", "HEADER_FORMAT",
} | set(_PEAK_FIELDS)

# SI units accepted when writing, with the factor to centimetres.
_SI_UNITS = {
    "M/S**2": ("acceleration", 100.0),
    "M/S": ("velocity", 100.0),
    "M": ("displacement", 100.0),
}

_KEY = re.compile(r"^([A-Z][A-Z0-9_/^]*):(.*)$")


class DynaWarning(UserWarning):
    """Warning about inconsistent DYNA header values."""


# ----------------------------------------------------------------------------
# High-level functions
# ----------------------------------------------------------------------------


def dynaread(source: Any) -> Record:
    """
    Read a DYNA 1.2 waveform file as a ShakeLab Record.

    Parameters
    ----------
    source : str, path-like, bytes or file object
        DYNA file, its content (a string with line breaks, or bytes),
        or an open file.

    Returns
    -------
    Record
        Waveform (float64 samples, in the units of the file) with its
        channel metadata and, in ``head.extra["dyna"]``, the other
        header fields.

    Raises
    ------
    ValueError
        If the file is invalid or holds a spectrum.
    """
    return Dyna.read(source).to_shakelab()


def dynawrite(
    target: Any,
    record: Record,
    fmt: str = DEFAULT_FORMAT,
    owrite: bool = False,
) -> None:
    """
    Write a ShakeLab Record to a DYNA 1.2 file.

    Parameters
    ----------
    target : str, path-like or text file object
        Output file, or a file opened in text mode.
    record : Record
        Waveform in cm/s**2, cm/s or cm (or in the SI units, converted
        to centimetres).
    fmt : str, optional
        ``%`` format of the samples. Default is ``"%.9g"`` (nine
        significant digits; ESM writes ``"%f"``).
    owrite : bool, optional
        If True, overwrite an existing file.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    ValueError
        If the units of the record are missing or not supported.
    """
    Dyna.from_shakelab(record).write(target, fmt=fmt, owrite=owrite)


def is_dyna(head: bytes) -> bool:
    """
    Tell whether the beginning of a file is a DYNA file.

    Parameters
    ----------
    head : bytes
        First bytes of the file.

    Returns
    -------
    bool
        True if the first line that is not blank starts with
        ``EVENT_NAME:``.

    Examples
    --------
    >>> with open("IV.MODE..HNZ.D.IT-2012-0008.ACC.MP.ASC", "rb") as fid:
    ...     is_dyna(fid.read(1024))
    True
    """
    if not isinstance(head, (bytes, bytearray, memoryview)):
        raise TypeError("head must be bytes")
    for line in bytes(head).decode("latin-1").splitlines():
        if line.strip():
            return line.lstrip().startswith("EVENT_NAME:")
    return False


# ----------------------------------------------------------------------------
# DYNA container
# ----------------------------------------------------------------------------


class Dyna:
    """
    Content of a DYNA 1.2 file.

    Parameters
    ----------
    header : dict, optional
        Header fields by name.
    data : array_like, optional
        Samples of a waveform, or ordinates of a spectrum.
    abscissa : array_like, optional
        Periods (s) of a response spectrum or frequencies (Hz) of a
        Fourier spectrum; None for a waveform.

    Attributes
    ----------
    header : dict
        Header fields in file order; empty fields are None, numeric
        fields are int or float, the others stripped strings (dates
        included, as written).
    data : numpy.ndarray
        Samples or spectral ordinates (float64).
    abscissa : numpy.ndarray or None
        Periods or frequencies of a spectrum.
    """

    def __init__(
        self,
        header: Optional[dict] = None,
        data: Any = None,
        abscissa: Any = None,
    ) -> None:
        if header is not None and not isinstance(header, dict):
            raise TypeError("header must be a dict")
        self.header: dict = dict(header or {})
        self.data = np.asarray(
            [] if data is None else data, dtype=np.float64
        )
        self.abscissa = (
            None if abscissa is None
            else np.asarray(abscissa, dtype=np.float64)
        )

    def __repr__(self) -> str:
        kind = self.header.get("DATA_TYPE") or "unknown"
        return (
            f"Dyna(sid={self.sid!r}, data_type={kind!r}, "
            f"ndata={self.data.size})"
        )

    # -- properties -----------------------------------------------------

    @property
    def sid(self) -> Optional[str]:
        """Source identifier ``NET.STA.LOC.CHA``, None if incomplete."""
        h = self.header
        net, sta, cha = (
            h.get("NETWORK"), h.get("STATION_CODE"), h.get("STREAM")
        )
        if not (net and sta and cha):
            return None
        return f"{net}.{sta}.{h.get('LOCATION') or ''}.{cha}"

    @property
    def is_timeseries(self) -> bool:
        """Whether the file holds a waveform rather than a spectrum."""
        return self.abscissa is None

    @property
    def start_ns(self) -> Optional[int]:
        """Time of the first sample, Unix ns (UTC); None if not given."""
        value = self.header.get("DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS")
        if not value:
            return None
        return _parse_time(value)

    # -- conversion -----------------------------------------------------

    @classmethod
    def read(cls, source: Any) -> Dyna:
        """
        Read a DYNA file.

        Parameters
        ----------
        source : str, path-like, bytes or file object
            DYNA file, its content (a string with line breaks, or
            bytes), or an open file.

        Raises
        ------
        ValueError
            If the header or the data are invalid.
        """
        return cls.from_text(read_text(source))

    @classmethod
    def from_text(cls, text: str) -> Dyna:
        """
        Parse the content of a DYNA file.

        Raises
        ------
        ValueError
            If the header or the data are invalid.
        """
        lines = text.splitlines()
        header: dict = {}
        k = 0
        for k, line in enumerate(lines):
            if not line.strip():
                continue
            match = _KEY.match(line.strip())
            if match is None:
                break
            _parse_field(header, match.group(1), match.group(2))
        else:
            k = len(lines)

        if "NDATA" not in header and "UNITS" not in header:
            raise ValueError("not a DYNA file: header fields not found")

        values = " ".join(lines[k:]).split()
        ncol = len(lines[k].split()) if k < len(lines) else 1
        try:
            data = np.array(values, dtype=np.float64)
        except ValueError as err:
            raise ValueError(f"invalid DYNA data: {err}") from None
        if ncol not in (1, 2) or data.size % ncol:
            raise ValueError(
                f"invalid DYNA data: {data.size} values in {ncol} columns"
            )

        dyna = cls(header)
        if ncol == 2:
            data = data.reshape(-1, 2)
            dyna.abscissa = data[:, 0].copy()
            dyna.data = data[:, 1].copy()
        else:
            dyna.data = data

        ndata = header.get("NDATA")
        if ndata is not None and ndata != dyna.data.size:
            warnings.warn(
                f"NDATA is {ndata}, {dyna.data.size} values read",
                DynaWarning,
                stacklevel=3,
            )
        return dyna

    def to_text(self, fmt: str = DEFAULT_FORMAT) -> str:
        """
        Return the content of the DYNA file.

        Parameters
        ----------
        fmt : str, optional
            ``%`` format of the data values.
        """
        lines = [
            f"{name}: {_format_value(value)}"
            for name, value in self._layout()
        ]
        if self.abscissa is None:
            lines.extend(fmt % v for v in self.data)
        else:
            lines.extend(
                f"{fmt % x} {fmt % y}"
                for x, y in zip(self.abscissa, self.data)
            )
        return "\n".join(lines) + "\n"

    def write(
        self,
        target: Any,
        fmt: str = DEFAULT_FORMAT,
        owrite: bool = False,
    ) -> None:
        """
        Write the DYNA file.

        Parameters
        ----------
        target : str, path-like or text file object
            Output file, or a file opened in text mode.
        fmt : str, optional
            ``%`` format of the data values.
        owrite : bool, optional
            If True, overwrite an existing file.

        Raises
        ------
        FileExistsError
            If the file exists and ``owrite`` is False.
        """
        write_text(target, self.to_text(fmt), owrite)

    def to_shakelab(self) -> Record:
        """
        Convert a waveform to a ShakeLab Record.

        Raises
        ------
        ValueError
            If the file holds a spectrum, or the sampling interval is
            missing.
        """
        from shakelab.signals.base import Record
        from shakelab.signals.metadata import ChannelEpoch

        if not self.is_timeseries:
            raise ValueError(
                f"the DYNA file holds a {self.header.get('DATA_TYPE')!r}"
                ", not a waveform; read it with Dyna"
            )
        h = self.header
        if h.get("SAMPLING_INTERVAL_S") is None:
            raise ValueError("SAMPLING_INTERVAL_S is missing")

        units = h.get("UNITS")
        for quantity in _QUANTITIES.values():
            if units == quantity[0]:
                units = quantity[1]
                break

        start = self.start_ns
        record = Record(
            data=self.data.copy(),
            delta=h["SAMPLING_INTERVAL_S"],
            time=None if start is None else Date.from_unix_ns(start),
            eid=h.get("EVENT_ID"),
            units=units,
        )
        record.head.meta = ChannelEpoch(self.sid, **self._channel())
        record.head.extra["dyna"] = {
            name: value for name, value in h.items()
            if name not in _RECORD_FIELDS
        }
        return record

    @classmethod
    def from_shakelab(cls, record: Record) -> Dyna:
        """
        Build the DYNA content of a ShakeLab Record.

        The channel metadata of the record and the fields in
        ``record.head.extra["dyna"]`` are written to the header.

        Raises
        ------
        ValueError
            If the units of the record are missing or not supported.
        """
        quantity, factor = _quantity(record.head.units)
        file_units, _, data_type, peak, peak_time = _QUANTITIES[quantity]

        header = {
            name: value
            for name, value in record.head.extra.get("dyna", {}).items()
            if name not in _RECORD_FIELDS
        }
        data = np.asarray(record.data, dtype=np.float64) * factor
        delta = record.head.delta

        sid = record.head.sid
        codes = sid.split(".") if sid else [None] * 4
        for name, code in zip(
            ("NETWORK", "STATION_CODE", "LOCATION", "STREAM"), codes
        ):
            header[name] = code or None
        header["EVENT_ID"] = record.head.eid

        meta = record.head.meta
        if meta.location is not None:
            header["STATION_LATITUDE_DEGREE"] = meta.location.latitude
            header["STATION_LONGITUDE_DEGREE"] = meta.location.longitude
            header["STATION_ELEVATION_M"] = meta.location.elevation
        if meta.depth is not None:
            header["SENSOR_DEPTH_M"] = meta.depth * 1000.0

        if record.head.time is not None:
            header["DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS"] = (
                _format_time(record.head.time)
            )
            header["DATE_TIME_FIRST_SAMPLE_PRECISION"] = "milliseconds"
        header["SAMPLING_INTERVAL_S"] = delta
        header["NDATA"] = data.size
        header["DURATION_S"] = data.size * delta
        header["UNITS"] = file_units
        header["DATA_TYPE"] = data_type
        header["HEADER_FORMAT"] = HEADER_FORMAT
        if data.size:
            k = int(np.argmax(np.abs(data)))
            header[peak] = float(data[k])
            header[peak_time] = k * delta
        return cls(header, data)

    def info(self) -> None:
        """Print the defined header fields."""
        print(repr(self))
        for name, value in self.header.items():
            if value is not None:
                print(f"{name:>40} = {value}")

    # -- helpers --------------------------------------------------------

    def _layout(self) -> list:
        """Header lines (name, value) in file order."""
        h = self.header
        known = set(HEADER_FIELDS)
        lines = []
        for field in _LAYOUT:
            if field == "PEAK":
                names = [n for n in _PEAK_FIELDS if n in h]
                if not names:
                    names = list(_peak_fields(h.get("UNITS")))
                lines.extend((n, h.get(n)) for n in names)
            elif field == "HEADER_FORMAT":
                lines.append((field, h.get(field) or HEADER_FORMAT))
            else:
                lines.append((field, h.get(field)))
        lines.extend(
            (name, value) for name, value in h.items() if name not in known
        )
        return lines

    def _channel(self) -> dict:
        """Channel metadata from the header."""
        h = self.header
        channel = {}
        lon = h.get("STATION_LONGITUDE_DEGREE")
        lat = h.get("STATION_LATITUDE_DEGREE")
        elev = h.get("STATION_ELEVATION_M")
        if lon is not None and lat is not None:
            channel["location"] = (
                (lon, lat) if elev is None else (lon, lat, elev)
            )
        if h.get("SENSOR_DEPTH_M") is not None:
            channel["depth"] = h["SENSOR_DEPTH_M"] / 1000.0
        return channel


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _parse_field(header: dict, name: str, value: str) -> None:
    """Store a header field, splitting fields written on one line."""
    value = value.strip()
    match = _KEY.match(value)
    while match is not None and match.group(1) in HEADER_FIELDS:
        # ESM 3.0: "DATABASE_VERSION: HEADER_FORMAT: DYNA 1.2".
        header[name] = None
        name, value = match.group(1), match.group(2).strip()
        match = _KEY.match(value)
    header[name] = _convert(name, value)


def _convert(name: str, value: str) -> Any:
    """Convert a header value to its type; None if empty."""
    if value == "" or value == "None":
        return None
    try:
        if name in _INT_FIELDS:
            return int(float(value))
        if name in _FLOAT_FIELDS:
            return float(value)
    except ValueError:
        warnings.warn(
            f"invalid value of {name}: {value!r} ignored",
            DynaWarning,
            stacklevel=5,
        )
        return None
    return value


def _format_value(value: Any) -> str:
    """Text of a header value."""
    if value is None:
        return ""
    if isinstance(value, float):
        text = f"{value:.12g}"
        if text.lstrip("-").isdigit():
            text += ".0"
        return text
    return str(value)


def _parse_time(value: str) -> int:
    """Unix ns of a ``YYYYMMDD_HHMMSS[.fff]`` string (UTC)."""
    try:
        day, clock = value.split("_")
        seconds, _, fraction = clock.partition(".")
        moment = datetime.strptime(day + seconds, "%Y%m%d%H%M%S").replace(
            tzinfo=timezone.utc
        )
    except ValueError:
        raise ValueError(f"invalid DYNA time: {value!r}") from None
    ns = int(round(float("0." + (fraction or "0")) * 1e9))
    return int(moment.timestamp()) * 10**9 + ns


def _format_time(time: Date) -> str:
    """``YYYYMMDD_HHMMSS.fff`` text of a date, rounded to the ms."""
    ms = round(time.to_unix_ns() / 1e6)
    moment = datetime.fromtimestamp(ms // 1000, tz=timezone.utc)
    return moment.strftime("%Y%m%d_%H%M%S") + f".{ms % 1000:03d}"


def _quantity(units: Optional[str]) -> tuple:
    """Quantity and factor to centimetres of the record units."""
    from shakelab.signals.response import same_units

    for quantity, fields in _QUANTITIES.items():
        if same_units(units, fields[1]) or same_units(units, fields[0]):
            return quantity, 1.0
    for name, (quantity, factor) in _SI_UNITS.items():
        if same_units(units, name):
            return quantity, factor
    raise ValueError(
        f"units {units!r} not supported by DYNA: the record must be in "
        "cm/s**2, cm/s, cm or in the corresponding SI units"
    )


def _peak_fields(units: Optional[str]) -> tuple:
    """Peak fields of the quantity of a file, from its units."""
    for fields in _QUANTITIES.values():
        if units == fields[0]:
            return fields[3:]
    return ()
