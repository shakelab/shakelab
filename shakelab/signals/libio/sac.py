# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
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
Lightweight SAC file reader and writer.

This module provides a compact SAC binary implementation used by ShakeLab.
It supports the standard SAC header layout and evenly sampled waveform data.
The public API is intentionally kept compatible with the previous module:

- sacread()
- sacwrite()
- Sac

Data blocks are read and written using NumPy arrays, avoiding slow
sample-by-sample Python loops.
"""

from os.path import isfile
from struct import pack, unpack

import numpy as np

from shakelab.libutils.time import Date


DEFAULT_BYTE_ORDER = "be"

HEADER_SIZE = 632
SAC_UNDEFINED = -12345


def sacread(input_file, byte_order=None):
    """
    Read a SAC file and return a ShakeLab Record.

    Parameters
    ----------
    input_file : str
        Input SAC file path.
    byte_order : str or None, optional
        Byte order, either ``"be"`` or ``"le"``.
        If ``None``, the byte order is detected automatically
        by trying big-endian first and then little-endian.

    Returns
    -------
    Record
        ShakeLab signal record.
    """
    from shakelab.signals.base import Record

    if byte_order is None:
        sac = None

        for order in ("be", "le"):
            try:
                sac = Sac(input_file, byte_order=order)
                break
            except ValueError:
                continue

        if sac is None:
            raise ValueError(
                "Unable to determine SAC byte order for file: {0}".format(
                    input_file
                )
            )

    else:
        sac = Sac(input_file, byte_order=byte_order)

    record = Record()
    record.head.delta = sac.delta
    record.head.time = Date(sac.time)
    record.data = np.asarray(sac.data[0]).copy()

    record.head.sid = "{0}.{1}..{2}".format(
        sac.head["KNETWK"],
        sac.head["KSTNM"],
        sac.head["KCMPNM"],
    )

    return record


def sacwrite(output_file, record, byte_order=DEFAULT_BYTE_ORDER,
             owrite=False):
    """
    Write a ShakeLab Record to a SAC file.

    Parameters
    ----------
    output_file : str
        Output SAC file path.
    record : Record
        ShakeLab signal record.
    byte_order : str, optional
        Byte order, either ``"be"`` or ``"le"``.
    owrite : bool, optional
        If ``True``, overwrite an existing file.
    """
    sac = Sac(byte_order=byte_order)

    sid = record.head.sid.split(".")
    data = np.asarray(record.data, dtype=np.float32)

    sac.head["DELTA"] = record.head.delta
    sac.head["NPTS"] = int(record.nsamp)
    sac.head["B"] = 0.0
    sac.head["E"] = (record.nsamp - 1) * record.head.delta
    sac.head["LEVEN"] = 1
    sac.head["IFTYPE"] = 1

    if data.size:
        sac.head["DEPMIN"] = float(np.min(data))
        sac.head["DEPMAX"] = float(np.max(data))
        sac.head["DEPMEN"] = float(np.mean(data))

    sac.head["NZYEAR"] = record.head.time.year
    sac.head["NZJDAY"] = record.head.time.ordinal_day
    sac.head["NZHOUR"] = record.head.time.hour
    sac.head["NZMIN"] = record.head.time.minute
    sac.head["NZSEC"] = int(record.head.time.second)
    sac.head["NZMSEC"] = int(round((record.head.time.second % 1) * 1000))

    if len(sid) > 0:
        sac.head["KNETWK"] = sid[0]

    if len(sid) > 1:
        sac.head["KSTNM"] = sid[1]

    if len(sid) > 2:
        sac.head["KHOLE"] = sid[2]

    if len(sid) > 3:
        sac.head["KCMPNM"] = sid[3]

    sac.data[0] = data
    sac.write(output_file, byte_order=byte_order, owrite=owrite)


class Sac:
    """
    SAC binary file container.

    Parameters
    ----------
    file : str, optional
        Input SAC file path.  If omitted, an empty SAC object is created.
    byte_order : str, optional
        Byte order, either ``"be"`` or ``"le"``.

    Attributes
    ----------
    head : dict
        SAC header fields.
    data : list
        Data blocks. ``data[0]`` is the primary waveform block, while
        ``data[1]`` is used for the optional second data block.
    byte : str
        Current byte order.
    """

    def __init__(self, file=None, byte_order=DEFAULT_BYTE_ORDER):
        self.head = _empty_header()
        self.data = [
            np.array([], dtype=np.float32),
            np.array([], dtype=np.float32),
        ]
        self.byte = byte_order

        if file is not None:
            self.read(file, byte_order=byte_order)

    def read(self, file, byte_order=DEFAULT_BYTE_ORDER):
        """
        Read a SAC file from disk.

        Parameters
        ----------
        file : str
            Input SAC file path.
        byte_order : str, optional
            Byte order, either ``"be"`` or ``"le"``.
        """
        self.byte = byte_order

        with open(file, "rb") as fid:
            self._read_header(fid)
            self._validate_header()

            npts = self.head["NPTS"]

            if npts < 0:
                raise ValueError("Invalid SAC NPTS value: {0}".format(npts))

            self.data[0] = _read_data_block(fid, npts, self.byte)

            if self._has_second_data_block():
                self.data[1] = _read_data_block(fid, npts, self.byte)
            else:
                self.data[1] = np.array([], dtype=np.float32)

    def write(self, file, byte_order=None, owrite=False):
        """
        Write a SAC file to disk.

        Parameters
        ----------
        file : str
            Output SAC file path.
        byte_order : str, optional
            Byte order, either ``"be"`` or ``"le"``.
        owrite : bool, optional
            If ``True``, overwrite an existing file.
        """
        if isfile(file) and not owrite:
            raise FileExistsError("File exists: {0}".format(file))

        if byte_order is not None:
            self.byte = byte_order

        self._prepare_header_for_write()

        with open(file, "wb") as fid:
            self._write_header(fid)
            _write_data_block(fid, self.data[0], self.byte)

            if len(self.data[1]):
                _write_data_block(fid, self.data[1], self.byte)

    def info(self):
        """
        Print non-default SAC header fields.
        """
        print("------------")

        for name, _, _, default in HDR_STRUCT:
            value = self.head[name]

            if value != default:
                print("{0:>12} = {1}".format(name, value))

    @property
    def delta(self):
        """Sampling interval in seconds."""
        return self.head["DELTA"]

    @property
    def time(self):
        """Reference time as a ShakeLab-compatible date string."""
        date = "{0:04d}-".format(self.head["NZYEAR"])
        date += "{0:03d}T".format(self.head["NZJDAY"])
        date += "{0:02d}:".format(self.head["NZHOUR"])
        date += "{0:02d}:".format(self.head["NZMIN"])
        date += "{0:02d}.".format(self.head["NZSEC"])
        date += "{0:04d}".format(self.head["NZMSEC"])

        return date

    def copy(self):
        """
        Return a deep copy of the SAC object.
        """
        new = Sac(byte_order=self.byte)
        new.head = self.head.copy()
        new.data = [self.data[0].copy(), self.data[1].copy()]

        return new

    def _read_header(self, fid):
        """Read SAC header fields."""
        self.head = {}

        for name, byte_num, dtype, _ in HDR_STRUCT:
            self.head[name] = _fread(fid, byte_num, dtype, self.byte)

    def _write_header(self, fid):
        """Write SAC header fields."""
        for name, byte_num, dtype, _ in HDR_STRUCT:
            _fwrite(fid, self.head[name], byte_num, dtype, self.byte)

    def _validate_header(self):
        """Validate basic SAC header consistency."""
        if self.head["NVHDR"] != 6:
            raise ValueError("Invalid SAC NVHDR value")

        if self.head["NPTS"] < 0:
            raise ValueError("Invalid SAC NPTS value")

        if self.head["DELTA"] <= 0:
            raise ValueError("Invalid SAC DELTA value")

        if self.head["LEVEN"] not in (0, 1):
            raise ValueError("Invalid SAC LEVEN value")

        if self.head["IFTYPE"] not in (1, 2, 3, 4):
            raise ValueError("Invalid SAC IFTYPE value")

    def _has_second_data_block(self):
        """Return True if the SAC file stores a second data block."""
        return self.head["LEVEN"] != 1 or self.head["IFTYPE"] in (2, 3)

    def _prepare_header_for_write(self):
        """Update basic header consistency before writing."""
        npts = len(self.data[0])
        self.head["NPTS"] = int(npts)

        if npts and self.head["DELTA"] != SAC_UNDEFINED:
            self.head["E"] = self.head["B"] + (npts - 1) * self.head["DELTA"]

        if npts:
            data = np.asarray(self.data[0], dtype=np.float32)
            self.head["DEPMIN"] = float(np.min(data))
            self.head["DEPMAX"] = float(np.max(data))
            self.head["DEPMEN"] = float(np.mean(data))


def _empty_header():
    """Return a SAC header dictionary initialized with default values."""
    return {name: default for name, _, _, default in HDR_STRUCT}


def _read_data_block(fid, npts, byte_order):
    """
    Read a SAC float32 data block.

    Parameters
    ----------
    fid : file object
        Open binary file.
    npts : int
        Number of samples.
    byte_order : str
        Byte order, either ``"be"`` or ``"le"``.

    Returns
    -------
    numpy.ndarray
        Float32 data array.
    """
    byte_count = npts * 4
    raw = fid.read(byte_count)

    if len(raw) != byte_count:
        raise EOFError("Unexpected end of SAC data block")

    dtype = np.dtype("f4").newbyteorder(_numpy_byte_order(byte_order))

    return np.frombuffer(raw, dtype=dtype, count=npts).copy()


def _write_data_block(fid, data, byte_order):
    """
    Write a SAC float32 data block.
    """
    dtype = np.dtype("f4").newbyteorder(_numpy_byte_order(byte_order))
    array = np.asarray(data, dtype=dtype)

    fid.write(array.tobytes())


def _fread(fid, byte_num, dtype, byte_order):
    """
    Read one SAC header field.
    """
    raw = fid.read(byte_num)

    if len(raw) != byte_num:
        raise EOFError("Unexpected end of SAC header")

    if dtype == "s":
        value = raw.decode("ascii", errors="replace")
        return value.split("\x00", 1)[0].strip()

    fmt = _struct_prefix(byte_order) + dtype

    return unpack(fmt, raw)[0]


def _fwrite(fid, value, byte_num, dtype, byte_order):
    """
    Write one SAC header field.
    """
    if dtype == "s":
        raw = _format_sac_string(value, byte_num)
        fid.write(raw)
        return

    fmt = _struct_prefix(byte_order) + dtype
    fid.write(pack(fmt, value))


def _format_sac_string(value, byte_num):
    """
    Format a SAC fixed-length string field.
    """
    if value is None:
        value = ""

    text = str(value)
    raw = text[:byte_num].ljust(byte_num).encode("ascii", errors="replace")

    return raw


def _struct_prefix(byte_order):
    """
    Return struct byte-order prefix.
    """
    if byte_order == "be":
        return ">"

    if byte_order == "le":
        return "<"

    raise ValueError("Unsupported byte order: {0}".format(byte_order))


def _numpy_byte_order(byte_order):
    """
    Return NumPy byte-order marker.
    """
    if byte_order == "be":
        return ">"

    if byte_order == "le":
        return "<"

    raise ValueError("Unsupported byte order: {0}".format(byte_order))


# INTERNAL: Header Structure
#   [0] - Field Name
#   [1] - Length in Bytes
#   [2] - Variable Type
#   [3] - Default Value

HDR_STRUCT = [
    ("DELTA", 4, "f", -12345),
    ("DEPMIN", 4, "f", -12345),
    ("DEPMAX", 4, "f", -12345),
    ("SCALE", 4, "f", -12345),
    ("ODELTA", 4, "f", -12345),
    ("B", 4, "f", -12345),
    ("E", 4, "f", -12345),
    ("O", 4, "f", -12345),
    ("A", 4, "f", -12345),
    ("INTERNAL1", 4, "f", -12345),
    ("T0", 4, "f", -12345),
    ("T1", 4, "f", -12345),
    ("T2", 4, "f", -12345),
    ("T3", 4, "f", -12345),
    ("T4", 4, "f", -12345),
    ("T5", 4, "f", -12345),
    ("T6", 4, "f", -12345),
    ("T7", 4, "f", -12345),
    ("T8", 4, "f", -12345),
    ("T9", 4, "f", -12345),
    ("F", 4, "f", -12345),
    ("RESP0", 4, "f", -12345),
    ("RESP1", 4, "f", -12345),
    ("RESP2", 4, "f", -12345),
    ("RESP3", 4, "f", -12345),
    ("RESP4", 4, "f", -12345),
    ("RESP5", 4, "f", -12345),
    ("RESP6", 4, "f", -12345),
    ("RESP7", 4, "f", -12345),
    ("RESP8", 4, "f", -12345),
    ("RESP9", 4, "f", -12345),
    ("STLA", 4, "f", -12345),
    ("STLO", 4, "f", -12345),
    ("STEL", 4, "f", -12345),
    ("STDP", 4, "f", -12345),
    ("EVLA", 4, "f", -12345),
    ("EVLO", 4, "f", -12345),
    ("EVEL", 4, "f", -12345),
    ("EVDP", 4, "f", -12345),
    ("MAG", 4, "f", -12345),
    ("USER0", 4, "f", -12345),
    ("USER1", 4, "f", -12345),
    ("USER2", 4, "f", -12345),
    ("USER3", 4, "f", -12345),
    ("USER4", 4, "f", -12345),
    ("USER5", 4, "f", -12345),
    ("USER6", 4, "f", -12345),
    ("USER7", 4, "f", -12345),
    ("USER8", 4, "f", -12345),
    ("USER9", 4, "f", -12345),
    ("DIST", 4, "f", -12345),
    ("AZ", 4, "f", -12345),
    ("BAZ", 4, "f", -12345),
    ("GCARC", 4, "f", -12345),
    ("INTERNAL2", 4, "f", -12345),
    ("INTERNAL3", 4, "f", -12345),
    ("DEPMEN", 4, "f", -12345),
    ("CMPAZ", 4, "f", -12345),
    ("CMPINC", 4, "f", -12345),
    ("XMINIMUM", 4, "f", -12345),
    ("XMAXIMUM", 4, "f", -12345),
    ("YMINIMUM", 4, "f", -12345),
    ("YMAXIMUM", 4, "f", -12345),
    ("UNUSED1", 4, "f", -12345),
    ("UNUSED2", 4, "f", -12345),
    ("UNUSED3", 4, "f", -12345),
    ("UNUSED4", 4, "f", -12345),
    ("UNUSED5", 4, "f", -12345),
    ("UNUSED6", 4, "f", -12345),
    ("UNUSED7", 4, "f", -12345),
    ("NZYEAR", 4, "i", -12345),
    ("NZJDAY", 4, "i", -12345),
    ("NZHOUR", 4, "i", -12345),
    ("NZMIN", 4, "i", -12345),
    ("NZSEC", 4, "i", -12345),
    ("NZMSEC", 4, "i", -12345),
    ("NVHDR", 4, "i", 6),
    ("NORID", 4, "i", -12345),
    ("NEVID", 4, "i", -12345),
    ("NPTS", 4, "i", -12345),
    ("INTERNAL4", 4, "i", -12345),
    ("NWFID", 4, "i", -12345),
    ("NXSIZE", 4, "i", -12345),
    ("NYSIZE", 4, "i", -12345),
    ("UNUSED8", 4, "i", -12345),
    ("IFTYPE", 4, "i", 1),
    ("IDEP", 4, "i", -12345),
    ("IZTYPE", 4, "i", -12345),
    ("UNUSED9", 4, "i", -12345),
    ("IINST", 4, "i", -12345),
    ("ISTREG", 4, "i", -12345),
    ("IEVREG", 4, "i", -12345),
    ("IEVTYP", 4, "i", -12345),
    ("IQUAL", 4, "i", -12345),
    ("ISYNTH", 4, "i", -12345),
    ("IMAGTYP", 4, "i", -12345),
    ("IMAGSRC", 4, "i", -12345),
    ("UNUSED10", 4, "i", -12345),
    ("UNUSED11", 4, "i", -12345),
    ("UNUSED12", 4, "i", -12345),
    ("UNUSED13", 4, "i", -12345),
    ("UNUSED14", 4, "i", -12345),
    ("UNUSED15", 4, "i", -12345),
    ("UNUSED16", 4, "i", -12345),
    ("UNUSED17", 4, "i", -12345),
    ("LEVEN", 4, "i", 1),
    ("LPSPOL", 4, "i", 0),
    ("LOVROK", 4, "i", 1),
    ("LCALDA", 4, "i", 1),
    ("UNUSED18", 4, "i", 0),
    ("KSTNM", 8, "s", "-12345"),
    ("KEVNM", 16, "s", "-12345"),
    ("KHOLE", 8, "s", "-12345"),
    ("KO", 8, "s", "-12345"),
    ("KA", 8, "s", "-12345"),
    ("KT0", 8, "s", "-12345"),
    ("KT1", 8, "s", "-12345"),
    ("KT2", 8, "s", "-12345"),
    ("KT3", 8, "s", "-12345"),
    ("KT4", 8, "s", "-12345"),
    ("KT5", 8, "s", "-12345"),
    ("KT6", 8, "s", "-12345"),
    ("KT7", 8, "s", "-12345"),
    ("KT8", 8, "s", "-12345"),
    ("KT9", 8, "s", "-12345"),
    ("KF", 8, "s", "-12345"),
    ("KUSER0", 8, "s", "-12345"),
    ("KUSER1", 8, "s", "-12345"),
    ("KUSER2", 8, "s", "-12345"),
    ("KCMPNM", 8, "s", "-12345"),
    ("KNETWK", 8, "s", "-12345"),
    ("KDATRD", 8, "s", "-12345"),
    ("KINST", 8, "s", "-12345"),
]