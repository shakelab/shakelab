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
Waveform input/output dispatch utilities.

This module provides high-level reader and writer functions for ShakeLab
waveform containers.  Format-specific implementations are delegated to the
modules under ``shakelab.signals.libio``.

Imports of format-specific modules are intentionally lazy, in order to avoid
circular imports with ``shakelab.signals.base`` and to keep module loading
lightweight.
"""

import os
from glob import glob


USE_LIBMSEED = True

MSEED_FORMATS = ("mseed", "miniseed", "ms", "seed")
SAC_FORMATS = ("sac",)
DYNA_FORMATS = ("dyna",)
TDMS_FORMATS = ("tdms",)

UNSUPPORTED_FORMATS = (
    "ascii",
    "seisan",
    "seg2",
    "dat",
    "gse",
    "reftek",
    "itaca",
)

FORMAT_ALIASES = {
    "miniseed": "mseed",
    "ms": "mseed",
    "seed": "mseed",
}

EXTENSION_FORMATS = {
    ".ms": "mseed",
    ".mseed": "mseed",
    ".miniseed": "mseed",
    ".seed": "mseed",
    ".sac": "sac",
    ".txt": "dyna",
    ".tdms": "tdms",
}


def get_mseed_module(use_libmseed=USE_LIBMSEED):
    """
    Return the selected MiniSEED backend module.

    Parameters
    ----------
    use_libmseed : bool, optional
        If ``True``, use the libmseed-backed module.  If ``False``, use the
        pure-Python MiniSEED implementation.

    Returns
    -------
    module
        MiniSEED I/O module.
    """
    if use_libmseed:
        from shakelab.signals.libio import cymseed as mseed
    else:
        from shakelab.signals.libio import mseed

    return mseed


def reader(file_path, stream_collection=None, format=None,
           byte_order=None, is_db=False):
    """
    Read waveform data from one or more files.

    Parameters
    ----------
    file_path : str or list
        Path to a waveform file, wildcard pattern, directory, list of files,
        or a text file containing a list of file paths if ``is_db=True``.
    stream_collection : StreamCollection, optional
        Existing collection to which loaded records are appended.  If omitted,
        a new collection is created.
    format : str, optional
        Input format.  If omitted, the format is inferred from the file
        extension.
    byte_order : str or None, optional
        Byte order used by formats that require it.  ``None`` lets supporting
        readers attempt automatic detection.
    is_db : bool, optional
        If ``True``, interpret ``file_path`` as a text file listing waveform
        paths.

    Returns
    -------
    StreamCollection
        Collection containing loaded records.

    Raises
    ------
    TypeError
        If ``file_path`` has an unsupported type.
    ValueError
        If a format cannot be inferred or is not recognized.
    NotImplementedError
        If the selected format is known but unsupported.
    """
    stream_collection = _ensure_stream_collection(stream_collection)
    file_list = _expand_file_list(file_path, is_db=is_db)

    for path in file_list:
        current_format = _resolve_format(path, format)
        stream_collection = _read_file(
            path,
            current_format,
            stream_collection,
            byte_order=byte_order,
        )

    return stream_collection


def writer(input_data, file_path, format=None, byte_order=None,
           owrite=False):
    """
    Write waveform data to disk.

    Parameters
    ----------
    input_data : Record, Stream or StreamCollection
        ShakeLab waveform container to write.
    file_path : str
        Output file path for single-file formats or output directory for
        multi-file formats such as SAC.
    format : str, optional
        Output format.  If omitted, the format is inferred from
        ``file_path``.
    byte_order : str or None, optional
        Byte order used by formats that require it.  ``None`` lets the target
        writer use its own default.
    owrite : bool, optional
        If ``True``, overwrite existing files.

    Raises
    ------
    TypeError
        If ``input_data`` is not a supported ShakeLab waveform container.
    ValueError
        If the format is not recognized.
    NotImplementedError
        If the selected format is known but unsupported.
    """
    stream_collection = _as_stream_collection(input_data)
    current_format = _resolve_format(file_path, format)

    if current_format == "mseed":
        _write_mseed(stream_collection, file_path, owrite=owrite)

    elif current_format == "sac":
        _write_sac(
            stream_collection,
            file_path,
            byte_order=byte_order,
            owrite=owrite,
        )

    elif current_format in UNSUPPORTED_FORMATS:
        raise NotImplementedError(
            "{0}: format not yet implemented".format(current_format)
        )

    else:
        raise ValueError("{0}: format not recognized".format(
            current_format
        ))


def _read_file(path, format, stream_collection, byte_order=None):
    """
    Read one file and append its content to a StreamCollection.
    """
    if format == "mseed":
        mseed = get_mseed_module(USE_LIBMSEED)
        return mseed.msread(path, stream_collection=stream_collection)

    if format == "sac":
        from shakelab.signals.libio import sac

        record = sac.sacread(path, byte_order=byte_order)
        stream_collection.append(record)

        return stream_collection

    if format == "dyna":
        from shakelab.signals.libio import dyna

        record = dyna.dynaread(path)
        stream_collection.append(record)

        return stream_collection

    if format == "tdms":
        from shakelab.signals.libio import tdms

        tdms_collection = tdms.tdms_stream_read(path)

        for stream in tdms_collection:
            stream_collection.append(stream)

        return stream_collection

    if format in UNSUPPORTED_FORMATS:
        raise NotImplementedError(
            "{0}: format not yet implemented".format(format)
        )

    raise ValueError("{0}: format not recognized".format(format))


def _write_mseed(stream_collection, file_path, owrite=False):
    """
    Write a StreamCollection to MiniSEED.
    """
    if os.path.exists(file_path) and not owrite:
        raise FileExistsError("File exists: {0}".format(file_path))

    mseed = get_mseed_module(USE_LIBMSEED)
    mseed.mswrite(file_path, stream_collection)


def _write_sac(stream_collection, directory, byte_order=None,
               owrite=False):
    """
    Write a StreamCollection to one SAC file per record.
    """
    from shakelab.signals.libio import sac

    os.makedirs(directory, exist_ok=True)

    for stream in stream_collection:
        for record in stream:
            sid = _safe_filename_part(record.head.sid)
            time = _safe_filename_part(record.head.time.iso8601)
            filename = "{0}_{1}.sac".format(sid, time)
            path = os.path.join(directory, filename)

            sac.sacwrite(
                path,
                record,
                byte_order=byte_order,
                owrite=owrite,
            )


def _ensure_stream_collection(stream_collection=None):
    """
    Return a StreamCollection instance.
    """
    if stream_collection is not None:
        return stream_collection

    from shakelab.signals.base import StreamCollection

    return StreamCollection()


def _as_stream_collection(input_data):
    """
    Convert a ShakeLab waveform container to a StreamCollection.
    """
    from shakelab.signals.base import Record, Stream, StreamCollection

    if isinstance(input_data, StreamCollection):
        return input_data

    if isinstance(input_data, (Record, Stream)):
        stream_collection = StreamCollection()
        stream_collection.append(input_data)

        return stream_collection

    raise TypeError(
        "input_data must be a Record, Stream or StreamCollection"
    )


def _expand_file_list(file_path, is_db=False):
    """
    Expand input paths to an ordered list of files.
    """
    if isinstance(file_path, (list, tuple)):
        return list(file_path)

    if not isinstance(file_path, str):
        raise TypeError("file_path must be str, list or tuple")

    if is_db:
        return _read_path_list(file_path)

    if _has_wildcards(file_path):
        return sorted(glob(file_path))

    if os.path.isdir(file_path):
        return sorted(glob(os.path.join(file_path, "*")))

    return [file_path]


def _read_path_list(path):
    """
    Read a text file containing one waveform path per line.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError("List file not found: {0}".format(path))

    with open(path, "r", encoding="utf-8") as fid:
        return [
            line.strip()
            for line in fid
            if line.strip() and not line.lstrip().startswith("#")
        ]


def _resolve_format(path, format=None):
    """
    Resolve and normalize a waveform format.
    """
    if format is not None:
        return _normalize_format(format)

    extension = os.path.splitext(path)[1].lower()

    if extension not in EXTENSION_FORMATS:
        raise ValueError(
            "Unable to infer waveform format from extension: {0}".format(
                path
            )
        )

    return EXTENSION_FORMATS[extension]


def _normalize_format(format):
    """
    Normalize format aliases.
    """
    if format is None:
        return None

    value = format.lower().strip()

    return FORMAT_ALIASES.get(value, value)


def _has_wildcards(path):
    """
    Return True if path contains shell-style wildcard characters.
    """
    return "*" in path or "?" in path or "[" in path


def _safe_filename_part(value):
    """
    Return a filesystem-friendly filename component.
    """
    text = str(value)

    for char in ("/", "\\", ":", " "):
        text = text.replace(char, "_")

    return text