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
waveform containers. Format-specific implementations are delegated to the
modules under ``shakelab.signals.libio``.

Imports of format-specific modules are intentionally lazy, in order to avoid
circular imports with ``shakelab.signals.base`` and to keep module loading
lightweight.

MiniSEED operations can use either the libmseed-backed ``cymseed`` module or
the pure-Python ``mseed`` module. The default backend is controlled by
``USE_LIBMSEED``. A backend can also be selected for an individual operation
by passing ``mseed_backend`` to :func:`reader` or :func:`writer`:

- ``"auto"`` uses the module-level ``USE_LIBMSEED`` preference;
- ``"libmseed"`` selects the ``cymseed`` backend;
- ``"python"`` selects the pure-Python backend.

The per-call option takes precedence over the module-level preference.
"""

import os
from glob import glob


USE_LIBMSEED = True

MSEED_BACKEND_AUTO = "auto"
MSEED_BACKEND_LIBMSEED = "libmseed"
MSEED_BACKEND_PYTHON = "python"

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


def get_mseed_module(backend=None):
    """
    Return the selected MiniSEED backend module.

    Parameters
    ----------
    backend : {"auto", "libmseed", "python"} or bool, optional
        MiniSEED backend selection. ``"auto"`` or ``None`` uses the
        module-level ``USE_LIBMSEED`` preference. ``"libmseed"`` selects
        the libmseed-backed ``cymseed`` module, while ``"python"`` selects
        the pure-Python ``mseed`` module.

        Boolean values are accepted for backward compatibility, with
        ``True`` selecting ``cymseed`` and ``False`` selecting ``mseed``.

    Returns
    -------
    module
        Selected MiniSEED backend module.

    Raises
    ------
    ValueError
        If the backend selector is not recognized.

    Notes
    -----
    Backend modules are imported lazily when this function is called.
    """
    if _use_libmseed_backend(backend):
        from shakelab.signals.libio import cymseed as mseed
        print("Using libmseed backend")
    else:
        from shakelab.signals.libio import mseed
        print("Using python backend")

    return mseed


def _use_libmseed_backend(backend=None):
    """
    Resolve a MiniSEED backend selector to a boolean preference.

    Parameters
    ----------
    backend : {"auto", "libmseed", "python"} or bool, optional
        Backend selector. ``"auto"`` or ``None`` resolves to the current
        value of ``USE_LIBMSEED``. Boolean values are returned unchanged.

    Returns
    -------
    bool
        ``True`` for the libmseed-backed ``cymseed`` module and ``False``
        for the pure-Python ``mseed`` module.

    Raises
    ------
    ValueError
        If the backend selector is not recognized.

    Notes
    -----
    This function resolves configuration only. It does not import a backend
    or check whether ``cymseed`` is installed.
    """
    if backend is None or backend == MSEED_BACKEND_AUTO:
        return USE_LIBMSEED

    if isinstance(backend, bool):
        return backend

    if backend == MSEED_BACKEND_LIBMSEED:
        return True

    if backend == MSEED_BACKEND_PYTHON:
        return False

    raise ValueError(
        "MiniSEED backend must be 'auto', 'libmseed' or 'python'"
    )


def reader(file_path=None, list_file=None, stream_collection=None,
           format=None, byte_order=None, **format_options):
    """
    Read waveform data from one or more files.

    Parameters
    ----------
    file_path : str, path-like or sequence, optional
        Waveform file path, wildcard pattern, directory, or sequence of
        waveform paths. Mutually exclusive with ``list_file``.
    list_file : str or path-like, optional
        Text file containing one waveform path per line. Mutually exclusive
        with ``file_path``. Relative paths are resolved against the directory
        containing the list file.
    stream_collection : StreamCollection, optional
        Existing collection to which loaded records are appended. If omitted,
        a new collection is created.
    format : str, optional
        Input format. If omitted, the format is inferred from each file
        extension.
    byte_order : str or None, optional
        Byte order used by formats that support explicit selection. ``None``
        lets the target reader use its default byte-order policy.
    **format_options
        Additional keyword arguments forwarded to the selected
        format-specific reader. For MiniSEED input, ``mseed_backend`` may
        be set to ``"auto"``, ``"libmseed"`` or ``"python"``. The default
        value, ``"auto"``, uses the module-level ``USE_LIBMSEED``
        preference.

    Returns
    -------
    StreamCollection
        Collection containing loaded records.

    Raises
    ------
    TypeError
        If an input source has an unsupported type or a format-specific
        option is invalid.
    ValueError
        If neither input source is provided, both are provided, or a format
        cannot be inferred.
    NotImplementedError
        If the selected format is known but unsupported.

    Examples
    --------
    Read MiniSEED using the configured default backend:

    >>> collection = reader("waveforms.mseed")

    Force the pure-Python MiniSEED backend for one operation:

    >>> collection = reader(
    ...     "waveforms.mseed",
    ...     mseed_backend="python",
    ... )

    Force the libmseed-backed implementation:

    >>> collection = reader(
    ...     "waveforms.mseed",
    ...     mseed_backend="libmseed",
    ... )
    """
    if file_path is not None and list_file is not None:
        raise ValueError(
            "file_path and list_file are mutually exclusive"
        )

    if file_path is None and list_file is None:
        raise ValueError(
            "Either file_path or list_file must be provided"
        )

    if list_file is not None:
        file_path = _read_file_list(list_file)

    stream_collection = _ensure_stream_collection(stream_collection)
    file_list = _expand_file_list(file_path)

    for path in file_list:
        current_format = _resolve_format(path, format)
        stream_collection = _read_file(
            path,
            current_format,
            stream_collection,
            byte_order=byte_order,
            **format_options,
        )

    stream_collection.sort()

    return stream_collection


def writer(input_data, file_path, format=None, byte_order=None,
           owrite=False, **format_options):
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
        Output format. If omitted, the format is inferred from
        ``file_path``.
    byte_order : str or None, optional
        Byte order used by formats that support explicit selection. ``None``
        lets the target writer use its default byte-order policy.
    owrite : bool, optional
        If ``True``, overwrite existing files.
    **format_options
        Additional keyword arguments forwarded to the selected
        format-specific writer. For MiniSEED output, ``mseed_backend`` may
        be set to ``"auto"``, ``"libmseed"`` or ``"python"``. The default
        value, ``"auto"``, uses the module-level ``USE_LIBMSEED``
        preference.

    Raises
    ------
    TypeError
        If ``input_data`` is not a supported ShakeLab waveform container or
        a format-specific option is invalid.
    ValueError
        If the format is not recognized.
    NotImplementedError
        If the selected format is known but unsupported.

    Examples
    --------
    Write MiniSEED using the pure-Python backend:

    >>> writer(
    ...     collection,
    ...     "waveforms.mseed",
    ...     mseed_backend="python",
    ...     encoding=11,
    ...     reclen=512,
    ... )
    """
    stream_collection = _as_stream_collection(input_data)
    current_format = _resolve_format(file_path, format)

    if current_format == "mseed":
        _write_mseed(
            stream_collection,
            file_path,
            owrite=owrite,
            **format_options,
        )

    elif current_format == "sac":
        _write_sac(
            stream_collection,
            file_path,
            byte_order=byte_order,
            owrite=owrite,
            **format_options,
        )

    elif current_format in UNSUPPORTED_FORMATS:
        raise NotImplementedError(
            "{0}: format not yet implemented".format(current_format)
        )

    else:
        raise ValueError("{0}: format not recognized".format(
            current_format
        ))


def _read_file(path, format, stream_collection, byte_order=None,
               **format_options):
    """
    Read one waveform file into a StreamCollection.

    Parameters
    ----------
    path : str
        Input waveform path.
    format : str
        Normalized waveform format identifier.
    stream_collection : StreamCollection
        Collection to which decoded records are appended.
    byte_order : str or None, optional
        Explicit byte order for readers that support it.
    **format_options
        Format-specific options. For MiniSEED, ``mseed_backend`` selects
        the backend and is consumed by this dispatcher rather than passed
        to the backend reader.

    Returns
    -------
    StreamCollection
        Updated waveform collection.
    """
    if format == "mseed":
        options = dict(format_options)
        backend = options.pop(
            "mseed_backend",
            MSEED_BACKEND_AUTO,
        )

        use_libmseed = _use_libmseed_backend(backend)
        mseed = get_mseed_module(use_libmseed)

        if byte_order is not None and not use_libmseed:
            options["byte_order"] = byte_order

        return mseed.msread(
            path,
            stream_collection=stream_collection,
            **options,
        )

    if format == "sac":
        from shakelab.signals.libio import sac

        options = dict(format_options)

        if byte_order is not None:
            options["byte_order"] = byte_order

        record = sac.sacread(path, **options)
        stream_collection.append(record)

        return stream_collection

    if format == "dyna":
        from shakelab.signals.libio import dyna

        record = dyna.dynaread(path, **format_options)
        stream_collection.append(record)

        return stream_collection

    if format == "tdms":
        from shakelab.signals.libio import tdms

        tdms_collection = tdms.tdms_stream_read(
            path,
            **format_options,
        )

        for stream in tdms_collection:
            stream_collection.append(stream)

        return stream_collection

    if format in UNSUPPORTED_FORMATS:
        raise NotImplementedError(
            "{0}: format not yet implemented".format(format)
        )

    raise ValueError("{0}: format not recognized".format(format))


def _write_mseed(stream_collection, file_path, owrite=False,
                 **format_options):
    """
    Write a StreamCollection using the selected MiniSEED backend.

    Parameters
    ----------
    stream_collection : StreamCollection
        Waveform collection to write.
    file_path : str or path-like
        Output MiniSEED path.
    owrite : bool, optional
        If ``True``, overwrite an existing output file.
    **format_options
        MiniSEED writer options. ``mseed_backend`` selects ``"auto"``,
        ``"libmseed"`` or ``"python"`` and is consumed by this dispatcher.
        All remaining options are forwarded to the selected ``mswrite()``
        implementation.

    Raises
    ------
    FileExistsError
        If the output file already exists and ``owrite`` is ``False``.
    ValueError
        If the MiniSEED backend selector is invalid.
    """
    if os.path.exists(file_path) and not owrite:
        raise FileExistsError("File exists: {0}".format(file_path))

    options = dict(format_options)
    backend = options.pop(
        "mseed_backend",
        MSEED_BACKEND_AUTO,
    )

    mseed = get_mseed_module(backend)

    mseed.mswrite(
        file_path,
        stream_collection,
        **options,
    )


def _write_sac(stream_collection, directory, byte_order=None,
               owrite=False, **format_options):
    """
    Write a StreamCollection to one SAC file per record.

    Additional format-specific options are forwarded unchanged to the SAC
    writer.
    """
    from shakelab.signals.libio import sac

    os.makedirs(directory, exist_ok=True)

    options = dict(format_options)
    options["owrite"] = owrite

    if byte_order is not None:
        options["byte_order"] = byte_order

    for stream in stream_collection:
        for record in stream:
            sid = _safe_filename_part(record.head.sid)
            time = _safe_filename_part(record.head.time.iso8601)
            filename = "{0}_{1}.sac".format(sid, time)
            path = os.path.join(directory, filename)

            sac.sacwrite(
                path,
                record,
                **options,
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


def _read_file_list(list_file):
    """
    Read waveform paths from a text file.

    Empty lines and lines starting with ``#`` are ignored. Relative paths
    are resolved against the directory containing the list file.

    Parameters
    ----------
    list_file : str or path-like
        Path to a text file containing one waveform path per line.

    Returns
    -------
    list of str
        Waveform paths contained in the list file.

    Raises
    ------
    TypeError
        If ``list_file`` is not a path-like object.
    """
    if not isinstance(list_file, (str, os.PathLike)):
        raise TypeError(
            "list_file must be a string or path-like object"
        )

    list_file = os.path.abspath(os.fspath(list_file))
    base_directory = os.path.dirname(list_file)
    file_list = []

    with open(list_file, "r", encoding="utf-8") as fid:
        for line in fid:
            path = line.strip()

            if not path or path.startswith("#"):
                continue

            path = os.path.expanduser(path)
            path = os.path.expandvars(path)

            if not os.path.isabs(path):
                path = os.path.join(base_directory, path)

            file_list.append(path)

    return file_list


def _expand_file_list(file_path):
    """
    Expand waveform input paths.

    Parameters
    ----------
    file_path : str, path-like or sequence
        Waveform file path, wildcard pattern, directory, or sequence of
        paths.

    Returns
    -------
    list of str
        Expanded waveform file paths.

    Raises
    ------
    TypeError
        If ``file_path`` or one of its elements has an unsupported type.
    """
    if isinstance(file_path, (str, os.PathLike)):
        path = os.path.expanduser(os.fspath(file_path))
        path = os.path.expandvars(path)

        if os.path.isdir(path):
            return sorted(
                os.path.join(path, name)
                for name in os.listdir(path)
                if os.path.isfile(os.path.join(path, name))
            )

        matches = sorted(glob(path))

        if matches:
            return matches

        return [path]

    if isinstance(file_path, (list, tuple)):
        file_list = []

        for path in file_path:
            file_list.extend(_expand_file_list(path))

        return file_list

    raise TypeError(
        "file_path must be a path or a sequence of paths"
    )


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


def _safe_filename_part(value):
    """
    Return a filesystem-friendly filename component.
    """
    text = str(value)

    for char in ("/", "\\", ":", " "):
        text = text.replace(char, "_")

    return text