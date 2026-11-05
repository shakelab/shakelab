# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Reading and writing of waveforms and channel metadata.

:func:`reader` reads waveform files of any supported format into a
:class:`~shakelab.signals.base.StreamCollection`, :func:`writer` writes
records, streams or collections, and :func:`metadata_reader` reads
channel metadata (StationXML, RESP or JSON) into a
:class:`~shakelab.signals.metadata.MetadataCollection`. The formats are
implemented in :mod:`shakelab.signals.libio`; :func:`formats` lists
them.

Files can be given as a path, a wildcard pattern, a directory, a
sequence of them, or a text file listing them (``list_file``). Unless
``format`` is given, the format of each file is inferred from its
extension or, for files without a known extension (e.g. the daily files
of an SDS archive) or with an ambiguous one (``.asc``), from its first
bytes; files of unknown format found in directories or patterns are
skipped with one :class:`IoWarning`. Pickle files are never recognized
from their content.

MiniSEED is read and written by the libmseed-based ``cymseed`` module
or by the pure-Python ``mseed`` module: :data:`USE_LIBMSEED` sets the
default, and the ``mseed_backend`` option (``"auto"``, ``"libmseed"``
or ``"python"``) selects one for a single call.

The ``pickle`` format (``.pkl``, ``.pickle``) saves records, streams
or collections as Python objects, with all their metadata. It is meant
for intermediate results: a pickle file can run code when it is read,
so only files from trusted sources should be read, and files written
by one version of ShakeLab may not be readable by another.

:func:`reader` can also select channels (``select``, with the patterns
of :mod:`shakelab.signals.codes`) and a time window (``starttime``,
``endtime``, as in :meth:`~shakelab.signals.base.Record.cut`), and
read an archive organized by channel and day (``layout="sds"``, the
SeisComP Data Structure), opening only the files of the channels and
days requested.

Format modules are imported only when a file of that format is read or
written.
"""

from __future__ import annotations

import importlib
import inspect
import os
import pickle
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from fnmatch import fnmatchcase
from glob import glob
from typing import TYPE_CHECKING, Any, Callable, Optional

from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import warn_external
from shakelab.signals.codes import SidPattern, parse_sid_pattern

if TYPE_CHECKING:
    from shakelab.signals.base import Record, StreamCollection
    from shakelab.signals.metadata import MetadataCollection

__all__ = [
    "IoWarning",
    "USE_LIBMSEED",
    "formats",
    "metadata_reader",
    "reader",
    "writer",
]

# Default MiniSEED backend: libmseed (cymseed) if True, else Python.
USE_LIBMSEED = True


class IoWarning(UserWarning):
    """Warning about files that are not read."""


# ----------------------------------------------------------------------------
# Formats
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class _Format:
    """Description of a waveform format and of its implementation."""

    description: str
    module: str
    read: Optional[str]
    write: Optional[str]
    extensions: tuple[str, ...] = ()
    suffix: str = ""
    per_record: bool = False
    detect: Optional[str] = None


_FORMATS = {
    "mseed": _Format(
        "MiniSEED 2 and 3", "shakelab.signals.libio.mseed",
        "msread", "mswrite", (".mseed", ".miniseed", ".ms", ".seed"),
        ".mseed", detect="is_mseed",
    ),
    "sac": _Format(
        "SAC binary", "shakelab.signals.libio.sac",
        "sacread", "sacwrite", (".sac",), ".sac", per_record=True,
        detect="is_sac",
    ),
    "dyna": _Format(
        "DYNA 1.2 (ESM, ITACA)", "shakelab.signals.libio.dyna",
        "dynaread", "dynawrite", (".asc",), ".asc", per_record=True,
        detect="is_dyna",
    ),
    "simpleascii": _Format(
        "Simple ASCII (EarthScope)", "shakelab.signals.libio.simpleascii",
        "simpleasciiread", "simpleasciiwrite", (), ".txt",
        detect="is_simpleascii",
    ),
    "ascii": _Format(
        "Text tables", "shakelab.signals.libio.ascii",
        "asciiread", "asciiwrite", (".csv",), ".csv",
    ),
    "tdms": _Format(
        "NI TDMS", "shakelab.signals.libio.tdms",
        "tdmsread", None, (".tdms",), detect="is_tdms",
    ),
    "pickle": _Format(
        "Python pickle (trusted files only)", __name__,
        "_read_pickle", "_write_pickle", (".pkl", ".pickle"), ".pkl",
    ),
}

_ALIASES = {
    "miniseed": "mseed",
    "ms": "mseed",
    "seed": "mseed",
    "pkl": "pickle",
}

_UNSUPPORTED = ("gse", "reftek", "seg2", "segy", "seisan")

# Extensions shared by several formats: the content decides.
_AMBIGUOUS = (".asc",)

# Bytes read to recognize a format from the content.
_HEAD_SIZE = 1024

_EXTENSIONS = {
    extension: name
    for name, fmt in _FORMATS.items()
    for extension in fmt.extensions
}

_MSEED_BACKENDS = {
    "libmseed": "shakelab.signals.libio.cymseed",
    "python": "shakelab.signals.libio.mseed",
}

# Fields of the file name templates of writer() (as in seedlink).
_NAME_FIELDS = ("net", "sta", "loc", "cha", "sid", "start", "time", "year",
                "jday")

# Archive layouts read by reader().
_LAYOUTS = ("sds",)

# Metadata types recognized from file names (see MetadataCollection.read).
_METADATA_TYPES = ("xml", "sxml", "resp", "json")


def formats() -> dict[str, dict]:
    """
    Return the supported waveform formats.

    Returns
    -------
    dict
        For each format name: ``description``, ``extensions`` (used to
        infer the format from file names), ``read`` and ``write``
        (whether it can be read and written), ``per_record`` (True
        for formats holding one record per file) and ``detect`` (True
        if files are recognized from their content).

    Examples
    --------
    >>> sorted(formats())
    ['ascii', 'dyna', 'mseed', 'pickle', 'sac', 'simpleascii', 'tdms']
    >>> formats()["sac"]["extensions"]
    ('.sac',)
    """
    return {
        name: {
            "description": fmt.description,
            "extensions": fmt.extensions,
            "read": fmt.read is not None,
            "write": fmt.write is not None,
            "per_record": fmt.per_record,
            "detect": fmt.detect is not None,
        }
        for name, fmt in _FORMATS.items()
    }


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


def reader(
    source: Any = None,
    *,
    list_file: Any = None,
    format: Optional[str] = None,
    layout: Optional[str] = None,
    select: Any = None,
    starttime: Any = None,
    endtime: Any = None,
    metadata: Any = None,
    stream_collection: Optional[StreamCollection] = None,
    **options: Any,
) -> StreamCollection:
    """
    Read waveform files into a StreamCollection.

    Parameters
    ----------
    source : str, path-like or sequence, optional
        File, wildcard pattern (``*``, ``?``, ``[...]``), directory
        (its files, not its subdirectories), or a list or tuple of
        them. ``~`` and environment variables are expanded.
    list_file : str or path-like, optional
        Text file listing the files to read, one per line, in any of
        the forms of ``source``. Empty lines and lines starting with
        ``#`` are ignored; relative paths refer to the folder of the
        list. The files listed are read after those of ``source``.
    format : str, optional
        Format of all files (see :func:`formats`). If omitted, it is
        inferred from the extension of each file or, without a known
        extension or with ``.asc``, from its first bytes.
    layout : {"sds"}, optional
        Organization of an archive whose root is ``source``: ``"sds"``
        is the SeisComP Data Structure, one MiniSEED file per channel
        and day in ``YEAR/NET/STA/CHA.D/NET.STA.LOC.CHA.D.YEAR.DAY``.
        Only the files of the channels and days requested are read
        (also the day before ``starttime``, since a record is stored in
        the file of the day it begins). ``starttime`` and ``endtime``
        are required; missing days are not an error.
    select : str or sequence of str, optional
        Channels to keep: patterns ``NET.STA.LOC.CHA`` with ``*``,
        ``?``, ``[...]`` and ``{A,B}``, codes missing at the end
        matching anything (see :mod:`shakelab.signals.codes`).
    starttime, endtime : Date or str, optional
        Time window to keep, bounds included: records outside it are
        dropped and the others cut to it (see
        :meth:`~shakelab.signals.base.Record.cut`). The names differ
        from the ``start`` option of some formats (the start time of
        records without one, e.g. text tables).
    metadata : MetadataCollection, ChannelEpoch, path or sequence, optional
        Channel metadata assigned to the records read, as with
        :meth:`~shakelab.signals.base.StreamCollection.set_metadata`.
        Paths are read with :func:`metadata_reader`, before the data;
        the metadata are assigned after the selection. The records
        already in ``stream_collection`` are not modified.
    stream_collection : StreamCollection, optional
        Collection to which the records are added; a new one if
        omitted.
    **options
        Options of the format reader, the same for all files, e.g.
        ``byte_order`` for SAC and for the Python MiniSEED backend,
        ``mseed_backend`` (``"auto"``, ``"libmseed"`` or ``"python"``)
        for MiniSEED, ``rate`` or ``columns`` for text tables (see the
        readers in :mod:`shakelab.signals.libio`).

    Returns
    -------
    StreamCollection
        Collection with the records read, sorted.

    Raises
    ------
    ValueError
        If no source is given, if the format of a file given by name
        is not recognized, if a pattern is invalid, if the time window
        is reversed, or if ``layout`` is unknown or given without
        ``starttime`` and ``endtime``.
    FileNotFoundError
        If a file, or the root of an archive, does not exist, or if a
        pattern matches no file.
    TypeError
        If a source has an unsupported type, or if an option is not
        accepted by the reader of a format.
    NotImplementedError
        If ``format`` is known but not supported yet.

    Warns
    -----
    IoWarning
        If files of unknown format, found in directories or patterns,
        are skipped (one warning listing them; hidden files are
        ignored), or if an archive has no data for the request.
    MetadataWarning
        If ``metadata`` is given and some channels read have no epoch
        in it (one warning listing them).

    Examples
    --------
    >>> sc = reader("waveforms.mseed")
    >>> sc = reader("event/", metadata="station.xml")
    >>> sc = reader(["a.sac", "b.sac"], byte_order="le")
    >>> sc = reader(list_file="files.txt")
    >>> sc = reader("data/*.mseed", mseed_backend="python")
    >>> sc = reader("SDS/2026/IV/MODE/HNZ.D/")   # files without extension
    >>> sc = reader("event/", select="IV.MODE..HN?",
    ...             starttime="2012-05-29T07:00:00",
    ...             endtime="2012-05-29T07:02:00")
    >>> sc = reader("/archive", layout="sds", select="IV.{MODE,MIRA}",
    ...             starttime="2026-05-03T23:50:00",
    ...             endtime="2026-05-04T00:20:00")
    """
    # Arguments and metadata are checked first, so that an error stops
    # before the data.
    patterns = _patterns(select)
    start, end = _window(starttime, endtime)
    metadata = _as_metadata(metadata)
    name = None if format is None else _format_name(format)
    if layout is None:
        files = _input_files(source, list_file)
    else:
        files = _archive_files(source, list_file, layout, patterns,
                               start, end)
        name = name or "mseed"
    selecting = patterns is not None or start is not None or (
        end is not None
    )

    # With a selection or with metadata, the data are read into a new
    # collection, so that the records already in stream_collection are
    # left as they are.
    if stream_collection is None or selecting or metadata is not None:
        from shakelab.signals.base import StreamCollection

        collection = StreamCollection()
    else:
        collection = stream_collection

    skipped = []
    for path, explicit in files:
        current = name or _format_of(path)
        if current is None:
            if explicit:
                raise ValueError(
                    f"{path}: format not recognized from the file name "
                    "or content; set format"
                )
            skipped.append(path)
            continue
        _read_file(path, current, collection, options)

    if skipped:
        _warn_skipped(skipped)

    if selecting:
        collection = _collection(
            _selected(collection.records(), patterns, start, end)
        )
    if layout is not None and not len(collection):
        warn_external(
            f"no data in the archive for {select or 'all channels'} "
            f"from {start} to {end}",
            IoWarning,
        )
    if metadata is not None:
        collection.set_metadata(metadata)
    if stream_collection is not None and collection is not stream_collection:
        for stream in collection:
            stream_collection.append(stream)
        collection = stream_collection

    collection.sort()

    return collection


def metadata_reader(
    source: Any = None,
    *,
    list_file: Any = None,
    format: Optional[str] = None,
) -> MetadataCollection:
    """
    Read channel metadata files into a MetadataCollection.

    Parameters
    ----------
    source : str, path-like or sequence, optional
        File, wildcard pattern, directory, or a list or tuple of them,
        as in :func:`reader`.
    list_file : str or path-like, optional
        Text file listing the files to read, as in :func:`reader`.
    format : str, optional
        Format of all files: ``"xml"`` (StationXML), ``"resp"`` or
        ``"json"``. If omitted, it is inferred from each file name
        (see :meth:`~shakelab.signals.metadata.MetadataCollection.read`).

    Returns
    -------
    MetadataCollection
        Collection with the content of all files.

    Raises
    ------
    ValueError
        If no source is given, or if the format of a file given by
        name is unknown.
    FileNotFoundError
        If a file does not exist or a pattern matches no file.
    TypeError
        If a source has an unsupported type.

    Warns
    -----
    IoWarning
        If files of unknown format, found in directories or patterns,
        are skipped (one warning listing them).
    MetadataWarning
        If epochs of the same channel from different files overlap.

    Notes
    -----
    Files are added in order; an epoch equal to one already read is
    ignored, so the same metadata can be listed twice. Text, bytes or
    open files are read with
    :class:`~shakelab.signals.metadata.MetadataCollection` directly.

    Examples
    --------
    >>> meta = metadata_reader(["station.xml", "RESP.IU.ANMO.00.BHZ"])
    >>> meta = metadata_reader("metadata/*.xml")
    """
    from shakelab.signals.metadata import MetadataCollection, _file_type

    collection = MetadataCollection()
    skipped = []
    for path, explicit in _input_files(source, list_file):
        if (
            format is None
            and not explicit
            and _file_type(path) not in _METADATA_TYPES
        ):
            skipped.append(path)
            continue
        collection.read(path, format=format, append=True)

    if skipped:
        _warn_skipped(skipped)

    return collection


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------


def writer(
    data: Any,
    target: Any,
    *,
    format: Optional[str] = None,
    owrite: bool = False,
    **options: Any,
) -> None:
    """
    Write waveforms to one or more files.

    Parameters
    ----------
    data : Record, Stream or StreamCollection
        Waveforms to write.
    target : str or path-like
        One of:

        - a file name: all records in one file. Formats holding one
          record per file (SAC, DYNA) accept a single record only;
        - a directory (existing, or ending with a separator): one file
          per channel, or per record for formats holding one record
          per file, named ``NET.STA.LOC.CHA`` (with ``.{start}`` for
          one record per file) and the extension of the format;
        - a file name template with the fields of the SeedLink
          archives (:mod:`shakelab.signals.seedlink.archiver`):
          ``{net}``, ``{sta}``, ``{loc}``, ``{cha}``, ``{sid}``, the
          start time ``{time}`` (a ``datetime``, with a ``strftime``
          format such as ``{time:%Y%m%dT%H%M%S}``), ``{year}``,
          ``{jday}`` (day of the year, three digits) and ``{start}``,
          short for ``{time:%Y%m%dT%H%M%S}``. The records with the same
          name go to the same file; missing folders are created.
    format : str, optional
        Output format (see :func:`formats`). If omitted, it is
        inferred from the extension of ``target``; it is required
        when ``target`` is a directory.
    owrite : bool, optional
        If True, overwrite existing files. Default is False.
    **options
        Options of the format writer, e.g. ``encoding``, ``reclen``
        and ``mseed_backend`` for MiniSEED, ``byte_order`` for SAC,
        ``fmt`` for DYNA (see the writers in
        :mod:`shakelab.signals.libio`).

    Raises
    ------
    TypeError
        If ``data`` is not a waveform container, ``target`` is not a
        path, or an option is not accepted by the format writer.
    ValueError
        If the format cannot be inferred or is unknown, if a format
        holding one record per file is given several records for one
        file, or if a template has unknown fields.
    FileExistsError
        If a file exists and ``owrite`` is False (no file is written).
    NotImplementedError
        If the format cannot be written.

    Examples
    --------
    >>> writer(sc, "waveforms.mseed", encoding=11)
    >>> writer(record, "record.sac")
    >>> writer(sc, "sac_out/", format="sac")
    >>> writer(sc, "{sta}/{sta}.{cha}.{start}.sac")
    >>> writer(sc, "{sta}.mseed")    # one MiniSEED file per station
    >>> writer(sc, "{year}/{jday}/{net}.{sta}.{loc}.{cha}.mseed")
    """
    from shakelab.signals.base import Record, Stream, StreamCollection

    if not isinstance(data, (Record, Stream, StreamCollection)):
        raise TypeError(
            "data must be a Record, a Stream or a StreamCollection, got "
            f"{type(data).__name__}"
        )
    if not isinstance(target, (str, os.PathLike)):
        raise TypeError("target must be a path or a file name template")

    target = os.path.expanduser(os.fspath(target))
    template = "{" in target
    directory = not template and (
        target.endswith(("/", os.sep)) or os.path.isdir(target)
    )

    if format is not None:
        name = _format_name(format)
    elif directory:
        raise ValueError(f"{target}: set format to write into a directory")
    else:
        name = _format_from_name(target)
        if name is None:
            raise ValueError(
                f"{target}: format not recognized from the file name; "
                "set format"
            )

    fmt = _FORMATS[name]
    if fmt.write is None:
        raise NotImplementedError(f"{name}: writing is not supported")

    options = dict(options)
    write = _function(name, "write", options)

    if not template and not directory:
        content = data
        if fmt.per_record:
            records = data.records()
            if len(records) != 1:
                raise ValueError(
                    f"{name}: one record per file, {len(records)} "
                    "records given; write into a directory or use a "
                    "name template such as "
                    f"'{{sta}}.{{cha}}.{{start}}{fmt.suffix}'"
                )
            content = records[0]
        write(target, content, owrite=owrite, **options)
        return

    if directory:
        fields = "{net}.{sta}.{loc}.{cha}"
        if fmt.per_record:
            fields += ".{start}"
        target = os.path.join(target, fields + fmt.suffix)

    files = _names(data.records(), target)

    if fmt.per_record:
        for path, records in files.items():
            if len(records) > 1:
                raise ValueError(
                    f"{path}: {len(records)} records with the same name; "
                    "add {start} to the name template"
                )
    if not owrite:
        for path in files:
            if os.path.exists(path):
                raise FileExistsError(f"File exists: {path}")

    for path, records in files.items():
        folder = os.path.dirname(path)
        if folder:
            os.makedirs(folder, exist_ok=True)
        content = records[0] if fmt.per_record else _collection(records)
        write(path, content, owrite=owrite, **options)


# ----------------------------------------------------------------------------
# Pickle
# ----------------------------------------------------------------------------


def _read_pickle(
    source: Any,
    stream_collection: StreamCollection,
) -> StreamCollection:
    """
    Read a pickle file written by :func:`writer` (trusted files only).

    The records of the Record, Stream or StreamCollection in the file
    are appended to ``stream_collection``.

    Raises
    ------
    TypeError
        If the file does not hold a ShakeLab waveform container.
    """
    from shakelab.signals.base import Record, Stream, StreamCollection

    with open(source, "rb") as fid:
        content = pickle.load(fid)

    if not isinstance(content, (Record, Stream, StreamCollection)):
        raise TypeError(
            f"{source}: the pickle file holds a "
            f"{type(content).__name__}, not a Record, Stream or "
            "StreamCollection"
        )
    for record in content.records():
        stream_collection.append(record)

    return stream_collection


def _write_pickle(target: Any, data: Any, owrite: bool = False) -> None:
    """
    Write a Record, Stream or StreamCollection as a pickle file.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    """
    if os.path.exists(target) and not owrite:
        raise FileExistsError(f"File exists: {target}")
    with open(target, "wb") as fid:
        pickle.dump(data, fid, protocol=pickle.HIGHEST_PROTOCOL)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _read_file(
    path: str,
    name: str,
    collection: StreamCollection,
    options: dict,
) -> None:
    """Read one file of a known format into a collection."""
    options = dict(options)
    read = _function(name, "read", options)
    if _FORMATS[name].per_record:
        collection.append(read(path, **options))
    else:
        read(path, stream_collection=collection, **options)


def _function(name: str, action: str, options: dict) -> Callable:
    """
    Return the reader or writer of a format, after checking options.

    ``mseed_backend`` is removed from ``options`` for MiniSEED.
    """
    fmt = _FORMATS[name]
    module = fmt.module
    if name == "mseed":
        module = _mseed_module(options.pop("mseed_backend", None))
    function = getattr(
        importlib.import_module(module),
        fmt.read if action == "read" else fmt.write,
    )

    fixed = ("source", "target", "stream_collection", "record", "data",
             "owrite")
    parameters = inspect.signature(function).parameters
    accepted = [key for key in parameters if key not in fixed]
    if name == "mseed":
        accepted.append("mseed_backend")
    unknown = sorted(set(options) - set(accepted))
    if unknown:
        raise TypeError(
            f"{name}: unknown option(s) {', '.join(unknown)}; accepted: "
            f"{', '.join(accepted) or 'none'}"
        )

    return function


def _mseed_module(backend: Optional[str]) -> str:
    """Module of a MiniSEED backend ("auto", "libmseed" or "python")."""
    if backend is None or backend == "auto":
        backend = "libmseed" if USE_LIBMSEED else "python"
    if backend not in _MSEED_BACKENDS:
        raise ValueError(
            "mseed_backend must be 'auto', 'libmseed' or 'python', got "
            f"{backend!r}"
        )
    return _MSEED_BACKENDS[backend]


def _format_name(format: Any) -> str:
    """Normalized name of a format given by the user."""
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    name = format.strip().lower()
    name = _ALIASES.get(name, name)
    if name in _UNSUPPORTED:
        raise NotImplementedError(f"{name}: format not supported yet")
    if name not in _FORMATS:
        raise ValueError(
            f"{format}: unknown format; supported: "
            f"{', '.join(sorted(_FORMATS))}"
        )
    return name


def _format_from_name(path: str) -> Optional[str]:
    """Format inferred from a file extension, or None."""
    return _EXTENSIONS.get(os.path.splitext(path)[1].lower())


def _format_of(path: str) -> Optional[str]:
    """
    Format of a file to read: from its extension, or from its first
    bytes for unknown and ambiguous extensions; None if unknown.
    """
    extension = os.path.splitext(path)[1].lower()
    name = _EXTENSIONS.get(extension)
    if name is not None and extension not in _AMBIGUOUS:
        return name

    with open(path, "rb") as fid:
        head = fid.read(_HEAD_SIZE)
    for candidate, fmt in _FORMATS.items():
        if fmt.detect is None:
            continue
        module = importlib.import_module(fmt.module)
        if getattr(module, fmt.detect)(head):
            return candidate

    return name


def _input_files(source: Any, list_file: Any) -> list[tuple[str, bool]]:
    """
    Files of ``source`` and ``list_file``, with a flag that is True for
    files given by name (False for those of directories and patterns).
    """
    if source is None and list_file is None:
        raise ValueError("source or list_file must be given")

    files = []
    if source is not None:
        files.extend(_expand(source))
    if list_file is not None:
        files.extend(_expand(_list_entries(list_file)))

    return files


def _expand(source: Any) -> list[tuple[str, bool]]:
    """Expand a path, pattern, directory or sequence (see _input_files)."""
    if isinstance(source, (list, tuple)):
        return [item for entry in source for item in _expand(entry)]

    if not isinstance(source, (str, os.PathLike)):
        raise TypeError(
            "source must be a path, a pattern, a directory or a list of "
            f"them, got {type(source).__name__}"
        )

    path = os.path.expandvars(os.path.expanduser(os.fspath(source)))

    if os.path.isdir(path):
        return [
            (os.path.join(path, entry), False)
            for entry in sorted(os.listdir(path))
            if not entry.startswith(".")
            and os.path.isfile(os.path.join(path, entry))
        ]

    if any(char in path for char in "*?["):
        matches = [match for match in sorted(glob(path))
                   if os.path.isfile(match)]
        if not matches:
            raise FileNotFoundError(f"No file matches {path}")
        return [(match, False) for match in matches]

    if not os.path.isfile(path):
        raise FileNotFoundError(f"No such file: {path}")

    return [(path, True)]


def _list_entries(list_file: Any) -> list[str]:
    """Entries of a list file, relative to its folder."""
    if not isinstance(list_file, (str, os.PathLike)):
        raise TypeError("list_file must be a path")

    list_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(os.fspath(list_file)))
    )
    folder = os.path.dirname(list_file)
    entries = []

    with open(list_file, "r", encoding="utf-8") as fid:
        for line in fid:
            entry = line.strip()
            if not entry or entry.startswith("#"):
                continue
            entry = os.path.expandvars(os.path.expanduser(entry))
            entries.append(os.path.join(folder, entry))

    return entries


def _warn_skipped(paths: list[str]) -> None:
    """One warning for the files of unknown format."""
    names = [os.path.basename(path) for path in paths]
    listed = ", ".join(names[:5]) + (", ..." if len(names) > 5 else "")
    warn_external(
        f"{len(paths)} file(s) of unknown format skipped: {listed}",
        IoWarning,
    )


def _as_metadata(metadata: Any) -> Any:
    """
    Metadata source of :func:`reader`: collections and channel epochs
    as they are; paths and sequences read with :func:`metadata_reader`.
    """
    if metadata is None:
        return None

    from shakelab.signals.metadata import ChannelEpoch, MetadataCollection

    if isinstance(metadata, (MetadataCollection, ChannelEpoch)):
        return metadata

    if isinstance(metadata, (str, os.PathLike, list, tuple)):
        return metadata_reader(metadata)

    raise TypeError(
        "metadata must be a MetadataCollection, a ChannelEpoch, a path "
        "or a sequence of paths"
    )


def _names(records: list[Record], template: str) -> dict[str, list]:
    """Records grouped by the file names given by a template."""
    files: dict[str, list] = {}
    for record in records:
        try:
            path = template.format_map(_name_fields(record))
        except (KeyError, IndexError, ValueError):
            raise ValueError(
                f"invalid file name template {template!r}; fields: "
                + ", ".join("{" + field + "}" for field in _NAME_FIELDS)
            ) from None
        files.setdefault(path, []).append(record)

    return files


def _name_fields(record: Record) -> dict[str, Any]:
    """Values of the template fields for a record."""
    sid = record.sid or ""
    codes = sid.split(".")
    if len(codes) != 4:
        codes = ["", "", "", sid]
    fields: dict[str, Any] = dict(zip(_NAME_FIELDS, codes + [sid]))
    if record.time is not None:
        nanoseconds = record.time.to_unix_ns()
        time = datetime.fromtimestamp(
            nanoseconds // 10 ** 9, tz=timezone.utc
        ) + timedelta(microseconds=nanoseconds % 10 ** 9 // 1000)
        fields.update(
            start=time.strftime("%Y%m%dT%H%M%S"),
            time=time,
            year=f"{time.year:04d}",
            jday=f"{time.timetuple().tm_yday:03d}",
        )
    return fields


def _collection(records: list[Record]) -> StreamCollection:
    """
    Collection holding the given records (not copies, nor merged),
    grouped into streams by consecutive identifiers.
    """
    from shakelab.signals.base import Stream, StreamCollection

    collection = StreamCollection()
    for record in records:
        if not collection.stream or collection.stream[-1].sid != record.sid:
            collection.stream.append(Stream(sid=record.sid))
        collection.stream[-1].record.append(record)

    return collection


def _patterns(select: Any) -> Optional[list[SidPattern]]:
    """Patterns of ``select`` (a string or a sequence), or None."""
    if select is None:
        return None
    if isinstance(select, str):
        select = [select]
    elif not isinstance(select, (list, tuple)):
        raise TypeError("select must be a pattern or a list of patterns")
    patterns = []
    for item in select:
        patterns.extend(parse_sid_pattern(item))
    return patterns


def _window(
    starttime: Any,
    endtime: Any,
) -> tuple[Optional[Date], Optional[Date]]:
    """Bounds of the time window, as dates."""
    bounds = []
    for name, value in (("starttime", starttime), ("endtime", endtime)):
        if value is not None and not isinstance(value, (Date, str)):
            raise TypeError(
                f"{name} must be a Date or an ISO 8601 string"
            )
        bounds.append(None if value is None else Date(value))
    if None not in bounds and bounds[0] > bounds[1]:
        raise ValueError("the end of the time window precedes its start")
    return bounds[0], bounds[1]


def _selected(
    records: list[Record],
    patterns: Optional[list[SidPattern]],
    start: Optional[Date],
    end: Optional[Date],
) -> list[Record]:
    """Records of the selected channels, cut to the time window."""
    kept = []
    for record in records:
        if patterns is not None and not any(
            pattern.matches(record.sid) for pattern in patterns
        ):
            continue
        if start is not None or end is not None:
            if record.time is None:
                raise ValueError(
                    f"{record.sid}: a record without start time cannot "
                    "be cut to a time window"
                )
            if (
                (start is not None and record.endtime < start)
                or (end is not None and record.time > end)
            ):
                continue
            try:
                record.cut(start, end)
            except ValueError:
                continue    # no sample within the window
        kept.append(record)
    return kept


def _archive_files(
    source: Any,
    list_file: Any,
    layout: Any,
    patterns: Optional[list[SidPattern]],
    start: Optional[Date],
    end: Optional[Date],
) -> list[tuple[str, bool]]:
    """Files of the archives in ``source`` needed for a request."""
    if layout not in _LAYOUTS:
        raise ValueError(
            f"unknown layout {layout!r}; supported: {', '.join(_LAYOUTS)}"
        )
    if start is None or end is None:
        raise ValueError(
            f"layout={layout!r} requires starttime and endtime"
        )
    if list_file is not None:
        raise ValueError("list_file cannot be used with layout")
    if source is None:
        raise ValueError("source must be the root of the archive")

    roots = source if isinstance(source, (list, tuple)) else [source]
    files = []
    for root in roots:
        if not isinstance(root, (str, os.PathLike)):
            raise TypeError("the root of an archive must be a path")
        root = os.path.expandvars(os.path.expanduser(os.fspath(root)))
        if not os.path.isdir(root):
            raise FileNotFoundError(f"No such archive: {root}")
        files.extend(
            (path, False)
            for path in _sds_files(root, patterns, start, end)
        )
    return files


def _sds_files(
    root: str,
    patterns: Optional[list[SidPattern]],
    start: Date,
    end: Date,
) -> list[str]:
    """
    Files of an SDS archive for the channels and days requested (the
    day before ``start`` included).
    """
    if patterns is None:
        patterns = [SidPattern("*", "*", "*", "*")]

    first = datetime(start.year, start.month, start.day) - timedelta(1)
    last = datetime(end.year, end.month, end.day)
    days: dict[str, set] = {}
    day = first
    while day <= last:
        days.setdefault(f"{day.year:04d}",
                        set()).add(f"{day.timetuple().tm_yday:03d}")
        day += timedelta(1)

    def folders(path: str, names: list[str]) -> list[str]:
        if not os.path.isdir(path):
            return []
        return sorted(
            entry for entry in os.listdir(path)
            if os.path.isdir(os.path.join(path, entry))
            and any(fnmatchcase(entry, name) for name in names)
        )

    files = []
    for year, jdays in days.items():
        top = os.path.join(root, year)
        for net in folders(top, [p.network for p in patterns]):
            current = [p for p in patterns if fnmatchcase(net, p.network)]
            for sta in folders(os.path.join(top, net),
                               [p.station for p in current]):
                here = [p for p in current if fnmatchcase(sta, p.station)]
                folder = os.path.join(top, net, sta)
                for cha in folders(folder,
                                   [p.channel + ".D" for p in here]):
                    path = os.path.join(folder, cha)
                    for entry in sorted(os.listdir(path)):
                        parts = entry.split(".")
                        if (
                            len(parts) == 7
                            and parts[4] == "D"
                            and parts[5] == year
                            and parts[6] in jdays
                            and any(p.matches(".".join(parts[:4]))
                                    for p in here)
                        ):
                            files.append(os.path.join(path, entry))
    return files

