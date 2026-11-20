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
File formats of YARRAY.

- Array files (JSON): the stations with their local coordinates and
  the waveform file of each component, for data without (or with
  unusable) metadata; read by :meth:`ArrayData.read
  <shakelab.exploration.yarray.data.ArrayData.read>` with
  :func:`read_array_file`. Data with metadata are given to
  :class:`~shakelab.exploration.yarray.data.ArrayData` as waveforms
  read by :func:`shakelab.signals.io.reader` instead.
- Result files (JSON Lines, with optional ``.npz`` arrays): written
  and read by the ``write`` and ``read`` methods of
  :class:`~shakelab.exploration.yarray.beamforming.BeamformingResult`
  and :class:`~shakelab.exploration.yarray.fdd.FDDResult`, with
  :func:`write_result` and :func:`read_result`; tables of the results
  in CSV with :func:`write_table`.

Array file
----------
Version 1, a JSON object::

    {
      "version": 1,
      "coordinate_system": {"type": "local_cartesian", "unit": "m"},
      "reference": [13.25, 46.06, 120.0],
      "data_root": "sac",
      "reader": {"format": "sac", "options": {}},
      "stations": [
        {"id": "S1001",
         "coordinates": {"x": 2048.0, "y": 1964.0, "z": 0.0},
         "coordinate_covariance": [[0.01, 0.0], [0.0, 0.01]],
         "components": {"Z": "S1001.Z.sac", "N": "S1001.N.sac",
                        "H1": {"file": "S1001.1.sac",
                               "orientation": [0.7, 0.7, 0.0],
                               "reader": {"format": "sac"}}}}
      ]
    }

Coordinates are East, North[, Up] in metres (a list or an object
``x``, ``y``, ``z``); a station without ``coordinates`` has no position
(enough for the frequency domain decomposition, not for the
wavenumber methods). The optional ``reference`` is the geographic
position (longitude, latitude[, elevation in m]) of their origin.
Component names Z, N, E (also UP, VERTICAL, NORTH, EAST, X, Y and, as
in the ArrayTool3C file lists, 1, 2, 3 for E, N, Z) imply the
orientation; other names need an explicit
``orientation`` (east, north, up). Files are relative to ``data_root``
(default: the folder of the array file); ``reader`` gives the
arguments of :func:`shakelab.signals.io.reader` (``format``,
``byte_order``, ``options``), globally or per component. Each file
must contain one record; the components of a station are cut to their
common time window.

Result file
-----------
Version 2, JSON Lines: one JSON object per line, in UTF-8. The first
line is the header: ``format`` (``"shakelab.yarray.result"``),
``version``, ``kind`` (``"beamforming"`` or ``"fdd"``), ``config``
(JSON representation of the configuration), ``frequencies`` (Hz) and
``units``; for the FDD also the channels, ``stations``,
``components`` and ``coordinates`` (East, North, Up in metres, null
when unknown). Each following line is one result written (one
analysis, or several joined); a file written with ``append=True``
holds the results of successive analyses with the same header, joined
when read. Lists follow the order of ``frequencies``; times are ISO
8601 with nanoseconds; missing values are null.

- Beamforming: ``block_times`` (by frequency), ``station_sets``
  (``stations`` and ``coordinates`` of each set of stations),
  ``block_stations`` (index of the set of each block, by frequency),
  ``maxima`` (table by column: ``block``, ``time``, ``frequency``,
  ``component``, ``rank``, ``wavenumber`` in rad/m, ``azimuth`` in
  degrees, ``velocity`` and ``slowness``, derived, ``power``,
  ``power_standard_error``, ``ellipticity``, ``ellipticity_angle``).
- FDD: ``block_times``, ``values`` (first singular value of each
  block) and ``vectors`` (``real`` and ``imag`` parts of the first
  singular vector of each block, by channel), by frequency.

``arrays`` is null or the name, relative to the folder of the file,
of a ``.npz`` file (folder ``<name>_arrays``, files named after the
first block time) with the arrays of the line: the spectra of the
beamforming (``spectra``, a JSON table, and ``spectrum<i>_<name>``),
or all the singular values (``values_<i>``, ``(n_blocks,
n_channels)``) and the kept singular vectors (``vectors_<i>``,
``(n_blocks, n_channels, n_vectors)``, complex) of the FDD, ``i``
being the index of the frequency. The ``.npz`` files are read without
pickle.
"""

from __future__ import annotations

import csv
import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional

import numpy as np

from shakelab.exploration.yarray.data import StationData
from shakelab.libutils.timeN import Date
from shakelab.libutils.utils import check_target
from shakelab.signals import io as waveform_io

__all__ = [
    "read_array_file",
    "read_result",
    "read_result_header",
    "result_path",
    "write_result",
    "write_table",
]

# Format name and version of the result files.
_RESULT_FORMAT = "shakelab.yarray.result"
_RESULT_VERSION = 2

# Version of the array file format.
_VERSION = 1

# Orientations implied by component names, and their order.
_CANONICAL = {
    **dict.fromkeys(("Z", "UP", "VERTICAL", "3"), (0.0, 0.0, 1.0)),
    **dict.fromkeys(("N", "NORTH", "Y", "2"), (0.0, 1.0, 0.0)),
    **dict.fromkeys(("E", "EAST", "X", "1"), (1.0, 0.0, 0.0)),
}
_ORDER = ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0))


# ----------------------------------------------------------------------------
# Array files
# ----------------------------------------------------------------------------


def read_array_file(source: Any) -> tuple[list[StationData], Any]:
    """
    Read the stations and the reference of an array file (JSON).

    Parameters
    ----------
    source : str or path-like
        Array file (see the module notes for the format).

    Returns
    -------
    stations : list of StationData
        Stations in local coordinates, each on the time window common
        to its components.
    reference : tuple of float or None
        Geographic position of the origin of the coordinates, if
        given.

    Raises
    ------
    FileNotFoundError
        If the array file or a waveform file does not exist.
    TypeError
        If a field has the wrong type.
    ValueError
        If the file is not valid JSON, the version or the coordinate
        system is not supported, a field is invalid, a file does not
        contain exactly one record, or the components of a station
        have no common time window.
    """
    path = Path(source).expanduser()

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"{path}: invalid JSON at line {error.lineno}, column "
            f"{error.colno}: {error.msg}"
        ) from error

    if not isinstance(data, Mapping):
        raise TypeError(f"{path}: the array file must be a JSON object")

    if data.get("version", _VERSION) != _VERSION:
        raise ValueError(f"{path}: unsupported version {data['version']}")

    _check_coordinate_system(data.get("coordinate_system"))
    root = Path(data.get("data_root", ".")).expanduser()

    if not root.is_absolute():
        root = path.parent / root

    default = _reader_options(data.get("reader"), "reader")
    entries = data.get("stations")

    if not isinstance(entries, list) or not entries:
        raise ValueError(f"{path}: 'stations' must be a non-empty list")

    reference = data.get("reference")

    if reference is not None:
        if (
            not isinstance(reference, list)
            or len(reference) not in (2, 3)
            or not all(isinstance(value, (int, float))
                       and not isinstance(value, bool)
                       for value in reference)
        ):
            raise ValueError(
                f"{path}: 'reference' must be [longitude, latitude"
                "[, elevation]]"
            )

        reference = tuple(float(value) for value in reference)

    stations = [
        _read_station(entry, f"stations[{index}]", root, default)
        for index, entry in enumerate(entries)
    ]
    return stations, reference


def _check_coordinate_system(value: Any) -> None:
    """Accept only local Cartesian coordinates in metres."""
    if value is None:
        return

    if not isinstance(value, Mapping):
        raise TypeError("'coordinate_system' must be an object")

    if value.get("type", "local_cartesian") != "local_cartesian":
        raise ValueError("only local Cartesian coordinates are supported")

    if value.get("unit", "m") != "m":
        raise ValueError("coordinates must be in metres")


def _reader_options(value: Any, context: str) -> dict:
    """Keyword arguments of the waveform reader."""
    if value is None:
        return {}

    if not isinstance(value, Mapping):
        raise TypeError(f"{context} must be an object")

    options = value.get("options", {})

    if not isinstance(options, Mapping):
        raise TypeError(f"{context}.options must be an object")

    arguments = dict(options)

    for key in ("format", "byte_order"):
        if value.get(key) is not None:
            arguments[key] = value[key]

    return arguments


def _read_station(
    entry: Any,
    context: str,
    root: Path,
    default: dict,
) -> StationData:
    """Read the components of one station."""
    if not isinstance(entry, Mapping):
        raise TypeError(f"{context} must be an object")

    sid = entry.get("id")

    if not isinstance(sid, str) or not sid.strip():
        raise ValueError(f"{context}.id must be a non-empty string")

    components = entry.get("components")

    if not isinstance(components, Mapping) or not 1 <= len(components) <= 3:
        raise ValueError(f"{context}.components must have 1 to 3 entries")

    items = []

    for name, spec in components.items():
        item = f"{context}.components.{name}"
        key = str(name).strip().upper()

        if isinstance(spec, str):
            spec = {"file": spec}

        if not isinstance(spec, Mapping) or not isinstance(
            spec.get("file"), str
        ):
            raise TypeError(f"{item} must be a file name or an object")

        orientation = spec.get("orientation", _CANONICAL.get(key))

        if orientation is None:
            raise ValueError(
                f"{item}: the orientation of component {name!r} must be "
                "given"
            )

        options = {**default, **_reader_options(spec.get("reader"), item)}
        record = _read_record(root / spec["file"], options, item)
        order = (
            _ORDER.index(tuple(orientation))
            if tuple(orientation) in _ORDER else len(_ORDER)
        )
        items.append((order, key, record, orientation))

    items.sort(key=lambda item: item[:2])
    records = [item[2] for item in items]
    start = max(rec.starttime for rec in records)
    end = min(rec.endtime for rec in records)

    if end < start:
        raise ValueError(f"{context}: the components have no common time")

    return StationData(
        sid=sid.strip(),
        records=[rec.extract(start, end) for rec in records],
        coordinates=(
            None if entry.get("coordinates") is None
            else _coordinates(entry["coordinates"], context)
        ),
        orientation=[item[3] for item in items],
        coordinate_covariance=entry.get("coordinate_covariance"),
    )


def _coordinates(value: Any, context: str) -> list:
    """Coordinates from a list or an object x, y[, z]."""
    if isinstance(value, Mapping):
        if "x" not in value or "y" not in value:
            raise ValueError(f"{context}.coordinates needs x and y")

        return [value["x"], value["y"]] + (
            [value["z"]] if "z" in value else []
        )

    if isinstance(value, list):
        return value

    raise TypeError(f"{context}.coordinates must be a list or an object")


def _read_record(path: Path, options: dict, context: str) -> Any:
    """Read the single record of a waveform file."""
    if not path.is_file():
        raise FileNotFoundError(f"{context}: file not found: {path}")

    records = waveform_io.reader(str(path), **options).records()

    if len(records) != 1:
        raise ValueError(
            f"{context}: {path} contains {len(records)} records, one is "
            "required"
        )

    return records[0]


# ----------------------------------------------------------------------------
# Result files
# ----------------------------------------------------------------------------


def write_result(
    header: dict,
    entry: dict,
    target: Any,
    arrays: Optional[Mapping] = None,
    times: Optional[Mapping] = None,
    append: bool = False,
    owrite: bool = False,
) -> Path:
    """
    Write one analysis to a result file (JSON Lines).

    Parameters
    ----------
    header : dict
        Header of the file (``kind``, ``config``, ``units``...),
        written as the first line of a new file; with ``append`` it
        must equal the header of the file.
    entry : dict
        Content of the analysis, written as one line.
    target : str or path-like
        File; folders are created. A name ending in ``.json`` is
        refused (the file is not a JSON document).
    arrays : mapping, optional
        Arrays written to a ``.npz`` file next to the result file
        (folder ``<name>_arrays``), referred to by the line.
    times : mapping, optional
        Block times (Date) by frequency: the ``.npz`` file is named
        after the first one (``20261116T101500.npz``; a number is
        added if it exists in a file being appended to).
    append : bool, optional
        Add the line to the file (created if missing). Default False.
    owrite : bool, optional
        Overwrite an existing file (ignored with ``append``). Default
        is False.

    Returns
    -------
    pathlib.Path
        Path of the result file.

    Raises
    ------
    ValueError
        If the name of the file ends in ``.json``, or ``append`` and
        the file has another header.
    FileExistsError
        If the file exists and neither ``append`` nor ``owrite`` is set.
    """
    path = result_path(target)

    check_target(path, owrite, append=append)

    exists = path.exists() and path.stat().st_size > 0

    if exists and append:
        current = read_result_header(path)

        if current.get("kind") != header.get("kind"):
            raise ValueError(
                f"{path}: the file holds {current.get('kind')} results"
            )

        if current != json.loads(_dumps(header)):
            raise ValueError(
                f"{path}: the file holds results with another "
                "configuration, other channels, frequencies or units"
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    line = dict(entry)
    line["arrays"] = None

    if arrays:
        folder = path.with_name(path.stem + "_arrays")
        folder.mkdir(exist_ok=True)
        name = _time_name(times or {})
        file = folder / f"{name}.npz"
        count = 1

        while append and file.exists():
            file = folder / f"{name}_{count}.npz"
            count += 1

        np.savez_compressed(file, **{
            key: np.asarray(value) for key, value in arrays.items()
        })
        line["arrays"] = f"{folder.name}/{file.name}"

    new = not (append and exists)
    text = _dumps(line)

    with path.open("w" if new else "a", encoding="utf-8") as stream:
        if new:
            stream.write(_dumps({
                "format": _RESULT_FORMAT, "version": _RESULT_VERSION,
                **header,
            }))

        stream.write(text)

    return path


def read_result_header(source: Any) -> dict:
    """
    Read the header of a result file.

    Parameters
    ----------
    source : str or path-like
        Result file (JSON Lines).

    Returns
    -------
    dict
        Header, without ``format`` and ``version``.

    Raises
    ------
    ValueError
        If the file is not a YARRAY result file or has another version
        of the format.
    """
    path = Path(source).expanduser()

    with path.open(encoding="utf-8") as stream:
        return _header(path, stream.readline())


def read_result(source: Any, arrays: bool = True) -> tuple[dict, list]:
    """
    Read a result file.

    Parameters
    ----------
    source : str or path-like
        Result file (JSON Lines) written by :func:`write_result`.
    arrays : bool, optional
        Load the ``.npz`` files of the analyses. Default True.

    Returns
    -------
    header : dict
        Header, without ``format`` and ``version``.
    entries : list of dict
        Analyses in the order they were written; ``arrays`` is a
        dictionary of arrays, or None if the analysis has none or
        ``arrays`` is False.

    Raises
    ------
    ValueError
        If the file is not a YARRAY result file, has another version of
        the format, contains no analysis or an invalid line.
    FileNotFoundError
        If a ``.npz`` file of an analysis is missing (read with
        ``arrays=False`` to ignore them).
    """
    path = Path(source).expanduser()

    with path.open(encoding="utf-8") as stream:
        header = _header(path, stream.readline())
        entries = []

        for number, text in enumerate(stream, start=2):
            if not text.strip():
                continue

            try:
                entry = json.loads(text)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path}: invalid line {number}: {error.msg}"
                ) from error

            if not isinstance(entry, dict):
                raise ValueError(f"{path}: invalid line {number}")

            entries.append(entry)

    if not entries:
        raise ValueError(f"{path}: the file holds no analysis")

    for entry in entries:
        name = entry.get("arrays")

        if not arrays or name is None:
            entry["arrays"] = None
            continue

        file = path.parent / name

        if not file.is_file():
            raise FileNotFoundError(
                f"{path}: the arrays {file} are missing (read with "
                "arrays=False to ignore them)"
            )

        with np.load(file, allow_pickle=False) as data:
            entry["arrays"] = {key: data[key] for key in data.files}

    return header, entries


def result_path(target: Any) -> Path:
    """
    Return the path of a result file, refusing the suffix ``.json``.

    Parameters
    ----------
    target : str or path-like
        Name of the file.

    Returns
    -------
    pathlib.Path
        Path, with ``~`` expanded.

    Raises
    ------
    ValueError
        If the name ends in ``.json``.
    """
    path = Path(target).expanduser()

    if path.suffix.lower() == ".json":
        raise ValueError(
            f"{path}: result files are JSON Lines, not JSON documents: "
            "use the suffix .jsonl"
        )

    return path


def write_table(target: Any, columns: Mapping, owrite: bool = False) -> None:
    """
    Write a table, given by columns, to a CSV file.

    Parameters
    ----------
    target : str or path-like
        CSV file; folders are created.
    columns : mapping
        Values of each column (name -> sequence), all of the same
        length. Missing values (None, NaN) are empty fields; times are
        written in ISO 8601 with nanoseconds.
    owrite : bool, optional
        Overwrite an existing file. Default False.

    Raises
    ------
    FileExistsError
        If the file exists and ``owrite`` is False.
    """
    path = Path(target).expanduser()

    check_target(path, owrite)

    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)

    def text(value: Any) -> str:
        value = _plain(value)

        if value is None:
            return ""

        return repr(value) if isinstance(value, float) else str(value)

    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)

        for row in zip(*(columns[name] for name in names)):
            writer.writerow([text(value) for value in row])


def _time_name(times: Mapping) -> str:
    """Name of a .npz file: the first block time (or ``result``)."""
    values = [time for items in times.values() for time in items]

    if not values:
        return "result"

    text = min(values).to_iso8601(precision=0)
    return text.replace("-", "").replace(":", "").rstrip("Z")


def _header(path: Path, text: str) -> dict:
    """Check and return the header line of a result file."""
    try:
        header = json.loads(text)
    except json.JSONDecodeError:
        header = None

    if not isinstance(header, dict) or (
        header.get("format") != _RESULT_FORMAT
    ):
        raise ValueError(f"{path}: not a YARRAY result file")

    if header.get("version") != _RESULT_VERSION:
        raise ValueError(
            f"{path}: unsupported result version {header.get('version')}"
        )

    return {key: value for key, value in header.items()
            if key not in ("format", "version")}


def _dumps(value: Any) -> str:
    """One line of JSON (NaN and infinities as null)."""
    return json.dumps(_plain(value), allow_nan=False,
                      separators=(",", ":")) + "\n"


def _plain(value: Any) -> Any:
    """Value with JSON types only."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}

    if isinstance(value, np.ndarray):
        return _plain(value.tolist())

    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]

    if isinstance(value, Date):
        return value.to_iso8601(precision=9)

    if isinstance(value, (bool, np.bool_)):
        return bool(value)

    if isinstance(value, (int, np.integer)):
        return int(value)

    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if math.isfinite(value) else None

    return value
