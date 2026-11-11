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
Input and output of YArray.

- :func:`read_array`: an :class:`~.data.ArrayData` from an array file
  (JSON), which lists the stations with their local coordinates and
  the waveform file of each component, for data without (or with
  unusable) metadata. Data with metadata are read with
  :func:`shakelab.signals.io.reader` and
  :meth:`~.data.ArrayData.from_collection` instead.
- :func:`write_result`, :func:`read_results`, :func:`read_result`:
  beamforming results in pickle files, also accumulated in one file
  (``append=True``).

Array file
----------
Version 1, a JSON object::

    {
      "version": 1,
      "coordinate_system": {"type": "local_cartesian", "unit": "m"},
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
``x``, ``y``, ``z``). Component names Z, N, E (also UP, VERTICAL,
NORTH, EAST, X, Y and, as in the ArrayTool3C file lists, 1, 2, 3 for
E, N, Z) imply the orientation; other names need an explicit
``orientation`` (east, north, up). Files are relative to ``data_root``
(default: the folder of the array file); ``reader`` gives the
arguments of :func:`shakelab.signals.io.reader` (``format``,
``byte_order``, ``options``), globally or per component. Each file
must contain one record; the records are cut to the time window common
to all of them.

Warnings
--------
Pickle files can run arbitrary code when read: read only result files
of trusted origin. Pickle is the provisional format of the results,
until the native ShakeLab format.
"""

from __future__ import annotations

import json
import pickle
from collections.abc import Mapping
from pathlib import Path
from typing import Any


from shakelab.exploration.yarray.beamforming import BeamformingResult
from shakelab.exploration.yarray.data import ArrayData, StationData
from shakelab.signals import io as waveform_io

__all__ = [
    "read_array",
    "read_result",
    "read_results",
    "write_result",
]

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


def read_array(source: Any) -> ArrayData:
    """
    Read an array from an array file (JSON).

    Parameters
    ----------
    source : str or path-like
        Array file (see the module notes for the format).

    Returns
    -------
    ArrayData
        Stations on the time window common to all records, in local
        coordinates (no geographic frame).

    Raises
    ------
    FileNotFoundError
        If the array file or a waveform file does not exist.
    TypeError
        If a field has the wrong type.
    ValueError
        If the file is not valid JSON, the version or the coordinate
        system is not supported, a field is invalid, a file does not
        contain exactly one record, or the records have no common
        time window.
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

    stations = [
        _read_station(entry, f"stations[{index}]", root, default)
        for index, entry in enumerate(entries)
    ]

    return _common_window(stations)


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

    return StationData(
        sid=sid.strip(),
        records=[item[2] for item in items],
        coordinates=_coordinates(entry.get("coordinates"), context),
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


def _common_window(stations: list) -> ArrayData:
    """Cut every record to the window common to all of them."""
    records = [rec for station in stations for rec in station.records]
    start = max(rec.starttime for rec in records)
    end = min(rec.endtime for rec in records)

    if end < start:
        raise ValueError("the records have no common time window")

    return ArrayData([
        StationData(
            station.sid,
            [rec.extract(start, end) for rec in station.records],
            station.coordinates,
            station.orientation,
            station.coordinate_covariance,
        )
        for station in stations
    ])


# ----------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------


def write_result(
    result: BeamformingResult,
    target: Any,
    append: bool = False,
    owrite: bool = False,
) -> None:
    """
    Write a beamforming result to a pickle file.

    Parameters
    ----------
    result : BeamformingResult
        Result to write.
    target : str or path-like
        File; folders are created.
    append : bool, optional
        Add the result to the end of the file (created if missing), to
        accumulate the results of successive analyses. Default False.
    owrite : bool, optional
        Overwrite an existing file (ignored with ``append``). Default
        is False.

    Raises
    ------
    TypeError
        If ``result`` is not a BeamformingResult.
    FileExistsError
        If the file exists and neither ``append`` nor ``owrite`` is set.
    """
    if not isinstance(result, BeamformingResult):
        raise TypeError("result must be a BeamformingResult")

    path = Path(target).expanduser()

    if path.exists() and not (append or owrite):
        raise FileExistsError(f"{path} exists (use owrite or append)")

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("ab" if append else "wb") as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)


def read_results(source: Any) -> list[BeamformingResult]:
    """
    Read all the results of a pickle file.

    Parameters
    ----------
    source : str or path-like
        File written by :func:`write_result` (trusted origin only).

    Returns
    -------
    list of BeamformingResult
        Results in the order they were written.

    Raises
    ------
    ValueError
        If the file contains no result or other objects.
    """
    results = []

    with Path(source).expanduser().open("rb") as stream:
        while True:
            try:
                item = pickle.load(stream)
            except EOFError:
                break

            if not isinstance(item, BeamformingResult):
                raise ValueError(
                    f"{source}: unexpected object {type(item).__name__}"
                )

            results.append(item)

    if not results:
        raise ValueError(f"{source}: no result")

    return results


def read_result(source: Any, concatenate: bool = False) -> BeamformingResult:
    """
    Read the result of a pickle file.

    Parameters
    ----------
    source : str or path-like
        File written by :func:`write_result` (trusted origin only).
    concatenate : bool, optional
        If the file holds several results, join them
        (:meth:`BeamformingResult.concatenate`) instead of raising.

    Returns
    -------
    BeamformingResult
        The result.

    Raises
    ------
    ValueError
        If the file holds several results and ``concatenate`` is False.
    """
    results = read_results(source)

    if len(results) == 1:
        return results[0]

    if not concatenate:
        raise ValueError(
            f"{source} holds {len(results)} results: use read_results or "
            "concatenate=True"
        )

    return BeamformingResult.concatenate(results)
