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
"""Read flat-layer velocity models from files and package resources."""

from __future__ import annotations

import csv
import json
from importlib import resources
from math import inf, isfinite
from pathlib import Path
from typing import Optional

from .ttsolver import LayeredModel1D

__all__ = ["list_models", "load_model", "model_info", "read_model"]

_CSV_FIELDS = ("thickness (m)", "vp (m/s)", "vs (m/s)")
_JSON_FIELDS = ("thickness_m", "vp_m_s", "vs_m_s")


def _models_dir():
    return resources.files(__package__).joinpath("models")


def list_models() -> tuple[str, ...]:
    """Return the packaged model identifiers in alphabetical order."""
    return tuple(sorted(
        entry.name[:-5] for entry in _models_dir().iterdir()
        if entry.name.endswith(".json") and entry.is_file()
    ))


def _read_catalog(name):
    if not isinstance(name, str):
        raise TypeError("Model name must be a string.")
    available = list_models()
    if name not in available:
        raise ValueError(
            f"Unknown model {name!r}. Available: {', '.join(available)}."
        )
    with _models_dir().joinpath(f"{name}.json").open(
        "r", encoding="utf-8",
    ) as stream:
        return _read_json(stream)


def load_model(name: str) -> LayeredModel1D:
    """Return a fresh model for an identifier from :func:`list_models`."""
    return _build_model(_json_columns(_read_catalog(name)))


def model_info(name: str) -> dict[str, str]:
    """Return the published metadata for a packaged model."""
    metadata = _read_catalog(name).get("metadata", {})
    if not isinstance(metadata, dict):
        raise ValueError("Model metadata must be a JSON object.")
    return metadata


def read_model(
    source: str | Path,
    *,
    format: Optional[str] = None,
) -> LayeredModel1D:
    """Read a user-supplied JSON or CSV file into a new layered model.

    Parameters
    ----------
    source : str or pathlib.Path
        Path to a UTF-8 model file. JSON uses a ``layers`` array with
        ``thickness_m``, ``vp_m_s`` and ``vs_m_s`` per layer. CSV uses
        the headers ``thickness (m)``, ``vp (m/s)`` and ``vs (m/s)``.
        Layers are ordered from the surface down, with thicknesses in
        metres and velocities in metres per second. The final layer
        may have infinite thickness, written as JSON ``null`` or CSV
        ``inf``. JSON metadata and additional CSV columns are optional.
    format : {'json', 'csv'}, optional
        File format. Defaults to the file extension.

    Returns
    -------
    LayeredModel1D
        Model ready for travel-time calculations.

    Raises
    ------
    ValueError
        If the file format or layer parameters are invalid.
    """
    path = Path(source)
    if format is not None and not isinstance(format, str):
        raise TypeError("Model format must be a string.")
    kind = (path.suffix[1:] if format is None else format).lower()
    if kind not in ("json", "csv"):
        raise ValueError(
            f"Unsupported model format {kind!r}. Use JSON or CSV."
        )
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        if kind == "json":
            columns = _json_columns(_read_json(stream))
        else:
            try:
                columns = _read_csv_columns(stream)
            except (TypeError, csv.Error) as exc:
                raise ValueError("Invalid CSV model row.") from exc
    return _build_model(columns)


def _read_csv_columns(stream, *, with_density=False):
    fields = _CSV_FIELDS + (("density (kg/m3)",) if with_density else ())
    reader = csv.DictReader(stream)
    columns = tuple([] for _ in fields)
    if not with_density:
        missing = set(fields).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Missing CSV model columns: {', '.join(sorted(missing))}."
            )
    for row in reader:
        for field, column in zip(fields, columns):
            column.append(float(row[field]))
    return columns


def _invalid_json_constant(value):
    raise ValueError(
        f"Invalid JSON number {value!r}. Use null for a half-space."
    )


def _read_json(stream):
    data = json.load(stream, parse_int=float,
                     parse_constant=_invalid_json_constant)
    if not isinstance(data, dict):
        raise ValueError("A JSON model must be an object.")
    return data


def _json_columns(data):
    layers = data.get("layers")
    if not isinstance(layers, list):
        raise ValueError("A JSON model must contain a 'layers' array.")

    columns = ([], [], [])
    for index, layer in enumerate(layers, 1):
        for field, column in zip(_JSON_FIELDS, columns):
            if not isinstance(layer, dict) or field not in layer:
                raise ValueError(f"Layer {index} is missing {field!r}.")
            value = layer[field]
            if field == "thickness_m" and value is None:
                column.append(inf)
            elif type(value) is float and isfinite(value):
                column.append(value)
            else:
                raise ValueError(f"Layer {index}: invalid {field!r}.")
    return columns


def _build_model(columns):
    if not columns[0]:
        raise ValueError("A model must contain at least one layer.")
    for index, (thickness, vp, vs) in enumerate(zip(*columns), 1):
        half_space = thickness == inf and index == len(columns[0])
        if not half_space and (not isfinite(thickness) or thickness <= 0):
            raise ValueError(f"Layer {index}: invalid thickness.")
        if not (isfinite(vp) and vp > 0 and isfinite(vs) and vs > 0):
            raise ValueError(
                f"Layer {index}: Vp and Vs must be finite and positive."
            )
    return LayeredModel1D(*columns)
