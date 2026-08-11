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
Input/output helpers for seismic-array processing.

The module separates external data access from numerical processing.

``read_from_files`` builds an :class:`ArrayData` snapshot from a JSON
configuration and waveform files. A future ``read_from_buffer`` helper
will build the same snapshot type from a circular live-data buffer.

Processing configurations support direct arbitrary-frequency Fourier
evaluation and fixed-window RFFT processing with nearest-bin or complex
interpolated frequency selection.

Beamforming results can currently be serialized with Python pickle.
This format is intended for development and trusted local workflows;
it is not yet the final public interchange format.
"""

from __future__ import annotations

import json
import pickle
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from shakelab.signals import io as waveform_io
from shakelab.signals.base import Record

from .beamforming import (
    BeamformingConfig,
    BeamformingMethod,
    BeamformingResult,
    ComponentMode,
    EllipticityGrid,
    GeometricUncertainty,
    GridQuantity,
    GridScale,
    OutputMode,
    SearchGrid,
    StatisticalUncertainty,
)
from .spectral import (
    FFTSelection,
    FrequencyConfig,
    FrequencyScale,
    TransformMethod,
    WindowConfig,
    WindowType,
)
from .data import ArrayData, StationData


FloatArray = NDArray[np.float64]


_CANONICAL_ORIENTATIONS = {
    "E": np.array([1.0, 0.0, 0.0]),
    "EAST": np.array([1.0, 0.0, 0.0]),
    "X": np.array([1.0, 0.0, 0.0]),
    "1": np.array([1.0, 0.0, 0.0]),
    "N": np.array([0.0, 1.0, 0.0]),
    "NORTH": np.array([0.0, 1.0, 0.0]),
    "Y": np.array([0.0, 1.0, 0.0]),
    "2": np.array([0.0, 1.0, 0.0]),
    "Z": np.array([0.0, 0.0, 1.0]),
    "VERTICAL": np.array([0.0, 0.0, 1.0]),
    "UP": np.array([0.0, 0.0, 1.0]),
    "3": np.array([0.0, 0.0, 1.0]),
}

_COMPONENT_ORDER = {
    "E": 0,
    "EAST": 0,
    "X": 0,
    "1": 0,
    "N": 1,
    "NORTH": 1,
    "Y": 1,
    "2": 1,
    "Z": 2,
    "VERTICAL": 2,
    "UP": 2,
    "3": 2,
}


def read_from_files(
    config_file: str | Path,
) -> ArrayData:
    """
    Read an array-processing snapshot from waveform files.

    Parameters
    ----------
    config_file
        Path to the array JSON configuration.

    Returns
    -------
    ArrayData
        Validated array-data snapshot.

    Raises
    ------
    FileNotFoundError
        If the configuration or a waveform file does not exist.
    TypeError
        If the JSON structure contains invalid value types.
    ValueError
        If the configuration is invalid or a waveform file does not
        resolve to exactly one record.
    """
    config_path = _resolve_config_path(config_file)
    config = read_array_config(config_path)

    _validate_version(config)

    data_root = _resolve_data_root(
        config,
        config_path,
    )

    reader_config = _parse_reader_config(
        config.get("reader")
    )

    station_entries = config.get("stations")

    if not isinstance(station_entries, list):
        raise TypeError(
            "Configuration field 'stations' must be a list."
        )

    if not station_entries:
        raise ValueError(
            "Configuration must contain at least one station."
        )

    stations = tuple(
        _load_station(
            station_entry,
            data_root=data_root,
            default_reader=reader_config,
            station_index=index,
        )
        for index, station_entry in enumerate(station_entries)
    )

    return ArrayData(stations)




def read_processing_config(
    file_path: str | Path,
) -> BeamformingConfig:
    """
    Read a beamforming-processing configuration from JSON.

    The JSON file is converted to the same :class:`BeamformingConfig`
    used by the Python API and command-line interface.
    """
    path = Path(file_path).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"Processing configuration not found: {path}"
        )

    try:
        with path.open(
            "r",
            encoding="utf-8",
        ) as stream:
            data = json.load(stream)

    except json.JSONDecodeError as error:
        raise ValueError(
            "Invalid processing JSON configuration at "
            f"line {error.lineno}, column {error.colno}: "
            f"{error.msg}"
        ) from error

    if not isinstance(data, Mapping):
        raise TypeError(
            "The processing configuration root must be a JSON object."
        )

    return processing_config_from_dict(
        data
    )


def write_processing_config(
    config: BeamformingConfig,
    file_path: str | Path,
) -> None:
    """
    Write a beamforming-processing configuration as readable JSON.
    """
    if not isinstance(config, BeamformingConfig):
        raise TypeError(
            "config must be a BeamformingConfig instance."
        )

    path = Path(file_path).expanduser().resolve()

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
    ) as stream:
        json.dump(
            processing_config_to_dict(
                config
            ),
            stream,
            indent=2,
        )

        stream.write("\n")


def processing_config_from_dict(
    data: Mapping[str, Any],
) -> BeamformingConfig:
    """
    Build :class:`BeamformingConfig` from a mapping.

    The JSON schema groups spectral-transform settings under the
    optional ``spectral`` object while retaining the established
    top-level ``frequencies`` and ``window`` sections.

    Existing configuration files without ``spectral`` remain valid
    and use direct Fourier evaluation.
    """
    if not isinstance(data, Mapping):
        raise TypeError(
            "Processing configuration must be a mapping."
        )

    version = data.get(
        "version",
        1,
    )

    if isinstance(version, bool) or not isinstance(version, int):
        raise TypeError(
            "Processing configuration field 'version' must be "
            "an integer."
        )

    if version != 1:
        raise ValueError(
            "Unsupported processing configuration version: "
            f"{version}"
        )

    frequencies = _parse_processing_frequencies(
        data.get("frequencies")
    )

    grid_data = _require_mapping(
        data,
        "grid",
    )

    window_data = _require_mapping(
        data,
        "window",
    )

    spectral_data = data.get(
        "spectral",
        {},
    )

    if not isinstance(spectral_data, Mapping):
        raise TypeError(
            "Processing field 'spectral' must be an object."
        )

    method_data = data.get(
        "method",
        {},
    )

    if isinstance(method_data, str):
        method_name = method_data
        method_options: Mapping[str, Any] = {}

    elif isinstance(method_data, Mapping):
        method_name = method_data.get(
            "name",
            BeamformingMethod.CAPON.value,
        )
        method_options = method_data

    else:
        raise TypeError(
            "Processing field 'method' must be a string or object."
        )

    output_data = data.get(
        "output",
        {},
    )

    if not isinstance(output_data, Mapping):
        raise TypeError(
            "Processing field 'output' must be an object."
        )

    uncertainty_data = data.get(
        "uncertainty",
        {},
    )

    if not isinstance(uncertainty_data, Mapping):
        raise TypeError(
            "Processing field 'uncertainty' must be an object."
        )

    ellipticity_data = data.get(
        "ellipticity",
        {},
    )

    if not isinstance(ellipticity_data, Mapping):
        raise TypeError(
            "Processing field 'ellipticity' must be an object."
        )

    peak_neighborhood = output_data.get(
        "neighborhood",
        [3, 3],
    )

    if (
        not isinstance(peak_neighborhood, Sequence)
        or isinstance(
            peak_neighborhood,
            (str, bytes),
        )
        or len(peak_neighborhood) != 2
    ):
        raise TypeError(
            "output.neighborhood must contain exactly two integers."
        )

    fft_padding = spectral_data.get(
        "fft_padding"
    )

    if fft_padding is not None:
        fft_padding = _as_integer(
            fft_padding,
            context="spectral.fft_padding",
        )

    return BeamformingConfig(
        frequencies=frequencies,
        grid=SearchGrid(
            minimum=_required_number(
                grid_data,
                "minimum",
            ),
            maximum=_required_number(
                grid_data,
                "maximum",
            ),
            size=_required_integer(
                grid_data,
                "size",
            ),
            quantity=GridQuantity(
                grid_data.get(
                    "quantity",
                    GridQuantity.VELOCITY.value,
                )
            ),
            scale=GridScale(
                grid_data.get(
                    "scale",
                    GridScale.LOG.value,
                )
            ),
            azimuth_size=_optional_integer(
                grid_data,
                "azimuth_size",
                180,
            ),
        ),
        window=WindowConfig(
            type=WindowType(
                window_data.get(
                    "type",
                    WindowType.CYCLES.value,
                )
            ),
            length=_optional_number(
                window_data,
                "length",
                20.0,
            ),
            taper=_optional_number(
                window_data,
                "taper",
                0.1,
            ),
            overlap=_optional_number(
                window_data,
                "overlap",
                0.5,
            ),
            windows_per_block=_optional_integer(
                window_data,
                "windows_per_block",
                40,
            ),
            block_overlap=_optional_number(
                window_data,
                "block_overlap",
                0.5,
            ),
        ),
        method=BeamformingMethod(
            method_name
        ),
        components=ComponentMode(
            data.get(
                "components",
                ComponentMode.VERTICAL.value,
            )
        ),
        diagonal_loading=_optional_number(
            method_options,
            "diagonal_loading",
            1e-3,
        ),
        music_sources=_optional_integer(
            method_options,
            "music_sources",
            1,
        ),
        ellipticity_grid=EllipticityGrid(
            minimum_angle=_optional_number(
                ellipticity_data,
                "minimum_angle",
                -85.0,
            ),
            maximum_angle=_optional_number(
                ellipticity_data,
                "maximum_angle",
                85.0,
            ),
            size=_optional_integer(
                ellipticity_data,
                "size",
                61,
            ),
        ),
        output=OutputMode(
            output_data.get(
                "mode",
                OutputMode.MAXIMA.value,
            )
        ),
        maxima=_optional_integer(
            output_data,
            "maxima",
            1,
        ),
        peak_neighborhood=(
            int(peak_neighborhood[0]),
            int(peak_neighborhood[1]),
        ),
        minimum_relative_power=_optional_number(
            output_data,
            "minimum_relative_power",
            0.0,
        ),
        statistical_uncertainty=StatisticalUncertainty(
            uncertainty_data.get(
                "statistical",
                StatisticalUncertainty.NONE.value,
            )
        ),
        geometric_uncertainty=GeometricUncertainty(
            uncertainty_data.get(
                "geometric",
                GeometricUncertainty.NONE.value,
            )
        ),
        strict=_optional_boolean(
            data,
            "strict",
            True,
        ),
        transform=TransformMethod(
            spectral_data.get(
                "transform",
                TransformMethod.DIRECT.value,
            )
        ),
        fft_selection=FFTSelection(
            spectral_data.get(
                "fft_selection",
                FFTSelection.NEAREST.value,
            )
        ),
        fft_padding=fft_padding,
    )


def processing_config_to_dict(
    config: BeamformingConfig,
) -> dict[str, Any]:
    """
    Convert :class:`BeamformingConfig` to the JSON-compatible schema.
    """
    if not isinstance(config, BeamformingConfig):
        raise TypeError(
            "config must be a BeamformingConfig instance."
        )

    return {
        "version": 1,
        "frequencies": _processing_frequencies_to_dict(
            config.frequencies
        ),
        "spectral": {
            "transform": config.transform.value,
            "fft_selection": config.fft_selection.value,
            "fft_padding": config.fft_padding,
        },
        "method": {
            "name": config.method.value,
            "diagonal_loading": config.diagonal_loading,
            "music_sources": config.music_sources,
        },
        "components": config.components.value,
        "grid": {
            "quantity": config.grid.quantity.value,
            "scale": config.grid.scale.value,
            "minimum": config.grid.minimum,
            "maximum": config.grid.maximum,
            "size": config.grid.size,
            "azimuth_size": config.grid.azimuth_size,
        },
        "window": {
            "type": config.window.type.value,
            "length": config.window.length,
            "taper": config.window.taper,
            "overlap": config.window.overlap,
            "windows_per_block": (
                config.window.windows_per_block
            ),
            "block_overlap": config.window.block_overlap,
        },
        "ellipticity": {
            "minimum_angle": (
                config.ellipticity_grid.minimum_angle
            ),
            "maximum_angle": (
                config.ellipticity_grid.maximum_angle
            ),
            "size": config.ellipticity_grid.size,
        },
        "output": {
            "mode": config.output.value,
            "maxima": config.maxima,
            "neighborhood": list(
                config.peak_neighborhood
            ),
            "minimum_relative_power": (
                config.minimum_relative_power
            ),
        },
        "uncertainty": {
            "statistical": (
                config.statistical_uncertainty.value
            ),
            "geometric": (
                config.geometric_uncertainty.value
            ),
        },
        "strict": config.strict,
    }


def _parse_processing_frequencies(
    value: Any,
) -> FrequencyConfig:
    """
    Parse explicit, sampled-range, or native-bin frequency settings.
    """
    if not isinstance(value, Mapping):
        raise TypeError(
            "Processing field 'frequencies' must be an object."
        )

    if "values" in value:
        values = value["values"]

        if (
            not isinstance(values, Sequence)
            or isinstance(values, (str, bytes))
        ):
            raise TypeError(
                "frequencies.values must be an array."
            )

        return FrequencyConfig.from_values(
            np.asarray(
                values,
                dtype=float,
            )
        )

    minimum = _required_number(
        value,
        "minimum",
    )

    maximum = _required_number(
        value,
        "maximum",
    )

    size_value = value.get(
        "size"
    )

    if size_value is None:
        size = None
    else:
        size = _as_integer(
            size_value,
            context="frequencies.size",
        )

    scale = FrequencyScale(
        value.get(
            "scale",
            FrequencyScale.LOG.value,
        )
    )

    return FrequencyConfig(
        minimum=minimum,
        maximum=maximum,
        size=size,
        scale=scale,
    )


def _processing_frequencies_to_dict(
    frequencies: Any,
) -> dict[str, Any]:
    """Convert a beamforming frequency definition to JSON data."""
    if isinstance(
        frequencies,
        FrequencyConfig,
    ):
        if frequencies.values is not None:
            return {
                "values": [
                    float(value)
                    for value in frequencies.values
                ],
            }

        return {
            "minimum": float(
                frequencies.minimum
            ),
            "maximum": float(
                frequencies.maximum
            ),
            "size": frequencies.size,
            "scale": frequencies.scale.value,
        }

    values = np.asarray(
        frequencies,
        dtype=float,
    )

    return {
        "values": [
            float(value)
            for value in values
        ],
    }


def _require_mapping(
    data: Mapping[str, Any],
    key: str,
) -> Mapping[str, Any]:
    """Return a required nested mapping."""
    value = data.get(
        key
    )

    if not isinstance(value, Mapping):
        raise TypeError(
            f"Processing field {key!r} must be an object."
        )

    return value


def _required_number(
    data: Mapping[str, Any],
    key: str,
) -> float:
    """Read a required finite numerical field."""
    if key not in data:
        raise ValueError(
            f"Missing required field: {key}"
        )

    return _as_finite_float(
        data[key],
        context=key,
    )


def _optional_number(
    data: Mapping[str, Any],
    key: str,
    default: float,
) -> float:
    """Read an optional finite numerical field."""
    if key not in data:
        return float(
            default
        )

    return _as_finite_float(
        data[key],
        context=key,
    )


def _required_integer(
    data: Mapping[str, Any],
    key: str,
) -> int:
    """Read a required integer field."""
    if key not in data:
        raise ValueError(
            f"Missing required field: {key}"
        )

    return _as_integer(
        data[key],
        context=key,
    )


def _optional_integer(
    data: Mapping[str, Any],
    key: str,
    default: int,
) -> int:
    """Read an optional integer field."""
    if key not in data:
        return int(
            default
        )

    return _as_integer(
        data[key],
        context=key,
    )


def _optional_boolean(
    data: Mapping[str, Any],
    key: str,
    default: bool,
) -> bool:
    """Read an optional Boolean field."""
    value = data.get(
        key,
        default,
    )

    if not isinstance(value, bool):
        raise TypeError(
            f"{key} must be a Boolean."
        )

    return value


def _as_integer(
    value: Any,
    *,
    context: str,
) -> int:
    """Convert a value to an integer without accepting booleans."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(
            f"{context} must be an integer."
        )

    return int(
        value
    )


def write_result(
    result: BeamformingResult,
    file_path: str | Path,
    *,
    append: bool = False,
) -> None:
    """
    Serialize one beamforming result to a trusted local pickle file.

    Parameters
    ----------
    result
        Beamforming result to serialize.
    file_path
        Destination file.
    append
        If True, append another serialized result to an existing file.
        Sequential results can be recovered with ``read_results``.

    Notes
    -----
    Pickle is used only as a provisional development format. Pickle
    files must never be loaded from untrusted sources.
    """
    if not isinstance(result, BeamformingResult):
        raise TypeError(
            "result must be a BeamformingResult instance."
        )

    path = Path(file_path).expanduser().resolve()
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    mode = "ab" if append else "wb"

    with path.open(mode) as stream:
        pickle.dump(
            result,
            stream,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


def read_results(
    file_path: str | Path,
) -> tuple[BeamformingResult, ...]:
    """
    Read all sequential beamforming results from a pickle file.

    Notes
    -----
    Only trusted local pickle files must be opened.
    """
    path = Path(file_path).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"Result file not found: {path}"
        )

    results: list[BeamformingResult] = []

    with path.open("rb") as stream:
        while True:
            try:
                result = pickle.load(stream)
            except EOFError:
                break

            if not isinstance(result, BeamformingResult):
                raise TypeError(
                    "The result file contains an unexpected object: "
                    f"{type(result).__name__}."
                )

            results.append(result)

    if not results:
        raise ValueError(
            f"Result file contains no beamforming results: {path}"
        )

    return tuple(results)


def read_result(
    file_path: str | Path,
) -> BeamformingResult:
    """
    Read exactly one beamforming result from a pickle file.

    Use ``read_results`` for incrementally appended result files.
    """
    results = read_results(
        file_path
    )

    if len(results) != 1:
        raise ValueError(
            "Result file contains more than one serialized result; "
            "use read_results() instead."
        )

    return results[0]


def read_array_config(
    config_file: str | Path,
) -> dict[str, Any]:
    """
    Read and decode an array JSON configuration.

    Parameters
    ----------
    config_file
        Path to the JSON configuration.

    Returns
    -------
    dict
        Decoded JSON object.
    """
    config_path = _resolve_config_path(config_file)

    try:
        with config_path.open(
            "r",
            encoding="utf-8",
        ) as stream:
            config = json.load(stream)

    except json.JSONDecodeError as error:
        raise ValueError(
            "Invalid JSON configuration at "
            f"line {error.lineno}, column {error.colno}: "
            f"{error.msg}"
        ) from error

    if not isinstance(config, dict):
        raise TypeError(
            "The array configuration root must be a JSON object."
        )

    return config


def _resolve_config_path(
    config_file: str | Path,
) -> Path:
    """Return an absolute existing configuration path."""
    if not isinstance(
        config_file,
        (str, Path),
    ):
        raise TypeError(
            "config_file must be a string or Path."
        )

    path = Path(config_file).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"Array configuration not found: {path}"
        )

    return path


def _validate_version(
    config: Mapping[str, Any],
) -> None:
    """Validate the configuration format version."""
    version = config.get("version", 1)

    if isinstance(version, bool) or not isinstance(version, int):
        raise TypeError(
            "Configuration field 'version' must be an integer."
        )

    if version != 1:
        raise ValueError(
            f"Unsupported array configuration version: {version}"
        )


def _resolve_data_root(
    config: Mapping[str, Any],
    config_path: Path,
) -> Path:
    """
    Resolve the waveform root directory.

    A missing ``data_root`` means the directory containing the JSON
    configuration.
    """
    value = config.get("data_root", ".")

    if not isinstance(value, str):
        raise TypeError(
            "Configuration field 'data_root' must be a string."
        )

    root = Path(value).expanduser()

    if not root.is_absolute():
        root = config_path.parent / root

    return root.resolve()


def _parse_reader_config(
    value: Any,
) -> dict[str, Any]:
    """
    Parse global or component-specific reader settings.

    Supported fields are:

    - ``format``
    - ``byte_order``
    - ``options``
    """
    if value is None:
        return {
            "format": None,
            "byte_order": None,
            "options": {},
        }

    if not isinstance(value, Mapping):
        raise TypeError(
            "Reader configuration must be a JSON object."
        )

    format_name = value.get("format")
    byte_order = value.get("byte_order")
    options = value.get("options", {})

    if format_name is not None and not isinstance(
        format_name,
        str,
    ):
        raise TypeError(
            "Reader field 'format' must be a string or null."
        )

    if byte_order is not None and not isinstance(
        byte_order,
        str,
    ):
        raise TypeError(
            "Reader field 'byte_order' must be a string or null."
        )

    if not isinstance(options, Mapping):
        raise TypeError(
            "Reader field 'options' must be a JSON object."
        )

    return {
        "format": format_name,
        "byte_order": byte_order,
        "options": dict(options),
    }


def _load_station(
    entry: Any,
    *,
    data_root: Path,
    default_reader: Mapping[str, Any],
    station_index: int,
) -> StationData:
    """Load one configured station."""
    context = f"stations[{station_index}]"

    if not isinstance(entry, Mapping):
        raise TypeError(
            f"{context} must be a JSON object."
        )

    station_id = entry.get("id")

    if not isinstance(station_id, str):
        raise TypeError(
            f"{context}.id must be a string."
        )

    station_id = station_id.strip()

    if not station_id:
        raise ValueError(
            f"{context}.id cannot be empty."
        )

    coordinates = _parse_coordinates(
        entry.get("coordinates"),
        context=f"{context}.coordinates",
    )

    covariance = _parse_covariance(
        entry.get("coordinate_covariance"),
        context=f"{context}.coordinate_covariance",
    )

    components = _parse_components(
        entry.get("components"),
        context=f"{context}.components",
    )

    component_entries = sorted(
        components,
        key=_component_sort_key,
    )

    records: list[Record] = []
    orientations: list[FloatArray] = []

    for component in component_entries:
        path = _resolve_waveform_path(
            component["file"],
            data_root,
        )

        reader_config = _merge_reader_config(
            default_reader,
            component.get("reader"),
        )

        record = _read_single_record(
            path,
            reader_config=reader_config,
            context=(
                f"station {station_id!r}, "
                f"component {component['name']!r}"
            ),
        )

        records.append(record)
        orientations.append(component["orientation"])

    return StationData(
        sid=station_id,
        records=tuple(records),
        coordinates=coordinates,
        orientation=np.stack(
            orientations,
            axis=0,
        ),
        coordinate_covariance=covariance,
    )


def _parse_coordinates(
    value: Any,
    *,
    context: str,
) -> tuple[float, ...]:
    """
    Parse station coordinates.

    Supported forms are:

    ``[x, y]``

    ``[x, y, z]``

    ``{"x": x, "y": y, "z": z}``
    """
    if isinstance(value, Mapping):
        if "x" not in value or "y" not in value:
            raise ValueError(
                f"{context} must define 'x' and 'y'."
            )

        coordinates = [
            value["x"],
            value["y"],
        ]

        if "z" in value:
            coordinates.append(value["z"])

    elif isinstance(value, Sequence) and not isinstance(
        value,
        (str, bytes),
    ):
        coordinates = list(value)

    else:
        raise TypeError(
            f"{context} must be an array or JSON object."
        )

    if len(coordinates) not in (2, 3):
        raise ValueError(
            f"{context} must contain two or three values."
        )

    parsed = tuple(
        _as_finite_float(
            item,
            context=f"{context}[{index}]",
        )
        for index, item in enumerate(coordinates)
    )

    return parsed


def _parse_covariance(
    value: Any,
    *,
    context: str,
) -> ArrayLike | None:
    """Parse an optional coordinate covariance matrix."""
    if value is None:
        return None

    try:
        covariance = np.asarray(
            value,
            dtype=float,
        )
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{context} must be a numerical matrix."
        ) from error

    if covariance.shape not in (
        (2, 2),
        (3, 3),
    ):
        raise ValueError(
            f"{context} must have shape (2, 2) or (3, 3)."
        )

    if not np.all(np.isfinite(covariance)):
        raise ValueError(
            f"{context} must contain only finite values."
        )

    return covariance


def _parse_components(
    value: Any,
    *,
    context: str,
) -> list[dict[str, Any]]:
    """
    Parse station component definitions.

    The compact form is:

    ``"components": {"E": "east.sac", "N": "north.sac"}``

    The extended form is:

    ``"components": {
        "H1": {
            "file": "axis_1.sac",
            "orientation": [0.7, 0.7, 0.0]
        }
    }``
    """
    if not isinstance(value, Mapping):
        raise TypeError(
            f"{context} must be a JSON object."
        )

    if not value:
        raise ValueError(
            f"{context} cannot be empty."
        )

    if len(value) > 3:
        raise ValueError(
            f"{context} cannot contain more than three components."
        )

    parsed: list[dict[str, Any]] = []

    for name, specification in value.items():
        if not isinstance(name, str):
            raise TypeError(
                f"{context} component names must be strings."
            )

        normalized_name = name.strip().upper()

        if not normalized_name:
            raise ValueError(
                f"{context} contains an empty component name."
            )

        parsed.append(
            _parse_component(
                normalized_name,
                specification,
                context=f"{context}.{name}",
            )
        )

    return parsed


def _parse_component(
    name: str,
    specification: Any,
    *,
    context: str,
) -> dict[str, Any]:
    """Parse one waveform-component definition."""
    if isinstance(specification, str):
        file_name = specification
        orientation = _canonical_orientation(
            name,
            context=context,
        )
        reader = None

    elif isinstance(specification, Mapping):
        file_name = specification.get("file")

        if not isinstance(file_name, str):
            raise TypeError(
                f"{context}.file must be a string."
            )

        orientation_value = specification.get("orientation")

        if orientation_value is None:
            orientation = _canonical_orientation(
                name,
                context=context,
            )
        else:
            orientation = _parse_orientation(
                orientation_value,
                context=f"{context}.orientation",
            )

        reader = specification.get("reader")

        if reader is not None:
            reader = _parse_reader_config(reader)

    else:
        raise TypeError(
            f"{context} must be a file path string or JSON object."
        )

    file_name = file_name.strip()

    if not file_name:
        raise ValueError(
            f"{context}.file cannot be empty."
        )

    return {
        "name": name,
        "file": file_name,
        "orientation": orientation,
        "reader": reader,
    }


def _canonical_orientation(
    name: str,
    *,
    context: str,
) -> FloatArray:
    """Return the default orientation associated with a component."""
    try:
        return _CANONICAL_ORIENTATIONS[name].copy()
    except KeyError as error:
        raise ValueError(
            f"{context} uses non-canonical component name {name!r}; "
            "an explicit orientation is required."
        ) from error


def _parse_orientation(
    value: Any,
    *,
    context: str,
) -> FloatArray:
    """Parse one component orientation vector."""
    try:
        orientation = np.asarray(
            value,
            dtype=float,
        )
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{context} must be a numerical array."
        ) from error

    if orientation.shape != (3,):
        raise ValueError(
            f"{context} must contain exactly three values."
        )

    if not np.all(np.isfinite(orientation)):
        raise ValueError(
            f"{context} must contain only finite values."
        )

    if np.linalg.norm(orientation) == 0.0:
        raise ValueError(
            f"{context} must define a non-zero direction."
        )

    return orientation


def _component_sort_key(
    component: Mapping[str, Any],
) -> tuple[int, str]:
    """
    Return deterministic component ordering.

    Canonical components are ordered East, North, Vertical. Components
    with explicit non-canonical names follow alphabetically.
    """
    name = component["name"]

    return (
        _COMPONENT_ORDER.get(name, 3),
        name,
    )


def _resolve_waveform_path(
    value: str,
    data_root: Path,
) -> Path:
    """Resolve one waveform file path."""
    path = Path(value).expanduser()

    if not path.is_absolute():
        path = data_root / path

    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"Waveform file not found: {path}"
        )

    return path


def _merge_reader_config(
    default: Mapping[str, Any],
    override: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Merge global and component reader settings."""
    merged = {
        "format": default.get("format"),
        "byte_order": default.get("byte_order"),
        "options": dict(
            default.get("options", {})
        ),
    }

    if override is None:
        return merged

    if override.get("format") is not None:
        merged["format"] = override["format"]

    if override.get("byte_order") is not None:
        merged["byte_order"] = override["byte_order"]

    merged["options"].update(
        override.get("options", {})
    )

    return merged


def _read_single_record(
    path: Path,
    *,
    reader_config: Mapping[str, Any],
    context: str,
) -> Record:
    """
    Read a waveform file and require exactly one record.

    The standard ShakeLab reader is used, preserving format inference
    and format-specific options.
    """
    collection = waveform_io.reader(
        file_path=str(path),
        format=reader_config.get("format"),
        byte_order=reader_config.get("byte_order"),
        **reader_config.get("options", {}),
    )

    records = tuple(
        record
        for stream in collection
        for record in stream
    )

    if not records:
        raise ValueError(
            f"No waveform record was read for {context}: {path}"
        )

    if len(records) != 1:
        raise ValueError(
            f"Expected one waveform record for {context}, "
            f"but {len(records)} records were read from {path}."
        )

    return records[0]


def _as_finite_float(
    value: Any,
    *,
    context: str,
) -> float:
    """Convert a JSON numerical value to a finite float."""
    if isinstance(value, bool) or not isinstance(
        value,
        (int, float),
    ):
        raise TypeError(
            f"{context} must be a number."
        )

    result = float(value)

    if not np.isfinite(result):
        raise ValueError(
            f"{context} must be finite."
        )

    return result