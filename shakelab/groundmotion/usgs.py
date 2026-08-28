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
ShakeMap readers and lightweight container classes.

This module provides generic access to standard USGS ShakeMap products,
independently of the ShakeScenario ground-motion provider layer.

Main components
---------------
ShakeMapGMData
    Container for ShakeMap ground-motion grid values and optional
    uncertainty values.

ShakeMapEvent
    Container for one ShakeMap event directory, including event metadata,
    ground-motion grid data, uncertainty data, and optional station data.

Design principles
-----------------
- This module performs file parsing, validation, interpolation, and export.
- It does not depend on GroundMotionProvider or ScenarioEvent.
- ShakeMap product fields are preserved using their native names
  (e.g. PGA, PGV, PSA03, STDPGA).
- Unit conversion is intentionally not performed here.
- Missing or inconsistent ShakeMap products fail explicitly rather than being
  silently filled with zero values.

Notes
-----
ShakeMap XML products are read from ``grid.xml`` and, when available,
``uncertainty.xml``. Both files are expected to use the same grid geometry.
"""

from __future__ import annotations

import csv
import json
import math
import xml.etree.ElementTree as ET

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata


_SHAKEMAP_NAMESPACE = "http://earthquake.usgs.gov/eqcenter/shakemap"
_DEFAULT_GRID_FILE = "grid.xml"
_DEFAULT_UNCERTAINTY_FILE = "uncertainty.xml"
_DEFAULT_INFO_FILE = "info.json"
_DEFAULT_STATIONS_FILE = "stationlist.json"


# ---------------------------------------------------------------------------
# Small immutable metadata structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ShakeMapGridField:
    """
    Metadata describing one ShakeMap grid field.

    Parameters
    ----------
    index
        One-based column index as stored by ShakeMap.
    name
        Field name, for example ``LON``, ``LAT``, ``PGA`` or ``STDPGA``.
    units
        Field units declared in the XML metadata.
    """

    index: int
    name: str
    units: str = ""


@dataclass(frozen=True)
class ShakeMapGridSpec:
    """
    ShakeMap grid geometry metadata.

    Parameters
    ----------
    lon_min, lon_max
        Longitude bounds in degrees.
    lat_min, lat_max
        Latitude bounds in degrees.
    lon_spacing, lat_spacing
        Grid spacing in degrees.
    nlon, nlat
        Number of grid nodes along longitude and latitude.
    """

    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float
    lon_spacing: float
    lat_spacing: float
    nlon: int
    nlat: int


@dataclass(frozen=True)
class ShakeMapGrid:
    """
    Parsed ShakeMap grid product.

    Parameters
    ----------
    rows
        Grid records, one dictionary per ShakeMap grid point.
    fields
        Ordered grid-field metadata.
    spec
        Optional grid geometry metadata from ``grid_specification``.
    """

    rows: tuple[Mapping[str, float], ...]
    fields: tuple[ShakeMapGridField, ...]
    spec: Optional[ShakeMapGridSpec] = None

    @property
    def field_names(self) -> tuple[str, ...]:
        """Return ordered field names."""
        return tuple(field.name for field in self.fields)

    @property
    def units(self) -> Dict[str, str]:
        """Return field-unit mapping."""
        return {
            field.name: field.units
            for field in self.fields
        }


# ---------------------------------------------------------------------------
# Ground-motion data container
# ---------------------------------------------------------------------------


class ShakeMapGMData:
    """
    Container for ShakeMap ground-motion values and optional uncertainties.

    ``grid`` and ``uncertainty`` preserve ShakeMap native field names.
    Coordinates are expected to be stored as ``LON`` and ``LAT``.

    When uncertainty data are supplied, the coordinate geometry is validated
    against the main grid.
    """

    def __init__(
        self,
        grid_data: Union[
            ShakeMapGrid,
            Sequence[Mapping[str, float]],
        ],
        uncertainty_data: Optional[
            Union[
                ShakeMapGrid,
                Sequence[Mapping[str, float]],
            ]
        ] = None,
    ) -> None:
        self._grid_product = self._normalize_product(grid_data)
        self._uncertainty_product = (
            self._normalize_product(uncertainty_data)
            if uncertainty_data is not None
            else None
        )

        self.grid: List[Dict[str, float]] = [
            dict(row) for row in self._grid_product.rows
        ]
        self.uncertainty: Optional[List[Dict[str, float]]] = None

        if self._uncertainty_product is not None:
            self.uncertainty = [
                dict(row)
                for row in self._uncertainty_product.rows
            ]
            self._validate_uncertainty_alignment()

        self._regular_axes: Optional[
            Tuple[np.ndarray, np.ndarray]
        ] = None
        self._regular_values: Dict[
            Tuple[str, str],
            np.ndarray
        ] = {}
        self._interpolators: Dict[
            Tuple[str, str, str],
            RegularGridInterpolator
        ] = {}

        self._prepare_regular_grid()

    @staticmethod
    def _normalize_product(
        product: Union[
            ShakeMapGrid,
            Sequence[Mapping[str, float]],
        ]
    ) -> ShakeMapGrid:
        """Normalize legacy row lists to a ShakeMapGrid object."""
        if isinstance(product, ShakeMapGrid):
            if not product.rows:
                raise ValueError("ShakeMap grid is empty.")
            return product

        if not isinstance(product, Sequence) or isinstance(
            product,
            (str, bytes),
        ):
            raise TypeError(
                "grid data must be a ShakeMapGrid or a sequence of mappings."
            )

        rows: List[Mapping[str, float]] = []
        for row in product:
            if not isinstance(row, Mapping):
                raise TypeError(
                    "ShakeMap grid rows must be mappings."
                )
            rows.append(dict(row))

        if not rows:
            raise ValueError("ShakeMap grid is empty.")

        fields = tuple(
            ShakeMapGridField(
                index=i + 1,
                name=str(name),
                units="",
            )
            for i, name in enumerate(rows[0].keys())
        )

        return ShakeMapGrid(
            rows=tuple(rows),
            fields=fields,
            spec=None,
        )

    @property
    def grid_product(self) -> ShakeMapGrid:
        """Return parsed main grid product."""
        return self._grid_product

    @property
    def uncertainty_product(self) -> Optional[ShakeMapGrid]:
        """Return parsed uncertainty product, if available."""
        return self._uncertainty_product

    @property
    def field_names(self) -> tuple[str, ...]:
        """Return fields available in the main grid."""
        return self._grid_product.field_names

    @property
    def uncertainty_field_names(self) -> tuple[str, ...]:
        """Return fields available in uncertainty data."""
        if self._uncertainty_product is None:
            return ()
        return self._uncertainty_product.field_names

    @property
    def units(self) -> Dict[str, str]:
        """
        Return combined field-unit mapping.

        Main-grid units are returned first and uncertainty-field units are
        added when present.
        """
        units = dict(self._grid_product.units)
        if self._uncertainty_product is not None:
            units.update(self._uncertainty_product.units)
        return units

    def __len__(self) -> int:
        return len(self.grid)

    def has_param(self, name: str) -> bool:
        """Return whether a parameter exists in grid or uncertainty data."""
        key = str(name).strip().upper()
        if key in self.field_names:
            return True
        if key in self.uncertainty_field_names:
            return True
        return False

    def get_param(self, name: str) -> List[float]:
        """
        Retrieve all values for a ShakeMap parameter.

        Parameters
        ----------
        name
            Native ShakeMap field name.

        Returns
        -------
        list of float
            Parameter values in native ShakeMap units.

        Raises
        ------
        KeyError
            If the field does not exist.
        """
        key = str(name).strip().upper()

        if (
            self._uncertainty_product is not None
            and key in self.uncertainty_field_names
        ):
            return [
                float(row[key])
                for row in self.uncertainty
            ]

        if key in self.field_names:
            return [
                float(row[key])
                for row in self.grid
            ]

        raise KeyError(
            f"ShakeMap parameter not found: {key!r}."
        )

    def get_point(self, index: int) -> Dict[str, float]:
        """
        Return all main-grid and uncertainty values for one point index.
        """
        record = self.grid[index].copy()
        if self.uncertainty is not None:
            for key, value in self.uncertainty[index].items():
                if key in ("LON", "LAT"):
                    continue
                record[key] = value
        return record

    def get_ground_motion(
        self,
        sites: Union[
            Tuple[float, float],
            Sequence[Tuple[float, float]],
        ],
        parameters: Optional[Sequence[str]] = None,
        method: str = "linear",
        outside: str = "nan",
    ) -> List[Dict[str, float]]:
        """
        Interpolate ShakeMap values at one or more site coordinates.

        Parameters
        ----------
        sites
            One ``(longitude, latitude)`` pair or a sequence of pairs.
        parameters
            Native ShakeMap parameters to interpolate. If omitted, all
            non-coordinate fields from grid and uncertainty products are used.
        method
            Interpolation method. Supported values are ``linear`` and
            ``nearest``.
        outside
            Behavior for points outside the ShakeMap grid:
            - ``"nan"``: return NaN;
            - ``"raise"``: raise ValueError.

        Returns
        -------
        list of dict
            One dictionary per requested site.
        """
        site_list = _normalize_site_tuples(sites)

        parameter_list = self._normalize_parameters(parameters)

        method_key = str(method).strip().lower()
        if method_key not in {"linear", "nearest"}:
            raise ValueError(
                "Interpolation method must be 'linear' or 'nearest'."
            )

        outside_key = str(outside).strip().lower()
        if outside_key not in {"nan", "raise"}:
            raise ValueError(
                "outside must be either 'nan' or 'raise'."
            )

        if self._regular_axes is not None:
            output = self._interpolate_regular(
                site_list,
                parameter_list,
                method=method_key,
            )
        else:
            output = self._interpolate_scattered(
                site_list,
                parameter_list,
                method=method_key,
            )

        if outside_key == "raise":
            for site, values in zip(site_list, output):
                missing = [
                    key
                    for key, value in values.items()
                    if not math.isfinite(value)
                ]
                if missing:
                    raise ValueError(
                        "Site is outside ShakeMap interpolation domain or "
                        "contains unavailable data: "
                        f"lon={site[0]}, lat={site[1]}, "
                        f"parameters={missing}."
                    )

        return output

    # ------------------------------------------------------------------
    # Grid validation / interpolation preparation
    # ------------------------------------------------------------------

    def _validate_uncertainty_alignment(self) -> None:
        """Validate grid/uncertainty size and point coordinates."""
        assert self.uncertainty is not None

        if len(self.grid) != len(self.uncertainty):
            raise ValueError(
                "Mismatch in grid and uncertainty lengths: "
                f"{len(self.grid)} != {len(self.uncertainty)}."
            )

        for index, (grid_row, unc_row) in enumerate(
            zip(self.grid, self.uncertainty)
        ):
            try:
                grid_lon = float(grid_row["LON"])
                grid_lat = float(grid_row["LAT"])
                unc_lon = float(unc_row["LON"])
                unc_lat = float(unc_row["LAT"])
            except KeyError as exc:
                raise ValueError(
                    "Both grid and uncertainty data must contain "
                    "LON and LAT fields."
                ) from exc

            if not (
                math.isclose(
                    grid_lon,
                    unc_lon,
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
                and math.isclose(
                    grid_lat,
                    unc_lat,
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    "Grid and uncertainty coordinates are not aligned "
                    f"at index {index}: "
                    f"({grid_lon}, {grid_lat}) != "
                    f"({unc_lon}, {unc_lat})."
                )

    def _prepare_regular_grid(self) -> None:
        """
        Detect a complete rectilinear grid and prepare cached arrays.

        If the data do not form a complete Cartesian longitude/latitude grid,
        interpolation falls back to ``scipy.interpolate.griddata``.
        """
        lons = np.asarray(
            [float(row["LON"]) for row in self.grid],
            dtype=float,
        )
        lats = np.asarray(
            [float(row["LAT"]) for row in self.grid],
            dtype=float,
        )

        unique_lons = np.unique(lons)
        unique_lats = np.unique(lats)

        if unique_lons.size * unique_lats.size != len(self.grid):
            return

        coord_to_index: Dict[Tuple[float, float], int] = {}

        for i, (lon, lat) in enumerate(zip(lons, lats)):
            key = (float(lon), float(lat))
            if key in coord_to_index:
                return
            coord_to_index[key] = i

        if len(coord_to_index) != len(self.grid):
            return

        for lat in unique_lats:
            for lon in unique_lons:
                if (float(lon), float(lat)) not in coord_to_index:
                    return

        self._regular_axes = (
            unique_lats.astype(float),
            unique_lons.astype(float),
        )

        for source_name, rows in self._parameter_sources():
            available = [
                key
                for key in rows[0].keys()
                if key not in ("LON", "LAT")
            ]

            for param in available:
                values = np.empty(
                    (unique_lats.size, unique_lons.size),
                    dtype=float,
                )

                for iy, lat in enumerate(unique_lats):
                    for ix, lon in enumerate(unique_lons):
                        row_index = coord_to_index[
                            (float(lon), float(lat))
                        ]
                        values[iy, ix] = float(
                            rows[row_index].get(
                                param,
                                float("nan"),
                            )
                        )

                self._regular_values[
                    (source_name, str(param).upper())
                ] = values

    def _parameter_sources(
        self,
    ) -> Iterable[Tuple[str, List[Dict[str, float]]]]:
        yield "grid", self.grid
        if self.uncertainty is not None:
            yield "uncertainty", self.uncertainty

    def _normalize_parameters(
        self,
        parameters: Optional[Sequence[str]],
    ) -> List[str]:
        if parameters is None:
            out = [
                key
                for key in self.field_names
                if key not in ("LON", "LAT")
            ]
            out.extend(
                key
                for key in self.uncertainty_field_names
                if key not in ("LON", "LAT")
                and key not in out
            )
            return list(out)

        out: List[str] = []
        for param in parameters:
            key = str(param).strip().upper()
            if not key:
                raise ValueError(
                    "ShakeMap parameter names must be non-empty."
                )
            if not self.has_param(key):
                raise KeyError(
                    f"ShakeMap parameter not found: {key!r}."
                )
            if key not in out:
                out.append(key)

        return out

    def _source_for_parameter(
        self,
        parameter: str,
    ) -> Tuple[str, List[Dict[str, float]]]:
        key = str(parameter).upper()

        if (
            self.uncertainty is not None
            and key in self.uncertainty_field_names
        ):
            return "uncertainty", self.uncertainty

        if key in self.field_names:
            return "grid", self.grid

        raise KeyError(
            f"ShakeMap parameter not found: {key!r}."
        )

    def _interpolate_regular(
        self,
        sites: List[Tuple[float, float]],
        parameters: Sequence[str],
        *,
        method: str,
    ) -> List[Dict[str, float]]:
        assert self._regular_axes is not None

        lat_axis, lon_axis = self._regular_axes

        xi = np.asarray(
            [(lat, lon) for lon, lat in sites],
            dtype=float,
        )

        output = [{} for _ in sites]

        for param in parameters:
            source_name, _rows = self._source_for_parameter(param)
            value_key = (source_name, param)
            values = self._regular_values[value_key]

            interp_key = (source_name, param, method)
            interpolator = self._interpolators.get(interp_key)

            if interpolator is None:
                interpolator = RegularGridInterpolator(
                    (lat_axis, lon_axis),
                    values,
                    method=method,
                    bounds_error=False,
                    fill_value=np.nan,
                )
                self._interpolators[interp_key] = interpolator

            interpolated = np.asarray(
                interpolator(xi),
                dtype=float,
            )

            for i, value in enumerate(interpolated):
                output[i][param] = float(value)

        return output

    def _interpolate_scattered(
        self,
        sites: List[Tuple[float, float]],
        parameters: Sequence[str],
        *,
        method: str,
    ) -> List[Dict[str, float]]:
        points = np.asarray(
            [
                (float(row["LON"]), float(row["LAT"]))
                for row in self.grid
            ],
            dtype=float,
        )
        xi = np.asarray(sites, dtype=float)

        output = [{} for _ in sites]

        for param in parameters:
            _source_name, rows = self._source_for_parameter(param)
            values = np.asarray(
                [
                    float(row.get(param, float("nan")))
                    for row in rows
                ],
                dtype=float,
            )

            interpolated = griddata(
                points=points,
                values=values,
                xi=xi,
                method=method,
                fill_value=np.nan,
            )

            for i, value in enumerate(interpolated):
                output[i][param] = float(value)

        return output

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------

    def export_csv(
        self,
        path: Union[str, Path],
        site_file: Optional[Union[str, Path]] = None,
        parameters: Optional[Sequence[str]] = None,
        format: str = "geojson",
        outside: str = "nan",
    ) -> None:
        """Export interpolated ShakeMap values to CSV."""
        path = Path(path)
        sites = self._get_sites(site_file, format)

        values = self.get_ground_motion(
            [
                (site["longitude"], site["latitude"])
                for site in sites
            ],
            parameters=parameters,
            outside=outside,
        )

        if not sites:
            raise ValueError("No sites available for export.")

        fieldnames = ["id", "longitude", "latitude"]
        fieldnames.extend(values[0].keys())

        with path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
            )
            writer.writeheader()

            for site, value in zip(sites, values):
                writer.writerow(
                    {
                        "id": site["id"],
                        "longitude": site["longitude"],
                        "latitude": site["latitude"],
                        **value,
                    }
                )

    def export_geojson(
        self,
        path: Union[str, Path],
        site_file: Optional[Union[str, Path]] = None,
        parameters: Optional[Sequence[str]] = None,
        format: str = "geojson",
        outside: str = "nan",
    ) -> None:
        """Export interpolated ShakeMap values to GeoJSON."""
        path = Path(path)
        sites = self._get_sites(site_file, format)

        values = self.get_ground_motion(
            [
                (site["longitude"], site["latitude"])
                for site in sites
            ],
            parameters=parameters,
            outside=outside,
        )

        features = []

        for site, value in zip(sites, values):
            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": [
                            site["longitude"],
                            site["latitude"],
                        ],
                    },
                    "properties": {
                        "id": site["id"],
                        **value,
                    },
                }
            )

        geojson = {
            "type": "FeatureCollection",
            "features": features,
        }

        path.write_text(
            json.dumps(
                geojson,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

    def _get_sites(
        self,
        site_file: Optional[Union[str, Path]],
        fmt: str,
    ) -> List[Dict[str, Any]]:
        """
        Load interpolation sites from CSV or GeoJSON.

        If no site file is supplied, the internal ShakeMap grid itself is
        returned as the site list.
        """
        if site_file is None:
            return [
                {
                    "id": str(i),
                    "longitude": float(row["LON"]),
                    "latitude": float(row["LAT"]),
                }
                for i, row in enumerate(self.grid)
            ]

        site_path = Path(site_file)
        fmt_key = str(fmt).strip().lower()

        id_keys = ("id", "site_id", "name")
        lon_keys = ("longitude", "lon")
        lat_keys = ("latitude", "lat")

        if fmt_key == "csv":
            with site_path.open(
                "r",
                encoding="utf-8",
            ) as f:
                reader = csv.DictReader(f)
                sites: List[Dict[str, Any]] = []

                for row in reader:
                    sites.append(
                        {
                            "id": _first_mapping_value(
                                row,
                                id_keys,
                                "id",
                            ),
                            "longitude": float(
                                _first_mapping_value(
                                    row,
                                    lon_keys,
                                    "longitude",
                                )
                            ),
                            "latitude": float(
                                _first_mapping_value(
                                    row,
                                    lat_keys,
                                    "latitude",
                                )
                            ),
                        }
                    )

                return sites

        if fmt_key == "geojson":
            data = _load_json(site_path)

            features = data.get("features")
            if not isinstance(features, list):
                raise ValueError(
                    "Invalid GeoJSON: missing 'features' array."
                )

            sites = []

            for feature in features:
                try:
                    properties = feature["properties"]
                    geometry = feature["geometry"]
                    coordinates = geometry["coordinates"]
                except (KeyError, TypeError) as exc:
                    raise ValueError(
                        "Invalid GeoJSON feature."
                    ) from exc

                if (
                    not isinstance(coordinates, Sequence)
                    or len(coordinates) < 2
                ):
                    raise ValueError(
                        "GeoJSON Point coordinates must contain "
                        "[longitude, latitude]."
                    )

                sites.append(
                    {
                        "id": _first_mapping_value(
                            properties,
                            id_keys,
                            "id",
                        ),
                        "longitude": float(coordinates[0]),
                        "latitude": float(coordinates[1]),
                    }
                )

            return sites

        raise ValueError(
            f"Unsupported site file format: {fmt!r}."
        )


# ---------------------------------------------------------------------------
# ShakeMap event container
# ---------------------------------------------------------------------------


class ShakeMapEvent:
    """
    Container for one standard ShakeMap event directory.

    The event can be built directly from in-memory objects or loaded from a
    directory containing standard ShakeMap products.
    """

    def __init__(
        self,
        info: Optional[Dict[str, Any]] = None,
        gm_data: Optional[ShakeMapGMData] = None,
        stations: Optional[Dict[str, Any]] = None,
        folder: Optional[Union[str, Path]] = None,
        load_stations: bool = True,
        load_uncertainty: bool = True,
        require_info: bool = False,
    ) -> None:
        self.info: Optional[Dict[str, Any]] = info
        self.gm_data: Optional[ShakeMapGMData] = gm_data
        self.stations: Optional[Dict[str, Any]] = stations
        self.folder: Optional[Path] = None

        if folder is not None:
            self.load_folder(
                folder,
                load_stations=load_stations,
                load_uncertainty=load_uncertainty,
                require_info=require_info,
            )

    @property
    def event_id(self) -> str:
        """Return ShakeMap event identifier when available."""
        return str(
            self._event_information().get(
                "event_id",
                "N/A",
            )
        )

    @property
    def magnitude(self) -> Any:
        """Return event magnitude metadata when available."""
        return self._event_information().get(
            "magnitude",
            "N/A",
        )

    @property
    def origin_time(self) -> Any:
        """Return event origin time metadata when available."""
        return self._event_information().get(
            "origin_time",
            "N/A",
        )

    @property
    def latitude(self) -> Any:
        """Return event latitude metadata when available."""
        return self._event_information().get(
            "latitude",
            "N/A",
        )

    @property
    def longitude(self) -> Any:
        """Return event longitude metadata when available."""
        return self._event_information().get(
            "longitude",
            "N/A",
        )

    @property
    def depth(self) -> Any:
        """Return event depth metadata when available."""
        return self._event_information().get(
            "depth",
            "N/A",
        )

    @property
    def units(self) -> Dict[str, str]:
        """
        Return ShakeMap field units.

        XML metadata is preferred because it describes the actual loaded
        products. ``info.json`` units are used only as an additional source.
        """
        units: Dict[str, str] = {}

        if self.gm_data is not None:
            units.update(self.gm_data.units)

        if self.info:
            ground_motions = (
                self.info.get("output", {})
                .get("ground_motions", {})
            )

            if isinstance(ground_motions, Mapping):
                for key, value in ground_motions.items():
                    if (
                        isinstance(value, Mapping)
                        and "units" in value
                    ):
                        units.setdefault(
                            str(key).upper(),
                            str(value["units"]),
                        )

        return units

    def _event_information(self) -> Mapping[str, Any]:
        if not self.info:
            return {}

        event_information = (
            self.info.get("input", {})
            .get("event_information", {})
        )

        if isinstance(event_information, Mapping):
            return event_information

        return {}

    def load_folder(
        self,
        folder: Union[str, Path],
        load_stations: bool = True,
        load_uncertainty: bool = True,
        require_info: bool = False,
    ) -> None:
        """
        Load a standard ShakeMap event directory.

        ``grid.xml`` is required. Other products are optional unless
        explicitly required.
        """
        event_dir = Path(folder).expanduser().resolve()

        if not event_dir.exists():
            raise FileNotFoundError(
                f"ShakeMap directory not found: {event_dir}"
            )

        if not event_dir.is_dir():
            raise ValueError(
                f"ShakeMap path is not a directory: {event_dir}"
            )

        self.folder = event_dir

        grid_path = event_dir / _DEFAULT_GRID_FILE
        if not grid_path.is_file():
            raise FileNotFoundError(
                f"ShakeMap grid file not found: {grid_path}"
            )

        info_path = event_dir / _DEFAULT_INFO_FILE

        if info_path.is_file():
            self.load_info(info_path)
        elif require_info:
            raise FileNotFoundError(
                f"ShakeMap info file not found: {info_path}"
            )
        else:
            self.info = None

        uncertainty_path = (
            event_dir / _DEFAULT_UNCERTAINTY_FILE
        )

        grid_product = _load_grid_product(grid_path)

        uncertainty_product = None

        if load_uncertainty and uncertainty_path.is_file():
            uncertainty_product = _load_grid_product(
                uncertainty_path
            )

        self.gm_data = ShakeMapGMData(
            grid_data=grid_product,
            uncertainty_data=uncertainty_product,
        )

        stations_path = event_dir / _DEFAULT_STATIONS_FILE

        if load_stations and stations_path.is_file():
            self.load_stations(stations_path)
        else:
            self.stations = None

    def load_info(
        self,
        path: Union[str, Path],
    ) -> None:
        """Load event metadata from ``info.json``."""
        self.info = _load_json(Path(path))

    def load_grid(
        self,
        path: Union[str, Path],
    ) -> None:
        """
        Load a main ShakeMap grid.

        Existing uncertainty data are discarded because their geometry may no
        longer correspond to the newly loaded grid.
        """
        product = _load_grid_product(Path(path))
        self.gm_data = ShakeMapGMData(
            grid_data=product,
        )

    def load_uncertainty(
        self,
        path: Union[str, Path],
    ) -> None:
        """
        Load uncertainty data and validate geometry against the main grid.
        """
        if self.gm_data is None:
            raise RuntimeError(
                "Grid must be loaded before uncertainty."
            )

        uncertainty = _load_grid_product(Path(path))

        self.gm_data = ShakeMapGMData(
            grid_data=self.gm_data.grid_product,
            uncertainty_data=uncertainty,
        )

    def load_stations(
        self,
        path: Union[str, Path],
    ) -> None:
        """Load ``stationlist.json``."""
        self.stations = _load_json(Path(path))

    def get_ground_motion(
        self,
        sites: Union[
            Tuple[float, float],
            Sequence[Tuple[float, float]],
        ],
        parameters: Optional[Sequence[str]] = None,
        method: str = "linear",
        outside: str = "nan",
    ) -> List[Dict[str, float]]:
        """Interpolate loaded ShakeMap fields at requested sites."""
        if self.gm_data is None:
            raise RuntimeError("Grid data not loaded.")

        return self.gm_data.get_ground_motion(
            sites=sites,
            parameters=parameters,
            method=method,
            outside=outside,
        )

    def export(
        self,
        path: Union[str, Path],
        site_file: Optional[Union[str, Path]] = None,
        parameters: Optional[Sequence[str]] = None,
        format: str = "geojson",
        outside: str = "nan",
    ) -> None:
        """Export ShakeMap values to CSV or GeoJSON."""
        if self.gm_data is None:
            raise RuntimeError("Grid data not loaded.")

        output_path = Path(path)
        suffix = output_path.suffix.lower()

        if suffix == ".csv":
            self.gm_data.export_csv(
                output_path,
                site_file=site_file,
                parameters=parameters,
                format=format,
                outside=outside,
            )
            return

        if suffix in (".geojson", ".json"):
            self.gm_data.export_geojson(
                output_path,
                site_file=site_file,
                parameters=parameters,
                format=format,
                outside=outside,
            )
            return

        raise ValueError(
            f"Unsupported output file extension: {suffix!r}."
        )

    def summary(self) -> None:
        """Print a concise summary of loaded ShakeMap products."""
        print("ShakeMap Summary")
        print("----------------")
        print(f"Event ID:     {self.event_id}")
        print(f"Magnitude:    {self.magnitude}")
        print(f"Origin time:  {self.origin_time}")
        print(
            "Epicenter:    "
            f"lat={self.latitude}, lon={self.longitude}"
        )
        print(f"Depth:        {self.depth}")

        if self.gm_data is None:
            print("Grid data:    Not loaded")
        else:
            print(f"Grid points:  {len(self.gm_data)}")
            print(
                "GM parameters: "
                + ", ".join(
                    key
                    for key in self.gm_data.field_names
                    if key not in ("LON", "LAT")
                )
            )

            if self.gm_data.uncertainty is None:
                print("Uncertainty:  Not loaded")
            else:
                print(
                    "STD parameters: "
                    + ", ".join(
                        key
                        for key in self.gm_data.uncertainty_field_names
                        if key not in ("LON", "LAT")
                    )
                )

            lons = self.gm_data.get_param("LON")
            lats = self.gm_data.get_param("LAT")

            print(
                "Longitude range: "
                f"{min(lons):.4f} - {max(lons):.4f}"
            )
            print(
                "Latitude range:  "
                f"{min(lats):.4f} - {max(lats):.4f}"
            )

        if self.stations:
            features = self.stations.get("features", [])
            count = len(features) if isinstance(features, list) else 0
            print(f"Stations:     {count}")
        else:
            print("Stations:     Not loaded")

        if self.units:
            print("Units:")
            for key, unit in sorted(self.units.items()):
                print(f"  {key}: {unit}")


# ---------------------------------------------------------------------------
# XML / JSON parsing helpers
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> Dict[str, Any]:
    """Load one JSON object from disk."""
    try:
        data = json.loads(
            path.read_text(encoding="utf-8")
        )
    except OSError as exc:
        raise ValueError(
            f"Cannot read JSON file: {path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON file: {path}"
        ) from exc

    if not isinstance(data, dict):
        raise ValueError(
            f"JSON root must be an object: {path}"
        )

    return data


def _load_grid(path: Path) -> List[Dict[str, float]]:
    """
    Backward-compatible helper returning only parsed ShakeMap rows.

    New code should prefer ``_load_grid_product`` because it preserves field
    metadata and grid geometry.
    """
    product = _load_grid_product(path)
    return [
        dict(row)
        for row in product.rows
    ]


def _load_grid_product(path: Path) -> ShakeMapGrid:
    """
    Parse a ShakeMap ``grid.xml`` or ``uncertainty.xml`` product.

    Returns
    -------
    ShakeMapGrid
        Parsed rows plus field and grid-geometry metadata.
    """
    try:
        tree = ET.parse(path)
    except (OSError, ET.ParseError) as exc:
        raise ValueError(
            f"Cannot parse ShakeMap XML file: {path}"
        ) from exc

    root = tree.getroot()

    fields = _parse_grid_fields(root)

    if not fields:
        raise ValueError(
            f"No grid_field definitions found in {path}."
        )

    field_names = [
        field.name
        for field in fields
    ]

    if len(set(field_names)) != len(field_names):
        raise ValueError(
            f"Duplicate grid_field names in {path}."
        )

    grid_data_element = _find_child(
        root,
        "grid_data",
    )

    if (
        grid_data_element is None
        or grid_data_element.text is None
        or not grid_data_element.text.strip()
    ):
        raise ValueError(
            f"No grid_data found in {path}."
        )

    rows: List[Mapping[str, float]] = []

    for line_number, line in enumerate(
        grid_data_element.text.strip().splitlines(),
        start=1,
    ):
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()

        if len(parts) != len(fields):
            raise ValueError(
                "ShakeMap grid_data column count mismatch "
                f"in {path}, line {line_number}: "
                f"{len(parts)} values for {len(fields)} fields."
            )

        try:
            values = [
                float(value)
                for value in parts
            ]
        except ValueError as exc:
            raise ValueError(
                "Non-numeric ShakeMap grid_data value "
                f"in {path}, line {line_number}."
            ) from exc

        rows.append(
            dict(
                zip(
                    field_names,
                    values,
                )
            )
        )

    if not rows:
        raise ValueError(
            f"ShakeMap grid contains no rows: {path}"
        )

    for required in ("LON", "LAT"):
        if required not in field_names:
            raise ValueError(
                f"ShakeMap grid missing required field {required!r}: "
                f"{path}"
            )

    spec = _parse_grid_specification(root)

    if spec is not None:
        expected_count = spec.nlon * spec.nlat
        if expected_count != len(rows):
            raise ValueError(
                "ShakeMap grid size does not match "
                "grid_specification: "
                f"expected {expected_count}, got {len(rows)} "
                f"in {path}."
            )

    return ShakeMapGrid(
        rows=tuple(rows),
        fields=tuple(fields),
        spec=spec,
    )


def _parse_grid_fields(
    root: ET.Element,
) -> List[ShakeMapGridField]:
    field_elements = _find_all_children(
        root,
        "grid_field",
    )

    fields: List[ShakeMapGridField] = []

    for order, element in enumerate(
        field_elements,
        start=1,
    ):
        name = element.attrib.get("name")
        if not name:
            raise ValueError(
                "ShakeMap grid_field missing 'name'."
            )

        raw_index = element.attrib.get(
            "index",
            str(order),
        )

        try:
            index = int(raw_index)
        except ValueError as exc:
            raise ValueError(
                "ShakeMap grid_field has invalid index: "
                f"{raw_index!r}."
            ) from exc

        fields.append(
            ShakeMapGridField(
                index=index,
                name=str(name).strip().upper(),
                units=str(
                    element.attrib.get("units", "")
                ).strip(),
            )
        )

    fields.sort(key=lambda field: field.index)

    return fields


def _parse_grid_specification(
    root: ET.Element,
) -> Optional[ShakeMapGridSpec]:
    element = _find_child(
        root,
        "grid_specification",
    )

    if element is None:
        return None

    required = {
        "lon_min": "lon_min",
        "lon_max": "lon_max",
        "lat_min": "lat_min",
        "lat_max": "lat_max",
        "lon_spacing": "nominal_lon_spacing",
        "lat_spacing": "nominal_lat_spacing",
        "nlon": "nlon",
        "nlat": "nlat",
    }

    values: Dict[str, Any] = {}

    for target, source in required.items():
        raw = element.attrib.get(source)
        if raw is None:
            return None

        try:
            if target in ("nlon", "nlat"):
                values[target] = int(raw)
            else:
                values[target] = float(raw)
        except ValueError:
            return None

    return ShakeMapGridSpec(**values)


def _find_child(
    root: ET.Element,
    local_name: str,
) -> Optional[ET.Element]:
    for element in root.iter():
        if _local_tag_name(element.tag) == local_name:
            return element
    return None


def _find_all_children(
    root: ET.Element,
    local_name: str,
) -> List[ET.Element]:
    return [
        element
        for element in root.iter()
        if _local_tag_name(element.tag) == local_name
    ]


def _local_tag_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", 1)[1]
    return tag


def _normalize_site_tuples(
    sites: Union[
        Tuple[float, float],
        Sequence[Tuple[float, float]],
    ],
) -> List[Tuple[float, float]]:
    if (
        isinstance(sites, tuple)
        and len(sites) == 2
        and all(
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            for value in sites
        )
    ):
        site_values = [sites]
    else:
        if not isinstance(sites, Sequence) or isinstance(
            sites,
            (str, bytes),
        ):
            raise TypeError(
                "sites must be a (lon, lat) tuple or a sequence "
                "of (lon, lat) tuples."
            )
        site_values = list(sites)

    out: List[Tuple[float, float]] = []

    for site in site_values:
        if (
            not isinstance(site, Sequence)
            or len(site) != 2
        ):
            raise ValueError(
                "Each site must contain longitude and latitude."
            )

        lon = float(site[0])
        lat = float(site[1])

        if not (
            math.isfinite(lon)
            and math.isfinite(lat)
        ):
            raise ValueError(
                "Site longitude and latitude must be finite."
            )

        out.append((lon, lat))

    return out


def _first_mapping_value(
    mapping: Mapping[str, Any],
    keys: Sequence[str],
    label: str,
) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]

    raise ValueError(
        f"Missing required field {label!r}; "
        f"expected one of {list(keys)}."
    )


__all__ = [
    "ShakeMapGridField",
    "ShakeMapGridSpec",
    "ShakeMapGrid",
    "ShakeMapGMData",
    "ShakeMapEvent",
]
