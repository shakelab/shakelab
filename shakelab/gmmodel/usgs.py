# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
ShakeMap reader and container classes for USGS ShakeMap products.

Includes:
- ShakeMapGMData: ground motion values and uncertainties
- ShakeMapEvent: metadata + shakemap data wrapper
"""
import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from scipy.interpolate import griddata


class ShakeMapGMData:
    """
    Container for ShakeMap ground motion values and uncertainties.

    Includes both grid.xml and uncertainty.xml contents, matched by index.
    """

    def __init__(self,
                 grid_data: List[Dict[str, float]],
                 uncertainty_data: Optional[List[Dict[str, float]]] = None):
        """
        Initialize the ShakeMapGMData object.

        Parameters
        ----------
        grid_data : list of dict
            Parsed records from 'grid.xml'.
        uncertainty_data : list of dict, optional
            Parsed records from 'uncertainty.xml', if available.
        """
        self.grid = grid_data
        self.uncertainty = uncertainty_data

        if uncertainty_data:
            if len(grid_data) != len(uncertainty_data):
                raise ValueError("Mismatch in grid and uncertainty lengths.")

    def get_param(self, name: str) -> List[float]:
        """
        Retrieve a list of values for a given parameter.

        Parameters
        ----------
        name : str
            Name of the parameter (e.g., 'PGA', 'STDPGA', 'LAT', 'LON').

        Returns
        -------
        list of float
            Values for the parameter across all grid points.
        """
        if name.startswith("STD") and self.uncertainty:
            return [row[name] for row in self.uncertainty if name in row]
        return [row[name] for row in self.grid if name in row]

    def get_point(self, index: int) -> Dict[str, float]:
        """
        Return a dictionary with all parameters at a given point index.

        Parameters
        ----------
        index : int
            Index of the point.

        Returns
        -------
        dict
            Combined parameters (including uncertainty if available).
        """
        record = self.grid[index].copy()
        if self.uncertainty:
            record.update(self.uncertainty[index])
        return record

    def __len__(self):
        return len(self.grid)

    def get_ground_motion(
        self,
        sites: Union[Tuple[float, float], List[Tuple[float, float]]],
        parameters: Optional[List[str]] = None,
        nan_fill: Optional[float] = None
        ) -> List[Dict[str, float]]:
        """
        Interpolate ground motion values at one or more site coordinates.
    
        Parameters
        ----------
        sites : tuple or list of tuples
            Coordinates (longitude, latitude) of the site(s).
        parameters : list of str, optional
            Parameters to retrieve. If None, all available are returned.
        nan_fill : float, optional
            Value to replace NaNs. If None, NaNs are preserved.
    
        Returns
        -------
        list of dict
            One dictionary per site with {parameter: value} pairs.
        """
        if isinstance(sites, tuple):
            sites = [sites]
    
        x = [row["LON"] for row in self.grid]
        y = [row["LAT"] for row in self.grid]
    
        if parameters is None:
            parameters = [k for k in self.grid[0] if k not in ("LON", "LAT")]
            if self.uncertainty:
                parameters += [
                    k for k in self.uncertainty[0]
                    if k.startswith("STD") and k not in ("LON", "LAT")
                ]
    
        output = [{} for _ in sites]
    
        for param in parameters:
            if param.startswith("STD") and self.uncertainty:
                z = [row.get(param, float("nan")) for row in self.uncertainty]
            else:
                z = [row.get(param, float("nan")) for row in self.grid]
    
            interp = griddata(
                points=list(zip(x, y)),
                values=z,
                xi=sites,
                method="linear"
            )
    
            for i, val in enumerate(interp):
                if nan_fill is not None and (val is None or str(val) == "nan"):
                    output[i][param] = nan_fill
                else:
                    output[i][param] = float(val)
    
        return output

    def export_csv(self,
                   path,
                   site_file: Optional[Path] = None,
                   parameters: Optional[List[str]] = None,
                   format: str = "geojson",
                   nan_fill: Optional[float] = None) -> None:
        """
        Export ground motion data to CSV.
    
        Parameters
        ----------
        path : str or Path
            Output CSV file path.
        site_file : str or Path, optional
            CSV or GeoJSON file with new site grid (id, longitude, latitude).
        parameters : list of str, optional
            Parameters to export (e.g., ['PGA', 'MMI']). If None, export all.
        format : str, default 'geojson'
            Format of site_file: 'csv' or 'geojson'.
        nan_fill : float, optional
            Value to substitute for NaN results. If None, NaNs are preserved.
        """
        path = Path(path)
        if site_file is not None:
            site_file = Path(site_file)
    
        sites = self._get_sites(site_file, format)
        values = self._interpolate(sites, parameters, nan_fill)
    
        sample_id = next(iter(values))
        fieldnames = ["id", "longitude", "latitude"]
        fieldnames += list(values[sample_id].keys())
    
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for site in sites:
                row = {
                    "id": site["id"],
                    "longitude": site["longitude"],
                    "latitude": site["latitude"],
                    **values[site["id"]]
                }
                writer.writerow(row)    

    def export_geojson(self,
                       path,
                       site_file: Optional[Path] = None,
                       parameters: Optional[List[str]] = None,
                       format: str = "geojson",
                       nan_fill: Optional[float] = None) -> None:
        """
        Export ground motion data to GeoJSON FeatureCollection.
    
        Parameters
        ----------
        path : str or Path
            Output GeoJSON file path.
        site_file : str or Path, optional
            CSV or GeoJSON file with new site grid (id, longitude, latitude).
        parameters : list of str, optional
            Parameters to export (e.g., ['PGA', 'MMI']). If None, export all.
        format : str, default 'geojson'
            Format of site_file: 'csv' or 'geojson'.
        nan_fill : float, optional
            Value to substitute for NaN results. If None, NaNs are preserved.
        """
        path = Path(path)
        if site_file is not None:
            site_file = Path(site_file)
    
        sites = self._get_sites(site_file, format)
        values = self._interpolate(sites, parameters, nan_fill)
    
        features = []
        for site in sites:
            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [site["longitude"], site["latitude"]]
                },
                "properties": {
                    "id": site["id"],
                    **values[site["id"]]
                }
            }
            features.append(feature)
    
        geojson = {
            "type": "FeatureCollection",
            "features": features
        }
    
        with open(path, "w", encoding="utf-8") as f:
            json.dump(geojson, f, indent=2)

    def _get_sites(self,
                   site_file: Optional[Path],
                   fmt: str) -> List[Dict[str, float]]:
        """
        Internal: Load site list from CSV or GeoJSON, or use grid as fallback.
    
        Supports flexible field names like 'lon', 'latitude', 'site_id'.
    
        Parameters
        ----------
        site_file : str or Path or None
            Path to site list file. If None, use internal ShakeMap grid.
        fmt : str
            Format: 'csv' or 'geojson'.
    
        Returns
        -------
        list of dict
            Each dict contains: id, longitude, latitude.
        """
        if site_file is None:
            return [
                {
                    "id": str(i),
                    "longitude": row["LON"],
                    "latitude": row["LAT"]
                }
                for i, row in enumerate(self.grid)
            ]
    
        site_file = Path(site_file)
    
        id_keys = ["id", "site_id", "name"]
        lon_keys = ["longitude", "lon"]
        lat_keys = ["latitude", "lat"]
    
        def _get_first_key(d: Dict[str, str],
                           keys: List[str],
                           label: str) -> str:
            for k in keys:
                if k in d:
                    return d[k]
            raise ValueError(f"Missing required field '{label}' "
                             f"(expected one of {keys})")
    
        if fmt == "csv":
            with open(site_file, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                sites = []
                for row in reader:
                    try:
                        site = {
                            "id": _get_first_key(row, id_keys, "id"),
                            "longitude": float(_get_first_key(
                                row, lon_keys, "longitude")),
                            "latitude": float(_get_first_key(
                                row, lat_keys, "latitude"))
                        }
                    except ValueError as e:
                        raise ValueError(f"Error parsing site file: {e}")
                    sites.append(site)
                return sites
    
        if fmt == "geojson":
            with open(site_file, "r", encoding="utf-8") as f:
                gj = json.load(f)
    
            if "features" not in gj:
                raise ValueError("Invalid GeoJSON: missing 'features' field")
    
            sites = []
            for feat in gj["features"]:
                try:
                    props = feat["properties"]
                    coords = feat["geometry"]["coordinates"]
                    site = {
                        "id": _get_first_key(props, id_keys, "id"),
                        "longitude": float(coords[0]),
                        "latitude": float(coords[1])
                    }
                except (KeyError, TypeError, IndexError, ValueError) as e:
                    raise ValueError(
                        "GeoJSON feature must include 'properties.id' and "
                        "Point geometry with [lon, lat] coordinates."
                    )
                sites.append(site)
            return sites
    
        raise ValueError(f"Unsupported site file format: {fmt}")
        
    def _interpolate(self,
                     sites: List[Dict],
                     parameters: Optional[List[str]],
                     nan_fill: Optional[float] = None
                     ) -> Dict[str, Dict[str, float]]:
        """
        Internal: interpolate specified parameters onto given site list.
    
        Parameters
        ----------
        sites : list of dict
            Sites with id, longitude, latitude.
        parameters : list of str, optional
            Parameters to interpolate. If None, use all available.
        nan_fill : float, optional
            Value to substitute for NaN results. If None, NaNs are preserved.
    
        Returns
        -------
        dict
            Dictionary keyed by site id, each with {param: value}.
        """
        x = [row["LON"] for row in self.grid]
        y = [row["LAT"] for row in self.grid]
    
        if parameters is None:
            parameters = [k for k in self.grid[0] if k not in ("LON", "LAT")]
            if self.uncertainty:
                parameters += [
                    k for k in self.uncertainty[0] if k not in ("LON", "LAT")
                ]
    
        output = {site["id"]: {} for site in sites}
        xi = [(s["longitude"], s["latitude"]) for s in sites]
    
        for param in parameters:
            if param.startswith("STD") and self.uncertainty:
                z = [row.get(param, float("nan")) for row in self.uncertainty]
            else:
                z = [row.get(param, float("nan")) for row in self.grid]
    
            interp = griddata(points=list(zip(x, y)), values=z, xi=xi,
                              method="linear")
    
            for val, site in zip(interp, sites):
                if nan_fill is not None and (val is None or str(val) == "nan"):
                    output[site["id"]][param] = nan_fill
                else:
                    output[site["id"]][param] = float(val)
    
        return output


class ShakeMapEvent:
    """
    Container class for a ShakeMap event, including metadata and grid data.

    Can be initialized directly from data dictionaries or by loading from
    a ShakeMap folder structure. Supports incremental loading of individual
    files after instantiation.
    """

    def __init__(self,
                 info: Optional[Dict] = None,
                 gm_data: Optional[ShakeMapGMData] = None,
                 stations: Optional[Dict] = None,
                 folder: Optional[Path] = None,
                 load_stations: bool = True,
                 load_uncertainty: bool = True):
        """
        Initialize the ShakeMapEvent object.

        Parameters
        ----------
        info : dict, optional
            Parsed contents of 'info.json'.
        gm_data : ShakeMapGMData, optional
            Ground motion data object.
        stations : dict, optional
            Parsed contents of 'stationlist.json'.
        folder : str or Path, optional
            If provided, load all files from this folder.
        load_stations : bool, default True
            Whether to load stationlist.json if available.
        load_uncertainty : bool, default True
            Whether to load uncertainty.xml if available.
        """
        self.info = None
        self.gm_data = None
        self.stations = None

        if folder is not None:
            self.load_folder(folder,
                             load_stations=load_stations,
                             load_uncertainty=load_uncertainty)
        else:
            self.info = info
            self.gm_data = gm_data
            self.stations = stations

    @property
    def event_id(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("event_id", "N/A")
    
    
    @property
    def magnitude(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("magnitude", "N/A")
    
    
    @property
    def origin_time(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("origin_time", "N/A")
    
    
    @property
    def latitude(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("latitude", "N/A")
    
    
    @property
    def longitude(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("longitude", "N/A")
    
    
    @property
    def depth(self) -> str:
        if not self.info:
            return "N/A"
        return self.info.get("input", {}).get(
            "event_information", {}
        ).get("depth", "N/A")

    @property
    def units(self) -> Dict[str, str]:
        if not self.info:
            return {}
        return {
            k: v.get("units", "unknown")
            for k, v in self.info.get("output", {}).get(
                "ground_motions", {}
                ).items()
        }

    def load_folder(self,
                    folder: Path,
                    load_stations: bool = True,
                    load_uncertainty: bool = True) -> None:
        """
        Load ShakeMap data from a standard folder.

        Parameters
        ----------
        folder : str or Path
            Path to the ShakeMap folder.
        load_stations : bool
            Whether to load stationlist.json if available.
        load_uncertainty : bool
            Whether to load uncertainty.xml if available.
        """
        folder = Path(folder)
        self.load_info(folder / "info.json")
        self.load_grid(folder / "grid.xml")
        if load_uncertainty and (folder / "uncertainty.xml").exists():
            self.load_uncertainty(folder / "uncertainty.xml")
        if load_stations and (folder / "stationlist.json").exists():
            self.load_stations(folder / "stationlist.json")

    def load_info(self, path: Path) -> None:
        """Load event metadata from an info.json file."""
        self.info = _load_json(path)

    def load_grid(self, path: Path) -> None:
        """Load ground motion data from a grid.xml file."""
        grid = _load_grid(path)
        self.gm_data = ShakeMapGMData(grid_data=grid)

    def load_uncertainty(self, path: Path) -> None:
        """Load uncertainty data from an uncertainty.xml file."""
        if not self.gm_data:
            raise RuntimeError("Grid must be loaded before uncertainty.")
        unc = _load_grid(path)
        self.gm_data.uncertainty = unc

    def load_stations(self, path: Path) -> None:
        """Load station observations from a stationlist.json file."""
        self.stations = _load_json(path)

    def get_ground_motion(
        self,
        sites: Union[Tuple[float, float], List[Tuple[float, float]]],
        parameters: Optional[List[str]] = None,
        nan_fill: Optional[float] = None
        ) -> List[Dict[str, float]]:
        """
        Get ground motion values at one or more sites via interpolation.
    
        Parameters
        ----------
        sites : tuple or list of tuples
            Coordinates (longitude, latitude) of one or more sites.
        parameters : list of str, optional
            Parameters to retrieve. If None, all available are returned.
        nan_fill : float, optional
            Value to replace NaNs. If None, NaNs are preserved.
    
        Returns
        -------
        list of dict
            One dictionary per site with {parameter: value} pairs.
        """
        if not self.gm_data:
            raise RuntimeError("Grid data not loaded.")
        return self.gm_data.get_ground_motion(
            sites=sites,
            parameters=parameters,
            nan_fill=nan_fill
        )

    def export(self,
               path,
               site_file: Optional[Path] = None,
               parameters: Optional[List[str]] = None,
               format: str = "geojson",
               nan_fill: Optional[float] = None) -> None:
        """
        Export the ground motion data to CSV or GeoJSON.

        Parameters
        ----------
        path : str or Path
            Output file path (.csv or .geojson).
        site_file : str or Path, optional
            Optional input file with site list (CSV or GeoJSON).
        parameters : list of str, optional
            Parameters to export. If None, all available will be used.
        format : str, default 'geojson'
            Format of site_file: 'csv' or 'geojson'.
        nan_fill : float, optional
            Value to assign for missing interpolated data (NaN).
        """
        path = Path(path)
        if site_file is not None:
            site_file = Path(site_file)

        suffix = path.suffix.lower()
        if suffix == ".csv":
            self.gm_data.export_csv(path,
                                    site_file=site_file,
                                    parameters=parameters,
                                    format=format,
                                    nan_fill=nan_fill)
        elif suffix in (".geojson", ".json"):
            self.gm_data.export_geojson(path,
                                        site_file=site_file,
                                        parameters=parameters,
                                        format=format,
                                        nan_fill=nan_fill)
        else:
            raise ValueError("Unsupported output file extension: " + suffix)

    def summary(self) -> None:
        """
        Print a summary of the ShakeMap event and its loaded components.
        """
        print("ShakeMap Summary")
        print("----------------")
    
        # Event metadata
        event_id = "N/A"
        magnitude = "N/A"
        origin_time = "N/A"
        lat = "N/A"
        lon = "N/A"
    
        if self.info:
            info_event = self.info.get("input", {}).get(
                "event_information", {}
            )
            event_id = info_event.get("event_id", "N/A")
            magnitude = info_event.get("magnitude", "N/A")
            origin_time = info_event.get("origin_time", "N/A")
            lat = info_event.get("latitude", "N/A")
            lon = info_event.get("longitude", "N/A")
    
        print(f"Event ID:     {event_id}")
        print(f"Magnitude:    {magnitude}")
        print(f"Origin time:  {origin_time}")
        print(f"Epicenter:    lat={lat}, lon={lon}")
    
        # Grid summary
        if self.gm_data:
            npts = len(self.gm_data)
            print(f"Grid points:  {npts}")
    
            sample_grid = self.gm_data.grid[0]
            grid_params = [
                k for k in sample_grid if k not in ("LON", "LAT")
            ]
            print("Available GM parameters:",
                  ", ".join(grid_params))
    
            if self.gm_data.uncertainty:
                sample_unc = self.gm_data.uncertainty[0]
                unc_params = [
                    k for k in sample_unc if k not in ("LON", "LAT")
                ]
                print("Available STD parameters:",
                      ", ".join(unc_params))
            else:
                print("Uncertainty:  Not loaded")
    
            lons = self.gm_data.get_param("LON")
            lats = self.gm_data.get_param("LAT")
            print(f"Longitude range:  {min(lons):.2f}° – {max(lons):.2f}°")
            print(f"Latitude range:   {min(lats):.2f}° – {max(lats):.2f}°")
        else:
            print("Grid data:    Not loaded")
    
        # Station info
        if self.stations:
            nsta = len(self.stations.get("features", []))
            print(f"Stations loaded: Yes ({nsta} entries)")
        else:
            print("Stations loaded: No")

        # Unit info
        if self.info:
            print("Units (from info.json):")
            units = self.get_units()
            for k, u in units.items():
                print(f"  {k}: {u}")

# -------------------------------------------------------------------------
# Internal helper functions
# -------------------------------------------------------------------------

def _load_json(path: Path) -> Dict:
    """Internal function to load a JSON file."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_grid(path: Path) -> List[Dict[str, float]]:
    """Internal function to load ShakeMap grid or uncertainty XML file."""
    tree = ET.parse(path)
    root = tree.getroot()
    ns = {"sm": "http://earthquake.usgs.gov/eqcenter/shakemap"}

    fields = root.findall("sm:grid_field", ns)
    columns = [field.attrib["name"] for field in fields]

    lines = root.find("sm:grid_data", ns).text.strip().split("\n")
    data = []
    for line in lines:
        values = list(map(float, line.strip().split()))
        data.append(dict(zip(columns, values)))
    return data