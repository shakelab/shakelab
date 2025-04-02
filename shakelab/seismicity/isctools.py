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
"""
import os
import requests
from urllib.parse import urlencode

from shakelab.libutils.time import Date

class ISCDownloader:
    """
    A downloader for retrieving ISC (International Seismological Centre)
    bulletins in ISF format.
    """
    BASE_URL = "https://www.isc.ac.uk/cgi-bin/web-db-run"

    params = {
        "CatalogueType" : "REVIEWED",
        "OutputFormat" : "ISF",
        "SearchAreaShape" : "RECT",
        "RectangleBottomLatitude" : "",
        "RectangleTopLatitude" : "",
        "RectangleLeftLongitude" : "",
        "RectangleRightLongitude" : "",
        "CircularLatitude" : "",
        "CircularLongitude" : "",
        "CircularRadius" : "",
        "MaxDistanceUnits" : "deg",
        "SeismicRegionNumber" : "",
        "GeogrephicRegionNumber" : "",
        "PolygonCoordinates" : "",
        "StartYear" : "",
        "StartMonth" : "",
        "StartDay" : "",
        "StartTime" : "00:00:00",
        "EndYear" : "",
        "EndMonth" : "",
        "EndDay" : "",
        "EndTime" : "23:59:59",
        "MinimumDepth" : "",
        "MaximumDepth" : "",
        "NoDepthEvents" : "on",
        "MinimumMagnitude" : "",
        "MaximumMagnitude" : "",
        "NoMagnitudeEvents" : "",
        "MagnitudeType" : "",
        "MagnitudeAgency" : "",
        "FocalMechanismAgency" : "Any",
        "IncludePhases" : "off",
        "MinimumPhaseNumber" : "",
        "MaximumPhaseNumber" : "",
        "NoKnownPhases" : "",
        "PrimeOnly" : "",
        "IncludeMagnitudes" : "on",
        "IncludeHeaders" : "on",
        "IncludeComments" : "off",
        "IncludeLinks" : "off"
        }

    def __init__(self, verbose=True):
        """
        Initializes the downloader with default parameters set to None.
        """
        self.content = ""
        self.verbose = verbose

    def help(self):
        """
        Print usage information and available parameters for ISCDownloader.
    
        This includes:
        - Description of mandatory and optional parameters
        - How to set spatial, temporal, and magnitude filters
        - How to download and save data
        - New features: error handling, verbose output, content validation
        """
        print("ISC Bulletin Downloader - Help")
        print("-----------------------------------")
        print("Mandatory parameters (must be set via specific methods):")
        print("  - Time window (start and end)")
        print("  - Spatial bounds (latitude and longitude)")
        print()
        print("Optional parameters (can be set via set_field or directly):")
        print("  - Magnitude bounds")
        print("  - Depth range")
        print("  - Magnitude type or agency")
        print("  - Include phases, comments, headers, etc.")
        print()
        print("New and improved features:")
        print("  - Robust error handling for HTTP requests and file saving")
        print("  - Automatic validation of downloaded ISF content")
        print("  - Informative messages during operations (if verbose=True)")
        print()
        print("Example usage:")
        print("  downloader = ISCDownloader()")
        print("  downloader.set_time_bounds(\"2022-01-01\", \"2022-01-31\")")
        print("  downloader.set_spatial_bounds((10.0, 15.0), (40.0, 45.0))")
        print("  downloader.set_magnitude_bounds(4.0, 6.0)")
        print("  downloader.download()")
        print("  downloader.save(\"output.isf\")")

    def print(self):

        print("\nCURRENT SETTINGS:\n")

        for key, value in self.params.items():
            if value == "": value = None
            print("\t{0} : {1}".format(key, value))

    def set_field(self, field_name, field_value):
        """
        """
        if field_name in self.params.keys():
            self.params[field_name] = str(field_value)
        else:
            print('Error: wrong field name')

    def set_time_bounds(self, start_time, end_time):
        """
        Set the start and end times for the ISC query.
        Format: "YYYY-MM-DDTHH:MM:SS.SSSS"
        """
        if not isinstance(start_time, Date):
            start_time = Date(start_time)

        if not isinstance(end_time, Date):
            end_time = Date(end_time)

        self.params["StartYear"] = str(start_time.year)
        self.params["StartMonth"] = str(start_time.month)
        self.params["StartDay"] = str(start_time.day)

        self.params["StartTime"] = f"{start_time.hour:02}:"
        self.params["StartTime"] += f"{start_time.minute:02}:"
        self.params["StartTime"] += f"{int(start_time.second):02}"

        self.params["EndYear"] = str(end_time.year)
        self.params["EndMonth"] = str(end_time.month)
        self.params["EndDay"] = str(end_time.day)

        self.params["EndTime"] = f"{end_time.hour:02}:"
        self.params["EndTime"] += f"{end_time.minute:02}:"
        self.params["EndTime"] += f"{int(end_time.second):02}"

    def set_spatial_bounds(self, longitude_bounds, latitude_bounds):
        """
        Set spatial bounds using tuples: (min, max) for latitude and longitude.
        """
        self.params["SearchAreaShape"] = "RECT"
        self.params["RectangleLeftLongitude"] = longitude_bounds[0]
        self.params["RectangleRightLongitude"] = longitude_bounds[1]
        self.params["RectangleBottomLatitude"] = latitude_bounds[0]
        self.params["RectangleTopLatitude"] = latitude_bounds[1]

    def set_magnitude_bounds(self, min_mag, max_mag):
        """
        Set minimum and maximum magnitude for the ISC query.
        """
        self.params["MinimumMagnitude"] = min_mag
        self.params["MinimumMagnitude"] = max_mag

    def download(self, timeout=60):
        """
        Download ISC bulletin data using the configured parameters.
        Parameters:
        timeout (int): Timeout for the HTTP request in seconds.
        """
        query_param_map = {
            "CatalogueType" : "request",
            "OutputFormat" : "out_format",
            "SearchAreaShape" : "searchshape",
            "RectangleBottomLatitude" : "bot_lat",
            "RectangleTopLatitude" : "top_lat",
            "RectangleLeftLongitude" : "left_lon",
            "RectangleRightLongitude" : "right_lon",
            "CircularLatitude" : "ctr_lat",
            "CircularLongitude" : "ctr_lon",
            "CircularRadius" : "radius",
            "MaxDistanceUnits" : "max_dist_units",
            "SeismicRegionNumber" : "srn",
            "GeogrephicRegionNumber" : "grn",
            "PolygonCoordinates" : "coordvals",
            "StartYear" : "start_year",
            "StartMonth" : "start_month",
            "StartDay" : "start_day",
            "StartTime" : "start_time",
            "EndYear" : "end_year",
            "EndMonth" : "end_month",
            "EndDay" : "end_day",
            "EndTime" : "end_time",
            "MinimumDepth" : "min_dep",
            "MaximumDepth" : "max_dep",
            "NoDepthEvents" : "null_dep",
            "MinimumMagnitude" : "min_mag",
            "MaximumMagnitude" : "max_mag",
            "NoMagnitudeEvents" : "null_mag",
            "MagnitudeType" : "req_mag_type",
            "MagnitudeAgency" : "req_mag_agcy",
            "FocalMechanismAgency" : "req_fm_agcy",
            "IncludePhases" : "include_phases",
            "MinimumPhaseNumber" : "min_def",
            "MaximumPhaseNumber" : "max_def",
            "NoKnownPhases" : "null_phs",
            "PrimeOnly" : "prime_only",
            "IncludeMagnitudes" : "include_magnitudes",
            "IncludeHeaders" : "include_headers",
            "IncludeComments" : "include_comments",
            "IncludeLinks" : "include_links"
            }

        query_params = {}
        for key, value in self.params.items():
            query_params[query_param_map[key]] = str(value)
    
        query_url = "{}?{}".format(self.BASE_URL, urlencode(query_params))
    
        if getattr(self, "verbose", True):
            print("Sending request to ISC server...")
    
        try:
            response = requests.get(query_url, timeout=timeout)
            response.raise_for_status()
            self.content = self._extract_isf_data(response.text)
    
            if getattr(self, "verbose", True):
                print("Download completed. Response size: {} characters"
                      .format(len(self.content)))
    
        except requests.exceptions.Timeout:
            print("Error: request timed out.")
        except requests.exceptions.HTTPError as e:
            print("HTTP error: {}".format(e))
        except requests.exceptions.RequestException as e:
            print("Connection error: {}".format(e))
        
        self.content = self._extract_isf_data(response.text)

    def save(self, output_file, overwrite=False):
        """
        Save the downloaded ISC bulletin to a file.
    
        Parameters:
            output_file (str): Path to the output file.
            overwrite (bool): If True, overwrite the file if it exists.
        """
        if not self.content:
            print("Warning: no data to save. Did you run download()?") 
            return
    
        if os.path.isfile(output_file) and not overwrite:
            print("Error: the file '{}' already exists.".format(output_file))
            print("Use overwrite=True to replace it.")
            return
    
        try:
            with open(output_file, "w", encoding="utf-8") as file:
                file.write(self.content)
            print("Bulletin successfully saved to '{}'".format(output_file))
        except OSError as e:
            print("Error saving file '{}': {}".format(output_file, e))

    def to_database(self):
        """
        """
        pass

    def _extract_isf_data(self, content):
        """
        Cleans the ISF content by removing everything before 'DATA_TYPE'
        and after 'STOP', excluding 'STOP' itself.
        """
        if not content or not isinstance(content, str):
            print("Warning: empty or invalid content received.")
            return ""
    
        lines = content.splitlines()
    
        try:
            start_idx = next(
                (i for i, line in enumerate(lines) if "DATA_TYPE" in line), None
            )
            end_idx = next(
                (i for i, line in enumerate(lines) if line.strip() == "STOP"), None
            )
        except Exception as e:
            print("Error while scanning lines: {}".format(e))
            return ""
    
        if start_idx is None or end_idx is None or start_idx >= end_idx:
            print("Warning: could not locate valid ISF section.")
            return ""
    
        cleaned_lines = lines[start_idx:end_idx]
        return "\n".join(cleaned_lines) + "\n"


# Example usage
if __name__ == "__main__":

    downloader = ISCDownloader()

    # Print detailed help
    # downloader.help()

    # downloader.print()

    # Set the required parameters and download
    downloader.set_time_bounds(
        "2022-02-28T00:00:00",
        "2023-03-01T23:59:59"
        )

    # Set the search bounds (longitude-latitude)
    downloader.set_spatial_bounds(
        (11.0, 14.0),
        (45.0, 47.0)
        )

    downloader.set_magnitude_bounds(4, 10)

    # Download the bulletin
    downloader.download()

    # Save data to file
    downloader.save("isc_reviewed_bulletin.isf", overwrite=True)

