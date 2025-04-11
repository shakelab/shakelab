#!/usr/bin/env python3

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
Example use:
fdsnclient.py INGV IV.MODE..HNN 2012-05-20T02:03:23.0Z 2012-05-20T02:06:23.0Z --filename waveform.mseed
"""
import argparse
import logging
from shakelab.signals.fdsnws import FDSNClient
from shakelab.libutils.time import Date

# Configure logging (optional)
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def download_waveform(data_center, fdsn_code, starttime, endtime, 
                      correct=False, output=None):
    """
    Download waveforms from the FDSN service using FDSNClient.

    Args:
        data_center (str): The FDSN data center to download waveforms from.
        fdsn_code (str): FDSN code identifying the station
            (e.g. "IU.ANMO..BHZ").
        starttime (str): Start time in the format "YYYY-MM-DDTHH:MM:SS".
        endtime (str): End time in the format "YYYY-MM-DDTHH:MM:SS".
        correct (bool): If True, apply instrument response correction.
        output (str): File path to save the waveform (optional).
    """
    client = FDSNClient(data_center=data_center)

    # Conversion to date object allows more flexibility
    # of the input time format
    starttime = Date(starttime)
    endtime = Date(endtime)

    try:
        logger.info(f"Requesting waveform for {fdsn_code} from {starttime} "
                    f"to {endtime} using data center {data_center}")

        resp_value = client.get_waveform(fdsn_code=fdsn_code, 
                                         starttime=starttime, 
                                         endtime=endtime, 
                                         correct=correct,
                                         output=output)
        return True

    except Exception as e:
        logger.error(f"Error while downloading waveform: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Download waveforms from an FDSN service"
    )

    # Changed order: data_center is now the first required argument
    parser.add_argument('data_center', type=str, 
                        help="FDSN data center to download from")
    parser.add_argument('fdsn_code', type=str, 
                        help='FDSN code for the request (e.g. "IU.ANMO..BHZ")')
    parser.add_argument('starttime', type=str, 
                        help="Start time in format YYYY-MM-DDTHH:MM:SS")
    parser.add_argument('endtime', type=str, 
                        help="End time in format YYYY-MM-DDTHH:MM:SS")

    # Optional parameters
    parser.add_argument('-c', '--correct', action='store_true', 
                        help='Correct the instrument response')
    parser.add_argument('-f', '--filename', type=str, 
                        help="Output file name (default: output.mseed)")

    args = parser.parse_args()

    # Call the function to download the waveforms
    download_waveform(data_center=args.data_center, 
                      fdsn_code=args.fdsn_code, 
                      starttime=args.starttime, 
                      endtime=args.endtime, 
                      correct=args.correct, 
                      output=args.filename)

if __name__ == "__main__":
    main()
