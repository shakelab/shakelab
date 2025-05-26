#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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

import logging
from shakelab.signals.fdsnws import FDSNClient
from shakelab.libutils.time import Date

# Configure logging (optional)
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def setup_parser(parser):
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
    parser.add_argument('-o', '--output', type=str, 
                        help="Output file name (default: output.mseed)")

def run(args):
    """
    """
    client = FDSNClient(data_center=args.data_center)

    # Conversion to date object allows more flexibility
    # of the input time format
    starttime = Date(args.starttime)
    endtime = Date(args.endtime)

    try:
        logger.info(f"Requesting waveform for {args.fdsn_code} from {args.starttime} "
                    f"to {args.endtime} using data center {args.data_center}")

        resp_value = client.get_waveform(fdsn_code=args.fdsn_code, 
                                         starttime=args.starttime, 
                                         endtime=args.endtime, 
                                         correct=args.correct,
                                         output=args.output)
        return True

    except Exception as e:
        logger.error(f"Error while downloading waveform: {e}")
        return False