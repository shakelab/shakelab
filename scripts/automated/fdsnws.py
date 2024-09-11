#!/usr/bin/env python3
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

def setup_parser(parser):
    parser.add_argument('network', type=str, help='Network name (e.g., INGV)')
    parser.add_argument('station', type=str, help='Station and channel (e.g., IV.MODE..HNN)')
    parser.add_argument('starttime', type=str, help='Start time in ISO8601 format (e.g., 2012-05-20T02:03:23.000000Z)')
    parser.add_argument('endtime', type=str, help='End time in ISO8601 format (e.g., 2012-05-20T02:06:23.000000Z)')
    parser.add_argument('--filename', type=str, required=True, help='Filename to save the downloaded waveform (e.g., waveform.mseed)')

def run(args):
    # Execute logic for fdsnws with provided arguments
    print(f"Downloading waveform for {args.network} {args.station} from {args.starttime} to {args.endtime}")
    # Actual implementation should interface with FDSN API to fetch data
    with open(args.filename, 'w') as f:
        f.write(f'Fake waveform for {args.station} saved in {args.filename}.')
    print(f"Waveform saved in {args.filename}")

