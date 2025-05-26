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

import argparse
import logging
from shakelab.cli import fdsnws
from shakelab.cli import convert

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def shakelab():
    """
    """
    parser = argparse.ArgumentParser(
        description='Shakelab: a python tool for engineering seismology.'
    )

    subparsers = parser.add_subparsers(
        dest='command', 
        help='Choose a functionality'
    )

    # Add subparser for fdsnws with detailed help
    fdsnws_parser = subparsers.add_parser(
        'fdsnws', 
        help='Download waveforms via FDSN.'
    )
    fdsnws.setup_parser(fdsnws_parser)

    # Add subparser for convert with detailed help
    convert_parser = subparsers.add_parser(
        'convert', 
        help='Convert a file to others formats.'
    )
    convert.setup_parser(convert_parser)

    args = parser.parse_args()

    if args.command == 'fdsnws':
        fdsnws.run(args)
    elif args.command == 'convert':
        convert.run(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    shakelab()