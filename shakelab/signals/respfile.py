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
RESP file parser for ShakeLab ResponseCollection.

This module provides functionality to parse RESP files or RESP strings
into a ShakeLab ResponseCollection, compatible with the internal
StageSet and StageResponse classes.
"""
import os
from shakelab.signals import response as rspm
from shakelab.libutils.time import Date


def parse_resp(resp_input):
    """
    Parse a RESP file or RESP string into a ResponseCollection.

    Args:
        resp_input (str): Path to RESP file.

    Returns:
        ResponseCollection: ShakeLab ResponseCollection object containing
                             parsed stages.
    """
    if os.path.isfile(resp_input):
        with open(resp_input, 'r') as f:
            resp_text = f.read()
    else:
        raise FileNotFoundError(f"RESP file not found: {resp_input}")

    lines = resp_text.splitlines()

    station, network, location, channel = None, None, "", None
    starttime, endtime = None, None
    stages = []
    current_stage = None
    rcoll = rspm.ResponseCollection()

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Header information
        if line.startswith('B050F03'):
            station = line.split(':', 1)[-1].strip()
        elif line.startswith('B050F16'):
            network = line.split(':', 1)[-1].strip()
        elif line.startswith('B052F03'):
            location = line.split(':', 1)[-1].strip()
        elif line.startswith('B052F04'):
            channel = line.split(':', 1)[-1].strip()
        elif line.startswith('B052F22'):
            starttime = parse_resp_date(line.split(':', 1)[-1].strip())
        elif line.startswith('B052F23'):
            endtime = parse_resp_date(line.split(':', 1)[-1].strip())

        # Start of new stage
        elif line.startswith(('B053F04', 'B054F04', 'B058F03')):
            if current_stage:
                stages.append(current_stage)
            current_stage = {}
            current_stage['stage_number'] = int(
                line.split(':', 1)[-1].strip()
            )

        # PAZ stage
        elif line.startswith('B053F07'):
            current_stage['normalization_factor'] = float(
                line.split(':', 1)[-1].strip()
            )
        elif line.startswith('B053F08'):
            current_stage['normalization_frequency'] = float(
                line.split(':', 1)[-1].strip()
            )
        elif line.startswith('B053F09'):
            current_stage['zeros'] = []
        elif line.startswith('B053F14'):
            current_stage['poles'] = []
        elif line.startswith('B053F15-18'):
            parts = line.split()
            real = float(parts[2])
            imag = float(parts[3])
            if 'poles' in current_stage:
                current_stage['poles'].append(complex(real, imag))
            elif 'zeros' in current_stage:
                current_stage['zeros'].append(complex(real, imag))

        # FIR stage
        elif line.startswith('B054F07'):
            current_stage['numerator'] = []
        elif line.startswith('B054F10'):
            current_stage['denominator'] = [1.0]
        elif line.startswith('B054F08-09'):
            parts = line.split()
            coeff = float(parts[2])
            if 'numerator' in current_stage:
                current_stage['numerator'].append(coeff)

        # Gain stage
        elif line.startswith('B058F04'):
            current_stage['sensitivity'] = float(
                line.split(':', 1)[-1].strip()
            )
            current_stage['frequency'] = 1.0

    if current_stage:
        stages.append(current_stage)

    # Assemble ResponseCollection
    if station and network and channel:
        sid = f"{network}.{station}.{location}.{channel}"
        stream_response = rspm.StreamResponse(sid)
        stageset = rspm.StageSet(starttime, endtime)

        for stage in stages:
            if 'poles' in stage or 'zeros' in stage:
                data = {
                    k: v for k, v in stage.items()
                    if k in rspm.StagePoleZero._KEYMAP
                }
                s = rspm.StagePoleZero(data)
            elif 'numerator' in stage:
                data = {
                    k: v for k, v in stage.items()
                    if k in rspm.StageFIR._KEYMAP
                }
                s = rspm.StageFIR(data)
            elif 'sensitivity' in stage:
                data = {
                    k: v for k, v in stage.items()
                    if k in rspm.StageGain._KEYMAP
                }
                s = rspm.StageGain(data)
            else:
                continue

            stageset.append(s)

        stream_response.append(stageset)
        rcoll.append(stream_response)

    return rcoll


def parse_resp_date(date_str):
    """
    Parse a RESP date string into a ShakeLab Date object.

    Args:
        date_str (str): Date string in the format 'YYYY,DOY,HH:MM:SS.sss'.

    Returns:
        Date: ShakeLab Date object representing the parsed date.
    """
    try:
        parts = date_str.replace(',', ' ').replace(':', ' ')\
                        .replace('.', ' ').split()
        year = int(parts[0])
        doy = int(parts[1])
        hour = int(parts[2])
        minute = int(parts[3])
        second = float(parts[4])

        return Date([year, doy, hour, minute, second])

    except Exception:
        return Date('now')
