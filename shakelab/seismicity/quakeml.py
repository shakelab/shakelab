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
Module to parse catalogues in quakeml format.
"""
import io, re
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime

from shakelab.seismicity.catalogue import EqDatabase, Event

# NOTE: CAMBIARE DATETIME

def parse_quakeml(xml):
    """
    Parse a QuakeML string or file path into a populated EqDatabase object.

    This parser extracts all origin and magnitude solutions, marks the
    preferred ones using preferredOriginID and preferredMagnitudeID,
    and populates each Event accordingly. XML namespaces are removed
    at the start for cleaner tag access.

    Parameters
    ----------
    xml : str or Path
        Either a QuakeML XML string or a path to a QuakeML file.

    Returns
    -------
    EqDatabase
        A database object containing Event instances with associated
        location and magnitude solutions.
    """
    if isinstance(xml, (str, Path)):
        if not xml.lstrip().startswith("<?xml"):
            path = Path(xml)
            if path.is_file():
                with open(path, "r", encoding="utf-8") as f:
                    xml = f.read()
            else:
                raise ValueError(f"File not found: {path}")
    else:
        raise TypeError("Input must be a QuakeML string or a file path.")

    xml = xml_strip(xml)
    it = ET.iterparse(io.StringIO(xml))

    # Strip namespaces for cleaner tag access
    for _, el in it:
        el.tag = re.sub(r"\{.*\}", "", el.tag)  # remove namespace
    root = it.root

    eqdb = EqDatabase()

    for event_elem in root.findall(".//event"):
        try:
            raw_id = event_elem.attrib.get("publicID", "unknown")
            eid = extract_event_id(raw_id)
            event = Event(eid)

            # Retrieve preferred IDs
            pref_oid = event_elem.findtext("preferredOriginID")
            pref_mid = event_elem.findtext("preferredMagnitudeID")

            # Parse all origins
            for origin in event_elem.findall("origin"):
                oid = origin.attrib.get("publicID", "")

                time_val = origin.findtext("time/value")
                time_unc = origin.findtext("time/uncertainty")
                lat_val = origin.findtext("latitude/value")
                lat_unc = origin.findtext("latitude/uncertainty")
                lon_val = origin.findtext("longitude/value")
                lon_unc = origin.findtext("longitude/uncertainty")
                dep_val = origin.findtext("depth/value")
                dep_unc = origin.findtext("depth/uncertainty")
                loc_agency = origin.findtext("creationInfo/agencyID")

                if not all([time_val, lat_val, lon_val, dep_val]):
                    continue

                dt = datetime.fromisoformat(time_val.replace("Z", "+00:00"))

                loc = {
                    'Year': dt.year,
                    'Month': dt.month,
                    'Day': dt.day,
                    'Hour': dt.hour,
                    'Minute': dt.minute,
                    'Second': dt.second + dt.microsecond / 1e6,
                    'Latitude': float(lat_val),
                    'Longitude': float(lon_val),
                    'Depth': float(dep_val) / 1000.0
                }

                if time_unc:
                    loc['SecError'] = float(time_unc)
                if lat_unc:
                    loc['LatError'] = float(lat_unc)
                if lon_unc:
                    loc['LonError'] = float(lon_unc)
                if dep_unc:
                    loc['DepError'] = float(dep_unc)
                if loc_agency:
                    loc['LocCode'] = loc_agency

                loc['LocPrime'] = (oid == pref_oid)
                event.location.add(loc, prime=loc['LocPrime'])

            # Parse all magnitudes
            for mag in event_elem.findall("magnitude"):
                mid = mag.attrib.get("publicID", "")

                mag_val = mag.findtext("mag/value")
                mag_type = mag.findtext("type")
                mag_err = mag.findtext("mag/uncertainty")
                mag_agency = mag.findtext("creationInfo/agencyID")

                if not mag_val:
                    continue

                m = {'MagSize': float(mag_val)}

                if mag_type:
                    m['MagType'] = mag_type
                if mag_err:
                    m['MagError'] = float(mag_err)
                if mag_agency:
                    m['MagCode'] = mag_agency

                m['MagPrime'] = (mid == pref_mid)
                event.magnitude.add(m, prime=m['MagPrime'])

            if len(event.location) > 0:
                eqdb.add(event)

        except Exception as err:
            print(f"[WARNING] Skipping event {eid}: {err}")
            continue

    return eqdb

def extract_event_id(public_id):
    """
    Extracts a clean event ID from a QuakeML publicID string.
    Returns the part after 'eventId=' if present,
    otherwise the last path component.
    """
    if "eventId=" in public_id:
        return public_id.split("eventId=")[-1]
    elif "/" in public_id:
        return public_id.rstrip("/").split("/")[-1]
    else:
        return public_id

def xml_strip(xml):
    """
    Remove white spaces and new lines from an XML string.

    Args:
        xml (str): XML content as a string.

    Returns:
        str: XML content with white spaces and new lines removed.
    """
    lines = xml.split('\n')
    buffer = []
    for x in lines:
        y = x.strip()
        if y: buffer.append(y)
    buffer = ''.join(buffer)
    return buffer