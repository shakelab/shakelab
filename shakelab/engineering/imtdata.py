# -*- coding: utf-8 -*-
# ****************************************************************************
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
# ****************************************************************************
"""
ShakeScenario IMT utilities.

This module converts ShakeLab waveform collections into the canonical
``ShakeScenarioIMT`` JSON structure.

The expected waveform source identifiers follow the FDSN convention:

    NETWORK.STATION.LOCATION.CHANNEL

For example:

    OX.30001..EHN
    OX.30001..EHE
    OX.30001..EHZ

Input waveforms are assumed to represent ground velocity in m/s.

Station geometry and channel metadata are supplied through a
``ShakeLabStation`` JSON document. Event metadata can be supplied through a
``ShakeLabEvent`` JSON document.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from shakelab.engineering.engpar import ResponseSpectrum
from shakelab.signals.base import StreamCollection


SCHEMA_TYPE = "ShakeScenarioIMT"
SCHEMA_VERSION = "1.0.0"

STATION_SCHEMA_TYPE = "ShakeLabStation"
EVENT_SCHEMA_TYPE = "ShakeLabEvent"

DEFAULT_PERIODS = (
    0.1,
    0.2,
    0.3,
    0.5,
    1.0,
    2.0,
)

DEFAULT_DAMPING = 0.05
DEFAULT_COMBINATION = "geometric_mean"


def build_imt_dict(
    stream_collection,
    stations,
    event,
    *,
    metadata=None,
    periods=DEFAULT_PERIODS,
    damping=DEFAULT_DAMPING,
    combination=DEFAULT_COMBINATION,
):
    """
    Build a ``ShakeScenarioIMT`` dictionary.

    Parameters
    ----------
    stream_collection : StreamCollection
        ShakeLab waveform collection containing horizontal ground-velocity
        components.
    stations : dict
        ``ShakeLabStation`` document containing station geometry and channel
        metadata.
    event : dict
        Event dictionary containing at least ``event_id``, ``origin_time``,
        ``magnitude``, ``longitude``, ``latitude``, and ``depth``.
    metadata : dict, optional
        Dataset metadata inserted into the output.
    periods : sequence of float, optional
        Oscillator periods in seconds for pseudo-spectral acceleration.
    damping : float, optional
        Fraction of critical damping. Default is 0.05.
    combination : str, optional
        Horizontal-component combination method. Currently only
        ``"geometric_mean"`` is supported.

    Returns
    -------
    dict
        JSON-serializable ``ShakeScenarioIMT`` dataset.

    Notes
    -----
    Input records are interpreted as ground velocity in m/s.

    PGV is calculated directly from the horizontal velocity components.
    PGA and SA require acceleration, obtained by differentiating independent
    copies of the velocity records.

    ``SA`` represents pseudo-spectral acceleration.

    The output station code is the logical station code from the
    ``ShakeLabStation`` metadata, not the complete waveform SID.
    """
    if not isinstance(stream_collection, StreamCollection):
        raise TypeError(
            "stream_collection must be a StreamCollection instance"
        )

    if combination != "geometric_mean":
        raise ValueError(
            "Only combination='geometric_mean' is currently supported"
        )

    periods = _validate_periods(periods)
    damping = _validate_damping(damping)

    station_index = _station_index(stations)
    waveform_index = _waveform_index(stream_collection)

    output_stations = []

    for station_code, station in station_index.items():
        components = waveform_index.get(station_code)

        if components is None:
            raise KeyError(
                f"No waveform data found for station {station_code!r}"
            )

        if "N" not in components or "E" not in components:
            raise KeyError(
                f"Station {station_code!r} requires both N and E "
                "horizontal components"
            )

        record_n = components["N"]
        record_e = components["E"]

        _check_horizontal_records(
            station_code,
            record_n,
            record_e,
        )

        values_n = _component_imts(
            record_n,
            periods,
            damping,
        )

        values_e = _component_imts(
            record_e,
            periods,
            damping,
        )

        pga = _geometric_mean(
            values_n["PGA"],
            values_e["PGA"],
        )

        pgv = _geometric_mean(
            values_n["PGV"],
            values_e["PGV"],
        )

        sa = np.sqrt(
            values_n["SA"] * values_e["SA"]
        )

        output_stations.append(
            {
                "station_code": station_code,
                "longitude": station["longitude"],
                "latitude": station["latitude"],
                "values": {
                    "PGA": float(pga),
                    "PGV": float(pgv),
                    "SA": [
                        float(value)
                        for value in sa
                    ],
                },
            }
        )

    return {
        "type": SCHEMA_TYPE,
        "schema_version": SCHEMA_VERSION,
        "metadata": dict(metadata or {}),
        "event": _normalize_event(event),
        "imts": {
            "PGA": {
                "units": "m/s2",
                "component": "horizontal",
                "combination": combination,
            },
            "PGV": {
                "units": "m/s",
                "component": "horizontal",
                "combination": combination,
            },
            "SA": {
                "quantity": "pseudo_acceleration",
                "units": "m/s2",
                "component": "horizontal",
                "combination": combination,
                "damping": damping,
                "period_units": "s",
                "periods": [
                    float(period)
                    for period in periods
                ],
            },
        },
        "stations": output_stations,
    }


def read_stream_collection(
    filepath,
    *,
    mseed_backend="auto",
):
    """
    Read waveform data into a ShakeLab StreamCollection.

    Parameters
    ----------
    filepath : str or Path
        Input waveform file.
    mseed_backend : {"auto", "libmseed", "python"}, optional
        MiniSEED backend passed to ShakeLab.

    Returns
    -------
    StreamCollection
        Loaded waveform collection.
    """
    collection = StreamCollection()

    collection.read(
        filepath,
        mseed_backend=mseed_backend,
    )

    return collection


def read_stations(filepath):
    """
    Read a ``ShakeLabStation`` JSON document.

    Parameters
    ----------
    filepath : str or Path
        Station metadata JSON file.

    Returns
    -------
    dict
        Parsed ``ShakeLabStation`` document.
    """
    data = _read_json(filepath)

    if data.get("type") != STATION_SCHEMA_TYPE:
        raise ValueError(
            "Station file must have type='ShakeLabStation'"
        )

    if "stations" not in data:
        raise ValueError(
            "ShakeLabStation file is missing 'stations'"
        )

    return data


def read_events(filepath):
    """
    Read a ``ShakeLabEvent`` JSON document.

    Parameters
    ----------
    filepath : str or Path
        Event catalogue JSON file.

    Returns
    -------
    dict
        Parsed ``ShakeLabEvent`` document.
    """
    data = _read_json(filepath)

    if data.get("type") != EVENT_SCHEMA_TYPE:
        raise ValueError(
            "Event file must have type='ShakeLabEvent'"
        )

    if "events" not in data:
        raise ValueError(
            "ShakeLabEvent file is missing 'events'"
        )

    return data


def get_event(events, event_id):
    """
    Return one event from a ``ShakeLabEvent`` document.

    Parameters
    ----------
    events : dict
        Parsed ``ShakeLabEvent`` document.
    event_id : str
        Event identifier.

    Returns
    -------
    dict
        Matching event dictionary.

    Raises
    ------
    KeyError
        If the event is not found or if the identifier is not unique.
    """
    if not isinstance(events, dict):
        raise TypeError("events must be a dictionary")

    if events.get("type") != EVENT_SCHEMA_TYPE:
        raise ValueError(
            "events must be a ShakeLabEvent document"
        )

    event_list = events.get("events")

    if not isinstance(event_list, list):
        raise ValueError(
            "ShakeLabEvent 'events' must be a list"
        )

    matches = [
        event
        for event in event_list
        if event.get("event_id") == event_id
    ]

    if not matches:
        raise KeyError(
            f"Event not found: {event_id}"
        )

    if len(matches) > 1:
        raise KeyError(
            f"Event identifier is not unique: {event_id}"
        )

    return matches[0]


def write_imt_json(
    data,
    filepath,
    *,
    indent=2,
):
    """
    Write a ``ShakeScenarioIMT`` dictionary to disk.

    Parameters
    ----------
    data : dict
        ``ShakeScenarioIMT`` dictionary.
    filepath : str or Path
        Output JSON file.
    indent : int, optional
        JSON indentation level. Default is 2.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with filepath.open(
        "w",
        encoding="utf-8",
    ) as fobj:
        json.dump(
            data,
            fobj,
            indent=indent,
            ensure_ascii=False,
        )
        fobj.write("\n")


def convert_file(
    waveform_file,
    stations_file,
    events_file,
    event_id,
    output_file,
    *,
    metadata=None,
    periods=DEFAULT_PERIODS,
    damping=DEFAULT_DAMPING,
    combination=DEFAULT_COMBINATION,
    mseed_backend="auto",
):
    """
    Convert waveform and metadata files to ``ShakeScenarioIMT`` JSON.

    Parameters
    ----------
    waveform_file : str or Path
        Input waveform file.
    stations_file : str or Path
        ``ShakeLabStation`` JSON file.
    events_file : str or Path
        ``ShakeLabEvent`` JSON file.
    event_id : str
        Identifier of the event to associate with the IMT dataset.
    output_file : str or Path
        Output JSON file.
    metadata : dict, optional
        Dataset metadata.
    periods : sequence of float, optional
        SA oscillator periods in seconds.
    damping : float, optional
        SA damping ratio.
    combination : str, optional
        Horizontal component combination method.
    mseed_backend : {"auto", "libmseed", "python"}, optional
        MiniSEED backend passed to ShakeLab.

    Returns
    -------
    dict
        Generated ``ShakeScenarioIMT`` dictionary.
    """
    collection = read_stream_collection(
        waveform_file,
        mseed_backend=mseed_backend,
    )

    stations = read_stations(
        stations_file
    )

    events = read_events(
        events_file
    )

    event = get_event(
        events,
        event_id,
    )

    data = build_imt_dict(
        collection,
        stations,
        event,
        metadata=metadata,
        periods=periods,
        damping=damping,
        combination=combination,
    )

    write_imt_json(
        data,
        output_file,
    )

    return data


def _read_json(filepath):
    """
    Read a JSON document from disk.
    """
    filepath = Path(filepath)

    if not filepath.is_file():
        raise FileNotFoundError(
            f"JSON file not found: {filepath}"
        )

    with filepath.open(
        "r",
        encoding="utf-8",
    ) as fobj:
        return json.load(fobj)


def _validate_periods(periods):
    """
    Validate SA oscillator periods.
    """
    periods = np.asarray(
        periods,
        dtype=float,
    )

    if periods.ndim != 1 or periods.size == 0:
        raise ValueError(
            "periods must be a non-empty one-dimensional array"
        )

    if (
        not np.all(np.isfinite(periods))
        or np.any(periods <= 0.0)
    ):
        raise ValueError(
            "periods must contain finite positive values"
        )

    return periods


def _validate_damping(damping):
    """
    Validate SA damping.
    """
    damping = float(damping)

    if not np.isfinite(damping) or damping <= 0.0:
        raise ValueError(
            "damping must be finite and strictly positive"
        )

    return damping


def _normalize_event(event):
    """
    Validate and normalize one event dictionary.
    """
    if not isinstance(event, dict):
        raise TypeError(
            "event must be a dictionary"
        )

    required = (
        "event_id",
        "origin_time",
        "magnitude",
        "longitude",
        "latitude",
        "depth",
    )

    missing = [
        key
        for key in required
        if key not in event
    ]

    if missing:
        raise ValueError(
            "Missing event keys: {0}".format(
                ", ".join(missing)
            )
        )

    output = {
        "event_id": str(
            event["event_id"]
        ),
        "origin_time": str(
            event["origin_time"]
        ),
        "magnitude": float(
            event["magnitude"]
        ),
        "longitude": float(
            event["longitude"]
        ),
        "latitude": float(
            event["latitude"]
        ),
        "depth": float(
            event["depth"]
        ),
    }

    for key in (
        "magnitude",
        "longitude",
        "latitude",
        "depth",
    ):
        if not np.isfinite(output[key]):
            raise ValueError(
                f"Event {key} must be finite"
            )

    return output


def _parse_sid(sid):
    """
    Parse an FDSN source identifier.

    Returns
    -------
    tuple
        Network, station, location, and channel codes.
    """
    if not isinstance(sid, str):
        raise TypeError(
            "Waveform SID must be a string"
        )

    parts = sid.split(".")

    if len(parts) != 4:
        raise ValueError(
            f"Invalid FDSN source identifier: {sid!r}"
        )

    network, station, location, channel = parts

    if not station:
        raise ValueError(
            f"Missing station code in SID: {sid!r}"
        )

    if not channel:
        raise ValueError(
            f"Missing channel code in SID: {sid!r}"
        )

    return (
        network,
        station,
        location,
        channel,
    )


def _station_index(stations):
    """
    Validate and index ShakeLabStation metadata by station code.
    """
    if not isinstance(stations, dict):
        raise TypeError(
            "stations must be a dictionary"
        )

    if stations.get("type") != STATION_SCHEMA_TYPE:
        raise ValueError(
            "stations must be a ShakeLabStation document"
        )

    station_list = stations.get("stations")

    if not isinstance(station_list, list):
        raise ValueError(
            "ShakeLabStation 'stations' must be a list"
        )

    index = {}

    for item in station_list:
        if not isinstance(item, dict):
            raise TypeError(
                "Each station definition must be a dictionary"
            )

        required = (
            "network_code",
            "station_code",
            "location_code",
            "longitude",
            "latitude",
            "channel_codes",
        )

        missing = [
            key
            for key in required
            if key not in item
        ]

        if missing:
            raise ValueError(
                "Missing station keys: {0}".format(
                    ", ".join(missing)
                )
            )

        station_code = str(
            item["station_code"]
        )

        if not station_code:
            raise ValueError(
                "station_code cannot be empty"
            )

        if station_code in index:
            raise ValueError(
                f"Duplicate station code: {station_code!r}"
            )

        longitude = float(
            item["longitude"]
        )
        latitude = float(
            item["latitude"]
        )

        if not np.isfinite(
            longitude
        ) or not np.isfinite(
            latitude
        ):
            raise ValueError(
                f"Invalid coordinates for station "
                f"{station_code!r}"
            )

        channel_codes = item["channel_codes"]

        if not isinstance(channel_codes, list):
            raise TypeError(
                f"channel_codes must be a list for station "
                f"{station_code!r}"
            )

        index[station_code] = {
            "network_code": str(
                item["network_code"]
            ),
            "location_code": str(
                item["location_code"]
            ),
            "longitude": longitude,
            "latitude": latitude,
            "channel_codes": list(
                channel_codes
            ),
        }

    return index


def _waveform_index(stream_collection):
    """
    Index waveform records by station code and orientation.
    """
    index = {}

    for stream in stream_collection:
        (
            _network,
            station_code,
            _location,
            channel,
        ) = _parse_sid(
            stream.sid
        )

        component = channel[-1].upper()

        if component not in (
            "N",
            "E",
            "Z",
        ):
            continue

        record = _stream_record(
            stream
        )

        station = index.setdefault(
            station_code,
            {},
        )

        if component in station:
            raise ValueError(
                f"Duplicate {component} component for station "
                f"{station_code!r}"
            )

        station[component] = record

    return index


def _stream_record(stream):
    """
    Return one continuous Record from a ShakeLab Stream.
    """
    if len(stream) == 0:
        raise ValueError(
            f"Empty stream: {stream.sid!r}"
        )

    if len(stream) == 1:
        return stream[0].copy()

    stream = stream.copy()
    stream.sort()

    record = stream[0].copy()

    for segment in stream[1:]:
        merged = record.append(
            segment.copy(),
            enforce=False,
        )

        if not merged:
            raise ValueError(
                f"Stream {stream.sid!r} contains non-contiguous "
                "records"
            )

    return record


def _check_horizontal_records(
    station_code,
    record_n,
    record_e,
):
    """
    Check consistency of the two horizontal records.
    """
    if not np.isclose(
        record_n.delta,
        record_e.delta,
        rtol=1e-9,
        atol=0.0,
    ):
        raise ValueError(
            f"Inconsistent sampling intervals for station "
            f"{station_code!r}"
        )

    if record_n.nsamp != record_e.nsamp:
        raise ValueError(
            f"Inconsistent sample counts for station "
            f"{station_code!r}"
        )

    if (
        record_n.time is not None
        and record_e.time is not None
        and record_n.time != record_e.time
    ):
        raise ValueError(
            f"Inconsistent start times for station "
            f"{station_code!r}"
        )


def _component_imts(
    velocity_record,
    periods,
    damping,
):
    """
    Compute PGV, PGA and pseudo-spectral acceleration for one component.
    """
    pgv = velocity_record.peak

    acceleration_record = velocity_record.copy()

    acceleration_record.differentiate(
        method="grad",
    )

    pga = acceleration_record.peak

    spectrum = ResponseSpectrum(
        acceleration_record,
        damping=damping,
    )

    spectrum.compute(
        periods
    )

    return {
        "PGA": float(pga),
        "PGV": float(pgv),
        "SA": np.asarray(
            spectrum.psa,
            dtype=float,
        ),
    }


def _geometric_mean(
    value_1,
    value_2,
):
    """
    Return the geometric mean of two non-negative values.
    """
    value_1 = float(value_1)
    value_2 = float(value_2)

    if value_1 < 0.0 or value_2 < 0.0:
        raise ValueError(
            "Geometric mean requires non-negative IMT values"
        )

    return np.sqrt(
        value_1 * value_2
    )
