from shakelab.engineering.shakescenario.imt import convert_file


WAVEFORM_FILE = "waveforms_friuli_76.mseed"
STATIONS_FILE = "stations.json"
EVENTS_FILE = "event_friuli_76.json"
OUTPUT_FILE = "imt_friuli_76.json"

EVENT_ID = "friuli_1976"

METADATA = {
    "name": "Friuli 1976 synthetic ground-motion field",
    "description": (
        "Ground-motion intensity measures derived from SPEED "
        "synthetic waveforms."
    ),
    "source": "SPEED",
    "crs": "EPSG:4326",
}

PERIODS = (
    0.1,
    0.2,
    0.3,
    0.5,
    1.0,
    2.0,
)

DAMPING = 0.05


imt_data = convert_file(
    waveform_file=WAVEFORM_FILE,
    stations_file=STATIONS_FILE,
    events_file=EVENTS_FILE,
    event_id=EVENT_ID,
    output_file=OUTPUT_FILE,
    metadata=METADATA,
    periods=PERIODS,
    damping=DAMPING,
)


print(
    f"Generated {OUTPUT_FILE} "
    f"with {len(imt_data['stations'])} stations."
)