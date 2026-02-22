# ShakeScenario

**ShakeScenario** is a lightweight TCP-based client--server service for
executing rapid seismic impact and damage scenarios within the ShakeLab
engineering framework.

The module lives inside:

    shakelab.engineering

and provides separation between:

-   Scenario submission (client side)
-   Parallel execution (server side)
-   Job persistence (sqlite3)
-   Model registry (model_id → directory resolution)
-   Artifact storage (filesystem)
-   Deterministic job outputs

------------------------------------------------------------------------

# Architecture Overview

ShakeScenario follows a job-oriented architecture.

Client (CLI): - Submits scenario jobs - Queries status and metadata -
Waits for completion - Cancels, deletes, or resets jobs - Lists
available model_id on the server

Server: - Accepts TCP connections - Validates and merges payload with
server defaults - Resolves model_id to model directories - Executes jobs
using a worker pool - Persists job state in sqlite3 - Writes per-job
artifacts to a dedicated directory - Stores debug artifacts (error.txt)

Communication protocol: - TCP sockets - Length-prefixed JSON messages
(8-byte big-endian header)

------------------------------------------------------------------------

# Repository Layout

    shakescenario/
    ├── shakescenario-server.py
    ├── shakescenario.py
    ├── protocol.py
    ├── database.py
    ├── models.py
    └── README.md

------------------------------------------------------------------------

# Server Configuration (config.json)

The server requires a JSON configuration file.

Example:

{ "schema_version": "1.0.0", "paths": { "db": "./shakescenario.db",
"workdir": "./runs", "model_root": "./models" }, "server": {},
"defaults": { "model_id": "ne_italy_default", "ground_motion": {
"provider": "gmpe", "gmpe_name": "BragatoSlejko2005", "distance_approx":
"ellipsoid" }, "impact_config": { "uncertainty_mode": "lognormal",
"typology_weighting": "count", "missing_taxonomy": "raise",
"no_damage_key": "D0", "tail_key": "GT_LAST" } } }

Precedence:

    payload > config.defaults > hardcoded

------------------------------------------------------------------------

# Model Registry (model_id)

Models are resolved using:

    model_root / model_id

Each model directory must contain:

    exposure.json
    fragility.json
    taxonomy_tree.json

Validation includes:

-   Strict model_id regex
-   Prevention of path traversal
-   Existence checks

List models:

    python3 shakescenario.py models list

------------------------------------------------------------------------

# Running the Server

From repository root:

    python3 shakelab/engineering/shakescenario/shakescenario-server.py         --config shakelab/engineering/shakescenario/config.json

Options:

-   --host
-   --port
-   --workers
-   --config (required)

------------------------------------------------------------------------

# Client Usage

Ping:

    python3 shakescenario.py ping

Submit:

    python3 shakescenario.py submit         --mag 5.2         --lon 12.34         --lat 45.67         --depth 10

Override model:

    python3 shakescenario.py submit         --mag 5.2         --lon 12.34         --lat 45.67         --depth 10         --model-id ne_italy_default

List jobs:

    python3 shakescenario.py list

Get job:

    python3 shakescenario.py get 1

Wait:

    python3 shakescenario.py wait 1 --timeout 60

Reset (dev):

    python3 shakescenario.py reset --yes-i-know

------------------------------------------------------------------------

# Job Lifecycle

States:

-   queued
-   running
-   completed
-   failed
-   canceled

Reset also resets AUTOINCREMENT.

------------------------------------------------------------------------

# Artifact Layout

    runs/
    └── job_000001/
        ├── request_resolved.json
        ├── meta.json
        ├── impact_result.json
        ├── result_manifest.json
        └── error.txt (if failed)

request_resolved.json: Full merged payload.

impact_result.json: Output of compute_impact_scenario +
save_impact_result.

result_manifest.json: List of artifacts.

error.txt: Debug info if job fails.

------------------------------------------------------------------------

# Protocol (v1)

Each message:

1.  8-byte big-endian length
2.  UTF-8 JSON payload

Request:

{ "v": 1, "op": "submit", "req_id": "uuid", "payload": { ... } }

Response (success):

{ "v": 1, "req_id": "uuid", "ok": true, "result": { ... } }

Response (error):

{ "v": 1, "req_id": "uuid", "ok": false, "error": { "code":
"ERROR_CODE", "message": "Description" } }

------------------------------------------------------------------------

# Integration Point

Scenario computation is implemented in:

    ShakeScenarioServer._run_job()

This integrates:

-   ExposureModel
-   FragilityCollection
-   TaxonomyTree
-   GroundMotionContext
-   ImpactConfig
-   compute_impact_scenario
-   save_impact_result

------------------------------------------------------------------------

# Design Goals

-   Minimal dependencies
-   Deterministic persistence
-   Reproducible artifacts
-   Safe model resolution
-   Explicit configuration
-   Developer-friendly debugging (error.txt)
