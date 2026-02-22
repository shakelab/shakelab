# ShakeScenario

ShakeScenario is a client--server application for running rapid seismic
damage and impact scenarios using the ShakeLab engineering modules.

It allows you to:

-   Start a background server that manages jobs
-   Submit earthquake scenarios from a command line client
-   Track execution status
-   Retrieve results and artifacts
-   Manage multiple model configurations on the server

This document explains how to install, configure, start, and use the
system from a user perspective.

------------------------------------------------------------------------

# 1. System Overview

ShakeScenario consists of:

SERVER shakescenario-server.py

CLIENT (CLI) shakescenario.py

The server performs the actual scenario computation. The client is used
to communicate with the server.

All computations are executed as jobs and stored persistently.

------------------------------------------------------------------------

# 2. Preparing the Environment

ShakeScenario is part of:

    shakelab.engineering

Make sure:

-   Python 3.10+ is available
-   ShakeLab is installed and importable
-   Required models are prepared (see Section 4)

No external dependencies beyond the standard library and ShakeLab are
required.

------------------------------------------------------------------------

# 3. Server Configuration

The server requires a JSON configuration file.

Example: config.json

{ "schema_version": "1.0.0", "paths": { "db": "./shakescenario.db",
"workdir": "./runs", "model_root": "./models" }, "defaults": {
"model_id": "ne_italy_default", "ground_motion": { "provider": "gmpe",
"gmpe_name": "BragatoSlejko2005", "distance_approx": "ellipsoid" },
"impact_config": { "uncertainty_mode": "lognormal",
"typology_weighting": "count", "missing_taxonomy": "raise",
"no_damage_key": "D0", "tail_key": "GT_LAST" } } }

Meaning:

db SQLite database file storing job history

workdir Directory where job results are written

model_root Directory containing model subdirectories

defaults Default parameters used if not specified by the client

------------------------------------------------------------------------

# 4. Preparing Models

Each model must be placed in:

    model_root / model_id

Example:

    models/
        ne_italy_default/
            exposure.json
            fragility.json
            taxonomy_tree.json

Each model directory MUST contain:

-   exposure.json
-   fragility.json
-   taxonomy_tree.json

The model_id is the directory name.

------------------------------------------------------------------------

# 5. Starting the Server

From the project root:

    python3 shakelab/engineering/shakescenario/shakescenario-server.py         --config config.json

Optional parameters:

    --host      Default: 127.0.0.1
    --port      Default: 6000
    --workers   Number of parallel jobs

The server will:

-   Validate configuration
-   Validate model directories
-   Open the database
-   Start listening for connections

------------------------------------------------------------------------

# 6. Using the Client

All commands are executed using:

    python3 shakescenario.py

------------------------------------------------------------------------

## 6.1 Check Server Connectivity

    python3 shakescenario.py ping

If successful, the server responds.

------------------------------------------------------------------------

## 6.2 Submit a Scenario (Minimal)

    python3 shakescenario.py submit         --mag 5.2         --lon 12.34         --lat 45.67         --depth 10

This uses default model and default ground motion configuration.

------------------------------------------------------------------------

## 6.3 Submit with Custom Model

    python3 shakescenario.py submit         --mag 5.2         --lon 12.34         --lat 45.67         --depth 10         --model-id ne_italy_default

------------------------------------------------------------------------

## 6.4 List Available Models

    python3 shakescenario.py models list

------------------------------------------------------------------------

## 6.5 List Jobs

    python3 shakescenario.py list

------------------------------------------------------------------------

## 6.6 Get Job Details

    python3 shakescenario.py get 1

Shows:

-   Status
-   Parameters
-   Error (if failed)
-   Result metadata

------------------------------------------------------------------------

## 6.7 Wait for Completion

    python3 shakescenario.py wait 1 --timeout 60

------------------------------------------------------------------------

## 6.8 Reset Database (Development)

WARNING: Deletes all jobs.

    python3 shakescenario.py reset --yes-i-know

------------------------------------------------------------------------

# 7. Understanding Job Results

Each job creates a directory:

    runs/job_000001/

Containing:

request_resolved.json Full parameters used after merging defaults

meta.json Execution metadata

impact_result.json Final computed impact results

result_manifest.json List of generated artifacts

error.txt Only present if the job failed

------------------------------------------------------------------------

# 8. Job Lifecycle

A job can be:

queued running completed failed canceled

Status is visible using:

    shakescenario.py list
    shakescenario.py get <id>

------------------------------------------------------------------------

# 9. Parameter Precedence

When submitting a job:

    client payload overrides
    server defaults override
    internal hardcoded defaults

------------------------------------------------------------------------

# 10. Typical Workflow

1.  Prepare models in model_root
2.  Configure config.json
3.  Start server
4.  Submit scenario
5.  Monitor job
6.  Retrieve results from runs/job_xxxxxx

------------------------------------------------------------------------

ShakeScenario is designed to provide reproducible, transparent, and
traceable scenario calculations with persistent storage and clear
artifact management.
