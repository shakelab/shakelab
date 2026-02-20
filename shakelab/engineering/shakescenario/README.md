# ShakeScenario

**ShakeScenario** is a lightweight TCP client--server service for
running job-based rapid seismic damage and impact scenarios within the
ShakeLab engineering framework.

The module is designed to live inside:

    shakelab.engineering

and provides a clean separation between:

-   scenario submission (client side),
-   parallel execution (server side),
-   job persistence (sqlite3),
-   artifact storage (filesystem).

------------------------------------------------------------------------

## Architecture

ShakeScenario follows a simple job-oriented architecture.

Client (CLI): - submits scenario jobs, - queries status and metadata, -
manages lifecycle (cancel, delete, reset).

Server: - accepts TCP connections, - persists jobs in sqlite3, -
executes jobs using a worker pool, - stores per-job artifacts in a
dedicated directory.

Communication is performed via:

-   TCP sockets
-   Length-prefixed JSON messages (8-byte big-endian header)

------------------------------------------------------------------------

## Repository Layout

    shakescenario/
    ├── shakescenario-server.py
    ├── shakescenario.py
    ├── protocol.py
    ├── database.py
    ├── models.py
    └── README.md

This simplified layout is intentional, since the module is contained
within `shakelab.engineering`.

------------------------------------------------------------------------

## Running the Server

Start the server:

``` bash
python3 shakescenario-server.py   --host 127.0.0.1   --port 6000   --db shakescenario.db   --workdir ./runs   --workers 4
```

Options:

-   `--host` : bind address
-   `--port` : TCP port
-   `--db` : sqlite database file
-   `--workdir` : directory for job artifacts
-   `--workers` : number of parallel workers

------------------------------------------------------------------------

## Using the Client

Ping the server:

``` bash
python3 shakescenario.py ping
```

Submit a scenario:

``` bash
python3 shakescenario.py submit   --mag 5.2   --lon 12.34   --lat 45.67   --depth 10
```

List jobs:

``` bash
python3 shakescenario.py list
```

Get job details:

``` bash
python3 shakescenario.py get 1
```

Wait for completion:

``` bash
python3 shakescenario.py wait 1 --timeout 60
```

------------------------------------------------------------------------

## Job Persistence

Jobs are stored in a sqlite3 database with the following states:

-   queued
-   running
-   completed
-   failed
-   canceled

Each job has:

-   creation timestamp
-   optional start/end timestamps
-   parameter JSON
-   result metadata JSON
-   error field (if any)

------------------------------------------------------------------------

## Artifact Storage

Each job produces a dedicated directory:

    runs/
    └── job_000001/
        ├── meta.json
        └── result_manifest.json

Future versions may include additional outputs such as:

-   impact CSV
-   GeoJSON
-   logs
-   ground motion summaries

------------------------------------------------------------------------

## Protocol (v1)

Each TCP message consists of:

1.  8-byte unsigned big-endian length
2.  UTF-8 JSON payload

Request structure:

``` json
{
  "v": 1,
  "op": "submit",
  "req_id": "uuid",
  "payload": { ... }
}
```

Response structure:

``` json
{
  "v": 1,
  "req_id": "uuid",
  "ok": true,
  "result": { ... }
}
```

Errors are returned with:

``` json
{
  "v": 1,
  "req_id": "uuid",
  "ok": false,
  "error": {
    "code": "ERROR_CODE",
    "message": "Description"
  }
}
```

------------------------------------------------------------------------

## Integration Point

The actual scenario computation must be implemented in:

    ShakeScenarioServer._run_job()

This is the location where ShakeLab ground motion, fragility, and impact
modules should be invoked.

------------------------------------------------------------------------

## Design Goals

-   Minimal dependencies (stdlib only)
-   Deterministic persistence (sqlite3)
-   Clean protocol contract
-   Parallel execution via worker pool
-   Reproducible artifact layout
-   Future extensibility (download, authentication, clustering)
