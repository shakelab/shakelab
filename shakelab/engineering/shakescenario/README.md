# ShakeScenario

ShakeScenario is a client-server application for running rapid seismic
damage and impact scenarios using the ShakeLab engineering modules.

It provides:

-   a background TCP server that manages persistent calculation jobs;
-   a Python/command-line client for submitting and managing scenarios;
-   support for multiple self-contained scenario models;
-   configurable ground-motion providers and asset-to-provider
    assignments;
-   persistent, reproducible job inputs and outputs;
-   integration with the ShakeScenario WebUI.

The scientific calculation is implemented by the underlying ShakeLab
engineering modules. ShakeScenario provides the service,
model-management, job-management, and persistence layers around those
modules.

## 1. Architecture

ShakeScenario is organized in three levels.

1.  **Engineering library**

    The scientific calculation is provided by modules under
    `shakelab.engineering`, in particular exposure, taxonomy, fragility,
    ground motion, and impact. These modules can also be used directly
    from Python without running ShakeScenario.

2.  **ShakeScenario service**

    `shakeserver.py` exposes scenario calculations as persistent jobs
    through a TCP service. `shakeclient.py` provides the corresponding
    Python client and command-line interface.

3.  **WebUI**

    The WebUI uses the ShakeScenario service rather than reimplementing
    the scientific calculation.

A scenario model is a self-contained bundle containing the exposure,
taxonomy tree, fragility models, ground-motion configuration, and
optional geometry used by the WebUI.

## 2. Requirements

ShakeScenario is part of ShakeLab.

Requirements:

-   Python 3.10 or newer;
-   ShakeLab installed and importable;
-   the scientific dependencies required by the selected ShakeLab
    engineering modules;
-   at least one valid ShakeScenario model.

The server and client communicate using a length-prefixed JSON TCP
protocol.

## 3. Server configuration

The server is started with a JSON configuration file.

Example:

``` json
{
  "schema_version": "1.0.0",
  "server": {
    "host": "127.0.0.1",
    "port": 6000,
    "workers": 4
  },
  "paths": {
    "db": "./db/shakescenario.db",
    "workdir": "./runs",
    "model_root": "./models"
  },
  "defaults": {
    "model_id": "fvg_20260630",
    "impact_config": {
      "uncertainty_mode": "lognormal",
      "output": "state",
      "typology_weighting": "count",
      "normalize_asset_probabilities": true,
      "missing_taxonomy": "raise",
      "no_damage_key": "D0",
      "include_typology_breakdown": true
    }
  }
}
```

### 3.1 Paths

`db`
:   SQLite database used to store persistent job information.

`workdir`
:   Root directory where job directories and calculation products are
    written.

`model_root`
:   Root directory containing the available scenario models.

### 3.2 Server settings

`host`
:   TCP interface on which the server listens.

`port`
:   TCP port.

`workers`
:   Number of background worker threads available for calculations.

Command-line values for `--host`, `--port`, and `--workers` override the
corresponding values in the configuration file.

### 3.3 Defaults

`model_id`
:   Scenario model used when the client does not explicitly select one.

`impact_config`
:   Default configuration passed to the engineering impact calculation.

Ground-motion configuration is **not** a server default. It belongs to
the selected scenario model. This keeps the service independent of the
specific ground-motion methodology used by a model.

## 4. Scenario models

Each model is stored in:

``` text
model_root/<model_id>/
```

For example:

``` text
models/
└── fvg_20260630/
    ├── manifest.json
    ├── exposure.json
    ├── taxonomy_tree.json
    ├── geometry.geojson
    ├── groundmotion.json
    └── fragility/
        ├── fragility_faravelli.json
        └── fragility_rosti.json
```

The directory name is the `model_id`.

### 4.1 Model manifest

Each model must contain `manifest.json`.

Example:

``` json
{
  "schema_version": "1.0.0",
  "files": {
    "exposure": "exposure.json",
    "taxonomy_tree": "taxonomy_tree.json",
    "fragility": [
      "fragility/fragility_faravelli.json",
      "fragility/fragility_rosti.json"
    ],
    "ground_motion": "groundmotion.json",
    "geometry": "geometry.geojson"
  }
}
```

The scientific files declared by the manifest are resolved relative to
the model directory.

The current service requires:

-   `exposure`;
-   `taxonomy_tree`;
-   one or more `fragility` files;
-   `ground_motion`.

`geometry` is used by higher-level applications such as the WebUI and is
not part of the core impact calculation.

## 5. Ground-motion configuration

Ground motion is configured at model level in `groundmotion.json`.

This layer is separate from the general-purpose ground-motion modelling
code under `shakelab.gmmodel`. The engineering ground-motion layer
configures providers and determines which configured provider is used
for each exposure asset.

A model may define several configured providers. Provider instances and
their contexts are shared by all assets assigned to them.

### 5.1 Single provider for all assets

The simplest configuration uses one default provider:

``` json
{
  "type": "ShakeLabGroundMotion",
  "schema_version": "1.0.0",
  "providers": [
    {
      "id": "regional_gmpe",
      "provider": "gmpe",
      "default": true,
      "config": {
        "gmpe_name": "BragatoSlejko2005",
        "distance_approx": "ellipsoid"
      }
    }
  ]
}
```

All assets without an explicit assignment use the provider marked as
`default`.

### 5.2 Multiple providers

Providers may contain explicit asset assignments. The same configured
provider can therefore be reused by many assets without duplicating its
configuration.

Assignments may also carry provider-specific asset parameters, such as a
station code for an observed-motion provider.

The intensity measure type (IMT) is not assigned in the ground-motion
file. It is determined by the fragility model used for the asset
taxonomy.

The current assignment model resolves one provider for each asset.
Combining multiple simultaneous providers for the same asset is
intentionally outside the current v1 semantics.

## 6. Damage probabilities

Fragility curves are interpreted as damage-state exceedance
probabilities:

``` text
P(DS >= D1)
P(DS >= D2)
...
P(DS >= Dn)
```

ShakeLab can convert these cumulative probabilities into mutually
exclusive damage-state probabilities:

``` text
P(D0) = 1 - P(DS >= D1)

P(Dk) = P(DS >= Dk) - P(DS >= D(k+1))
        for k = 1, ..., n-1

P(Dn) = P(DS >= Dn)
```

`D0` represents no damage.

The impact configuration supports:

`output = "exceed"`
:   Damage probabilities are reported as exceedance probabilities
    `D1..Dn`.

`output = "state"`
:   Damage probabilities are reported as mutually exclusive states
    `D0..Dn`.

Expected counts are always calculated using the mutually exclusive state
representation `D0..Dn`.

There is no additional `GT_LAST` damage state in the current convention.
The last configured damage state `Dn` is the most severe state
represented by the damage scale.

## 7. Starting the server

From a working directory containing the server configuration:

``` bash
python3 /path/to/shakelab/engineering/shakescenario/shakeserver.py \
    --config config.json
```

Optional command-line overrides include:

``` text
--host
--port
--workers
```

At startup the server validates its configuration and the default model,
initializes the database, and reports the listening address and
available models.

A typical interactive startup looks like:

``` text
ShakeScenario server
  host:       127.0.0.1
  port:       6000
  workers:    4
  model root: /path/to/models
  default:    fvg_20260630
  models:     fvg_20260630

Server ready. Press Ctrl+C to stop.
```

Press `Ctrl+C` to stop an interactively running server.

## 8. Using the client

The command-line client is:

``` bash
python3 shakeclient.py
```

The default connection is `127.0.0.1:6000`. Use `--host` and `--port` to
connect elsewhere.

### 8.1 Check server connectivity

``` bash
python3 shakeclient.py ping
```

### 8.2 Submit a scenario

``` bash
python3 shakeclient.py submit \
    --mag 5.2 \
    --lon 13.25 \
    --lat 46.25 \
    --depth 10
```

Optional event metadata include:

``` text
--origin-time
--tag
```

For example:

``` bash
python3 shakeclient.py submit \
    --tag test_fvg_m52 \
    --origin-time "2026-07-08T12:32:15Z" \
    --mag 5.2 \
    --lon 13.25 \
    --lat 46.25 \
    --depth 10
```

The selected model determines the exposure, taxonomy, fragility, and
ground-motion configuration. The client therefore does not select a GMPE
or other ground-motion provider directly.

### 8.3 Submit with a specific model

``` bash
python3 shakeclient.py submit \
    --model-id fvg_20260630 \
    --mag 5.2 \
    --lon 13.25 \
    --lat 46.25 \
    --depth 10
```

### 8.4 Submit a JSON request

A request can also be supplied using:

``` bash
python3 shakeclient.py submit --request-json request.json
```

A standard request has the form:

``` json
{
  "models": {
    "model_id": "fvg_20260630"
  },
  "scenario": {
    "event": {
      "magnitude": 5.2,
      "hypocentre": {
        "longitude": 13.25,
        "latitude": 46.25,
        "elevation": -10000.0
      }
    }
  }
}
```

`scenario.ground_motion` is not part of the current standard request.
Ground motion is defined by the selected model.

### 8.5 List available models

``` bash
python3 shakeclient.py models list
```

### 8.6 List jobs

``` bash
python3 shakeclient.py list
```

### 8.7 Get job details

``` bash
python3 shakeclient.py get 1
```

### 8.8 Wait for completion

``` bash
python3 shakeclient.py wait 1 --timeout 120 --interval 2
```

### 8.9 Delete a job

``` bash
python3 shakeclient.py delete 1
```

To remove associated job files as well:

``` bash
python3 shakeclient.py delete 1 --purge
```

### 8.10 Reset the database

Development/administrative operation:

``` bash
python3 shakeclient.py reset --yes-i-know
```

To remove job directories as well:

``` bash
python3 shakeclient.py reset --yes-i-know --purge
```

## 9. Job lifecycle

A job can have one of the following states:

``` text
queued
running
completed
failed
canceled
```

Jobs are persisted in the SQLite database. Calculations are executed by
the server worker pool.

## 10. Job directories and artifacts

Each submitted job is assigned a directory such as:

``` text
runs/job_000001/
```

The current job layout is:

``` text
job_000001/
├── manifest.json
├── request.json
├── results/
│   ├── impact_assets.json
│   └── impact_summary.json
└── logs/
    ├── job.log
    ├── warnings.json
    └── error.log
```

`manifest.json`
:   Job metadata, status, model identifier, ground-motion provider
    provenance, and relative artifact paths.

`request.json`
:   Resolved request used by the calculation.

`results/impact_assets.json`
:   Per-asset impact results, including ground motion used by the
    calculation and damage results.

`results/impact_summary.json`
:   Lightweight aggregate statistics.

`logs/job.log`
:   Job execution log.

`logs/warnings.json`
:   Structured warnings.

`logs/error.log`
:   Error information for failed jobs.

Artifact paths stored in the manifest are relative to the job directory.

## 11. Ground-motion provenance

The impact result records the configured `provider_id` associated with
each ground-motion value. This makes it possible to identify which
configured provider produced the intensity measure used for a given
asset.

The job manifest also records lightweight model-level ground-motion
provenance, including the configured provider identifiers and the
default provider identifier.

The complete ground-motion configuration remains part of the selected
model rather than being duplicated into every job request.

## 12. Parameter precedence

For request parameters that support defaults, precedence is:

``` text
client payload
    overrides
server defaults
```

The server default `model_id` is used if no model is selected by the
client.

Ground-motion provider selection is different: it is resolved by the
ground-motion configuration of the selected model and its asset
assignments, not by client/server parameter precedence.

## 13. Typical workflow

1.  Prepare a scenario model under `model_root`.
2.  Add its `manifest.json` and `groundmotion.json`.
3.  Configure the ShakeScenario server.
4.  Start the server.
5.  Check connectivity with `ping`.
6.  Submit an earthquake scenario.
7.  Monitor the job with `list`, `get`, or `wait`.
8.  Inspect or retrieve the generated artifacts.
9.  Use the same service through the WebUI when required.

## 14. Design principles

ShakeScenario follows a few deliberate architectural rules:

-   scientific functionality remains in reusable `shakelab.engineering`
    modules;
-   the TCP service builds on those library functions rather than
    duplicating them;
-   the WebUI builds on the service layer;
-   scenario models are self-contained and portable;
-   exposure data are not modified merely to encode ground-motion
    routing;
-   configured ground-motion providers are shared across assigned
    assets;
-   provider-specific asset parameters belong to ground-motion
    assignments;
-   the IMT is determined by fragility, not duplicated in ground-motion
    configuration;
-   job outputs use relative artifact paths and preserve calculation
    provenance;
-   model and request formats are versioned explicitly.

This separation allows the same engineering functionality to be used
directly as a Python library, through the ShakeScenario service, or
through the WebUI.
