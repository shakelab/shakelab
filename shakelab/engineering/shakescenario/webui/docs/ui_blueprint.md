# ShakeScenario WebUI Blueprint

Version: 0.1  
Scope: first design reference for the ShakeScenario WebUI  
Status: working draft

---

## 1. Purpose

The ShakeScenario WebUI is a browser-based working interface for
seismologists and civil-protection operators.

It must be treated as a technical application, not as a generic website.
Its purpose is to support operational and scientific work on earthquake
damage scenarios, including:

- inspecting calculated ShakeScenario runs;
- selecting and comparing scenarios;
- visualizing geospatial layers;
- reviewing event metadata;
- reviewing damage and impact summaries;
- accessing output files, logs, and reports;
- launching or preparing new scenarios in future versions.

The interface should remain lightweight, stable, understandable, and easy to
maintain over the long term.

The design priority is operational clarity, not visual decoration.

---

## 2. General design philosophy

The WebUI should behave more like a desktop GIS-style application delivered
through a browser than like a conventional website.

The user should not feel that they are navigating between independent pages.
The main working context should remain visible as much as possible, especially
the map and the selected scenario.

The old Protezione Civile regionale dashboard is not a graphical reference.
It is only a functional reference for understanding which information may be
important in an operational context, especially the list of municipalities with
damage.

The internal ShakeScenario WebUI draft is the preferred visual reference:
clean, technical, modern, neutral colors, light cards, compact panels, and
minimal visual noise.

---

## 3. Technology stack

Preferred stack:

- Apache as public web server;
- Gunicorn as WSGI server;
- Flask as lightweight Python backend;
- HTML for structure;
- CSS for layout and style;
- vanilla JavaScript for frontend logic;
- Leaflet for map visualization.

Avoid unless strictly necessary:

- Django;
- React;
- Vue;
- Angular;
- Node.js build systems;
- Webpack/Vite;
- SQL ORM layers;
- unnecessary frontend component libraries;
- Docker as a mandatory runtime dependency.

The guiding principle is to keep the number of moving parts small.

---

## 4. Deployment concept

Recommended production architecture:

```text
Browser
  |
  | HTTPS / HTTP
  v
Apache
  |
  | reverse proxy for /api/*
  v
Gunicorn
  |
  | WSGI
  v
Flask WebUI backend
  |
  v
ShakeScenario data directory
```

Apache should be responsible for:

- public network access;
- TLS/HTTPS;
- static files, if useful;
- reverse proxying API requests to Gunicorn.

Gunicorn should run locally, for example:

```text
127.0.0.1:8000
```

and should be managed by `systemd`.

Flask should provide a small REST API and serve the main HTML template.

---

## 5. Data root

The WebUI reads ShakeScenario data from a configurable root directory.

Recommended environment variable:

```bash
SHAKESCENARIO_HOME=/path/to/shakescenario_data
```

This allows the same code to run both on a local development machine and on
the production server.

Current expected structure:

```text
shakescenario_data/
  db/
    shakescenario.db

  logs/
    server.log
    ... [still to implement]

  models/
    model_id/
      exposure.json
      taxonomy_tree.json
      fragility/
        fragility_1.json
        fragility_2.json
        ...
      polygons.json

  runs/
    job_000001/
      meta.json
      impact_result.json
      request_resolved.json
      result_manifest.json
      error.txt
```

This structure is still subject to future refinement.

---

## 6. Main layout

Preferred layout:

```text
+------------------------------------------------------------------+
| Header                                                           |
+------------------------------------------------------------------+

+--------------------------------------+---------------------------+
|                                      | Event details             |
|                                      |                           |
|                                      +---------------------------+
|                Map                   | Municipalities / impact   |
|                                      | summary                   |
|                                      |                           |
+--------------------------------------+---------------------------+

+------------------------------------------------------------------+
| Scenario search filters                                          |
+------------------------------------------------------------------+

+------------------------------------------------------------------+
| Calculated scenario archive                                      |
+------------------------------------------------------------------+

+------------------------------------------------------------------+
| Status footer                                                    |
+------------------------------------------------------------------+
```

The map is the operational core of the interface.

The right-side column is dedicated to information about the currently selected
scenario. It is divided into:

1. event/scenario details;
2. municipalities with damage and impact summary.

Below the map area, the user finds:

1. filters for searching calculated scenarios;
2. the archive table of calculated runs.

This keeps the main visual analysis area at the top while still allowing a
complete archive of scenarios to be available without changing page.

---

## 7. Main UI components

### 7.1 Header

Purpose:

- identify the application;
- expose global actions;
- provide compact system-level controls.

Initial elements:

- application name: ShakeScenario;
- short subtitle;
- refresh button;
- new scenario button.

Possible future elements:

- user/session information;
- settings;
- language selector;
- system status indicator;
- link to documentation;
- link to logs or administration page.

The header should remain compact and should not reduce the available vertical
space for the map.

---

### 7.2 MapPanel

Purpose:

- display the selected scenario spatially;
- visualize damage, shaking, exposure, and future result layers;
- provide GIS-style controls.

Initial elements:

- Leaflet map;
- layer selector;
- map subtitle/status;
- collapsible legend;
- optional map tools.

Initial thematic layers:

- damage;
- shaking;
- exposure.

Possible future thematic layers:

- buildings;
- population;
- fatalities;
- injuries;
- homeless;
- economic losses;
- lifelines;
- landslides;
- liquefaction;
- critical facilities;
- road or infrastructure disruption.

Possible overlays:

- municipal boundaries;
- epicentre;
- stations;
- labels;
- administrative boundaries;
- model polygons;
- exposure assets;
- selected municipality highlight.

The layer selector should become one of the central controls of the interface.
It should not be treated as a secondary option.

The legend should be floating and collapsible, more similar to a GIS control
than to a static dashboard legend.

Future interaction:

- click polygon to display local values;
- click municipality in the impact table to zoom/highlight on map;
- switch layer without reloading the whole page;
- enable/disable overlays;
- fit map to selected scenario extent.

---

### 7.3 EventDetailsPanel

Purpose:

- summarize the selected scenario and earthquake source information;
- provide direct access to scenario outputs.

Initial fields:

- job ID;
- status;
- origin time;
- magnitude;
- depth;
- latitude;
- longitude;
- epicentral area or place;
- model ID;
- calculation time, if available.

Possible future fields:

- event ID from seismic catalogue;
- agency/source;
- event type: observed, manual, exercise, synthetic;
- calculation mode;
- ground-motion model;
- exposure model;
- fragility model;
- operator/user;
- creation time;
- last update time;
- warning/error summary.

Actions:

- download JSON;
- download GeoJSON;
- open or generate report;
- open log;
- open error file, if present.

Preferred rendering:

- compact definition-list/table style;
- labels on the left;
- values on the right;
- no large empty text blocks.

---

### 7.4 ImpactPanel

Purpose:

- summarize the main consequences of the selected scenario;
- provide the list of municipalities with expected damage.

Initial content:

- table of municipalities with damage.

Preferred table fields:

- municipality name;
- ISTAT code;
- expected damaged buildings;
- D4 + D5, if available;
- PGA or intensity, if available;
- dominant damage class, if available.

Possible future fields:

- affected population;
- expected fatalities;
- expected injuries;
- homeless;
- economic losses;
- critical facilities;
- number of unusable buildings;
- schools/hospitals affected;
- bridges or lifelines affected.

Interaction:

- click a municipality to zoom to it on the map;
- hover row to highlight municipality;
- sort by damage indicator;
- filter by minimum damage threshold.

The title may initially be “Municipalities with damage”, but the component
should be designed as a more general impact summary, so that it can later host
additional tabs or sections.

Possible future tabs:

```text
Municipalities | Buildings | Population | Casualties | Lifelines
```

---

### 7.5 FiltersPanel

Purpose:

- filter calculated scenarios in the archive.

Initial filters:

- job ID;
- magnitude range;
- status;
- model ID, if available;
- date range, if available.

Possible future filters:

- event type;
- operator/user;
- region;
- municipality;
- minimum impact level;
- completed/failed/running;
- exercise vs real event;
- model version;
- fragility model;
- exposure model.

The filters should not reload the page. They should filter the already loaded
archive where possible. Server-side filtering may be added later if the archive
becomes very large.

---

### 7.6 ArchivePanel

Purpose:

- list calculated scenarios;
- allow scenario selection;
- provide the main navigation mechanism across runs.

Initial columns:

- job ID;
- date;
- magnitude;
- depth;
- epicentral area/place;
- model ID;
- status.

Possible future columns:

- event ID;
- source agency;
- calculation time;
- number of municipalities affected;
- maximum damage level;
- operator;
- scenario type;
- creation time;
- notes;
- warning/error indicator.

Interaction:

- clicking a row selects the scenario;
- selected row is highlighted;
- selected scenario updates:
  - map;
  - event details;
  - impact panel;
  - footer status.

Future interaction:

- double click opens report/details modal;
- sortable columns;
- pagination if archive grows;
- quick search;
- right-click or action menu for downloads/logs.

---

### 7.7 Footer

Purpose:

- provide compact system feedback.

Initial content:

- server status;
- selected scenario;
- last update time, if available.

Possible future content:

- backend version;
- ShakeScenario version;
- database status;
- number of runs loaded;
- current data root;
- warning/error indicator.

The footer should remain compact and non-invasive.

---

## 8. Frontend JavaScript architecture

Recommended structure:

```text
static/js/
  api.js
  state.js
  map.js
  details.js
  impact.js
  archive.js
  filters.js
  app.js
```

Each module should have a single responsibility.

### 8.1 api.js

Responsibility:

- communicate with Flask API;
- wrap `fetch`;
- provide typed/structured helper functions;
- centralize error handling for API calls.

Possible functions:

```javascript
getStatus()
getModels()
getRuns()
getRun(jobId)
getRunLayers(jobId)
getRunLayer(jobId, layerName)
getRunMunicipalities(jobId)
```

### 8.2 state.js

Responsibility:

- maintain shared application state.

Minimum state:

```text
runs
selectedRunId
selectedRun
selectedLayer
filters
serverStatus
```

Possible future state:

```text
models
availableLayers
selectedMunicipality
mapOverlays
lastUpdated
errors
```

State changes should be explicit and predictable.

### 8.3 map.js

Responsibility:

- initialize Leaflet;
- manage base layers;
- manage thematic layers;
- manage overlays;
- render legend;
- handle map interactions.

Possible functions:

```javascript
initMap()
setScenarioLayer(layerData)
setLegend(legendData)
setSelectedLayer(layerName)
highlightMunicipality(code)
zoomToMunicipality(code)
clearScenario()
```

### 8.4 details.js

Responsibility:

- render event/scenario metadata.

Possible functions:

```javascript
renderDetails(run)
clearDetails()
setActionButtonsEnabled(enabled)
```

### 8.5 impact.js

Responsibility:

- render municipality and impact summaries.

Possible functions:

```javascript
renderMunicipalities(items)
clearImpact()
bindMunicipalitySelection(callback)
```

### 8.6 archive.js

Responsibility:

- render scenario archive table;
- handle row selection;
- visually mark selected row.

Possible functions:

```javascript
renderArchive(runs)
setSelectedRun(jobId)
bindRunSelection(callback)
```

### 8.7 filters.js

Responsibility:

- read filter form;
- apply local filters;
- reset filters.

Possible functions:

```javascript
getFilters()
applyFilters(runs, filters)
clearFilters()
bindFilterEvents(callback)
```

### 8.8 app.js

Responsibility:

- initialize the application;
- coordinate modules;
- implement application workflow.

Startup flow:

```text
init application
  -> initialize map
  -> load server status
  -> load runs
  -> render archive
  -> select most recent completed run, if available
```

Scenario selection flow:

```text
user selects run
  -> update state
  -> load run details
  -> load default layer
  -> load municipalities
  -> update map
  -> update event details
  -> update impact panel
  -> update footer
```

Layer change flow:

```text
user changes layer
  -> update selectedLayer
  -> load layer data
  -> update map
  -> update legend
```

---

## 9. Backend API draft

Initial API endpoints:

```text
GET /api/v1/status
GET /api/v1/models
GET /api/v1/runs
GET /api/v1/runs/<job_id>
GET /api/v1/runs/<job_id>/layers
GET /api/v1/runs/<job_id>/layer/<layer_name>
GET /api/v1/runs/<job_id>/municipalities
```

Possible future endpoints:

```text
POST /api/v1/runs
GET  /api/v1/runs/<job_id>/log
GET  /api/v1/runs/<job_id>/report
GET  /api/v1/runs/<job_id>/download/<resource>
GET  /api/v1/runs/<job_id>/error
GET  /api/v1/models/<model_id>
GET  /api/v1/models/<model_id>/polygons
```

The API should be simple, explicit, and stable.

The frontend should not need to know the internal filesystem layout in detail.
That knowledge belongs mostly to the backend.

---

## 10. Expected run files

A calculated scenario currently lives in:

```text
runs/job_000001/
```

Expected files:

```text
meta.json
impact_result.json
request_resolved.json
result_manifest.json
error.txt
```

### 10.1 meta.json

Expected role:

- identify the run;
- describe status;
- describe event metadata;
- describe model references;
- provide creation/calculation timestamps.

Possible fields:

```json
{
  "job_id": "job_000001",
  "status": "completed",
  "created_at": "2026-07-02T10:00:00Z",
  "completed_at": "2026-07-02T10:00:10Z",
  "event": {
    "origin_time": "2026-07-02T09:58:00Z",
    "magnitude": 5.2,
    "depth_km": 8.0,
    "latitude": 46.2,
    "longitude": 13.1,
    "place": "Friuli Venezia Giulia"
  },
  "model_id": "fvg_2024"
}
```

This is only an indicative contract and can be adjusted to the actual
ShakeScenario output.

### 10.2 impact_result.json

Expected role:

- store damage and impact results;
- provide municipality-level summaries;
- provide aggregate statistics.

Possible contents:

```text
summary
municipalities
assets/buildings
population impact
damage classes
```

### 10.3 request_resolved.json

Expected role:

- store the full resolved input request used for the run;
- make the calculation reproducible;
- preserve all defaults expanded by the backend.

### 10.4 result_manifest.json

Expected role:

- list available result resources;
- describe available layers;
- describe downloadable files.

Possible contents:

```json
{
  "layers": [
    {
      "name": "damage",
      "type": "geojson",
      "path": "damage.geojson"
    },
    {
      "name": "shaking",
      "type": "geojson",
      "path": "shaking.geojson"
    }
  ],
  "downloads": [
    {
      "name": "impact_result",
      "type": "json",
      "path": "impact_result.json"
    }
  ]
}
```

### 10.5 error.txt

Expected role:

- store run failure information;
- remain absent or empty for successful jobs.

For failed jobs, the archive should display a failed status and provide a way
to inspect the error.

---

## 11. Application state

Minimum frontend state:

```text
runs
selectedRunId
selectedRun
selectedLayer
filters
serverStatus
```

Suggested startup state transition:

```text
application starts
  -> load server status
  -> load runs
  -> render archive
  -> select most recent completed run, if available
```

Suggested scenario selection transition:

```text
user selects run
  -> update selectedRunId
  -> load run details
  -> load default layer
  -> load municipalities
  -> update map
  -> update event details
  -> update impact panel
  -> update footer
```

Suggested layer transition:

```text
user changes layer
  -> update selectedLayer
  -> load layer data
  -> update map
  -> update legend
```

Suggested refresh transition:

```text
user clicks refresh
  -> reload server status
  -> reload runs
  -> preserve selectedRunId if still available
  -> otherwise select most recent completed run
```

---

## 12. Styling principles

General style:

- neutral technical appearance;
- clean white panels;
- light borders;
- limited shadows;
- compact spacing;
- clear typography;
- no heavy decorative elements.

Color use:

- neutral background for the page;
- white panels;
- blue or similar accent for primary actions;
- semantic colors for status:
  - completed;
  - running;
  - failed;
  - warning.

Typography:

- system UI font stack;
- compact but readable;
- clear difference between labels and values.

Tables:

- compact row height;
- visible hover state;
- selected row state;
- no heavy gridlines;
- sticky header may be useful later.

Map:

- large visual priority;
- floating collapsible legend;
- layer selector always visible;
- overlays should not obscure important map content.

---

## 13. Interaction principles

- Selecting a scenario should update all dependent components.
- Changing map layer should not change the selected scenario.
- Filters should not destroy current data; they only change archive visibility.
- Refresh should preserve current selection if possible.
- Failed scenarios should remain visible and inspectable.
- Errors should be explicit but not disruptive.
- The interface should be usable during an operational event without requiring
  page navigation.

---

## 14. Error handling

Expected error cases:

- backend not reachable;
- no `SHAKESCENARIO_HOME` configured;
- missing `runs/` directory;
- empty archive;
- malformed `meta.json`;
- missing `impact_result.json`;
- missing layer listed in manifest;
- failed scenario with `error.txt`.

Frontend behavior:

- show clear empty/error states;
- do not crash the whole interface;
- preserve the ability to view the archive if one run is malformed;
- display failed runs in the archive;
- provide access to error text when available.

Backend behavior:

- return structured JSON errors;
- avoid exposing unnecessary internal tracebacks to the user interface;
- log full details in server logs.

---

## 15. Design decisions already agreed

- Use the visual style of the internal WebUI draft.
- Use the layout logic of the `draft_webpage` concept.
- Do not use the PCR dashboard as visual reference.
- Use the PCR dashboard only as a functional hint for key operational
  information, especially municipalities with damage.
- Treat the application as a GIS-like technical working tool.
- Keep dependencies minimal.
- Keep frontend JavaScript modular.
- Use `SHAKESCENARIO_HOME` to switch between local testing and production.
- Use Apache + Gunicorn + Flask for deployment.
- Prefer filesystem-based run discovery at this stage.
- Avoid a heavy database-driven design until there is a clear need.

---

## 16. Implementation roadmap

### Phase 1: Static layout

- finalize `templates/index.html`;
- complete base `static/css/style.css`;
- verify layout visually with placeholder content;
- ensure the interface scales reasonably on a desktop monitor.

### Phase 2: Frontend module skeleton

- create JavaScript module files;
- initialize app state;
- initialize Leaflet map;
- render static placeholder archive;
- test component update functions.

### Phase 3: Backend API alignment

- configure `SHAKESCENARIO_HOME`;
- implement or update `/api/v1/status`;
- implement or update `/api/v1/runs`;
- implement run detail endpoint;
- implement layer and municipality endpoints.

### Phase 4: Scenario selection

- load scenario archive from backend;
- select most recent completed scenario;
- update event details;
- update map layer;
- update municipality table;
- update footer.

### Phase 5: Layer management

- load available layers from manifest;
- implement layer selector;
- render damage layer;
- render shaking layer;
- render exposure/polygon layer;
- implement legend update.

### Phase 6: Operational functions

- implement downloads;
- implement log/error viewing;
- implement report access;
- implement new scenario form/workflow;
- implement job status refresh.

### Phase 7: Refinement

- improve responsive behavior;
- improve table usability;
- add sorting/filtering;
- add municipality highlight;
- add map overlay controls;
- refine status and error messages.

---

## 17. Open questions

These points are intentionally left open and should be resolved during
implementation:

- exact schema of `meta.json`;
- exact schema of `impact_result.json`;
- exact structure of `result_manifest.json`;
- whether municipality polygons come from model `polygons.json` or from run
  outputs;
- how to represent multiple fragility models in the UI;
- whether the archive should eventually read from SQLite instead of only from
  the filesystem;
- how to expose logs safely;
- how to implement user authentication, if required;
- how to distinguish real events, exercises, and synthetic scenarios;
- how to handle multiple models and model versions;
- how to structure report generation.

---

## 18. Immediate next steps

Current working order:

1. save this blueprint in `WebUI/docs/ui_blueprint.md`;
2. complete `static/css/style.css`;
3. create JavaScript skeleton modules;
4. update `index.html` only if needed by the component design;
5. align Flask API with the current data directory structure;
6. test with one or more sample `runs/job_*` directories.

