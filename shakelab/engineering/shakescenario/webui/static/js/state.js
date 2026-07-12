"use strict";

const AppState = {
    runs: [],
    selectedRunId: null,
    selectedRun: null,
    selectedLayer: "damage",
    filters: {
        jobId: "",
        magnitudeMin: null,
        magnitudeMax: null,
        status: ""
    },
    serverStatus: null,
    impactAssets: [],
    availableLayers: []
};

function setRuns(runs) {
    AppState.runs = Array.isArray(runs) ? runs : [];
}

function setSelectedRun(run) {
    AppState.selectedRun = run || null;
    AppState.selectedRunId = run && run.job_id ? run.job_id : null;
}

function setSelectedLayer(layerName) {
    AppState.selectedLayer = layerName || "damage";
}

function setFilters(filters) {
    AppState.filters = {
        jobId: filters.jobId || "",
        magnitudeMin: filters.magnitudeMin,
        magnitudeMax: filters.magnitudeMax,
        status: filters.status || ""
    };
}

function setServerStatus(status) {
    AppState.serverStatus = status || null;
}

function setImpactAssets(items) {
    AppState.impactAssets = Array.isArray(items) ? items : [];
}

function setAvailableLayers(layers) {
    AppState.availableLayers = Array.isArray(layers) ? layers : [];
}