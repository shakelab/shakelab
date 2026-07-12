"use strict";

const API_BASE = "/api/v1";

async function requestJson(url, options = {}) {
    const response = await fetch(url, {
        headers: {
            "Accept": "application/json",
            ...(options.headers || {})
        },
        ...options
    });

    let payload = null;

    try {
        payload = await response.json();
    } catch (error) {
        payload = null;
    }

    if (!response.ok) {
        const message = payload && payload.message
            ? payload.message
            : payload && payload.error
                ? payload.error
                : `HTTP error ${response.status}`;

        throw new Error(message);
    }

    return payload;
}

async function getStatus() {
    return requestJson(`${API_BASE}/status`);
}

async function getServer() {
    return requestJson(`${API_BASE}/server`);
}

async function getModels() {
    return requestJson(`${API_BASE}/models`);
}

async function submitScenario(data) {
    return requestJson(`${API_BASE}/scenarios`, {
        method: "POST",
        headers: {
            "Content-Type": "application/json"
        },
        body: JSON.stringify(data)
    });
}

async function getRuns() {
    return requestJson(`${API_BASE}/runs`);
}

async function getRun(jobId) {
    return requestJson(`${API_BASE}/runs/${encodeURIComponent(jobId)}`);
}

async function deleteRun(jobId) {
    return requestJson(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}`,
        {
            method: "DELETE"
        }
    );
}

async function getRunLayers(jobId) {
    return requestJson(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}/layers`
    );
}

async function getRunLayer(jobId, layerName) {
    return requestJson(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}` +
        `/layer/${encodeURIComponent(layerName)}`
    );
}

async function getRunImpactAssets(jobId) {
    return requestJson(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}/impact-assets`
    );
}

async function requestBlob(url) {
    const response = await fetch(url);

    if (!response.ok) {
        throw new Error(`HTTP error ${response.status}`);
    }

    return response.blob();
}

async function getRunJsonBlob(jobId) {
    return requestBlob(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}/download/json`
    );
}

async function getRunGeoJsonBlob(jobId) {
    return requestBlob(
        `${API_BASE}/runs/${encodeURIComponent(jobId)}/download/geojson`
    );
}