"use strict";

/**
 * Return a safe string representation.
 */
function safeValue(value, fallback = "-") {
    if (value === null || value === undefined || value === "") {
        return fallback;
    }

    return String(value);
}

/**
 * Format a floating-point coordinate.
 */
function formatCoordinate(value) {
    if (value === null || value === undefined) {
        return "-";
    }

    return Number(value).toFixed(4);
}

/**
 * Format depth.
 */
function formatDepth(value) {
    if (value === null || value === undefined) {
        return "-";
    }

    return `${Number(value).toFixed(1)} km`;
}

/**
 * Render one row.
 */
function detailRow(label, value) {
    return `
        <div class="detail-item">
            <div class="detail-label">${label}</div>
            <div class="detail-value">${safeValue(value)}</div>
        </div>
    `;
}

/**
 * Clear panel.
 */
function clearDetails() {

    const panel = document.getElementById("event-details");

    panel.innerHTML = `
        <div class="empty-state">
            No scenario selected.
        </div>
    `;

    setActionButtonsEnabled(false);
}

/**
 * Enable or disable action buttons.
 */
function setActionButtonsEnabled(enabled) {

    [
        "btn-download-json",
        "btn-download-geojson",
        "btn-report"
    ].forEach((id) => {

        const button = document.getElementById(id);

        if (button) {
            button.disabled = !enabled;
        }

    });
}

/**
 * Render the selected run.
 */
function renderDetails(run) {

    if (!run) {
        clearDetails();
        return;
    }

    const event = run.event || {};

    const html = [
        detailRow("Scenario ID", run.job_id),
        detailRow("Model", run.model_id),
        detailRow("Date and Time", event.origin_time),
        detailRow("Magnitude", event.magnitude),
        detailRow("Longitude", formatCoordinate(event.longitude)),
        detailRow("Latitude", formatCoordinate(event.latitude)),
        detailRow("Depth", formatDepth(event.depth_km)),
        detailRow("Status", run.status)
    ].join("");

    document.getElementById("event-details").innerHTML = html;

    setActionButtonsEnabled(true);
}

function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");

    link.href = url;
    link.download = filename;

    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    URL.revokeObjectURL(url);
}

async function downloadSelectedJson() {
    if (!AppState.selectedRunId) {
        return;
    }

    const blob = await getRunJsonBlob(AppState.selectedRunId);
    downloadBlob(blob, `${AppState.selectedRunId}_impact_assets.json`);
}

async function downloadSelectedGeoJson() {
    if (!AppState.selectedRunId) {
        return;
    }

    const blob = await getRunGeoJsonBlob(AppState.selectedRunId);
    downloadBlob(blob, `${AppState.selectedRunId}_damage.geojson`);
}

function bindDownloadButtons() {
    const jsonButton = document.getElementById("btn-download-json");
    const geojsonButton = document.getElementById("btn-download-geojson");
    const reportButton = document.getElementById("btn-report");

    if (jsonButton) {
        jsonButton.addEventListener("click", downloadSelectedJson);
    }

    if (geojsonButton) {
        geojsonButton.addEventListener("click", downloadSelectedGeoJson);
    }

    if (reportButton) {
        reportButton.disabled = true;
    }
}