"use strict";

/**
 * Update footer information.
 */
function updateFooter() {
    const serverStatus = document.getElementById("server-status");
    const selectedScenario = document.getElementById("selected-scenario");

    if (serverStatus) {
        const status = AppState.serverStatus &&
            AppState.serverStatus.status
            ? AppState.serverStatus.status
            : "unknown";

        serverStatus.textContent = `Server status: ${status}`;
    }

    if (selectedScenario) {
        selectedScenario.textContent = AppState.selectedRunId
            ? `Selected scenario: ${AppState.selectedRunId}`
            : "No selected scenario";
    }
}

/**
 * Load server status.
 */
async function loadStatus() {
    try {
        const status = await getStatus();
        setServerStatus(status);
    } catch (error) {
        console.error("Unable to load server status:", error);
        setServerStatus({
            status: "offline",
            error: error.message
        });
    }

    updateFooter();
}

/**
 * Return true for completed runs.
 */
function isCompletedRun(run) {
    return String(run.status || "").toLowerCase() === "completed";
}

/**
 * Select default run after archive loading.
 */
function getDefaultRun(runs) {
    if (!Array.isArray(runs) || runs.length === 0) {
        return null;
    }

    if (!CONFIG.autoSelectLatest) {
        return null;
    }

    const completed = runs.filter(isCompletedRun);

    if (completed.length > 0) {
        return completed[0];
    }

    return runs[0];
}

/**
 * Load and render archive.
 */
async function loadRuns() {
    try {
        const payload = await getRuns();
        const runs = Array.isArray(payload) ? payload : payload.runs;

        setRuns(runs || []);

        const filtered = applyFilters(AppState.runs, AppState.filters);
        renderArchive(filtered);

        const defaultRun = getDefaultRun(filtered);

        if (defaultRun) {
            await selectRun(defaultRun.job_id);
        }
    } catch (error) {
        console.error("Unable to load runs:", error);
        setRuns([]);
        renderArchive([]);
        clearDetails();
        clearImpact();
    }
}

/**
 * Load available layers for a run.
 */
async function loadLayers(jobId) {
    try {
        const payload = await getRunLayers(jobId);
        const layers = Array.isArray(payload) ? payload : payload.layers;

        setAvailableLayers(layers || []);
    } catch (error) {
        console.warn("Unable to load run layers:", error);
        setAvailableLayers([]);
        renderLayerSelector();
    }
}

/**
 * Render available map layers in the layer selector.
 */
function renderLayerSelector() {
    const selector = document.getElementById("layer-selector");

    if (!selector) {
        return;
    }

    const layers = AppState.availableLayers;

    if (!Array.isArray(layers) || layers.length === 0) {
        selector.innerHTML = `
            <option value="damage">Damage</option>
        `;
        setSelectedLayer("damage");
        return;
    }

    selector.innerHTML = layers.map((layer) => {
        const name = layer.name;
        const title = layer.title || layer.name;

        return `
            <option value="${name}">
                ${title}
            </option>
        `;
    }).join("");

    const selectedExists = layers.some((layer) => {
        return layer.name === AppState.selectedLayer;
    });

    if (!selectedExists) {
        setSelectedLayer(layers[0].name);
    }

    selector.value = AppState.selectedLayer;
}

/**
 * Load and render selected layer.
 */
async function loadSelectedLayer(jobId) {
    try {
        const layerName = AppState.selectedLayer || CONFIG.defaultLayer;
        const layerData = await getRunLayer(jobId, layerName);

        setScenarioLayer(layerData);
        setDefaultDamageLegend();
    } catch (error) {
        console.warn("Unable to load selected layer:", error);
        clearScenarioLayer();
    }
}

/**
 * Load and render asset-level impact summary.
 */
async function loadImpactAssets(jobId) {
    try {
        const payload = await getRunImpactAssets(jobId);
        const items = Array.isArray(payload)
            ? payload
            : payload.assets;

        setImpactAssets(items || []);
        renderImpactAssets(AppState.impactAssets);
    } catch (error) {
        console.warn("Unable to load impact assets:", error);
        setImpactAssets([]);
        clearImpact();
    }
}

/**
 * Select a scenario run.
 */
async function selectRun(jobId) {
    if (!jobId) {
        return;
    }

    try {
        const run = await getRun(jobId);

        setSelectedRun(run);
        highlightSelectedRun(jobId);
        renderDetails(run);
        setEpicenter(run.event);

        const subtitle = document.getElementById("map-subtitle");
        if (subtitle) {
            subtitle.textContent = `Scenario ${jobId}`;
        }

        await loadLayers(jobId);
        renderLayerSelector();
        await loadSelectedLayer(jobId);
        await loadImpactAssets(jobId);

        updateFooter();
    } catch (error) {
        console.error(`Unable to select run ${jobId}:`, error);
        clearDetails();
        clearImpact();
        clearScenarioLayer();
    }
}

/**
 * Apply filters and re-render archive.
 */
function handleFilters(filters) {
    setFilters(filters);

    const filtered = applyFilters(AppState.runs, AppState.filters);
    renderArchive(filtered);

    if (
        AppState.selectedRunId &&
        !filtered.some((run) => run.job_id === AppState.selectedRunId)
    ) {
        setSelectedRun(null);
        highlightSelectedRun(null);
        clearDetails();
        clearImpact();
        clearScenarioLayer();
        updateFooter();
    } else {
        highlightSelectedRun(AppState.selectedRunId);
    }
}

/**
 * Handle layer selector changes.
 */
async function handleLayerChange(event) {
    const layerName = event.target.value;

    setSelectedLayer(layerName);

    if (AppState.selectedRunId) {
        await loadSelectedLayer(AppState.selectedRunId);
    }
}

async function handleRunDeletion(jobId) {
    const confirmed = window.confirm(
        `Delete scenario ${jobId}?\n\n` +
        "This operation cannot be undone."
    );

    if (!confirmed) {
        return;
    }

    try {
        await deleteRun(jobId);

        if (AppState.selectedRunId === jobId) {
            setSelectedRun(null);
            clearDetails();
            clearImpact();
            clearScenarioLayer();

            const subtitle = document.getElementById("map-subtitle");
            if (subtitle) {
                subtitle.textContent = "No scenario selected";
            }
        }

        await loadRuns();

    } catch (error) {
        window.alert(`Unable to delete scenario ${jobId}:\n${error.message}`);
    }
}

/**
 * Bind global UI events.
 */
function bindAppEvents() {
    const refreshButton = document.getElementById("btn-refresh");
    const layerSelector = document.getElementById("layer-selector");

    if (refreshButton) {
        refreshButton.addEventListener("click", async () => {
            await loadStatus();
            await loadRuns();
        });
    }

    if (layerSelector) {
        layerSelector.addEventListener("change", handleLayerChange);
    }

    bindRunSelection(selectRun);
    bindRunDeletion(handleRunDeletion);
    bindFilterEvents(handleFilters);

    bindImpactAssetSelection((assetId) => {
        console.log("Impact asset selected:", assetId);
        /* Future: highlightImpactAsset(assetId); */
    });

    bindMapControls();

    bindArchiveSorting(() => {
        const filtered = applyFilters(AppState.runs, AppState.filters);
        renderArchive(filtered);
    });

    bindDownloadButtons();
    bindNewScenarioDialog();
}

/**
 * Initialize application.
 */
async function initApp() {
    initMap();
    setDefaultDamageLegend();

    bindAppEvents();

    clearDetails();
    clearImpact();
    updateFooter();

    await loadStatus();
    await loadRuns();
}

document.addEventListener("DOMContentLoaded", () => {
    initApp();
});