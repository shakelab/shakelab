"use strict";

function formatUtcNow() {
    return new Date().toISOString().replace(/\.\d{3}Z$/, "Z");
}

function getScenarioModelIds(payload) {
    if (!payload) {
        return [];
    }

    if (Array.isArray(payload)) {
        return payload;
    }

    if (Array.isArray(payload.model_ids)) {
        return payload.model_ids;
    }

    if (
        payload.models &&
        Array.isArray(payload.models.model_ids)
    ) {
        return payload.models.model_ids;
    }

    return [];
}

function setScenarioMessage(text, cssClass = "") {
    const message = document.getElementById("new-scenario-message");

    if (!message) {
        return;
    }

    message.textContent = text || "";
    message.className = `form-message ${cssClass}`.trim();
}

async function populateScenarioModels() {
    const selector = document.getElementById("scenario-model");

    if (!selector) {
        return;
    }

    selector.innerHTML = `
        <option value="">Loading models...</option>
    `;

    try {
        const payload = await getModels();
        const modelIds = getScenarioModelIds(payload);

        if (modelIds.length === 0) {
            selector.innerHTML = `
                <option value="">No models available</option>
            `;
            return;
        }

        selector.innerHTML = modelIds.map((modelId) => {
            return `
                <option value="${modelId}">
                    ${modelId}
                </option>
            `;
        }).join("");

    } catch (error) {
        selector.innerHTML = `
            <option value="">Unable to load models</option>
        `;

        setScenarioMessage(error.message, "error");
    }
}

function resetScenarioForm() {
    const originTime = document.getElementById("scenario-origin-time");

    if (originTime && !originTime.value) {
        originTime.value = formatUtcNow();
    }

    setScenarioMessage("");
}

async function openNewScenarioDialog() {
    const modal = document.getElementById("new-scenario-modal");

    if (!modal) {
        return;
    }

    resetScenarioForm();
    modal.classList.remove("hidden");

    await populateScenarioModels();
}

function closeNewScenarioDialog() {
    const modal = document.getElementById("new-scenario-modal");

    if (!modal) {
        return;
    }

    modal.classList.add("hidden");
    setScenarioMessage("");
}

function readScenarioForm() {
    return {
        tag: document.getElementById("scenario-tag").value.trim(),
        model_id: document.getElementById("scenario-model").value,
        origin_time: document
            .getElementById("scenario-origin-time")
            .value.trim(),
        magnitude: Number(document.getElementById(
            "scenario-magnitude"
        ).value),
        longitude: Number(document.getElementById(
            "scenario-longitude"
        ).value),
        latitude: Number(document.getElementById(
            "scenario-latitude"
        ).value),
        depth: Number(document.getElementById("scenario-depth").value)
    };
}

function validateScenarioPayload(payload) {
    const required = [
        "tag",
        "model_id",
        "origin_time",
        "magnitude",
        "longitude",
        "latitude",
        "depth"
    ];

    for (const key of required) {
        if (
            payload[key] === null ||
            payload[key] === undefined ||
            payload[key] === "" ||
            Number.isNaN(payload[key])
        ) {
            return `Missing or invalid field: ${key}`;
        }
    }

    return null;
}

async function handleScenarioSubmit(event) {
    event.preventDefault();

    const payload = readScenarioForm();
    const error = validateScenarioPayload(payload);

    if (error) {
        setScenarioMessage(error, "error");
        return;
    }

    setScenarioMessage("Submitting scenario...");

    try {
        const response = await submitScenario(payload);
        const jobId = response.job_id || "unknown";

        setScenarioMessage(
            `Scenario submitted successfully: ${jobId}`,
            "success"
        );

        await loadRuns();

    } catch (submitError) {
        setScenarioMessage(submitError.message, "error");
    }
}

function bindNewScenarioDialog() {
    const newButton = document.getElementById("btn-new-run");
    const closeButton = document.getElementById("btn-close-new-scenario");
    const cancelButton = document.getElementById("btn-cancel-new-scenario");
    const form = document.getElementById("new-scenario-form");
    const modal = document.getElementById("new-scenario-modal");

    if (newButton) {
        newButton.addEventListener("click", openNewScenarioDialog);
    }

    if (closeButton) {
        closeButton.addEventListener("click", closeNewScenarioDialog);
    }

    if (cancelButton) {
        cancelButton.addEventListener("click", closeNewScenarioDialog);
    }

    if (form) {
        form.addEventListener("submit", handleScenarioSubmit);
    }

    if (modal) {
        modal.addEventListener("click", (event) => {
            if (event.target === modal) {
                closeNewScenarioDialog();
            }
        });
    }
}