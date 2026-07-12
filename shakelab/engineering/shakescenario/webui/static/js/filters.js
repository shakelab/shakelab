"use strict";

/**
 * Read filters from the form.
 */
function getFilters() {
    const jobId = document.getElementById("filter-job-id");
    const magMin = document.getElementById("filter-mag-min");
    const magMax = document.getElementById("filter-mag-max");
    const status = document.getElementById("filter-status");

    return {
        jobId: jobId ? jobId.value.trim() : "",
        magnitudeMin: magMin && magMin.value !== ""
            ? Number(magMin.value)
            : null,
        magnitudeMax: magMax && magMax.value !== ""
            ? Number(magMax.value)
            : null,
        status: status ? status.value : ""
    };
}

/**
 * Apply filters locally to the loaded runs.
 */
function applyFilters(runs, filters) {
    if (!Array.isArray(runs)) {
        return [];
    }

    return runs.filter((run) => {
        const event = run.event || {};
        const magnitude = Number(event.magnitude);

        if (filters.jobId) {
            const jobId = String(run.job_id || "").toLowerCase();
            const query = filters.jobId.toLowerCase();

            if (!jobId.includes(query)) {
                return false;
            }
        }

        if (filters.status) {
            const status = String(run.status || "").toLowerCase();

            if (status !== filters.status.toLowerCase()) {
                return false;
            }
        }

        if (
            filters.magnitudeMin !== null &&
            !Number.isNaN(filters.magnitudeMin)
        ) {
            if (Number.isNaN(magnitude) ||
                magnitude < filters.magnitudeMin) {
                return false;
            }
        }

        if (
            filters.magnitudeMax !== null &&
            !Number.isNaN(filters.magnitudeMax)
        ) {
            if (Number.isNaN(magnitude) ||
                magnitude > filters.magnitudeMax) {
                return false;
            }
        }

        return true;
    });
}

/**
 * Clear filter controls.
 */
function clearFilters() {
    const jobId = document.getElementById("filter-job-id");
    const magMin = document.getElementById("filter-mag-min");
    const magMax = document.getElementById("filter-mag-max");
    const status = document.getElementById("filter-status");

    if (jobId) {
        jobId.value = "";
    }

    if (magMin) {
        magMin.value = "2";
    }

    if (magMax) {
        magMax.value = "10";
    }

    if (status) {
        status.value = "";
    }
}

/**
 * Bind filter form events.
 */
function bindFilterEvents(callback) {
    const form = document.getElementById("scenario-filters");
    const clearButton = document.getElementById("btn-clear-filters");

    if (form) {
        form.addEventListener("submit", (event) => {
            event.preventDefault();

            if (callback) {
                callback(getFilters());
            }
        });
    }

    if (clearButton) {
        clearButton.addEventListener("click", () => {
            clearFilters();

            if (callback) {
                callback(getFilters());
            }
        });
    }
}