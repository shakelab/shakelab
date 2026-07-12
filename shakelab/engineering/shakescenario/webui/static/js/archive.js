"use strict";

/**
 * Render the archive table.
 */
const ArchiveSortState = {
    key: "origin_time",
    direction: "desc"
};

function getRunSortValue(run, key) {
    const event = run.event || {};

    if (key === "origin_time") {
        return event.origin_time || "";
    }

    if (key === "magnitude") {
        return Number(event.magnitude);
    }

    if (key === "longitude") {
        return Number(event.longitude);
    }

    if (key === "latitude") {
        return Number(event.latitude);
    }

    if (key === "depth_km") {
        return Number(event.depth_km);
    }

    if (key === "model_id") {
        return run.model_id || "";
    }

    if (key === "status") {
        return run.status || "";
    }

    return run.job_id || "";
}

function compareRuns(a, b, key) {
    const valueA = getRunSortValue(a, key);
    const valueB = getRunSortValue(b, key);

    const numA = Number(valueA);
    const numB = Number(valueB);

    if (!Number.isNaN(numA) && !Number.isNaN(numB)) {
        return numA - numB;
    }

    return String(valueA).localeCompare(String(valueB));
}

function sortRuns(runs) {
    const sorted = [...runs];

    sorted.sort((a, b) => {
        const result = compareRuns(a, b, ArchiveSortState.key);

        if (ArchiveSortState.direction === "desc") {
            return -result;
        }

        return result;
    });

    return sorted;
}

function formatTableNumber(value, decimals = 3) {
    const number = Number(value);

    if (Number.isNaN(number)) {
        return "-";
    }

    return number.toFixed(decimals);
}

function renderArchive(runs) {
    const body = document.getElementById("scenario-table-body");
    const counter = document.getElementById("archive-count");

    if (!body) {
        return;
    }

    if (!Array.isArray(runs) || runs.length === 0) {
        body.innerHTML = `
            <tr>
                <td colspan="9" class="empty-state">
                    No scenarios found.
                </td>
            </tr>
        `;

        if (counter) {
            counter.textContent = "0 scenarios";
        }

        return;
    }

    const sortedRuns = sortRuns(runs);

    if (counter) {
        counter.textContent = `${sortedRuns.length} scenario(s)`;
    }

    const rows = sortedRuns.map((run) => {
        const event = run.event || {};

        return `
            <tr data-job="${run.job_id}">
                <td>${safeValue(run.job_id)}</td>
                <td>${safeValue(event.origin_time)}</td>
                <td>${safeValue(event.magnitude)}</td>
                <td>${formatTableNumber(event.longitude, 4)}</td>
                <td>${formatTableNumber(event.latitude, 4)}</td>
                <td>${formatTableNumber(event.depth_km, 1)}</td>
                <td>${safeValue(run.model_id)}</td>
                <td>${renderStatus(run.status)}</td>
                <td class="scenario-actions">
                    <button
                        type="button"
                        class="delete-run-button"
                        data-delete-job="${run.job_id}"
                        title="Delete scenario">
                        🗑
                    </button>
                </td>
            </tr>
        `;
    }).join("");

    body.innerHTML = rows;

    updateArchiveSortHeader();
}

/**
 * Render status pill.
 */
function renderStatus(status) {

    const value = (status || "").toLowerCase();

    let css = "status-running";

    if (value === "completed") {
        css = "status-completed";
    }

    if (value === "failed") {
        css = "status-failed";
    }

    return `
        <span class="status-pill ${css}">
            ${safeValue(status)}
        </span>
    `;
}

/**
 * Highlight selected row.
 */
function highlightSelectedRun(jobId) {

    document
        .querySelectorAll("#scenario-table-body tr")
        .forEach((row) => {

            row.classList.remove("selected");

            if (
                row.dataset.job === jobId
            ) {
                row.classList.add("selected");
            }

        });
}

/**
 * Install click callback.
 */
function bindRunSelection(callback) {

    const body = document.getElementById(
        "scenario-table-body"
    );

    if (!body) {
        return;
    }

    body.addEventListener("click", (event) => {

        if (event.target.closest("[data-delete-job]")) {
            return;
        }

        const row = event.target.closest(
            "tr[data-job]"
        );

        if (!row) {
            return;
        }

        const jobId = row.dataset.job;

        if (callback) {
            callback(jobId);
        }

    });

}

function updateArchiveSortHeader() {
    document
        .querySelectorAll("#scenario-table th[data-sort-key]")
        .forEach((header) => {
            const key = header.dataset.sortKey;
            const baseLabel = header.textContent
                .replace(" ▲", "")
                .replace(" ▼", "");

            if (key === ArchiveSortState.key) {
                const arrow = ArchiveSortState.direction === "asc"
                    ? " ▲"
                    : " ▼";

                header.textContent = baseLabel + arrow;
            } else {
                header.textContent = baseLabel;
            }
        });
}

function bindArchiveSorting(callback) {
    document
        .querySelectorAll("#scenario-table th[data-sort-key]")
        .forEach((header) => {
            header.addEventListener("click", () => {
                const key = header.dataset.sortKey;

                if (ArchiveSortState.key === key) {
                    ArchiveSortState.direction =
                        ArchiveSortState.direction === "asc"
                            ? "desc"
                            : "asc";
                } else {
                    ArchiveSortState.key = key;
                    ArchiveSortState.direction = "asc";
                }

                if (callback) {
                    callback();
                }
            });
        });
}

function bindRunDeletion(callback) {
    const body = document.getElementById("scenario-table-body");

    if (!body) {
        return;
    }

    body.addEventListener("click", async (event) => {
        const button = event.target.closest("[data-delete-job]");

        if (!button) {
            return;
        }

        const jobId = button.dataset.deleteJob;

        if (!jobId) {
            return;
        }

        if (callback) {
            await callback(jobId);
        }
    });
}