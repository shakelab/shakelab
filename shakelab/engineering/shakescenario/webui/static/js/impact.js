"use strict";

/**
 * Clear the impact panel.
 */
function clearImpact() {
    const panel = document.getElementById("impact-assets-summary");

    if (!panel) {
        return;
    }

    panel.innerHTML = `
        <div class="empty-state">
            No impact results available.
        </div>
    `;
}

function formatNumber(value, decimals = 1) {
    const number = Number(value);

    if (Number.isNaN(number)) {
        return "-";
    }

    return number.toFixed(decimals);
}

/**
 * Return the first available property from an object.
 */
function getProperty(object, keys, fallback = "-") {
    if (!object) {
        return fallback;
    }

    for (const key of keys) {
        if (object[key] !== undefined && object[key] !== null) {
            return object[key];
        }
    }

    return fallback;
}

/**
 * Render asset-level impact summary.
 */
function renderImpactAssets(items) {
    const panel = document.getElementById("impact-assets-summary");

    if (!panel) {
        return;
    }

    if (!Array.isArray(items) || items.length === 0) {
        clearImpact();
        return;
    }

    const rows = items.map((item) => {
        const assetId = getProperty(item, [
            "asset_id",
            "id"
        ]);
        
        const name = getProperty(item, [
            "name"
        ]);

        const totalBuildings = getProperty(item, [
            "n_units",
            "total_units",
            "count"
        ]);

        const severe = getProperty(item, [
            "damage_d4_d5",
            "d4_d5",
            "D4_D5"
        ]);

        const pga = getProperty(item, [
            "pga",
            "PGA"
        ]);

        return `
            <tr data-asset-id="${assetId}">
                <td>${assetId}</td>
                <td>${name}</td>
                <td>${formatNumber(totalBuildings, 1)}</td>
                <td>${formatNumber(severe, 1)}</td>
                <td>${formatNumber(pga, 3)}</td>
            </tr>
        `;
    }).join("");

    panel.innerHTML = `
        <div class="table-wrapper">
            <table class="impact-assets-table">
                <thead>
                    <tr>
                        <th>Asset ID</th>
                        <th>Name</th>
                        <th>Buildings</th>
                        <th>D4+D5</th>
                        <th>PGA (g)</th>
                    </tr>
                </thead>
                <tbody>
                    ${rows}
                </tbody>
            </table>
        </div>
    `;
}

/**
 * Bind click events on impact asset rows.
 */
function bindImpactAssetSelection(callback) {
    const panel = document.getElementById("impact-assets-summary");

    if (!panel) {
        return;
    }

    panel.addEventListener("click", (event) => {
        const row = event.target.closest("tr[data-asset-id]");

        if (!row) {
            return;
        }

        const assetId = row.getAttribute("data-asset-id");
        
        if (callback) {
            callback(assetId);
        }
    });
}