"use strict";

let map = null;
let scenarioLayer = null;
let assetHighlightLayer = null;
let epicenterLayer = null;

const DEFAULT_MAP_CENTER = [46.15, 13.10];
const DEFAULT_MAP_ZOOM = 8;

function initMap() {
    map = L.map("map", {
        zoomControl: true
    }).setView(DEFAULT_MAP_CENTER, DEFAULT_MAP_ZOOM);

    L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
        maxZoom: 19,
        attribution: "&copy; OpenStreetMap contributors"
    }).addTo(map);

    setTimeout(refreshMapSize, 100);
    setTimeout(refreshMapSize, 500);

    return map;
}

function clearScenarioLayer() {
    if (scenarioLayer !== null) {
        map.removeLayer(scenarioLayer);
        scenarioLayer = null;
    }
}

function clearAssetHighlight() {
    if (assetHighlightLayer !== null) {
        map.removeLayer(assetHighlightLayer);
        assetHighlightLayer = null;
    }
}

function clearEpicenter() {
    if (epicenterLayer !== null) {
        map.removeLayer(epicenterLayer);
        epicenterLayer = null;
    }
}

function setEpicenter(event) {
    clearEpicenter();

    if (!event) {
        return;
    }

    const lat = Number(event.latitude);
    const lon = Number(event.longitude);

    if (Number.isNaN(lat) || Number.isNaN(lon)) {
        return;
    }

    const icon = L.divIcon({
        className: "epicenter-marker",
        html: "★",
        iconSize: [28, 28],
        iconAnchor: [14, 14]
    });

    epicenterLayer = L.marker([lat, lon], {
        icon: icon,
        zIndexOffset: 1000
    }).addTo(map);
}

function getFeatureValue(feature, keys) {
    if (!feature || !feature.properties) {
        return null;
    }

    for (const key of keys) {
        if (feature.properties[key] !== undefined) {
            return feature.properties[key];
        }
    }

    return null;
}

function formatPopupNumber(value, decimals = 1) {
    const number = Number(value);

    if (Number.isNaN(number)) {
        return "-";
    }

    return number.toFixed(decimals);
}

function getDamageColor(value) {
    if (value === null || value === undefined || Number.isNaN(value)) {
        return "#ffffff";
    }

    if (value < 5) {
        return "#ffffff";
    }
    if (value < 20) {
        return "#fee8a8";
    }
    if (value < 50) {
        return "#fdd17a";
    }
    if (value < 100) {
        return "#fdae61";
    }
    if (value < 200) {
        return "#f46d43";
    }
    if (value < 500) {
        return "#d73027";
    }

    return "#a50026";
}

function styleScenarioFeature(feature) {
    const damage = getFeatureValue(feature, [
        "damage_d4_d5",
        "D4_D5",
        "d4_d5"
    ]);

    return {
        color: "#374151",
        weight: 0.7,
        opacity: 0.8,
        fillColor: getDamageColor(Number(damage)),
        fillOpacity: 0.65
    };
}

function bindScenarioPopup(feature, layer) {
    const props = feature.properties || {};

    const assetId = getFeatureValue(feature, [
        "asset_id",
        "id"
    ]);

    const name = getFeatureValue(feature, [
        "name",
        "asset_name"
    ]);

    const totalBuildings = getFeatureValue(feature, [
        "n_units",
        "total_units",
        "count"
    ]);

    const severe = getFeatureValue(feature, [
        "damage_d4_d5",
        "d4_d5",
        "D4_D5"
    ]);

    const pga = getFeatureValue(feature, [
        "pga",
        "PGA"
    ]);

    const html = `
        <div class="popup-title">${safeValue(name)}</div>
        <div class="popup-row">
            <span class="popup-label">Asset ID</span>
            <span>${safeValue(assetId)}</span>
        </div>
        <div class="popup-row">
            <span class="popup-label">Buildings</span>
            <span>${formatPopupNumber(totalBuildings, 1)}</span>
        </div>
        <div class="popup-row">
            <span class="popup-label">D4+D5</span>
            <span>${formatPopupNumber(severe, 1)}</span>
        </div>
        <div class="popup-row">
            <span class="popup-label">PGA (g)</span>
            <span>${formatPopupNumber(pga, 3)}</span>
        </div>
    `;

    layer.bindPopup(html);
}

function pointToLayer(feature, latlng) {
    const damage = getFeatureValue(feature, [
        "damage_d4_d5",
        "damage_total",
        "damage",
        "D4_D5",
        "d4_d5"
    ]);

    const units = getFeatureValue(feature, [
        "n_units",
        "count",
        "asset_count"
    ]);

    let radius = 5;

    if (units !== null && !Number.isNaN(Number(units))) {
        radius = Math.max(4, Math.min(14, Math.sqrt(Number(units)) * 0.35));
    }

    return L.circleMarker(latlng, {
        radius: radius,
        color: "#374151",
        weight: 0.7,
        opacity: 0.9,
        fillColor: getDamageColor(Number(damage)),
        fillOpacity: 0.75
    });
}

function setScenarioLayer(geojson) {
    clearScenarioLayer();
    clearAssetHighlight();

    if (!geojson) {
        return;
    }

    scenarioLayer = L.geoJSON(geojson, {
        style: styleScenarioFeature,
        pointToLayer: pointToLayer,
        onEachFeature: bindScenarioPopup
    }).addTo(map);

    try {
        const bounds = scenarioLayer.getBounds();
        if (bounds.isValid()) {
            map.fitBounds(bounds, {
                padding: [20, 20]
            });
        }
    } catch (error) {
        console.warn("Unable to fit map bounds:", error);
    }

    setTimeout(refreshMapSize, 100);
}

function setLegend(items) {
    const legend = document.getElementById("legend-content");

    if (!legend) {
        return;
    }

    if (!items || items.length === 0) {
        legend.innerHTML = "<div class=\"empty-state\">No legend</div>";
        return;
    }

    legend.innerHTML = items.map((item) => {
        return `
            <div class="legend-row">
                <span class="legend-color"
                      style="background:${item.color};"></span>
                <span>${item.label}</span>
            </div>
        `;
    }).join("");
}

function setDefaultDamageLegend() {
    setLegend([
        {label: "0–5", color: "#ffffff"},
        {label: "5–20", color: "#fee8a8"},
        {label: "20–50", color: "#fdd17a"},
        {label: "50–100", color: "#fdae61"},
        {label: "100–200", color: "#f46d43"},
        {label: "200–500", color: "#d73027"},
        {label: "500+", color: "#a50026"}
    ]);
}

function toggleLegend() {
    const legend = document.getElementById("map-legend");

    if (legend) {
        legend.classList.toggle("collapsed");
    }
}

function bindMapControls() {
    const legendButton = document.getElementById("btn-toggle-legend");

    if (legendButton) {
        legendButton.addEventListener("click", toggleLegend);
    }
}

function refreshMapSize() {
    if (map) {
        map.invalidateSize();
    }
}