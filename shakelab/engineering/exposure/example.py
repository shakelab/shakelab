from exposure import Exposure

# ------------------------------------------------------------
# 1) Load exposure model from JSON
# ------------------------------------------------------------

exposure = Exposure.from_json(
    "exposure_example.json",
    validate=True,
)

print("Exposure metadata:")
print(exposure.metadata)

print("\nNumber of assets:", len(exposure.assets))

# ------------------------------------------------------------
# 2) List taxonomies present in the exposure
# ------------------------------------------------------------

taxonomies = exposure.list_taxonomies()

print("\nTaxonomies present:")
for tax in taxonomies:
    print(" ", tax)

# ------------------------------------------------------------
# 3) Inspect one asset and its typologies
# ------------------------------------------------------------

asset = exposure.assets[0]

print("\nExample asset:")
print("  id:", asset.id)
print("  aggregated:", asset.aggregated)

loc = asset.reference_location
lon = loc.get("longitude")
lat = loc.get("latitude")

print("  reference location:", lon, lat)

print("\nTypologies for this asset:")
for typ in asset.typologies:
    print(
        f"  taxonomy={typ.taxonomy} | "
        f"count={typ.count}"
    )

# ------------------------------------------------------------
# 4) Typical loop over assets (downstream usage)
# ------------------------------------------------------------

print("\nAssets summary:")
for asset in exposure.assets:
    loc = asset.reference_location
    lon = loc.get("longitude")
    lat = loc.get("latitude")

    for typ in asset.typologies:
        print(
            f"  asset={asset.id} | "
            f"taxonomy={typ.taxonomy} | "
            f"count={typ.count} | "
            f"lon={lon} | lat={lat}"
        )

# ------------------------------------------------------------
# 5) Save exposure back to JSON
# ------------------------------------------------------------

exposure.to_json(
    "exposure_out.json",
    validate=True,
    indent=2,
)

print("\nSaved to exposure_out.json")
