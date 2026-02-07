from fragility import (
    LognormalDamageModel,
    DiscreteDamageModel,
    DamageScale,
    IMBounds,
    FragilityModel,
    FragilityCollection,
)

# ------------------------------------------------------------
# 1) Define a damage scale and IM bounds
# ------------------------------------------------------------

damage_scale = DamageScale(
    id="D1D3",
    levels=("D1", "D2", "D3"),
)

im_bounds = IMBounds(min=0.01, max=5.0)

# ------------------------------------------------------------
# 2) Define single-damage-state models
# ------------------------------------------------------------

d1 = LognormalDamageModel(theta=0.05, beta=0.60)
d2 = LognormalDamageModel(theta=0.15, beta=0.60)
d3 = LognormalDamageModel(theta=0.30, beta=0.60)

# ------------------------------------------------------------
# 3) Build a multi-state FragilityModel (lognormal)
# ------------------------------------------------------------

frag_model = FragilityModel(
    id="GENERIC_RC_PGA",
    taxonomy="RC_GENERIC",
    imt="PGA",
    model_type="lognormal_continuous",
    damage_scale=damage_scale,
    im_bounds=im_bounds,
    damage_models={
        "D1": d1,
        "D2": d2,
        "D3": d3,
    },
    metadata={
        "name": "Generic RC fragility (example)",
        "year": 2026,
    },
)

# ------------------------------------------------------------
# 4) Use the model (deterministic IM)
# ------------------------------------------------------------

im = 0.20

print("PoE at IM =", im)
for level, poe in frag_model.poe_all(im).items():
    print(f"  {level}: {poe[0]:.3f}")

# ------------------------------------------------------------
# 5) Convert exceedance to state probabilities
# ------------------------------------------------------------

print("\nState probabilities")
state_p = frag_model.state_probabilities(im)
for level, p in state_p.items():
    print(f"  {level}: {p[0]:.3f}")

# ------------------------------------------------------------
# 6) Add the model to a collection
# ------------------------------------------------------------

collection = FragilityCollection(
    metadata={
        "name": "Example fragility collection",
        "date": "2026-02-07",
    },
    models={
        frag_model.id: frag_model,
    },
)

print("\nAvailable models:", collection.list_ids())

# Access the model by id
model = collection.get("GENERIC_RC_PGA")
print("PoE D2 @ IM=0.2:", model.poe("D2", 0.20)[0])

# ------------------------------------------------------------
# 7) Load fragility models from JSON and use them (with IM
#    uncertainty)
# ------------------------------------------------------------

from fragility import FragilityCollection

# Load the fragility database from JSON
collection_json = FragilityCollection.from_json(
    "fragility_example.json"
)

print("\nModels loaded from JSON:")
print(collection_json.list_ids())

# Select a model by id
model = collection_json.get("GENERIC_RC_LR_PGA_LN_V1")

# ------------------------------------------------------------
# 7a) Deterministic IM (reference case)
# ------------------------------------------------------------

im = 0.30
print("\nDeterministic IM =", im)

for level, poe in model.poe_all(im).items():
    print(f"  {level}: {poe[0]:.3f}")

# ------------------------------------------------------------
# 7b) Lognormal IM uncertainty (GMPE-like)
# ------------------------------------------------------------

im_median = 0.30
im_beta = 0.40  # sigma ln(IM) from GMPE

print("\nLognormal IM uncertainty")
print(f"Median IM = {im_median}, sigma ln(IM) = {im_beta}")

poe_ln = model.poe_lognormal_im_all(im_median, im_beta)

for level, poe in poe_ln.items():
    print(f"  {level}: {poe:.3f}")

# ------------------------------------------------------------
# 7c) Discrete IM distribution (e.g. scenario ensemble)
# ------------------------------------------------------------

im_values = [0.10, 0.20, 0.40]
probabilities = [0.2, 0.5, 0.3]

print("\nDiscrete IM distribution")
print("IM values:", im_values)
print("Probabilities:", probabilities)

poe_disc = model.poe_distribution_all(im_values, probabilities)

for level, poe in poe_disc.items():
    print(f"  {level}: {poe:.3f}")

# ------------------------------------------------------------
# 7d) State probabilities with deterministic IM
# ------------------------------------------------------------

print("\nState probabilities (deterministic IM)")
state_p = model.state_probabilities(im)

for level, p in state_p.items():
    print(f"  {level}: {p[0]:.3f}")
