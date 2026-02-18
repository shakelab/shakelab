from shakelab.engineering.taxonomy.taxonomy_tree import TaxonomyTree
from shakelab.engineering.fragility.fragility import FragilityCollection

# Load taxonomy → fragility mapping
tree = TaxonomyTree.from_json("taxonomy_tree_example.json")

# Load fragility models
fragilities = FragilityCollection.from_json(
    "../fragility/fragility_example.json"
)

# Resolve taxonomy into fragility models and weights
models = tree.resolve("RC-LRS", fragilities)

print(models)