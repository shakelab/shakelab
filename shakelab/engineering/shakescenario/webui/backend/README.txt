1. config.py      → legge SHAKESCENARIO_HOME
2. models.py      → dataclass per Run, Event, Layer, AssetImpact
3. repository.py  → legge runs/, meta.json, manifest, impact
4. api.py         → definisce le route Flask /api/v1/*
5. app.py         → crea Flask app e registra API