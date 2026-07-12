"""Main entry point for the ShakeScenario WebUI."""

from __future__ import annotations

from flask import Flask, render_template

from backend.api import create_api_blueprint
from backend.config import (
    ensure_data_root,
    load_config,
)


def create_app() -> Flask:
    """Create and configure the ShakeScenario WebUI application."""
    config = load_config()
    ensure_data_root(config)

    app = Flask(__name__)
    app.config["WEBUI_CONFIG"] = config

    app.register_blueprint(
        create_api_blueprint(config)
    )

    @app.get("/")
    def index():
        """Serve the main application page."""
        return render_template("index.html")

    return app


app = create_app()


if __name__ == "__main__":
    app.run(
        host=app.config["WEBUI_CONFIG"].web_host,
        port=app.config["WEBUI_CONFIG"].web_port,
        debug=app.config["WEBUI_CONFIG"].debug,
    )