"""Main entry point for the ShakeScenario WebUI."""

from __future__ import annotations

from flask import Flask, render_template

from shakelab.engineering.shakescenario.webui.backend.api import (
    create_api_blueprint,
)
from shakelab.engineering.shakescenario.webui.backend.config import (
    ensure_data_root,
    load_config,
    WebUIConfig,
)


def create_app(config: WebUIConfig | None = None) -> Flask:
    """Create the WebUI application using an optional explicit config."""
    if config is None:
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


def main() -> None:
    """Run the WebUI development server."""
    app = create_app()
    app.run(
        host=app.config["WEBUI_CONFIG"].web_host,
        port=app.config["WEBUI_CONFIG"].web_port,
        debug=app.config["WEBUI_CONFIG"].debug,
    )

if __name__ == "__main__":
    main()
