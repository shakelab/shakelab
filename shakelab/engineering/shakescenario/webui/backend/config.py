"""Configuration utilities for the ShakeScenario WebUI."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_HOME = "shakescenario_data"
DEFAULT_SERVER_HOST = "127.0.0.1"
DEFAULT_SERVER_PORT = 5001
DEFAULT_SERVER_TIMEOUT = 10.0
DEFAULT_WEB_HOST = "127.0.0.1"
DEFAULT_WEB_PORT = 8000
DEFAULT_DEBUG = True
DEFAULT_CONFIG_FILE = "webui_config.json"


@dataclass(frozen=True)
class WebUIConfig:
    """Runtime configuration for the WebUI backend."""

    home: Path
    server_host: str
    server_port: int
    server_timeout: float
    web_host: str
    web_port: int
    debug: bool

    @property
    def db_dir(self) -> Path:
        """Return the database directory."""
        return self.home / "db"

    @property
    def logs_dir(self) -> Path:
        """Return the logs directory."""
        return self.home / "logs"

    @property
    def models_dir(self) -> Path:
        """Return the models directory."""
        return self.home / "models"

    @property
    def runs_dir(self) -> Path:
        """Return the calculated runs directory."""
        return self.home / "runs"

    @property
    def database_path(self) -> Path:
        """Return the SQLite database path."""
        return self.db_dir / "shakescenario.db"


def _read_config_file() -> dict[str, Any]:
    """Read the optional WebUI JSON configuration file."""
    path = Path(os.environ.get(
        "SHAKESCENARIO_WEBUI_CONFIG",
        DEFAULT_CONFIG_FILE,
    ))

    if not path.exists():
        return {}

    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)

    if not isinstance(data, dict):
        raise ValueError(
            f"WebUI configuration file must contain a JSON object: {path}"
        )

    return data


def _get_value(
    file_config: dict[str, Any],
    env_name: str,
    file_name: str,
    default: Any,
) -> Any:
    """Return configuration value using env > file > default priority."""
    if env_name in os.environ:
        return os.environ[env_name]

    if file_name in file_config:
        return file_config[file_name]

    return default


def _as_bool(value: Any) -> bool:
    """Convert common configuration values to bool."""
    if isinstance(value, bool):
        return value

    if isinstance(value, str):
        return value.strip().lower() in (
            "1",
            "true",
            "yes",
            "on",
        )

    return bool(value)


def load_config() -> WebUIConfig:
    """Load the WebUI runtime configuration.

    Configuration values are resolved using the following precedence:

    1. Environment variables.
    2. Optional ``webui_config.json`` configuration file.
    3. Built-in defaults.

    The configuration includes both the ShakeScenario data location and the
    network settings required by the WebUI and the ShakeScenario server.
    """
    file_config = _read_config_file()

    home = _get_value(
        file_config,
        "SHAKESCENARIO_HOME",
        "home",
        DEFAULT_HOME,
    )

    server_host = _get_value(
        file_config,
        "SHAKESCENARIO_SERVER_HOST",
        "server_host",
        DEFAULT_SERVER_HOST,
    )

    server_port = int(_get_value(
        file_config,
        "SHAKESCENARIO_SERVER_PORT",
        "server_port",
        DEFAULT_SERVER_PORT,
    ))

    server_timeout = float(_get_value(
        file_config,
        "SHAKESCENARIO_SERVER_TIMEOUT",
        "server_timeout",
        DEFAULT_SERVER_TIMEOUT,
    ))

    web_host = _get_value(
        file_config,
        "SHAKESCENARIO_WEB_HOST",
        "web_host",
        DEFAULT_WEB_HOST,
    )

    web_port = int(_get_value(
        file_config,
        "SHAKESCENARIO_WEB_PORT",
        "web_port",
        DEFAULT_WEB_PORT,
    ))

    debug = _as_bool(_get_value(
        file_config,
        "SHAKESCENARIO_WEB_DEBUG",
        "debug",
        DEFAULT_DEBUG,
    ))

    return WebUIConfig(
        home=Path(home).expanduser().resolve(),
        server_host=str(server_host),
        server_port=server_port,
        server_timeout=server_timeout,
        web_host=str(web_host),
        web_port=web_port,
        debug=debug,
    )


def ensure_data_root(config: WebUIConfig) -> None:
    """Check that the configured data root exists.

    Raises
    ------
    FileNotFoundError
        If the configured ShakeScenario data root does not exist.
    """
    if not config.home.exists():
        raise FileNotFoundError(
            f"ShakeScenario data root does not exist: {config.home}"
        )

    if not config.runs_dir.exists():
        raise FileNotFoundError(
            f"ShakeScenario runs directory does not exist: "
            f"{config.runs_dir}"
        )