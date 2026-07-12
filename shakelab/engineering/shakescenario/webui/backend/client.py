# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# ****************************************************************************
"""
ShakeScenario client wrapper.

This module provides a lightweight adapter between the WebUI backend and the
ShakeScenario client library.
"""

from __future__ import annotations

from typing import Any

from shakelab.engineering.shakescenario.shakeclient import ShakeScenarioClient

from .config import WebUIConfig


class Client:
    """Wrapper around the ShakeScenario client."""

    def __init__(self, config: WebUIConfig) -> None:
        self._client = ShakeScenarioClient(
            host=config.server_host,
            port=config.server_port,
            timeout=config.server_timeout,
        )

    def ping(self) -> dict[str, Any]:
        """Ping the ShakeScenario server."""
        return self._client.ping()

    def list_models(self) -> dict[str, Any]:
        """Return the registered exposure models."""
        return self._client.models_list()

    def submit(
        self,
        *,
        tag: str | None,
        model_id: str,
        origin_time: str,
        magnitude: float,
        longitude: float,
        latitude: float,
        depth: float,
        gmpe: str | None = None,
        distance_approx: str | None = None,
    ) -> dict[str, Any]:
        """
        Submit a new scenario.

        Parameters
        ----------
        depth
            Hypocentral depth in kilometres.
        """
        payload: dict[str, Any] = {
            "models": {
                "model_id": model_id,
            },
            "scenario": {
                "event": {
                    "origin_time": origin_time,
                    "magnitude": float(magnitude),
                    "hypocentre": {
                        "longitude": float(longitude),
                        "latitude": float(latitude),
                        "elevation": -1000.0 * float(depth),
                    },
                },
            },
        }

        if tag:
            payload["tag"] = tag

        if gmpe is not None or distance_approx is not None:
            ground_motion: dict[str, Any] = {
                "provider": "gmpe",
            }

            if gmpe is not None:
                ground_motion["gmpe_name"] = gmpe

            if distance_approx is not None:
                ground_motion["distance_approx"] = distance_approx

            payload["scenario"]["ground_motion"] = ground_motion

        return self._client.submit(payload)
