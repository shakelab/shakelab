# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
High-level Python interface for YArray processing.

This module defines the configurable :class:`Yarray` processing engine.
It provides the application-level bridge between the YArray data model,
the numerical beamforming core, and the input/output layer.

A processing configuration may be supplied when creating the object,
changed later with :meth:`Yarray.set_config`, or overridden for a single
processing call.

The class stores only the default processing configuration. Processing
results are returned explicitly and are never retained as mutable object
state.

Examples
--------
Configure at construction time:

    processor = Yarray("rtbf_processing.json")
    result = processor.run_from_files("array_config.json")

Configure later:

    processor = Yarray()
    processor.set_config(config)
    result = processor.run(array)

Override the default configuration for one call:

    result = processor.run(
        array,
        config="alternative_processing.json",
    )
"""

from __future__ import annotations

from pathlib import Path

from .beamforming import BeamformingConfig, BeamformingResult, beamform
from .data import ArrayData
from .io import read_from_files, read_processing_config, write_result


ConfigInput = BeamformingConfig | str | Path


class Yarray:
    """
    Configurable high-level YArray processing engine.

    Parameters
    ----------
    config : BeamformingConfig or path-like, optional
        Default processing configuration. A path is interpreted as a
        JSON processing-configuration file and may select direct or
        RFFT spectral processing.

    Notes
    -----
    A configuration supplied directly to :meth:`run` or
    :meth:`run_from_files` overrides the stored default configuration
    for that call only.
    """

    def __init__(self, config: ConfigInput | None = None) -> None:
        self._config: BeamformingConfig | None = None

        if config is not None:
            self.set_config(config)

    @property
    def config(self) -> BeamformingConfig | None:
        """
        Return the stored default processing configuration.

        Returns
        -------
        BeamformingConfig or None
            Current default configuration, or ``None`` when the
            processor has not yet been configured.
        """
        return self._config

    @property
    def is_configured(self) -> bool:
        """
        Return whether a default processing configuration is available.

        Returns
        -------
        bool
            ``True`` when a default configuration has been set.
        """
        return self._config is not None

    def set_config(self, config: ConfigInput) -> BeamformingConfig:
        """
        Set the default processing configuration.

        Parameters
        ----------
        config : BeamformingConfig or path-like
            Processing configuration or JSON configuration file.
            Spectral-transform options are resolved by the I/O layer.

        Returns
        -------
        BeamformingConfig
            Resolved configuration stored by the processor.
        """
        resolved = self._resolve_config(config)
        self._config = resolved

        return resolved

    def clear_config(self) -> None:
        """
        Remove the stored default processing configuration.

        After this call, processing methods require an explicit
        ``config`` argument until a new default configuration is set.
        """
        self._config = None

    def run(
        self,
        array: ArrayData,
        *,
        config: ConfigInput | None = None,
        output_file: str | Path | None = None,
        append: bool = False,
    ) -> BeamformingResult:
        """
        Process an in-memory array snapshot.

        Parameters
        ----------
        array : ArrayData
            Validated array-data snapshot.
        config : BeamformingConfig or path-like, optional
            Per-call processing configuration. When omitted, the
            default configuration stored in the instance is used.
        output_file : path-like, optional
            File where the result is serialized. When omitted, the
            result is returned without writing to disk.
        append : bool, optional
            Append the serialized result to an existing result file.
            This option requires ``output_file``.

        Returns
        -------
        BeamformingResult
            Processing result.
        """
        if not isinstance(array, ArrayData):
            raise TypeError(
                "array must be an ArrayData instance."
            )

        resolved_config = self._config_for_run(config)

        if output_file is None and append:
            raise ValueError(
                "append=True requires output_file."
            )

        result = beamform(
            array,
            resolved_config,
        )

        if output_file is not None:
            write_result(
                result,
                output_file,
                append=append,
            )

        return result

    def run_from_files(
        self,
        input_file: str | Path,
        *,
        config: ConfigInput | None = None,
        output_file: str | Path | None = None,
        append: bool = False,
    ) -> BeamformingResult:
        """
        Read an array snapshot from files and process it.

        Parameters
        ----------
        input_file : path-like
            Array input JSON configuration passed to
            :func:`yarray.io.read_from_files`.
        config : BeamformingConfig or path-like, optional
            Per-call processing configuration. When omitted, the
            default configuration stored in the instance is used.
        output_file : path-like, optional
            File where the result is serialized.
        append : bool, optional
            Append the serialized result to an existing result file.

        Returns
        -------
        BeamformingResult
            Processing result.
        """
        array = read_from_files(input_file)

        return self.run(
            array,
            config=config,
            output_file=output_file,
            append=append,
        )

    def _config_for_run(
        self,
        config: ConfigInput | None,
    ) -> BeamformingConfig:
        """
        Resolve a per-call configuration or the stored default.

        Parameters
        ----------
        config : BeamformingConfig or path-like or None
            Optional per-call configuration.

        Returns
        -------
        BeamformingConfig
            Configuration to use for the processing call.
        """
        if config is not None:
            return self._resolve_config(config)

        if self._config is None:
            raise ValueError(
                "No processing configuration is available. "
                "Provide config=... or call set_config() first."
            )

        return self._config

    @staticmethod
    def _resolve_config(config: ConfigInput) -> BeamformingConfig:
        """
        Resolve one supported configuration representation.

        Parameters
        ----------
        config : BeamformingConfig or path-like
            Processing configuration or JSON file.

        Returns
        -------
        BeamformingConfig
            Resolved processing configuration.
        """
        if isinstance(config, BeamformingConfig):
            return config

        if isinstance(config, (str, Path)):
            return read_processing_config(config)

        raise TypeError(
            "config must be a BeamformingConfig or JSON file path."
        )
