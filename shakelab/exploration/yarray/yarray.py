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
Command-line interface for YArray processing and visualization.

This module is the user-facing command-line wrapper around the
:class:`shakelab.exploration.yarray.runner.Yarray` processing engine
and the reusable plotting functions defined in ``graphics.py``.

Two workflows are supported:

1. process waveform files and optionally visualize the new result;
2. load an existing serialized result and visualize it without
   repeating the numerical processing.

The JSON configuration used for processing may additionally contain a
top-level ``graphics`` section. Processing options, including direct or
RFFT spectral transforms, are resolved by ``io.py``; this module reads
only the optional graphics section. The same configuration file
therefore remains directly usable by the Python ``Yarray`` API.

Examples
--------
Process and save a result, applying graphics settings from the JSON
configuration:

    python3 -m shakelab.exploration.yarray.yarray \
        --input array_config.json \
        --config yarray_config.json \
        --output result.pkl

Visualize an existing result using the same graphics configuration:

    python3 -m shakelab.exploration.yarray.yarray \
        --result result.pkl \
        --config yarray_config.json

Command-line graphics options override values read from the JSON file.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from shakelab.exploration.yarray.graphics import (
    plot_dispersion,
    plot_ellipticity,
    plot_spectrum,
)
from shakelab.exploration.yarray.io import read_result
from shakelab.exploration.yarray.runner import Yarray


# Absolute package imports are intentional here. Unlike the library modules,
# this file is also designed to work as a directly executed CLI script.
_DEFAULT_FIGURE_FORMAT = "png"
_DEFAULT_FIGURE_DPI = 170


def _build_parser() -> argparse.ArgumentParser:
    """
    Build the YArray command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured command-line parser.
    """
    parser = argparse.ArgumentParser(
        prog="yarray",
        description=(
            "Run seismic-array processing or visualize saved results."
        ),
    )

    source = parser.add_mutually_exclusive_group(
        required=True
    )

    source.add_argument(
        "--input",
        type=Path,
        help=(
            "Array input JSON configuration. Requires --config and "
            "runs the numerical processing."
        ),
    )

    source.add_argument(
        "--result",
        type=Path,
        help=(
            "Existing serialized BeamformingResult to visualize "
            "without repeating the processing."
        ),
    )

    parser.add_argument(
        "--config",
        type=Path,
        help=(
            "YArray JSON configuration. Required with --input; "
            "it controls processing, spectral transforms, and optional "
            "graphics settings. With --result it is optional when all "
            "graphics options are supplied directly."
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Optional serialized result file. Valid only with --input."
        ),
    )

    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to an existing result file.",
    )

    parser.add_argument(
        "--plot-dispersion",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Enable or disable the frequency-velocity plot.",
    )

    parser.add_argument(
        "--plot-ellipticity",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Enable or disable the frequency-ellipticity plot.",
    )

    parser.add_argument(
        "--plot-spectrum",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Enable or disable one detailed FK spectrum plot.",
    )

    parser.add_argument(
        "--component",
        type=str,
        help="Component used by dispersion and spectrum plots.",
    )

    parser.add_argument(
        "--ranks",
        type=int,
        nargs="+",
        help="One-based spectral-maximum ranks to display.",
    )

    parser.add_argument(
        "--frequency",
        type=float,
        help="Requested frequency for the detailed spectrum.",
    )

    parser.add_argument(
        "--block",
        type=int,
        help="Covariance-block index for the detailed spectrum.",
    )

    parser.add_argument(
        "--spectrum-view",
        choices=(
            "kxky",
            "azimuth",
        ),
        help="Detailed spectrum representation.",
    )

    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Show generated figures interactively.",
    )

    parser.add_argument(
        "--save-figures",
        type=Path,
        help="Directory where generated figures are saved.",
    )

    parser.add_argument(
        "--figure-format",
        type=str,
        help="Image format used when saving figures.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        help="Resolution used when saving figures.",
    )

    return parser


def _read_graphics_config(
    file_path: str | Path | None,
) -> dict[str, Any]:
    """
    Read the optional top-level graphics section from a JSON file.

    Parameters
    ----------
    file_path : path-like or None
        YArray configuration file.

    Returns
    -------
    dict
        Graphics configuration. An empty dictionary is returned when
        no configuration file or no ``graphics`` section is present.
    """
    if file_path is None:
        return {}

    path = Path(file_path).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"YArray configuration not found: {path}"
        )

    try:
        with path.open(
            "r",
            encoding="utf-8",
        ) as stream:
            data = json.load(stream)

    except json.JSONDecodeError as error:
        raise ValueError(
            "Invalid YArray JSON configuration at "
            f"line {error.lineno}, column {error.colno}: "
            f"{error.msg}"
        ) from error

    if not isinstance(data, Mapping):
        raise TypeError(
            "The YArray configuration root must be a JSON object."
        )

    graphics = data.get(
        "graphics",
        {},
    )

    if not isinstance(graphics, Mapping):
        raise TypeError(
            "Configuration field 'graphics' must be a JSON object."
        )

    return dict(graphics)


def _graphics_options(
    args: argparse.Namespace,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """
    Merge graphics JSON settings with command-line overrides.

    Command-line values have precedence over configuration-file values.
    """
    spectrum = config.get(
        "spectrum",
        {},
    )

    if isinstance(spectrum, bool):
        spectrum = {
            "enabled": spectrum,
        }

    if not isinstance(spectrum, Mapping):
        raise TypeError(
            "graphics.spectrum must be a Boolean or JSON object."
        )

    save = config.get(
        "save",
        {},
    )

    if isinstance(save, bool):
        save = {
            "enabled": save,
        }

    if not isinstance(save, Mapping):
        raise TypeError(
            "graphics.save must be a Boolean or JSON object."
        )

    ranks = (
        args.ranks
        if args.ranks is not None
        else config.get("ranks")
    )

    if ranks is not None:
        if (
            not isinstance(ranks, Sequence)
            or isinstance(ranks, (str, bytes))
        ):
            raise TypeError(
                "graphics.ranks must be an array of positive integers."
            )

        ranks = tuple(
            int(rank)
            for rank in ranks
        )

        if not ranks or any(
            rank < 1
            for rank in ranks
        ):
            raise ValueError(
                "Graphics ranks must contain positive integers."
            )

    save_directory = (
        args.save_figures
        if args.save_figures is not None
        else save.get("directory")
    )

    save_enabled = bool(
        save.get(
            "enabled",
            save_directory is not None,
        )
    )

    if args.save_figures is not None:
        save_enabled = True

    return {
        "dispersion": (
            args.plot_dispersion
            if args.plot_dispersion is not None
            else bool(config.get("dispersion", False))
        ),
        "ellipticity": (
            args.plot_ellipticity
            if args.plot_ellipticity is not None
            else bool(config.get("ellipticity", False))
        ),
        "spectrum": (
            args.plot_spectrum
            if args.plot_spectrum is not None
            else bool(spectrum.get("enabled", False))
        ),
        "component": (
            args.component
            if args.component is not None
            else config.get("component")
        ),
        "ranks": ranks,
        "frequency": (
            args.frequency
            if args.frequency is not None
            else spectrum.get("frequency")
        ),
        "block": (
            args.block
            if args.block is not None
            else int(spectrum.get("block", 0))
        ),
        "spectrum_view": (
            args.spectrum_view
            if args.spectrum_view is not None
            else spectrum.get("view", "kxky")
        ),
        "show": (
            args.show
            if args.show is not None
            else bool(config.get("show", False))
        ),
        "save": save_enabled,
        "save_directory": save_directory,
        "figure_format": (
            args.figure_format
            if args.figure_format is not None
            else save.get(
                "format",
                _DEFAULT_FIGURE_FORMAT,
            )
        ),
        "dpi": (
            args.dpi
            if args.dpi is not None
            else int(
                save.get(
                    "dpi",
                    _DEFAULT_FIGURE_DPI,
                )
            )
        ),
    }


def _validate_graphics_options(
    options: Mapping[str, Any],
) -> None:
    """Validate merged graphics settings."""
    if options["spectrum"]:
        frequency = options["frequency"]

        if frequency is None:
            raise ValueError(
                "A spectrum plot requires a frequency."
            )

        if not np.isfinite(frequency) or frequency <= 0.0:
            raise ValueError(
                "Spectrum frequency must be finite and positive."
            )

        if options["block"] < 0:
            raise ValueError(
                "Spectrum block must be non-negative."
            )

        if options["spectrum_view"] not in {
            "kxky",
            "azimuth",
        }:
            raise ValueError(
                "Spectrum view must be 'kxky' or 'azimuth'."
            )

    if options["save"]:
        if options["save_directory"] is None:
            raise ValueError(
                "Saving figures requires a destination directory."
            )

        if options["dpi"] <= 0:
            raise ValueError(
                "Figure DPI must be positive."
            )


def _has_graphics(
    options: Mapping[str, Any],
) -> bool:
    """Return whether at least one plot has been requested."""
    return bool(
        options["dispersion"]
        or options["ellipticity"]
        or options["spectrum"]
    )


def _save_figure(
    figure,
    *,
    directory: Path,
    stem: str,
    file_format: str,
    dpi: int,
) -> Path:
    """Save one Matplotlib figure and return the written path."""
    directory = directory.expanduser().resolve()
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    file_format = file_format.lstrip(".")
    path = directory / f"{stem}.{file_format}"

    figure.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
    )

    return path


def _generate_graphics(
    result,
    options: Mapping[str, Any],
) -> tuple[Path, ...]:
    """
    Generate requested figures from a BeamformingResult.

    Returns
    -------
    tuple of Path
        Paths of figures written to disk.
    """
    written: list[Path] = []
    figures = []

    component = options["component"]
    ranks = options["ranks"]

    if options["dispersion"]:
        figure, _ = plot_dispersion(
            result,
            component=component,
            ranks=ranks,
        )

        figures.append(
            (
                figure,
                (
                    "dispersion"
                    if component is None
                    else f"dispersion_{component}"
                ),
            )
        )

    if options["ellipticity"]:
        figure, _ = plot_ellipticity(
            result,
            component=(
                "rayleigh_joint"
                if component is None
                else component
            ),
            ranks=ranks,
        )

        figures.append(
            (
                figure,
                (
                    "ellipticity"
                    if component is None
                    else f"ellipticity_{component}"
                ),
            )
        )

    if options["spectrum"]:
        figure, _ = plot_spectrum(
            result,
            frequency=options["frequency"],
            block=options["block"],
            component=component,
            view=options["spectrum_view"],
            ranks=ranks,
        )

        frequency = float(
            options["frequency"]
        )

        stem = (
            "spectrum"
            f"_f{frequency:g}"
            f"_b{options['block']}"
            f"_{options['spectrum_view']}"
        )

        if component is not None:
            stem += f"_{component}"

        figures.append(
            (
                figure,
                stem,
            )
        )

    if options["save"]:
        directory = Path(
            options["save_directory"]
        )

        for figure, stem in figures:
            written.append(
                _save_figure(
                    figure,
                    directory=directory,
                    stem=stem,
                    file_format=options["figure_format"],
                    dpi=options["dpi"],
                )
            )

    if options["show"]:
        plt.show()
    else:
        for figure, _ in figures:
            plt.close(
                figure
            )

    return tuple(
        written
    )


def _load_or_process(
    args: argparse.Namespace,
):
    """Return a BeamformingResult from processing or serialized input."""
    if args.input is not None:
        if args.config is None:
            raise ValueError(
                "--input requires --config."
            )

        processor = Yarray(
            args.config
        )

        return processor.run_from_files(
            args.input,
            output_file=args.output,
            append=args.append,
        )

    if args.output is not None:
        raise ValueError(
            "--output is valid only when processing --input."
        )

    if args.append:
        raise ValueError(
            "--append is valid only when processing --input."
        )

    return read_result(
        args.result
    )


def main(argv: Sequence[str] | None = None) -> int:
    """
    Run the YArray command-line interface.

    Parameters
    ----------
    argv : sequence of str, optional
        Command-line arguments without the executable name. When
        omitted, arguments are read from ``sys.argv``.

    Returns
    -------
    int
        Process exit code.
    """
    parser = _build_parser()
    args = parser.parse_args(
        argv
    )

    try:
        graphics_config = _read_graphics_config(
            args.config
        )

        options = _graphics_options(
            args,
            graphics_config,
        )

        _validate_graphics_options(
            options
        )

        result = _load_or_process(
            args
        )

        written = ()

        if _has_graphics(options):
            written = _generate_graphics(
                result,
                options,
            )

    except (
        FileNotFoundError,
        TypeError,
        ValueError,
        IndexError,
        np.linalg.LinAlgError,
    ) as error:
        parser.error(
            str(error)
        )

    if args.input is not None:
        print(
            "YArray processing completed."
        )
    else:
        print(
            "YArray result loaded."
        )

    print(
        f"Maxima: {len(result.maxima)}"
    )

    if result.spectra is None:
        print(
            "Spectra: not retained"
        )
    else:
        print(
            f"Spectra: {len(result.spectra)}"
        )

    if args.output is not None:
        print(
            f"Result written to: {args.output}"
        )

    for path in written:
        print(
            f"Figure written to: {path}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
