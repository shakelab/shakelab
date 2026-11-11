# ****************************************************************************
#
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
Command-line array processing.

The data come from an array file (``--input``, see
:func:`~.io.read_array`), from waveform files read by ShakeLab with
their metadata (``--data``, ``--metadata``), or a saved result is only
plotted (``--result``). The configuration is a JSON file (``--config``,
see :mod:`.config`), whose values the processing options of the
command line replace; without it the options must define at least the
frequencies and the grid. The figures (dispersion, ellipticity, spectrum,
density of the maxima) are chosen with the options or the section
``graphics`` of the configuration file.

Examples
--------
.. code-block:: shell

    slab-yarray --input array.json --config capon.json -o result.pkl \\
        --plot-dispersion --save-figures figures
    slab-yarray --input array.json --config capon.json \\
        --fmin 2 --fmax 8 --nfreq 20 --method music -o music.pkl
    slab-yarray --data "data/*.mseed" --metadata array.xml \\
        --select "XX.*..HH?" --config rtbf.json -o rtbf.pkl
    slab-yarray --result rtbf.pkl --plot-spectrum --frequency 4.5 --show

The command ``slab-yarray`` is installed with ShakeLab; the same
interface runs as ``python -m shakelab.exploration.yarray``.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional

import numpy as np

from shakelab.exploration.yarray.beamforming import Beamformer
from shakelab.exploration.yarray.config import BeamformingConfig
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.io import (
    read_array,
    read_result,
    write_result,
)

__all__ = ["main"]

# Processing options: argument name -> (JSON section, key).
_OPTIONS = {
    "frequency_scale": ("frequencies", "scale"),
    "method": ("method", "name"),
    "diagonal_loading": ("method", "diagonal_loading"),
    "music_sources": ("method", "music_sources"),
    "grid_quantity": ("grid", "quantity"),
    "grid_scale": ("grid", "scale"),
    "grid_min": ("grid", "minimum"),
    "grid_max": ("grid", "maximum"),
    "grid_size": ("grid", "size"),
    "azimuth_size": ("grid", "azimuth_size"),
    "window_type": ("window", "type"),
    "window_length": ("window", "length"),
    "window_taper": ("window", "taper"),
    "window_overlap": ("window", "overlap"),
    "windows_per_block": ("window", "windows_per_block"),
    "block_overlap": ("window", "block_overlap"),
    "ellipticity_angle_min": ("ellipticity", "minimum_angle"),
    "ellipticity_angle_max": ("ellipticity", "maximum_angle"),
    "ellipticity_angle_size": ("ellipticity", "size"),
    "output_mode": ("output", "mode"),
    "maxima": ("output", "maxima"),
    "minimum_relative_power": ("output", "minimum_relative_power"),
    "statistical_uncertainty": ("uncertainty", "statistical"),
    "geometric_uncertainty": ("uncertainty", "geometric"),
    "transform": ("spectral", "transform"),
    "fft_selection": ("spectral", "fft_selection"),
    "fft_padding": ("spectral", "fft_padding"),
}


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------


def main(argv: Optional[list] = None, prog: Optional[str] = None) -> int:
    """
    Run the command-line array processing.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; ``sys.argv[1:]`` if omitted.
    prog : str, optional
        Program name shown in the help.

    Returns
    -------
    int
        Exit status: 0 on success, 1 on errors, 130 when interrupted.
    """
    parser = _parser(prog)
    args = parser.parse_args(argv)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            warnings.showwarning = _show_warning
            return _run(args)
    except KeyboardInterrupt:
        return 130
    except (OSError, ValueError, TypeError, np.linalg.LinAlgError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


def _run(args: argparse.Namespace) -> int:
    """Process or load a result, save it and draw the figures."""
    graphics = _graphics_options(args)

    if args.result is not None:
        for name in ("config", "output", "write_config"):
            if getattr(args, name) is not None:
                raise ValueError(f"--{name.replace('_', '-')} needs data")

        result = read_result(args.result, concatenate=True)
    else:
        config = BeamformingConfig.from_dict(_config_dict(args))

        if args.write_config is not None:
            config.write(args.write_config, owrite=args.owrite)

        array = _array(args)
        print(array.summary())
        result = Beamformer(config).process(array)

        if args.output is not None:
            write_result(result, args.output, append=args.append,
                         owrite=args.owrite)
            print(f"result written to {args.output}")

    print(result.summary())

    for path in _figures(result, graphics):
        print(f"figure written to {path}")

    return 0


# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------


def _parser(prog: Optional[str]) -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog=prog or "slab-yarray",
        description="Frequency-wavenumber processing of seismic arrays.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path,
                        help="array file (JSON) with coordinates and files")
    source.add_argument("--data", nargs="+",
                        help="waveform files, patterns or folders")
    source.add_argument("--result", type=Path,
                        help="saved result to plot (no processing)")

    data = parser.add_argument_group("data (with --data)")
    data.add_argument("--metadata", nargs="+",
                      help="metadata files (positions, orientations)")
    data.add_argument("--select", help="channel pattern NET.STA.LOC.CHA")
    data.add_argument("--start", help="start time (ISO 8601)")
    data.add_argument("--end", help="end time (ISO 8601)")

    output = parser.add_argument_group("configuration and output")
    output.add_argument("--config", type=Path,
                        help="configuration file (JSON)")
    output.add_argument("--write-config", type=Path,
                        help="write the configuration used (JSON)")
    output.add_argument("-o", "--output", type=Path,
                        help="result file (pickle)")
    output.add_argument("--append", action="store_true",
                        help="append the result to the file")
    output.add_argument("--owrite", action="store_true",
                        help="overwrite existing files")

    processing = parser.add_argument_group(
        "processing (replace the values of --config)"
    )
    processing.add_argument("--frequencies", type=float, nargs="+",
                            help="explicit frequencies (Hz)")
    processing.add_argument("--fmin", type=float, help="minimum frequency")
    processing.add_argument("--fmax", type=float, help="maximum frequency")
    processing.add_argument(
        "--nfreq", type=int,
        help="number of frequencies (omitted with --fmin/--fmax and "
        "--transform rfft: every FFT bin)",
    )
    processing.add_argument("--frequency-scale", choices=("log", "linear"))
    processing.add_argument("--method",
                            choices=("bartlett", "capon", "music"))
    processing.add_argument("--components", choices=(
        "vertical", "radial", "transverse", "three_component_separate",
        "rayleigh_joint",
    ))
    processing.add_argument("--diagonal-loading", type=float)
    processing.add_argument("--music-sources", type=int)
    processing.add_argument("--grid-quantity",
                            choices=("velocity", "slowness", "wavenumber"))
    processing.add_argument("--grid-scale", choices=("log", "linear"))
    processing.add_argument("--grid-min", type=float)
    processing.add_argument("--grid-max", type=float)
    processing.add_argument("--grid-size", type=int)
    processing.add_argument("--azimuth-size", type=int)
    processing.add_argument("--window-type", choices=("cycles", "fixed"))
    processing.add_argument("--window-length", type=float,
                            help="cycles, or seconds with fixed windows")
    processing.add_argument("--window-taper", type=float)
    processing.add_argument("--window-overlap", type=float)
    processing.add_argument("--windows-per-block", type=int)
    processing.add_argument("--block-overlap", type=float)
    processing.add_argument("--ellipticity-angle-min", type=float)
    processing.add_argument("--ellipticity-angle-max", type=float)
    processing.add_argument("--ellipticity-angle-size", type=int)
    processing.add_argument("--output-mode",
                            choices=("maxima", "spectra", "both"))
    processing.add_argument("--maxima", type=int)
    processing.add_argument("--minimum-relative-power", type=float)
    processing.add_argument("--statistical-uncertainty",
                            choices=("none", "empirical"))
    processing.add_argument("--geometric-uncertainty",
                            choices=("none", "expected_coherence"))
    processing.add_argument("--strict", action=argparse.BooleanOptionalAction)
    processing.add_argument("--transform", choices=("direct", "rfft"))
    processing.add_argument("--fft-selection", choices=("nearest", "interp"))
    processing.add_argument("--fft-padding", type=int)

    figures = parser.add_argument_group(
        "figures (replace the section 'graphics' of --config)"
    )
    figures.add_argument("--plot-dispersion",
                         action=argparse.BooleanOptionalAction)
    figures.add_argument("--plot-ellipticity",
                         action=argparse.BooleanOptionalAction)
    figures.add_argument("--plot-spectrum",
                         action=argparse.BooleanOptionalAction)
    figures.add_argument("--plot-density",
                         action=argparse.BooleanOptionalAction,
                         help="density of the maxima versus frequency")
    figures.add_argument(
        "--density-coordinate",
        choices=("log_slowness", "log_velocity", "wavenumber"),
        help="coordinate of the density (log_slowness)",
    )
    figures.add_argument("--component", help="component of the figures")
    figures.add_argument("--ranks", type=int, nargs="+",
                         help="ranks of the maxima")
    figures.add_argument("--frequency", type=float,
                         help="frequency of the spectrum (Hz)")
    figures.add_argument("--block", type=int,
                         help="covariance block of the spectrum")
    figures.add_argument("--spectrum-view", choices=("kxky", "azimuth"))
    figures.add_argument(
        "--spectrum-quantity",
        choices=("power", "standard_error", "ellipticity"),
    )
    figures.add_argument("--show", action=argparse.BooleanOptionalAction,
                         help="show the figures")
    figures.add_argument("--save-figures", type=Path,
                         help="folder of the figures")
    figures.add_argument("--figure-format", help="image format (png)")
    figures.add_argument("--dpi", type=int, help="resolution (150)")

    return parser


def _config_dict(args: argparse.Namespace) -> dict:
    """Configuration from the file, with the options replacing values."""
    data: dict = {}

    if args.config is not None:
        data = json.loads(args.config.read_text(encoding="utf-8"))

        if not isinstance(data, Mapping):
            raise ValueError(f"{args.config}: not a JSON object")

        data = {key: value for key, value in data.items()
                if key != "graphics"}

    if isinstance(data.get("method"), str):
        data["method"] = {"name": data["method"]}

    if args.frequencies is not None:
        data["frequencies"] = {"values": args.frequencies}
    elif args.fmin is not None or args.fmax is not None:
        if args.fmin is None or args.fmax is None:
            raise ValueError("--fmin and --fmax go together")

        data["frequencies"] = {
            "minimum": args.fmin, "maximum": args.fmax, "size": args.nfreq,
        }
    elif args.nfreq is not None:
        data.setdefault("frequencies", {})["size"] = args.nfreq

    for name, (section, key) in _OPTIONS.items():
        value = getattr(args, name)

        if value is not None:
            data.setdefault(section, {})[key] = value

    for name in ("components", "strict"):
        if getattr(args, name) is not None:
            data[name] = getattr(args, name)

    if "frequencies" not in data or "grid" not in data:
        raise ValueError(
            "the configuration needs frequencies and grid: give --config "
            "or --frequencies/--fmin --fmax and --grid-min --grid-max "
            "--grid-size"
        )

    return data


def _array(args: argparse.Namespace) -> ArrayData:
    """Array data from an array file or waveform files."""
    if args.input is not None:
        return read_array(args.input)

    from shakelab.signals.io import metadata_reader, reader

    metadata = None if args.metadata is None else metadata_reader(
        args.metadata
    )
    collection = reader(
        args.data, metadata=metadata, select=args.select,
        starttime=args.start, endtime=args.end,
    )
    return ArrayData.from_collection(collection, select=args.select)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def _graphics_options(args: argparse.Namespace) -> dict:
    """Figure options: section 'graphics' of the config, then options."""
    section: dict = {}

    if args.config is not None and args.config.is_file():
        data = json.loads(args.config.read_text(encoding="utf-8"))
        section = dict(data.get("graphics", {})) if isinstance(
            data, Mapping
        ) else {}

    spectrum = section.get("spectrum", {})
    spectrum = {"enabled": spectrum} if isinstance(spectrum, bool) else dict(
        spectrum
    )
    density = section.get("density", {})
    density = {"enabled": density} if isinstance(density, bool) else dict(
        density
    )
    save = section.get("save", {})
    save = {"enabled": save} if isinstance(save, bool) else dict(save)

    def pick(value: Any, default: Any) -> Any:
        return default if value is None else value

    directory = pick(args.save_figures, save.get("directory"))

    return {
        "dispersion": pick(args.plot_dispersion,
                           bool(section.get("dispersion", False))),
        "ellipticity": pick(args.plot_ellipticity,
                            bool(section.get("ellipticity", False))),
        "spectrum": pick(args.plot_spectrum,
                         bool(spectrum.get("enabled", False))),
        "density": pick(args.plot_density,
                        bool(density.get("enabled", False))),
        "coordinate": pick(args.density_coordinate,
                           density.get("coordinate", "log_slowness")),
        "component": pick(args.component, section.get("component")),
        "ranks": pick(args.ranks, section.get("ranks")),
        "frequency": pick(args.frequency, spectrum.get("frequency")),
        "block": pick(args.block, spectrum.get("block", 0)),
        "view": pick(args.spectrum_view, spectrum.get("view", "kxky")),
        "quantity": pick(args.spectrum_quantity,
                         spectrum.get("quantity", "power")),
        "show": pick(args.show, bool(section.get("show", False))),
        "directory": None if directory is None else Path(directory),
        "save": args.save_figures is not None or bool(
            save.get("enabled", directory is not None)
        ),
        "format": pick(args.figure_format, save.get("format", "png")),
        "dpi": pick(args.dpi, save.get("dpi", 150)),
    }


def _figures(result: Any, options: dict) -> list:
    """Draw the requested figures; return the files written."""
    wanted = [
        name for name in ("dispersion", "ellipticity", "spectrum", "density")
        if options[name]
    ]

    if not wanted:
        return []

    if options["save"] and options["directory"] is None:
        raise ValueError("saving figures needs --save-figures DIR")

    if "spectrum" in wanted and options["frequency"] is None:
        raise ValueError("the spectrum figure needs --frequency")

    from shakelab.graphics import yarray as graphics

    plt = graphics._pyplot()
    figures = []
    component = options["component"]

    if "dispersion" in wanted:
        figure, _ = graphics.plot_dispersion(
            result, component=component, ranks=options["ranks"]
        )
        figures.append((figure, f"dispersion_{component or 'all'}"))

    if "ellipticity" in wanted:
        figure, _ = graphics.plot_ellipticity(result, ranks=options["ranks"])
        figures.append((figure, "ellipticity"))

    if "spectrum" in wanted:
        figure, _ = graphics.plot_spectrum(
            result, options["frequency"], block=options["block"],
            component=component, view=options["view"],
            quantity=options["quantity"], ranks=options["ranks"],
        )
        figures.append((
            figure,
            f"spectrum_f{options['frequency']:g}_b{options['block']}_"
            f"{options['view']}_{options['quantity']}",
        ))

    if "density" in wanted:
        from shakelab.exploration.yarray.statistics import maximum_density

        density = maximum_density(
            result, component=component, ranks=options["ranks"],
            coordinate=options["coordinate"],
        )
        figure, _ = graphics.plot_density(
            density, peaks=density.peaks(3, minimum_relative_density=0.2)
        )
        figures.append((figure, f"density_{density.component}"))

    written = []

    if options["save"]:
        options["directory"].mkdir(parents=True, exist_ok=True)

        for figure, stem in figures:
            path = options["directory"] / f"{stem}.{options['format']}"
            figure.savefig(path, dpi=options["dpi"], bbox_inches="tight")
            written.append(path)

    if options["show"]:
        plt.show()

    for figure, _ in figures:
        plt.close(figure)

    return written


def _show_warning(message: Any, category: Any, *args: Any,
                  **kwargs: Any) -> None:
    """Print warnings as one line."""
    print(f"warning: {message}", file=sys.stderr)
