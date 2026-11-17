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

The command ``slab-yarray`` has three subcommands, with the same
options for data, configuration and figures where they apply:

- ``beamforming``: frequency-wavenumber processing (dispersion,
  ellipticity, spectra, density of the maxima);
- ``fdd``: frequency domain decomposition (singular values, peaks,
  mode shapes) and, with ``--efdd``, the enhanced FDD (natural
  frequencies and damping ratios);
- ``response``: geometric response and limits of the array.

The data come from an array file (``--input``, see
:meth:`~.data.ArrayData.read`) or from waveform files read by ShakeLab
with their metadata (``--data``, ``--metadata``); a saved result can be
only plotted (``--result``), and the response needs only the
coordinates (``--coordinates``). The configuration is a JSON file
(``--config``, see :mod:`.config`), whose values the processing options
of the command line replace; without it the options must define the
required parameters. The figures are chosen with the options or with
the section ``graphics`` of the configuration file.

Examples
--------
.. code-block:: shell

    slab-yarray beamforming --input array.json --config capon.json \\
        -o result.jsonl --plot-dispersion --save-figures figures
    slab-yarray beamforming --data "data/*.mseed" --metadata array.xml \\
        --select "XX.*..HH?" --config rtbf.json -o rtbf.jsonl
    slab-yarray beamforming --result rtbf.jsonl --plot-spectrum \\
        --frequency 4.5 --show

    slab-yarray fdd --input array.json --config fdd.json -o fdd.jsonl \\
        --plot-singular-values --plot-mode-shape --show
    slab-yarray fdd --data "building/*.mseed" --metadata building.xml \\
        --fmin 0.5 --fmax 15 --window-type fixed --window-length 20 \\
        --components east --efdd --modes-csv modes.csv --plot-efdd \\
        --save-figures figures
    slab-yarray fdd --result fdd.jsonl --peak-range 0.2 0.8 --plot-mac

    slab-yarray response --input array.json --plot-response --show
    slab-yarray response --coordinates stations.txt --fmin 1 --fmax 20 \\
        --search-max 0.2 --plot-all --save-figures figures

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

from shakelab.exploration.yarray.beamforming import (
    Beamformer,
    BeamformingResult,
)
from shakelab.exploration.yarray.config import BeamformingConfig, FDDConfig
from shakelab.exploration.yarray.data import ArrayData
from shakelab.exploration.yarray.fdd import FDD, FDDResult, modes_to_csv
from shakelab.libutils.utils import command_message

__all__ = ["main"]

# Processing options of the beamforming: names of the BeamformingConfig
# parameters.
_BEAMFORMING_OPTIONS = (
    "frequencies", "fmin", "fmax", "nfreq", "frequency_scale", "method",
    "components", "diagonal_loading", "music_sources", "grid_quantity",
    "grid_scale", "grid_min", "grid_max", "grid_size", "azimuth_size",
    "window_type", "window_length", "window_taper", "window_overlap",
    "windows_per_block", "block_overlap", "ellipticity_angle_min",
    "ellipticity_angle_max", "ellipticity_angle_size", "spectra",
    "maxima", "minimum_relative_power", "statistical_uncertainty",
    "geometric_uncertainty", "strict", "transform", "fft_selection",
    "fft_padding",
)

# Processing options of the FDD: names of the FDDConfig parameters.
_FDD_OPTIONS = (
    "frequencies", "fmin", "fmax", "nfreq", "frequency_scale",
    "window_type", "window_length", "window_taper", "window_overlap",
    "windows_per_block", "block_overlap", "subfrequencies", "bandwidth",
    "normalization", "components", "azimuth", "vectors",
)

# Figures of the response analysis.
_RESPONSE_FIGURES = ("geometry", "response", "resolution", "velocity",
                     "search_domain")


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------


def main(argv: Optional[list] = None, prog: Optional[str] = None) -> int:
    """
    Run the command-line array processing.

    Parameters
    ----------
    argv : list of str, optional
        Arguments, starting with the subcommand (``beamforming``,
        ``fdd`` or ``response``); ``sys.argv[1:]`` if omitted.
    prog : str, optional
        Program name shown in the help.

    Returns
    -------
    int
        Exit status: 0 on success, 1 on errors, 130 when interrupted.
    """
    argv = sys.argv[1:] if argv is None else list(argv)
    args = _parser(prog).parse_args(argv)
    run = {
        "beamforming": _run_beamforming,
        "fdd": _run_fdd,
        "response": _run_response,
    }[args.command]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            warnings.showwarning = _show_warning
            return run(args)
    except KeyboardInterrupt:
        return 130
    except (OSError, ValueError, TypeError, np.linalg.LinAlgError) as error:
        print(f"error: {command_message(error)}", file=sys.stderr)
        return 1


# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------


def _parser(prog: Optional[str]) -> argparse.ArgumentParser:
    """Build the argument parser with its subcommands."""
    parser = argparse.ArgumentParser(
        prog=prog or "slab-yarray",
        description="Seismic-array processing: beamforming, frequency "
        "domain decomposition and array response.",
    )
    commands = parser.add_subparsers(dest="command", required=True,
                                     metavar="command")
    _beamforming_parser(commands)
    _fdd_parser(commands)
    _response_parser(commands)
    return parser


def _add_input(parser: argparse.ArgumentParser, flag: str, kind: Any,
               text: str) -> None:
    """Add the sources (one of --input, --data, flag) and data options."""
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path,
                        help="array file (JSON) with coordinates and files")
    source.add_argument("--data", nargs="+",
                        help="waveform files, patterns or folders")
    source.add_argument(flag, type=kind, help=text)

    data = parser.add_argument_group("data (with --data)")
    data.add_argument("--metadata", nargs="+",
                      help="metadata files (positions, orientations)")
    data.add_argument("--select", help="channel pattern NET.STA.LOC.CHA")
    data.add_argument("--start", help="start time (ISO 8601)")
    data.add_argument("--end", help="end time (ISO 8601)")


def _add_output(parser: argparse.ArgumentParser, arrays: str,
                table: str) -> Any:
    """Add the configuration and output options; return the group."""
    output = parser.add_argument_group("configuration and output")
    output.add_argument("--config", type=Path,
                        help="configuration file (JSON)")
    output.add_argument("--write-config", type=Path,
                        help="write the configuration used (JSON)")
    output.add_argument("-o", "--output", type=Path,
                        help="result file (JSON Lines, .jsonl)")
    output.add_argument("--append", action="store_true",
                        help="append the result to the file")
    output.add_argument("--no-arrays", action="store_true",
                        help=f"do not write the {arrays} (.npz) of the "
                        "result")
    output.add_argument("--csv", type=Path,
                        help=f"write the table of {table} (CSV), also of "
                        "a --result")
    output.add_argument("--owrite", action="store_true",
                        help="overwrite existing files")
    return output


def _add_frequencies_and_windows(group: Any, nfreq: str) -> None:
    """Add the options of frequencies and windows."""
    group.add_argument("--frequencies", type=float, nargs="+",
                       help="explicit frequencies (Hz)")
    group.add_argument("--fmin", type=float, help="minimum frequency")
    group.add_argument("--fmax", type=float, help="maximum frequency")
    group.add_argument("--nfreq", type=int, help=nfreq)
    group.add_argument("--frequency-scale", choices=("log", "linear"))
    group.add_argument("--window-type", choices=("cycles", "fixed"))
    group.add_argument("--window-length", type=float,
                       help="cycles, or seconds with fixed windows")
    group.add_argument("--window-taper", type=float)
    group.add_argument("--window-overlap", type=float)
    group.add_argument("--windows-per-block", type=int)
    group.add_argument("--block-overlap", type=float)


def _add_figure_files(group: Any) -> None:
    """Add the options to show and save the figures."""
    group.add_argument("--show", action=argparse.BooleanOptionalAction,
                       help="show the figures")
    group.add_argument("--save-figures", type=Path,
                       help="folder of the figures")
    group.add_argument("--figure-format", help="image format (png)")
    group.add_argument("--dpi", type=int, help="resolution (150)")


def _beamforming_parser(commands: Any) -> None:
    """Add the subcommand of the beamforming."""
    parser = commands.add_parser(
        "beamforming", help="frequency-wavenumber processing",
        description="Frequency-wavenumber processing of seismic arrays.",
    )
    _add_input(parser, "--result", Path,
               "saved result to plot (no processing)")
    _add_output(parser, "spectra", "the maxima")

    processing = parser.add_argument_group(
        "processing (replace the values of --config)"
    )
    _add_frequencies_and_windows(
        processing,
        "number of frequencies (omitted with --fmin/--fmax and "
        "--transform rfft: every FFT bin)",
    )
    processing.add_argument("--method",
                            choices=("bartlett", "capon", "music"))
    processing.add_argument("--components", choices=(
        "vertical", "radial", "transverse", "three_component_separate",
        "rayleigh_joint", "rayleigh_angular", "rayleigh_direct",
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
    processing.add_argument("--ellipticity-angle-min", type=float)
    processing.add_argument("--ellipticity-angle-max", type=float)
    processing.add_argument("--ellipticity-angle-size", type=int)
    processing.add_argument("--spectra",
                            action=argparse.BooleanOptionalAction,
                            help="keep the complete spectra")
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
        choices=("log_slowness", "log_velocity", "wavenumber",
                 "ellipticity_angle"),
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
    _add_figure_files(figures)


def _fdd_parser(commands: Any) -> None:
    """Add the subcommand of the frequency domain decomposition."""
    parser = commands.add_parser(
        "fdd", help="frequency domain decomposition",
        description="Frequency domain decomposition of seismic arrays: "
        "singular values, peaks and mode shapes; enhanced FDD for "
        "natural frequencies and damping ratios.",
    )
    _add_input(parser, "--result", Path,
               "saved result to analyse and plot (no processing)")
    output = _add_output(parser, "singular values and vectors",
                         "the first singular vectors")
    output.add_argument("--modes-csv", type=Path,
                        help="write the table of the modes (CSV)")

    processing = parser.add_argument_group(
        "processing (replace the values of --config)"
    )
    _add_frequencies_and_windows(
        processing,
        "number of frequencies (omitted with --fmin/--fmax and fixed "
        "windows: every FFT frequency)",
    )
    processing.add_argument("--subfrequencies", type=int,
                            help="sub-frequencies averaged (1)")
    processing.add_argument("--bandwidth", type=float,
                            help="relative band of the sub-frequencies")
    processing.add_argument("--normalization",
                            choices=("none", "coherence"))
    processing.add_argument("--components", choices=(
        "vertical", "north", "east", "radial", "transverse", "horizontal",
        "three_component",
    ))
    processing.add_argument("--azimuth", type=float,
                            help="azimuth of radial and transverse "
                            "(degrees from north)")
    processing.add_argument("--vectors", type=int,
                            help="singular vectors kept per block (3)")

    modes = parser.add_argument_group("peaks and modes")
    modes.add_argument("--peak-range", type=float, nargs=2,
                       metavar=("FMIN", "FMAX"),
                       help="band of the peaks (Hz; default: all)")
    modes.add_argument("--prominence", type=float,
                       help="smallest prominence of a peak (2)")
    modes.add_argument("--peak-count", type=int,
                       help="largest number of peaks (all)")
    modes.add_argument("--modes", type=float, nargs="+",
                       help="frequencies of the modes (Hz), instead of "
                       "the peaks")
    modes.add_argument("--efdd", action="store_true",
                       help="enhanced FDD of the modes (needs data)")
    modes.add_argument("--mac", type=float,
                       help="smallest MAC of the bell (0.8)")
    modes.add_argument("--cycles", type=float,
                       help="windows of the EFDD in cycles of the mode "
                       "(100)")
    modes.add_argument("--duration", type=float,
                       help="windows of the EFDD in seconds")
    modes.add_argument("--levels", type=float, nargs=2,
                       metavar=("LOW", "HIGH"),
                       help="amplitudes of the extrema of the decay "
                       "(0.1 0.9)")

    figures = parser.add_argument_group(
        "figures (replace the section 'graphics' of --config)"
    )
    figures.add_argument("--plot-singular-values",
                         action=argparse.BooleanOptionalAction)
    figures.add_argument("--count", type=int,
                         help="singular values drawn (3)")
    figures.add_argument("--plot-mode-shape",
                         action=argparse.BooleanOptionalAction,
                         help="shape of every mode")
    figures.add_argument("--along",
                         choices=("profile", "x", "y", "z", "index"),
                         help="abscissa of the shapes (profile, or index "
                         "without positions)")
    figures.add_argument("--plot-mac", action=argparse.BooleanOptionalAction,
                         help="MAC between the modes")
    figures.add_argument("--plot-efdd", action=argparse.BooleanOptionalAction,
                         help="bell and decay of every mode (with --efdd)")
    _add_figure_files(figures)


def _response_parser(commands: Any) -> None:
    """Add the subcommand of the response analysis."""
    parser = commands.add_parser(
        "response", help="geometric response and limits",
        description="Geometric response and limits of a seismic array.",
    )
    _add_input(
        parser, "--coordinates", Path,
        "text file of station coordinates: lines 'x y [z]' or "
        "'name x y [z]' (comma or space separated, '#' comments)",
    )

    analysis = parser.add_argument_group("analysis")
    analysis.add_argument("--resolution-level", type=float, default=0.5,
                          help="relative power of the main lobe (0.5)")
    analysis.add_argument("--ambiguity-level", type=float, default=0.5,
                          help="relative power of the ambiguities (0.5)")
    analysis.add_argument("--expected", action="store_true",
                          help="response expected under the coordinate "
                          "uncertainties of the array")
    analysis.add_argument("--grid-max", type=float,
                          help="extent of the wavenumber grid (automatic)")
    analysis.add_argument("--grid-size", type=int,
                          help="odd samples per axis of the grid "
                          "(automatic)")
    analysis.add_argument("--frequencies", type=float, nargs="+",
                          help="frequencies of the velocity limits (Hz)")
    analysis.add_argument("--fmin", type=float, help="minimum frequency")
    analysis.add_argument("--fmax", type=float, help="maximum frequency")
    analysis.add_argument("--nfreq", type=int, default=30,
                          help="number of log-spaced frequencies (30)")
    analysis.add_argument("--search-max", type=float,
                          help="radius of a circular search domain to check")

    figures = parser.add_argument_group("figures")

    for name in _RESPONSE_FIGURES:
        flag = name.replace("_", "-")
        figures.add_argument(f"--plot-{flag}", action="store_true",
                             help=f"{name.replace('_', ' ')} figure")

    figures.add_argument("--plot-all", action="store_true",
                         help="all the figures available")
    _add_figure_files(figures)


# ----------------------------------------------------------------------------
# Data and configuration
# ----------------------------------------------------------------------------


def _configuration(kind: type, names: tuple,
                   args: argparse.Namespace, required: str) -> Any:
    """Configuration from the file, with the options replacing values."""
    config = kind(
        args.config,
        **{name: getattr(args, name) for name in names
           if getattr(args, name) is not None},
    )

    if not config.complete:
        raise ValueError(f"the configuration needs {required}")

    return config


def _array(args: argparse.Namespace) -> ArrayData:
    """Array data from an array file or waveform files."""
    if args.input is not None:
        return ArrayData(args.input)

    from shakelab.signals.io import metadata_reader, reader

    metadata = None if args.metadata is None else metadata_reader(
        args.metadata
    )
    collection = reader(
        args.data, metadata=metadata, select=args.select,
        starttime=args.start, endtime=args.end,
    )
    return ArrayData(collection, select=args.select)


def _check_result_only(args: argparse.Namespace) -> None:
    """Refuse the options that need data with --result."""
    for name in ("config", "output", "write_config"):
        if getattr(args, name) is not None:
            raise ValueError(f"--{name.replace('_', '-')} needs data")


# ----------------------------------------------------------------------------
# Beamforming
# ----------------------------------------------------------------------------


def _run_beamforming(args: argparse.Namespace) -> int:
    """Process or load a result, save it and draw the figures."""
    graphics = _beamforming_graphics(args)

    if args.result is not None:
        _check_result_only(args)
        result = BeamformingResult(args.result)
    else:
        config = _configuration(
            BeamformingConfig, _BEAMFORMING_OPTIONS, args,
            "frequencies and grid: give --config or --frequencies/--fmin "
            "--fmax and --grid-min --grid-max --grid-size",
        )

        if args.write_config is not None:
            config.write(args.write_config, owrite=args.owrite)

        array = _array(args)
        print(array.summary())
        result = Beamformer(array, config).process()

        if args.output is not None:
            result.write(args.output, append=args.append,
                         owrite=args.owrite, arrays=not args.no_arrays)
            print(f"result written to {args.output}")

    if args.csv is not None:
        result.to_csv(args.csv, owrite=args.owrite)
        print(f"maxima written to {args.csv}")

    print(result.summary())

    for path in _beamforming_figures(result, graphics):
        print(f"figure written to {path}")

    return 0


def _beamforming_graphics(args: argparse.Namespace) -> dict:
    """Figure options: section 'graphics' of the config, then options."""
    section = _graphics_section(args.config)
    spectrum = _toggle(section.get("spectrum", {}))
    density = _toggle(section.get("density", {}))

    return {
        "dispersion": _pick(args.plot_dispersion,
                            bool(section.get("dispersion", False))),
        "ellipticity": _pick(args.plot_ellipticity,
                             bool(section.get("ellipticity", False))),
        "spectrum": _pick(args.plot_spectrum,
                          bool(spectrum.get("enabled", False))),
        "density": _pick(args.plot_density,
                         bool(density.get("enabled", False))),
        "coordinate": _pick(args.density_coordinate,
                            density.get("coordinate", "log_slowness")),
        "component": _pick(args.component, section.get("component")),
        "ranks": _pick(args.ranks, section.get("ranks")),
        "frequency": _pick(args.frequency, spectrum.get("frequency")),
        "block": _pick(args.block, spectrum.get("block", 0)),
        "view": _pick(args.spectrum_view, spectrum.get("view", "kxky")),
        "quantity": _pick(args.spectrum_quantity,
                          spectrum.get("quantity", "power")),
        **_saving(args, section),
    }


def _beamforming_figures(result: Any, options: dict) -> list:
    """Draw the requested figures; return the files written."""
    wanted = [
        name for name in ("dispersion", "ellipticity", "spectrum", "density")
        if options[name]
    ]

    if not wanted:
        return []

    if "spectrum" in wanted and options["frequency"] is None:
        raise ValueError("the spectrum figure needs --frequency")

    _check_saving(options)

    from shakelab.graphics import yarray as graphics

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
        density = result.density(
            component=component, ranks=options["ranks"],
            coordinate=options["coordinate"],
        )
        figure, _ = graphics.plot_density(
            density, peaks=density.peaks(3, minimum_relative_density=0.2)
        )
        figures.append((figure, f"density_{density.component}"))

    return _output_figures(figures, options)


# ----------------------------------------------------------------------------
# Frequency domain decomposition
# ----------------------------------------------------------------------------


def _run_fdd(args: argparse.Namespace) -> int:
    """Process or load an FDD result, find its modes, draw the figures."""
    graphics = _fdd_graphics(args)
    fdd = None

    if args.result is not None:
        _check_result_only(args)

        if args.efdd:
            raise ValueError("--efdd needs data")

        result = FDDResult(args.result)
    else:
        config = _configuration(
            FDDConfig, _FDD_OPTIONS, args,
            "frequencies: give --config or --frequencies/--fmin --fmax",
        )

        if args.write_config is not None:
            config.write(args.write_config, owrite=args.owrite)

        array = _array(args)
        print(array.summary())
        fdd = FDD(array, config)
        result = fdd.process()

        if args.output is not None:
            result.write(args.output, append=args.append,
                         owrite=args.owrite, arrays=not args.no_arrays)
            print(f"result written to {args.output}")

    if args.csv is not None:
        result.to_csv(args.csv, owrite=args.owrite)
        print(f"singular vectors written to {args.csv}")

    print(result.summary())
    modes = _fdd_modes(args, result, fdd)

    for mode in modes:
        print(mode.summary())

    if args.modes_csv is not None:
        modes_to_csv(args.modes_csv, modes, owrite=args.owrite)
        print(f"modes written to {args.modes_csv}")

    for path in _fdd_figures(result, modes, graphics):
        print(f"figure written to {path}")

    return 0


def _fdd_modes(args: argparse.Namespace, result: FDDResult,
               fdd: Optional[FDD]) -> list:
    """Modes at the given frequencies or at the peaks (FDD or EFDD)."""
    if args.modes is not None:
        frequencies = sorted(args.modes)
    else:
        fmin, fmax = args.peak_range or (None, None)
        options = {"fmin": fmin, "fmax": fmax, "count": args.peak_count}

        if args.prominence is not None:
            options["prominence"] = args.prominence

        frequencies = [float(value) for value in result.peaks(**options)]

    if not frequencies:
        print("no peaks: no modes")
        return []

    if not args.efdd:
        return [result.mode(value) for value in frequencies]

    options = {
        name: getattr(args, name) for name in ("mac", "cycles", "duration")
        if getattr(args, name) is not None
    }

    if args.levels is not None:
        options["levels"] = tuple(args.levels)

    return fdd.efdd(frequencies, **options)


def _fdd_graphics(args: argparse.Namespace) -> dict:
    """Figure options: section 'graphics' of the config, then options."""
    section = _graphics_section(args.config)
    values = _toggle(section.get("singular_values", {}))
    shape = _toggle(section.get("mode_shape", {}))

    return {
        "singular_values": _pick(args.plot_singular_values,
                                 bool(values.get("enabled", False))),
        "count": _pick(args.count, values.get("count", 3)),
        "mode_shape": _pick(args.plot_mode_shape,
                            bool(shape.get("enabled", False))),
        "along": _pick(args.along, shape.get("along")),
        "mac": _pick(args.plot_mac, bool(section.get("mac", False))),
        "efdd": _pick(args.plot_efdd, bool(section.get("efdd", False))),
        **_saving(args, section),
    }


def _fdd_figures(result: FDDResult, modes: list, options: dict) -> list:
    """Draw the requested figures of the FDD; return the files written."""
    wanted = [
        name for name in ("singular_values", "mode_shape", "mac", "efdd")
        if options[name]
    ]

    if not wanted:
        return []

    if set(wanted) & {"mode_shape", "mac", "efdd"} and not modes:
        raise ValueError("the figures of the modes need peaks or --modes")

    if "efdd" in wanted and modes[0].efdd is None:
        raise ValueError("the EFDD figure needs --efdd")

    _check_saving(options)

    from shakelab.graphics import yarray as graphics

    figures = []

    if "singular_values" in wanted:
        figure, _ = graphics.plot_singular_values(result,
                                                  count=options["count"])
        figures.append((figure, "singular_values"))

    if "mode_shape" in wanted:
        for mode in modes:
            figure, _ = graphics.plot_mode_shape(mode,
                                                 along=options["along"])
            figures.append((figure, f"mode_shape_f{mode.frequency:.4g}"))

    if "mac" in wanted:
        figures.append((graphics.plot_mac(modes)[0], "mac"))

    if "efdd" in wanted:
        for mode in modes:
            figure, _ = graphics.plot_efdd(mode)
            figures.append((figure, f"efdd_f{mode.frequency:.4g}"))

    return _output_figures(figures, options)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def _graphics_section(path: Optional[Path]) -> dict:
    """Section 'graphics' of a configuration file (empty if none)."""
    if path is None or not path.is_file():
        return {}

    data = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(data, Mapping):
        return {}

    return dict(data.get("graphics", {}))


def _toggle(value: Any) -> dict:
    """A figure entry of the section: a flag or a dictionary."""
    return {"enabled": value} if isinstance(value, bool) else dict(value)


def _pick(value: Any, default: Any) -> Any:
    """The command-line value, or the default when it is not given."""
    return default if value is None else value


def _saving(args: argparse.Namespace, section: Mapping) -> dict:
    """Options to show and save the figures."""
    save = _toggle(section.get("save", {}))
    directory = _pick(args.save_figures, save.get("directory"))

    return {
        "show": _pick(args.show, bool(section.get("show", False))),
        "directory": None if directory is None else Path(directory),
        "save": args.save_figures is not None or bool(
            save.get("enabled", directory is not None)
        ),
        "format": _pick(args.figure_format, save.get("format", "png")),
        "dpi": _pick(args.dpi, save.get("dpi", 150)),
    }


def _check_saving(options: dict) -> None:
    """Refuse saving without a folder."""
    if options["save"] and options["directory"] is None:
        raise ValueError("saving figures needs --save-figures DIR")


def _output_figures(figures: list, options: dict) -> list:
    """Save, show and close figures ``(figure, stem)``; return the files."""
    from shakelab.graphics import yarray as graphics

    plt = graphics._pyplot()
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


# ----------------------------------------------------------------------------
# Array response
# ----------------------------------------------------------------------------


def _read_coordinates(path: Path) -> tuple:
    """Return names (or None) and coordinates of a text file."""
    names, rows = [], []

    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("#", 1)[0].replace(",", " ").split()

        if not fields:
            continue

        try:
            float(fields[0])
        except ValueError:
            names.append(fields[0])
            fields = fields[1:]

        if len(fields) not in (2, 3):
            raise ValueError(f"{path}: invalid line {line.strip()!r}")

        rows.append([float(value) for value in fields] + [0.0] * (
            3 - len(fields)
        ))

    if names and len(names) != len(rows):
        raise ValueError(f"{path}: names on some lines only")

    return (names or None), np.array(rows)


def _response_frequencies(args: argparse.Namespace) -> Optional[np.ndarray]:
    """Frequencies of the velocity limits, if requested."""
    if args.frequencies is not None:
        return np.array(sorted(args.frequencies))

    if args.fmin is None and args.fmax is None:
        return None

    if args.fmin is None or args.fmax is None:
        raise ValueError("--fmin and --fmax go together")

    return np.geomspace(args.fmin, args.fmax, args.nfreq)


def _run_response(args: argparse.Namespace) -> int:
    """Analyse the response of an array and draw the figures."""
    from shakelab.exploration.yarray.response import ArrayResponse

    names = None

    if args.coordinates is not None:
        names, array = _read_coordinates(args.coordinates)
        source = array
    else:
        source = _array(args)
        names = list(source.sids)
        print(source.summary())

    response = ArrayResponse(source, grid_max=args.grid_max,
                             grid_size=args.grid_size)
    response.info()
    analysis = response.analyze(args.resolution_level,
                                args.ambiguity_level,
                                expected=args.expected)
    analysis.info()
    frequencies = _response_frequencies(args)
    limits = None

    if frequencies is not None:
        if analysis.limits.aliasing_radius is None:
            print("no upper wavenumber limit: velocity limits not "
                  "available")
        else:
            limits = analysis.velocity_limits(frequencies)
            print("frequency (Hz)  minimum velocity  maximum velocity")

            for f, low, high in zip(limits.frequencies, limits.minimum,
                                    limits.maximum):
                print(f"{f:14.4g}  {low:16.6g}  {high:16.6g}")

    domain = None

    if args.search_max is not None:
        domain = response.search_domain(args.search_max,
                                        args.ambiguity_level,
                                        expected=args.expected)
        domain.info()

    wanted = [
        name for name in _RESPONSE_FIGURES
        if args.plot_all or getattr(args, f"plot_{name}")
    ]

    if args.plot_all:
        wanted = [
            name for name in wanted
            if (name != "velocity" or limits is not None)
            and (name != "search_domain" or domain is not None)
        ]
    elif "velocity" in wanted and limits is None:
        raise ValueError("the velocity figure needs frequencies and an "
                         "upper wavenumber limit")
    elif "search_domain" in wanted and domain is None:
        raise ValueError("the search-domain figure needs --search-max")

    for path in _response_figures(args, wanted, analysis, names, limits,
                                  domain):
        print(f"figure written to {path}")

    return 0


def _response_figures(args: argparse.Namespace, wanted: list,
                      analysis: Any, names: Any, limits: Any,
                      domain: Any) -> list:
    """Draw the figures of the response; return the files written."""
    if not wanted:
        return []

    options = _saving(args, {})
    _check_saving(options)

    from shakelab.graphics import yarray as graphics

    figures = []

    if "geometry" in wanted:
        figure, _ = graphics.plot_array_geometry(
            analysis.response, labels=names
        )
        figures.append((figure, "array_geometry"))

    if "response" in wanted:
        figures.append((graphics.plot_array_response(analysis)[0],
                        "array_response"))

    if "resolution" in wanted:
        figures.append((graphics.plot_resolution(analysis, polar=True)[0],
                        "resolution"))

    if "velocity" in wanted:
        figures.append((graphics.plot_velocity_limits(limits)[0],
                        "velocity_limits"))

    if "search_domain" in wanted:
        figures.append((graphics.plot_search_domain(analysis, domain)[0],
                        "search_domain"))

    return _output_figures(figures, options)


# ----------------------------------------------------------------------------
# Warnings
# ----------------------------------------------------------------------------


def _show_warning(message: Any, category: Any, *args: Any,
                  **kwargs: Any) -> None:
    """Print warnings as one line."""
    print(f"warning: {message}", file=sys.stderr)
