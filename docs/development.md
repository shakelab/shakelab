# ShakeLab development and architecture guidelines

*Version 2.10, 18 November 2026*

[ShakeLab documentation](README.md) › Development guidelines

This document defines the design philosophy of ShakeLab and the
principles that all contributions follow. Parts I to III state
principles. The concrete scientific and technical conventions (time,
coordinates, units, numerical representation, supported versions) are
normative and are listed in [conventions.md](conventions.md).

The guidelines apply to human developers and to AI assistants working
on the code base alike.

## Contents

- [Part I. Architecture](#part-i-architecture)
- [Part II. Scientific semantics](#part-ii-scientific-semantics)
- [Part III. Engineering practice](#part-iii-engineering-practice)

---

## Part I. Architecture

### 1. Project identity

ShakeLab is a modular scientific Python framework for seismology,
engineering seismology, monitoring, geophysical data processing and
related scientific and operational applications. It is expected to
expand into further scientific and technical domains.

The current modules, classes, formats and applications do not define
the boundaries of ShakeLab. They are implementations of the
architecture, not the architecture itself. ShakeLab must evolve as a
coherent scientific framework, not as a collection of scripts or
project-specific solutions.

Its design objectives are:

- scientific correctness and transparency;
- modularity, composability and reusability;
- clear separation of responsibilities;
- intentionally designed, reasonably stable public APIs;
- limited and justified dependencies;
- computational efficiency;
- applicability to both research and operational workflows.

ShakeLab is a library of functions for seismologists, geophysicists
and engineers. A feature should not merely solve its immediate use
case: whenever reasonable, it should become a reusable scientific
capability of the library. Modules that are not used by other parts of
the code base are still library functionality, and are maintained and
improved as such.

### 2. Scientific software first

Scientific meaning takes precedence over programming convenience.
Abstractions and APIs should reflect scientific concepts, not the
incidental structure of input files, user interfaces or individual
research workflows.

Important scientific assumptions must remain visible in the code and
in the API. Do not hide scientifically relevant choices behind
implementation details or unexplained defaults.

Scientific proposals are evaluated critically before implementation.
Check, as relevant: dimensional and mathematical consistency, physical
interpretation, limiting behavior, numerical stability, sensitivity to
sampling, validity range and compatibility with established
definitions. If a formulation is questionable, resolve the scientific
issue before embedding it in the library. Scientific coding is not
mere translation of an idea into Python.

### 3. Layers and dependency direction

Organize functionality by conceptual responsibility. A recurring
layering, from the top down, is:

```
user-facing interfaces and applications
        |
orchestration
        |
scientific algorithms
        |
domain abstractions
        |
low-level operations
        |
fundamental representations (primitives)
```

Dependencies flow downward: each layer may use the layers below it,
never those above. A lower-level component must remain usable without
importing the components that consume it. The scientific core must not
depend on any application, user interface, service or deployment.

Not every subsystem needs every layer. Do not create artificial layers
to reproduce the diagram; place each responsibility that genuinely
exists at the appropriate level.

Avoid circular dependencies. When one appears, first check whether a
responsibility sits at the wrong level. Local or dynamic imports are
acceptable for genuinely optional or deferred dependencies (for
example a convenience method on a primitive that delegates to a higher
module), never to conceal a wrong structure.

### 4. One canonical implementation

This is the fundamental principle of ShakeLab. Each scientific
operation has one authoritative implementation, placed at the lowest
level at which its meaning is complete and reusable. Higher levels
expose convenience interfaces, coordinate or apply the operation to
collections, and delegate to the canonical implementation:

```
canonical implementation         e.g. great_circle_distance()
        |
domain-object interface          e.g. WgsPoint.epicentral_distance_to()
        |
collection-level delegation
        |
high-level workflow, CLI, service, GUI
```

When a concept exists for individual objects and for collections, the
operation belongs to the element. Collections coordinate and delegate,
and add only behavior that genuinely belongs to the collection:
selection, grouping, synchronization, aggregation and cross-element
validation.

Scientific functionality first exists as a reusable Python API.
Scripts, command-line tools, notebooks, services, web and graphical
applications consume that API; they never reimplement it. One
scientific core, multiple access mechanisms.

### 5. Representation, computation and orchestration

Keep three concerns distinguishable. *Representation* defines what the
data or scientific objects are. *Computation* implements
transformations, measurements, models and algorithms. *Orchestration*
decides which components are used and in which order.

Do not embed workflow decisions inside scientific algorithms, and do
not duplicate scientific calculations inside orchestration code.

More generally, each class, function and module has one clear
responsibility: parsers parse, converters convert, calculators
calculate, aggregators aggregate, renderers render, clients
communicate, orchestrators orchestrate. These are not naming rules;
the point is that distinct responsibilities remain separable. A class
that manages representation, computation, configuration, persistence
and visualization at once is crossing boundaries.

### 6. Domain objects and data representation

Represent meaningful scientific concepts explicitly rather than
propagating loose dictionaries, tuples or unrelated arrays. A domain
object represents a coherent concept, maintains its invariants,
exposes the operations that naturally belong to it and separates
internal representation from external use.

Domain objects are not universal service containers. An object
representing scientific data does not need to read every format,
implement every analysis, render itself, manage configuration or
orchestrate workflows. Place such operations in the appropriate
conversion, scientific, I/O or orchestration layer.

Conversely, do not introduce custom classes when a scalar, a NumPy
array or a standard collection already expresses the concept clearly.
Use NumPy arrays where array semantics are natural, preserve dtype
intentionally and state expected shapes where ambiguity could produce
wrong results.

### 7. Generality, composition and extensibility

Functionality should be reusable beyond the dataset or project that
motivated it. Avoid assumptions about one experiment, station,
network, region, dataset, file location, deployment or workflow,
unless they are intrinsic to the scientific concept.

Avoid premature abstraction as well. Generalize when there is a
genuine conceptual commonality, not because two pieces of code look
similar. Prefer small abstractions that can evolve over frameworks
designed for hypothetical requirements. The target is simple, general,
explicit and composable, not maximally abstract.

Build complex functionality by composing small operations with clear
semantics, each independently useful, rather than monolithic functions
that mix loading, preprocessing, calculation, visualization and
export. High-level convenience workflows compose canonical operations;
they are not alternative implementations.

When several implementations play the same conceptual role (formats,
ground-motion models, providers), define a small contract around that
role, so that consumers do not grow chains of implementation-specific
conditionals. Choose the simplest mechanism that expresses the
contract: a documented convention, duck typing, a `typing.Protocol`, a
base class or composition. Do not create abstract base classes
preemptively.

Before adding functionality, classify it as general scientific,
domain-specific but reusable, application-specific, deployment
configuration or one-off analysis, and place it accordingly. When an
application reveals a missing general capability, extract the
scientific concept rather than moving the application workflow into
the core.

New scientific domains are not forced into existing abstractions
because those abstractions exist. Decide whether the new functionality
uses existing primitives, extends an existing domain, belongs as a
sibling domain or needs a new abstraction. Reuse a concept when its
semantics match, not when its data structure looks similar. The
architecture grows by adding coherent components, not by overloading
existing ones.

---

## Part II. Scientific semantics

### 8. Scientific explicitness and metadata

Make scientifically relevant conventions explicit: units, coordinate
systems, reference frames, time references, sign conventions,
component orientation, normalization, sampling conventions, model
definitions and reference quantities. Where several conventions exist
in practice, state which one ShakeLab implements and apply it
consistently. The conventions adopted library-wide are listed in
[conventions.md](conventions.md); a module that needs a different one
must declare it.

Metadata read from a source may provide defaults, but explicit user
input takes precedence over metadata-derived or inferred values, and
the precedence must be predictable. Distinguish mandatory, optional,
inferred and explicitly provided values. Do not silently override
explicit parameters with metadata, and do not invent missing mandatory
information. Inference is acceptable only when it is scientifically
well defined and documented.

### 9. Units and physical quantities

Physical quantities keep clear semantics. Units are stated in the API
and documentation, never left to be inferred from variable names.
ShakeLab uses SI units internally, with a short list of declared
exceptions where established practice requires it
([conventions, Units](conventions.md#3-units)).

Conversions happen at identifiable boundaries (I/O, a named conversion
function, a clearly named property such as `depth_km`), never repeated
or hidden across several layers. Scientifically different quantities
remain distinct even when their arithmetic is similar: organize code
by physical meaning, not by similarity of formulas.

### 10. Time semantics

Time is scientific data, not a generic scalar. Distinguish absolute
time, relative time, duration, sampling interval and sample index, and
make conversions between them explicit. Absolute times have a
well-defined reference system and deterministic semantics; ambiguous
timestamps are not accepted. Operations that change temporal extent,
sampling or reference time update or preserve the related metadata
consistently. The concrete time conventions are in
[conventions, Time](conventions.md#1-time).

### 11. Spatial and coordinate semantics

Coordinates are scientific data, not generic vectors. Distinguish
geographic, projected and local metric coordinates, distances,
directions, orientations and reference frames. Coordinate
transformations are explicit and centralized in the geodetic module.
Quantities in incompatible reference systems are never combined
silently. The concrete conventions, including the order of coordinates
and the vertical reference, are in
[conventions, Space and coordinates](conventions.md#2-space-and-coordinates).

### 12. Scientific algorithms and numerics

When implementing a scientific method: understand its mathematical
formulation, assumptions and validity domain; define its inputs and
outputs; map mathematical quantities clearly to code; keep it separate
from I/O and visualization; handle meaningful boundary conditions. Do
not alter an established formulation because another form is easier to
implement. If ShakeLab deliberately implements a particular variant,
say so and cite the reference.

Do not assume that mathematically equivalent formulations are
numerically equivalent. Pay attention to precision, conditioning,
overflow and underflow, boundary effects, missing values, array shapes
and unnecessary copies. Vectorize with NumPy when it improves clarity
and efficiency, not when it obscures the algorithm without a
measurable benefit.

### 13. Reproducibility

Results must be reproducible. Stochastic methods (sampling,
simulations, Monte Carlo, random meshes) take an explicit seed,
through a parameter named `seed` (an integer, `None` for a
non-reproducible run, or a `numpy.random.Generator`), and never use
the global NumPy or Python random state
([conventions, Random numbers](conventions.md#5-random-numbers-and-reproducibility)).
Library code does not depend on wall-clock time, environment variables
or the iteration order of unordered containers in ways that change
scientific results.

---

## Part III. Engineering practice

### 14. Input/output, formats and serialization

Scientific I/O is layered:

```
storage format
        |  parser / decoder
format-specific representation
        |  conversion
ShakeLab domain representation
        |
public I/O API (read / write)
```

Not every format needs every layer explicitly, but storage-format
details must not leak into scientific processing. Downstream
operations work on ShakeLab domain representations, independently of
the original format. Parsing, scientific interpretation and
serialization are separate responsibilities, and a domain object is
not designed around any one persistence format.

When implementing a format: keep binary layout logic close to the
specification and away from scientific processing; preserve
scientifically relevant metadata and numerical types; distinguish
malformed data from unsupported features, and missing optional
information from invalid mandatory data; avoid unnecessary copies; aim
at scientifically meaningful round trips.

Interoperability with external ecosystems (for example ObsPy) is
provided through optional adapters (`to_obspy`, `from_obspy` or
equivalent) placed in separate modules. The core never depends on
them.

### 15. Performance, memory and mutability

Correctness and clarity come first, but performance is a design
constraint. Optimize roughly in this order: algorithm; data
representation; avoiding unnecessary computation, copies, parsing and
conversion; efficient array operations and indexing; lower-level
optimization only when measured and justified. For large datasets
memory is part of performance: prefer subset or streaming processing
and avoid large intermediate arrays.

Mutation is deliberate and visible. Never modify input data
unexpectedly. When an operation conceptually produces a new scientific
object, returning a new object is usually clearer; for large arrays
in-place operation may be justified. The choice is made per operation,
documented in the docstring, and consistent across similar operations.

### 16. Validation, errors, warnings and logging

Validate each invariant at the layer that owns it: structure in the
low-level representation, concept invariants in the domain object,
method assumptions in the algorithm, relationships between components
in the workflow. Do not repeat the same check at every layer.

Errors are meaningful at the level where they occur. Distinguish
invalid scientific input, malformed data, unsupported features,
missing mandatory information, unavailable optional information,
invalid configuration, external service failures and internal
programming errors. Use the most specific built-in exception
(`ValueError`, `TypeError`, `NotImplementedError`, ...) or a ShakeLab
exception when callers need to distinguish the case. A ShakeLab
exception derives from the built-in exception it specializes (a file
without waveforms raises `NotWaveformError`, a `ValueError`), so that
code handling the general case keeps working.

Never catch exceptions broadly to let execution continue, and never
replace an invalid operation with a scientifically different fallback.
A fallback is acceptable only when it preserves the intended meaning
and is explicit.

Never return a wrong or incomplete result silently. Where the data do
not allow the intended operation (margins missing around a window,
channels without metadata, files without waveforms in a folder), the
operation either raises or does what it can and says what it left out.
Warnings raised deep in the library are attributed to the caller's
code (`shakelab.libutils.utils.warn_external`); one call gives one
warning that lists all the cases (one per channel in real-time
streams), not one per record.

Raise an exception when the result would be wrong. Issue a warning
(`warnings.warn`, with an appropriate category) when the result is
valid but the user should know something: an approximation, a
deprecated call, a value outside the calibrated range of a model. Use
the `logging` module, with per-module loggers, for operational
information. Library code does not use `print` for debugging or error
reporting; the only exceptions are the methods whose purpose is to
show a summary to the user (`info()`, see Section 17) and the
command-line modules (`cli`).

### 17. Public API

The public API is designed, not the sum of whatever exists in a
module. Keep it small, coherent, predictable, scientifically
meaningful and reasonably stable. A name is public if it is listed in
the module's `__all__`, or, where `__all__` is absent, if it does not
begin with an underscore. Every module should therefore define
`__all__`, and package-level imports should export only the intended
API. Private helpers use the `_` prefix.

Similar operations follow similar conventions (naming, argument order,
optional parameters, mutation, metadata precedence, return values,
exceptions) when their semantics are genuinely similar. Semantic
correctness takes precedence over superficial uniformity. The
conventions adopted so far are:

- reading and writing: `source` (path, pattern, content or open file)
  and `target` (file or folder) as argument names, `format` for the
  format, `owrite` to overwrite existing files (False by default); the
  data are the first argument of a writer; the arguments after the
  first (after `target` in writers) are keyword-only; an existing file
  without `owrite` gives the same `FileExistsError` everywhere
  (`shakelab.libutils.utils.check_target`: `<file> exists (use
  owrite=True)`), and the commands have `--owrite` and show it in the
  message;
- return values: methods that change an object in place (`read`,
  `append`, `cut`, `filter`...) return None, and raise when they cannot
  do what was asked rather than returning a flag (`Record.append`
  raises for records that cannot be joined); functions that download
  to a `target` instead of returning the data return the number of
  bytes written (0: no data, no file);
- time windows: `starttime` and `endtime`, `Date` objects or ISO 8601
  strings, bounds included, keyword-only; the time of the first sample
  of data that do not have one (text tables, TDMS) is
  `first_sample_time`, a different concept;
- servers: one address, read by `shakelab.libutils.address`
  (`host[:port]`, IPv6 in brackets, and for web services
  `[scheme://]host[:port][/prefix]` with `https` if the scheme is
  omitted, never falling back to plain HTTP); a `port` given separately
  takes precedence over the one in the address; the same form on the
  command line, with `--port`;
- metadata attached to data: `metadata=` accepts metadata at hand (a
  `MetadataCollection` or files) wherever data are read or received
  (`io`, `fdsnws`, `seedlink`), and `True` in the services to download
  them;
- channel selection: `select`, with the patterns of
  `shakelab.signals.codes` (`NET.STA.LOC.CHA`, `*`, `?`, `[...]`), a
  pattern or a list of patterns, with the same meaning while reading
  (`io.reader`) and on collections in memory (`select` methods, which
  return copies);
- containers: a container is a sequence of its elements (a record of
  samples, a stream of records, a collection of streams); an integer
  gives one element, a slice a new object of the same class, and `in`
  follows the indexing, with an operand of the wrong kind raising
  `TypeError` instead of answering False;
- options with a fixed set of values are enumerations, and their
  string values are accepted wherever the enumeration is;
- objects that can describe themselves have `summary()`, returning
  the text, and `info()`, printing it: a description of the object
  (for a client, its server, configuration and state), not of
  something else;
- an argument of the wrong type raises `TypeError`, a value out of its
  domain `ValueError`.

ShakeLab is a library used outside this repository. A public name is
never removed because it is unused inside the repository: it is fixed,
kept, or deprecated as described in the next section.

### 18. Compatibility, deprecation and transitional modules

Preserve established public behavior unless it obstructs a significant
scientific or architectural improvement. Compatibility protects users,
not technical debt: accidental internal behavior is not preserved
because it exists.

A breaking change follows this procedure:

- identify the affected public behavior and its users;
- keep the old name working for at least one minor release, emitting a
  `DeprecationWarning` that names the replacement;
- record the deprecation and the planned removal version in the
  changelog;
- remove the old name only after the deprecation period.

Until the first stable release, public API changes may be made
directly, without a deprecation period, provided they are recorded in
the changelog.

A behavior that is scientifically wrong (for example swapped
coordinates or a wrong unit) is a bug: it is fixed directly, with a
changelog entry, and is not preserved for compatibility.

A transitional module that rewrites an existing one (the `N` suffix,
as in `timeN` and `geodeticN`) is temporary. Its life cycle is:
develop and test the new module; migrate all internal users; replace
the old module by a deprecation shim that re-exports or wraps the new
implementation; rename the new module to the definitive name. Two
parallel architectures are never maintained indefinitely.

### 19. Dependencies and supported versions

Keep dependencies limited and justified. In order of preference: the
standard library, NumPy, SciPy, existing ShakeLab functionality, other
packages only when they provide substantial value. Do not add a
dependency to avoid a small amount of clear code. External libraries
may provide algorithms, formats or infrastructure, but ShakeLab owns
the domain abstractions that define its public API. In particular,
ObsPy is not an internal dependency.

Dependencies are of two kinds. Core dependencies are required by the
scientific core and listed in `dependencies`. Optional dependencies
(GUI, web interfaces, plotting backends, external adapters) are
declared as extras in `pyproject.toml`, used only in dedicated
modules, imported lazily, and raise an `ImportError` with an
installation hint when missing. Importing `shakelab` or any scientific
module never requires an optional dependency.

`pyproject.toml` is the single source of the package metadata, the
version and the dependency requirements; the supported versions are
listed in [conventions.md](conventions.md#6-supported-versions).

### 20. Applications, configuration and visualization

Applications consume reusable components. They may configure, select
implementations, orchestrate, expose services, manage jobs and
persistence, communicate with external systems and present results.
They do not duplicate algorithms provided by the core.

Command-line tools are named `slab-<name>` (`slab-seedlink`,
`slab-fdsn`, `slab-yarray`), live in the `cli` module of their
subpackage, also run as `python -m shakelab.<subpackage>`, and only
parse options and call the library. Where a command mirrors a Python
interface, its options have the names of the parameters (`--grid-max`
for `grid_max`).

Configuration determines how components are assembled and deployed,
never the scientific meaning of core objects. File-system locations,
host names, endpoints, credentials and operational defaults stay out
of scientific modules, which receive the scientific parameters they
need explicitly. Configuration precedence is explicit and predictable.

Visualization consumes results; it is never the canonical place where
they are computed. Results remain accessible without plotting,
rendering transformations do not alter scientific meaning, and GUI
event handling stays separate from scientific logic.

### 21. Module organization

Organize modules by coherent concepts, not by the chronology of
development. Before creating a module, determine the responsibility it
owns, its layer, what it depends on, what may depend on it and whether
an existing module already owns the concept. Avoid generic dumping
grounds such as an ever-growing `utils.py`: a helper that acquires
domain meaning moves to the concept that owns it.

Modules of a subpackage do not import the subpackage's own
`__init__`: definitions shared by its modules live in a dedicated
module (`signals/libio/common.py`).

Dated or numbered copies of a module (`providers_20260825.py`) and
example scripts do not belong in the package tree: history lives in
git, examples in `examples/` at the repository root.

A module that is experimental or no longer maintained is marked as
such in its module docstring and, when imported, emits a
`FutureWarning` (experimental) or `DeprecationWarning` (scheduled for
removal). Its retirement follows the deprecation procedure.

### 22. Testing and scientific validation

Scientific functionality must be testable without a GUI, a running
server, network access, a specific deployment or unrelated file-system
state. Separate deterministic calculations from external effects.

- Tests use pytest and live under `tests/`, mirroring the package
  structure. They follow the coding style of the library (file header,
  a docstring on every test). The test folders are not packages, so
  file names are unique across `tests/` (for example
  `test_yarray_data.py` in `tests/exploration/yarray/`).
- Every corrected defect gets a regression test that fails before the
  fix.
- Scientific tests compare against published reference values,
  analytical solutions or independent implementations, with the source
  cited in the test docstring.
- Tolerances are justified by the method (numerical precision,
  published rounding), never chosen merely to make a test pass. Use
  `numpy.testing` or `pytest.approx` with explicit `rtol`/`atol`.
- Tests do not need network access: services are replaced by local
  fake servers (SeedLink and FDSN in `tests/signals/`). A test that
  truly needs network access or large data is marked and skipped by
  default.
- Tests are also run with `-W error::UserWarning`, and those using
  sockets or files with `python -X dev -W error::ResourceWarning`, to
  find unexpected warnings and resources left open.
- The scripts in `examples/` are references for the user and are not
  run by the test suite (some need a server or user data); they are
  checked by hand when the API they use changes.
- The complete suite (`python -m pytest -q`) passes before any merge
  into the main branches. GitHub Actions runs it automatically on
  every push and pull request, on Linux and macOS, for all supported
  Python versions.

### 23. Documentation

Public scientific functionality documents, as appropriate: purpose,
scientific meaning, inputs and outputs, units, conventions,
assumptions, important parameters, exceptions and a short example
using the public API. Implemented methods cite their scientific
reference.

Each module starts with a docstring describing its purpose, and, where
relevant, its "Design choices" and "Conventions". Documentation does
not compensate for a confusing API: make the API clear first.

Repository documentation is organized as follows, without duplicating
information between places:

- `README.md` (repository root): overview, installation, license, and
  links to the documents below;
- `README.md` in every folder of the package: an introduction to that
  part of ShakeLab for its users (model below);
- `docs/`: user-facing documentation, organized by subpackage as the
  package (below), with the library-wide [conventions](conventions.md)
  and these guidelines at the top level;
- `examples/` (repository root): short scripts organized by subpackage
  as the package, with an index in `examples/README.md`. They
  illustrate interfaces, configurations and functionality, with
  generic file and server names (`array.json`, `data/*.mseed`,
  `rtserve.example.org`): they do not run on the test data and are
  not run by the test suite. Demonstrations on real data live outside
  the repository;
- `tests/`: organized by subpackage as the package;
- [`CHANGELOG.md`](../CHANGELOG.md) (repository root), recording
  user-visible changes, deprecations and scientific corrections.

The development history and the state of the work (blocks, commits,
open points) belong to the changelog and to the project notes, never
to READMEs or guides.

#### Organization of `docs/`

- [`docs/README.md`](README.md) is the single entry point: what
  ShakeLab is and the principles that guide it, told for its users, a
  short presentation of the subpackages, and the map of the whole
  documentation, linking every document. The root `README.md` links
  to it.
- Every documented subpackage has a folder with an overview named
  after it (`docs/signals/signals.md`,
  `docs/exploration/exploration.md`): a discursive introduction to the
  subpackage, how its parts work together, and links to its guides,
  to the READMEs of its folders and to its examples. The quick start
  stays in the README of the package; the overview links to it.
- Guides, one per topic, sit in the folder of their subpackage
  (`docs/signals/waveforms.md`, `docs/exploration/yarray.md`): theory
  where needed, complete use with a concrete example for every method,
  design choices, validation and limitations.
- Every document starts, above its title, with a navigation line back
  to the entry point and to the overview of its subpackage:
  `[ShakeLab documentation](../README.md) › [Signals](signals.md) ›
  Waveforms`.
- No other `README.md` is used inside `docs/`.

#### Figures

Figures support the explanations that need a picture (waveforms,
filters, spectra, spectrograms, array responses). They are didactic:
each shows the user what a concept, a method or an option does, or
when to choose one over another; comparisons that a sentence or a
table states as well (for example a numerical agreement with a
reference implementation) stay in the text. They are reproducible:

- the images are in `docs/figures/`, organized as the guides
  (`docs/figures/signals/waveforms_filters.png`), and are committed;
- the scripts that make them are in `docs/figures/scripts/`, one per
  guide (`signals_waveforms.py`, one function per figure), with a
  common style module (`style.py`) and `make_figures.py`, which
  regenerates all of them or those of one guide;
- the scripts run on the test data of the repository (an exception to
  the rule of generic names, which holds for examples), use
  `shakelab.graphics` where it offers the plot and matplotlib
  otherwise, and are deterministic (`seed`);
- PNG on an opaque white background, of moderate size (a few hundred
  kB at most, 1600 px wide for a full-width figure); SVG only for
  schematic drawings;
- schematic drawings (how the parts of a module work together, the
  flow of the data) are block diagrams with few blocks, each with a
  name and at most one short line, arrows that follow the flows, and
  shaded areas that group what belongs together (a thread, a
  module); the details stay in the text;
- each figure has a caption in the guide that states the data and
  their source (`*Figure 3. ...*`, below the image);
- figures are regenerated by hand when the code they show changes;
  the test suite does not run the scripts.

#### Equations and names in the documentation

Equations in the guides (and in READMEs, where needed) are written in
LaTeX, which GitHub renders with MathJax: displayed equations in a
fenced block with the language `math`, symbols in the text between
`` $` `` and `` `$ `` (for example `` $`L_{\mathrm{eff}}`$ ``), a form
protected from the Markdown syntax. Vectors and matrices are bold,
written with the Unicode mathematical bold characters (𝐚, 𝐂, 𝐤):
GitHub does not show `\mathbf` or `\boldsymbol` in bold. Some macros
are refused by GitHub (`\operatorname`: use `\mathrm{tr}`) or do not
exist in MathJax (`\bm`); matrices (`bmatrix`) work in `math` blocks.
A formula is checked in the preview of a GitHub issue before it is
committed when it uses a construct not already in the guides. Long
expressions go in a displayed block rather than in the text, so that
lines stay within 79 characters. Names of the code
(parameters, functions, classes) stay in backticks, next to the symbol
they correspond to (`` $`k_R`$ (`resolution_radius`) ``). Docstrings
stay in plain text.

Package names written as words are in capitals, for example YARRAY,
in the documentation (titles and text) and in the docstrings, comments
and messages of the code; in the code itself, in paths and in commands
they keep the lower-case name of the module (`yarray`,
`shakelab.exploration.yarray`, `slab-yarray`).

#### Folder READMEs

The README of a folder introduces its functionality to a user who
meets it for the first time: what it is for, how it is used, how it is
organized and where to read more. It is informative and didactic, in
plain descriptive prose, with a short paragraph of explanation before
and after each example; complete descriptions of options, theory and
formats belong to the guide in `docs/`. Every README has the same
sections, in the same order, numbered consecutively: the sections
marked optional are left out when they do not apply, and the
following ones take the next number.

```markdown
# Name: function in a few words

Introduction (two or three paragraphs): what the package does, for
whom and for which problems, where it fits in ShakeLab, and the ways
it is used (library, command line, real time...).

## 1. Features
One line per capability, five to ten items.

## 2. Quick start
The main objects, one complete minimal example (15-25 lines) and what
its result contains.

## 3. Common tasks
### 3.1 ... (a task as title: "Build the array data", "Plan a survey")
Three to six tasks, each with an explanation, a short example (5-15
lines) and the essential options.

## 4. Command line            (optional: packages with a command)
## 5. Format                  (optional: packages defining a format)
## 6. Structure
Table "Module | Main objects | Role", and a few lines on layers,
data flow and generic tools used from other packages.

## 7. Conventions             (optional: only those of the package)
## 8. See also
The guide in docs/, the examples, related packages.

## 9. References              (optional: scientific packages)
```

Rules of style:

- Present tense, no dates, versions, commit references, development
  states or comparisons with earlier implementations; numbers that
  change with the code (processing times) belong to the guide.
- The examples follow the rules of `examples/` (generic names, one
  object per step, `info()`, short comments on the options) and use
  only the public API; every example is checked against the code
  before it is committed.
- Lines of at most 79 characters, except table rows; code blocks with
  their language (`python`, `sh`, `json`); relative links to the
  guides and examples.
- Library-wide conventions are not repeated: link to
  [conventions](conventions.md).
- When a folder contains subfolders with their own README, its README
  presents them in "Structure" and links to them.

The model README is `shakelab/exploration/yarray/README.md`.

### 24. Coding conventions and naming

ShakeLab is written in Python and follows PEP 8. English is used for
identifiers, comments, docstrings, documentation, commit messages and
all repository-facing text. The style below was set on
`engineering/impact.py`; for structure, names and sections, a new or
revised module is compared with the most recently revised module of
the same kind (`signals/libio/sac.py` for formats, `signals/seedlink`
for network clients, `signals/spectrum.py` for data derived from a
record):

- the standard ShakeLab license header at the top of every file,
  followed by the module docstring (see below);
- lines of at most 79 characters;
- type hints on public functions and methods, with
  `from __future__ import annotations`; `Optional[X]` for values that
  may be None (not `X | None`), built-in generics (`list[str]`,
  `dict[str, Any]`, `tuple[int, ...]`) rather than `typing.List` and
  similar;
- NumPy-style (numpydoc) docstrings: a summary line ending with a
  period; multi-line docstrings opened on the line after the quotes;
  sections in numpydoc order (Parameters, Attributes, Returns or
  Yields, Raises, Warns, See Also, Notes, References, Examples) with
  underlines as long as their titles; every parameter entry written
  `name : type`, with `, optional` when it has a default (alternatives
  as `str, path-like or file object`; `**kwargs` without a type); a
  Returns (or Yields) section in functions that return a value;
  parameters grouped on one line (`x, y : float`) only when they share
  type and description; module docstrings may add free sections
  (Conventions, Design choices, Reading, Writing);
- comment separators: module sections as three lines of 78 characters
  (`# ` followed by 76 dashes, the title, the dashes again); sections
  inside a class as one line of 73 characters (`    # -- title ---`);
- imports grouped as standard library, third-party packages and
  ShakeLab, each group in alphabetical order; explicit imports, never
  wildcard imports.

The license header is:

```python
# ****************************************************************************
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
# ****************************************************************************
```

Files of third parties included in the repository keep their own
license header.

Prefer descriptive names, explicit control flow, small coherent
functions, readable numerical expressions, limited side effects and
deterministic behavior. Avoid unnecessary metaprogramming, obscure
constructs, deep inheritance, deeply nested logic, global mutable
state and broad exception handling. Comments explain why, not what.

Names describe concepts, not implementation accidents. Use
well-defined domain terminology and conventional abbreviations only
(PGA, SA, Vs30, Rjb). Use one name per concept across the library, and
never the same name for different quantities. Where a quantity is not
in its canonical unit, the unit appears in the name (`depth_km`,
`distance_km`). Naming is part of the scientific API.

### 25. Designing new components

Before writing a substantial component, answer:

- **Scientific concept**: what does it represent, measure, transform
  or calculate?
- **Inputs and outputs**: what are its scientific boundaries, units
  and conventions?
- **Layer**: where does the responsibility belong?
- **Reuse**: which existing representations or operations apply?
- **Canonical implementation**: where does the authoritative behavior
  live?
- **Public API**: what is the smallest useful interface?
- **Extensibility**: which variations are real and which hypothetical?
- **Dependencies**: what external functionality is actually needed?
- **Performance**: what data sizes and access patterns are expected?
- **Validation**: how will correctness be established and tested?

Derive the structure from these answers, not from an arbitrary class
hierarchy.

### 26. Refactoring and existing code

Existing code is the current implementation, not the architecture. It
may contain historical decisions, prototypes, transitional APIs and
technical debt. Infer the architecture from the scientific concepts
and these principles.

When existing code conflicts with them: identify the inconsistency,
check whether it is intentional, assess compatibility constraints
(Section 18), decide the direction and improve incrementally. A
refactor must improve something concrete: responsibility, duplication,
dependency direction, API coherence, transparency, testability or
performance. Do not preserve poor architecture because it exists, and
do not rewrite mature, coherent code because another style is also
valid.

### 27. Decision rule and core principle

When several designs are reasonable, prefer the one that preserves
scientific meaning most explicitly, gives each responsibility one
location and each operation one implementation, keeps dependencies
directed, exposes the smallest coherent API, avoids assumptions of the
immediate use case and remains understandable without extensive
framework knowledge. Optimize for conceptual clarity, scientific
reliability and long-term maintainability, not for the smallest number
of lines.

```
  a small set of coherent scientific concepts
+ one authoritative implementation of each operation
+ explicit scientific semantics
+ composable domain abstractions
+ directed dependencies
= a framework that can grow without losing coherence
```

A user should not need to understand an application to use a
scientific algorithm, and an application should not need to reproduce
an algorithm to use it.

