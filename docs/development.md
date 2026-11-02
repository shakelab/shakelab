# ShakeLab development and architecture guidelines

*Version 2.1, 31 October 2026*

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
simulations, Monte Carlo, random meshes) accept an explicit
`numpy.random.Generator` or a seed, through a parameter conventionally
named `rng`. They never use the global NumPy or Python random state.
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
exception when callers need to distinguish the case.

Never catch exceptions broadly to let execution continue, and never
replace an invalid operation with a scientifically different fallback.
A fallback is acceptable only when it preserves the intended meaning
and is explicit.

Raise an exception when the result would be wrong. Issue a warning
(`warnings.warn`, with an appropriate category) when the result is
valid but the user should know something: an approximation, a
deprecated call, a value outside the calibrated range of a model. Use
the `logging` module, with per-module loggers, for operational
information. Library code does not use `print` for debugging or error
reporting; the only exception are methods whose purpose is to show a
summary to the user, such as the `info()` methods of the waveform
classes.

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
correctness takes precedence over superficial uniformity.

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
  a docstring on every test).
- Every corrected defect gets a regression test that fails before the
  fix.
- Scientific tests compare against published reference values,
  analytical solutions or independent implementations, with the source
  cited in the test docstring.
- Tolerances are justified by the method (numerical precision,
  published rounding), never chosen merely to make a test pass. Use
  `numpy.testing` or `pytest.approx` with explicit `rtol`/`atol`.
- Tests needing network access or large data are marked and skipped by
  default.
- Every script in `examples/` is run by `tests/test_examples.py` with a
  non-interactive plotting backend, so that the examples stay in line
  with the library.
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
- `docs/`: user-facing guides by topic, the library-wide
  [conventions](conventions.md) and these guidelines;
- `examples/` (repository root): runnable scripts, organized by
  subpackage as the package and the tests (`examples/signals/`, ...),
  using the data in `tests/<subpackage>/data/`, with an index in
  `examples/README.md`;
- `README.md` in each subpackage: its modules, design choices,
  supported formats and specific conventions, with links to the
  guides in `docs/` and the scripts in `examples/`;
- [`CHANGELOG.md`](../CHANGELOG.md) (repository root), recording
  user-visible changes, deprecations and scientific corrections.

### 24. Coding conventions and naming

ShakeLab is written in Python and follows PEP 8. English is used for
identifiers, comments, docstrings, documentation, commit messages and
all repository-facing text. The reference style is that of
`engineering/impact.py`:

- the standard ShakeLab license header at the top of every file,
  followed by the module docstring (see below);
- lines of at most 79 characters;
- type hints on public functions and methods, with
  `from __future__ import annotations`;
- NumPy-style (numpydoc) docstrings; when type hints are present, the
  Parameters section may list names and descriptions only;
- comment separators of the form `# ----` (or `# ====` for major
  sections) to structure long modules;
- explicit imports, never wildcard imports.

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

