# Baysor native implementation roadmap

Date: 2026-09-01

## Purpose

This document isolates the work from the cross-repository
[scalable Baysor integration roadmap](roadmap.md) that is owned by the
[`vibspatial/Baysor`](https://github.com/vibspatial/Baysor) C++ repository.
It is the implementation sequence for the reusable native engine and the Baysor
CLI. The combined roadmap remains the source for product-wide architecture,
Python packaging, SpatialData integration, tiled execution, reconciliation, and
scientific validation.

The repository boundary is:

> `vibspatial/Baysor` owns reusable scientific C++ operations and the CLI;
> `baysor-python` owns bindings, Python APIs, packaging, SpatialData, and Dask
> orchestration.

Native work is completed, reviewed, and tested here. `baysor-python` consumes an
exact reviewed commit and must not carry private patches against its Baysor
checkout.

## Scope

This repository owns:

- the source-level C++ segmentation request, result, cancellation, and execution
  contracts;
- extraction of the scientific run from the CLI into `baysor_lib`;
- reusable native serializers and reporting facilities over the structured
  result;
- an embeddable CMake library target with no Python dependency;
- the CLI frontend and its parity with the immutable upstream reference;
- a reusable native boundary-estimation operation;
- stable transcript-identity preservation in native results and Parquet output;
- optional native output selection for tiled workloads;
- globally reusable confidence/noise calibration needed for scientifically
  supported tiled execution; and
- focused native tests for correctness, ownership, cancellation, repeated calls,
  serialization, and CLI parity.

This repository does **not** own:

- nanobind code or Python types;
- `scikit-build-core`, Python wheels, or Python Stable-ABI decisions;
- `baysor_python.segment(...)`, `boundaries(...)`, or SpatialData adapters;
- tile planning, Dask workers, retries, manifests, or worker recycling;
- cross-tile cell reconciliation or Python-side global product construction; or
- Harpy integration.

These exclusions keep the native library independently buildable and testable.
No native public header may require Python, nanobind, SpatialData, Dask, Arrow
objects supplied by PyArrow, or other Python-layer types.

## Target architecture

The authoritative workflow will be one coarse-grained C++ operation:

```text
                   C++ run_segmentation(...)
                              |
                +-------------+-------------+
                |                           |
                v                           v
         Baysor CLI                  nanobind module
    parse CLI arguments              convert Python input
    call shared function             release Python GIL
    serialize outputs                call shared function
    map errors to exit codes         translate exceptions
```

Only the left-hand frontend is implemented in this repository. The right-hand
frontend is shown to define the consumer boundary; it is implemented in
`baysor-python`.

The public native API is a normal C++17 source-level API. It is not a frozen
binary SDK and must not expose internal BMM types merely to make them bindable.
The request and result should be coarse enough that both frontends can use them
without duplicating scientific orchestration or parameter resolution.

## Implementation sequence

The native work is divided into small, reviewable slices. The mapping to the
combined roadmap is explicit so that the two documents cannot silently diverge.

| Native slice | Combined-roadmap responsibility | Required handoff |
| --- | --- | --- |
| N0: Native baseline and regression fixture | Foundation of Slice 1A.1 | None |
| N1: Public segmentation contracts | Slice 1A.1 | None |
| N2: Scientific orchestration extraction | Slice 1A.1 | None |
| N3: Native lifecycle and embeddability gate | Slice 1A.1 | `baysor-python` pins the reviewed commit in Slice 1A.2 |
| N4: CLI frontend and parity | Slice 1B | `baysor-python` advances to the reviewed parity commit |
| N5: Public boundary operation | Native portion of Phase 2 | `baysor-python` binds the reviewed operation |
| N6: Stable transcript identity | Tiled correctness prerequisite | `baysor-python` removes any row-order compatibility path |
| N7: Selective native products | Tiled performance optimization | Tile calls request only authoritative products |
| N8: Reusable confidence/noise calibration | Tiled production prerequisite | Tile calls apply one global calibration |
| N9: Reusable graph calibration, if required | Conditional tiled-quality mitigation | Implement only when validation supplies evidence |

The dependency order is:

```text
N0 -> N1 -> N2 -> N3 -> N4
                         |
                         +-> N5 -> native boundary handoff
                         +-> N6 -> stable-identity handoff
                         +-> N7 -> selective-output handoff
                         +-> N8 -> global-calibration handoff
                                      |
                                      +-> N9 only if validation requires it
```

N5 through N8 may be developed as separate branches after N4. They must each
retain the N4 CLI-parity gate and be handed to `baysor-python` as explicit commit
updates.

## Immediate next work

The repository development environment, CMake test preset, and individual
GoogleTest discovery are in place. The next product work is Native Slice N0:

1. verify the complete native test suite from a clean configure and build;
2. select the small segmentation regression input and locked CLI configuration;
3. record the immutable reference CLI's semantic outputs; and
4. review that fixture as the behavioural gate for N1 through N4.

Extraction of `run_segmentation(...)` begins only after this reference is
reproducible. The deferred actual-UCB experiment is not a prerequisite for this
work.

### Native Slice N0: Establish the native baseline and regression fixture

This slice makes native behavior measurable before orchestration is moved.

Deliverables:

- preserve the immutable upstream `cpp-0.8.3` baseline at commit
  `d7077a7ded6f4b941915badc894f767532d39fd2`;
- perform maintained work on a dedicated feature branch rather than changing the
  baseline branch;
- retain the reproducible CMake test preset and documented user-space macOS
  development environment;
- keep individual GoogleTest cases discoverable through CTest and development
  tools;
- select a small, deterministic-enough native regression fixture that exercises
  molecule loading, an optional prior, parameter resolution, confidence fitting,
  segmentation, and the output serializers; and
- capture the reference CLI invocation, resolved parameters, native dependency
  identity, OpenMP settings, and semantic outputs for that fixture.

The regression fixture must be small enough for routine native CI. It is not the
deferred actual-UCB reference experiment. The reference comparison must account
for harmless cell relabelling and polygon representation differences while still
checking the scientific result.

Exit criterion: a clean checkout can configure, build, and run its native tests,
and the untouched reference CLI can produce the recorded fixture outputs under
locked settings.

### Native Slice N1: Define the public segmentation contracts

Define the source-level API before moving implementation out of the CLI.

Deliverables:

- `SegmentationRequest`, containing the molecule-input specification, optional
  prior specification, filtering, segmentation, clustering, confidence/noise,
  requested-product, and execution options;
- `SegmentationResult`, owning retained molecule identity and data, cell/noise
  assignments, confidence fields, optional clusters, cell statistics,
  boundaries, counts, convergence diagnostics, resolved options, and native
  provenance;
- a thread-safe C++17 cancellation token with a documented ownership and
  cross-thread use contract;
- a distinct cancelled outcome that cannot be mistaken for a complete result;
- documented invalid-request, input/output, native-processing, and cancellation
  error behaviour; and
- public headers that use ordinary C++ types and do not expose `argv`, CLI
  aliases, process exit policy, presentation strings, Python objects, or internal
  BMM ownership.

Data-dependent defaults such as inferred scale or initial component count belong
to the shared operation and must be returned as resolved values. A frontend may
parse user input into the request, but it may not independently reproduce this
resolution.

Exit criterion: the request, result, cancellation, error, and ownership contracts
are reviewable as a coherent native API and can be compiled without the CLI.

### Native Slice N2: Extract the scientific orchestration

Implement `run_segmentation(request, cancellation)` in `baysor_lib` by moving the
existing orchestration out of `cmd_run(...)`.

The shared operation owns this sequence:

```text
validate typed request
        |
load and filter molecules; load optional prior
        |
resolve data-dependent options and the confidence/noise model
        |
build the molecule graph and perform optional molecule clustering
        |
initialize and run the existing 2D or 3D BMM
        |
derive requested assignments, statistics, boundaries, counts, and diagnostics
        |
return an owned SegmentationResult
```

Deliverables:

- move scientific orchestration and canonical option resolution into the public
  operation without changing the underlying algorithms;
- reuse the existing loaders, prior handling, confidence estimator, graph and
  clustering implementations, BMM, boundary estimator, and scientific
  calculations;
- separate serialization into reusable native functions over
  `SegmentationResult`;
- use moves and normal C++ lifetime management so returned products do not refer
  to destroyed working objects;
- keep output-directory policy, filename selection, command-line presentation,
  and process exit codes outside the operation; and
- replace process-wide logging manipulation with an explicitly supplied progress
  or logging interface where library reporting is required.

This is an orchestration and ownership refactor. It is not an algorithm rewrite,
a Python port, or an opportunity to change scientific defaults.

Exit criterion: a focused C++ test invokes a complete segmentation and inspects
its owned result without calling CLI or Python code.

### Native Slice N3: Complete cancellation, lifecycle, and embeddability

Harden the extracted operation for repeated in-process use and hand it off as a
clean CMake consumer target.

Deliverables:

- add cooperative cancellation checks after major phases and at safe BMM
  iteration boundaries;
- ensure cancellation publishes no complete result and does not run completed
  result serializers;
- test invalid requests and structured native failures;
- test result ownership after temporary working state has been destroyed;
- test repeated same-process calls for hidden global-state leakage;
- verify that no call replaces a process-wide logger, calls `std::exit`, or
  relies on CLI-owned mutable state;
- make `baysor_lib` consumable through `add_subdirectory(...)`, including
  position-independent code and transitive build requirements needed by an
  embedding consumer; and
- retain standalone CLI and native-test builds without introducing nanobind or a
  Python dependency.

The repeated-call test should distinguish an allocator's bounded retained-memory
plateau from linear growth in live allocations. Hard process reclamation remains
a supervisor policy outside this repository.

Exit criterion: direct native tests cover success, invalid input, cancellation,
owned-result lifetime, and repeated calls; a clean CMake consumer can link the
library without modifying the Baysor checkout.

Handoff: select and review the green N3 commit. `baysor-python` pins that exact
commit as its Slice 1A.2 source dependency and records both the upstream baseline
and fork commit.

### Native Slice N4: Make the CLI a frontend over the shared operation

Refactor the CLI only after the direct native operation is established.

Deliverables:

- reduce `cmd_run(...)` to CLI/config parsing, construction of
  `SegmentationRequest`, invocation of `run_segmentation(...)`, shared
  serialization, user-facing reporting, and conversion of outcomes to exit
  codes;
- remove every second implementation of data-dependent parameter resolution or
  scientific orchestration from the CLI;
- preserve existing supported command-line arguments, configuration behaviour,
  output naming, and default full-output behaviour;
- compare the refactored CLI with the immutable reference CLI on the N0 fixture
  using the same inputs, options, OpenMP settings, and output mode; and
- retain reference and candidate logs and resolved options as parity evidence.

Parity includes retained transcript identity, cell/noise assignments, molecule
confidence, count matrices, cell statistics, resolved parameters, and Baysor
revision. Cell identifiers and polygons are compared semantically after
normalizing harmless relabelling, polygon orientation, and starting-vertex
differences. Because the reference CLI has no general segmentation seed, parity
runs should be repeated when measuring parallel variability.

Exit criterion: the refactored CLI is semantically equivalent to the reference
CLI and contains no scientific orchestration path independent of
`run_segmentation(...)`.

Handoff: `baysor-python` advances its pinned source to the reviewed N4 commit
before implementing or releasing the native Python binding.

### Native Slice N5: Establish a reusable boundary operation

Expose the existing boundary semantics as a stable, independently callable C++
library operation. The Python array conversion belongs in `baysor-python`, not
here.

Deliverables:

- define a coarse boundary request and owned result around the existing Baysor
  estimator;
- accept the complete molecule coordinates and final cell/noise assignments;
- support sparse global cell identifiers;
- support optional target-cell selection for bounded-memory batches while still
  accepting all contextual molecules required by those target cells;
- accept or compute one global boundary-distance parameter and allow the same
  value to be reused across batches;
- define explicit behaviour for empty, one-molecule, two-molecule, collinear,
  and ordinary cells;
- preserve existing supported 2D and 3D semantics without introducing a second
  boundary algorithm; and
- add focused determinism, thread-count, edge-case, target-selection, and
  CLI-geometry parity tests.

The operation must not merge tile-local polygons. Its purpose is to estimate one
new boundary from the final reconciled molecule assignments supplied by a
caller.

Exit criterion: identical coordinates, assignments, parameters, and context
produce geometry semantically equivalent to the CLI boundary output, both for a
complete call and for supported target-cell batches.

Handoff: `baysor-python` binds this reviewed operation and owns packed-array
conversion, GIL release, exception translation, Shapely conversion, and
SpatialData shapes.

### Native Slice N6: Round-trip stable transcript identity

Make transcript identity an explicit native correctness contract rather than an
assumption about output row order.

Deliverables:

- accept a configured signed 64-bit transcript identifier from supported table
  inputs;
- preserve it through filtering, compaction, segmentation, and result ownership;
- emit it in the authoritative Parquet molecule output without lossy conversion;
- validate length, type, and uniqueness where the native contract requires it;
- retain existing input compatibility when no identifier is supplied; and
- test Parquet round trips, filtered rows, large and sparse identifiers, and
  invalid or duplicate identifiers.

The current source-specific transcript field is useful groundwork but is not the
complete generic Parquet round-trip contract. A Python compatibility adapter may
temporarily validate row identity, but it must not remain the production tiled
identity mechanism after this slice.

Exit criterion: every retained input transcript carrying an identifier has the
same identifier in the structured native result and Parquet molecule output,
independently of filtering or cell assignment.

### Native Slice N7: Add selective scientific products

Avoid computing and writing tile-local products that the tiled workflow will
discard.

Deliverables:

- add typed request flags for authoritative molecule assignments and optional
  boundaries, cell statistics, count matrices, diagnostic colour embeddings,
  HTML reports, and other expensive products;
- ensure disabled products are neither calculated nor serialized;
- preserve the CLI's current full-output defaults;
- record the requested and produced product set in the resolved result and
  provenance; and
- prove that a molecule-focused run returns the same assignments and confidence
  values as a full-product run under the same scientific settings.

This slice is a runtime, memory, and storage optimization. It does not relax the
stable-identity or scientific-parity gates, and it is not required to establish
the untiled native API.

Exit criterion: a caller can request the molecule products needed for tiled
reconciliation without paying for unused tile boundaries, matrices, or reports,
while the default CLI result remains unchanged.

### Native Slice N8: Expose reusable confidence/noise calibration

Remove the most important invocation-wide statistical inconsistency between
independent tile runs.

Deliverables:

- select and document at least one supported native contract:
  - accept and preserve globally precomputed per-molecule confidence; or
  - fit, export, and apply a reusable signal/noise calibration;
- include the KNN definition, fitted mixture parameters, convergence information,
  and summary signal/noise proportions in the calibration result and provenance;
- let `SegmentationRequest` distinguish the existing fit-for-this-run behaviour
  from the supported global-calibration path;
- retain the existing untiled and CLI default behaviour unless the caller
  explicitly selects the reusable path;
- validate calibration inputs and fail closed on incompatible KNN definitions or
  parameters; and
- test full-data fit/export/apply, application to partitions, repeatability, and
  parity of the unchanged default path.

Engineering tiled runs may use tile-local fitting only as an explicitly recorded
non-production validation mode. Production-supported tiled execution requires
this slice and the corresponding `baysor-python` application path.

Exit criterion: separate segmentation calls can apply one compatible global
confidence/noise model and report enough provenance to prove that they did so.

### Native Slice N9: Add reusable molecule-graph calibration only if required

This is a conditional mitigation, not scheduled native work. Per-invocation
molecule-graph edge quantiles must first be measured during tiled-versus-untiled
validation.

Trigger: begin this slice only if tile-wide graph-threshold variation produces
material assignment differences after global gene filtering, fixed scale,
disabled or globally consistent clustering, sufficient halos, and N8 global
confidence handling are in place.

If triggered, deliverables are:

- an explicit fit/export/apply contract for graph thresholds or equivalent
  globally calibrated graph parameters;
- parameter validation and provenance;
- unchanged default untiled behaviour; and
- partitioned-versus-global application tests plus CLI-parity coverage.

Exit criterion: all tile calls can apply the same validated graph calibration
without changing the reference untiled default path.

## Cross-repository handoff contract

Every native handoff follows the same sequence:

1. implement and run focused tests in `vibspatial/Baysor`;
2. pass the applicable direct-native and CLI-parity gates;
3. review and commit the native change here;
4. select the exact green integration commit;
5. advance `baysor-python/vendor/Baysor` to that commit; and
6. record the upstream baseline and integration commit in Python build
   provenance.

The consumer repository must not patch a dirty submodule. If Python integration
reveals a native defect or missing generic CMake capability, the correction is
made and tested here, followed by a new commit-based handoff.

## Native test strategy

Each slice runs only the focused native tests needed for its change during
development. The complete native suite and reference CLI parity fixture are
release and handoff gates.

Required coverage across the sequence includes:

- request validation and data-dependent option resolution;
- complete direct segmentation through the shared operation;
- 2D and supported 3D dispatch;
- prior and no-prior paths;
- result ownership after working state destruction;
- cancellation before execution, between phases, and at a safe BMM iteration
  boundary;
- repeated same-process execution and global-state isolation;
- serializer consistency between direct and CLI paths;
- semantic parity with the immutable reference CLI;
- boundary edge cases, sparse cell identifiers, target-cell batches, and
  thread-count stability;
- stable transcript-ID round trips through filtering and Parquet;
- selective-output equivalence for authoritative assignments; and
- confidence/noise calibration fit/export/apply consistency.

Tests must use explicit OpenMP settings and record them in parity evidence.
Bitwise identity is not required for values with legitimate parallel or geometric
representation variability, but every normalization and tolerance must be
scientifically justified and documented.

## Native production gates

A native slice is complete only when:

- its public contract and supported error behaviour are documented;
- focused tests pass on the supported native platforms relevant to the slice;
- the change introduces no Python or application-layer dependency;
- default CLI scientific behaviour remains unchanged unless the roadmap
  explicitly authorizes and validates a change;
- the implementation has one authoritative scientific path;
- ownership, cancellation, logging, and process-global behaviour are safe for an
  embedding process; and
- the reviewed commit is suitable for an auditable, exact consumer pin.

The actual-data reference, Dask lifecycle, reconciliation quality, seam metrics,
grid-shift experiments, SpatialData round trips, and wheel release matrix remain
product gates in the combined roadmap and in `baysor-python`; they are not
duplicated as native implementation work here.

## Explicitly deferred native redesign

This roadmap does not implement tiling as a C++ loop over independent tiles.
Python will orchestrate coarse native tile calls, and the C++ repository will
provide only the scientific operations and globally reusable parameters those
calls require.

A genuinely native tiled Baysor algorithm would require domain decomposition of
the molecule graph and exchange of assignments, component statistics, and cell
lifecycle decisions during BMM iterations. That is a separate scientific
redesign. It should be considered only if the validated Python-orchestrated
core-plus-halo approach cannot meet the production quality gates after N8 and,
when evidence requires it, N9.
