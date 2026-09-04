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
- the CLI frontend and its parity with the recorded pre-extraction reference;
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

### Data transport and object ownership

The selected segmentation architecture is deliberately hybrid. Paths are used
at the large-data transport boundaries, while the scientific implementation
operates on owning C++ objects:

```text
Python or SpatialData objects
            |
            | prepare stable tabular input
            v
CSV or Parquet path + typed SegmentationRequest
            |
            | nanobind calls the C++ library directly; no CLI subprocess
            v
owning C++ MoleculeData and algorithm working state
            |
            v
owned C++ SegmentationResult
            |
            | reusable native serializers
            v
output artifacts on disk + paths and provenance
            |
            v
Python result objects or coordinated SpatialData elements
```

The initial `SegmentationRequest` is therefore path-oriented: its molecule input
is a typed file-source specification for supported CSV or Parquet data, not a
Python object or an unstructured CLI argument. `run_segmentation(...)` loads that
source and then operates entirely on native objects. A path-oriented nanobind
call is still an in-process native-library call with the GIL released; it does
not invoke the Baysor executable, parse CLI arguments, or start a subprocess.

`SegmentationResult` means the owned, in-memory C++ scientific result. It is
distinct from a Python-layer artifact descriptor, here called
`SegmentationArtifacts`, which contains validated output paths, status, and
provenance after native serialization. The CLI may serialize the native result
and map errors to exit codes. A Python worker may serialize the same native
result and return `SegmentationArtifacts` without exposing or transferring the
complete native object graph to Python. A higher-level SpatialData adapter may
accept and return Python objects while keeping this staging and serialization
boundary as an implementation detail.

This choice supports large datasets, bounded worker memory, retryable and
inspectable outputs, deterministic worker recycling, and future Dask execution.
In particular, tile tasks exchange only immutable descriptors and artifact
paths; complete transcript clouds and native results do not pass through the
Dask scheduler. The costs are staging I/O, temporary-storage requirements, and
some additional latency for small interactive runs.

A direct in-memory segmentation input from NumPy, PyArrow, or another Python
object is not part of the initial contract. Such an API is not automatically
zero-copy: it requires explicit dtype, layout, string/category, mutation,
lifetime, and GIL-release rules. It may be added as a second input-source form
only when profiling demonstrates a material benefit and its ownership and
copying semantics have been designed and tested. The path-oriented form remains
the required scalable and Dask-compatible contract.

This is also a scope decision for the native implementation, not a limitation of
nanobind. Nanobind could accept NumPy buffers, but the current canonical
`MoleculeData` owns its columns in `std::vector` instances and the processing
pipeline constructs additional Eigen and BMM working storage. A straightforward
array binding would therefore still copy the eager Python data into Baysor-owned
memory. Making that route genuinely zero-copy would require a broader refactor of
loading, filtering, compaction, gene encoding, matrix construction, and borrowed
buffer lifetimes. N1 through N4 must not combine the segmentation extraction with
that data-model refactor.

The initial file route instead reuses the existing native Arrow/Parquet loader as
the one canonical implementation of projection, filtering, gene encoding, prior
handling, and `MoleculeData` construction. This reduces scientific, memory, and
ownership risk while the public operation is extracted. The native roadmap does
not introduce a second array-input route speculatively.

The Python integration must measure prepared-input write time, native input-read
time, segmentation time, peak RSS, and temporary-storage volume on representative
untiled and tile-sized inputs. An in-memory input may be proposed as a separate
reviewed slice only if those measurements show that staging is material. Its
first version may still use an explicit copy into `MoleculeData`; borrowed or
Arrow-backed zero-copy processing requires its own correctness, lifetime, and
memory benchmarks.

The reusable boundary operation is intentionally different. It accepts packed
native coordinate and assignment arrays because it is a focused operation with
a small, explicit buffer contract and is needed after Python has reconciled
tile-local assignments. This exception does not require the complete
segmentation API to expose Python-owned arrays or internal BMM objects.

### Randomness and reproducibility contract

The normative caller-facing guarantee lives in the root-level
[native segmentation API contract](../contract.md). This roadmap records its
design rationale and the implementation and verification work needed to sustain
it.

The current feature branch, before extraction begins, is the behavioural
starting point. N0 records the exact committed source revision used to build the
reference CLI; the historical upstream commit remains source provenance rather
than a second behavioural oracle. The reference CLI uses fixed internal seeds
but exposes no general seed option. Its reference is generated in a fresh
process with one OpenMP thread, records the implementation's implicit
segmentation seed of `1`, and is repeated to verify semantic stability.

The maintained native API will make randomness explicit after that baseline is
captured:

- `SegmentationRequest` carries a 64-bit `random_seed`, defaulting to `1` for
  compatibility with the reference segmentation stream;
- `run_segmentation(...)` creates run-local random state from that value rather
  than relying on resetting mutable process-global state;
- segmentation-affecting code receives explicit state from a versioned run-level
  stream and named substream contract, with the default seed mapped to the N0
  compatibility behaviour;
- diagnostic-only operations use separate derived streams so enabling them cannot
  consume the scientific stream and change assignments;
- the resolved seed and effective OpenMP settings are returned in provenance;
  and
- the refactored CLI exposes the same setting as `--seed`, while the Python
  frontend forwards that native setting rather than defining another seed.

For the same input, options, build, platform, seed, and single-thread execution,
repeated calls must produce the same semantic result, including repeated calls
inside one process. A seed alone does not guarantee bitwise-identical
multi-threaded results: the current stochastic BMM uses per-thread random streams
with dynamic OpenMP scheduling, and parallel floating-point work may also vary in
its final bits. Multi-threaded comparisons must therefore lock and record the
thread configuration and characterize residual variability. Deterministic
random draws independent of worker scheduling, for example streams derived from
the run seed, iteration, and molecule identity, require a separately reviewed
change and are not an N0 gate.

## Implementation sequence

The native work is divided into small, reviewable slices. The mapping to the
combined roadmap is explicit so that the two documents cannot silently diverge.

| Native slice | Combined-roadmap responsibility | Required handoff |
| --- | --- | --- |
| N0: Native baseline and regression fixture | Foundation of Slice 1A.1 | None |
| N1: Public segmentation contracts | Slice 1A.1 | None |
| N2: Scientific orchestration extraction | Slice 1A.1 | None |
| N2b: Scientific parity matrix | N2 review and merge hardening | None |
| N3: Native lifecycle and embeddability gate | Slice 1A.1 | `baysor-python` pins the reviewed commit in Slice 1A.2 |
| N4: CLI frontend and parity | Slice 1B | `baysor-python` advances to the reviewed parity commit |
| N5: Public boundary operation | Native portion of Phase 2 | `baysor-python` binds the reviewed operation |
| N6: Stable transcript identity | Tiled correctness prerequisite | `baysor-python` removes any row-order compatibility path |
| N7: Selective native products | Tiled performance optimization | Tile calls request only authoritative products |
| N8: Reusable confidence/noise calibration | Tiled production prerequisite | Tile calls apply one global calibration |
| N9: Reusable graph calibration, if required | Conditional tiled-quality mitigation | Implement only when validation supplies evidence |

The dependency order is:

```text
N0 -> N1 -> N2 -> N2b -> N3 -> N4
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

Native Slice N0 is implemented and verified. It records commit
`f46a1e1dce1606d0ea644f4f8f1cf682597ba65c` as the pre-extraction behavioural
baseline, uses the small legacy-output fixture under
`tests/fixtures/native_baseline`, and adds one CTest-discovered semantic CLI
regression. The N0 verification build passed all 115 tests then present,
including that regression.

Native Slice N1 is complete. Commit
`f74417481feb988952eedc22c688e1aefb5e8eb7` adds the public owning request and
result types, produced-product tracking, the versioned random-stream contract,
structured native errors, the cross-thread cancellation source/token pair, and
the durable native API contract. Its focused contract tests and the N0 semantic
CLI regression pass.

Native Slice N2 is complete in the current working tree as of 2026-09-03. The
public `run_segmentation(...)` operation now performs canonical option
resolution, molecule and prior loading, confidence estimation, graph and
clustering setup, 2D/3D BMM execution, and owned scientific-result
materialization. Random state is run-local and derived from the N1 substream
contract, result-based serializers cover the locked legacy products, and the
direct operation plus unchanged CLI both pass the N0 semantic fixture.

Native Slice N2b is complete in the current working tree as of 2026-09-04. Its
versioned pre-N2 references and generation recipe live under
`tests/fixtures/scientific_parity`. The current CLI and direct operation both
match the duplicate-coordinate, MRF, Louvain, Leiden, and 3D references; focused
tests also cover legacy-seed equivalence and the NCV-colour product boundary.

Native Slice N3 is complete in the current working tree as of 2026-09-04. The
operation now checks cancellation between major phases and at complete BMM
iteration boundaries, restores caller OpenMP state on every exit path, rejects
nested OpenMP and overlapping same-process calls, and records unambiguous
OpenMP and Arrow execution provenance. Focused lifecycle tests cover structured
failures, cancellation, repeated calls, thread restoration, and unsupported
entry contexts. `baysor::baysor` is a position-independent CMake consumer target,
and the independent `tests/cmake_consumer` project configures, links, and runs
without building the CLI.

The next implementation slice is N4, which makes the CLI a frontend over the
shared operation and repeats the semantic parity gates.

The deferred actual-UCB experiment remains outside this sequence.

### Native Slice N0: Establish the native baseline and regression fixture

Implementation status: complete as of 2026-09-01.

This slice creates the trustworthy "before" result for the later extraction of
`run_segmentation(...)`. The existing native tests exercise many individual
components, but they do not by themselves prove that moving the complete
workflow out of `cmd_run(...)` preserves the assembled scientific result. N0
therefore adds one small end-to-end CLI fixture and an automated semantic
comparison before orchestration is moved.

N0 does not define the new public segmentation API, move code out of the CLI,
change an algorithm or default, add Python integration, or run the actual-UCB
reference experiment.

Deliverables:

- record the exact committed pre-extraction revision of the current feature
  branch as the behavioural baseline; the upstream `cpp-0.8.3` commit remains
  ancestry and source provenance, not a separate N0 oracle;
- use that revision directly without reconstructing a second historical build or
  maintaining a separate baseline worktree; all native behaviour already
  committed there is accepted as the baseline without a separate deviation
  catalogue;
- retain the reproducible CMake test preset and documented user-space macOS
  development environment;
- keep individual GoogleTest cases discoverable through CTest and development
  tools;
- select a small, deterministic-enough native regression input that exercises
  molecule loading, stable transcript identity, an optional prior, parameter
  resolution, confidence fitting, segmentation, cell products, boundaries, and
  the legacy output serializers;
- record one locked CLI configuration, including its input and output formats,
  prior, scale and initialization behaviour, clustering, iteration and
  convergence settings, and a one-thread OpenMP configuration; use the legacy
  output style so the segmented-molecule CSV retains `transcript_id`;
- run the recorded pre-extraction CLI and retain its parsed scientific result as
  the reference output;
- record the reference commit, resolved parameters, effective OpenMP settings,
  the implicit segmentation seed of `1`, and the locked configuration and
  generation command with that output; compiler and dependency inventories
  belong to CI and release-build provenance rather than this portable fixture;
  and
- add one focused automated semantic comparison that can later evaluate both the
  direct C++ operation and the refactored CLI against the same reference. It does
  not need to be a general-purpose comparison framework.

The reference comparison covers at least:

- retained transcript identifiers and molecule rows;
- cell/noise assignments and molecule confidence;
- resolved scientific options;
- count matrices and cell statistics; and
- boundary geometry.

The comparator must inspect parsed scientific content rather than require raw
output files to be byte-identical. It may normalize harmless cell relabelling,
polygon orientation, and polygon starting vertices, and may use explicitly
justified numerical tolerances. It must not normalize away different molecule
partitions, counts, confidence values, or other scientific changes.

The regression fixture must be small enough for routine native CI. It is not the
deferred actual-UCB reference experiment. Because the pre-extraction CLI exposes
no general segmentation seed, each reference attempt starts in a fresh process
and uses one OpenMP thread. Repeating the small run under those locked settings
is sufficient to confirm semantic stability. The fixture should be chosen so
that its meaningful assignments are stable; any accepted tolerance or semantic
normalization must be recorded with the fixture rather than introduced later to
make a regression pass. Stable transcript identity is validated through the
legacy segmented-molecule CSV; adding it to Parquet remains Native Slice N6.

Exit criterion: a clean checkout can configure, build, and run its native tests,
the small input, locked configuration, reference output, and portable generation
recipe are versioned, the recorded pre-extraction CLI can reproduce the semantic
result under those settings, and one focused automated comparator is ready to
gate N1 through N4.

### Native Slice N1: Define the public segmentation contracts

Implementation status: complete as of 2026-09-03 in commit
`f74417481feb988952eedc22c688e1aefb5e8eb7`.

Define the stable source-level API that the Baysor CLI and the future nanobind
module will share. This slice defines a coherent, compilable native boundary;
it does not yet move the scientific workflow out of `cmd_run(...)` or change an
algorithm. That orchestration extraction belongs to N2.

The current `cmd_run(...)` accepts a mixture of scientific configuration, input
and output paths, CLI presentation choices, and serialization choices. It runs
the complete workflow, writes files immediately, and reduces the outcome to an
integer exit code. Its local lambda named `run_segmentation` is only a
dimension-dispatch closure over captured CLI state; it is not an independently
callable library operation. N1 replaces this implicit interface with an explicit
contract, conceptually:

```cpp
SegmentationOutcome run_segmentation(
    const SegmentationRequest& request,
    const CancellationToken& cancellation
);
```

N1 defines the participating types and the declaration of this operation. N2
implements it by moving the existing orchestration behind the new boundary.

#### Consumer entry point and API exposure

The code-aligned diagram and field-level reference now live in the root-level
[native segmentation API contract](../contract.md). It identifies
`baysor_python.segment(...)` as the future public Python entry point, the
nanobind module as a private adapter, and `baysor::run_segmentation(...)` as the
public native source API. It also enumerates the request, result, cancellation,
and error contracts and separates them from native implementation details.

At the end of N1, the public native types and function declaration exist, but
the function has no scientific implementation yet. N2 makes the native call
operational, N4 changes the CLI to use it, and Python Slices 1C and 1D add the
private binding and public Python entry point respectively. No Python-callable
API is implemented in this repository.

#### Request contract

`SegmentationRequest` is an owning, typed description of one run. It contains:

- a path-oriented molecule-input specification, including column mapping and
  filtering options;
- an optional prior-segmentation specification;
- segmentation, clustering, and confidence/noise options;
- an explicit requested-product selection, defaulting to the CLI-compatible
  full scientific result;
- a 64-bit `random_seed`, defaulting to `1`; and
- execution settings, including the requested native thread configuration.

The request may compose or refine the existing native option structures, but it
must also own the run inputs that are currently passed separately to
`cmd_run(...)`. The supported initial contract remains path-oriented: C++ loads
the prepared CSV or Parquet input. An in-memory Python-array contract is not part
of this slice.

The request records unresolved user intent where resolution requires the loaded
data. For example, an automatically selected scale or initial component count
is not independently calculated by the CLI or Python frontend. The shared
operation resolves it once and returns the effective value in the result.

#### Result and ownership contract

`SegmentationResult` owns the scientific products needed by the native
serializers and the future binding:

- retained molecule data and stable source transcript identity;
- cell/noise assignments, molecule confidence, and assignment confidence;
- optional molecule-cluster assignments;
- stable cell identifiers and cell statistics;
- 2D boundaries and, where applicable, 3D boundary stacks;
- cell-by-gene counts;
- convergence and other requested diagnostic data;
- resolved options and the actual produced-product set; and
- native run provenance, including the resolved seed and effective thread
  settings.

These values use normal C++ ownership and move semantics. The result must not
contain references, pointers, or views into a destroyed `BmmData` or another
temporary working object. It also does not encode output directories, filenames,
console presentation, or process exit policy. Serialization operates separately
over a completed result.

The requested-product field is established here with a backward-compatible
full-result default. N7 later makes selective materialization a measured runtime
and memory optimization, so callers can omit unused boundaries, matrices, or
reports without changing assignments.

#### Cancellation contract

N1 introduces a small C++17-compatible cancellation facility backed by
thread-safe shared state. The public contract separates the two capabilities:

- a `CancellationSource` owns the shared state and may request cancellation;
- a copyable, read-only `CancellationToken` observes that same state and is the
  value passed to `run_segmentation(...)`.

This source/token split lets an embedding caller retain the control capability
on one thread while the segmentation thread receives only the observation
capability. In `baysor-python`, for example, a worker-side registry may retain
the source while the GIL-released native call receives its token. The registry
and Dask control protocol remain outside this repository; the native contract
must not depend on Python, nanobind, or Dask.

The documented contract must state:

- how the source creates and owns the shared state, and the copy and lifetime
  semantics of the token;
- that a token may safely be observed by the segmentation thread while another
  thread uses its source to request cancellation;
- that cancellation is cooperative rather than forced thread termination; and
- that cancellation produces a distinct outcome and never a partially valid
  `SegmentationResult`.

The outcome type must therefore make success and cancellation structurally
distinct, for example through a variant of `SegmentationResult` and a dedicated
`SegmentationCancelled` value. N2 begins using this contract in the extracted
operation. N3 completes the safe major-phase and BMM-iteration checkpoints that
act on it.

#### Randomness and errors

N1 defines one run-local randomness contract. The request supplies a master
64-bit seed, scientific substreams and diagnostic-only substreams are separate,
and any named-substream derivation is documented and versioned. Seed `1` maps to
the N0 one-thread compatibility behaviour. The resolved seed, substream-contract
version, and effective native thread settings are retained in result provenance.
N2 performs the mechanical work of passing this state through every stochastic
scientific operation.

The public API also documents distinct behaviour for invalid requests, molecule
or prior input failures, native-processing failures, serializer/output failures,
and cancellation. The library preserves structured native outcomes; the CLI may
later translate them into messages and exit codes, and nanobind may translate
them into documented Python exceptions.

Deliverables:

- `SegmentationRequest`, containing the molecule-input specification, optional
  prior specification, filtering, segmentation, clustering, confidence/noise,
  requested-product, 64-bit `random_seed` with default `1`, and execution
  options;
- `SegmentationResult`, owning retained molecule identity and data, cell/noise
  assignments, confidence fields, optional clusters, cell statistics,
  boundaries, counts, convergence diagnostics, resolved options, and native
  provenance;
- a thread-safe C++17 cancellation source/token pair with a documented
  ownership, lifetime, and cross-thread use contract;
- a distinct cancelled outcome that cannot be mistaken for a complete result;
- a run-local random-state contract covering every segmentation-affecting source
  of randomness and separating scientific streams from diagnostic-only streams;
- a documented, versioned derivation of any named substreams, including the
  compatibility mapping for the default seed;
- documented invalid-request, input/output, native-processing, and cancellation
  error behaviour; and
- public headers that use ordinary C++ types and do not expose `argv`, CLI
  aliases, process exit policy, presentation strings, Python objects, or internal
  BMM ownership.

Focused contract tests compile and link only against `baysor_lib`. They verify
that a request can be constructed without CLI code, its defaults include seed
`1`, result values have safe owning/move semantics, a source can request
cancellation while another thread observes its token, and a cancelled outcome
cannot be treated as a successful result. They do not duplicate the N0
scientific regression or prematurely test the N2 implementation.

Exit criterion: the request, result, cancellation, error, and ownership contracts
are documented and reviewable as a coherent native API; their focused tests can
compile and link without the CLI executable; and the existing CLI scientific
path remains unchanged pending N2.

### Native Slice N2: Extract the scientific orchestration

Implementation status: complete as of 2026-09-03.

Implement `run_segmentation(request, cancellation)` in `baysor_lib` by moving the
existing orchestration out of `cmd_run(...)`.

#### Starting point and target boundary

At the start of N2, [`src/segmentation/segmentation.cpp`](../src/segmentation/segmentation.cpp)
implements the seed-derivation and structured-error support established by N1,
but it does not define the public run operation. The assembled scientific
workflow still lives inside
[`cmd_run(...)`](../src/cli/main.cpp), where input preparation, parameter
resolution, segmentation, product construction, file writing, logging, and exit
codes are interleaved.

N2 turns the N1 declaration into the directly callable native engine:

```text
SegmentationRequest + CancellationToken
                    |
                    v
        baysor::run_segmentation(...)
                    |
                    v
          SegmentationOutcome
             /             \
            v               v
SegmentationResult   SegmentationCancelled
```

The operation accepts the input path and typed scientific intent through
`SegmentationRequest`. It does not accept an output directory, serialization
format, filename, CLI command, Python object, logger configuration, or process
exit policy. It copies request options into run-local state because the public
request is const, performs canonical data-dependent resolution there, and
returns the effective values through `SegmentationResult::resolved_options`.

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

Concretely, N2 moves or extracts the following orchestration from
`cmd_run(...)`:

1. validate the typed molecule, prior, segmentation, neighborhood-composition,
   requested-product, seed, and execution settings;
2. load and filter the prepared CSV or Parquet molecule table and load the
   optional prior;
3. infer prior-derived scale, initial cell count, and other values that require
   the loaded molecule cloud;
4. estimate molecule confidence/noise, build the molecule graph, and perform
   optional MRF, Louvain, or Leiden molecule clustering;
5. dispatch to the existing `initialize_bmm_data<2>`/`<3>` and `bmm<2>`/`<3>`
   implementations according to the loaded dimensionality; and
6. derive the requested scientific products and move their durable values into
   one owning `SegmentationResult`.

The existing loaders, prior code, confidence estimator, graph construction,
clustering implementations, BMM, boundary estimator, and scientific calculations
remain authoritative. N2 changes their coordination and ownership, not their
mathematics, defaults, or update order.

#### Result materialization and serialization boundary

The current CLI calculates cell statistics, boundaries, and the count matrix
next to the code that immediately writes those products. N2 separates
calculation from representation on disk:

```text
BmmData and temporary working objects
                    |
                    | copy or move durable scientific values
                    v
          owning SegmentationResult
                    |
                    | separate reusable serializers
                    v
       CSV / Parquet / geometry / matrix artifacts
```

The returned result contains retained molecule data, cell/noise assignments,
assignment confidence, optional molecule clusters and neighborhood-composition
colours, stable cell identifiers, cell statistics, 2D or 3D boundaries, the
cell-by-gene count matrix, report-neutral diagnostics, resolved options,
produced-product flags, and native provenance. No pointer, reference, or view
into `BmmData`, a clustering model, a report object, or CLI-owned state may
escape the operation.

Native serializers become reusable functions over a completed result. The
scientific operation itself writes no output artifacts. HTML generation and
other presentation code consume report-neutral result data outside the run
operation. The default request still produces the complete scientific result;
N7, rather than N2, is responsible for proving and optimizing selective
materialization.

#### Run-local randomness

N2 must replace the current mixture of process-global and hard-coded random
state with one run-local context derived from `SegmentationRequest::random_seed`.
The relevant paths include duplicate-point jitter, stochastic BMM assignment,
molecule clustering, and neighborhood-composition or diagnostic embeddings.

The default seed `1` preserves the N0 single-thread scientific stream. Other
master seeds use the versioned substream derivation established by N1. Scientific
and diagnostic-only streams remain separate, so requesting a diagnostic or
colour product cannot consume random values that would otherwise affect cell
assignments. Any internal function-signature changes needed to pass a generator
or derived seed are part of N2; the random source must not become another public
Python or CLI contract.

The result records the master seed, derived substreams, substream-contract
version, and effective native thread setting in provenance. Same-seed
single-thread repeatability, repeated-call isolation, and multithreaded
characterization are hardened and tested in N3.

#### Errors, cancellation, and CLI overlap

Scientific validation or execution failures that are currently reduced to a log
message and integer return code become categorized `SegmentationError`
exceptions. Cooperative cancellation remains the distinct
`SegmentationCancelled` outcome and must never publish a partially valid result.
Comprehensive cancellation checkpoints, repeated-call lifecycle behaviour, and
embeddability hardening belong to N3.

N2 establishes the authoritative direct library path without claiming that the
CLI refactor is complete. The existing N0 CLI regression must remain green while
the operation is extracted. N4 then reduces `cmd_run(...)` to a frontend over
the operation, removes any remaining duplicate orchestration, and performs the
full pre-extraction CLI-parity gate.

Deliverables:

- move scientific orchestration and canonical option resolution into the public
  operation without changing the underlying algorithms;
- initialize run-local random state from `SegmentationRequest::random_seed` and
  pass it explicitly to duplicate-point jitter, stochastic BMM assignment, and
  every clustering step that can affect segmentation;
- preserve the N0 one-thread stream for seed `1`, including existing fixed
  scientific subsystem seeds where required by parity, while defining how other
  master seeds deterministically derive those subsystem streams;
- inventory fixed diagnostic seeds and ensure report or colour generation cannot
  perturb the segmentation stream;
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

#### Verification and exit criterion

A focused C++ test uses the N0 fixture to construct a `SegmentationRequest`,
invokes a complete segmentation through `baysor::run_segmentation(...)`, and
inspects the resulting scientific values without invoking the CLI or Python. It
checks at least the retained molecule identity, cell/noise assignments,
confidence, resolved options, cell identifiers and statistics, boundaries,
count matrix, produced-product flags, seed, and effective-thread provenance.
Where the N0 reference is applicable, comparison is semantic rather than based
on byte-identical serialization.

Exit criterion: the direct native test passes, the returned data remains valid
after all temporary algorithm state has been destroyed, the existing N0 CLI
regression remains green, and no algorithm or scientific default has changed.
Cancellation-depth, repeated-call, embeddable-consumer, and final CLI-parity
gates remain explicitly assigned to N3 and N4.

### Native Slice N2b: Harden extraction with a scientific parity matrix

Implementation status: complete in the current working tree as of 2026-09-04.

The implemented fixture set and historical generation recipe are stored in
`tests/fixtures/scientific_parity`. The reference scientific artifacts were
generated twice from the pre-extraction commit and were byte-identical between
runs; only the machine-specific invocation comment in the resolved-parameter
dump was removed. CTest discovers separate current-CLI and direct-operation
cases for every matrix entry, so a failure identifies both the scientific branch
and the frontend that diverged.

N2b tests the algorithm branches whose random-state plumbing changed during N2
but which are not exercised by the original N0 fixture. It is a focused
equivalence gate, not an expansion of the segmentation algorithm or a general
benchmark suite.

N2b must protect against two different failure modes:

1. the new direct operation may orchestrate unchanged scientific components
   differently from the legacy CLI; and
2. the seed-plumbing changes inside shared scientific components may alter both
   the current CLI and the direct operation in the same way.

A comparison between only the current CLI and the direct operation detects the
first failure mode but cannot exclude the second, because both link the same
modified `baysor_lib`. The scientific oracle for every end-to-end N2b case is
therefore a recorded result from the pre-extraction behavioural baseline,
commit `f46a1e1dce1606d0ea644f4f8f1cf682597ba65c`:

```text
                 recorded pre-N2 result
                    from f46a1e1...
                       /         \
                      v           v
          current legacy CLI     direct native operation
```

Generate each reference twice from fresh single-threaded processes with
`OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, and the pre-extraction CLI's implicit
seed `1`. Retain the small input, locked configuration, parsed scientific
outputs, baseline revision, and generation recipe. Routine CI consumes these
versioned references; it does not rebuild the historical checkout. Both the
current legacy CLI and `run_segmentation(...)` must independently match the same
reference. Comparing the two current paths remains a useful diagnostic, but it
is not the historical oracle.

Add the following small parity matrix:

1. a 2D fixture containing duplicate coordinates, exercising the jittered graph
   path through the pre-N2 reference, current legacy CLI, and direct native
   operation;
2. a fixture exercising MRF molecule clustering through all three paths;
3. small parameterized Louvain and Leiden cases covering their seeded NCV basis,
   clustering dispatch, molecule-cluster assignments, and downstream
   segmentation through all three paths;
4. a small 3D case exercising dimensional dispatch, segmentation, joined
   boundaries, and per-Z boundary stacks through all three paths; and
5. a focused unit test proving that a legacy default call and the equivalent
   explicitly supplied legacy seed produce identical random sequences and
   scientific results.

#### NCV scientific and representation boundary

N2b distinguishes the neighborhood-composition features used by the scientific
workflow from the hexadecimal colours derived from them:

```text
neighborhood-composition vectors
            |
            +-- Louvain/Leiden clustering -> molecule types -> BMM
            |
            +-- UMAP/interpolation -> Lab/RGB -> per-molecule hex colours
```

The first path is scientific. Louvain and Leiden molecule-cluster assignments
and their downstream cell/noise partition, confidence, counts, statistics, and
boundaries must match the pre-N2 reference. MRF clustering does not consume the
NCV representation and is tested independently.

The second path is a post-segmentation visualization product. Its absolute UMAP
orientation and resulting `#RRGGBB` values have no stable biological identity
and may be sensitive to numerical library, compiler, or platform differences.
Exact NCV colour strings and UMAP coordinates are therefore excluded from the
historical scientific parity oracle.

Exclusion from the historical oracle does not remove product coverage. Focused
tests must verify that requesting NCV colours returns one well-formed colour per
retained molecule, is repeatable for the same build, input, seed, and one-thread
configuration, and does not change molecule clusters, cell/noise assignments,
confidence, counts, statistics, or boundaries compared with an otherwise
identical request that omits colours. These tests validate the colour-product
contract and the separation of diagnostic randomness without treating an
absolute colour as a scientific result.

Reference comparisons remain semantic: molecule identity, noise/cell partition,
confidence, molecule clustering where applicable, statistics, boundaries,
counts, and resolved settings are compared without requiring byte-identical
serialization or identical arbitrary cell-label numbering. The existing N0
comparators should be reused or minimally generalized; N2b does not create a
general-purpose comparison framework. Use seed `1` and one native thread for the
compatibility comparisons. Tests of other seeds should assert repeatability and
invariants rather than equality to a historical CLI result that had no
configurable master seed.

Exit criterion: every new pre-N2 reference is reproducible under its locked
settings, both the current legacy CLI and direct operation match it semantically,
the supporting legacy-seed test passes, and the existing N0 regression remains
green. Review finds no unexplained scientific difference in the covered 2D, 3D,
duplicate-coordinate, MRF, Louvain, Leiden, or NCV-based clustering paths; the
NCV colour product also passes its shape, repeatability, and scientific
non-interference tests. Any difference must be explained and accepted before N2
is merged; it must not be hidden by broad numerical tolerances or fixture
regeneration.

### Native Slice N3: Complete cancellation, lifecycle, and embeddability

Implementation status: complete as of 2026-09-04.

Harden the extracted operation for repeated in-process use and hand it off as a
clean CMake consumer target.

#### Thread-resource and concurrency contract

N3 must define the thread setting precisely before an embedding layer relies on
it. `SegmentationExecutionOptions::native_threads` is a call-scoped
request for Baysor's OpenMP team-size limit; it is not a process-wide limit on
every thread that native dependencies may create. A value of zero inherits the
already configured OpenMP runtime policy. A positive value applies only for the
segmentation call, and the caller's prior OpenMP setting must be restored on
success, cancellation, and exception paths.

`omp_get_max_threads()` reports the configured maximum for a future OpenMP team,
not the number of threads that a parallel region actually used. OpenMP dynamic
team adjustment may use fewer threads. N3 must therefore make provenance
unambiguous: retain the requested value, record the configured OpenMP maximum
and dynamic-adjustment setting, and do not describe either as an observed thread
count unless a parallel region is explicitly instrumented. The library should
inherit `OMP_DYNAMIC` rather than silently changing that process policy.

Arrow CSV and Parquet loading can use Arrow's own executor independently of
OpenMP. N3 must inventory this path and either expose a separate, deliberately
named I/O-thread policy or document and test the selected fixed behaviour. The
OpenMP request and its provenance must not imply that Arrow threads are included
in the same budget.

The supported native entry context is a top-level caller that is not already
inside an active OpenMP parallel region. Calling Baysor from inside such a region
can invoke nested-parallelism and active-level policies that change or serialize
the requested teams. N3 must detect this condition before applying a thread
override, fail with a documented structured error, and test that behaviour. This
restriction does not affect the normal Python/Dask worker entry path or the
OpenMP regions created internally by Baysor.

Sequential repeated calls in one process are supported. Concurrent segmentation
calls in one process remain unsupported unless N3 proves that their OpenMP
configuration, random state, logging, and dependencies are isolated. This
limitation must be explicit. The corresponding `baysor-python` execution model
permits at most one active Baysor call per worker process: one tile uses the
worker's allocated OpenMP budget, while multiple tiles run concurrently only in
separate, appropriately budgeted worker processes. This avoids oversubscription
from multiplying the worker count by the OpenMP team size without preventing
later benchmark-driven parallelism.

The execution contract deliberately supports three configurations:

1. **OpenMP-focused:** one worker processes tiles sequentially and assigns many
   OpenMP threads to the active tile;
2. **tile-focused:** multiple worker processes segment different tiles
   concurrently with one OpenMP thread per tile; and
3. **hybrid:** multiple worker processes segment different tiles concurrently,
   with a bounded multi-threaded OpenMP team in each worker.

For CPU compute threads, the coordinator must maintain the budget

```text
concurrent Baysor worker processes * OpenMP threads per worker
    <= allocated CPU cores
```

For example, a 16-core allocation may use `1 * 16`, `16 * 1`, `4 * 4`, or
`8 * 2`, subject to available memory and measured performance. These are all
supported configurations; benchmarks determine the appropriate default for a
dataset and deployment. Concurrent tiles multiply the per-tile working-memory
requirement, and Arrow I/O threads are accounted for separately from the formula
above.

Cancellation remains cooperative. Python, Dask, or another embedding caller may
request cancellation from a control thread, but the native call returns only
after it reaches a safe checkpoint; the token does not forcibly terminate an
OpenMP team. Hard termination and complete allocator-memory reclamation remain
worker-supervisor policies based on process isolation, outside this native API.

N3 implements only the native half of that cancellation path: the shared atomic
state, source/token ownership, token propagation, safe polling checkpoints, and
distinct cancelled outcome. Baysor C++ does not listen to Python, Dask, operating
system signals, or a Dask `Future` directly. In particular, calling
`Future.cancel()` does not by itself invoke
`CancellationSource::request_cancellation()` or set the token observed by the
running native operation.

`baysor-python` must provide the other half. Its nanobind adapter retains the
`CancellationSource` in a worker-side active-attempt registry, releases the GIL
around `run_segmentation(...)`, and exposes a control path that reaches the same
worker process and invokes `request_cancellation()` from a separate control
thread. The Dask coordinator must use that explicit bridge and wait for the
distinct cancelled outcome for a bounded grace period. If the request cannot
reach the worker or Baysor does not reach a safe checkpoint in time, the Nanny
terminates and replaces the dedicated worker process. Implementing and testing
that registry, control protocol, timeout, and supervisor fallback belongs to
`baysor-python`, not N3 in this repository.

Deliverables:

- add cooperative cancellation checks after major phases and at safe BMM
  iteration boundaries;
- ensure cancellation publishes no complete result and does not run completed
  result serializers;
- test invalid requests and structured native failures;
- test result ownership after temporary working state has been destroyed;
- test repeated same-process calls for hidden global-state leakage;
- test that repeated same-process calls with the same seed and one thread produce
  the same semantic result, without one call advancing another call's stream;
- define `native_threads` as an OpenMP request rather than a total native-process
  thread cap, and document the supported same-process concurrency model;
- preserve an inherited OpenMP configuration when `native_threads` is zero and
  restore the prior configuration after positive requests on successful,
  cancelled, and exceptional exits;
- replace or precisely document the misleading `effective_native_threads`
  terminology: provenance must distinguish the requested value, configured
  OpenMP maximum, and dynamic-team setting from any actually observed team size;
- inventory Arrow CSV and Parquet reader threading and make its relationship to
  the OpenMP budget explicit, without implying that one setting controls both;
- require a top-level native entry context and reject a call made from an active
  OpenMP parallel region before changing its thread configuration;
- add focused tests for zero/inherited configuration, positive thread requests,
  restoration across every exit path, nested-entry rejection, and thread
  provenance;
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
owned-result lifetime, repeated calls, and OpenMP setting/restoration semantics;
thread provenance cannot be mistaken for an observed team size; Arrow I/O
threading, the top-level entry requirement, and the supported same-process
concurrency model are documented; and a clean CMake consumer can link the
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
- expose `--seed` as the CLI spelling of the shared 64-bit `random_seed`, default
  it to `1`, and record its resolved value in run provenance;
- compare the refactored CLI with the recorded pre-extraction CLI on the N0
  fixture using seed `1`, the same inputs, options, one-thread OpenMP
  configuration, and output mode; and
- retain reference and candidate logs and resolved options as parity evidence.

Parity includes retained transcript identity, cell/noise assignments, molecule
confidence, count matrices, cell statistics, resolved parameters, and Baysor
revision. Cell identifiers and polygons are compared semantically after
normalizing harmless relabelling, polygon orientation, and starting-vertex
differences. Because the pre-extraction CLI has no general segmentation seed,
parity runs compare its implicit seed-`1` stream with the refactored CLI's
explicit `--seed 1` stream in one-thread mode. Multi-threaded runs remain a
separate repeatability characterization and must not be assumed bitwise
deterministic merely because they use the same seed.

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
development. The complete native suite and pre-extraction CLI parity fixture are
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
- same-seed one-thread repeatability in fresh and repeated same-process calls;
- proof that diagnostic-product selection does not change scientific assignments
  by consuming the run's scientific random stream;
- resolved-seed and effective-thread provenance;
- serializer consistency between direct and CLI paths;
- semantic parity with the recorded pre-extraction CLI;
- boundary edge cases, sparse cell identifiers, target-cell batches, and
  thread-count stability;
- stable transcript-ID round trips through filtering and Parquet;
- selective-output equivalence for authoritative assignments; and
- confidence/noise calibration fit/export/apply consistency.

Tests must use an explicit random seed and explicit OpenMP settings and record
both in parity evidence.
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
