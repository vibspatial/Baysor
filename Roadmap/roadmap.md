# Scalable Baysor integration roadmap

Date: 2026-08-28

## Executive summary

The separate `baysor_python` package should integrate the current native C++
Baysor as a points-first workflow with two execution modes:

1. an untiled mode that establishes correctness and provides a reference
   result; and
2. an optional tiled mode for datasets or deployment environments that cannot
   run Baysor comfortably as one process.

Both modes should share the same input preparation, `baysor_python` segmentation
runner, output contract, provenance model, and SpatialData importer. Tiling should
be an execution strategy within this integration rather than a separate raster
segmentation implementation.

The benchmark in [benchmark.md](benchmark.md) showed that an approximately
47-million-transcript synthetic CosMx mosaic completed a one-iteration resource
run with 19.81 GB peak resident memory. This does not establish segmentation
quality or the runtime of a converged run, but it demonstrates that the UCB-sized
sample does not require tiling solely to fit on the available 32 GB machine.
Consequently, the production reference experiment should be untiled on the
actual sample. The untiled result is also needed as the reference against which
the tiled implementation will be validated.

The tiled design must reconcile molecule-to-cell assignments using stable
transcript identities in overlapping tiles. It should not rasterize every tile
and rely primarily on pixel-IoU stitching. Baysor's molecule assignments are the
authoritative result; cell polygons, count matrices, cell statistics, and
optional raster labels should be derived globally after reconciliation.

The separate `baysor-python` repository will own the pinned Baysor source and
native build. Its segmentation and boundary APIs will call a shared C++ library
through a thin native Python extension. It will also own the generic
Python-orchestrated core-plus-halo tiling, its Dask Distributed execution layer,
assignment reconciliation, rescue logic, and a thin SpatialData adapter. Harpy
may provide a convenience API and pass a caller-managed Dask client, but will own
neither Baysor-derived C++ code nor Baysor-specific execution or tiling semantics.

## Goals

The integration should:

- run a pinned C++ Baysor implementation through a versioned `baysor_python`
  package without building or vendoring Baysor inside downstream packages;
- accept a SpatialData points element and an optional labels element as a
  transcript-native prior;
- support large backed SpatialData objects without materializing the complete
  transcript table in pandas;
- preserve one stable identity for every retained transcript;
- produce globally consistent molecule assignments, cell identifiers, shapes,
  counts, and cell statistics;
- reproduce Baysor's boundary semantics through the `baysor_python` native API,
  with explicit upstream provenance and parity tests;
- expose reusable Python-orchestrated tiled segmentation from `baysor_python`,
  independently of Harpy, using Dask Distributed for execution while keeping the
  scientific tile and result contracts independent of Dask;
- expose all scale-sensitive and resource-sensitive Baysor parameters;
- make runs reproducible, inspectable, resumable, and safe to retry; and
- demonstrate that tiled results are not materially dependent on tile seams or
  grid placement.

The supported scope through the tiled-validation milestone is 2D. Arbitrary
affine transformations, 3D Baysor, and automatic distributed deployment are
outside that supported scope and require separate production design and
validation before release.

### Product quality policy

This roadmap describes staged delivery of a professional product, not an MVP or
a sequence of disposable prototypes. Every completed phase must meet production
engineering standards within its declared scope: stable public contracts,
reproducible and offline-capable release builds, supported error behaviour,
provenance, focused automated tests, documentation, licence compliance, and a
defined upgrade path.

Phasing limits which capabilities are supported; it does not lower their quality
bar. A capability that has not yet met its scientific acceptance gates may be
available only as an explicitly non-production validation mode. It must not be a
default, be described as production-ready, or silently weaken a quality gate.

## Why this should be a new points-first integration

Harpy's existing `baysor_callable` targets the older Julia CLI contract. It
writes per-chunk CSV and TIFF inputs, sets `JULIA_NUM_THREADS`, reads GeoJSON,
and immediately rasterizes the polygons. The generic `segment_points` path then
reconciles independently generated raster labels across Dask chunks.

That design is not a good fit for current Baysor C++:

- Baysor already performs its expensive work in native C++ and exposes a CLI;
- it accepts and emits Parquet;
- its main result is a molecule-to-cell assignment with confidence and noise
  information;
- its Parquet bundle includes cell statistics and GeoParquet boundaries; and
- converting each local result to a raster before stitching loses the strongest
  cross-tile evidence: the assignments of the same transcripts in both tiles.

A dedicated high-level SpatialData operation in `baysor_python` should present
the complete workflow. A convenience operation such as
`harpy.pt.segment_baysor` may delegate to it without reimplementing preparation,
execution, or import. The exact public names can be settled during API review.
Harpy's existing Julia-era callable can be deprecated independently after the
new integration is established.

## `baysor_python` integration boundary

`baysor-python` will be a separately versioned repository and Python
distribution, exposing the import package `baysor_python`. It will own:

- the pinned Baysor source and native dependency build;
- one coarse-grained public C++ segmentation operation shared by the CLI and
  Python binding;
- `segment(...)`, implemented over a thin native binding to that operation;
- `segment_tiled(...)`, implementing core-plus-halo planning, staging,
  Dask-backed tile execution, assignment reconciliation, and rescue;
- the supported local Dask cluster configuration and the contract for accepting
  a caller-managed Dask client;
- `boundaries(...)`, implemented as a direct array-oriented native binding;
- a thin SpatialData adapter for input preparation, prior sampling, stable
  transcript identity, and construction of coordinated output elements;
- translation of native failures into documented Python exceptions and result
  objects; and
- Baysor-version reporting, build provenance, wheels, and CLI-versus-Python
  parity tests.

Downstream packages such as Harpy may own:

- a thin convenience API over the public `baysor_python` SpatialData adapter;
- Harpy-specific defaults, metadata, visualization, and optional rasterization;
- application-specific Dask cluster provisioning and passing an existing client
  into the `baysor_python` execution contract; and
- end-to-end scientific and seam validation.

The dependency direction is intentionally one-way:

```text
Harpy or another application -> baysor-python -> SpatialData
                                      |
                                      +-> Dask Distributed (`tiled` extra only)
```

`baysor-python` must not import or depend on Harpy. Its only direct Python runtime
dependency in the Phase 1 release is SpatialData. The low-level `segment(...)`
and `boundaries(...)` contracts should nevertheless remain independent of
SpatialData so they can be tested and reused without the high-level adapter.
Dask Distributed enters as an optional `baysor-python` dependency for the tiled
execution phases, for example through a `tiled` package extra. Importing and
using the Phase 1 untiled API must not require that extra.

The authoritative scientific workflow must live in one C++ library operation,
conceptually `run_segmentation(request) -> result`. The existing CLI becomes a
thin frontend that parses arguments, calls that operation, and serializes the
result. The native Python extension calls the same operation, releases the GIL
during the long-running calculation, and translates C++ failures into documented
Python exceptions. Neither frontend may duplicate parameter resolution or
scientific orchestration.

The request includes one 64-bit native `random_seed`, defaulting to `1`. The
shared operation owns run-local random state and returns the resolved seed in
provenance. The refactored CLI exposes it as `--seed`; the Python API forwards
the same field and must not maintain a separate Python random stream for native
segmentation. Enabling reporting or other optional diagnostics must not consume
the scientific stream and change assignments.

```text
                         C++ run_segmentation(...)
                                  |
                    +-------------+-------------+
                    |                           |
                    v                           v
             Baysor CLI                  nanobind module
        parse CLI arguments             convert Python input
        call shared function            release Python GIL
        serialize outputs               call shared function
                                        translate exceptions
                                        return structured result
```

Both arrows therefore lead to the same compiled scientific operation. The CLI
and Python module are adapters around that operation, not separate segmentation
implementations.

The direct native call is the low-level backend for untiled segmentation. Tiled
execution invokes that same extension inside Dask worker processes supervised by
Dask Nannies. Worker isolation and recycling are execution policies, not reasons
to make the CLI the Python API backend.

During development, the pinned Baysor source may be included as a Git submodule.
Published source distributions must bundle the exact source snapshot and must not
require Git, submodule initialization, or network access during installation.

### Incorporating the pinned Baysor source and native extension

Phase 1 will incorporate the native C++ Baysor implementation without porting its
segmentation algorithm to Python. The upstream source-provenance baseline is tag
[`cpp-0.8.3`](https://github.com/kharchenkolab/Baysor/tree/cpp-0.8.3), commit
`d7077a7ded6f4b941915badc894f767532d39fd2`. This is an immutable revision; builds
must not follow a moving upstream branch.

Native integration development is owned by the maintained
[`vibspatial/Baysor`](https://github.com/vibspatial/Baysor) fork. Its `cpp`
branch preserves the immutable upstream baseline. Work must be performed on a
dedicated feature branch, such as `feature/segmentation-api`, rather than by
modifying that baseline branch. The original `kharchenkolab/Baysor` revision is
the upstream ancestry and source-provenance baseline. N0 records the current
feature branch's exact pre-extraction revision as its behavioural baseline; a
later reviewed commit in the `vibspatial/Baysor` fork is the integration revision
consumed by `baysor-python`.

The development checkout should use this layout:

```text
baysor-python/
├── CMakeLists.txt
├── pyproject.toml
├── vendor/
│   └── Baysor/                 # vibspatial/Baysor at an exact integration commit
├── licenses/
│   └── Baysor-LICENSE          # upstream MIT notice
├── src/baysor_python/
│   ├── _baysor_core.cpp        # thin nanobind module
│   ├── _segment.py             # public Python segmentation API
│   ├── _models.py              # request and result types
│   ├── _exceptions.py
│   └── spatialdata.py          # high-level SpatialData adapter
└── tests/
```

The submodule is a developer source-management mechanism, not a runtime fetch.
It points to `https://github.com/vibspatial/Baysor.git`, while its commit history
retains the relationship to the immutable upstream baseline. The fork carries
the C++ entry-point extraction as a small, reviewable series of commits that can
be proposed upstream. The submodule must pin the exact reviewed integration
commit; the build manifest must record both the upstream baseline and the fork
integration commit. Native changes must first be committed and tested in
`vibspatial/Baysor`; `baysor-python` must never carry a dirty submodule or local
patches as an implicit part of its build.

Release source distributions must contain an exported copy of that exact tree.
Baysor's CMake currently obtains several pinned header-only dependencies through
`FetchContent`; an offline-buildable source distribution must export those exact
dependency sources too and direct CMake to the local copies. Neither wheel nor
source-distribution installation may clone repositories at build time.

### Repository ownership and integration handoff

The repository boundary follows the product architecture:

The native implementation slices owned by this repository are maintained in
[roadmap_Baysor.md](roadmap_Baysor.md). This section retains the cross-repository
ownership and handoff contract.

> `vibspatial/Baysor` owns the reusable segmentation engine;
> `baysor-python` owns the Python binding, packaging, and higher-level Python
> integration.

The responsibilities are:

| Responsibility | Owning repository |
| --- | --- |
| Preserve the upstream baseline and native regression fixture | `vibspatial/Baysor` |
| Define and implement `SegmentationRequest`, `SegmentationResult`, cancellation, and `run_segmentation(...)` | `vibspatial/Baysor` |
| Extract the scientific pipeline from the CLI and keep the CLI as an adapter | `vibspatial/Baysor` |
| Provide a reusable CMake library target suitable for embedding | `vibspatial/Baysor` |
| Test native results, errors, ownership, repeated calls, cancellation, and CLI parity | `vibspatial/Baysor` |
| Pin the reviewed fork commit and record source/build provenance | `baysor-python` |
| Build, package, and test the nanobind extension | `baysor-python` |
| Provide the public Python, SpatialData, and later Dask APIs | `baysor-python` |

The fork must remain independent of nanobind, Python packaging, SpatialData, and
Dask. Its public operation is a normal source-level C++ API with no Python types.
Conversely, `baysor-python` consumes that API rather than duplicating scientific
orchestration or applying private patches to the vendored checkout.

The handoff is commit based. Native work is reviewed and tested in
`vibspatial/Baysor`; a green integration commit is selected; and
`baysor-python` advances `vendor/Baysor` to that exact commit. Any later native
change repeats the same sequence. This makes a submodule-pointer update an
explicit dependency update with auditable provenance.

The implementation sequence across the repositories is:

1. preserve the fork's `cpp` baseline and create the native feature branch and
   regression fixture;
2. define, extract, and test the coarse-grained C++ operation in
   `vibspatial/Baysor`;
3. select the reviewed Slice 1A native commit, pin it in `baysor-python`, and
   pass the consumer build smoke check;
4. refactor and parity-test the CLI in `vibspatial/Baysor`, then advance the
   `baysor-python` submodule to that reviewed Slice 1B commit; and
5. build the nanobind extension in `baysor-python` against the final pinned
   native contract.

### Packaging and installation model

The current pure-Python Hatchling configuration supports only the repository's
pre-native packaging state. Phase 1 native packaging replaces it with
`scikit-build-core`, which invokes the repository's top-level CMake project. CMake
will include the pinned Baysor project, build its existing C++17 `baysor_lib` as
position-independent code, build a small nanobind extension such as
`baysor_python._baysor_core`, and link the library implementation into that
extension. The CMake install rule places the extension in the `baysor_python`
package. Importing it loads compiled machine code into Python; it does not launch
an executable or compile Baysor at segmentation time.

The source and release flow is:

```text
development checkout
vendor/Baysor/ Git submodule at an exact integration commit
             |
             | export checked-out files, licences, and pinned source dependencies
             v
published source distribution (.tar.gz)
ordinary source files; no nested Git repository and no runtime submodule fetch
             |
             | scikit-build-core -> CMake -> C++ compiler in release CI
             v
platform wheel (.whl)
Python modules + _baysor_core native extension + required redistributable libraries
             |
             | pip selects and unpacks the matching wheel
             v
installed package
no Git, submodule initialization, CMake, compiler, or build-time network access
```

The artifacts have distinct responsibilities:

- **Developer checkout:** `vendor/Baysor` is an active Git submodule. Developers
  initialize it and CMake compiles directly from that pinned checkout.
- **Source distribution:** the release process materializes the checked-out
  Baysor files as ordinary files under `vendor/Baysor`; submodule metadata is not
  part of the archive. It also contains the nanobind source, CMake files,
  licences, and exact sources for dependencies that Baysor otherwise obtains via
  `FetchContent`. Building this archive must not clone a repository or select a
  moving dependency revision.
- **Platform wheel:** release CI compiles ahead of time. The wheel contains the
  Python code, the importable native extension, the build manifest, licence
  notices, and every required redistributable shared library not guaranteed by
  the platform. It need not contain the Baysor C++ source or a standalone Baysor
  executable.
- **Reference CLI:** the `baysor` executable remains a build target and parity
  oracle over the same `run_segmentation(...)` operation. It is built and tested
  in development and CI, but is not the backend of the Python API and is not a
  required wheel payload. It may be distributed separately. If a Python-facing
  command is later useful, its entry point should call the supported Python API.

Self-contained source distribution means that all project-controlled and
`FetchContent` source is present and the project does not download code while it
builds. It does not embed a compiler or every system development library. A
source build may still require the native prerequisites documented in Baysor's
[installation guide](https://github.com/kharchenkolab/Baysor/blob/cpp/docs/installation.md).
Normal supported installations must select a prebuilt wheel and therefore must
not require those build tools.

Phase 1 supports ordinary, GIL-enabled CPython 3.12 through 3.14. CPython is the
interpreter and binary-API boundary for the native extension; it is not another
implementation layer beside nanobind. Nanobind supplies the C++ binding code that
calls the CPython C API and can compile that code against either a version-specific
CPython ABI or CPython's Limited API and Stable ABI (`abi3`). It does not require
one wheel per CPython minor version.

The preferred release artifact is an `abi3` extension with a Python 3.12 ABI
floor, configured through nanobind's stable-ABI support and the corresponding
`scikit-build-core` wheel tag. One such wheel per platform and architecture can
serve the complete declared CPython 3.12--3.14 range. This reduces artifact count,
not validation scope: the same finished wheel must still be installed and tested
under every supported CPython minor version.

The binding's coarse request/result boundary makes stable-ABI dispatch overhead
unlikely to be material relative to a complete Baysor run. Phase 1 must verify
that assumption and all required nanobind features, including GIL release,
cancellation, exception translation, and array conversion. Conventional
version-specific CPython wheels remain the supported fallback only if the
stable-ABI build has a measured compatibility, correctness, or material
performance problem. The reason for taking that fallback must be recorded in the
release documentation.

PyPy and free-threaded CPython builds are outside the Phase 1 support matrix;
either would require its own compatibility and native-wheel validation. Nanobind
split mode is likewise not selected: this package has one coarse native extension,
so introducing a separate compiled nanobind backend dependency is not justified
without a measured need.

This is an established open-source integration pattern rather than a
project-specific packaging invention. The closest precedent is
[PyValhalla](https://github.com/valhalla/valhalla/blob/master/src/bindings/python/CMakeLists.txt),
which links a substantial C++ routing engine into a nanobind `STABLE_ABI` module
and publishes one Python-3.12-floor `abi3` wheel per platform. Comparable
production uses include
[Tesseract Robotics](https://github.com/tesseract-robotics/tesseract_nanobind/blob/main/CMakeLists.txt),
[ONNX Optimizer](https://github.com/onnx/optimizer/blob/main/CMakeLists.txt), and
[SHAP](https://github.com/shap/shap/blob/master/CMakeLists.txt). These precedents
validate the general architecture, but they do not replace Baysor-specific tests
for native dependencies, Arrow/PyArrow coexistence, cancellation, or repeated
in-process execution.

The intended build configuration is conceptually:

```text
nanobind_add_module(_baysor_core NB_STATIC STABLE_ABI ...)
target_link_libraries(_baysor_core PRIVATE baysor_lib ...)

[tool.scikit-build]
wheel.py-api = "cp312"

result: baysor_python-<version>-cp312-abi3-<platform>.whl
```

`STABLE_ABI` constrains the compiled extension to CPython's Limited API;
`wheel.py-api = "cp312"` declares the corresponding compatibility floor in the
wheel tag. Both sides must agree. An `abi3` filename alone is not evidence that
the binary is compliant, so every repaired wheel must pass
[`abi3audit --strict`](https://github.com/pypa/abi3audit). This ABI audit is
separate from `auditwheel` or `delocate`, which inspect and repair platform-native
library dependencies.

`cibuildwheel` should drive the release matrix. Linux wheels must be repaired and
inspected with `auditwheel`; macOS wheels must likewise be repaired and inspected
with `delocate`. The release workflow must run the strict Stable-ABI audit after
platform repair and before clean-environment tests. Those tests must install the
finished wheel itself, not a working tree or an editable build.

Native dependency packaging is part of the product contract. Header-only code is
compiled into the extension. Baysor library code should be linked into the
extension rather than installed as a separately discoverable public Baysor
library. Other native libraries should be linked statically when technically and
licence compatible; eligible shared libraries must otherwise be bundled with
correct run paths. Libraries guaranteed by the wheel's platform policy remain
system provided.

Arrow and Parquet require an explicit packaging decision because SpatialData also
uses PyArrow. The preferred outcome is a self-contained wheel with no dependency
on an arbitrary system Arrow installation, but Phase 1 must compare a private
bundled Arrow/Parquet build with linking against the libraries supplied by a
strictly compatible PyArrow release. The chosen design must pass same-process
SpatialData/PyArrow and Baysor segmentation tests, native dependency inspection,
clean-install tests, and a wheel-size assessment. The public path-oriented API
must not exchange Arrow C++ objects across the extension boundary; this limits
ABI coupling whichever packaging strategy is selected.

The package must include a machine-readable build manifest containing at least:

- the `baysor-python` version;
- the upstream Baysor baseline tag and full commit;
- the `vibspatial/Baysor` integration commit;
- the native extension checksum;
- target platform and architecture;
- compiler and relevant native dependency versions; and
- an inventory of the maintained integration commits relative to the immutable
  upstream baseline.

The production Python API does not select an arbitrary system Baysor executable;
it uses the C++ revision compiled into its extension. The CLI is retained for
parity tests, diagnostics, and standalone Baysor use. Any necessary Baysor source
change should be an auditable commit in `vibspatial/Baysor` and proposed
upstream; release builds must never apply an undocumented local edit.

The binding must be deliberately coarse-grained. It should expose the shared run
operation and the boundary operation, not Baysor's internal BMM classes or full
object graph. The Phase 1 binding may remain path-oriented: C++ reads the prepared
Parquet input, runs segmentation, uses shared serializers, and returns structured
output paths and provenance. An in-memory array API should be added only when
profiling demonstrates that it provides a material benefit and its ownership and
copying semantics are defined.

### Why a native binding, and why nanobind rather than Cython

Cython could call Baysor's C++ code, but it would not remove the need to first
define a stable public C++ operation. It would add a generated-C++ interface
layer around essentially the same coarse request/result boundary. Cython remains
a valid option when substantial Python-like code is being compiled or when a
project already standardizes on Cython, but neither condition applies here.

A small nanobind module maps more directly onto Baysor's existing C++17 and CMake
build. Its responsibility is limited to converting the stable request and result,
releasing the GIL, defining ownership, and translating exceptions. The
scientific operation and serializers remain ordinary C++ and are independently
testable without Python. This keeps the integration native without committing to
binding Baysor's internal object graph.

The choice is a maintainability decision, not a claim that Cython is incapable
or that every call must be in-memory. It should be revisited only if a concrete
platform, ABI, or interface limitation appears. The CLI remains supported as a
standalone frontend and an independent parity oracle, but it is not the Python
runtime backend.

Baysor's MIT licence and original copyright notice must accompany redistributed
source and binaries. The surrounding `baysor-python` code remains under its own
BSD-3-Clause licence.

### Why tiling is Python-orchestrated

The selected production architecture is Python-orchestrated tiling: Python plans
coarse overlapping tasks and calls native Baysor once per expanded tile. The
segmentation, molecular graph, BMM, and boundary kernels remain C++. Python does
not process individual molecules inside the iterative segmentation algorithm, so
scheduler overhead is amortized over substantial native jobs.

Moving the same independent-tile loop into C++ would not restore global Baysor
context or improve seam semantics. It would reproduce the same approximation
while moving manifests, scheduling, checkpointing, retries, and distributed I/O
into a less suitable layer.

A genuinely native tiled Baysor algorithm would be a different project: it would
partition the molecule graph and exchange assignments, component statistics, and
cell lifecycle decisions between tiles during BMM iterations. Such domain
decomposition could be scientifically superior, but it is a major Baysor
algorithm redesign and is outside the supported tiled architecture. It should be
reconsidered only if halo, rescue, seam, and grid-shift validation show that
independent overlapping runs cannot meet the acceptance gates.

`segment_tiled(...)` will use Dask Distributed as its supported execution engine.
When no client is supplied, `baysor-python` creates and owns a configured local
cluster. A caller such as Harpy may instead pass a compatible Dask client; it
does not implement a second Baysor executor. Tile descriptions, manifests,
transactional outputs, reconciliation, and scientific acceptance remain ordinary
`baysor-python` contracts rather than Dask collection semantics.

### Dask execution and worker-lifecycle contract

The tiled implementation should submit one coarse Future per expanded tile. A
task receives a small immutable tile descriptor, reads its staged Parquet input,
calls the nanobind operation, writes to an attempt-specific output directory, and
returns only structured paths, status, resource measurements, and provenance.
Large transcript clouds and native result objects must not pass through the Dask
scheduler.

Dask does not need a Baysor-specific file abstraction and does not inspect the
contents of these files. `baysor-python` teaches Dask the handoff through ordinary
task functions and two small, serializable contracts:

- `PreparedTile` identifies the tile, its durable prepared-input URI, core and
  halo bounds, input fingerprint, schema version, and preparation provenance;
- `TileArtifacts` identifies the tile and attempt, the validated native output
  URIs, status, resource measurements, resolved parameters, and native build
  provenance.

The corresponding dependency chain is:

```text
plan immutable core-plus-halo descriptors
                    |
                    v
stage_tile(source, descriptor)
  read/select with Dask; write durable prepared Parquet; release eager partition
                    |
                    | Future[PreparedTile]
                    v
segment_tile(prepared_tile, options)
  create unique attempt directory; call nanobind with paths
  validate outputs; atomically publish completion manifest
                    |
                    | Future[TileArtifacts]
                    v
coordinator validates completed manifest
                    |
                    v
reconcile committed artifact paths
```

Passing the `PreparedTile` Future to `segment_tile(...)` creates the scheduling
dependency: segmentation cannot start until staging has completed successfully.
The worker receives the small resolved descriptor and C++ opens the referenced
Parquet input itself. A pandas DataFrame, Dask partition, NumPy buffer, or complete
native result is not handed across the native segmentation boundary. Task
submissions that create attempt files must use non-pure/unique task semantics so
the scheduler does not incorrectly substitute a prior side effect for a new
attempt.

Writing files from a Dask task, including through an in-process C++ call, is an
intentional controlled side effect rather than an unsupported Dask pattern. It
is permitted here only because inputs are immutable, every execution has a
unique attempt identity and directory, retries never reuse partial files, and a
small atomic completion manifest is the sole commit point. The task is submitted
with non-pure semantics and returns only the committed `TileArtifacts`; no later
stage treats the mere presence of output files as success.

Every attempt writes under a unique location and publishes no completion marker
until the required outputs and schemas have been validated. A failed or cancelled
attempt remains uncommitted; an automatic or manual retry creates a new attempt
rather than trusting partial files. The small completed manifest, not directory
existence alone, is the durable success record and the value consumed by later
reconciliation tasks.

Paths are an explicit deployment contract. In the supported local cluster, all
worker processes share the configured run directory. A multi-host client must
provide a shared filesystem or supported object store visible to the relevant
workers, or a reviewed localization/upload layer that copies an input into worker
scratch and publishes outputs to durable shared storage before returning. Local
worker scratch is never the sole location of a prepared input or completed result
when the worker will be recycled. Accessibility, write permissions, atomic-publish
capability, and recovery after worker restart are checked during preflight.

The initial local production configuration is deliberately sequential:

```text
Dask workers             = 1
Dask threads per worker  = 1
simultaneous tiles       = 1
OpenMP threads per tile  = all allocated cores
worker lifetime          = one successful, failed, or cancelled tile attempt
```

The Dask thread executes the bound call and nanobind releases the GIL. OpenMP,
not the Dask thread pool, supplies parallelism within Baysor. If later benchmarks
use `W` simultaneous workers with `T` OpenMP threads each, the executor must
enforce `W * T <= C`, where `C` is the CPU allocation, and must also satisfy the
measured per-worker memory budget. A Python thread pool must not be used to run
multiple Baysor segmentations inside one worker process.

The bound C++ operation runs in the Dask worker process; it does not implicitly
create a child process. Dask continues to track the Future as processing and the
Nanny can monitor or replace the worker process, but the scheduler cannot inspect
the native call's internal phases or forcibly pre-empt its execution thread.
Releasing the GIL is therefore mandatory: it keeps the worker's Python control
thread responsive enough to service status and cooperative-cancellation
requests. GIL release does not itself make native code interruptible, so wall
timeouts and the hard-restart fallback remain part of the execution contract.

The process and supervision topology is:

```text
Python coordinator
        |
        | submit tile / inspect Future
        v
Dask scheduler ------------------------> Dask worker process
        |                                     |
        | targeted restart                    +-- Python task thread
        v                                           |
Dask Nanny                                          +-- nanobind
        |                                                 |
        +---- supervises/restarts worker                  +-- Baysor C++
                                                               |
                                                               +-- OpenMP threads
```

The scheduler communicates with the worker for normal task execution. The Nanny
is a separate supervisor, not part of that task-data path. Cooperative
cancellation reaches the registered native token through the still-responsive
worker process; hard cancellation asks the Nanny to replace that entire worker.

The one-worker configuration is the first validated operating mode, not an
architectural limit. Tile attempts already have immutable descriptors, isolated
attempt directories, and small structured results, so increasing `W` must be an
execution-policy change rather than a change to segmentation or reconciliation.
The implementation must keep worker count and native threads per worker as
independent explicit settings.

Every local worker must be Nanny-supervised. OpenMP and any other discovered
native thread-pool settings must be installed in the worker environment before
the native extension is imported and recorded in the run manifest. Dask resource
annotations reserve a Baysor execution slot, but they are not a substitute for
explicit native thread configuration or host-level CPU allocation.

Native C++ allocations are largely unmanaged memory from Dask's perspective.
Dask RSS limits and Nanny termination are safety guards; spilling cannot reduce
an active Baysor working set and allocator retention must not be mistaken for
reusable Dask-managed data. After a successful tile attempt has validated its
outputs and atomically published its completion manifest, the coordinator
validates the returned `TileArtifacts`; failed and cancelled attempts remain
uncommitted. The Nanny-supervised worker is proactively restarted after every
outcome and before another tile is submitted to it. This provides deterministic
reclamation of the native heap, OpenMP thread stacks, and library caches. The
worker identity before and after recycling is retained in provenance.

Cancellation uses an explicit worker control channel; constructing a token only
inside the task would not let the coordinator reach an executing native call.
Each worker owns a thread-safe active-attempt registry keyed by the unique
attempt identifier. Immediately before entering the binding, `segment_tile(...)`
creates and registers a native `CancellationSource`, passes its read-only token
to C++, and unregisters it with RAII/`finally` cleanup on every exit path. The
registry stores at most the worker's single active Baysor attempt and no native
result data.

The two-stage cancellation sequence is:

1. the coordinator locates the worker currently processing the attempt and sends
   a targeted request through a supported Dask worker control mechanism, such as
   a reviewed `Client.run(...)` or WorkerPlugin protocol;
2. the worker-side handler finds the attempt and sets the source's thread-safe
   cancellation flag, while C++ checks its token at safe phase and iteration
   boundaries;
3. the coordinator keeps the Future-to-worker association and waits for the
   distinct cancelled outcome for a configured grace period; and
4. if the native call does not return, the coordinator restarts that dedicated
   worker through its Nanny and only then releases or cancels the Future.

`Future.cancel()` alone is not treated as a hard interruption of executing C++
code and must not be issued first, because doing so can discard task-tracking
information needed for the targeted shutdown. If an attempt has already
finished and unregistered, the control handler reports it as inactive and the
coordinator rechecks the Future instead of treating that race as an error. The
control path must be integration-tested while the worker's sole task thread is
inside the GIL-released binding.

A native exception becomes a failed Future. A fatal native crash or OOM loses
the worker and is recovered by its Nanny. An unresponsive call reaches its wall
timeout and triggers the same dedicated-worker restart. In all cases, output
promotion is independent of task termination: incomplete attempt directories
have no completion manifest, remain uncommitted, and are never consumed by
reconciliation. A child subprocess per tile remains a possible stronger
isolation mode if production evidence requires it, but is not part of the first
supported design because the one-task worker is already disposable.

The locally owned cluster is dedicated to Baysor. A caller-managed client is
accepted only after preflight proves that its Baysor workers are exclusive,
Nanny-supervised or equivalently restartable, configured with one Dask execution
thread, equipped with the required native build and thread settings, and able to
access the staged inputs and output store. `baysor-python` must not restart a
shared worker containing unrelated tasks or irreplaceable data. Unsupported
clients fail preflight rather than silently weakening memory or cancellation
guarantees.

#### Controlled multi-worker extension

After the sequential lifecycle is correct, the local cluster may run multiple
worker processes. Every worker still has exactly one Dask execution thread and
may run at most one Baysor tile attempt. A supported configuration is therefore
described by:

```text
W = simultaneous Dask worker processes
T = OpenMP threads used by each Baysor call
C = compute cores available after coordinator and system headroom
M = memory available to the tiled run

W * T <= C
W * measured_peak_tile_memory + coordinator_memory + safety_margin <= M
```

In particular, `T = 1` permits several single-threaded tiles to run concurrently,
up to the CPU, memory, and storage-throughput limits. Each worker advertises one
Baysor execution slot, and the coordinator binds no more than one active tile to
that worker. Dask resource annotations express the scheduling constraint; worker
provisioning and, where required, OS or cluster-manager CPU affinity enforce the
actual CPU allocation.

The coordinator uses a bounded sliding window rather than submitting the entire
tile set. It initially assigns at most one tile to each eligible worker. As each
Future completes, the coordinator gathers and validates its small result, commits
it as the sole run-manifest writer, restarts that specific Nanny-supervised
worker, waits for its fresh replacement, and only then assigns the next tile to
that worker. Other workers continue their current tiles during this recycling.
Explicit worker assignment prevents a queued tile from starting in a worker
before its required restart.

`OMP_NUM_THREADS` and `OMP_THREAD_LIMIT` must be set consistently with `T` before
native import. The implementation must also inventory and constrain other
native thread pools used by the complete task, including input and output code;
setting the Dask thread count and OpenMP alone is not accepted as proof that the
host cannot be oversubscribed. Effective process and thread counts are measured
and retained in provenance.

The required fixed-core benchmark matrix includes at least `W=1, T=C`,
`W=2, T=floor(C/2)`, and a feasible single-threaded configuration with `T=1`.
Comparisons use end-to-end throughput, per-tile latency, aggregate peak RSS,
storage throughput, worker-restart overhead, failure rate, and scientific repeat
agreement. The sequential configuration remains the default until another point
in this matrix passes the same resource and scientific acceptance gates.

## Output contract

A completed run should add coordinated elements to SpatialData.

### Molecule assignments

The authoritative result is a points element containing the original retained
transcript columns plus at least:

- a stable transcript identifier;
- the final global cell identifier;
- `is_noise`;
- Baysor's molecule confidence, when available;
- assignment confidence, when available;
- the sampled global prior label; and
- optional globally meaningful molecule-cluster information.

There must be exactly one final row per retained input transcript. Original
transcript attributes and feature-panel metadata should be preserved.

### Cell shapes

A shapes element should contain one final polygon per retained global cell. Its
index must use the same instance identifiers as the table. Tile-local polygons
are intermediate data and must not be exposed as the final segmentation.

### Cell table

An AnnData table should contain:

- the global cell-by-gene count matrix;
- cell centroids;
- transcript counts;
- noise and assignment-confidence summaries where meaningful;
- dominant prior nucleus and nucleus-ownership QC fields;
- area and shape statistics; and
- SpatialData region and instance annotations linking it to the shapes element.

Counts must be rebuilt from the final molecule assignments. Tile-level feature
matrices cannot be concatenated because halo transcripts are duplicated and
tile-local cell identifiers are not globally meaningful.

### Optional labels raster

Raster labels may be generated from the final global shapes when explicitly
requested. Shapes and molecule assignments should remain the primary result;
the workflow must not require a full-mosaic raster merely to represent Baysor's
output.

### Provenance

SpatialData output metadata should record:

- the `baysor_python` version, embedded Baysor revision, native-extension path,
  and extension checksum;
- all resolved Baysor and adapter-side parameters;
- source points, prior labels, and coordinate system;
- global filtering decisions and feature-panel identity;
- the master native seed and, for tiled runs, the deterministic effective seed
  assigned to each stable tile identity;
- tile core and halo bounds;
- per-tile input and output checksums and status;
- Dask and Distributed versions, cluster ownership, worker configuration,
  effective OpenMP settings, task keys, and worker identities;
- per-tile wall time, peak process memory, cancellation state, and recycling
  outcome;
- reconciliation thresholds and ambiguity counts; and
- final coverage and seam-QC metrics.

## Common untiled foundation

Phase 1 establishes a low-level
`baysor_python.segment(...)` API and then places a thin SpatialData adapter on
top. The low-level Python API accepts prepared inputs and a structured request,
calls the nanobind extension, and returns paths, logs, resource measurements,
resolved native versions, and status in a structured result. The extension calls
the same coarse C++ operation used by the CLI. The adapter prepares and validates
Parquet, validates the output schemas, and constructs coordinated SpatialData
elements. A Harpy convenience wrapper may delegate to this adapter without being
a dependency of `baysor-python`.

The low-level interface is deliberately path-oriented and conceptually has this
shape:

```python
result = baysor_python.segment(
    molecules="molecules.parquet",
    output_dir="segmentation",
    prior_segmentation=":prior_cell",
    scale=53.0,
    scale_std="25%",
    n_cells_init=50_000,
    cluster_method="none",
    iters=100,
    tol=0.005,
    omp_threads=8,
)
```

It should return a structured, path-first result rather than eagerly loading a
large output into memory. The result includes output paths, native logs, return
status, wall time, resource measurements, resolved parameters, and the
package/Baysor build manifest. Expected output files and required columns must be
validated before a successful result is returned.

The corresponding high-level operation,
`segment_spatialdata(...)`, validates the selected points and coordinate system,
creates or preserves stable transcript identifiers, samples an optional labels
prior into a transcript column, writes the prepared Parquet input, calls
`segment(...)`, validates its outputs, and constructs coordinated SpatialData
points, shapes, and table elements. It must not contain another implementation of
the segmentation algorithm.

The operational flow is:

```text
SpatialData points + optional labels prior
                    |
                    v
prepare stable transcript-native Parquet input
                    |
                    v
low-level segment(...) -> nanobind -> shared C++ run operation
                    |
                    v
validate native Parquet/geometry/count outputs
                    |
                    v
construct coordinated SpatialData elements and provenance
```

The public API contract must expose:

- `scale` and `scale_std`;
- `n_cells_init`;
- `cluster_method` and related cluster parameters;
- `prior_segmentation_confidence`;
- `min_molecules_per_cell`;
- `iters` and `tol`;
- a 64-bit native `random_seed`, defaulting to `1`;
- native build identity and supported execution settings;
- OpenMP thread count;
- work directory and intermediate-retention policy; and
- overwrite and resume behavior.

The integration must always set an explicit initial cell count for large
nuclei-prior runs. Automatic initialization produced approximately 1.43 million
components in the benchmark and led to pathologically large polygon/output
work. For the UCB sample, the candidate range remains 40,000 to 60,000 initial
components.

The actual-sample reference study should evaluate:

- the Cellpose 4 nuclei sampled into a transcript column such as `prior_cell`;
- `cluster_method=none`;
- sample-specific `scale` estimated from the nuclei or cell labels;
- `n_cells_init` values of 40,000, 50,000, and 60,000;
- prior confidence 0.5, with 0.2 as a sensitivity check;
- a 100-iteration run before committing to 500 iterations; and
- `tol=0.005` unless convergence diagnostics justify another value.

## Tiled architecture

The following algorithm belongs to `baysor_python.segment_tiled(...)`. The
`baysor_python` SpatialData adapter prepares Dask-independent tile descriptors and
imports the reconciled result. `baysor-python` owns the Dask task implementation;
Harpy may supply an existing compatible client or a convenience entry point, but
it does not implement a second Harpy-specific execution, tiling, or reconciliation
algorithm.

### Statistical non-equivalence of tiled Baysor

Independent overlapping Baysor runs are not mathematically equivalent to one
untiled run. The implementation must therefore treat tiling as a validated
approximation, not merely as an execution optimization that is assumed to preserve
the result.

At pinned C++ revision `d7077a7ded6f4b941915badc894f767532d39fd2`, the
algorithm has both local and invocation-wide dependencies:

- gene filtering and encoding are computed over the invocation;
- `scale` and `scale_std` may be estimated from all prior-cell centres in the
  invocation;
- molecule confidence is obtained by fitting a two-component signal/noise model
  to KNN-distance statistics over the invocation;
- molecule-graph edge filtering and weights use invocation-wide edge-length
  quantiles;
- initial component placement and automatic component counts depend on the
  molecules in the invocation;
- optional molecule clustering is fitted over the invocation; and
- the BMM assignment update is predominantly local because a molecule considers
  only components represented among its molecule-graph neighbours.

The final boundary estimator is less globally statistical than its whole-cloud
signature suggests. For an ordinary cell with at least three assigned molecules,
it triangulates that cell's assigned molecules and uses non-cell molecules inside
the cell's bounding box to remove admixture. A dataset-wide mean nearest-neighbour
distance supplies the offset used for degenerate one- and two-molecule cells. A
global boundary rebuild can therefore reproduce the untiled boundary semantics
when final assignments and required local context are the same; the larger source
of tiled-versus-untiled difference is the assignment calculation that precedes
the polygon.

A naive implementation that estimates all parameters independently per tile,
uses no halo, or merges tile polygons is expected to be materially worse. With
large tiles, a nuclei prior, one global feature panel and scale, disabled molecule
clustering, sufficient halos, assignment reconciliation, rescue runs, and global
boundary rebuilding, the working expectation is more limited:

- dense-tissue cells well inside a tile core should often be close to the untiled
  result;
- differences should be enriched near seams, sparse or background regions,
  tissue transitions, unusually large cells, cells without prior nuclei, and
  ambiguous neighbours; and
- per-tile confidence/noise fits can produce differences throughout a tile, not
  only near a seam, particularly when tiles cover biologically different density
  regimes.

The proposed full expanded UCB tile contains approximately 4.42 million
transcripts at average density, so sampling variance in tile-wide estimates may
be small. Spatial non-identical-distribution and algorithmic stochasticity remain
the concerns. No claim that the difference is scientifically negligible may be
made until actual-data validation is complete.

### How tiling artefacts are avoided

The tiled workflow does **not** treat boundary polygons as the objects that must
be stitched. It treats molecule assignments as the authoritative intermediate
result. Cell identities are reconciled through the stable transcripts shared by
overlapping tiles, and one new global polygon is generated only after that
reconciliation is complete.

Consider one biological cell crossing the boundary between two tile cores. The
halo makes the cell and its surrounding molecule context visible in both Baysor
runs. Baysor may call it `cell_42` in tile A and `cell_17` in tile B, but many of
the same stable transcript IDs in the shared halo will be assigned to both local
cells. Sufficient agreement on those transcripts, supported by assignment
confidence, spatial proximity, and the prior nucleus, establishes that the two
local identifiers represent one global cell:

```text
tile_A/cell_42 ---\
                   +--- global_cell_9001
tile_B/cell_17 ---/
```

After all compatible local identifiers have been mapped to global cells,
`baysor_python.segment_tiled(...)` selects one final assignment for every
transcript. It then gathers the complete molecule cloud of `global_cell_9001`
from both sides of the seam and invokes the native boundary API to estimate a
single new boundary from that global cloud. The SpatialData adapter converts and
imports the resulting assignments and geometry.

The intended flow is therefore:

```text
overlapping tile context
          |
          v
match local cells using shared transcript identities
          |
          v
select one global cell assignment per transcript
          |
          v
estimate one new global polygon per cell
```

Tile-local polygons may be retained for diagnostics, but they are not clipped,
glued together, or unioned to form the final boundary. Direct polygon stitching
can preserve straight cuts at tile edges, introduce gaps or overlaps, and create
a shape that disagrees with the final molecule assignments. Re-estimating the
boundary after transcript reconciliation removes the tile seam from the geometry
construction altogether.

The halo is the first defence against artefacts because it prevents a core seam
from also being a Baysor context boundary. Shared-transcript matching is the
second defence because it resolves duplicate local identities. Global boundary
re-estimation is the third defence because it prevents tile edges from appearing
in the final cell geometry.

### 1. Global preflight and filtering

Input preparation must happen globally before tiling:

- require backed SpatialData for large runs;
- validate finite coordinates and required columns;
- establish or create a stable signed 64-bit transcript identifier;
- validate uniqueness of that identifier;
- resolve the target coordinate system;
- define identity and translation as the supported coordinate relationships
  between points and prior labels for this release scope;
- apply QV, feature-class, excluded-gene, and rare-gene filtering once;
- establish one global feature panel and gene encoding; and
- estimate one global `scale` that is reused by every tile.

Per-tile gene filtering is not allowed because it would give different gene
models to different tiles. Per-tile scale estimation is likewise not allowed.
After global filtering, tile runs must disable any second rare-gene filter that
could drop globally retained genes in low-count tiles.

Confidence/noise calibration is the largest unresolved global-statistics issue.
The pinned CLI currently recalculates molecule confidence for every invocation,
even when a confidence column is present. `baysor-python` should pursue an
upstream-compatible API that can either:

1. accept and preserve globally precomputed per-molecule confidence; or
2. accept one globally fitted signal/noise calibration and apply it consistently
   in every tile.

The global calibration may be fitted exactly or from a reproducible,
spatially-stratified sample, but sampled calibration must first be compared with
the untiled fit. If neither capability is available during engineering
validation, per-tile confidence fitting may be used only in a non-production
validation mode and must be recorded explicitly as an approximation. Its output
is not eligible for production promotion. Tile confidence distributions, fitted
signal/noise parameters, and noise fractions must be included in QC.

The molecule graph is locally constructed, but its edge filtering and weight
reference are based on per-invocation edge-length distributions. The tiled
implementation must record those values per tile and test whether they vary with
tissue region. If this variation produces material assignment differences, the
native API must support globally calibrated graph thresholds before production
promotion.

### 2. Sample the prior at transcript positions

The nuclei label must be sampled at each transcript coordinate to create a
transcript-native prior column. Transcripts outside a nucleus must remain in the
table with prior label `0`; they must not be discarded.

The existing blockwise point-to-label aggregation code can provide useful
building blocks, but the Baysor preparation path needs a variant that retains
background rows and preserves all required transcript attributes.

If global label identifiers are very large or sparse, each tile may compact its
prior labels and keep a sidecar mapping from tile-local prior labels to global
labels. This avoids allocations based on a large maximum label while preserving
global nucleus identity for reconciliation.

### 3. Plan tile cores and halos

Tile cores form a disjoint half-open partition of the data extent. Every
transcript therefore has exactly one core owner. Each core is expanded by a halo
before Baysor is run so cells and molecule neighborhoods near a core edge have
context from the adjacent tile.

For the UCB dimensions, the first configuration to validate is:

- core size: 10,000 by 10,000 coordinate units;
- halo: `4 * scale`, or 212 units when `scale=53`;
- six columns by three rows, for 18 tiles;
- approximately 4.42 million transcripts in a full expanded tile at the average
  dataset density; and
- approximately 1.081 times the original transcript volume after halo
  duplication.

The halo multiplier remains a validation parameter until the acceptance study
selects a supported default. Values of `2 * scale`, `4 * scale`, and `6 * scale`
must be compared.

Tile planning should use a configurable maximum transcript budget. A coarse
spatial count pass should estimate expanded-tile counts before writing inputs.
Tiles exceeding the budget should be recursively split. The planner must handle
irregular edge tiles and T-junctions without changing core ownership semantics.

### 4. Stage tile inputs in one pass

`baysor_python.segment_tiled(...)` should not run Baysor's crop flags repeatedly
against the full source Parquet file, because each run would still need to scan
the source. Instead, one staging pass, executed locally or through the owned Dask
workflow, should route transcripts to all expanded tiles that contain them.

The staging pass should:

- retain global coordinates;
- include only the columns Baysor and reconciliation need;
- duplicate halo transcripts while retaining their stable identifier;
- attach the sampled prior and any approved global filter/cluster columns;
- write one Parquet input per expanded tile; and
- write a manifest containing counts, bounds, checksums, local prior mappings,
  resolved parameters, and run status.

The manifest makes the workflow resumable. A completed tile may be reused only
when its input checksum, `baysor_python` version, embedded Baysor/native-extension
identity, Dask execution-policy version, and resolved parameters all match.

### 5. Run Baysor through Dask Distributed

`baysor_python.segment_tiled(...)` submits one native `segment(...)` call per tile
as a coarse Dask Future. By default it owns a local Nanny-supervised cluster; a
validated caller-managed client uses the identical tile task. The initial local
configuration submits one tile at a time to one worker with one Dask execution
thread, while the native call uses the explicitly configured OpenMP thread budget.

The orchestrator submits only as many tile attempts as the active worker and
memory budgets permit. Completed task results are gathered as small structured
metadata, validated, atomically committed in the run manifest, and followed by
proactive worker recycling. Automatic Dask rescheduling or retries must use a new
attempt directory and may reuse no partial output from a failed worker.

On the 32 GB benchmark machine, the required resource comparison is:

- one Baysor process using eight OpenMP threads; versus
- two Baysor processes using four OpenMP threads each.

Each tile must use:

- the same fixed global scale and scale standard deviation;
- the globally selected feature panel, with tile-local rare-gene filtering
  disabled;
- globally consistent QV and feature-class filtering;
- the same global confidence/noise calibration for production runs; an explicitly
  recorded per-tile fit is permitted only in non-production validation mode;
- `cluster_method=none` for the production baseline unless validation supports a
  different globally consistent method;
- `--skip-ncv-color`;
- the same iteration and convergence parameters; and
- an explicit tile-level initial component count.

A practical tile-level initial count is approximately two to three times the
number of active prior nuclei in the expanded tile, with a global cell-density
fallback for tissue without prior nuclei. Automatic initialization must not be
used.

Tile logs, resolved parameters, wall time, peak process memory, Dask task and
worker identities, native status or exception, and worker-recycling outcome must
be retained. A failed tile should not invalidate successful tiles, but the
workflow must not reconcile a partial tile set.

### 6. Reconcile local cells through shared transcripts in `baysor_python`

The same stable transcripts occur in both runs on either side of an internal
seam. They provide direct evidence for whether two tile-local cell identifiers
represent the same biological cell.

For each pair of spatially adjacent tiles:

1. Join the two molecule outputs by stable transcript identifier within their
   shared halo.
2. For every pair of non-noise local cell identifiers, count transcripts assigned
   to both cells.
3. Compute evidence including shared-transcript count, overlap coefficient,
   Jaccard similarity, mutual-best status, centroid distance, assignment
   confidence, and agreement on the dominant global prior nucleus.
4. Generate candidate match edges only when minimum evidence criteria are met.
5. Process candidate edges from strongest to weakest with constrained
   union-find.

The union operation must maintain the invariant that a global component contains
at most one local cell from any individual tile. This prevents a chain of
pairwise matches from merging two cells that the same Baysor run considered
distinct.

Many-to-one matches, conflicting nucleus evidence, low-overlap matches, and
components that would violate the one-cell-per-tile invariant are ambiguous.
They should be recorded for QC and resolved conservatively rather than merged
silently.

For example, if tile A splits a region into `cell_A1` and `cell_A2` while tile B
calls the same region `cell_B1`, the reconciler must not merge all three cells:
doing so would merge two cells that one Baysor run explicitly kept separate. It
should accept at most the strongest compatible match and flag the remaining
conflict. An ambiguous or halo-touching group can be rerun as a rescue region
centred on the seam with a larger halo, or on a grid shifted so that the seam is
inside a tile core. If the conflict remains, retaining a conservative separation
with an explicit QC flag is preferable to an unsupported merge.

Thresholds should be learned from tiled-versus-untiled comparisons. They must
not be chosen solely from the synthetic resource benchmark.

### 7. Select one final assignment per transcript

After local cell identifiers have been mapped to global components, a transcript
may still have several tile predictions. The preferred prediction should come
from the tile in which that transcript is farthest from the outer boundary of
the expanded tile, because this prediction had the most spatial context.
Assignment confidence can be used as a secondary criterion. The unique core
owner provides the deterministic fallback.

Noise is a valid competing assignment and should not be overwritten merely
because another tile assigned the transcript to a low-confidence cell.

The result of this stage is a single partitioned molecule table with exactly one
row per retained transcript and one globally unique cell identifier or noise.

### 8. Build global cell products

All final products must be rebuilt from the reconciled molecule table:

- cell-by-gene counts;
- centroids and transcript counts;
- confidence and nucleus-ownership summaries;
- cell-level QC fields; and
- cell boundaries.

`baysor_python` returns the authoritative reconciled molecule assignments and
native boundary result. Its SpatialData adapter derives the tables, QC fields,
shapes, and optional labels from that result rather than from tile-local
products. Harpy may add presentation or application-specific metadata without
changing the authoritative assignments.

Tile count matrices must be ignored. Tile polygons must not simply be clipped
and unioned because doing so can create visible seams and can disagree with the
final transcript assignments.

The selected boundary path is to estimate each final cell polygon once from the
globally reconciled molecule assignments. Thus, a cell crossing a seam receives
one boundary calculation over the molecules on both sides, rather than two
polygon fragments followed by a geometric merge.

"Global molecule cloud" does not mean that the estimator may see only the
molecules assigned to the target cell. Baysor's boundary algorithm also uses
nearby molecules assigned to other cells or to noise to reject admixture around
the Delaunay boundary. An exact whole-dataset call therefore receives the complete
reconciled coordinates and cell labels. A bounded-memory batched call may receive
only a subset of target cells, but it must also receive every contextual molecule
intersecting those cells' required bounding regions. The global boundary-distance
parameter must be computed once or passed unchanged to every batch.

### `baysor_python` boundary-estimator decision

Neither Harpy nor `baysor_python` will reimplement the estimator in Python, and
neither will merge tile-local polygons. The separate `baysor-python` repository
will build the pinned Baysor source and expose its boundary estimator through a
thin `pybind11` or `nanobind` extension. Harpy will depend on this public Python
API rather than owning Baysor-derived C++ code.

The native extension boundary must remain array-oriented and independent of
SpatialData, GeoPandas, Shapely, Arrow, and Parquet. Its conceptual contract is:

- input: contiguous molecule coordinates, final global cell labels, optional
  target cell identifiers, and an optional precomputed global boundary-distance
  parameter;
- output: packed polygon vertices, polygon offsets, and the corresponding global
  cell identifiers; and
- execution: release the Python GIL, support sparse global cell identifiers, and
  handle empty, one-molecule, two-molecule, and collinear cells explicitly.

The `baysor_python` SpatialData adapter will prepare the arrays, call
`baysor_python.boundaries(...)`, convert the packed result to Shapely polygons
and a GeoDataFrame, and construct the shapes element. A Harpy wrapper may expose
the same result without duplicating this conversion.

`baysor_python` should be built with CMake and `scikit-build-core`. The required
release matrix covers macOS ARM64 and Linux x86-64 for the package's supported
Python versions. Its releases must include Baysor's MIT notice, the exact
upstream commit and original source paths, an inventory of binding-related Baysor
changes, and an update procedure. The Python package must report both its own
version and the embedded Baysor revision.

Parity tests must compare the native API geometry with the pinned Baysor CLI
output after normalizing irrelevant polygon orientation and starting-vertex
differences. Repeated-call tests must also check memory release, exception
translation, and OpenMP behavior.

## Baysor upstream requirements and optimizations

### Stable transcript identity in Parquet output

This is the most important correctness requirement. At revision
`d7077a7ded6f4b941915badc894f767532d39fd2`, Baysor reads numeric
`transcript_id`, and legacy CSV can write it, but `molecules.parquet` does not
include it. Tiled reconciliation should not depend permanently on undocumented
row-order preservation.

The production solution is an upstream Baysor change that round-trips
`transcript_id` in Parquet output. Until that exists, a pinned compatibility
adapter may attach input IDs to output rows only after strict validation of row
count and exact gene/x/y equality. The adapter must fail closed when validation
does not hold, remain explicitly documented as a compatibility constraint, and
must not become the permanent identity contract.

### Selective outputs

Tile-level count matrices, cell statistics, and polygons are not authoritative,
yet Parquet output currently generates them. Upstream switches to skip these
products would reduce tile runtime, memory, and disk use. A molecules-only mode
is desirable.

### Stable native library entry points

Phase 1 must add one stable, coarse-grained C++ segmentation entry point before
adding Python integration. Conceptually it accepts a `SegmentationRequest` and
returns a `SegmentationResult`. The request carries the molecule input and
resolved run options; the result owns the authoritative molecule assignments,
confidence values, cell statistics, boundaries, counts, resolved options, and
run provenance needed by the shared serializers.

The operation must contain the scientific orchestration that currently resides
in the CLI's `cmd_run(...)`. CLI parsing, presentation, and process exit codes
remain in the CLI frontend. File serialization should be a shared library
facility over `SegmentationResult`, not interleaved differently in the CLI and
binding. Both frontends must call the same operation so parameter resolution and
scientific behaviour cannot drift.

The Python binding should expose only the coarse operation and stable result
views needed by `baysor_python`; it must not expose internal BMM implementation
types as public Python API. The binding releases the GIL around the native run,
defines ownership for every returned view, and translates C++ exceptions at the
module boundary.

The operation must also accept a thread-safe cooperative cancellation token. A
separate source retains the ability to request cancellation from another thread;
the token passed to the run operation provides read-only access to their shared
state. This source/token contract remains ordinary C++ and has no Dask or Python
dependency. It checks that token only at defined safe phase and iteration
boundaries and returns a distinct cancelled outcome without publishing a
complete result. The CLI and nanobind frontends use the same cancellation path;
process termination remains a supervisor fallback when native cancellation does
not complete within its grace period.

The Phase 2 boundary API will bind the existing boundary functions, but
bounded-memory use requires explicit target-cell selection while retaining all
required contextual molecules. This capability should preferably be implemented
in Baysor's public C++ library rather than as a divergent algorithm in the
binding layer.

### Global confidence and noise calibration

At the pinned revision, each CLI segmentation invocation recomputes molecule
confidence from its own KNN-distance distribution and signal/noise mixture fit;
an input confidence column is not preserved. Consequently, overlapping tile
runs can acquire tile-wide statistical differences even far from a seam. This
is the most important remaining difference after gene filtering, scale, and
clustering have been made globally consistent.

Before tiled mode is production-supported, Baysor must expose either:

- an option to accept and preserve globally precomputed per-molecule confidence;
  or
- a fit/export/apply interface for the signal/noise calibration, so one global
  calibration can be reused by every tile.

The fit result and provenance should include the KNN definition, mixture-model
parameters, convergence information, and summary signal/noise proportions. A
reproducible spatially stratified sample may be used if an exact global fit is
not practical, but its result must first be compared with the untiled fit.

Engineering validation may use tile-local confidence fitting to exercise the
workflow, but those results are non-production and must retain each tile's fitted
parameters and confidence distribution and test for spatial discontinuities.
Global or consistently applied calibration is a scientific-comparability
requirement before tiled mode is production-supported, not merely a runtime
optimization.

## Implementation sequence

### Phase 1: Establish `baysor_python` segmentation and the untiled integration

Phase 1 is divided into independently reviewable implementation slices. Each
slice must complete its focused tests and contract before the next layer is
added. The full actual-data reference is deliberately deferred; it is not an
entry gate for these slices. Native engine changes are implemented and reviewed
in `vibspatial/Baysor`; Python binding and product-integration changes are
implemented in `baysor-python`. The submodule commit is the only source handoff
between them.

#### Slice 1A: Extract and hand off the coarse-grained public C++ operation

This slice contains no Python API, nanobind, or SpatialData implementation. It is
split into native-engine work in `vibspatial/Baysor` and a source-integration
handoff in `baysor-python`.

##### Slice 1A.1: Native API extraction in `vibspatial/Baysor`

Deliverables:

- preserve the immutable `cpp` baseline at upstream commit
  `d7077a7ded6f4b941915badc894f767532d39fd2` and perform the extraction on a
  dedicated feature branch;
- establish a small native regression fixture before changing orchestration,
  recording the exact committed pre-extraction revision of the current feature
  branch as its behavioural baseline; the upstream commit remains ancestry and
  source provenance rather than a second regression oracle;
- keep the entry-point extraction as a small, auditable commit series that can
  be proposed upstream;
- define `SegmentationRequest`, `SegmentationResult`, and
  `run_segmentation(...)` in the Baysor C++ library;
- include a 64-bit `random_seed` with default `1`, run-local random state, and
  resolved-seed provenance in that native contract;
- define a C++17-compatible thread-safe cancellation token and safe cancellation
  checkpoints for the complete run operation;
- move scientific orchestration and parameter resolution out of CLI
  `cmd_run(...)` into that operation without changing algorithm behaviour;
- keep file writers and reporting as reusable library facilities over the result;
- remove CLI-only process concerns, global logger replacement, and process exit
  codes from the scientific operation; and
- add focused C++ tests that call, cancel, and repeat the operation directly on a
  small fixture.

The fork owns every modification required to make `baysor_lib` reusable,
including any generic CMake option or position-independent-code support needed
to embed the library. It must not add nanobind or Python-specific types. Those
belong to `baysor-python` after the native contract is stable.

At the pinned revision,
[`cmd_run(...)`](https://github.com/kharchenkolab/Baysor/blob/d7077a7ded6f4b941915badc894f767532d39fd2/src/cli/main.cpp#L44)
combines typed-option validation, data loading, prior handling, data-dependent
parameter inference, confidence fitting, graph construction, molecule clustering,
2D/3D BMM dispatch, cell statistics, boundaries, count matrices, serialization,
HTML reporting, output-directory creation, process-wide logger replacement, and
integer exit codes. Its local lambda named `run_segmentation` is only a
dimension-dispatch closure over that CLI state; it is not the public operation
defined by this slice.

The target C++ boundary is:

```text
SegmentationRequest + CancellationToken
                    |
                    v
        C++ run_segmentation(...)
                    |
          +---------+---------+
          |                   |
          v                   v
SegmentationResult     distinct cancelled outcome
```

`SegmentationRequest` is a typed native contract. It contains the molecule-input
specification, optional prior-segmentation specification, molecule/filtering
options, segmentation and clustering options, confidence/noise options, requested
scientific products, a 64-bit native random seed, and execution settings. It
contains no `argv`, CLI aliases, presentation strings, process exit policy, or
Python objects. The frontend parses its configuration into this request. The
operation performs the canonical data-dependent resolution that requires loaded
molecules or a prior, such as inferred scale and initial cell count, and records
the final resolved values and seed in the result. No frontend may independently
reproduce that resolution or substitute another random stream.

`SegmentationResult` owns the complete scientific products required by the
serializers: retained molecule identity and data, cell/noise assignments,
molecule and assignment confidence, optional molecule clusters, cell statistics,
boundaries, counts, convergence and diagnostic data, final resolved options, and
run provenance. Ownership must use normal C++ lifetime management and moves where
appropriate; it must not expose dangling views into a destroyed `BmmData` or
require callers to manage Baysor internals.

The extracted operation owns this complete scientific sequence:

```text
validate typed request
        |
load and filter molecules; load optional prior
        |
resolve data-dependent parameters and confidence/noise model
        |
build molecule graph and perform optional molecule clustering
        |
initialize and run the existing 2D or 3D BMM implementation
        |
derive assignments, statistics, boundaries, counts, and diagnostics
        |
return an owned SegmentationResult
```

This is an orchestration and ownership refactor, not an algorithm port or rewrite.
The existing loaders, confidence estimator, graph and clustering functions, BMM,
boundary estimator, and scientific calculations remain authoritative. The public
API is one complete source-level C++ operation in `baysor_lib`; it is not a frozen
binary SDK and does not expose BMM components or the internal object graph.

Serialization remains a reusable library layer over the completed result:

```text
SegmentationResult
        |
        +-- molecule-table serializer
        +-- cell-statistics serializer
        +-- boundary serializer
        +-- count-matrix serializer
        +-- parameter/provenance serializer
        +-- diagnostic/report generator
```

Raw CLI parsing, output-path and filename selection, console presentation,
process-wide logger configuration, and conversion of failures to executable exit
codes stay outside `run_segmentation(...)`. The operation may report progress
through an explicitly supplied interface, but it must not replace application
global state, call `std::exit`, create directories through shell commands, or
reduce structured failures to integer return codes.

Cancellation is cooperative and has one shared native implementation. The token
is safe to set from another thread and is checked after major phases and at safe
BMM iteration boundaries. Adding those BMM checkpoints must not alter the update
sequence when cancellation is not requested. Cancellation returns a distinct
outcome without a complete `SegmentationResult`; completed serializers must not
run for that attempt. Process termination remains a later supervisor fallback,
not the primary native cancellation mechanism.

Randomness likewise has one shared native implementation. The recorded N0
pre-extraction CLI is captured in a fresh process with one OpenMP thread; its
implicit segmentation seed is `1`. The maintained operation creates
run-local random state from `SegmentationRequest::random_seed` rather than
resetting mutable process-global state. Scientific and diagnostic-only random
streams are separated so requesting a report cannot change segmentation. The
same input, build, platform, seed, options, and one-thread configuration must be
repeatable across fresh and repeated same-process calls. A seed alone does not
promise bitwise-identical multi-threaded execution because dynamic OpenMP
scheduling and parallel floating-point work can still vary; those runs retain a
locked thread configuration and a measured repeatability baseline.

Focused C++ tests must cover a complete direct run, inspection of core result and
resolved-option fields, invalid-request error handling, cancellation before a
completed result is published, result ownership after working objects are
destroyed, and repeated same-process calls without hidden global-state leakage.
The same-process coverage includes repeated one-thread calls with the same seed
and verifies that one call does not advance another call's random stream.
Full pinned-CLI parity is the Slice 1B gate.

Native exit criterion: a C++ test can execute a complete segmentation through
`run_segmentation(...)`, inspect its owned structured scientific result, request
cooperative cancellation, and repeat the operation without invoking CLI or Python
code.

##### Slice 1A.2: Pin the native handoff in `baysor-python`

Deliverables:

- add `vendor/Baysor` as a submodule of
  `https://github.com/vibspatial/Baysor.git` at the reviewed Slice 1A.1
  integration commit;
- retain Baysor's MIT licence and record both the immutable upstream baseline
  and selected fork integration commit in build provenance;
- add a CI guard that rejects an absent, dirty, or unexpectedly advanced
  submodule checkout; and
- add a focused configuration/build smoke check proving that the pinned native
  library target can be consumed without modifying its checkout.

This handoff does not yet build the nanobind module. Its purpose is to establish
one reproducible native source dependency and prove that `baysor-python` consumes
the public library target cleanly.

Slice 1A exit criterion: the native operation passes its focused C++ tests in
`vibspatial/Baysor`, and `baysor-python` pins and verifies the exact reviewed
commit without local submodule modifications.

#### Slice 1B: Make the CLI a frontend over the shared operation

This slice is implemented in `vibspatial/Baysor`.

Deliverables:

- reduce the CLI to argument/configuration parsing, invocation of
  `run_segmentation(...)`, serialization, and user-facing reporting;
- expose `--seed` as the CLI spelling of the shared native seed, with default
  `1`, and record it in resolved options and provenance;
- ensure the CLI and direct C++ path use the same parameter resolution and
  serializers;
- compare the refactored CLI with the untouched pinned CLI on the same fixture;
  and
- after the parity gate passes, advance `baysor-python/vendor/Baysor` to the
  reviewed Slice 1B integration commit and update its recorded provenance.

Exit criterion: the refactored CLI is semantically equivalent to the pinned
reference and contains no second segmentation orchestration path.

#### Slice 1C: Add the thin native Python binding

This slice is implemented in `baysor-python`. Generic build-system changes needed
by an embeddable C++ consumer remain owned by `vibspatial/Baysor`; the nanobind
module, wheel build, and Python contracts do not move into the fork.

Deliverables:

- replace the pure-Python build backend with `scikit-build-core` and build the
  pinned Baysor library as position-independent code plus a small nanobind
  module;
- add CMake install rules that place the linked native extension in the Python
  package without making the standalone CLI a Python runtime dependency;
- configure the binding for CPython's Stable ABI with a Python 3.12 floor and
  prove that every required nanobind feature compiles against the Limited API;
- expose the coarse run operation as `baysor_python._baysor_core`;
- retain a path-oriented request and result boundary for large backed datasets;
- release the GIL for the native calculation;
- bind the native cancellation source/token contract and maintain a thread-safe
  active-attempt registry so a separate worker control thread can request
  cooperative shutdown while the run has released the GIL;
- translate C++ exceptions into documented Python exceptions;
- define ownership and lifetime for every returned object or view; and
- emit the package/Baysor build manifest.

This slice does not bind Baysor's internal object graph or design a general
zero-copy API. Those capabilities require measured need and a separate contract.
If focused compatibility or performance tests reject the Stable ABI, this slice
must record the concrete failure and switch to conventional version-specific
CPython extension builds. It must not make that switch merely because they are
the default nanobind build mode.

Exit criterion: Python calls the shared C++ operation directly and receives a
validated structured result without launching the CLI.

#### Slice 1D: Establish the public `baysor_python.segment(...)` API

Deliverables:

- typed request, options, result, and exception contracts;
- validation of paths, columns, supported parameters, output location, and
  OpenMP settings;
- stable version and native-build reporting;
- native log, resource, and resolved-parameter capture; and
- focused Python tests for validation, exception translation, cancellation,
  output-schema failures, and repeated-call memory behaviour.

Repeated direct calls must show no linear growth in live native allocations or
RSS beyond an explicitly measured allocator plateau. The API does not promise
that an in-process allocator returns all freed pages to the operating system;
full-run supervisors may select a disposable worker when hard reclamation is a
requirement.

The public API calls the bound C++ operation; it does not offer a separate CLI
backend or executable override.

Exit criterion: `baysor_python.segment(...)` reproduces the refactored CLI on the
small parity fixture and exposes a documented, stable Python contract.

#### Slice 1E: Add the untiled SpatialData adapter

Deliverables:

- global prior sampling and stable transcript-ID preparation;
- backed Parquet preparation without full pandas materialization;
- invocation of `baysor_python.segment(...)`;
- output-schema and transcript-identity validation;
- construction of coordinated SpatialData points, shapes, table, and optional
  labels elements; and
- complete provenance and a focused SpatialData write/read integration test.

Exit criterion: the SpatialData path reproduces the same native assignments and
coordinated products without depending on Harpy.

#### Slice 1F: Complete supported release packaging

Deliverables:

- assemble a source distribution that materializes the exact Baysor submodule
  tree, pinned `FetchContent` sources, build files, build manifest, and licence
  notices as ordinary archive files;
- prove that the source distribution configures without Git, submodule
  initialization, or build-time network access when its documented system native
  prerequisites are present;
- build the preferred Python-3.12-floor `abi3` wheels from the produced source
  distribution with `cibuildwheel` for macOS ARM64 and Linux x86-64, or build the
  documented conventional per-minor CPython fallback if Slice 1C rejected the
  Stable ABI;
- verify that the preferred artifacts carry the expected `cp312-abi3` wheel tag
  and run `abi3audit --strict` after platform repair, rejecting any non-stable
  CPython symbol or ABI-floor mismatch;
- complete and document the Arrow/Parquet packaging decision after testing both
  private bundling and compatible PyArrow linkage in a SpatialData process;
- repair and inspect the wheels with `delocate` or `auditwheel`, including their
  native dependency inventory and run paths;
- verify that wheels contain the importable native extension and required
  redistributable libraries, but do not depend on a system Baysor executable or
  submodule checkout;
- install the same `abi3` wheel under CPython 3.12, 3.13, and 3.14 on each
  platform, or install the matching per-minor fallback wheel, and run import,
  small segmentation, cancellation, exception, array-conversion, and
  SpatialData/PyArrow coexistence tests; and
- documented upstream-update, ABI, provenance, and release procedures.

Exit criterion: the source distribution is complete and network-independent as
defined above; every platform in the release matrix installs the finished wheel
without a compiler or system Baysor installation, imports the native module,
reports the expected Baysor revision, has no unresolved or unintended native
dependencies, coexists with SpatialData/PyArrow, and passes the small native
parity test. An `abi3` release must prove that the same wheel passes this gate on
every supported CPython minor version and has a clean strict Stable-ABI audit.

#### Phase 1 parity contract

All parity paths must use the exact same input files, pinned Baysor revision,
configuration, OpenMP thread count, output mode, and controllable execution
settings. The recorded N0 pre-extraction CLI does not expose a seed flag and is
therefore run as a fresh one-thread process using its implicit seed of `1`. The
shared native operation, refactored CLI, and Python API expose that same seed
explicitly; baseline parity uses `random_seed=1`. Multi-threaded comparisons lock
the thread count and are repeated to characterize scheduling and floating-point
variability rather than assuming that a seed guarantees bitwise identity.

Comparison includes retained transcript identity, cell/noise assignments,
molecule confidence, count matrices, cell statistics, resolved parameters, and
the Baysor revision. Cell identifiers and polygon output are compared
semantically after normalizing harmless relabelling, polygon orientation, and
starting-vertex differences.

Overall Phase 1 exit criterion: the CLI and native Python API share one C++
operation, the Python result is semantically equivalent to the pinned CLI result,
the SpatialData adapter preserves that result, and supported wheels pass their
release tests.

### Phase 2: Establish the `baysor_python` native boundary API

Extend the established native binding with an independently callable boundary
operation and resolve its semantic-parity questions before building the tiled
workflow around it.

Deliverables:

- an array-oriented binding to the pinned Baysor C++ boundary implementation;
- `baysor_python.boundaries(...)` and a thin SpatialData conversion adapter;
- explicit full-context and batched-target input contracts;
- packed polygon output and conversion to SpatialData shapes;
- focused edge-case, repeated-call, and CLI-parity tests;
- macOS ARM64 and Linux x86-64 build and wheel smoke tests; and
- dependency, licence, provenance, source-bundling, and upstream-update
  documentation.

Exit criterion: for identical molecule coordinates, assignments, and parameters,
the Python API produces geometry equivalent to pinned Baysor and its wheel can be
installed on every platform in the required release matrix.

### Phase 3: Implement `baysor_python.segment_tiled(...)` planning and staging

Deliverables:

- Dask-independent tile descriptor, attempt-output, and result contracts;
- a master run seed and a deterministic effective native seed derived from the
  master seed and stable tile identity, both retained in each tile descriptor;
- SpatialData-adapter stable transcript-ID and prior preparation against that
  contract;
- one globally filtered feature panel and a contract that disables a second
  rare-gene filter inside tile runs;
- a global confidence/noise-calibration input contract, with an explicit
  unsupported-capability result from backends that do not implement it;
- half-open core and halo planner;
- density/count preflight and adaptive tile splitting;
- single-pass routing to tile Parquet files;
- local-to-global prior-label mappings;
- checksummed run manifest; and
- resume and stale-output detection.

Exit criterion: every input transcript has exactly one core owner and the
expected halo memberships, with no uncovered or unexpectedly duplicated rows.

### Phase 4: Implement resource-aware Dask tile execution

Deliverables:

- an optional `tiled` dependency extra containing a supported Dask Distributed
  version range;
- a locally owned `LocalCluster` using Nanny-supervised worker processes;
- coarse Future submission over the Phase 3 tile descriptors and small structured
  task results rather than molecule-cloud transfer through the scheduler;
- a caller-supplied Dask-client contract with worker, native-build, storage, and
  restart-capability preflight;
- the initial one-worker, one-Dask-thread, one-tile-at-a-time policy;
- a separately reviewed bounded multi-worker extension with independent
  `worker_count` and `openmp_threads_per_worker` settings;
- explicit worker assignment, one Baysor execution slot per worker, and a sliding
  submission window that never queues a second tile ahead of worker recycling;
- explicit pre-import OpenMP configuration, GIL release, native-thread reporting,
  an audit of other native thread pools, and enforcement of the total CPU budget;
- attempt-specific output directories, validation, atomic manifest commit, and
  proactive worker recycling after every tile attempt;
- a per-worker active-attempt cancellation registry, a targeted coordinator-to-
  worker control request, cooperative native cancellation before Future release,
  and Nanny worker restart after a grace period, with no promotion of partial
  outputs;
- explicit wall timeouts and documented handling for Python/native exceptions,
  worker loss, fatal native failure, OOM, and an unresponsive native call;
- RSS and unmanaged-memory monitoring with memory-budget admission control;
- explicit per-tile `n_cells_init` calculation;
- retained logs, resource traces, confidence/noise-fit summaries, and
  molecule-graph threshold summaries, including Dask task and worker identities;
- application of the common global confidence/noise calibration once supported
  by the native backend;
- workflow-controlled retry/resume behavior that never reuses an incomplete
  attempt directory;
- reuse of the same effective tile seed for every retry or resumed execution of
  a given compatible tile descriptor; and
- validation that every expected tile completed with compatible schemas and
  parameters.

Exit criterion: a complete set of reproducible tile-local molecule assignments
can be generated without exceeding the configured memory or CPU budget; each
successful, failed, or cancelled attempt leaves a valid manifest state and the
next tile runs in a fresh worker process. The multi-worker slice additionally
proves that each worker runs at most one tile, recycling one worker does not
interrupt the others, and observed native concurrency satisfies the declared
`W * T <= C` budget.

### Phase 5: Implement reconciliation in `baysor_python`

Deliverables:

- overlap joins by transcript ID;
- candidate cell-match metrics;
- constrained union-find;
- ambiguity reporting;
- final per-transcript assignment selection; and
- deterministic global cell relabeling.

Exit criterion: one final row exists per retained transcript, reconciliation is
order-independent, and all graph invariants pass.

### Phase 6: Build and import global products

Deliverables:

- reconciled assignments and global boundaries returned by
  `baysor_python.segment_tiled(...)`, with complete contextual molecules for each
  exact or batched boundary call;
- SpatialData-adapter global sparse cell-by-gene aggregation;
- cell statistics and nucleus-ownership QC;
- SpatialData points, shapes, and table elements constructed by
  `baysor_python`;
- optional rasterization; and
- complete provenance metadata.

Exit criterion: assignments, shapes, table instances, and optional labels are
mutually consistent and survive a SpatialData write/read round trip.

### Deferred validation milestone: Establish the actual-data reference

This milestone was formerly Phase 0. It is intentionally deferred until Phase 1
has produced a working native API and small CLI-parity test. It must be completed
before Phase 7 and before tiled mode can be assigned production defaults.

Run the actual UCB mosaic untiled with the Cellpose nuclei prior and explicit
initial cell counts. Compare `cluster_method=none` with Louvain on representative
crops before deciding whether clustering materially improves segmentation.

Deliverables:

- selected baseline parameters, including the native random seed;
- repeated untiled runs under the locked seed, build, and execution settings to
  quantify residual stochastic variability;
- an untiled 100-iteration full-mosaic reference result;
- resource measurements;
- visual overlays and biological QC; and
- frozen reference molecule-assignment datasets and an untiled-versus-untiled
  disagreement baseline for tiled comparisons.

Exit criterion: a scientifically plausible untiled result and a parameter set
worth reproducing in tiled mode before tiled-quality validation begins.

### Phase 7: Validate tiled quality and choose defaults

Deliverables:

- tiled-versus-untiled comparison on the actual UCB mosaic;
- untiled-versus-untiled baseline comparison under the same locked settings;
- halo-size comparison;
- tile-size and fixed-core concurrency benchmarks covering at least one worker
  using all compute cores, two workers dividing the cores, and a feasible
  one-native-thread-per-worker configuration;
- a half-tile grid-shift experiment;
- seam-specific QC plots; and
- documented validated defaults and failure thresholds.

Exit criterion: tiled results meet the agreed quality gates and do not show
material seam or grid-placement dependence.

### Phase 8: Remote and multi-node Dask hardening

Only after the local tiled workflow is correct:

- harden `baysor-python` support for caller-managed Dask clients and dedicated
  Baysor worker resources;
- document and validate Nanny or deployment-equivalent worker restart semantics,
  CPU allocation, native environment propagation, and supported cluster versions;
- support shared worker-local staging and remote-backed inputs where practical;
- improve adaptive splitting for heterogeneous density;
- define intermediate cleanup and retention policies; and
- test worker loss, partial retries, and manifest recovery.

Exit criterion: distributed execution changes throughput and capacity without
changing the scientific result or reconciliation semantics.

## Test strategy

### Focused unit tests

- tile cores cover the domain exactly once;
- half-open bounds behave correctly on seams and dataset edges;
- halo membership is correct at faces, corners, and T-junctions;
- stable IDs remain unique after routing;
- tile initial-cell estimates are explicit and bounded;
- candidate scores are invariant to input ordering;
- constrained union-find rejects two cells from the same tile in one component;
- ambiguous split, merge, and noise cases remain flagged;
- final transcript selection is deterministic;
- tile plans and reconciled results are invariant between the locally owned Dask
  cluster and a compatible caller-managed client;
- caller-managed clients without exclusive restartable Baysor workers fail
  capability preflight; and
- manifests reject parameter, input, or native-build mismatches.

### Small integration tests

Use a tiny synthetic dataset with known cells crossing tile seams. Run both
untiled and tiled Baysor and verify coverage, matching, counts, shapes, and
round-trip SpatialData integrity. The test should include duplicate coordinates,
background transcripts, cells without a prior nucleus, and nuclei close to a
tile boundary.

At the `baysor_python` layer, verify that `segment(...)` matches the direct pinned
CLI for identical inputs, parameters, output mode, and controllable execution
settings. The recorded pre-extraction CLI reference runs in a fresh one-thread
process with its implicit seed `1`; the refactored CLI and Python path use
explicit seed `1`.
Native tests additionally verify same-seed repeated same-process calls and that
diagnostic selection does not perturb segmentation. Multi-threaded comparisons
are repeated when characterizing parallel variability. Verify separately that
`boundaries(...)` matches the CLI-generated geometry for identical assignments.
Native tests must also cover exception translation, sparse cell identifiers,
noise, degenerate cells, and supported wheel installation.

Local Dask integration tests must additionally verify that only one tile is
admitted per worker, the worker has one Dask execution thread, the effective
OpenMP budget is recorded, the worker PID changes after each tile, and
allocator-held RSS cannot accumulate across tile attempts. Cancellation tests
must prove that the worker control path remains responsive while the sole task
thread is inside the GIL-released binding, reaches the correct active-attempt
source, and yields no committed partial result. Failure-injection fixtures must
cover a normal exception, injected worker loss, cooperative cancellation, and a
simulated non-responsive native operation that exceeds its grace period and is
replaced by a fresh Nanny worker. Each retry must use a new attempt identity and
directory. The same fixture must produce identical tile outputs through the
locally owned cluster and a compatible externally supplied client.

Multi-worker integration tests must use heterogeneous tile durations and verify
that the bounded sliding window keeps other workers productive while one worker
is recycled, never admits two simultaneous tiles to one worker, respects both
CPU and aggregate-memory admission limits, and produces the same committed
manifest and scientific result regardless of task-completion order.

### Actual-data validation

Run the untiled reference at least twice under a locked native seed, parameters,
build, and execution settings before judging tiled quality. This measures
Baysor's ordinary stochastic or parallel-execution variability.
Tiled-versus-untiled disagreement must be interpreted relative to this
untiled-versus-untiled baseline rather than against an unrealistic requirement
of bitwise identity.

Compare the tiled result to the untiled UCB reference using:

- retained and assigned transcript counts;
- noise fraction;
- molecule-assignment agreement after cell matching;
- per-cell transcript-set Jaccard similarity;
- total and substantial cell counts;
- molecules per cell and cell-by-gene count correlation;
- one-nucleus-per-cell ownership;
- split nuclei and multi-nucleus cells;
- cell area and shape distributions;
- known-marker coherence;
- unmatched cells and substantial split/merge event rates;
- per-tile confidence distributions, signal/noise fit parameters, graph
  thresholds, and noise fractions;
- disagreement as a function of distance to the nearest seam; and
- visual overlays at seams and in representative tissue regions.

Run at least the proposed `2 * scale`, `4 * scale`, and `6 * scale` halo
experiments, and repeat the selected tiled run after shifting the grid by half a
core in x and y. A scientifically reliable tiled implementation should be
nearly invariant to grid placement, and its residual disagreement in well
contextualized interiors should be close to the untiled repeat baseline.

## Production acceptance gates

The exact biological thresholds should be finalized after the untiled
actual-data run. The following engineering gates are mandatory:

- 100% of retained transcript IDs have exactly one final output row;
- no unknown, duplicated, or uncovered transcript IDs;
- no global cell contains two local cells from the same tile;
- all shapes and table instances refer to existing global cell IDs;
- matched-cell count profiles correlate with the untiled reference at greater
  than 0.99;
- tiled-versus-untiled assignment disagreement in well-contextualized interiors
  does not materially exceed the untiled-versus-untiled baseline;
- confidence, graph-threshold, and noise summaries show no unexplained
  tile-boundary discontinuities or tile-wide shifts;
- nucleus-ownership metrics do not materially degrade relative to untiled;
- seam-proximal disagreement is not materially worse than interior
  disagreement; and
- shifting the tile grid does not materially change the result;
- the observed Dask-worker and OpenMP configuration stays within the declared
  CPU and memory budgets; and
- cancelled or failed tile attempts leave no committed result, while every
  subsequent attempt starts in a fresh compatible worker process.

Failure of a quality gate must make the result ineligible for production use; it
must not be reduced to a warning hidden in logs.

## Principal risks

### Loss of global statistical context

Global prefiltering, one feature encoding, and fixed global `scale` and
`scale_std` remove several avoidable sources of drift. Optional molecule
clustering is disabled in the production baseline because fitting it
independently per tile would introduce another invocation-wide model.

The principal remaining risk is the current CLI's per-invocation confidence and
noise fit. Halos do not make that model global, so it can shift assignments
throughout a biologically atypical tile rather than only near its boundary.
Per-invocation graph thresholds and stochastic component initialization are
secondary sources of difference. The tiled implementation must retain and
compare these statistics and preserve global confidence or apply one exported
global calibration. If controlled tile runs still exceed the untiled
repeat baseline, the next steps are globally fixed graph calibration and larger
or adaptive tiles. Failure after those mitigations is the trigger to investigate
native domain decomposition rather than silently weakening the quality gates.

Global or consistently transferred clustering may enter the supported scope only
if untiled crop experiments show a meaningful scientific benefit and the
resulting implementation passes the same production gates.

### Native memory and Dask worker lifecycle

Baysor's native working set is unmanaged memory from Dask's perspective. Dask
spilling therefore cannot make room inside an executing segmentation, and an
allocator may retain freed pages after a successful call. Relying only on Nanny
memory termination would turn ordinary allocator retention into emergency worker
loss and task rescheduling.

Admission control must use measured peak tile memory with a safety margin. The
normal path proactively recycles a dedicated Nanny-supervised worker after each
attempt; the Nanny threshold is a last-resort guard. Tests must distinguish a
bounded in-process allocator plateau from linear growth, and the run manifest
must make worker loss, cancellation, retry, and output promotion auditable.

The worker remains under Dask's lifecycle control while C++ executes in process:
the Future is visible as processing and the Nanny observes process RSS and
liveness, but Dask cannot inspect or forcibly stop the native stack. The released
GIL permits a targeted worker control request to set the registered native
cancellation source. If cooperative checkpoints do not return the call before
the grace deadline, replacing the exclusive worker process is the hard stop.

### Oversubscription and unsafe shared-cluster cancellation

Dask counts its own task-execution threads, not the OpenMP team created inside a
Baysor call. Incorrect worker provisioning can therefore oversubscribe a host
even when Dask reports no excess tasks. The initial one-worker policy and the
`W * T <= C` contract are mandatory until resource benchmarks validate another
configuration. Other native libraries may maintain their own thread pools, so
the effective process and thread count must be measured rather than inferred
only from Dask and OpenMP settings. Deployment-level CPU affinity is the hard
enforcement mechanism where strict isolation is required.

Hard cancellation restarts a worker process and is unsafe when that worker owns
unrelated tasks or irreplaceable data. Caller-managed clients must provide
exclusive restartable Baysor workers. A shared cluster that cannot satisfy this
contract is unsupported rather than accepted with weaker cancellation or memory
guarantees.

### Incorrect cell merging

Aggressive transitive reconciliation can merge neighboring cells. Stable shared
transcript evidence, the one-cell-per-tile component invariant, prior-nucleus
evidence, and explicit ambiguity reporting are required safeguards.

### Insufficient halo

Cells or relevant molecule neighborhoods larger than the halo can remain
truncated. Halo-size experiments, boundary-touch flags, and grid-shift
validation are required before selecting a default.

### Excessive intermediate data

Halo duplication is modest for the proposed UCB grid, but tile inputs and full
tile output bundles can still consume substantial disk. Checksummed resume,
selective upstream output, a documented retention policy, and capacity
preflight are required.

### Undocumented Baysor output assumptions

Row-order preservation is visible in the pinned implementation but is not an
adequate long-term identity contract. Parquet transcript-ID round-tripping is a
production requirement.

### `baysor_python` and Baysor version drift

A binding release can silently diverge from its reference CLI if source revisions,
parameter defaults, serializers, or orchestration are not coupled. Every release
must compile and report one Baysor revision, require the CLI and binding to call
the same C++ run operation, and test Python-versus-CLI parity on that revision.

## Decision

Proceed in the following order:

1. pin Baysor and extract one coarse-grained public C++ segmentation operation;
2. make the Baysor CLI a thin frontend over that operation and prove parity with
   the untouched pinned CLI;
3. expose the shared C++ operation through a thin nanobind module;
4. establish the typed `baysor_python.segment(...)` API and prove native
   Python-versus-CLI parity on the small fixture;
5. add the SpatialData preparation and import adapter and repeat parity through
   that complete path;
6. complete the supported source distribution and platform-wheel release matrix;
7. expose Baysor boundary estimation as an array-oriented native
   `baysor_python` API and verify CLI parity;
8. implement optional Python-orchestrated core-plus-halo tiling as the reusable
   Dask-backed `baysor_python.segment_tiled(...)` workflow, while keeping its tile
   and result contracts independent of Dask collections;
9. implement shared-transcript reconciliation, conservative rescue, and final
   boundary estimation inside that `baysor_python` workflow;
10. provide a locally owned, Nanny-supervised Dask cluster for native tile calls,
    first validate one worker and one tile at a time, then add bounded multi-worker
    execution with explicit CPU and memory admission, accept a compatible
    caller-managed client, recycle workers after tile attempts, and construct
    reconciled SpatialData results;
11. establish the deferred untiled actual-UCB reference after Phase 1 and no later
   than immediately before tiled-quality validation; and
12. promote tiled mode only after it matches the untiled reference and passes
   untiled-repeat, halo, seam, and grid-shift validation.

Tiling is therefore a planned scalability capability, not the default for the
current 47-million-transcript sample. The architecture must nevertheless make
it possible to scale beyond a single machine without changing the authoritative
data model or scientific interpretation of the result.

Until globally consistent confidence/noise handling and the actual-data quality
gates are satisfied, tiled output is not production-supported and untiled Baysor
remains the scientific reference.

The selected implementation is explicitly Python-orchestrated tiling with native
Baysor jobs, not a C++ loop around independent tiles. Native domain decomposition
within Baysor is the defined contingency only if the validated
overlap-and-reconciliation approach cannot meet the scientific quality gates.
