# Native segmentation API contract

This document describes the public C++ segmentation boundary in this repository
and how [`baysor-python`](https://github.com/vibspatial/baysor-python) consumes
it. The declarations in
[`include/baysor/segmentation`](include/baysor/segmentation) are the source of
truth; this document explains their intended use and ownership.

## Status

The public native types and `baysor::run_segmentation(...)` operation below are
implemented and tested. The operation executes the complete scientific workflow
directly through `baysor_lib`; it does not invoke the CLI or write output files.
This repository does not provide a Python extension; the nanobind adapter and
public Python API belong to the separate `baysor-python` distribution.

This is a versioned C++17 source-level API designed to be shared by the CLI and
native-language bindings. It is not a frozen binary ABI.

## Entry point and exposed layers

The natural entry point for a Python caller is the public
`baysor_python.segment(...)` function. Users, Harpy adapters, and Dask
tile tasks do not invoke the Baysor executable, import the private nanobind
module directly, or construct Baysor's internal BMM objects.

```text
Python caller: user / Harpy adapter / Dask tile task
                              |
                              | public Python API
                              v
              baysor_python.segment(...)
                              |
                              | private nanobind adapter
                              | - convert typed Python values
                              | - retain the CancellationSource
                              | - release the GIL around native work
                              | - translate outcomes and errors
                              v
     baysor_python._baysor_core.run_segmentation(...)
                              |
                              | construct ordinary C++ contract values
                              v
+---------------- public native API in vibspatial/Baysor ----------------+
|                                                                        |
|          SegmentationRequest + CancellationToken                       |
|                              |                                         |
|                              v                                         |
|              baysor::run_segmentation(...)                             |
|                    /             |              \                       |
|                   v              v               v                      |
|       SegmentationResult  SegmentationCancelled  throws                |
|       (complete, owned)   (no partial result)     SegmentationError     |
|                                                                        |
+------------------------------------------------------------------------+
                              |
                              | reusable native serializers on success
                              v
       validated files on disk + Python SegmentationArtifacts

The current Baysor CLI retains its locked compatibility orchestration while its
frontend is migrated to this API. It is not the implementation used by the
Python caller.
```

The exposed layers are:

| Layer | Contract | Intended audience |
| --- | --- | --- |
| Public Python product API | `baysor_python.segment(...)`, typed Python options, `SegmentationArtifacts`, and documented Python exceptions | Supported entry point for Python callers; owned by `baysor-python` |
| Private Python/native adapter | `baysor_python._baysor_core` and its nanobind conversions, GIL handling, cancellation registry, and exception translation | Implementation detail of `baysor-python`; callers must not depend on its exact signature or binding types |
| Public native source API | `SegmentationRequest`, `CancellationToken`, `SegmentationOutcome`, `SegmentationResult`, `SegmentationCancelled`, `SegmentationError`, and `run_segmentation(...)` in namespace `baysor` | Supported C++17 integration surface for native consumers and the nanobind module; the compatibility CLI is being migrated to the same surface |
| Native internals | `BmmData`, graphs, clustering and estimator working objects, CLI arguments, logger configuration, filenames, and process exit codes | Private implementation details that may evolve behind the public operation |

## Native operation

The sole coarse-grained scientific entry point is:

```cpp
[[nodiscard]] baysor::SegmentationOutcome baysor::run_segmentation(
    const baysor::SegmentationRequest& request,
    const baysor::CancellationToken& cancellation
);
```

It accepts exactly two parameters:

1. `request` is an owning, typed description of the input, scientific options,
   requested products, seed, and native execution settings.
2. `cancellation` is a read-only view of shared cancellation state. A caller
   that needs cancellation retains the corresponding `CancellationSource` and
   may request cancellation from another thread. A default-constructed token is
   never cancelled.

The operation does not accept CLI arguments, Python objects, an output
directory, filenames, logger configuration, or process-exit policy.

It performs the following scientific sequence using the same native loaders and
algorithms as the CLI: validate and resolve options, load molecules and an
optional prior, estimate molecule confidence, build the molecule graph, perform
optional molecule clustering, run the 2D or 3D BMM, and materialize the requested
scientific products. Values that survive the call are copied or moved into the
owning result; internal graphs, clustering models, and BMM state do not escape.

## `SegmentationRequest`

`SegmentationRequest` contains the following top-level fields:

| Field | Type | Default | Purpose |
| --- | --- | --- | --- |
| `molecules` | `MoleculeInputSpecification` | Empty path plus default molecule options | Prepared CSV or Parquet input and its read-time interpretation |
| `prior` | `PriorInputOptions` | No prior | Optional column, image, or boundary prior |
| `segmentation` | `SegmentationOptions` | See below | Segmentation, clustering, and confidence settings |
| `neighborhood_composition` | `NeighborhoodCompositionOptions` | Size `0`, method `"ri"` | Scientific neighborhood-composition settings |
| `requested_products` | `SegmentationProducts` | Every product requested | Scientific products the caller wants materialized |
| `random_seed` | `std::uint64_t` | `1` | Master seed for the versioned native random-stream contract |
| `execution` | `SegmentationExecutionOptions` | OpenMP threads `0`, Arrow threads enabled | Native resource request |

Sentinel values such as automatic scale or cell-count selection remain
unresolved in the request. `run_segmentation(...)` resolves them once, using the
loaded data, and records the effective values in the result.

### Molecule input

`MoleculeInputSpecification` has two fields:

| Field | Type | Purpose |
| --- | --- | --- |
| `path` | `std::string` | Path to the prepared CSV or Parquet molecule table |
| `options` | `MoleculeInputOptions` | Column mapping, filtering, dimensionality, and confidence-estimation options |

`MoleculeInputOptions` exposes:

| Field | Default |
| --- | --- |
| `x_col`, `y_col`, `z_col` | `"x"`, `"y"`, `"z"` |
| `gene_col`, `qv_col` | `"gene"`, `"qv"` |
| `force_2d` | `false` |
| `min_molecules_per_gene` | `1` |
| `exclude_genes` | Empty comma-separated pattern expression |
| `min_molecules_per_cell` | `0` |
| `confidence_nn_id` | `0` |
| `min_qv` | `-1.0` |
| `x_min`, `y_min`, `z_min` | Negative infinity |
| `x_max`, `y_max`, `z_max` | Positive infinity |

The initial native boundary is deliberately path-oriented. It does not accept a
NumPy array, pandas DataFrame, PyArrow table, or Dask collection.

### Prior input

`PriorInputOptions` exposes:

| Field | Type | Default |
| --- | --- | --- |
| `type` | `PriorInputType` | `None`; alternatives are `Column`, `Image`, and `Boundary` |
| `path` | `std::string` | Empty |
| `column_name` | `std::string` | Empty |
| `unassigned_label` | `std::string` | `"0"` |
| `min_molecules_per_segment` | `int` | `0` |
| `estimate_scale_from_prior` | `bool` | `true` |

### Segmentation options

`SegmentationOptions` exposes:

| Field | Default |
| --- | --- |
| `scale` | `-1.0`, meaning automatic resolution |
| `scale_std` | `"25%"` |
| `cluster_method` | `ClusterMethod::Mrf`; alternatives are `None`, `Louvain`, and `Leiden` |
| `n_clusters` | `0` |
| `cluster_resolution` | `1.0` |
| `cluster_graph_k` | `15` |
| `cluster_n_dims` | `20` |
| `cluster_basis_sample_size` | `100000` |
| `prior_segmentation_confidence` | `0.2` |
| `iters` | `500` |
| `tol` | `0.0`, meaning run all requested iterations |
| `n_cells_init` | `0` |
| `nuclei_genes` | Empty |
| `cyto_genes` | Empty |

### Neighborhood composition and execution

`NeighborhoodCompositionOptions` contains `neighborhood_size` with default `0`
and `method` with default `"ri"`.

`SegmentationExecutionOptions::native_threads` controls Baysor's call-scoped
OpenMP maximum. A positive value requests that number of threads, `0` inherits
the caller's configured OpenMP maximum, and a negative value is invalid. It is
not a total process thread cap: the call applies this request through the OpenMP
runtime and does not resize Arrow's separate process-global thread pools.

`use_arrow_threads` independently enables or disables Arrow's parallel CSV and
Parquet decoding and defaults to `true`. It is a Boolean decoding policy, not an
Arrow thread-count setting. When enabled, the available concurrency is governed
by Arrow's process-level configuration; when disabled, the covered readers use
their serial decoding paths. Neither setting promises that the process contains
no other library or background threads.

### Requested products

`SegmentationProducts` contains these Boolean fields, all `true` in a default
request:

- `molecule_assignments`
- `molecule_confidence`
- `assignment_confidence`
- `molecule_clusters`
- `neighborhood_composition_colors`
- `cell_statistics`
- `boundaries`
- `count_matrix`
- `diagnostics`

Until the selective-materialization work is implemented, the engine may produce
more than requested. The result's `produced_products` field is authoritative
about what is actually present.

## Randomness and native threads

Each call constructs its core scientific random generator from `random_seed`.
The generator is passed explicitly through duplicate-coordinate jitter, graph
construction, single-thread stochastic BMM assignment, and boundary estimation.
Molecule clustering, neighborhood-composition colours, and diagnostic embeddings
receive deterministic seeds derived from the same master seed through the
versioned `RandomSubstream` contract. The default master seed `1` maps the core
stream to seed `1` and the other established subsystem streams to seed `42`,
preserving the one-thread compatibility baseline.

A positive `native_threads` value configures the OpenMP maximum for the duration
of the call; zero inherits the caller's current value. The previous OpenMP
maximum and dynamic-team policy are restored on success, cancellation, and
failure. Baysor inherits dynamic-team adjustment rather than disabling it.

Provenance deliberately does not claim to record an effective or observed team
size. `requested_native_threads` is the request (`0` still means inherit),
`configured_openmp_max_threads` is the maximum reported by OpenMP after applying
that request, and `openmp_dynamic_enabled` records whether OpenMP may use a
smaller team. `arrow_threads_enabled` records the independent I/O policy.

The public reproducibility contract is:

| Execution | Guarantee |
| --- | --- |
| One native thread, with the same prepared input, options, Baysor and dependency builds, platform, and master seed | Repeated calls, including calls made in the same process, produce the same semantic scientific result. |
| More than one native thread, with the same seed and thread count | Baysor initializes the same versioned random streams, but does not guarantee an identical result. Dynamic scheduling can associate work with different streams, and parallel floating-point evaluation order can vary. |
| Different native thread counts | Results may differ. |
| Different builds or platforms | Bitwise identity is not guaranteed; compatibility must be assessed using the semantic output contract and declared numerical tolerances. |

Here, the same *semantic scientific result* means the same retained molecule
identities and the same noise/cell partition, allowing a consistent relabelling
of otherwise equivalent cell identifiers, with floating-point products compared
using their declared tolerances. It does not mean byte-identical output files,
identical serialization order, or cross-platform bitwise equality.

Callers that require deterministic replay must explicitly set
`native_threads = 1`; leaving it at `0` does not provide that guarantee because
the runtime may select multiple threads. Every completed run records the master
seed, random-stream contract version, derived stream seeds, requested and
configured OpenMP settings, dynamic-team policy, and Arrow I/O policy in
`SegmentationResult::provenance`. Consumers comparing tiled or distributed runs
should lock and record both the seed and execution configuration.

### Entry context and concurrency

`run_segmentation(...)` must be entered from outside an active OpenMP parallel
region. A nested call is rejected before Baysor changes the caller's OpenMP
configuration. Sequential calls in one process are supported; concurrent calls
in the same process are rejected immediately with
`UnsupportedExecutionContext`. This fail-fast rule prevents concurrent calls
from racing over OpenMP configuration or silently multiplying native teams.

An embedding worker may therefore execute at most one Baysor segmentation at a
time. Distributed orchestration may choose one process with a large OpenMP team,
many processes with one OpenMP thread each, or a bounded hybrid. For CPU compute
threads it should maintain:

```text
concurrent worker processes * OpenMP threads per worker <= allocated CPU cores
```

Arrow reader concurrency and per-tile memory are budgeted separately. In
particular, changing `native_threads` for a call does not resize an Arrow pool
that the worker process has already initialized. Process-level deployment must
therefore configure Arrow separately or set `use_arrow_threads = false` when
parallel decoding is not part of its resource budget. This model is compatible
with a Dask nanny supervising one worker process: cancellation is cooperative
within the process, while forced termination and guaranteed memory reclamation
remain supervisor responsibilities.

## Return and failure contract

`SegmentationOutcome` is:

```cpp
using SegmentationOutcome =
    std::variant<SegmentationResult, SegmentationCancelled>;
```

There are three mutually exclusive outcomes:

- success returns one complete, owned `SegmentationResult`;
- cooperative cancellation returns `SegmentationCancelled`; and
- failure throws `SegmentationError` with a structured error code and message.

Cancellation is not an exception and is never represented as a partially valid
result. The operation observes cancellation before work, between major
scientific phases, and at complete BMM iteration boundaries. It never interrupts
an E-step/M-step transition halfway through. Cancellation is cooperative, so a
request made while one indivisible native primitive is running takes effect when
that primitive reaches its next checkpoint. Callers must not assume an
instantaneous response or a fixed maximum cancellation latency. Complete BMM
iterations are the current safe checkpoint granularity; coordinated checkpoints
may be added inside long-running phases if measurements show that this latency
is unacceptable, provided all participating native workers leave and join
normally. Serializers accept only a completed `SegmentationResult`, so a
cancelled operation cannot serialize a partial run.

### `SegmentationResult`

The success value owns all data needed after the algorithm's working objects
have been destroyed:

| Field | Type | Contents |
| --- | --- | --- |
| `molecules` | `MoleculeData` | Retained coordinates, gene encoding, confidence and source metadata |
| `cell_assignments` | `std::vector<int>` | One assignment per retained molecule; `0` is noise and positive value `i` maps to `cell_ids[i - 1]` |
| `assignment_confidence` | `std::vector<double>` | Per-molecule assignment confidence when produced |
| `molecule_clusters` | `std::vector<int>` | One-based per-molecule cluster identifiers when produced |
| `neighborhood_composition_colors` | `std::vector<std::string>` | Per-molecule neighborhood-composition colors when produced |
| `cell_ids` | `std::vector<std::string>` | Stable cell identifiers ordered by positive assignment value |
| `cell_statistics` | `std::optional<CellStatistics>` | Named dense per-cell statistics table |
| `boundaries_2d` | `std::optional<PolygonCollection>` | Cell ID to a `2 x N` boundary-vertex matrix |
| `boundaries_3d` | `std::optional<PolygonStack>` | Per-slice collections of cell boundary polygons |
| `count_matrix` | `std::optional<CellByGeneCounts>` | Sparse cell-by-gene count matrix and owned axis labels |
| `diagnostics` | `std::optional<SegmentationDiagnostics>` | Requested report-neutral convergence and model diagnostics |
| `resolved_options` | `ResolvedSegmentationOptions` | Effective input, prior, scientific, neighborhood, and execution settings after resolution |
| `produced_products` | `SegmentationProducts` | Flags that distinguish an omitted product from a valid empty product |
| `provenance` | `NativeRunProvenance` | Native version, build, random-stream, seed, and thread provenance |

`MoleculeData` owns:

- coordinate arrays `x`, `y`, and optional `z`;
- one-based encoded `gene` values and their `gene_names` dictionary;
- molecule `confidence`, internal `cluster`, `prior_segmentation`, and
  `nuclei_probs` arrays; and
- preserved `source_transcript_id` values when present in the input.

The nested result objects are:

- `CellStatistics`: `cell_ids`, column names in `columns`, and a dense
  row-by-column `values` matrix;
- `CellByGeneCounts`: `cell_ids`, `gene_names`, and a sparse `values` matrix
  whose rows are cells and columns are genes;
- `SegmentationDiagnostics`: optional confidence-estimation diagnostics, a
  component-count trace, optional molecule-clustering diagnostics, and optional
  neighborhood-composition diagnostics;
- `ResolvedSegmentationOptions`: resolved molecule, prior, segmentation,
  neighborhood-composition, and execution options; and
- `NativeRunProvenance`: `baysor_version`, `build_revision`, master
  `random_seed`, random-substream contract version and derived streams, and
  requested/configured OpenMP, dynamic-team, and Arrow I/O settings.

For completeness, the diagnostic and provenance structures expose these
fields:

| Type | Fields |
| --- | --- |
| `ConfidenceEstimationDiagnostics` | `edge_lengths`, `fit_differences`, `neighbor_index`, `signal_mean`, `signal_standard_deviation`, `noise_mean`, `noise_standard_deviation` |
| `ComponentCountSnapshot` | `cells_by_minimum_molecule_count` |
| `MoleculeClusteringDiagnostics` | `max_differences`, `assignment_change_fractions` |
| `NeighborhoodCompositionDiagnostics` | `sample_molecule_indices`, `sample_embedding_x`, `sample_embedding_y`, `chosen_confidence_threshold`, `anchor_count` |
| `SegmentationDiagnostics` | `confidence_estimation`, `component_count_trace`, `molecule_clustering`, `neighborhood_composition` |
| `RandomSubstreamProvenance` | `stream`, `seed` |
| `NativeRunProvenance` | `baysor_version`, `build_revision`, `random_seed`, `random_substream_contract_version`, `random_substreams`, `requested_native_threads`, `configured_openmp_max_threads`, `openmp_dynamic_enabled`, `arrow_threads_enabled` |

No result field points into a temporary `BmmData` or CLI-owned object. Normal
C++ ownership and move semantics apply.

### `SegmentationCancelled`

`SegmentationCancelled` is intentionally an empty marker type. It has no result
attributes because cancellation must not publish incomplete scientific data.
The caller can distinguish it from success with the normal `std::variant`
facilities.

### `SegmentationError`

`SegmentationError` derives from `std::runtime_error`. Its `what()` value carries
the message and `code()` returns one of:

| Error code | Category |
| --- | --- |
| `InvalidRequest` | Invalid or inconsistent typed input |
| `UnsupportedExecutionContext` | Nested OpenMP entry or another active in-process segmentation call |
| `MoleculeInput` | Molecule input loading or validation failure |
| `PriorInput` | Prior input loading or validation failure |
| `NativeProcessing` | Failure in native scientific processing |
| `Serialization` | Failure while serializing a completed result |

The CLI maps structured errors to messages and exit codes. The nanobind
adapter maps them to documented Python exceptions without reducing them to an
integer process status.

## CMake embedding

The native library is position-independent and exported to build-tree consumers
as `baysor::baysor`. A parent CMake project can consume an exact Baysor checkout
without building the CLI:

```cmake
add_subdirectory(path/to/Baysor baysor EXCLUDE_FROM_ALL)
target_link_libraries(my_extension PRIVATE baysor::baysor)
```

The target publishes the C++17 requirement, public include directory, and
transitive link requirements needed by a static consumer. The standalone CLI is
built by default only when Baysor is the top-level project; an embedding parent
can opt in with `BAYSOR_BUILD_CLI=ON`. The clean consumer under
`tests/cmake_consumer` is the executable link check for this contract.

Returned results own their values and release them through normal C++ lifetime
rules. Repeated calls retain no Baysor run object or random state. A native
allocator may retain freed arenas, so stable resident memory is not evidence of
a live-object leak; use a sanitizer or platform leak checker to detect linear
growth in live allocations. A caller that requires hard memory reclamation must
isolate the call in a supervised process.

## Paths, native objects, and Python artifacts

The input path belongs to `SegmentationRequest` because it is the scalable
transport into the native operation. The output directory does not belong to
that request. On success, `run_segmentation(...)` returns an in-memory native
`SegmentationResult`; a frontend then passes it to reusable native serializers.

[`include/baysor/reporting/output.h`](include/baysor/reporting/output.h) provides
result-based overloads for the molecule table (CSV or Parquet), cell statistics
(CSV or Parquet), count matrix (TSV), and 2D boundaries (GeoJSON). They require
the corresponding product to be present and report failures with
`SegmentationErrorCode::Serialization`. File naming and output-directory policy
remain the frontend's responsibility.

The public Python call returns a small `SegmentationArtifacts`
descriptor containing validated output paths, status, resource measurements,
resolved settings, and provenance. It does not send the complete native result
object graph through Python or the Dask scheduler.

This separation gives each layer one responsibility:

```text
prepared CSV or Parquet path
            |
            v
SegmentationRequest -> run_segmentation(...) -> SegmentationResult
                                                  |
                                                  v
                                      reusable native serializers
                                                  |
                                                  v
                                  files + SegmentationArtifacts
```
