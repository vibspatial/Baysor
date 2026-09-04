#pragma once

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/segmentation/cancellation.h"
#include "baysor/utils/options.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace baysor {

inline constexpr std::uint64_t kDefaultSegmentationSeed = 1;
inline constexpr std::uint32_t kRandomSubstreamContractVersion = 1;

/// Stable random-stream categories used by the public segmentation contract.
///
/// Contract version 1 preserves the pre-extraction, one-thread defaults when
/// the master seed is 1: the shared core scientific stream uses seed 1 and the
/// formerly fixed clustering, neighborhood-composition, and diagnostic streams
/// use seed 42. Other master seeds are deterministically expanded with the
/// SplitMix64 finalizer over
///
///   master_seed + 0x9e3779b97f4a7c15 * (stream_id + 1)
///
/// The numeric stream identifiers and derivation must not change without
/// incrementing kRandomSubstreamContractVersion.
enum class RandomSubstream : std::uint64_t {
    CoreScientific = 0,
    MoleculeClustering = 1,
    NeighborhoodComposition = 2,
    Diagnostics = 3
};

[[nodiscard]] std::uint64_t derive_random_substream_seed(
    std::uint64_t master_seed,
    RandomSubstream stream
) noexcept;

[[nodiscard]] const char* random_substream_name(RandomSubstream stream) noexcept;

/// A path-oriented molecule input and its read-time interpretation.
struct MoleculeInputSpecification {
    std::string path;
    MoleculeInputOptions options;
};

/// Scientific neighborhood-composition settings currently mixed into the CLI's
/// plotting options. Presentation-only settings intentionally do not enter the
/// segmentation request.
struct NeighborhoodCompositionOptions {
    int neighborhood_size = 0;
    std::string method = "ri";
};

/// Requested native execution resources.
///
/// native_threads controls only Baysor's OpenMP team-size request. A positive
/// value requests that many threads and zero inherits the caller's configured
/// OpenMP maximum. use_arrow_threads independently controls parallel CSV and
/// Parquet decoding; it does not size Arrow's process-global CPU pool. Callers
/// requiring deterministic semantic replay must request exactly one OpenMP
/// thread; parallel execution is not guaranteed to be deterministic, even with
/// an unchanged seed and thread count.
struct SegmentationExecutionOptions {
    int native_threads = 0;
    bool use_arrow_threads = true;
};

/// Products requested from one segmentation run.
///
/// The default requests the complete scientific result. Selective
/// materialization is a later optimization; until then an implementation may
/// produce more fields and must describe them in produced_products.
struct SegmentationProducts {
    bool molecule_assignments = true;
    bool molecule_confidence = true;
    bool assignment_confidence = true;
    bool molecule_clusters = true;
    bool neighborhood_composition_colors = true;
    bool cell_statistics = true;
    bool boundaries = true;
    bool count_matrix = true;
    bool diagnostics = true;

    [[nodiscard]] static SegmentationProducts none() noexcept {
        SegmentationProducts products;
        products.molecule_assignments = false;
        products.molecule_confidence = false;
        products.assignment_confidence = false;
        products.molecule_clusters = false;
        products.neighborhood_composition_colors = false;
        products.cell_statistics = false;
        products.boundaries = false;
        products.count_matrix = false;
        products.diagnostics = false;
        return products;
    }
};

/// Owning description of a single segmentation run.
///
/// Values that depend on loaded data, such as an automatic scale or component
/// count, remain unresolved here and are returned through resolved_options.
struct SegmentationRequest {
    MoleculeInputSpecification molecules;
    PriorInputOptions prior;
    SegmentationOptions segmentation;
    NeighborhoodCompositionOptions neighborhood_composition;
    SegmentationProducts requested_products;
    /// Master seed for the versioned, run-local random streams. Together with
    /// one-thread execution, unchanged input, options, build, and platform, this
    /// guarantees repeatable semantic results. Parallel runs are not guaranteed
    /// deterministic.
    std::uint64_t random_seed = kDefaultSegmentationSeed;
    SegmentationExecutionOptions execution;
};

/// Rectangular, named per-cell statistics table.
struct CellStatistics {
    std::vector<std::string> cell_ids;
    std::vector<std::string> columns;
    Eigen::MatrixXd values; ///< Rows follow cell_ids; columns follow columns.
};

/// Sparse cell-by-gene count matrix with owned axis labels.
struct CellByGeneCounts {
    std::vector<std::string> cell_ids;
    std::vector<std::string> gene_names;
    Eigen::SparseMatrix<float> values; ///< Rows are cells; columns are genes.
};

/// Report-neutral convergence values retained for molecule clustering.
struct MoleculeClusteringDiagnostics {
    std::vector<double> max_differences;
    std::vector<double> assignment_change_fractions;
};

/// Fitted signal/noise-model values required by native diagnostics.
struct ConfidenceEstimationDiagnostics {
    std::vector<double> edge_lengths;
    std::vector<double> fit_differences;
    int neighbor_index = 0;
    double signal_mean = 0.0;
    double signal_standard_deviation = 0.0;
    double noise_mean = 0.0;
    double noise_standard_deviation = 0.0;
};

/// Number of retained cells at configured minimum-molecule thresholds for one
/// BMM iteration.
struct ComponentCountSnapshot {
    std::unordered_map<int, int> cells_by_minimum_molecule_count;
};

/// Report-neutral neighborhood-composition embedding diagnostics.
struct NeighborhoodCompositionDiagnostics {
    std::vector<int> sample_molecule_indices;
    std::vector<double> sample_embedding_x;
    std::vector<double> sample_embedding_y;
    double chosen_confidence_threshold = 0.95;
    int anchor_count = 0;
};

/// Owned diagnostics needed to reproduce native reports without retaining
/// BmmData or another working object.
struct SegmentationDiagnostics {
    std::optional<ConfidenceEstimationDiagnostics> confidence_estimation;
    std::vector<ComponentCountSnapshot> component_count_trace;
    std::optional<MoleculeClusteringDiagnostics> molecule_clustering;
    std::optional<NeighborhoodCompositionDiagnostics> neighborhood_composition;
};

/// Effective scientific and execution options after data-dependent resolution.
struct ResolvedSegmentationOptions {
    MoleculeInputOptions molecules;
    PriorInputOptions prior;
    SegmentationOptions segmentation;
    NeighborhoodCompositionOptions neighborhood_composition;
    SegmentationExecutionOptions execution;
};

struct RandomSubstreamProvenance {
    RandomSubstream stream = RandomSubstream::CoreScientific;
    std::uint64_t seed = kDefaultSegmentationSeed;
};

/// Native implementation and random-state provenance for one completed run.
struct NativeRunProvenance {
    std::string baysor_version;
    std::string build_revision;
    std::uint64_t random_seed = kDefaultSegmentationSeed;
    std::uint32_t random_substream_contract_version = kRandomSubstreamContractVersion;
    std::vector<RandomSubstreamProvenance> random_substreams;
    int requested_native_threads = 0;
    int configured_openmp_max_threads = 0;
    bool openmp_dynamic_enabled = false;
    bool arrow_threads_enabled = true;
};

/// Fully owned scientific result.
///
/// Empty optional products and produced_products distinguish an omitted product
/// from a valid product that happens to contain no rows. No member refers to
/// BmmData or to CLI-owned state.
struct SegmentationResult {
    MoleculeData molecules;
    std::vector<int> cell_assignments; ///< 0 is noise; positive value i maps to cell_ids[i - 1].
    std::vector<double> assignment_confidence; ///< One value per retained molecule when produced.
    std::vector<int> molecule_clusters; ///< One-based cluster IDs, one per molecule when produced.
    std::vector<std::string> neighborhood_composition_colors; ///< One value per molecule when produced.
    std::vector<std::string> cell_ids; ///< Stable IDs ordered by positive assignment value.
    std::optional<CellStatistics> cell_statistics;
    std::optional<PolygonCollection> boundaries_2d;
    std::optional<PolygonStack> boundaries_3d;
    std::optional<CellByGeneCounts> count_matrix;
    std::optional<SegmentationDiagnostics> diagnostics;
    ResolvedSegmentationOptions resolved_options;
    SegmentationProducts produced_products = SegmentationProducts::none();
    NativeRunProvenance provenance;
};

/// Distinct cooperative-cancellation outcome. A cancelled run never returns a
/// partially valid SegmentationResult.
struct SegmentationCancelled {};

using SegmentationOutcome = std::variant<SegmentationResult, SegmentationCancelled>;

enum class SegmentationErrorCode {
    InvalidRequest,
    UnsupportedExecutionContext,
    MoleculeInput,
    PriorInput,
    NativeProcessing,
    Serialization
};

/// Structured native failure. Cancellation is represented by
/// SegmentationCancelled rather than this exception type.
class SegmentationError : public std::runtime_error {
public:
    SegmentationError(SegmentationErrorCode code, std::string message);

    [[nodiscard]] SegmentationErrorCode code() const noexcept;

private:
    SegmentationErrorCode code_;
};

/// Execute one complete native segmentation.
///
/// The operation loads the prepared molecule input, resolves data-dependent
/// options, runs the scientific workflow, and returns owned result values. It
/// does not choose filenames or write output artifacts. Failures throw
/// SegmentationError; cooperative cancellation returns SegmentationCancelled.
/// Calls must originate outside an active OpenMP parallel region. Sequential
/// calls in one process are supported; overlapping calls in one process are
/// rejected with UnsupportedExecutionContext.
[[nodiscard]] SegmentationOutcome run_segmentation(
    const SegmentationRequest& request,
    const CancellationToken& cancellation
);

} // namespace baysor
